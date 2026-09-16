#!/usr/bin/env node
// Register this repo's MCP servers with an installed Loom.
//
// Loom has no user-extension mechanism -- its Pi extensions are hardcoded as -e
// args in loom/bin/loom.js. But it does merge MCP servers: on every launch it
// reads $PI_CODING_AGENT_DIR/mcp.json, sets its own `galaxy` / `brc-analytics`
// keys, and writes the whole object back, leaving foreign keys alone. So we just
// add ours there.
//
// Usage: node register.mjs [--remove] [--dry-run] [--self-test]

import { readdirSync, existsSync, readFileSync, writeFileSync, mkdirSync, chmodSync, statSync, rmSync } from "node:fs";
import { join, dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { homedir, tmpdir } from "node:os";
import assert from "node:assert/strict";

const HERE = dirname(fileURLToPath(import.meta.url));

/** Same resolution as loom/bin/loom.js:157 -- must not drift from it. */
function mcpConfigPath() {
  return join(process.env.PI_CODING_AGENT_DIR || join(homedir(), ".pi", "agent"), "mcp.json");
}

/** Every sibling directory holding a server.py. The directory name becomes the
 * mcp.json key, which pi-mcp-adapter turns into the tool prefix (hyphens ->
 * underscores). New servers in this repo need no edit here. */
function discover(root = HERE) {
  return readdirSync(root, { withFileTypes: true })
    .filter((e) => e.isDirectory() && existsSync(join(root, e.name, "server.py")))
    .map((e) => e.name)
    .sort();
}

/** --directory is load-bearing: without it uv resolves server.py against whatever
 * cwd the agent spawned in. directTools mirrors loom's own galaxy entry. */
function entry(root, name) {
  return {
    command: "uv",
    args: ["run", "--directory", resolve(root, name), "server.py"],
    directTools: true,
  };
}

function apply(cfg, names, root, remove) {
  cfg.mcpServers ??= {};
  for (const name of names) {
    if (remove) delete cfg.mcpServers[name];
    else cfg.mcpServers[name] = entry(root, name);
  }
  return cfg;
}

function read(path) {
  return existsSync(path) ? JSON.parse(readFileSync(path, "utf-8")) : {};
}

/** mcp.json can hold a Galaxy API key and loom keeps it 0600, so we must too. The
 * writeFileSync mode only applies when the file is created -- hence the chmod. */
function write(path, cfg) {
  mkdirSync(dirname(path), { recursive: true });
  writeFileSync(path, JSON.stringify(cfg, null, 2) + "\n", { mode: 0o600 });
  try {
    chmodSync(path, 0o600);
  } catch {
    // best-effort: filesystems without POSIX modes (e.g. some Windows mounts)
  }
}

function main(argv) {
  const remove = argv.includes("--remove");
  const dryRun = argv.includes("--dry-run");
  const path = mcpConfigPath();
  const names = discover();

  if (!names.length) {
    console.error(`No server.py found under ${HERE} -- nothing to register.`);
    process.exit(1);
  }

  const cfg = apply(read(path), names, HERE, remove);

  if (dryRun) {
    console.log(`# would write ${path}`);
    console.log(JSON.stringify(cfg, null, 2));
    return;
  }

  write(path, cfg);
  console.log(`${remove ? "Removed from" : "Registered in"} ${path}:`);
  for (const n of names) console.log(`  ${n}  ->  ${n.replace(/-/g, "_")}_*`);
  if (!remove) console.log("\nRestart loom (or Orbit), then run /mcp to confirm.");
}

// The one thing worth checking: that we merge into loom's file rather than stomp
// it, and that install/remove round-trips exactly.
function selfTest() {
  const dir = join(tmpdir(), `h2o-mcp-test-${process.pid}`);
  mkdirSync(dir, { recursive: true });
  const path = join(dir, "mcp.json");
  const seed = {
    mcpServers: {
      galaxy: { command: "uvx", args: ["galaxy-mcp>=1.9.0"], env: { GALAXY_API_KEY: "secret" } },
      "brc-analytics": { url: "https://dev.brc-analytics.org/api/v1/mcp/", directTools: true },
    },
  };
  write(path, seed);
  const before = readFileSync(path, "utf-8");

  const names = discover();
  assert.ok(names.length > 0, "discover() found no servers");

  write(path, apply(read(path), names, HERE, false));
  const after = read(path);
  assert.deepEqual(after.mcpServers.galaxy, seed.mcpServers.galaxy, "loom's galaxy entry was clobbered");
  assert.deepEqual(after.mcpServers["brc-analytics"], seed.mcpServers["brc-analytics"], "brc-analytics clobbered");
  for (const n of names) {
    assert.equal(after.mcpServers[n].command, "uv");
    assert.ok(existsSync(join(after.mcpServers[n].args[2], "server.py")), `${n}: --directory points nowhere`);
  }

  const once = readFileSync(path, "utf-8");
  write(path, apply(read(path), names, HERE, false));
  assert.equal(readFileSync(path, "utf-8"), once, "not idempotent");

  write(path, apply(read(path), names, HERE, true));
  assert.equal(readFileSync(path, "utf-8"), before, "--remove did not restore the original");

  if (process.platform !== "win32") {
    assert.equal(statSync(path).mode & 0o777, 0o600, "mcp.json is not 0600");
  }
  rmSync(dir, { recursive: true, force: true });
  console.log(`ok  ${names.length} server(s): ${names.join(", ")}`);
}

if (process.argv.includes("--self-test")) selfTest();
else main(process.argv.slice(2));
