#!/usr/bin/env node
// Register this repo's MCP servers with an installed Loom.
//
// Loom has no user-extension mechanism -- its Pi extensions are hardcoded as -e
// args in loom/bin/loom.js. But it does merge MCP servers: on every launch it
// reads $PI_CODING_AGENT_DIR/mcp.json, sets its own `galaxy` / `brc-analytics`
// keys, and writes the whole object back, leaving foreign keys alone. So we just
// add ours there.
//
// A server opts in by dropping an `mcp-server.json` next to its code. That file
// IS the mcp.json entry -- no invented schema -- with two placeholders:
//
//   ${dir}   the directory holding the manifest
//   ${root}  the repo root
//
// plus an optional "name" (defaults to the manifest's directory name) that sets
// the mcp.json key, which pi-mcp-adapter turns into the tool prefix. Any other
// ${VAR} is left alone for pi-mcp-adapter to interpolate from the environment at
// spawn time. Servers therefore live wherever they like -- a self-contained
// PEP-723 script and a package sharing the repo's pyproject both work.
//
// Usage: node register.mjs [--remove] [--dry-run] [--self-test]

import { readdirSync, existsSync, readFileSync, writeFileSync, mkdirSync, chmodSync, statSync, rmSync } from "node:fs";
import { join, dirname, basename, resolve, relative } from "node:path";
import { fileURLToPath } from "node:url";
import { homedir, tmpdir } from "node:os";
import assert from "node:assert/strict";

const HERE = dirname(fileURLToPath(import.meta.url));
const ROOT = dirname(HERE); // mcp/register.mjs -> repo root
const MANIFEST = "mcp-server.json";
const SKIP = new Set([".git", "node_modules", "__pycache__", ".venv", "venv", ".mypy_cache", ".pytest_cache"]);

/** Same resolution as loom/bin/loom.js:157 -- must not drift from it. */
function mcpConfigPath() {
  return join(process.env.PI_CODING_AGENT_DIR || join(homedir(), ".pi", "agent"), "mcp.json");
}

/** Every mcp-server.json in the repo, as {name, dir, entry}. */
function discover(root = ROOT) {
  const found = [];
  const walk = (dir) => {
    for (const e of readdirSync(dir, { withFileTypes: true })) {
      if (e.isDirectory()) {
        if (!SKIP.has(e.name) && !e.name.startsWith(".")) walk(join(dir, e.name));
      } else if (e.name === MANIFEST) {
        found.push(load(dir, root));
      }
    }
  };
  walk(root);
  return found.sort((a, b) => a.name.localeCompare(b.name));
}

function load(dir, root) {
  const path = join(dir, MANIFEST);
  let raw;
  try {
    raw = JSON.parse(readFileSync(path, "utf-8"));
  } catch (err) {
    throw new Error(`${relative(root, path)}: ${err.message}`);
  }
  const { name = basename(dir), ...entry } = raw;
  if (!entry.command && !entry.url) {
    throw new Error(`${relative(root, path)}: needs a "command" (stdio) or "url" (http)`);
  }
  return { name, dir, entry: substitute(entry, { dir, root }) };
}

/** Replace only ${dir} and ${root}; any other ${VAR} belongs to pi-mcp-adapter,
 * which interpolates it from the live process env when it spawns the server. */
function substitute(value, vars) {
  if (typeof value === "string") {
    return value.replace(/\$\{(dir|root)\}/g, (_, k) => resolve(vars[k]));
  }
  if (Array.isArray(value)) return value.map((v) => substitute(v, vars));
  if (value && typeof value === "object") {
    return Object.fromEntries(Object.entries(value).map(([k, v]) => [k, substitute(v, vars)]));
  }
  return value;
}

function apply(cfg, servers, remove) {
  cfg.mcpServers ??= {};
  for (const s of servers) {
    if (remove) delete cfg.mcpServers[s.name];
    else cfg.mcpServers[s.name] = s.entry;
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
  const servers = discover();

  if (!servers.length) {
    console.error(`No ${MANIFEST} found under ${ROOT} -- nothing to register.`);
    process.exit(1);
  }

  const cfg = apply(read(path), servers, remove);

  if (dryRun) {
    console.log(`# would write ${path}`);
    console.log(JSON.stringify(cfg, null, 2));
    return;
  }

  write(path, cfg);
  console.log(`${remove ? "Removed from" : "Registered in"} ${path}:`);
  for (const s of servers) {
    console.log(`  ${s.name.padEnd(20)} ${s.name.replace(/-/g, "_")}_*   (${relative(ROOT, s.dir) || "."})`);
  }
  if (remove) return;

  console.log("\nRestart loom (or Orbit), then run /mcp to confirm.");
  // A server sharing the repo's pyproject pulls the full scientific stack on its
  // first launch. That easily outruns pi-mcp-adapter's startup window, and the
  // server just looks broken -- so pre-warm the environment out of band.
  if (existsSync(join(ROOT, "pyproject.toml"))) {
    console.log(`\nFirst launch installs this repo's Python environment, which is large enough`);
    console.log(`to time out the agent's MCP startup. Warm it once, ahead of time:`);
    console.log(`\n  uv sync --directory ${ROOT}\n`);
  }
}

// The things worth checking: that we merge into loom's file rather than stomp it,
// that install/remove round-trips, and that placeholders resolve the way the two
// very different servers in this repo need them to.
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

  const servers = discover();
  assert.ok(servers.length >= 2, `expected >=2 servers, found ${servers.length}`);

  write(path, apply(read(path), servers, false));
  const after = read(path);
  assert.deepEqual(after.mcpServers.galaxy, seed.mcpServers.galaxy, "loom's galaxy entry was clobbered");
  assert.deepEqual(after.mcpServers["brc-analytics"], seed.mcpServers["brc-analytics"], "brc-analytics clobbered");

  for (const s of servers) {
    const e = after.mcpServers[s.name];
    assert.ok(e, `${s.name} missing from mcp.json`);
    assert.ok(e.command || e.url, `${s.name}: neither command nor url`);
    for (const a of e.args ?? []) {
      assert.ok(!/\$\{(dir|root)\}/.test(a), `${s.name}: unresolved placeholder in ${a}`);
      if (a.startsWith("/")) assert.ok(existsSync(a), `${s.name}: path arg does not exist: ${a}`);
    }
  }

  // ${IMMPORT_API_KEY} must survive untouched -- pi-mcp-adapter resolves it at spawn.
  const h2o = after.mcpServers["hypothesis2omics"];
  assert.equal(h2o?.env?.IMMPORT_API_KEY, "${IMMPORT_API_KEY}", "env placeholder was substituted too early");

  const once = readFileSync(path, "utf-8");
  write(path, apply(read(path), servers, false));
  assert.equal(readFileSync(path, "utf-8"), once, "not idempotent");

  write(path, apply(read(path), servers, true));
  assert.equal(readFileSync(path, "utf-8"), before, "--remove did not restore the original");

  if (process.platform !== "win32") {
    assert.equal(statSync(path).mode & 0o777, 0o600, "mcp.json is not 0600");
  }
  rmSync(dir, { recursive: true, force: true });
  console.log(`ok  ${servers.length} server(s): ${servers.map((s) => s.name).join(", ")}`);
}

if (process.argv.includes("--self-test")) selfTest();
else main(process.argv.slice(2));
