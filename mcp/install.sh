#!/usr/bin/env bash
# Add this repo's MCP servers to an already-installed Loom.
#
#   curl -fsSL https://raw.githubusercontent.com/NIAID-BRC-Codeathons/hypothesis2omics/main/mcp/install.sh | bash
#
# Clones (or updates) the repo, then hands off to mcp/register.mjs, which merges
# the servers into loom's ~/.pi/agent/mcp.json. Flags are passed through:
#   ... | bash -s -- --dry-run
#   ... | bash -s -- --remove
#
# Override the clone location with H2O_MCP_DIR, or the branch with H2O_MCP_BRANCH.
set -euo pipefail

REPO="https://github.com/NIAID-BRC-Codeathons/hypothesis2omics.git"
BRANCH="${H2O_MCP_BRANCH:-main}"
DIR="${H2O_MCP_DIR:-$HOME/.loom/mcp/hypothesis2omics}"

need() {
  command -v "$1" >/dev/null 2>&1 && return 0
  echo "error: \`$1\` not found on your PATH." >&2
  echo "  $2" >&2
  exit 1
}

need git "Install git, then re-run this script."
need node "Node is what runs the registration step. It comes with Loom
  (npm install -g @galaxyproject/loom) -- if Loom isn't installed yet, start there."

# uv is only needed when the agent launches a server, not to register one. Warn,
# don't block: registering now and installing uv later works fine.
if ! command -v uv >/dev/null 2>&1; then
  cat >&2 <<'NOTE'
warning: `uv` was not found on your PATH.

The MCP servers run as `uv run server.py`, so their tools will fail to start
until you install uv:

  * Install script:  curl -LsSf https://astral.sh/uv/install.sh | sh
  * Homebrew:        brew install uv

Registering anyway -- install uv before you use the tools.

NOTE
fi

if [ ! -e "$DIR" ]; then
  echo "Cloning $REPO ($BRANCH) -> $DIR"
  mkdir -p "$(dirname "$DIR")"
  git clone --depth 1 --branch "$BRANCH" "$REPO" "$DIR"
elif [ -d "$DIR/.git" ]; then
  echo "Updating $DIR"
  # --ff-only so local edits are surfaced rather than silently discarded.
  if ! git -C "$DIR" pull --ff-only; then
    echo "error: could not fast-forward $DIR (local commits or edits?)." >&2
    echo "  Resolve them, or delete the directory and re-run this script." >&2
    exit 1
  fi
else
  echo "error: $DIR exists but is not a git clone. Move it aside and re-run." >&2
  exit 1
fi

exec node "$DIR/mcp/register.mjs" "$@"
