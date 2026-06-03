#!/usr/bin/env bash
# Launch the high-lift tool with ONE command: the RANS backend (rans_server.py, GPU)
# plus the Vite frontend, together. Ctrl-C stops both.
#
#   ./run-highlift.sh [vite-port]      # default 5173
#
# The backend is one-per-machine (GPU, port 8077). If one is already running (e.g. your
# dev copy's), this reuses it instead of starting a second — so a demo worktree on its own
# port shares the same backend. Override the compute venv with FLOW360_VENV_PY.
# Must be run in a normal (non-sandboxed) shell — the GPU solver won't run otherwise.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PORT="${1:-5173}"
PY="${FLOW360_VENV_PY:-/home/qiqi/flexcompute/compute/.venv/bin/python}"

# Start the backend only if nothing is already listening on 8077; stop it on exit.
if (exec 3<>/dev/tcp/127.0.0.1/8077) 2>/dev/null; then
  echo "[run-highlift] backend already up on :8077 — reusing it"
else
  echo "[run-highlift] starting rans_server.py on :8077"
  "$PY" "$HERE/rans/rans_server.py" &
  SERVER_PID=$!
  trap 'kill $SERVER_PID 2>/dev/null' EXIT   # not exec'd below, so this fires on Ctrl-C
fi

# Frontend in the foreground (Ctrl-C ends it, then the trap stops the backend).
cd "$HERE/flexfoil-ui"
echo "[run-highlift] open http://localhost:$PORT/highlift-preview.html"
npm run dev -- --port "$PORT"
