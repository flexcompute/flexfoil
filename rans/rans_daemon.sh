#!/bin/bash
# RANS runner daemon — run ONCE in a normal (non-sandboxed) interactive shell:
#
#   nohup bash /home/qiqi/flexcompute/flexfoil/rans/rans_daemon.sh \
#        > /tmp/rans_daemon.log 2>&1 &
#
# It watches a queue directory and executes the full pipeline (mesh → case →
# GPU solve) for each submitted job, so the agent can drive Flow360 without you
# launching each run. The GPU solver only works in this shell context, not in the
# agent's sandbox — hence this split.
#
# Job protocol (the agent writes these):
#   $QUEUE/<job>/case.json     the contours JSON (rans CaseConfig)
#   $QUEUE/<job>/SUBMITTED     marker file → daemon picks the job up
# Daemon writes back:
#   $QUEUE/<job>/out/          run_case.py output (summary.json, *.pvtu, forces…)
#   $QUEUE/<job>/run.log       combined solver/pipeline log
#   $QUEUE/<job>/DONE          marker with the exit code (agent polls for this)
set -u
QUEUE="${RANS_QUEUE:-/tmp/rans_queue}"
RANS_DIR="$(cd "$(dirname "$0")" && pwd)"
PY="${FLOW360_VENV_PY:-/home/qiqi/flexcompute/compute/.venv/bin/python}"
mkdir -p "$QUEUE"
echo "[rans-daemon] watching $QUEUE  (pipeline: $RANS_DIR/run_case.py, py: $PY)"

while true; do
  for sub in "$QUEUE"/*/SUBMITTED; do
    [ -e "$sub" ] || continue
    job="$(dirname "$sub")"
    [ -e "$job/DONE" ] && continue
    [ -e "$job/case.json" ] || { echo "$(date -Is) [skip] $job: no case.json"; continue; }
    echo "$(date -Is) [run ] $job"
    "$PY" "$RANS_DIR/run_case.py" "$job/case.json" -o "$job/out" --solve > "$job/run.log" 2>&1
    rc=$?
    echo "$rc" > "$job/DONE"
    echo "$(date -Is) [done] $job  rc=$rc"
  done
  sleep 2
done
