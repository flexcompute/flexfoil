#!/bin/bash
# Build the standalone 2D anisotropic mesher (mesh2d) against the Flow360 compute
# install. The prebuilt binary at ../bin/mesh2d already works on this machine;
# rebuild only if the compute libraries change.
#
# Approach: borrow the include (-I/-isystem) flags from any CavityBasedMesher
# translation unit in compile_commands.json (dropping the per-target *_EXPORTS
# defines), then link the in-house mesh libraries the mesher uses.
set -euo pipefail
ROOT="${FLOW360_COMPUTE_ROOT:-/home/qiqi/flexcompute/compute}"
SRC="$ROOT/src/Flow360Core"
BUILD="$ROOT/build/release/Flow360Core"
RELEASE="$ROOT/install/release"
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${1:-$HERE/../bin/mesh2d}"

CC="$BUILD/compile_commands.json"
[ -f "$CC" ] || { echo "compile_commands.json not found at $CC (configure the build first)"; exit 1; }

# Pull include flags from a representative cavity-mesher compile command.
INCS=$(python3 - "$CC" <<'PY'
import json, sys, shlex
cmds = json.load(open(sys.argv[1]))
entry = next(c for c in cmds if "CavityBasedMesher" in c["file"])
flags = [t for t in shlex.split(entry["command"]) if t.startswith(("-I", "-isystem", "-std="))]
print(" ".join(flags))
PY
)

set -x
g++ -O2 -std=c++17 $INCS \
    "$HERE/mesh2d.cpp" -o "$OUT" \
    -L"$RELEASE/lib" \
    -lflow360cavitybasedmesher -lflow360meshreader -lflow360meshwriter \
    -lflow360meshstructs -lflow360meshprimitives -lflow360algorithms -lflow360logging \
    -lmpi -Wl,--copy-dt-needed-entries -Wl,-rpath,"$RELEASE/lib"
set +x
echo "built $OUT"
