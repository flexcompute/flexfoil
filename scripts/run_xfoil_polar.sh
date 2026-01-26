#!/bin/bash
# Run XFOIL polar for comparison

XFOIL=/Applications/ESP128/EngSketchPad/bin/xfoil
AIRFOIL=testdata/naca0012.dat

cd /Users/harry/flexfoil-boundary-layer

# Create XFOIL script
cat > /tmp/xfoil_script.txt << 'EOF'
LOAD testdata/naca0012.dat
OPER
VISC 1e6
ITER 100
PACC
xfoil_polar.txt

ASEQ -4 10 1

PACC

QUIT
EOF

# Run XFOIL
$XFOIL < /tmp/xfoil_script.txt

echo "=== XFOIL Polar Results ==="
cat xfoil_polar.txt
