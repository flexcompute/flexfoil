#!/bin/bash
# Generate XFOIL buffer coordinates for comparison

cd /Users/harry/flexfoil-boundary-layer

XFOIL=./Xfoil/xfoil

for NACA in 0012 2412 4412; do
    echo "=== Generating NACA $NACA buffer ==="
    
    # Create XFOIL script
    cat > /tmp/xfoil_script_$NACA.txt << EOF
NACA $NACA
PSAV naca${NACA}_xfoil_buffer.dat
PANE
PSAV naca${NACA}_xfoil_paneled_new.dat

QUIT
EOF
    
    # Run XFOIL
    $XFOIL < /tmp/xfoil_script_$NACA.txt > /tmp/xfoil_output_$NACA.txt 2>&1
    
    echo "  Buffer saved to naca${NACA}_xfoil_buffer.dat"
    echo "  Paneled saved to naca${NACA}_xfoil_paneled_new.dat"
done

echo ""
echo "Done! Files generated:"
ls -la naca*_xfoil_*.dat
