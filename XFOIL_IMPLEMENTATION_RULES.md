# XFOIL Implementation Rules

Match XFOIL exactly. Read the Fortran source, not approximations.

## Source Files

XFOIL source at `/Xfoil/src/`:
- `xfoil.f` - PANGEN (paneling), geometry
- `xpanel.f` - Panel method (PSILIN)
- `xoper.f` - Operating point (CLCALC)
- `spline.f` - Spline interpolation
- `xbl.f` - Boundary layer

## Rules

1. **Read the source** - Don't guess. Port the Fortran directly.
2. **Match exactly** - Numerical results must be identical. If not, diff line-by-line and fix.
3. **Use XFOIL's splines** - Hermite cubic (1st derivative unknowns), not natural cubic (2nd derivative).
4. **Match defaults** - NPAN=160, CVPAR=1.0, CTERAT=0.15, RDSTE=0.667, IPFAC=5.
5. **Test against XFOIL output** - Run XFOIL, compare, document any differences.

```bash
# Generate XFOIL reference
echo "NACA 0012
PANE
PSAV naca0012_xfoil.dat
QUIT" | ./Xfoil/bin/xfoil
```

## Implementation Status

### NACA Generator ✓
`naca::naca4()` - exact port of XFOIL's NACA4. TE bunching AN=1.5, NSIDE=123.

### XFOIL Workflow
NACA → buffer (XB/YB) → spline fit → PANGEN → working coords (X/Y).

Buffer is never used directly. Frontend auto-applies PANGEN after NACA generation.

### Paneling (PANGEN) ✓
RMS error vs XFOIL: **2.1e-7** (machine precision). Perfectly symmetric.

Key details:
- `XfoilSpline` for curvature (CURV)
- `LEFIND` for exact LE arc-length
- `CVAVG` from 7 samples around LE
- `Spline1D` with zero-third-derivative BCs

### Inviscid Solver ✓
PSILIN formulation with blunt TE handling. **Cl at α=0: 0.000000** (perfect symmetry).

Key details:
- Both TE nodes (upper/lower) for blunt TE
- TE panel angle: `ATAN2(-SX, SY) + PI`
- ANTE/ASTE for TE gap

### Not Implemented
- Boundary layer
- Viscous-inviscid coupling

## Reference

| Routine | File | Purpose |
|---------|------|---------|
| PANGEN | xfoil.f | Panel distribution |
| LEFIND | xfoil.f | Find LE |
| PSILIN | xpanel.f | Influence coefficients |
| CLCALC | xoper.f | Cl/Cm integration |
| SPLINE/SEVAL/CURV | spline.f | Spline routines |
