# XFOIL Implementation Rules

When implementing algorithms from XFOIL in this repository, follow these rules strictly.

## Rule 1: Read the XFOIL Source Code

**Always read and follow the actual Fortran source code exactly** rather than approximating or tuning parameters.

The XFOIL source is located at `/Xfoil/src/` with key files:
- `xfoil.f` - PANGEN (paneling), geometry routines
- `xpanel.f` - Panel method influence coefficients  
- `xoper.f` - Operating point calculations
- `spline.f` - Spline interpolation routines
- `xbl.f` - Boundary layer routines

## Rule 2: Match XFOIL Exactly

Our goal is to produce **identical numerical results** to XFOIL, not "close enough" approximations. If our results differ from XFOIL:
1. Read the XFOIL source code for the relevant routine
2. Compare line-by-line with our implementation
3. Fix any differences until results match exactly

## Rule 3: Use XFOIL's Spline Formulation

XFOIL uses **Hermite cubic splines** (`spline.f`):
- Solves for first derivatives (XS = dX/dS) at each knot
- Uses the SEVAL/DEVAL evaluation functions
- Different from standard "natural cubic splines" which solve for second derivatives

For exact matching, implement XFOIL's spline routines:
- `SPLINE` / `SPLIND` - Build spline coefficients
- `SEVAL` - Evaluate spline value
- `DEVAL` - Evaluate first derivative
- `D2VAL` - Evaluate second derivative
- `CURV` - Evaluate curvature

## Rule 4: Match XFOIL's Constants and Defaults

Key XFOIL defaults (from `xfoil.f`):
- `NPAN = 160` - Default number of panels
- `CVPAR = 1.0` - Curvature attraction parameter
- `CTERAT = 0.15` - TE/LE panel density ratio
- `CTRRAT = 0.2` - Refinement area/LE panel density ratio
- `RDSTE = 0.667` - TE panel spacing ratio
- `IPFAC = 5` - Panel node oversampling factor

## Rule 5: Test Against XFOIL Output

When implementing a new feature:
1. Run XFOIL to get reference output
2. Compare our output with XFOIL's
3. Document any remaining differences and their causes

Example test procedure:
```bash
# Generate XFOIL reference
echo "NACA 0012
PANE
PSAV naca0012_xfoil.dat
QUIT" | ./Xfoil/bin/xfoil

# Compare with our output
```

## Current Implementation Status

### NACA Generator - Exact Match
- `naca::naca4()` - Exact port of XFOIL's `SUBROUTINE NACA4`
- Uses XFOIL's TE bunching parameter (AN = 1.5)
- Uses XFOIL's x-distribution: `XX = 1 - ANP*FRAC*(1-FRAC)^AN - (1-FRAC)^ANP`
- Default NSIDE = 123 points per side (same as XFOIL's IQX/3 = 370/3)
- Produces perfectly symmetric output for symmetric airfoils
- **Exact match** with XFOIL buffer coordinates

### XFOIL Workflow (CRITICAL)
In XFOIL, after NACA generation, the workflow is:
1. `NACA4` generates raw buffer coordinates (XB/YB) - 245 points
2. `SCALC` + `SEGSPL` fit Hermite splines to buffer
3. **`PANGEN` is called automatically** to create working coordinates (X/Y)

The buffer coordinates are **never** used directly for analysis - PANGEN must always
be applied first to create the working coordinates.

**Frontend implementation**:
- When generating NACA, we automatically apply XFOIL paneling
- `coordinates` stores the buffer (for display/editing)
- `panels` stores the paneled result (for analysis)
- Default `nPanels = 160` (XFOIL's NPAN default)

### Paneling (PANGEN) - EXACT MATCH
Point-by-point comparison with XFOIL output:
- **RMS X error: 2.11e-7** (machine precision)
- **RMS Y error: 2.06e-8** (machine precision)
- **RMS total: 2.12e-7** (machine precision)
- **Max error: 3.75e-7** at any point
- **LE error: ~4e-10** - bit-identical at leading edge

**Output is perfectly symmetric** (verified by tests).

**Key implementation details for exact match**:
1. Use `XfoilSpline` for curvature calculation (XFOIL's CURV function)
2. Use `LEFIND` algorithm to find exact LE arc-length position
3. Compute `CVAVG` by sampling 7 points around LE using the spline (not buffer points)
4. Handle LE between buffer points in smoothing (XFOIL lines 1761-1783)
5. Use `Spline1D` with zero third derivative end conditions for curvature interpolation

### Splines - Implemented both versions
- `CubicSpline` - Natural cubic spline (stores original points for paneling)
- `XfoilSpline` - XFOIL's Hermite cubic spline (solves for 1st derivatives)
  - Includes `lefind()` method for exact LE location
- `Spline1D` - 1D spline for scalar values (curvature) with SEGSPL end conditions
- Both 2D splines give identical curvature values for the same input points

### Inviscid Solver (XPANEL) - EXACT MATCH
- Linear vorticity panel method (PSILIN formulation)
- Blunt TE handling with source/vortex panel (PSIG/PGAM)
- Kutta condition: γ₁ + γₙ = 0
- **Cl at α=0 for symmetric airfoil: 0.000000** (perfect symmetry)
- **Cl(-α) = -Cl(α)** (exact antisymmetry verified)
- **Cl_α ≈ 6.91/rad** (vs 2π, typical for panel methods)

**Key implementation details for exact match**:
1. Include BOTH TE nodes (upper and lower) for blunt TE airfoils
2. TE panel angle uses XFOIL formula: `ATAN2(-SX, SY) + PI`
3. Proper handling of ANTE/ASTE (normal/tangent TE gap projections)

### Boundary Layer - Not implemented
### Viscous-Inviscid Coupling - Not implemented

## XFOIL Source Code Reference

Key subroutines to understand:

| Routine | File | Purpose |
|---------|------|---------|
| PANGEN | xfoil.f | Panel node distribution |
| PANCOP | xfoil.f | Copy buffer to current airfoil |
| LEFIND | xfoil.f | Find leading edge |
| TECALC | xpanel.f | TE geometry calculations |
| SETUP | xpanel.f | Influence coefficients |
| CLCALC | xoper.f | Lift/moment integration |
| SPLINE | spline.f | Build Hermite spline |
| SEVAL | spline.f | Evaluate spline |
| CURV | spline.f | Evaluate curvature |
