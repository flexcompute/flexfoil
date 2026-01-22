# XFOIL Instrumented - Debug Version

This is a modified version of XFOIL that outputs JSON-formatted debug data for comparison with RustFoil.

## Building

```bash
cd bin
make
```

This will produce `xfoil_instrumented` binary.

## Usage

```bash
# Run with standard XFOIL commands
echo "NACA 0012
OPER
VISC 3000000
ALFA 4

QUIT" | ./xfoil_instrumented
```

## Output

The instrumented version produces `xfoil_debug.json` containing detailed snapshots of:

### Inviscid Solver
- **NCALC**: Normal vectors (NX, NY) at each node
- **APCALC**: Panel angles
- **TECALC**: Trailing edge geometry (ANTE, ASTE, DSTE, SHARP)
- **PSILIN**: Influence coefficients (PSIS, PSID, DZDG) for panel pairs
- **GGCALC**: Influence matrix AIJ and base solutions GAMU(1), GAMU(2)
- **SPECAL**: Combined gamma and QINV for given alpha
- **CLCALC**: Lift, moment coefficients and Cp distribution
- **STFIND**: Stagnation point location (IST, SST)

### Viscous Solver
- **VISCAL**: Iteration-level summary (CL, CD, residuals)
- **BLVAR**: Boundary layer secondary variables at each station
- **BLDIF**: Jacobian matrices (VS1, VS2) and residuals
- **MRCHUE**: Initial march state at each station
- **UPDATE**: Newton iteration deltas applied
- **QDCALC**: DIJ influence matrix
- **BLSOLV**: Newton solve system info
- **Closure functions**: HSL, CFL, DAMPL outputs

## JSON Format

```json
{
  "events": [
    {
      "call_id": 1,
      "subroutine": "VISCAL",
      "iteration": 1,
      "alpha_rad": 0.0698,
      "reynolds": 3000000.0
    },
    {
      "call_id": 2,
      "subroutine": "BLVAR",
      "side": 1,
      "ibl": 3,
      "input": {...},
      "output": {...}
    }
  ]
}
```

## Comparison with RustFoil

Use this output to compare with RustFoil's viscous solver at each stage:

1. Compare initial march states (MRCHUE vs RustFoil march_fixed_ue)
2. Compare Jacobian matrices (BLDIF vs RustFoil bldif)
3. Compare closure functions (HSL, CFL vs RustFoil closures)
4. Compare Newton iteration convergence

## Files Modified

- `xblsys.f`: BLVAR, BLDIF, DAMPL, HSL, CFL, HKIN instrumented
- `xbl.f`: MRCHUE, UPDATE instrumented
- `xoper.f`: VISCAL instrumented
- `xpanel.f`: QDCALC instrumented
- `xsolve.f`: BLSOLV instrumented

## New Files

- `xfoil_debug.f`: Debug output module with JSON helpers
