# RustFoil vs XFOIL Systematic Comparison Results

## Changes Made

### 1. Fixed CD_f Integration (Ue² factor)
**Location**: `crates/rustfoil-solver/src/viscous/forces.rs`

XFOIL computes friction drag as:
```
CD_f = ∫ TAU * dx * 2/QINF²
```
where `TAU = 0.5 * ρ * Ue² * Cf`.

Our original code was missing the Ue² factor:
```rust
// WRONG: cf_avg * dx
// CORRECT: ue_sq_avg * cf_avg * dx
```

### 2. Fixed Far Wake Extrapolation
**Location**: `crates/rustfoil-solver/src/viscous/forces.rs`

XFOIL computes total CD using Squire-Young at the far wake (x≈2c), not at TE.
Analysis of XFOIL wake data shows θ follows `θ*Ue² ≈ const` (not `θ*Ue = const`).

Updated extrapolation:
```rust
// Old: θ_far = θ_TE * Ue_TE / Ue_far  (1% error)
// New: θ_far = θ_TE * (Ue_TE / Ue_far)²  (matches XFOIL wake evolution)
```

## Results

### Before Fixes
| Alpha | CL Error | CD Error |
|-------|----------|----------|
| 0° | N/A | -25.2% |
| 2° | -21.6% | -30.2% |
| 4° | -5.4% | -36.4% |
| 6° | +0.2% | -54.9% |
| 8° | +10.3% | +8.7% |

### After Fixes
| Alpha | CL Error | CD Error | Notes |
|-------|----------|----------|-------|
| 0° | N/A | **-0.7%** | Near perfect |
| 2° | -21.6% | **-4.5%** | Good |
| 4° | -5.4% | **-4.9%** | Good |
| 6° | +0.2% | -50.4% | High H at TE (separation) |
| 8° | +10.3% | +32.2% | Different flow regime |

**Key improvement**: CD error at moderate angles (0-4°) reduced from 25-36% to within 5%.

## Remaining Issues

### 1. High H at TE (α=6-8°)
At α=6°, upper surface TE has H=7.96 (should be ~2.0 for attached turbulent flow).
This indicates near-separation or improper turbulent closure.

**Debug output at α=6°:**
```
Upper TE: θ=2.76e-3, H=7.96, Ue=0.77
Lower TE: θ=1.30e-3, H=3.62, Ue=0.78
```

The H spike prevents proper wake extrapolation and causes the CD to be underestimated.

### 2. CL Errors at Low Alpha (0-2°)
- α=0°: Non-zero CL for symmetric airfoil (stagnation discretization)
- α=2°: 22% CL error (Newton not fully converging)

These are related to the Newton iteration not reaching tight tolerance (residual plateaus at ~1.6).

## Key Intermediate Quantities Comparison

### θ at TE (α=4°)
| Surface | RustFoil | XFOIL | Error |
|---------|----------|-------|-------|
| Upper | 3.99e-3 | 4.76e-3 | -16% |
| Lower | 1.43e-3 | 1.42e-3 | +0.7% |
| Combined | 5.42e-3 | 6.18e-3 | -12% |

### Ue at Key Locations (α=4°)
| Location | RustFoil | XFOIL | Error |
|----------|----------|-------|-------|
| Upper x=0.5 | 1.168 | 1.187 | -1.6% |
| Upper TE | 0.805 | 0.765 | +5.2% |
| Lower x=0.5 | 1.042 | 1.031 | +1.1% |
| Lower TE | 0.802 | 0.765 | +4.7% |

Ue values match within 5% - good.

### H at Key Locations (α=4°)
| Location | RustFoil | XFOIL | Error |
|----------|----------|-------|-------|
| Upper x=0.3 | 3.15 | 1.69 | +86% ← Issue |
| Upper x=0.5 | 1.45 | 1.47 | -1% |
| Upper TE | 1.78 | 1.91 | -7% |
| Lower x=0.9 | 4.58 | 5.22 | -12% |
| Lower TE | 3.77 | 1.47 | +157% ← Issue |

The H values show significant discrepancies near transition, indicating
issues with the turbulent closure immediately after transition.

## Recommendations

1. **Investigate turbulent H evolution after transition** - H should drop from ~2.6 (laminar) to ~1.5 (turbulent) quickly after transition, but our code shows H staying high.

2. **Add wake stations** - Computing a proper wake (x=1 to x=2) instead of extrapolation would improve accuracy at all angles.

3. **Fix Newton convergence at low alpha** - The residual plateaus suggest an issue with the system conditioning or update limiting.
