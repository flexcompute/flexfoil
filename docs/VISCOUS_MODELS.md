# Viscous Boundary Layer Models in RustFoil

## Overview

RustFoil implements integral boundary layer methods for viscous flow analysis. This document describes the mathematical formulations, closure relations, and numerical methods used.

## Table of Contents

1. [Integral Boundary Layer Equations](#integral-boundary-layer-equations)
2. [Laminar Models](#laminar-models)
3. [Turbulent Models](#turbulent-models)
4. [Transition Prediction](#transition-prediction)
5. [Viscous-Inviscid Coupling](#viscous-inviscid-coupling)
6. [Model Comparison](#model-comparison)
7. [References](#references)

---

## Integral Boundary Layer Equations

All methods start from the **von Kármán momentum integral equation**:

```
dθ/ds + (H + 2 - M²)(θ/Uₑ)(dUₑ/ds) = Cf/2
```

where:
- **θ** = momentum thickness
- **δ*** = displacement thickness  
- **H = δ*/θ** = shape factor
- **Uₑ** = boundary layer edge velocity
- **Cf** = skin friction coefficient
- **M** = Mach number (≈0 for incompressible)
- **s** = arc length along surface

This equation alone is insufficient—we need **closure relations** to relate H, Cf, and other quantities. Different methods provide different second equations and closures.

---

## Laminar Models

### Thwaites' Method (1949)

A one-equation method that solves directly for θ using an empirical correlation.

**Key Equation:**
```
θ²(s) = (0.45ν/Uₑ⁶) ∫₀ˢ Uₑ⁵ ds'
```

**Closure Relations:**

The pressure gradient parameter:
```
λ = (θ²/ν)(dUₑ/ds)
```

Shape factor H(λ):
- λ = 0 (flat plate): H ≈ 2.59 (Blasius)
- λ > 0 (favorable): H decreases toward 2.0
- λ < 0 (adverse): H increases toward 3.5 (separation at λ ≈ -0.09)

```
H(λ) = 2.61 - 3.75λ + 5.24λ²    for λ ≥ 0
H(λ) = 2.61 + 10.3|λ|            for λ < 0
```

Skin friction via Thwaites' L(λ):
```
Cf = 2L(λ)/Reθ
L(λ) = 0.22 + 1.57λ - 1.8λ²     for λ ≥ 0
```

**Advantages:** Simple, fast, no differential equation to march
**Limitations:** Only valid for attached laminar flow

---

## Turbulent Models

RustFoil provides three turbulent boundary layer models:

### 1. Head's Entrainment Method (1958)

A two-equation method using momentum thickness θ and entrainment shape factor H₁.

**Governing Equations:**

1. Momentum equation:
```
dθ/ds + (H + 2)(θ/Uₑ)(dUₑ/ds) = Cf/2
```

2. Entrainment equation:
```
d(H₁θ)/ds = Cₑ(H₁)
```

where H₁ = (δ - δ*)/θ is the entrainment shape factor.

**Closure Relations:**

Shape factor relation (Ludwieg-Tillmann):
```
H₁ = 3.3 + 0.8234(H - 1.1)⁻¹·²⁸⁷   for H ≤ 1.6
H₁ = 3.3 + 1.5501(H - 0.6778)⁻³·⁰⁶⁴  for H > 1.6
```

Entrainment coefficient:
```
Cₑ = 0.0306(H₁ - 3.0)⁻⁰·⁶¹⁶⁹
```

Skin friction (Ludwieg-Tillmann):
```
Cf = 0.246 × 10⁻⁰·⁶⁷⁸ᴴ × Reθ⁻⁰·²⁶⁸
```

**Advantages:** Well-validated, stable marching
**Limitations:** Limited history effects, struggles with strong separation

---

### 2. Green's Lag-Entrainment Method (1977)

An extension of Head's method that adds a lag equation for history effects.

**Governing Equations:**

1. Momentum equation (same as Head)
2. Entrainment equation (same as Head)  
3. Lag equation:
```
d(Cτ)/ds = (Cτ_eq - Cτ)/L
```

where:
- Cτ = τw/(ρUₑ²) is the shear stress coefficient
- Cτ_eq is the equilibrium value from entrainment
- L is a lag length scale (typically ~10θ)

**Advantages:** Better history effects than Head
**Limitations:** More complex, additional tuning constants

---

### 3. XFOIL's Cτ Lag-Dissipation Method (Drela, 1987)

The method used in XFOIL/ISES codes. Uses momentum thickness and shear lag coefficient as primary variables.

**Governing Equations:**

1. Momentum equation:
```
dθ/ds + (H + 2 - M²)(θ/Uₑ)(dUₑ/ds) = Cf/2
```

2. Shear lag equation:
```
d(Cτ)/ds = (Cτ_eq - Cτ)/L
```

**Key Constants (from XFOIL):**
```
SCCON = 5.6    // Shear coefficient lag constant
GACON = 6.70   // G-β equilibrium constant A
GBCON = 0.75   // G-β equilibrium constant B
GCCON = 18.0   // G-β wall term constant
DLCON = 0.9    // Dissipation length ratio
```

**Equilibrium Cτ from G-β locus:**
```
G = GACON × √(1 + GBCON × β) + GCCON/(H × Reθ × √(Cf/2))
β = (θ/Uₑ)(dUₑ/ds) × 2/(Cf × Hs)
Cτ_eq = (Hs × Cf/2) × (G/6.7)²
```

**Skin Friction (Coles wall law):**
```
Cf = CFFAC × 0.3 × exp(-1.33H) × (ln(Reθ)/2.3)^(-1.74 - 0.31H)
```

**Advantages:** 
- Most accurate for separation prediction
- Handles strong adverse pressure gradients
- Validated extensively in XFOIL

**Limitations:**
- More complex implementation
- Requires careful tuning of constants

---

## Transition Prediction

### The eᴺ Method (Smith-Gamberoni, van Ingen)

Transition is predicted when the integrated amplification factor N reaches a critical value Ncrit.

**Amplification Rate:**
```
dN/ds = σ(Reθ, H)
```

where σ is the spatial growth rate of the most unstable Tollmien-Schlichting wave.

**Critical Reθ (Drela correlation):**
```
Reθ_crit = 155 + exp(6.54 - 0.75 × ln(Reθ_crit) + H)
```

This implicit equation gives Reθ_crit ≈ 200-400 depending on H.

**Amplification rate correlation:**
```
σ = f(Reθ/Reθ_crit - 1, H)
```

**Transition Criteria:**
- Default: Ncrit = 9 (clean wind tunnel)
- Low turbulence: Ncrit = 11-14
- High turbulence: Ncrit = 4-8

**Alternative: Michel's Criterion**
A simpler empirical criterion:
```
Reθ_tr = 1.174 × (1 + 22400/Rex) × Rex^0.46
```

---

## Viscous-Inviscid Coupling

### Transpiration Velocity Model

The boundary layer displaces the inviscid flow by an amount δ*. Rather than modifying the geometry, we apply a **transpiration velocity**:

```
Vn = d(Uₑδ*)/ds
```

This normal velocity is added to the inviscid boundary condition:
```
ψ(surface) = ψ₀ + ∫ Vn ds
```

### Global Newton-Raphson Iteration

The coupled system is solved iteratively:

**Algorithm:**
```
1. Initialize: Solve inviscid flow for Uₑ(s)
2. Loop until converged:
   a. March boundary layer with current Uₑ → get δ*(s)
   b. Compute transpiration: Vn = d(Uₑδ*)/ds  
   c. Re-solve inviscid with transpiration BC → get new Uₑ
   d. Compute residual: R = ||δ*_new - δ*_old||
   e. Apply under-relaxation: δ* = ω × δ*_new + (1-ω) × δ*_old
   f. If R < tolerance, converged
```

**Convergence Parameters:**
- Relaxation factor: ω = 0.5-0.9 (lower for separated flows)
- Tolerance: ||Δδ*/δ*|| < 10⁻⁴
- Maximum iterations: 50-100

**Jacobian Handling:**
XFOIL uses a full Newton method with analytical Jacobians. We use a quasi-Newton approach with numerical differentiation for simplicity.

---

## Model Comparison

### Accuracy vs. Complexity Trade-off

| Model | Equations | Separation | History | Speed |
|-------|-----------|------------|---------|-------|
| Thwaites (laminar) | 1 | Poor | None | ★★★★★ |
| Head | 2 | Fair | Limited | ★★★★ |
| Green Lag | 3 | Good | Good | ★★★ |
| XFOIL Cτ | 2 | Excellent | Excellent | ★★★ |

### When to Use Each Model

- **Thwaites**: Quick laminar-only estimates, educational purposes
- **Head**: General attached turbulent flows, fast polar generation
- **Green Lag**: Flows with mild separation, better accuracy than Head
- **XFOIL Cτ**: High accuracy requirements, strong pressure gradients, validation cases

### Validation Data

All models validated against:
1. Blasius flat plate (laminar): θ/x = 0.664/√Rex
2. Turbulent flat plate: Cf = 0.0592/Rex^(1/5)
3. NACA 0012 at Re = 10⁶ (comparison with XFOIL)
4. NACA 4412 with separation (comparison with experiment)

---

## Drag Calculation

### Friction Drag (XFOIL Method)
Wall shear stress integrated along both surfaces:
```
τ = 0.5 × ρ × Uₑ² × Cf

Cd_f = Σ 0.5×(τᵢ + τᵢ₋₁) × Δx × 2/U∞²
```
where Δx is projected onto the freestream direction:
```
Δx = (xᵢ - xᵢ₋₁)×cos(α) + (yᵢ - yᵢ₋₁)×sin(α)
```

### Total Drag (Squire-Young Formula)
From momentum deficit at wake end:
```
Cd = 2 × (θ_wake/c) × (Uₑ_wake/U∞)^((5+H_wake)/2)
```

The wake is marched downstream from the trailing edge with:
- Cf = 0 (no wall)
- H decreasing toward equilibrium (~1.0)
- θ approximately conserved (θ×Uₑ ≈ const)

### Pressure Drag
```
Cd_p = Cd_total - Cd_f
```

---

## References

1. **Thwaites, B.** (1949). "Approximate Calculation of the Laminar Boundary Layer." *Aeronautical Quarterly*, Vol. 1.

2. **Head, M.R.** (1958). "Entrainment in the Turbulent Boundary Layer." *ARC R&M 3152*.

3. **Green, J.E., Weeks, D.J., Brooman, J.W.F.** (1977). "Prediction of Turbulent Boundary Layers and Wakes in Compressible Flow by a Lag-Entrainment Method." *ARC R&M 3791*.

4. **Drela, M.** (1987). "Two-Dimensional Transonic Aerodynamic Design and Analysis Using the Euler Equations." *Ph.D. Thesis, MIT*.

5. **Drela, M.** (1989). "XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils." *Conference on Low Reynolds Number Airfoil Aerodynamics*, Notre Dame.

6. **White, F.M.** (2006). *Viscous Fluid Flow*, 3rd Ed. McGraw-Hill.

7. **Cebeci, T., Bradshaw, P.** (1977). *Momentum Transfer in Boundary Layers*. Hemisphere.

---

## Implementation Notes

### Numerical Stability

1. **θ floor**: Prevent θ < 10⁻¹⁰ to avoid division by zero
2. **H clamp**: Keep 1.0 < H < 10.0 for turbulent, 2.0 < H < 4.0 for laminar
3. **Cf floor**: Cf > 10⁻⁶ to prevent numerical issues
4. **Step size**: Δs < 0.1θ for accurate marching

### Stagnation Point Treatment

At s = 0:
- θ(0) = 0 causes singularities
- Use Thwaites similarity: θ ∝ s^0.5 near stagnation
- Start marching from s > 0 with initialized values

### Wake Treatment

In the wake (x > trailing edge):
- No wall: Cf = 0
- H increases (wake spreading)
- Use wake-specific closure relations
