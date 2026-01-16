# Mark Drela's Linear Vorticity Panel Method (XFOIL)

## Overview

This document describes the inviscid panel method used in XFOIL, based on direct analysis of the XFOIL 6.99 source code (`xpanel.f`, `xoper.f`, `xfoil.f`). The method uses **linear vorticity panels** with a **stream function formulation**.

## Key Concepts

### 1. Panel Representation

XFOIL uses **N nodes** defining **N-1 panels** on the airfoil surface, plus a special trailing edge (TE) panel connecting node N back to node 1.

- Nodes are numbered 1 to N, traversing the airfoil **counterclockwise** (starting from TE upper surface, around LE, back to TE lower surface)
- Each panel has a **linearly varying vorticity distribution** between its two endpoint nodes
- The vorticity at node j is denoted γⱼ (stored in `GAM(J)`)

### 2. Stream Function Formulation

The boundary condition is that the **stream function ψ is constant** on the airfoil surface:

```
ψ(node_i) = ψ₀   for all i = 1, ..., N
```

where ψ₀ is the (unknown) internal stream function value.

### 3. Linear Vorticity Influence Coefficients

For a panel from node JO to node JP with linearly varying vorticity, the stream function at field point (XI, YI) is computed using the **sum/difference formulation**:

```fortran
GSUM = GAM(JP) + GAM(JO)    ! Sum of endpoint vorticities
GDIF = GAM(JP) - GAM(JO)    ! Difference of endpoint vorticities

PSI = PSI + QOPI*(PSIS*GSUM + PSID*GDIF)
```

where `QOPI = 1/(4π)` and:

```fortran
PSIS = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2)
PSID = ((X1+X2)*PSIS + 0.5*(RS2*G2-RS1*G1 + X1*X1-X2*X2))*DXINV
```

The geometric quantities are:
- `(X1, YY)` = position of field point relative to node JO in panel-aligned coordinates
- `(X2, YY)` = position of field point relative to node JP in panel-aligned coordinates  
- `G1 = ln(r₁²)`, `G2 = ln(r₂²)` where r₁, r₂ are distances to panel endpoints
- `T1 = atan2(X1, YY)`, `T2 = atan2(X2, YY)` - angles to panel endpoints
- `DXINV = 1/(X1-X2)` = inverse panel length in local coordinates

### 4. Influence Coefficient Matrix Assembly

The influence matrix `AIJ(I,J)` stores dψ/dγⱼ at node i:

```fortran
C------ dPsi/dGam
DZDG(JO) = DZDG(JO) + QOPI*(PSIS-PSID)
DZDG(JP) = DZDG(JP) + QOPI*(PSIS+PSID)
```

This comes from the chain rule:
- Contribution to node JO: `QOPI*(PSIS - PSID)` (coefficient of γ_JO in GSUM - GDIF)
- Contribution to node JP: `QOPI*(PSIS + PSID)` (coefficient of γ_JP in GSUM + GDIF)

### 5. System Matrix Structure

The system has **N+1 unknowns**: γ₁, γ₂, ..., γₙ, ψ₀

The matrix equation is:

```
| A₁₁  A₁₂  ... A₁ₙ  -1 | | γ₁  |   | -ψ∞(node 1) |
| A₂₁  A₂₂  ... A₂ₙ  -1 | | γ₂  |   | -ψ∞(node 2) |
|  :    :   ...  :    : | |  :  | = |      :       |
| Aₙ₁  Aₙ₂  ... Aₙₙ  -1 | | γₙ  |   | -ψ∞(node N) |
|  1    0   ...  1    0 | | ψ₀  |   |      0       |
```

where:
- `Aᵢⱼ = dψ/dγⱼ` at node i (the influence coefficient)
- The `-1` column enforces ψ = ψ₀ on the surface
- The last row is the **Kutta condition**: γ₁ + γₙ = 0

### 6. Freestream Contribution

The freestream stream function at point (x, y) is:

```fortran
PSIINF = QINF*(COS(ALFA)*Y - SIN(ALFA)*X)
```

The RHS of the system is `-ψ∞(node_i)`:
```fortran
RES1 =  QINF*Y(I)    ! For alpha = 0°
RES2 = -QINF*X(I)    ! For alpha = 90°
```

### 7. Two-Solution Superposition

XFOIL solves the system **twice** - once for α=0° and once for α=90° - storing the results in `GAMU(I,1)` and `GAMU(I,2)`. The solution for any angle of attack is then:

```fortran
GAM(I) = COS(ALFA)*GAMU(I,1) + SIN(ALFA)*GAMU(I,2)
```

This is the key efficiency of XFOIL - the expensive matrix factorization is done once, then any α can be computed instantly.

### 8. Surface Velocity

The surface velocity equals the local vorticity strength:

```
V_surface(i) = γᵢ
```

This is because:
- Inside the body: V = 0
- Outside: V = γ (jump across vortex sheet)
- On the surface: V = γ (taking the outside limit)

### 9. Trailing Edge Treatment

For a **sharp trailing edge** (`SHARP = .TRUE.`):
- The standard Kutta condition γ₁ + γₙ = 0 is used
- Additionally, the velocity along the TE bisector is set to zero at a point just inside the TE

For a **blunt trailing edge**:
- A TE panel is added with both source and vortex components
- `SIGTE = 0.5*SCS*(GAM(N) - GAM(1))` - source strength
- `GAMTE = -0.5*SDS*(GAM(N) - GAM(1))` - vortex strength
- where SCS, SDS are geometric factors based on TE geometry

### 10. Lift and Moment Calculation

Lift is computed by integrating the pressure coefficient around the airfoil:

```fortran
CL = ∫ Cp * dx  (integrated around airfoil)
```

where the pressure coefficient uses the Karman-Tsien compressibility correction:

```fortran
CGINC = 1.0 - (GAM(I)/QINF)**2    ! Incompressible Cp
CPG = CGINC/(BETA + BFAC*CGINC)   ! Compressible Cp
```

with `BETA = sqrt(1 - M∞²)`.

## Key Constants

```fortran
PI   = 3.14159265...
HOPI = 1/(2π) = 0.159154943...
QOPI = 1/(4π) = 0.079577471...
```

## Implementation Notes

### Panel Local Coordinate System

For each panel from node JO to JP:
1. Compute panel length: `DSO = |P_JP - P_JO|`
2. Compute unit tangent: `(SX, SY) = (X_JP - X_JO, Y_JP - Y_JO) / DSO`
3. Transform field point to local coords:
   - `X1 = SX*(XI-X_JO) + SY*(YI-Y_JO)` (tangential distance to JO)
   - `X2 = SX*(XI-X_JP) + SY*(YI-Y_JP)` (tangential distance to JP)
   - `YY = SX*(YI-Y_JO) - SY*(XI-X_JO)` (normal distance, same for both)

### Self-Influence

When the field point is at a panel endpoint (IO = JO or IO = JP):
- The logarithm terms G1 or G2 are set to 0
- The arctangent terms T1 or T2 are set to 0
- This avoids the singularity and gives the correct principal value

### Sign Conventions

- Positive γ induces **counterclockwise** circulation
- For a lifting airfoil at positive α, γ > 0 on the upper surface, γ < 0 on the lower surface
- The Kutta condition γ₁ + γₙ = 0 ensures smooth flow departure at the TE

## Comparison with Constant-Strength Methods

| Aspect | Constant Vortex | Linear Vortex (XFOIL) |
|--------|-----------------|----------------------|
| Unknowns | N (panel strengths) | N (node strengths) |
| Continuity | Discontinuous γ | Continuous γ |
| Accuracy | O(1/N) | O(1/N²) |
| LE resolution | Poor | Excellent |
| Influence calc | Simpler | More complex |

## References

1. Drela, M. "XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils", Conference on Low Reynolds Number Airfoil Aerodynamics, 1989.
2. XFOIL 6.99 Source Code, MIT (2013)
3. Katz, J. and Plotkin, A. "Low-Speed Aerodynamics", Cambridge University Press, 2001.
