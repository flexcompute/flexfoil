## Session Log: flexfoil-boundary-layer

### 14:30 - Project: flexfoil-boundary-layer

#### Viscous Boundary Layer Implementation - Complete

**Experiment:** Implement multiple turbulent BL models with documentation and tuned Newton-Raphson solver

**Result:** ✅ Success - All 44 tests passing

---

### What We Built

#### 1. Multiple Turbulent BL Models
Created `TurbulentModel` enum with three options:

| Model | Method | Best For |
|-------|--------|----------|
| **Head** (default) | θ + H₁ entrainment | Fast, attached flows |
| **XfoilCtau** | θ + Cτ shear lag | Separation, accuracy |
| **GreenLag** | θ + H₁ + Cτ | History effects |

**Files:**
- `crates/rustfoil-solver/src/boundary_layer/xfoil_turb.rs` - New XFOIL Cτ model
- `crates/rustfoil-solver/src/boundary_layer/mod.rs` - Model selection enum

#### 2. Comprehensive Documentation
Created `docs/VISCOUS_MODELS.md` with:
- Full mathematical formulations (von Kármán, Thwaites, Head, XFOIL Cτ)
- Closure relations (Ludwieg-Tillmann, Coles wall law)
- eᴺ transition prediction
- V-I coupling algorithm
- Model comparison table
- 7 academic references

#### 3. Tuned Newton-Raphson V-I Coupling
Updated `ViscousConfig` with:
```rust
relaxation_initial: 0.7,
relaxation_min: 0.3,
adaptive_relaxation: true,  // Auto-reduce on divergence
tolerance: 1e-4,
max_iterations: 100,
```

Features:
- Adaptive under-relaxation
- RMS + max residual monitoring
- Stall detection & recovery

#### 4. WASM Bindings
New exports for frontend:
```javascript
TurbModel.Head / XfoilCtau / GreenLag
analyze_viscous_extended(coords, alpha, Re, model, n_crit)
get_turbulent_model_info()  // JSON for "About" section
```

---

### Validation Tests

| Test | Status |
|------|--------|
| NACA 0012 Cl vs XFOIL | ✅ Within 0.5% |
| Lift curve slope (2π correction) | ✅ |
| Thwaites H(λ=0) vs Blasius | ✅ 2.61 vs 2.59 |
| Coles Cf formula | ✅ |
| XFOIL Cτ marching | ✅ |
| Newton-Raphson convergence | ✅ |

**Total:** 44 tests passing

---

### Files Changed
```
docs/VISCOUS_MODELS.md                          (new)
docs/session_log_2026-01-18.md                  (new)
crates/rustfoil-solver/src/lib.rs               (exports)
crates/rustfoil-solver/src/boundary_layer/
  ├── mod.rs                                    (TurbulentModel enum)
  ├── xfoil_turb.rs                             (new - XFOIL Cτ model)
  ├── closure.rs                                (Thwaites fix)
  └── transition.rs                             (Michel test fix)
crates/rustfoil-solver/src/viscous/coupling.rs  (Newton-Raphson tuning)
crates/rustfoil-solver/tests/validation.rs      (new tests)
crates/rustfoil-wasm/src/lib.rs                 (model selection API)
crates/rustfoil-wasm/Cargo.toml                 (serde_json dep)
```

---

### Key Learnings

1. **XFOIL uses Cτ (shear lag) not H₁ (entrainment)** - The lag equation `dCτ/ds = (Cτ_eq - Cτ)/L` captures history effects better than Head's method

2. **Coles Cf formula** - Uses `log₁₀(Reθ)` not `ln(Reθ)`: `GRT/2.3026 = log₁₀(Reθ)`

3. **Adaptive relaxation critical** - Fixed relaxation fails on separated flows; need to reduce ω when residual increases

4. **Panel methods give ~10% higher lift slope** - NACA 0012 at 12% thickness gives dCl/dα ≈ 6.9/rad vs theory 2π ≈ 6.28/rad

---

### Next Steps
- [ ] Add model selection UI controls in `SolvePanel.tsx`
- [ ] Test viscous polars against XFOIL reference data
- [ ] Add wake closure for Cd calculation
- [ ] Consider compressibility corrections (Karman-Tsien)

---

### Related Notes
- [[Boundary_Layer_Theory]]
- [[XFOIL_Implementation]]
- [[Panel_Methods]]
