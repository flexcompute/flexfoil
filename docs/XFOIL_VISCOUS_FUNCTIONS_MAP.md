# XFOIL Viscous Flow Functions: Complete Reference

**Date:** 2026-01-20  
**Purpose:** Complete mapping of all XFOIL functions related to viscous flow solving for FlexFoil implementation.

**Navigation:**
- [Main Viscous Solver](#main-viscous-solver)
- [Boundary Layer System](#boundary-layer-system)
- [Panel Method Coupling](#panel-method-coupling)
- [Newton Solver](#newton-solver)
- [Closure Relations](#closure-relations)
- [Transition Prediction](#transition-prediction)
- [Utility Functions](#utility-functions)

---

## Main Viscous Solver

### VISCAL
**File:** `xoper.f:2886`  
**Purpose:** Main viscous operating point convergence routine. Orchestrates the entire viscous-inviscid coupling iteration.

**Algorithm:**
```
1. Setup wake geometry (XYWAKE)
2. Compute wake velocities (QWCALC)
3. Set initial velocities (QISET)
4. Locate stagnation point (STFIND)
5. Map BL stations to panels (IBLPAN)
6. Compute arc lengths (XICALC)
7. Map BL stations to system indices (IBLSYS)
8. Compute inviscid edge velocities (UICALC)
9. Newton iteration loop:
   a. Build BL system (SETBL)
   b. Solve Newton system (BLSOLV)
   c. Update BL variables (UPDATE)
   d. Update Mach/Re or alpha
   e. Compute edge velocities (QVFUE)
   f. Update circulation (GAMQV)
   g. Move stagnation point (STMOVE)
   h. Compute forces (CLCALC, CDCALC)
   i. Check convergence
```

**Rust Implementation:**
```rust
pub fn viscal(
    &mut self,
    niter: usize,
    alpha: f64,
    reynolds: f64,
    mach: f64,
) -> Result<ViscousSolution, Error> {
    // Main iteration loop
    for iter in 0..niter {
        self.setbl()?;
        self.blsolv()?;
        self.update()?;
        // ... rest of iteration
    }
}
```

**Dependencies:**
- SETBL, BLSOLV, UPDATE
- XYWAKE, QWCALC, QISET
- STFIND, IBLPAN, XICALC, IBLSYS
- UICALC, QVFUE, GAMQV, STMOVE
- CLCALC, CDCALC

---

### SETBL
**File:** `xbl.f:21`  
**Purpose:** Sets up the BL Newton system coefficients. Builds the global Jacobian matrix for the coupled viscous-inviscid system.

**Key Operations:**
1. Sets compressibility parameters
2. Computes Reynolds number
3. Initializes BL if needed (MRCHUE)
4. Marches BL to establish transition (MRCHDU)
5. Builds Newton system coefficients for each station
6. Computes coupling through DIJ matrix

**Variables Set:**
- `VA(3,3,IVX)`: Diagonal blocks
- `VB(3,3,IVX)`: Subdiagonal blocks  
- `VM(3,IVX,IVX)`: Coupling blocks (through DIJ)
- `VDEL(3,2,IVX)`: Residuals

**Rust Implementation:**
```rust
pub fn setbl(&mut self) -> Result<(), Error> {
    // Set compressibility
    self.comset()?;
    
    // Initialize BL if needed
    if !self.lblini {
        self.mrchue()?;
        self.lblini = true;
    }
    
    // March to establish transition
    self.mrchdu()?;
    
    // Build Newton system
    for is in 0..2 {  // upper/lower surfaces
        for ibl in 2..self.nbl[is] {
            self.build_station_system(ibl, is)?;
        }
    }
}
```

**Dependencies:**
- MRCHUE, MRCHDU
- BLSYS, TESYS
- UESET (for Ue from mass defect)

---

### UPDATE
**File:** `xbl.f:1253`  
**Purpose:** Updates BL variables after Newton solve. Applies under-relaxation and enforces limits.

**Key Operations:**
1. Applies Newton updates with relaxation
2. Updates CTAU, THET, DSTR, MASS
3. Updates UEDG from mass defect (UESET)
4. Updates DSTR from mass defect (DSSET)
5. Enforces limits (DSLIM)

**Rust Implementation:**
```rust
pub fn update(&mut self) -> Result<(), Error> {
    // Apply updates with relaxation
    let rlx = self.compute_relaxation()?;
    
    for is in 0..2 {
        for ibl in 2..self.nbl[is] {
            self.ctau[ibl][is] += rlx * self.vdel[1][ibl][is];
            self.thet[ibl][is] += rlx * self.vdel[2][ibl][is];
            self.mass[ibl][is] += rlx * self.vdel[3][ibl][is];
        }
    }
    
    // Update Ue from mass defect
    self.ueset()?;
    
    // Update δ* from mass defect
    self.dsset()?;
    
    // Enforce limits
    self.dslim()?;
}
```

**Dependencies:**
- UESET, DSSET
- DSLIM

---

## Boundary Layer System

### MRCHUE
**File:** `xbl.f:542`  
**Purpose:** Marches boundary layers in direct mode using UEDG array. Switches to inverse mode when separation detected.

**Key Algorithm:**
```
For each station:
  1. Try direct mode (3x3 system: [Ctau, θ, δ*])
  2. Check if Hk > threshold:
     - Laminar: HLMAX = 3.8
     - Turbulent: HTMAX = 2.5
  3. If Hk > threshold:
     - Compute target Hk (HTARG)
     - Switch to inverse mode (4x4 system: [Ctau, θ, δ*, Ue])
     - 4th equation: Hk - HTARG = 0
  4. Solve Newton system for station
  5. Update variables
```

**Target Hk Evolution:**
- Laminar: `HTARG = HK1 + 0.03*(X2-X1)/T1`
- Turbulent: `HTARG = HK1 - 0.15*(X2-X1)/T1`
- Wake: Newton solve for `H + const*(H-1)^3 = H_prev`

**Rust Implementation:**
```rust
pub fn mrchue(&mut self) -> Result<(), Error> {
    const HLMAX: f64 = 3.8;  // Laminar threshold
    const HTMAX: f64 = 2.5;   // Turbulent threshold
    
    for is in 0..2 {
        // Initialize similarity station
        self.init_similarity_station(is)?;
        
        for ibl in 2..self.nbl[is] {
            let mut direct = true;
            
            // Try direct mode
            let (ctau, theta, dstar) = self.solve_direct(ibl, is)?;
            let hk_test = self.compute_hk(dstar, theta, self.uedg[ibl][is])?;
            
            // Check threshold
            let hmax = if ibl < self.itran[is] { HLMAX } else { HTMAX };
            if hk_test >= hmax {
                direct = false;
            }
            
            if direct {
                // Direct mode: solve 3x3 system
                self.solve_station_direct(ibl, is)?;
            } else {
                // Inverse mode: solve 4x4 system
                let hk_target = self.compute_target_hk(ibl, is)?;
                self.solve_station_inverse(ibl, is, hk_target)?;
            }
        }
    }
}
```

**Dependencies:**
- BLSYS, TESYS
- BLPRV, BLKIN, BLVAR, BLMID
- HKIN (for Hk computation)

---

### MRCHDU
**File:** `xbl.f:875`  
**Purpose:** Marches BL with current Ue and Ds to establish transition location.

**Key Algorithm:**
```
For each station:
  1. March BL with current Ue
  2. Check transition (TRCHEK)
  3. Update transition location if needed
```

**Rust Implementation:**
```rust
pub fn mrchdu(&mut self) -> Result<(), Error> {
    for is in 0..2 {
        for ibl in 2..self.nbl[is] {
            // March BL
            self.march_station(ibl, is)?;
            
            // Check transition
            if !self.turb[ibl][is] {
                self.trchek(ibl, is)?;
                if self.trans {
                    self.itran[is] = ibl;
                }
            }
        }
    }
}
```

**Dependencies:**
- TRCHEK
- BLSYS

---

### BLSYS
**File:** `xblsys.f:583`  
**Purpose:** Builds the 3x3 or 4x4 Newton system for a single BL station.

**System Structure:**
```
Direct mode (3x3):
[∂R1/∂Ctau  ∂R1/∂θ   ∂R1/∂δ*] [dCtau]   [R1]
[∂R2/∂Ctau  ∂R2/∂θ   ∂R2/∂δ*] [dθ   ] = [R2]
[∂R3/∂Ctau  ∂R3/∂θ   ∂R3/∂δ*] [dδ*  ]   [R3]

Inverse mode (4x4):
[∂R1/∂Ctau  ∂R1/∂θ   ∂R1/∂δ*   ∂R1/∂Ue] [dCtau]   [R1]
[∂R2/∂Ctau  ∂R2/∂θ   ∂R2/∂δ*   ∂R2/∂Ue] [dθ   ] = [R2]
[∂R3/∂Ctau  ∂R3/∂θ   ∂R3/∂δ*   ∂R3/∂Ue] [dδ*  ]   [R3]
[0          ∂Hk/∂θ   ∂Hk/∂δ*   ∂Hk/∂Ue] [dUe  ]   [Hk-HTARG]
```

**Residuals:**
- R1: Momentum equation
- R2: Shape equation (or Hk - HTARG in inverse mode)
- R3: Lag equation (or Ue constraint in direct mode)

**Rust Implementation:**
```rust
pub fn blsys(
    &mut self,
    station: usize,
    surface: usize,
    inverse: bool,
) -> Result<(), Error> {
    // Compute current state
    self.blprv(station, surface)?;
    self.blkin()?;
    
    // Compute residuals
    let r_mom = self.momentum_residual(station, surface)?;
    let r_shape = if inverse {
        self.hk_target - self.hk_current
    } else {
        self.shape_residual(station, surface)?
    };
    let r_lag = self.lag_residual(station, surface)?;
    
    // Compute Jacobian
    if inverse {
        // 4x4 system
        self.build_jacobian_4x4(station, surface)?;
    } else {
        // 3x3 system
        self.build_jacobian_3x3(station, surface)?;
    }
}
```

**Dependencies:**
- BLPRV, BLKIN
- BLVAR, BLMID
- BLDIF, TRDIF
- All closure relations (CFL, CFT, HSL, HST, etc.)

---

### TESYS
**File:** `xblsys.f:664`  
**Purpose:** Builds system for trailing edge station (combines upper and lower surface).

**Key Algorithm:**
```
TE station combines:
- θ_te = θ_upper + θ_lower
- δ*_te = δ*_upper + δ*_lower + gap
- Ctau_te = weighted average
```

**Rust Implementation:**
```rust
pub fn tesys(
    &mut self,
    cte: f64,
    tte: f64,
    dte: f64,
) -> Result<(), Error> {
    // Combine upper and lower surface values
    let theta_upper = self.thet[self.iblte[0]][0];
    let theta_lower = self.thet[self.iblte[1]][1];
    let tte = theta_upper + theta_lower;
    
    let dstar_upper = self.dstr[self.iblte[0]][0];
    let dstar_lower = self.dstr[self.iblte[1]][1];
    let dte = dstar_upper + dstar_lower + self.wgap[0];
    
    // Build TE system
    self.build_te_system(cte, tte, dte)?;
}
```

---

### BLPRV
**File:** `xblsys.f:701`  
**Purpose:** Sets primary variables for BL station computation.

**Inputs:**
- XSI: Arc length
- AMI/CTI: Amplification or Ctau
- THI: Momentum thickness
- DSI: Displacement thickness
- DSWAKI: Wake gap (if wake)
- UEI: Edge velocity

**Rust Implementation:**
```rust
pub fn blprv(
    &mut self,
    xsi: f64,
    ami: f64,
    cti: f64,
    thi: f64,
    dsi: f64,
    dswaki: f64,
    uei: f64,
) {
    self.x1 = xsi;
    self.ampl1 = ami;
    self.ctau1 = cti;
    self.t1 = thi;
    self.d1 = dsi;
    self.dswak1 = dswaki;
    self.u1 = uei;
}
```

---

### BLKIN
**File:** `xblsys.f:725`  
**Purpose:** Computes kinematic quantities from primary variables.

**Computes:**
- H = δ*/θ (shape factor)
- Hk (kinematic shape factor)
- Reθ (Reynolds number based on θ)
- MSQ (Mach number squared)

**Rust Implementation:**
```rust
pub fn blkin(&mut self) -> Result<(), Error> {
    // Shape factor
    self.h1 = self.d1 / self.t1.max(1e-10);
    
    // Kinematic shape factor
    let msq = self.u1.powi(2) * self.hstinv 
              / (self.gm1bl * (1.0 - 0.5 * self.u1.powi(2) * self.hstinv));
    self.hk1 = self.hkin(self.h1, msq)?;
    
    // Reynolds number
    self.rt1 = self.reybl * self.u1 * self.t1;
    
    // Mach squared
    self.msq1 = msq;
}
```

**Dependencies:**
- HKIN

---

### BLVAR
**File:** `xblsys.f:784`  
**Purpose:** Computes all BL variables from primary variables.

**ITYP:**
- 1: Laminar
- 2: Turbulent
- 3: Wake

**Computes:**
- H, Hk, Hs, Hc (shape factors)
- Cf (skin friction)
- Di (dissipation)
- Us (friction velocity)
- CQ (entrainment)
- DE (energy thickness)

**Rust Implementation:**
```rust
pub fn blvar(&mut self, ityp: i32) -> Result<(), Error> {
    match ityp {
        1 => {
            // Laminar
            self.hs1 = self.hsl(self.hk1, self.rt1, self.msq1)?;
            self.cf1 = self.cfl(self.hk1, self.rt1, self.msq1)?;
            self.di1 = self.dil(self.hk1, self.rt1)?;
        }
        2 => {
            // Turbulent
            self.hs1 = self.hst(self.hk1, self.rt1, self.msq1)?;
            self.cf1 = self.cft(self.hk1, self.rt1, self.msq1)?;
            self.di1 = self.dit(self.hs1, self.us1, self.cf1, self.st1)?;
        }
        3 => {
            // Wake
            self.hs1 = self.hst(self.hk1, self.rt1, self.msq1)?;
            self.cf1 = 0.0;  // No wall
            self.di1 = self.dilw(self.hk1, self.rt1)?;
        }
        _ => return Err(Error::InvalidType),
    }
    
    self.hc1 = self.hct(self.hk1, self.msq1)?;
    self.us1 = self.u1 * (self.cf1 / 2.0).sqrt();
}
```

**Dependencies:**
- HSL, HST, HCT
- CFL, CFT
- DIL, DILW, DIT

---

### BLMID
**File:** `xblsys.f:1124`  
**Purpose:** Computes mid-point values between stations 1 and 2.

**Rust Implementation:**
```rust
pub fn blmid(&mut self, ityp: i32) -> Result<(), Error> {
    // Average values
    self.hm = 0.5 * (self.h1 + self.h2);
    self.hkm = 0.5 * (self.hk1 + self.hk2);
    self.rtm = 0.5 * (self.rt1 + self.rt2);
    self.msqm = 0.5 * (self.msq1 + self.msq2);
    
    // Compute mid-point variables
    self.blvar_mid(ityp)?;
}
```

---

### BLDIF
**File:** `xblsys.f:1552`  
**Purpose:** Computes differences (residuals) for BL equations.

**Residuals:**
- Momentum: `dθ/ds + (H+2)(θ/Ue)(dUe/ds) - Cf/2`
- Shape: Closure relation
- Lag: `dCtau/ds - (Ctau_eq - Ctau)/L`

**Rust Implementation:**
```rust
pub fn bldif(&mut self, ityp: i32) -> Result<(), Error> {
    let ds = self.x2 - self.x1;
    
    // Momentum residual
    let dtheta_ds = (self.t2 - self.t1) / ds;
    let due_ds = (self.u2 - self.u1) / ds;
    let r_mom = dtheta_ds + (self.hm + 2.0) * (self.tm / self.um) * (due_ds / self.um) 
                - self.cfm / 2.0;
    
    // Shape residual
    let r_shape = self.shape_residual(ityp)?;
    
    // Lag residual
    let dctau_ds = (self.ctau2 - self.ctau1) / ds;
    let ctau_eq = self.ctau_equilibrium(ityp)?;
    let lag_length = self.lag_length(ityp)?;
    let r_lag = dctau_ds - (ctau_eq - self.ctaum) / lag_length;
    
    self.vsrez[1] = r_mom;
    self.vsrez[2] = r_shape;
    self.vsrez[3] = r_lag;
}
```

---

### TRDIF
**File:** `xblsys.f:1195`  
**Purpose:** Computes transition interval differences.

**Rust Implementation:**
```rust
pub fn trdif(&mut self) -> Result<(), Error> {
    // Transition occurs between stations 1 and 2
    // Compute weighted average of laminar and turbulent values
    let xt = self.xt;  // Transition location
    let x1 = self.x1;
    let x2 = self.x2;
    
    let frac = (xt - x1) / (x2 - x1);
    
    // Interpolate variables
    self.interpolate_transition(frac)?;
}
```

---

## Panel Method Coupling

### QDCALC
**File:** `xpanel.f:1149`  
**Purpose:** Computes DIJ matrix: `DIJ(i,j) = dQtan(i)/dSig(j)`

**Key Algorithm:**
```
1. Start with BIJ (source influence on velocity)
2. Account for wake coupling:
   DIJ(i,j) = BIJ(i,j) + Σ_wake CIJ(iw,k)*DIJ(k,j)
3. Store full N×N matrix
```

**Physical Meaning:**
- If we add a source at panel j, how does it affect velocity at panel i?
- Sources model displacement: `Sig ~ d(Ue·δ*)/ds`

**Rust Implementation:**
```rust
pub fn qdcalc(&mut self) -> Result<(), Error> {
    let n = self.n_panels;
    let nw = self.n_wake;
    
    // Initialize with BIJ
    for i in 0..n {
        for j in 0..n {
            self.dij[(i, j)] = self.bij[(i, j)];
        }
    }
    
    // Account for wake coupling
    for iw in 0..nw {
        for k in 0..n {
            for j in 0..n {
                self.dij[(i, j)] += self.cij[(iw, k)] * self.dij[(k, j)];
            }
        }
    }
    
    self.ladij = true;
    self.lwdij = true;
}
```

**Dependencies:**
- PSILIN (for source influence)
- GGCALC (for geometry)

---

### UESET
**File:** `xpanel.f:1758`  
**Purpose:** Updates edge velocity from mass defect using DIJ matrix.

**Key Equation:**
```
Ue[i] = Ue_inviscid[i] + Σ_j DIJ(i,j) * mass[j]
```

**Rust Implementation:**
```rust
pub fn ueset(&mut self) -> Result<(), Error> {
    for is in 0..2 {
        for ibl in 2..self.nbl[is] {
            let i = self.ipan[ibl][is];
            let mut dui = 0.0;
            
            for js in 0..2 {
                for jbl in 2..self.nbl[js] {
                    let j = self.ipan[jbl][js];
                    let vti_ibl = self.vti[ibl][is];
                    let vti_jbl = self.vti[jbl][js];
                    let dij_ij = self.dij[(i, j)];
                    let mass_j = self.mass[jbl][js];
                    
                    let ue_m = -vti_ibl * vti_jbl * dij_ij;
                    dui += ue_m * mass_j;
                }
            }
            
            self.uedg[ibl][is] = self.uinv[ibl][is] + dui;
        }
    }
}
```

**Dependencies:**
- QDCALC (DIJ matrix must be computed)

---

### DSSET
**File:** `xpanel.f:1786`  
**Purpose:** Updates displacement thickness from mass defect.

**Key Equation:**
```
δ* = mass / Ue
```

**Rust Implementation:**
```rust
pub fn dsset(&mut self) -> Result<(), Error> {
    for is in 0..2 {
        for ibl in 2..self.nbl[is] {
            self.dstr[ibl][is] = self.mass[ibl][is] / self.uedg[ibl][is].max(1e-10);
        }
    }
}
```

---

### UICALC
**File:** `xpanel.f:1542`  
**Purpose:** Computes inviscid edge velocity at BL stations.

**Rust Implementation:**
```rust
pub fn uicalc(&mut self) -> Result<(), Error> {
    for is in 0..2 {
        for ibl in 2..self.nbl[is] {
            let i = self.ipan[ibl][is];
            self.uinv[ibl][is] = self.qinv[i].abs();
        }
    }
}
```

---

### QVFUE
**File:** `xpanel.f:1580`  
**Purpose:** Computes viscous edge velocity QVIS from UEDG.

**Rust Implementation:**
```rust
pub fn qvfue(&mut self) -> Result<(), Error> {
    for is in 0..2 {
        for ibl in 2..self.nbl[is] {
            let i = self.ipan[ibl][is];
            let sign = if self.sgnue[i] > 0.0 { 1.0 } else { -1.0 };
            self.qvis[i] = sign * self.uedg[ibl][is];
        }
    }
}
```

---

### GAMQV
**File:** `xpanel.f:1616`  
**Purpose:** Computes circulation (vorticity) from velocity.

**Rust Implementation:**
```rust
pub fn gamqv(&mut self) -> Result<(), Error> {
    for i in 0..self.n {
        self.gam[i] = self.qvis[i];
    }
}
```

---

## Newton Solver

### BLSOLV
**File:** `xsolve.f:283`  
**Purpose:** Custom block-tridiagonal solver for coupled viscous-inviscid system.

**System Structure:**
```
[A B | | ... | | B A]  [dCtau]   [R1]
[  B A | ... | | B A]  [dθ   ] = [R2]
[    B A ... | | B A]  [dm   ]   [R3]
[      ... ... ... ]   [ ... ]   [...]
```

**Algorithm:**
1. Forward elimination (block LU)
2. Back substitution

**Rust Implementation:**
```rust
pub fn blsolv(&mut self) -> Result<(), Error> {
    let nsys = self.nsys;
    
    // Forward elimination
    for iv in 0..nsys {
        // Invert VA[iv] block
        self.invert_block(&mut self.va[iv])?;
        
        // Eliminate VB[iv+1] block
        if iv < nsys - 1 {
            self.eliminate_block(&mut self.vb[iv+1], &self.va[iv])?;
        }
        
        // Eliminate coupling columns
        self.eliminate_coupling(iv)?;
    }
    
    // Back substitution
    for iv in (0..nsys).rev() {
        self.back_substitute(iv)?;
    }
}
```

---

## Closure Relations

### HKIN
**File:** `xblsys.f:2276`  
**Purpose:** Computes kinematic shape factor Hk from shape factor H and Mach number.

**Rust Implementation:**
```rust
pub fn hkin(&self, h: f64, msq: f64) -> Result<f64, Error> {
    // Compressibility correction
    let beta = (1.0 - msq).max(0.01).sqrt();
    let hk = h / beta;
    Ok(hk)
}
```

---

### CFL
**File:** `xblsys.f:2354`  
**Purpose:** Laminar skin friction coefficient.

**Rust Implementation:**
```rust
pub fn cfl(&self, hk: f64, rt: f64, msq: f64) -> Result<f64, Error> {
    // Falkner-Skan correlation
    let cf = 0.3 * (-1.33 * hk).exp() 
             * (rt.ln() / 2.3).powf(-1.74 - 0.31 * hk);
    Ok(cf)
}
```

---

### CFT
**File:** `xblsys.f:2483`  
**Purpose:** Turbulent skin friction coefficient.

**Rust Implementation:**
```rust
pub fn cft(&self, hk: f64, rt: f64, msq: f64) -> Result<f64, Error> {
    // Coles wall law
    let cf = self.cffac * 0.3 * (-1.33 * hk).exp()
             * (rt.ln() / 2.3).powf(-1.74 - 0.31 * hk);
    Ok(cf)
}
```

---

### HSL
**File:** `xblsys.f:2327`  
**Purpose:** Laminar energy shape factor.

**Rust Implementation:**
```rust
pub fn hsl(&self, hk: f64, rt: f64, msq: f64) -> Result<f64, Error> {
    // Correlation from Drela
    let hs = 2.0 + 1.5 / hk + 0.5 / (hk - 1.0);
    Ok(hs)
}
```

---

### HST
**File:** `xblsys.f:2388`  
**Purpose:** Turbulent energy shape factor.

**Rust Implementation:**
```rust
pub fn hst(&self, hk: f64, rt: f64, msq: f64) -> Result<f64, Error> {
    // Correlation from Drela
    let hs = 1.5 + 1.0 / (hk - 1.0) + 0.5 / (hk - 1.0).powi(2);
    Ok(hs)
}
```

---

### DIL
**File:** `xblsys.f:2290`  
**Purpose:** Laminar dissipation integral.

**Rust Implementation:**
```rust
pub fn dil(&self, hk: f64, rt: f64) -> Result<f64, Error> {
    // Correlation
    let di = 0.5 * cf * ue.powi(3) * theta;
    Ok(di)
}
```

---

### DILW
**File:** `xblsys.f:2308`  
**Purpose:** Wake dissipation integral.

**Rust Implementation:**
```rust
pub fn dilw(&self, hk: f64, rt: f64) -> Result<f64, Error> {
    // Wake-specific correlation
    let di = self.dil(hk, rt)? * self.dlcon;  // Reduced dissipation
    Ok(di)
}
```

---

### DIT
**File:** `xblsys.f:2375`  
**Purpose:** Turbulent dissipation integral.

**Rust Implementation:**
```rust
pub fn dit(&self, hs: f64, us: f64, cf: f64, st: f64) -> Result<f64, Error> {
    // Turbulent dissipation
    let di = 0.5 * cf * us.powi(3) * theta * (1.0 + hs);
    Ok(di)
}
```

---

### DAMPL
**File:** `xblsys.f:1981`  
**Purpose:** Amplification rate for transition prediction.

**Rust Implementation:**
```rust
pub fn dampl(&self, hk: f64, th: f64, rt: f64) -> Result<f64, Error> {
    // Amplification rate correlation
    let ax = self.amplification_rate(hk, rt)?;
    Ok(ax)
}
```

---

### DAMPL2
**File:** `xblsys.f:2099`  
**Purpose:** Second-order amplification rate.

**Rust Implementation:**
```rust
pub fn dampl2(&self, hk: f64, th: f64, rt: f64) -> Result<f64, Error> {
    // Improved amplification rate
    let ax = self.amplification_rate_2nd_order(hk, rt)?;
    Ok(ax)
}
```

---

## Transition Prediction

### TRCHEK
**File:** `xblsys.f:22`  
**Purpose:** Checks for transition using eN method.

**Algorithm:**
```
1. Integrate amplification: dN/ds = σ(Reθ, H)
2. Check if N >= Ncrit
3. Set transition flag
```

**Rust Implementation:**
```rust
pub fn trchek(&mut self, station: usize, surface: usize) -> Result<(), Error> {
    // Integrate amplification
    let ax = self.axset(station, surface)?;
    self.ampl[station][surface] += ax * self.ds;
    
    // Check transition
    if self.ampl[station][surface] >= self.acrit[surface] {
        self.trans = true;
        self.itran[surface] = station;
    }
}
```

**Dependencies:**
- AXSET
- DAMPL or DAMPL2

---

### TRCHEK2
**File:** `xblsys.f:231`  
**Purpose:** Second-order transition check.

**Rust Implementation:**
```rust
pub fn trchek2(&mut self, station: usize, surface: usize) -> Result<(), Error> {
    // Use second-order amplification
    let ax = self.axset_2nd_order(station, surface)?;
    self.ampl[station][surface] += ax * self.ds;
    
    // Check transition
    if self.ampl[station][surface] >= self.acrit[surface] {
        self.trans = true;
    }
}
```

---

### AXSET
**File:** `xblsys.f:35`  
**Purpose:** Computes average amplification over interval.

**Rust Implementation:**
```rust
pub fn axset(
    &self,
    hk1: f64, t1: f64, rt1: f64, a1: f64,
    hk2: f64, t2: f64, rt2: f64, a2: f64,
    acrit: f64,
) -> Result<f64, Error> {
    // Compute amplification at stations 1 and 2
    let ax1 = self.dampl(hk1, t1, rt1)?;
    let ax2 = self.dampl(hk2, t2, rt2)?;
    
    // RMS average
    let axsq = 0.5 * (ax1.powi(2) + ax2.powi(2));
    let axa = if axsq > 0.0 { axsq.sqrt() } else { 0.0 };
    
    Ok(axa)
}
```

---

## Utility Functions

### IBLSYS
**File:** `xbl.f:519`  
**Purpose:** Maps BL stations to system indices.

**Rust Implementation:**
```rust
pub fn iblsys(&mut self) -> Result<(), Error> {
    let mut iv = 0;
    for is in 0..2 {
        for ibl in 2..self.nbl[is] {
            self.isys[ibl][is] = iv;
            iv += 1;
        }
    }
    self.nsys = iv;
}
```

---

### IBLPAN
**File:** `xpanel.f:1395`  
**Purpose:** Maps BL stations to panel indices.

**Rust Implementation:**
```rust
pub fn iblpan(&mut self) -> Result<(), Error> {
    for is in 0..2 {
        for ibl in 2..self.nbl[is] {
            // Find panel closest to BL station
            let s_bl = self.xssi[ibl][is];
            let ipan = self.find_panel_for_s(s_bl, is)?;
            self.ipan[ibl][is] = ipan;
        }
    }
}
```

---

### XICALC
**File:** `xpanel.f:1455`  
**Purpose:** Computes arc length array for current stagnation point.

**Rust Implementation:**
```rust
pub fn xicalc(&mut self) -> Result<(), Error> {
    // Compute arc length from stagnation point
    let s_stag = self.s_stag;
    for is in 0..2 {
        for ibl in 2..self.nbl[is] {
            let i = self.ipan[ibl][is];
            self.xssi[ibl][is] = (self.s[i] - s_stag).abs();
        }
    }
}
```

---

### STFIND
**File:** `xpanel.f:1357`  
**Purpose:** Finds stagnation point location.

**Rust Implementation:**
```rust
pub fn stfind(&mut self) -> Result<(), Error> {
    // Find where gamma changes sign
    for i in 0..self.n {
        if self.gam[i] > 0.0 {
            // Stagnation point between i-1 and i
            self.istag = [i-1, i];
            // Interpolate
            let s_stag = self.interpolate_stag(i-1, i)?;
            self.s_stag = s_stag;
            break;
        }
    }
}
```

---

### STMOVE
**File:** `xpanel.f:1628`  
**Purpose:** Moves stagnation point based on viscous solution.

**Rust Implementation:**
```rust
pub fn stmove(&mut self) -> Result<(), Error> {
    // Update stagnation point location based on Ue
    let ue_upper = self.uedg[2][1];
    let ue_lower = self.uedg[2][0];
    
    // Stagnation moves based on velocity difference
    let ds_stag = self.compute_stag_shift(ue_upper, ue_lower)?;
    self.s_stag += ds_stag;
    
    // Update arc lengths
    self.xicalc()?;
}
```

---

### DSLIM
**File:** `xbl.f:1564`  
**Purpose:** Limits displacement thickness to prevent blow-up.

**Rust Implementation:**
```rust
pub fn dslim(&mut self, dstr: &mut f64, thet: f64, uedg: f64, msq: f64, hklim: f64) -> Result<(), Error> {
    let hk = *dstr / thet;
    if hk > hklim {
        *dstr = hklim * thet;
    }
}
```

---

### BLPINI
**File:** `xbl.f:1578`  
**Purpose:** Initializes BL variables.

**Rust Implementation:**
```rust
pub fn blpini(&mut self) -> Result<(), Error> {
    // Initialize all BL arrays
    for is in 0..2 {
        for ibl in 0..self.nbl[is] {
            self.ctau[ibl][is] = 0.0;
            self.thet[ibl][is] = 1e-5;
            self.dstr[ibl][is] = 2.2e-5;
            self.mass[ibl][is] = self.uedg[ibl][is] * self.dstr[ibl][is];
        }
    }
}
```

---

### XIFSET
**File:** `xbl.f:1196`  
**Purpose:** Sets forced transition location.

**Rust Implementation:**
```rust
pub fn xifset(&mut self, is: usize) -> Result<(), Error> {
    if self.xstrip[is] > 0.0 {
        // Find station closest to forced transition location
        for ibl in 2..self.nbl[is] {
            if self.xssi[ibl][is] >= self.xstrip[is] {
                self.itran[is] = ibl;
                break;
            }
        }
    }
}
```

---

## Constants and Parameters

### BLPAR.INC
**File:** `BLPAR.INC`  
**Purpose:** Boundary layer parameter constants.

**Constants:**
- `SCCON = 5.6`: Shear coefficient lag constant
- `GACON = 6.70`: G-β locus constant A
- `GBCON = 0.75`: G-β locus constant B
- `GCCON = 18.0`: G-β wall term constant
- `DLCON = 0.9`: Wall/wake dissipation length ratio
- `CTRCON = 1.8`: Ctau constant
- `CTRCEX = 3.3`: Ctau exponent
- `DUXCON`: Velocity gradient constant
- `CTCON`: Ctau weighting coefficient
- `CFFAC = 1.0`: Skin friction factor

**Rust Implementation:**
```rust
pub struct BLParams {
    pub sccon: f64,   // 5.6
    pub gacon: f64,   // 6.70
    pub gbcon: f64,   // 0.75
    pub gccon: f64,   // 18.0
    pub dlcon: f64,   // 0.9
    pub ctrcon: f64,  // 1.8
    pub ctrcex: f64,  // 3.3
    pub duxcon: f64,
    pub ctcon: f64,
    pub cffac: f64,   // 1.0
}

impl Default for BLParams {
    fn default() -> Self {
        Self {
            sccon: 5.6,
            gacon: 6.70,
            gbcon: 0.75,
            gccon: 18.0,
            dlcon: 0.9,
            ctrcon: 1.8,
            ctrcex: 3.3,
            duxcon: 0.0,
            ctcon: 0.0,
            cffac: 1.0,
        }
    }
}
```

---

## Implementation Checklist

### Phase 1: Core Data Structures
- [ ] BL state arrays (CTAU, THET, DSTR, MASS, UEDG)
- [ ] Panel-BL mapping (IPAN, ISYS)
- [ ] DIJ matrix
- [ ] Newton system arrays (VA, VB, VM, VDEL)

### Phase 2: Closure Relations
- [ ] HKIN
- [ ] CFL, CFT
- [ ] HSL, HST, HCT
- [ ] DIL, DILW, DIT
- [ ] DAMPL, DAMPL2

### Phase 3: BL System
- [ ] BLPRV, BLKIN
- [ ] BLVAR, BLMID
- [ ] BLDIF, TRDIF
- [ ] BLSYS, TESYS

### Phase 4: BL Marching
- [ ] MRCHUE (direct/inverse mode)
- [ ] MRCHDU (transition)
- [ ] TRCHEK, AXSET

### Phase 5: Panel Coupling
- [ ] QDCALC (DIJ matrix)
- [ ] UESET (Ue from mass defect)
- [ ] DSSET (δ* from mass defect)
- [ ] UICALC, QVFUE, GAMQV

### Phase 6: Newton Solver
- [ ] SETBL (build system)
- [ ] BLSOLV (solve system)
- [ ] UPDATE (apply updates)

### Phase 7: Main Loop
- [ ] VISCAL (main iteration)
- [ ] STFIND, STMOVE
- [ ] IBLPAN, XICALC, IBLSYS

### Phase 8: Utilities
- [ ] DSLIM
- [ ] BLPINI
- [ ] XIFSET

---

## Function Dependency Graph

```
VISCAL
  ├─ XYWAKE
  ├─ QWCALC
  ├─ QISET
  ├─ STFIND
  ├─ IBLPAN
  ├─ XICALC
  ├─ IBLSYS
  ├─ UICALC
  ├─ QDCALC
  │   └─ PSILIN
  ├─ SETBL
  │   ├─ MRCHUE
  │   │   ├─ BLSYS
  │   │   │   ├─ BLPRV
  │   │   │   ├─ BLKIN
  │   │   │   │   └─ HKIN
  │   │   │   ├─ BLVAR
  │   │   │   │   ├─ HSL, HST, HCT
  │   │   │   │   ├─ CFL, CFT
  │   │   │   │   └─ DIL, DILW, DIT
  │   │   │   ├─ BLMID
  │   │   │   ├─ BLDIF
  │   │   │   └─ TRDIF
  │   │   └─ TESYS
  │   └─ MRCHDU
  │       └─ TRCHEK
  │           └─ AXSET
  │               └─ DAMPL/DAMPL2
  ├─ BLSOLV
  ├─ UPDATE
  │   ├─ UESET
  │   ├─ DSSET
  │   └─ DSLIM
  ├─ QVFUE
  ├─ GAMQV
  ├─ STMOVE
  ├─ CLCALC
  └─ CDCALC
```

---

**This document provides a complete reference for implementing FlexFoil