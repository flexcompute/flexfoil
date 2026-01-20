# FlexFoil Implementation Specification

## Purpose
This document provides exact Rust function signatures and implementation templates for recreating XFOIL's viscous solver. Each function corresponds to an XFOIL Fortran routine.

---

## Module: `closures.rs` - Boundary Layer Closure Relations

### Constants

```rust
/// Boundary layer model parameters (BLPAR.INC equivalent)
pub struct BlParams {
    pub sccon: f64,   // 5.6 - Shear lag constant
    pub gacon: f64,   // 6.70 - G-beta constant
    pub gbcon: f64,   // 0.75 - G-beta constant  
    pub gccon: f64,   // 18.0 - G-beta wall term
    pub dlcon: f64,   // 0.9 - Wake dissipation ratio
    pub ctrcon: f64,  // 1.8 - Ctau transition constant
    pub ctrcex: f64,  // 3.3 - Ctau transition exponent
    pub duxcon: f64,  // 1.0 - dUe/dx weighting
    pub ctcon: f64,   // 0.5/(gacon² × gbcon) - Ctau equilibrium
    pub cffac: f64,   // 1.0 - Cf correction factor
}

impl Default for BlParams {
    fn default() -> Self {
        let gacon = 6.70;
        let gbcon = 0.75;
        Self {
            sccon: 5.6,
            gacon,
            gbcon,
            gccon: 18.0,
            dlcon: 0.9,
            ctrcon: 1.8,
            ctrcex: 3.3,
            duxcon: 1.0,
            ctcon: 0.5 / (gacon * gacon * gbcon),
            cffac: 1.0,
        }
    }
}
```

### `hkin` - Kinematic Shape Parameter (HKIN)

```rust
/// Kinematic shape parameter Hk from H and M² (Whitfield)
/// 
/// XFOIL: xblsys.f:2276
/// 
/// # Arguments
/// * `h` - Shape parameter δ*/θ
/// * `msq` - Edge Mach number squared
/// 
/// # Returns
/// (Hk, dHk/dH, dHk/dM²)
pub fn hkin(h: f64, msq: f64) -> (f64, f64, f64) {
    let hk = (h - 0.29 * msq) / (1.0 + 0.113 * msq);
    let hk_h = 1.0 / (1.0 + 0.113 * msq);
    let hk_msq = (-0.29 - 0.113 * hk) / (1.0 + 0.113 * msq);
    (hk, hk_h, hk_msq)
}
```

### `hsl` - Laminar Energy Shape Parameter (HSL)

```rust
/// Laminar H* correlation
/// 
/// XFOIL: xblsys.f:2327
/// 
/// # Arguments
/// * `hk` - Kinematic shape parameter
/// * `ret` - Momentum thickness Reynolds number (unused in laminar)
/// * `msq` - Edge Mach number squared (unused in laminar)
/// 
/// # Returns
/// (H*, dH*/dHk, dH*/dRet, dH*/dM²)
pub fn hsl(hk: f64, _ret: f64, _msq: f64) -> (f64, f64, f64, f64) {
    let (hs, hs_hk) = if hk < 4.35 {
        let tmp = hk - 4.35;
        let hs = 0.0111 * tmp * tmp / (hk + 1.0)
               - 0.0278 * tmp.powi(3) / (hk + 1.0)
               + 1.528
               - 0.0002 * (tmp * hk).powi(2);
        let hs_hk = 0.0111 * (2.0 * tmp - tmp * tmp / (hk + 1.0)) / (hk + 1.0)
                  - 0.0278 * (3.0 * tmp * tmp - tmp.powi(3) / (hk + 1.0)) / (hk + 1.0)
                  - 0.0002 * 2.0 * tmp * hk * (tmp + hk);
        (hs, hs_hk)
    } else {
        let hs2 = 0.015;
        let hs = hs2 * (hk - 4.35).powi(2) / hk + 1.528;
        let hs_hk = hs2 * 2.0 * (hk - 4.35) / hk
                  - hs2 * (hk - 4.35).powi(2) / (hk * hk);
        (hs, hs_hk)
    };
    
    (hs, hs_hk, 0.0, 0.0)
}
```

### `hst` - Turbulent Energy Shape Parameter (HST)

```rust
/// Turbulent H* correlation with Re_θ dependence
/// 
/// XFOIL: xblsys.f:2388
/// 
/// # Arguments
/// * `hk` - Kinematic shape parameter
/// * `ret` - Momentum thickness Reynolds number
/// * `msq` - Edge Mach number squared
/// 
/// # Returns
/// (H*, dH*/dHk, dH*/dRet, dH*/dM²)
pub fn hst(hk: f64, ret: f64, msq: f64) -> (f64, f64, f64, f64) {
    const HSMIN: f64 = 1.500;
    const DHSINF: f64 = 0.015;
    
    // Limited Re_θ dependence for Re_θ < 200
    let (rtz, rtz_rt) = if ret > 200.0 {
        (ret, 1.0)
    } else {
        (200.0, 0.0)
    };
    
    let (ho, ho_rt) = if ret > 400.0 {
        (3.0 + 400.0 / ret, -400.0 / (ret * ret))
    } else {
        (4.0, 0.0)
    };
    
    let (hs, hs_hk, hs_rt) = if hk < ho {
        // Attached branch
        let hr = (ho - hk) / (ho - 1.0);
        let hr_hk = -1.0 / (ho - 1.0);
        let hr_rt = (1.0 - hr) / (ho - 1.0) * ho_rt;
        
        let hs = (2.0 - HSMIN - 4.0 / rtz) * hr * hr * 1.5 / (hk + 0.5) + HSMIN + 4.0 / rtz;
        let hs_hk = -(2.0 - HSMIN - 4.0 / rtz) * hr * hr * 1.5 / (hk + 0.5).powi(2)
                  + (2.0 - HSMIN - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hr_hk;
        let hs_rt = (2.0 - HSMIN - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hr_rt
                  + (hr * hr * 1.5 / (hk + 0.5) - 1.0) * 4.0 / (rtz * rtz) * rtz_rt;
        
        (hs, hs_hk, hs_rt)
    } else {
        // Separated branch
        let grt = rtz.ln();
        let hdif = hk - ho;
        let rtmp = hk - ho + 4.0 / grt;
        let htmp = 0.007 * grt / (rtmp * rtmp) + DHSINF / hk;
        let htmp_hk = -0.014 * grt / rtmp.powi(3) - DHSINF / (hk * hk);
        let htmp_rt = -0.014 * grt / rtmp.powi(3) * (-ho_rt - 4.0 / (grt * grt * rtz) * rtz_rt)
                    + 0.007 / (rtmp * rtmp * rtz) * rtz_rt;
        
        let hs = hdif * hdif * htmp + HSMIN + 4.0 / rtz;
        let hs_hk = hdif * 2.0 * htmp + hdif * hdif * htmp_hk;
        let hs_rt = hdif * hdif * htmp_rt - 4.0 / (rtz * rtz) * rtz_rt
                  + hdif * 2.0 * htmp * (-ho_rt);
        
        (hs, hs_hk, hs_rt)
    };
    
    // Compressibility correction (Whitfield)
    let fm = 1.0 + 0.014 * msq;
    let hs_out = (hs + 0.028 * msq) / fm;
    let hs_hk_out = hs_hk / fm;
    let hs_rt_out = hs_rt / fm;
    let hs_msq = 0.028 / fm - 0.014 * hs_out / fm;
    
    (hs_out, hs_hk_out, hs_rt_out, hs_msq)
}
```

### `cfl` - Laminar Skin Friction (CFL)

```rust
/// Laminar skin friction coefficient (Falkner-Skan)
/// 
/// XFOIL: xblsys.f:2354
/// 
/// # Arguments
/// * `hk` - Kinematic shape parameter
/// * `ret` - Momentum thickness Reynolds number
/// * `msq` - Edge Mach number squared (unused)
/// 
/// # Returns
/// (Cf, dCf/dHk, dCf/dRet, dCf/dM²)
pub fn cfl(hk: f64, ret: f64, _msq: f64) -> (f64, f64, f64, f64) {
    let (cf, cf_hk) = if hk < 5.5 {
        let tmp = (5.5 - hk).powi(3) / (hk + 1.0);
        let cf = (0.0727 * tmp - 0.07) / ret;
        let cf_hk = (-0.0727 * tmp * 3.0 / (5.5 - hk) - 0.0727 * tmp / (hk + 1.0)) / ret;
        (cf, cf_hk)
    } else {
        let tmp = 1.0 - 1.0 / (hk - 4.5);
        let cf = (0.015 * tmp * tmp - 0.07) / ret;
        let cf_hk = (0.015 * tmp * 2.0 / (hk - 4.5).powi(2)) / ret;
        (cf, cf_hk)
    };
    
    let cf_rt = -cf / ret;
    (cf, cf_hk, cf_rt, 0.0)
}
```

### `cft` - Turbulent Skin Friction (CFT)

```rust
/// Turbulent skin friction coefficient (Coles)
/// 
/// XFOIL: xblsys.f:2483
/// 
/// # Arguments
/// * `hk` - Kinematic shape parameter
/// * `ret` - Momentum thickness Reynolds number
/// * `msq` - Edge Mach number squared
/// * `params` - BL parameters (for cffac)
/// 
/// # Returns
/// (Cf, dCf/dHk, dCf/dRet, dCf/dM²)
pub fn cft(hk: f64, ret: f64, msq: f64, cffac: f64) -> (f64, f64, f64, f64) {
    const GAM: f64 = 1.4;
    let gm1 = GAM - 1.0;
    
    let fc = (1.0 + 0.5 * gm1 * msq).sqrt();
    let grt = (ret / fc).ln().max(3.0);
    
    let gex = -1.74 - 0.31 * hk;
    let arg = (-1.33 * hk).max(-20.0);
    let thk = (4.0 - hk / 0.875).tanh();
    
    let cfo = cffac * 0.3 * arg.exp() * (grt / 2.3026_f64.ln()).powf(gex);
    let cf = (cfo + 1.1e-4 * (thk - 1.0)) / fc;
    
    let cf_hk = (-1.33 * cfo - 0.31 * (grt / 2.3026_f64.ln()).ln() * cfo
               - 1.1e-4 * (1.0 - thk * thk) / 0.875) / fc;
    let cf_rt = gex * cfo / (fc * grt) / ret;
    let cf_msq = gex * cfo / (fc * grt) * (-0.25 * gm1 / (fc * fc)) 
               - 0.25 * gm1 * cf / (fc * fc);
    
    (cf, cf_hk, cf_rt, cf_msq)
}
```

### `dil` - Laminar Dissipation (DIL)

```rust
/// Laminar dissipation function 2CD/H* (Falkner-Skan)
/// 
/// XFOIL: xblsys.f:2290
/// 
/// # Arguments
/// * `hk` - Kinematic shape parameter
/// * `ret` - Momentum thickness Reynolds number
/// 
/// # Returns
/// (2CD/H*, d(2CD/H*)/dHk, d(2CD/H*)/dRet)
pub fn dil(hk: f64, ret: f64) -> (f64, f64, f64) {
    let (di, di_hk) = if hk < 4.0 {
        let di = (0.00205 * (4.0 - hk).powf(5.5) + 0.207) / ret;
        let di_hk = (-0.00205 * 5.5 * (4.0 - hk).powf(4.5)) / ret;
        (di, di_hk)
    } else {
        let hkb = hk - 4.0;
        let den = 1.0 + 0.02 * hkb * hkb;
        let di = (-0.0016 * hkb * hkb / den + 0.207) / ret;
        let di_hk = (-0.0016 * 2.0 * hkb * (1.0 / den - 0.02 * hkb * hkb / (den * den))) / ret;
        (di, di_hk)
    };
    
    let di_rt = -di / ret;
    (di, di_hk, di_rt)
}
```

### `hct` - Density Shape Parameter (HCT)

```rust
/// Density shape parameter H** (Whitfield)
/// 
/// XFOIL: xblsys.f:2514
/// 
/// # Arguments
/// * `hk` - Kinematic shape parameter
/// * `msq` - Edge Mach number squared
/// 
/// # Returns
/// (H**, dH**/dHk, dH**/dM²)
pub fn hct(hk: f64, msq: f64) -> (f64, f64, f64) {
    let hc = msq * (0.064 / (hk - 0.8) + 0.251);
    let hc_hk = msq * (-0.064 / (hk - 0.8).powi(2));
    let hc_msq = 0.064 / (hk - 0.8) + 0.251;
    (hc, hc_hk, hc_msq)
}
```

---

## Module: `transition.rs` - Transition Prediction

### `dampl` - Amplification Rate (DAMPL)

```rust
/// Envelope amplification rate dN/dx for e^N method
/// 
/// XFOIL: xblsys.f:1981
/// 
/// Reference: Drela & Giles, AIAA J. Oct 1987
/// 
/// # Arguments
/// * `hk` - Kinematic shape parameter
/// * `theta` - Momentum thickness
/// * `ret` - Momentum thickness Reynolds number
/// 
/// # Returns
/// (AX, dAX/dHk, dAX/dθ, dAX/dRet)
pub fn dampl(hk: f64, theta: f64, ret: f64) -> (f64, f64, f64, f64) {
    const DGR: f64 = 0.08;
    
    let hmi = 1.0 / (hk - 1.0);
    let hmi_hk = -hmi * hmi;
    
    // Critical Re_θ correlation
    let aa = 2.492 * hmi.powf(0.43);
    let aa_hk = (aa / hmi) * 0.43 * hmi_hk;
    
    let bb = (14.0 * hmi - 9.24).tanh();
    let bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;
    
    let grcrit = aa + 0.7 * (bb + 1.0);
    let grc_hk = aa_hk + 0.7 * bb_hk;
    
    let gr = ret.log10();
    let gr_rt = 1.0 / (2.302585 * ret);
    
    if gr < grcrit - DGR {
        // Below critical: no amplification
        return (0.0, 0.0, 0.0, 0.0);
    }
    
    // Ramp function for smooth turn-on
    let rnorm = (gr - (grcrit - DGR)) / (2.0 * DGR);
    let rn_hk = -grc_hk / (2.0 * DGR);
    let rn_rt = gr_rt / (2.0 * DGR);
    
    let (rfac, rfac_hk, rfac_rt) = if rnorm >= 1.0 {
        (1.0, 0.0, 0.0)
    } else {
        let rfac = 3.0 * rnorm * rnorm - 2.0 * rnorm.powi(3);
        let rfac_rn = 6.0 * rnorm - 6.0 * rnorm * rnorm;
        (rfac, rfac_rn * rn_hk, rfac_rn * rn_rt)
    };
    
    // Amplification envelope slope
    let arg = 3.87 * hmi - 2.52;
    let arg_hk = 3.87 * hmi_hk;
    
    let ex = (-arg * arg).exp();
    let ex_hk = ex * (-2.0 * arg * arg_hk);
    
    let dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
    let dadr_hk = 0.028 - 0.0345 * ex_hk;
    
    // m(H) correlation
    let af = -0.05 + 2.7 * hmi - 5.5 * hmi * hmi + 3.0 * hmi.powi(3);
    let af_hmi = 2.7 - 11.0 * hmi + 9.0 * hmi * hmi;
    let af_hk = af_hmi * hmi_hk;
    
    // Final amplification rate
    let ax = (af * dadr / theta) * rfac;
    let ax_hk = (af_hk * dadr / theta + af * dadr_hk / theta) * rfac
              + (af * dadr / theta) * rfac_hk;
    let ax_th = -ax / theta;
    let ax_rt = (af * dadr / theta) * rfac_rt;
    
    (ax, ax_hk, ax_th, ax_rt)
}
```

### `trchek2` - Transition Check (TRCHEK2)

```rust
/// Second-order transition check
/// 
/// XFOIL: xblsys.f:231
/// 
/// Solves implicit amplification equation and locates transition point.
/// 
/// # Arguments
/// * `st1` - Upstream station
/// * `st2` - Downstream station  
/// * `ncrit` - Critical amplification (typically 9.0)
/// 
/// # Returns
/// Transition result with location and sensitivities
pub struct TransitionResult {
    pub occurs: bool,         // Did transition occur?
    pub forced: bool,         // Was it forced (trip)?
    pub xt: f64,             // Transition location
    pub ampl2: f64,          // Amplification at station 2
    // Sensitivities...
    pub xt_a1: f64,
    pub xt_t1: f64,
    pub xt_d1: f64,
    pub xt_u1: f64,
    // ... etc
}

pub fn trchek2(
    st1: &BlStation,
    st2: &BlStation,
    ncrit: f64,
    xiforc: f64,  // Forced transition location
) -> TransitionResult {
    const DAEPS: f64 = 5.0e-5;
    
    // Calculate average amplification rate
    let (ax, ax_hk1, ax_t1, ax_rt1, ax_a1,
         ax_hk2, ax_t2, ax_rt2, ax_a2) = axset(
        st1.hk, st1.theta, st1.ret, st1.ampl,
        st2.hk, st2.theta, st2.ret, st2.ampl,
        ncrit
    );
    
    // Initial guess for N2
    let mut ampl2 = st1.ampl + ax * (st2.x - st1.x);
    
    // Newton iteration for implicit N2
    for _ in 0..30 {
        // Determine weighting factors
        let (amplt, amplt_a2, sfa, sfa_a1, sfa_a2) = if ampl2 <= ncrit {
            // No transition yet
            (ampl2, 1.0, 1.0, 0.0, 0.0)
        } else {
            // Transition in interval
            let amplt = ncrit;
            let sfa = (amplt - st1.ampl) / (ampl2 - st1.ampl);
            (amplt, 0.0, sfa, (sfa - 1.0) / (ampl2 - st1.ampl), -sfa / (ampl2 - st1.ampl))
        };
        
        // ... interpolation and iteration logic
        
        // Check convergence
        let res = ampl2 - st1.ampl - ax * (st2.x - st1.x);
        if res.abs() < DAEPS {
            break;
        }
    }
    
    // Determine if transition occurred
    let trfree = ampl2 >= ncrit;
    let trforc = xiforc > st1.x && xiforc <= st2.x;
    
    TransitionResult {
        occurs: trfree || trforc,
        forced: trforc && (!trfree || xiforc < /* xt from free */),
        xt: if trforc { xiforc } else { /* interpolated xt */ },
        ampl2,
        // ... sensitivities
        ..Default::default()
    }
}
```

---

## Module: `equations.rs` - BL Equations

### Equation System Structure

```rust
/// Residuals and Jacobians for one BL interval
pub struct BlEquations {
    /// Residuals: [shear_lag/amp, momentum, shape]
    pub res: [f64; 3],
    
    /// Jacobian wrt upstream station: d(res)/d[ctau, theta, dstar, ue, xi]
    pub jac1: [[f64; 5]; 3],
    
    /// Jacobian wrt downstream station
    pub jac2: [[f64; 5]; 3],
    
    /// Mach sensitivity
    pub res_msq: [f64; 3],
    
    /// Reynolds sensitivity  
    pub res_re: [f64; 3],
}

/// Equation type indicator
#[derive(Clone, Copy, PartialEq)]
pub enum BlType {
    Similarity,  // ITYP = 0
    Laminar,     // ITYP = 1
    Turbulent,   // ITYP = 2
    Wake,        // ITYP = 3
}
```

### `bldif` - Finite Difference Equations (BLDIF)

```rust
/// Set up finite-difference Newton coefficients
/// 
/// XFOIL: xblsys.f:1552
/// 
/// # Equations
/// 1. Laminar: dN/dx = AX(Hk, θ, Rθ)
///    Turbulent: Shear lag equation
/// 2. Momentum: d(ln θ)/dx + (H+2-M)·d(ln Ue)/dx = Cf·x/(2θ)
/// 3. Shape: d(ln H*)/dx + (2H**/H*+1-H)·d(ln Ue)/dx = (Cf/2 - 2CD/H*)·x/θ
pub fn bldif(
    st1: &BlStation,
    st2: &BlStation,
    bl_type: BlType,
    params: &BlParams,
) -> BlEquations {
    let mut eq = BlEquations::default();
    
    let dx = st2.x - st1.x;
    
    // Logarithmic differences (or similarity special case)
    let (xlog, ulog, tlog, hlog) = if bl_type == BlType::Similarity {
        (1.0, params.bule, 0.5 * (1.0 - params.bule), 0.0)
    } else {
        ((st2.x / st1.x).ln(),
         (st2.ue / st1.ue).ln(),
         (st2.theta / st1.theta).ln(),
         (st2.hs / st1.hs).ln())
    };
    
    // Local upwinding based on Hk change
    let upw = compute_upwinding(st1.hk, st2.hk, bl_type);
    
    match bl_type {
        BlType::Similarity => {
            // Equation 1: Zero amplification
            eq.res[0] = -st2.ampl;
            eq.jac2[0][0] = 1.0;
        }
        
        BlType::Laminar => {
            // Equation 1: Amplification
            let (ax, ax_hk1, ax_t1, ax_rt1, ax_a1,
                 ax_hk2, ax_t2, ax_rt2, ax_a2) = axset(
                st1.hk, st1.theta, st1.ret, st1.ampl,
                st2.hk, st2.theta, st2.ret, st2.ampl,
                params.ncrit,
            );
            
            eq.res[0] = -(st2.ampl - st1.ampl - ax * dx);
            
            let z_ax = -dx;
            eq.jac1[0][0] = z_ax * ax_a1 - 1.0;
            eq.jac1[0][1] = z_ax * (ax_hk1 * st1.hk_t + ax_t1 + ax_rt1 * st1.rt_t);
            eq.jac1[0][2] = z_ax * ax_hk1 * st1.hk_d;
            eq.jac1[0][3] = z_ax * (ax_hk1 * st1.hk_u + ax_rt1 * st1.rt_u);
            eq.jac1[0][4] = ax;
            
            eq.jac2[0][0] = z_ax * ax_a2 + 1.0;
            eq.jac2[0][1] = z_ax * (ax_hk2 * st2.hk_t + ax_t2 + ax_rt2 * st2.rt_t);
            eq.jac2[0][2] = z_ax * ax_hk2 * st2.hk_d;
            eq.jac2[0][3] = z_ax * (ax_hk2 * st2.hk_u + ax_rt2 * st2.rt_u);
            eq.jac2[0][4] = -ax;
        }
        
        BlType::Turbulent | BlType::Wake => {
            // Equation 1: Shear lag
            let ald = if bl_type == BlType::Wake { params.dlcon } else { 1.0 };
            
            // Upwind-averaged quantities
            let sa = (1.0 - upw) * st1.ctau + upw * st2.ctau;
            let cqa = (1.0 - upw) * st1.cq + upw * st2.cq;
            let cfa = (1.0 - upw) * st1.cf + upw * st2.cf;
            let hka = (1.0 - upw) * st1.hk + upw * st2.hk;
            
            // Simple averages
            let usa = 0.5 * (st1.us + st2.us);
            let rta = 0.5 * (st1.ret + st2.ret);
            let dea = 0.5 * (st1.de + st2.de);
            let da = 0.5 * (st1.dstar + st2.dstar);
            
            // Equilibrium dUe/dx (NEW correlation 12 Oct 94)
            let gcc = if bl_type == BlType::Turbulent { params.gccon } else { 0.0 };
            let hkc = (hka - 1.0 - gcc / rta).max(0.01);
            let hr = hkc / (params.gacon * ald * hka);
            let uq = (0.5 * cfa - hr * hr) / (params.gbcon * da);
            
            let scc = params.sccon * 1.333 / (1.0 + usa);
            let slog = (st2.ctau / st1.ctau).ln();
            
            eq.res[0] = -(scc * (cqa - sa * ald) * dx
                       - dea * 2.0 * slog
                       + dea * 2.0 * (uq * dx - ulog) * params.duxcon);
            
            // ... Jacobian entries (extensive)
        }
    }
    
    // Equation 2: Momentum integral
    {
        let ha = 0.5 * (st1.h + st2.h);
        let ma = 0.5 * (st1.mach + st2.mach);
        let xa = 0.5 * (st1.x + st2.x);
        let ta = 0.5 * (st1.theta + st2.theta);
        
        // Cf term with central + endpoint weighting for accuracy
        let cfm = compute_midpoint_cf(st1, st2, bl_type, params);
        let cfx = 0.50 * cfm * xa / ta 
                + 0.25 * (st1.cf * st1.x / st1.theta + st2.cf * st2.x / st2.theta);
        
        let btmp = ha + 2.0 - ma;
        eq.res[1] = -(tlog + btmp * ulog - xlog * 0.5 * cfx);
        
        // ... Jacobian entries
    }
    
    // Equation 3: Shape parameter (kinetic energy)
    {
        let ha = 0.5 * (st1.h + st2.h);
        let hsa = 0.5 * (st1.hs + st2.hs);
        let hca = 0.5 * (st1.hc + st2.hc);
        
        let xot1 = st1.x / st1.theta;
        let xot2 = st2.x / st2.theta;
        
        let dix = (1.0 - upw) * st1.di * xot1 + upw * st2.di * xot2;
        let cfx = (1.0 - upw) * st1.cf * xot1 + upw * st2.cf * xot2;
        
        let btmp = 2.0 * hca / hsa + 1.0 - ha;
        eq.res[2] = -(hlog + btmp * ulog + xlog * (0.5 * cfx - dix));
        
        // ... Jacobian entries
    }
    
    eq
}
```

---

## Module: `inverse.rs` - Inverse Mode for Separation

### `solve_inverse` - Inverse Mode System

```rust
/// Separation thresholds
pub const HLMAX: f64 = 3.8;  // Laminar
pub const HTMAX: f64 = 2.5;  // Turbulent

/// Solve BL equations in inverse mode (prescribed Hk, solve for Ue)
/// 
/// XFOIL: xbl.f MRCHUE lines 730-758
/// 
/// This is the KEY function for stall prediction.
/// When Hk exceeds threshold, we:
/// 1. Prescribe a target Hk (extrapolated from upstream)
/// 2. Add equation: Hk = Hk_target
/// 3. Solve for Ue as an unknown (instead of prescribing it)
/// 
/// # Arguments
/// * `st1` - Upstream station (known)
/// * `st2_init` - Initial guess for downstream station
/// * `htarg` - Target Hk value
/// * `is_laminar` - True if laminar, false if turbulent
/// * `is_wake` - True if in wake region
pub fn solve_inverse(
    st1: &BlStation,
    st2_init: &BlStation,
    htarg: f64,
    is_laminar: bool,
    is_wake: bool,
) -> Result<BlStation, SolverError> {
    let mut st2 = st2_init.clone();
    let mut ue = st2.ue;
    
    for _ in 0..25 {
        // Build 3 BL equations
        let bl_type = if is_wake {
            BlType::Wake
        } else if is_laminar {
            BlType::Laminar
        } else {
            BlType::Turbulent
        };
        
        let eq = bldif(st1, &st2, bl_type, &Default::default());
        
        // Add 4th equation: Hk = Hk_target (INVERSE MODE)
        // VS2(4,1) = 0
        // VS2(4,2) = HK2_T2
        // VS2(4,3) = HK2_D2
        // VS2(4,4) = HK2_U2
        // VSREZ(4) = HTARG - HK2
        
        let mut jac = [[0.0; 4]; 4];
        let mut res = [0.0; 4];
        
        // Copy BL equations (rows 0-2)
        for i in 0..3 {
            res[i] = eq.res[i];
            jac[i][0] = eq.jac2[i][0];  // d/dCtau
            jac[i][1] = eq.jac2[i][1];  // d/dTheta
            jac[i][2] = eq.jac2[i][2];  // d/dDstar
            jac[i][3] = eq.jac2[i][3];  // d/dUe
        }
        
        // Row 3: Hk constraint (INVERSE MODE EQUATION)
        jac[3][0] = 0.0;
        jac[3][1] = st2.hk_t;  // dHk/dTheta
        jac[3][2] = st2.hk_d;  // dHk/dDstar
        jac[3][3] = st2.hk_u;  // dHk/dUe
        res[3] = htarg - st2.hk;
        
        // Solve 4x4 system
        let delta = solve_4x4(&jac, &res)?;
        
        // Relaxation
        let dmax = delta[1].abs() / st2.theta
            .max(delta[2].abs() / st2.dstar)
            .max(delta[3].abs() / ue);
        let rlx = if dmax > 0.3 { 0.3 / dmax } else { 1.0 };
        
        // Update - INCLUDING Ue (this is the key difference from direct mode)
        st2.ctau += rlx * delta[0];
        st2.theta += rlx * delta[1];
        st2.dstar += rlx * delta[2];
        ue += rlx * delta[3];  // ← UE IS UPDATED IN INVERSE MODE
        
        // Clamp
        if !is_laminar {
            st2.ctau = st2.ctau.clamp(1e-7, 0.30);
        }
        
        // Recompute secondary variables
        st2 = BlStation::from_primary(st2.x, st2.theta, st2.dstar, ue, st2.ctau);
        
        // Check convergence
        if dmax < 1e-5 {
            return Ok(st2);
        }
    }
    
    Err(SolverError::InverseModeNotConverged)
}

/// Compute target Hk for inverse mode
/// 
/// XFOIL: xbl.f lines 692-721
pub fn compute_htarg(
    st1: &BlStation,
    dx: f64,
    is_laminar: bool,
    is_transition: bool,
    is_wake: bool,
    xt: f64,  // Transition location (if is_transition)
) -> f64 {
    let htarg = if is_laminar {
        // Laminar: relatively slow Hk increase downstream
        st1.hk + 0.03 * dx / st1.theta
    } else if is_transition {
        // Transition interval: weighted laminar and turbulent
        let dx_lam = xt - st1.x;
        let dx_turb = dx - dx_lam;
        st1.hk + (0.03 * dx_lam - 0.15 * dx_turb) / st1.theta
    } else if is_wake {
        // Wake: asymptotic behavior with backward Euler
        let const_val = 0.03 * dx / st1.theta;
        let mut hk2 = st1.hk;
        for _ in 0..3 {
            hk2 = hk2 - (hk2 + const_val * (hk2 - 1.0).powi(3) - st1.hk)
                      / (1.0 + 3.0 * const_val * (hk2 - 1.0).powi(2));
        }
        hk2
    } else {
        // Turbulent: relatively fast Hk decrease downstream
        st1.hk - 0.15 * dx / st1.theta
    };
    
    // Limit to reasonable values
    let hmax = if is_laminar { HLMAX } else { HTMAX };
    if is_wake {
        htarg.max(1.01)
    } else {
        htarg.max(hmax)
    }
}
```

---

## Module: `solver.rs` - Newton Solver

### `blsolv` - Block Elimination Solver (BLSOLV)

```rust
/// Custom block-elimination solver for the coupled BL system
/// 
/// XFOIL: xsolve.f:283
/// 
/// System structure:
/// ```text
/// [A  |  |  .  |  |  .  |][d1]   [R1]
/// [B  A  |  .  |  |  .  |][d2]   [R2]
/// [|  B  A  .  |  |  .  |][d3] = [R3]
/// [.  .  .  .  |  |  .  |][. ]   [. ]
/// [|  Z  |  |  B  A  .  |][dn]   [Rn]
/// ```
/// 
/// A, B = 3×2 BL equation Jacobians
/// | = 3×1 mass defect influence vectors
/// Z = 3×2 TE coupling block
pub struct BlNewtonSystem {
    /// Diagonal blocks VA(3,2,NSYS)
    pub va: Vec<[[f64; 2]; 3]>,
    /// Off-diagonal blocks VB(3,2,NSYS)
    pub vb: Vec<[[f64; 2]; 3]>,
    /// Mass influence VM(3,NSYS,NSYS)
    pub vm: Vec<Vec<[f64; 3]>>,
    /// TE coupling VZ(3,2)
    pub vz: [[f64; 2]; 3],
    /// Residuals VDEL(3,2,NSYS) - col 0 = residual, col 1 = sensitivity
    pub vdel: Vec<[[f64; 2]; 3]>,
    /// System size
    pub nsys: usize,
    /// TE station indices
    pub ivte1: usize,
    pub ivte2: usize,
}

impl BlNewtonSystem {
    pub fn solve(&mut self, vaccel: f64) -> Vec<[f64; 3]> {
        let nsys = self.nsys;
        
        // Forward elimination
        for iv in 0..nsys {
            let ivp = iv + 1;
            
            // Invert VA(iv) block (3x3 with mass column)
            
            // Normalize first row
            let pivot = 1.0 / self.va[iv][0][0];
            self.va[iv][0][1] *= pivot;
            for l in iv..nsys {
                self.vm[iv][l][0] *= pivot;
            }
            self.vdel[iv][0][0] *= pivot;
            self.vdel[iv][0][1] *= pivot;
            
            // Eliminate lower first column
            for k in 1..3 {
                let vtmp = self.va[iv][k][0];
                self.va[iv][k][1] -= vtmp * self.va[iv][0][1];
                for l in iv..nsys {
                    self.vm[iv][l][k] -= vtmp * self.vm[iv][l][0];
                }
                self.vdel[iv][k][0] -= vtmp * self.vdel[iv][0][0];
                self.vdel[iv][k][1] -= vtmp * self.vdel[iv][0][1];
            }
            
            // Continue with rows 2, 3...
            // (similar pattern for full 3x3 inversion)
            
            // Eliminate VB(iv+1) block
            if ivp < nsys {
                for k in 0..3 {
                    let vtmp1 = self.vb[ivp][k][0];
                    let vtmp2 = self.vb[ivp][k][1];
                    let vtmp3 = self.vm[ivp][iv][k];
                    
                    for l in ivp..nsys {
                        self.vm[ivp][l][k] -= vtmp1 * self.vm[iv][l][0]
                                            + vtmp2 * self.vm[iv][l][1]
                                            + vtmp3 * self.vm[iv][l][2];
                    }
                    // Similar for vdel...
                }
            }
            
            // Handle VZ at TE station
            if iv == self.ivte1 {
                let ivz = self.ivte2 + 1;
                // Eliminate VZ coupling...
            }
            
            // Eliminate lower VM columns (with acceleration)
            for kv in iv+2..nsys {
                for row in 0..3 {
                    let vtmp = self.vm[kv][iv][row];
                    if vtmp.abs() > vaccel {
                        for l in ivp..nsys {
                            self.vm[kv][l][row] -= vtmp * self.vm[iv][l][2];
                        }
                        self.vdel[kv][row][0] -= vtmp * self.vdel[iv][2][0];
                        self.vdel[kv][row][1] -= vtmp * self.vdel[iv][2][1];
                    }
                }
            }
        }
        
        // Backward substitution
        let mut solution = vec![[0.0; 3]; nsys];
        for iv in (1..nsys).rev() {
            let vtmp = self.vdel[iv][2][0];
            for kv in 0..iv {
                self.vdel[kv][0][0] -= self.vm[kv][iv][0] * vtmp;
                self.vdel[kv][1][0] -= self.vm[kv][iv][1] * vtmp;
                self.vdel[kv][2][0] -= self.vm[kv][iv][2] * vtmp;
            }
        }
        
        // Extract solution
        for iv in 0..nsys {
            solution[iv] = [
                self.vdel[iv][0][0],
                self.vdel[iv][1][0],
                self.vdel[iv][2][0],
            ];
        }
        
        solution
    }
}
```

---

## Module: `viscal.rs` - Main Solver Loop

```rust
/// Main viscous solution convergence loop
/// 
/// XFOIL: xoper.f:2886 VISCAL
pub fn viscal(
    state: &mut FlexFoilState,
    max_iter: usize,
) -> Result<ViscousSolution, SolverError> {
    const EPS: f64 = 1.0e-4;
    
    // 1. Initialize if needed
    if !state.wake_initialized {
        state.compute_wake()?;
    }
    if !state.bl_initialized {
        state.init_bl_from_inviscid()?;
    }
    
    // 2. Set up influence matrix
    if !state.dij_valid {
        state.compute_dij()?;
    }
    
    // 3. Newton iteration
    for iter in 0..max_iter {
        // Build Newton system (SETBL)
        let mut system = state.build_bl_system()?;
        
        // Solve (BLSOLV)
        let deltas = system.solve(state.params.vaccel);
        
        // Update with relaxation (UPDATE)
        let (rlx, rms, rmx) = state.apply_update(&deltas)?;
        
        // Update flow state
        if state.alpha_specified {
            state.update_mach_re_from_cl()?;
        } else {
            state.update_inviscid_from_alpha()?;
        }
        
        // Compute viscous velocities
        state.compute_qvis_from_uedg()?;
        state.compute_gam_from_qvis()?;
        
        // Move stagnation point
        state.move_stagnation_point()?;
        
        // Compute forces
        state.compute_cl_cd()?;
        
        // Log progress
        log::info!(
            "Iter {}: rms={:.4e}, max={:.4e}, CL={:.4f}, CD={:.5f}",
            iter, rms, rmx, state.cl, state.cd
        );
        
        // Check convergence
        if rms < EPS {
            return Ok(ViscousSolution {
                converged: true,
                iterations: iter + 1,
                cl: state.cl,
                cd: state.cd,
                cm: state.cm,
                alpha: state.alpha,
                xtr: [state.xtr[0], state.xtr[1]],
            });
        }
    }
    
    Err(SolverError::NotConverged {
        iterations: max_iter,
        rms: state.last_rms,
    })
}
```

---

## Testing Strategy

### Unit Tests for Closures

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_hkin() {
        // Test case from XFOIL (M=0)
        let (hk, hk_h, hk_msq) = hkin(2.5, 0.0);
        assert_relative_eq!(hk, 2.5, epsilon = 1e-10);
        assert_relative_eq!(hk_h, 1.0, epsilon = 1e-10);
        assert_relative_eq!(hk_msq, -0.29, epsilon = 1e-10);
        
        // Test case with compressibility
        let (hk, _, _) = hkin(2.5, 0.25);
        assert_relative_eq!(hk, 2.144, epsilon = 0.001);
    }
    
    #[test]
    fn test_cfl_falkner_skan() {
        // Flat plate: Hk ≈ 2.59, Cf = 0.664/sqrt(Re_x)
        // At Re_θ = 1000: Cf ≈ 0.00066
        let (cf, _, _, _) = cfl(2.59, 1000.0, 0.0);
        assert_relative_eq!(cf, 0.00066, epsilon = 0.0001);
    }
    
    #[test]
    fn test_dampl_below_critical() {
        // Below critical Re_θ: no amplification
        let (ax, _, _, _) = dampl(2.5, 0.001, 100.0);
        assert_eq!(ax, 0.0);
    }
    
    #[test]
    fn test_dampl_above_critical() {
        // Above critical Re_θ: positive amplification
        let (ax, _, _, _) = dampl(2.5, 0.001, 1000.0);
        assert!(ax > 0.0);
    }
}
```

### Integration Tests

```rust
#[test]
fn test_naca0012_polar() {
    let mut solver = FlexFoilSolver::new();
    solver.load_airfoil("NACA0012");
    solver.set_reynolds(1e6);
    solver.set_ncrit(9.0);
    
    // Test several angles
    let test_cases = [
        (0.0, 0.0, 0.006),      // (alpha, Cl_expected, Cd_expected)
        (4.0, 0.44, 0.0065),
        (8.0, 0.88, 0.010),
        (12.0, 1.18, 0.020),
        (16.0, 1.05, 0.080),    // Post-stall: Cl should drop
    ];
    
    for (alpha, cl_exp, cd_exp) in test_cases {
        solver.set_alpha(alpha.to_radians());
        let result = solver.solve_viscous(100).unwrap();
        
        assert_relative_eq!(result.cl, cl_exp, epsilon = 0.05);
        assert_relative_eq!(result.cd, cd_exp, epsilon = 0.003);
    }
}

#[test]
fn test_stall_prediction() {
    // This is the critical test: ensure Cl_max exists
    let mut solver = FlexFoilSolver::new();
    solver.load_airfoil("NACA0012");
    solver.set_reynolds(1e6);
    
    let mut cl_values = Vec::new();
    for alpha_deg in 0..=20 {
        solver.set_alpha((alpha_deg as f64).to_radians());
        if let Ok(result) = solver.solve_viscous(100) {
            cl_values.push((alpha_deg, result.cl));
        }
    }
    
    // Find Cl_max
    let cl_max = cl_values.iter().map(|(_, cl)| *cl).fold(0.0, f64::max);
    let alpha_clmax = cl_values.iter().find(|(_, cl)| *cl == cl_max).unwrap().0;
    
    // Verify stall occurs (Cl decreases after max)
    let post_stall: Vec<_> = cl_values.iter()
        .filter(|(a, _)| *a > alpha_clmax)
        .collect();
    
    assert!(!post_stall.is_empty(), "Should have post-stall data");
    assert!(
        post_stall.iter().any(|(_, cl)| *cl < cl_max),
        "Cl should decrease after stall"
    );
    
    // Verify Cl_max is reasonable (around 1.1-1.3 for NACA0012)
    assert!(cl_max > 1.0 && cl_max < 1.5);
}
```

---

*This specification provides the foundation for implementing FlexFoil with exact XFOIL equivalence.*
