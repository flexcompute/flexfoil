# FlexFoil

A 2D panel code with integral boundary layer based on Mark Drela's tools - modern, high-performance airfoil analysis engine written in Rust, compiling to WebAssembly for real-time (60 Hz) web-based interaction.

## Overview

FlexFoil (RustFoil core) is a ground-up rewrite of Mark Drela's XFOIL, architected from the start to support:

- **Multi-body configurations** (slats, main wing, flaps)
- **Real-time feedback** during geometry manipulation
- **WebAssembly output** for browser-based UIs
- **Modern Rust idioms** with strong type safety

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                         Applications                             │
├─────────────────────┬─────────────────────┬─────────────────────┤
│   rustfoil-wasm     │   rustfoil-cli      │   flexfoil-ui       │
│   WASM/JS bindings  │   Dev/testing CLI   │   React frontend    │
├─────────────────────┴─────────────────────┴─────────────────────┤
│                       rustfoil-solver                            │
│   • Inviscid panel method (Phase 2)                             │
│   • Boundary layer equations (Phase 3)                          │
│   • Viscous-inviscid interaction (Phase 4)                      │
├─────────────────────────────────────────────────────────────────┤
│                        rustfoil-core                             │
│   • Geometry types (Point, Panel, Body)                         │
│   • Cubic spline interpolation                                  │
│   • Coordinate I/O                                              │
└─────────────────────────────────────────────────────────────────┘
```

## Project Structure

```
flexfoil/
├── Cargo.toml              # Workspace root
├── crates/
│   ├── rustfoil-core/      # Geometry, splines, core types
│   ├── rustfoil-solver/    # Panel method, BL, VII
│   ├── rustfoil-wasm/      # WebAssembly bindings
│   └── rustfoil-cli/       # Command-line tool
├── flexfoil-ui/            # React + TypeScript frontend
└── testdata/               # Sample airfoil files
```

## Quick Start

### Build

```bash
# Build all crates
cargo build --release

# Run tests
cargo test

# Build WASM (requires wasm-pack)
cd crates/rustfoil-wasm
wasm-pack build --target web
```

### Frontend

```bash
cd flexfoil-ui
npm install
npm run dev
```

### CLI Usage

```bash
# Analyze an airfoil at 5° angle of attack
cargo run --bin rustfoil -- analyze testdata/naca0012.dat --alpha 5.0

# Generate a polar
cargo run --bin rustfoil -- polar testdata/naca0012.dat --alpha-start -5 --alpha-end 15

# Repanel with cosine spacing
cargo run --bin rustfoil -- repanel testdata/naca0012.dat --panels 100
```

### WASM Usage (JavaScript)

```javascript
import init, { RustFoil, analyze_airfoil } from 'rustfoil-wasm';

await init();

// Quick one-shot analysis
const coords = [1.0, 0.0, 0.5, 0.1, 0.0, 0.0, 0.5, -0.1, 1.0, 0.0];
const result = analyze_airfoil(coords, 5.0);
console.log(`Cl = ${result.cl}`);

// Interactive usage
const foil = new RustFoil();
foil.set_coordinates(coords);
foil.set_alpha(5.0);
const solution = foil.solve();
```

## Development Phases

### Phase 1: Foundation
- [x] Geometry types (`Point`, `Panel`, `Body`)
- [x] Cubic spline interpolation
- [x] XFOIL-exact NACA generator
- [x] XFOIL-exact paneling (PANGEN)
- [x] WASM bridge
- [ ] Coordinate file I/O (Selig, Lednicer formats)

### Phase 2: Inviscid Solver (Current)
- [x] Linear vorticity panel method (PSILIN)
- [x] Kutta condition
- [x] Blunt TE handling
- [x] XFOIL-exact Cl/Cm
- [ ] Wake modeling
- [ ] Multi-body interactions
- [ ] Real-time Cp display

### Phase 3: Boundary Layer
- [ ] Thwaites method (laminar)
- [ ] Head/Green method (turbulent)
- [ ] Transition prediction (eN)
- [ ] Separation detection

### Phase 4: Viscous-Inviscid Interaction
- [ ] Global Newton-Raphson
- [ ] Transpiration velocity
- [ ] Convergence acceleration

### Phase 5: Advanced Features
- [ ] Multi-body viscous
- [ ] Wake interaction
- [ ] Inverse design

## Technical Details

### Coordinate Convention
- **X-axis:** Downstream (freestream direction). LE at x≈0, TE at x≈1.
- **Y-axis:** Upward. Upper surface has positive y.
- **Panel ordering:** Counter-clockwise from TE lower to TE upper.

### Dependencies
- `nalgebra` - Linear algebra (SIMD-optimized)
- `wasm-bindgen` - JavaScript interop
- `serde` - Serialization for WASM

### Performance Targets
- Inviscid solve: < 16ms for 200 panels (60 Hz capable)
- WASM bundle: < 500 KB gzipped

## References

1. Drela, M. "XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils" (1989)
2. Katz, J. & Plotkin, A. "Low-Speed Aerodynamics" (2001)
3. Cebeci, T. & Bradshaw, P. "Momentum Transfer in Boundary Layers" (1977)
4. Fidkowski, K. J. "A Coupled Inviscid-Viscous Airfoil Analysis Solver, Revisited," *AIAA Journal*, Vol. 60, No. 5, 2022, pp. 2961–2971. [doi:10.2514/1.J061341](https://doi.org/10.2514/1.J061341)

## License

FlexFoil is licensed under the [MIT License](LICENSE).

Copyright (c) 2026 Flexcompute, Inc. and Harry Smith.

This repository includes third-party components and references commercial
dependencies. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for details
on AG Grid Enterprise, the TWK Everett typeface, and other attributions.

### Feature requests and usage analytics

The desktop and mobile **Feedback & requests** widget links to the public
[feature request tracker](https://github.com/flexcompute/flexfoil/issues?q=is%3Aissue+label%3Aenhancement).
Without `VITE_FEEDBACK_SHEET_URL`, the form opens a GitHub issue draft with the
chosen category and message using an explicit issue template. Labels come from
template metadata, so visitors do not need permission to assign labels. Users must sign in and submit on GitHub; a draft is
not a saved request. Feature requests use the existing `enhancement` label and
GitHub's open/closed status and discussion. A configured feedback service retains
the existing POST flow; its opaque response cannot confirm that a request was saved.

Local development emits no custom usage events by default. Set
`VITE_SOLVE_RUN_ANALYTICS=true` to enable them. The tracked `.env.production`
enables this flag for production builds, including deployed and packaged UIs;
set `VITE_SOLVE_RUN_ANALYTICS=false` to override it. Each visitor must also accept analytics. **Analytics preferences**
in Help (desktop) or the mobile footer can reopen the consent choice. Rejection
or unavailable consent storage suppresses custom events. The initial Google tag
configuration is deferred until acceptance, with advertising signals disabled.

| Event | Meaning | Coarse parameters |
| --- | --- | --- |
| `solve_run` | Solve requested, regardless of convergence | `solve_mode`, `solver_mode`, `n_panels` |
| `feature_use` | Panel selected, DAT imported/exported, SVG exported, documentation or feedback/tracker opened | `feature`, optional `panel_id` |
| `feedback_handoff` | GitHub draft opened; not a submitted issue | `feedback_type` |
| `feedback_sent` | Configured feedback POST completed; persistence unconfirmed | `feedback_type` |

Event page metadata uses the application base URL and referrer origin only.
Custom payloads do not include geometry, filenames, feedback text, contacts or
precise location. Solving, importing and exporting still work with analytics off.
Google Analytics provides approximate country/city dimensions from network
location; no browser geolocation permission is requested. City availability and
accuracy depend on property settings, network routing and privacy thresholds.

For reports, open the GA4 property containing measurement ID `G-065GK6XBSR`.
Register event-scoped custom dimensions `feature`, `panel_id`, `solve_mode`,
`solver_mode`, and `feedback_type`; register `n_panels` as a custom metric if needed.
Use an Exploration with event count and total users, broken down by feature/panel
or solve mode and Country/City. Filter for the events above and the Flexfoil
application base URL. Reports describe consenting visitors, not every user.
Property access and dimension registration are account-side setup; this code does
not configure the property or expose admin reports to public app users. Check
Enhanced Measurement settings so automatic events do not collect shared route
state or download names separately from these sanitized custom events.

The regression checks run locally without contacting GA:

```sh
(cd flexfoil-ui && npx vitest run --maxWorkers=1 && npm run build)
npx playwright test --config e2e/analytics/playwright.config.mts
```

The browser fixture mounts the real SolvePanel, feedback widget and consent UI
with an offline solver and in-memory run persistence. It verifies UI/event
contracts, not solver numerics or delivery to GA. Chrome must be installed.
