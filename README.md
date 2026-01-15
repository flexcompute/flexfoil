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

### Phase 1: Foundation (Current)
- [x] Geometry types (`Point`, `Panel`, `Body`)
- [x] Cubic spline interpolation
- [x] Cosine repaneling
- [x] WASM bridge skeleton
- [ ] Coordinate file I/O (Selig, Lednicer formats)

### Phase 2: Inviscid Solver
- [x] Influence coefficient matrix
- [x] Kutta condition
- [x] Linear system solve
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

## License

MIT
