# Third-Party Notices

This file documents third-party software, assets, and commercial dependencies
included in or required by FlexFoil.

---

## Commercial Dependencies

### AG Grid Enterprise

FlexFoil's web UI (`flexfoil-ui`) uses [AG Grid Enterprise](https://www.ag-grid.com/)
for its data grid component. AG Grid Enterprise is proprietary software that
requires a separate commercial license from AG Grid Ltd.

- **License required:** AG Grid Enterprise License
- **Licensor:** AG Grid Ltd (<https://www.ag-grid.com/license-pricing/>)
- **Impact:** Without a valid license key (provided via the
  `VITE_AG_GRID_LICENSE_KEY` environment variable), the grid will display
  evaluation watermarks. All other FlexFoil functionality is unaffected.
- **Alternative:** The open-source AG Grid Community edition (MIT license) is
  also listed as a dependency and covers core grid functionality.

### TWK Everett Typeface

The FlexFoil web UI loads the TWK Everett typeface from Flexcompute's CDN. TWK
Everett is a commercial typeface designed by Teo Tuominen and published by
TypeWithKnife.

- **License required:** TWK Everett Commercial License
- **Licensor:** TypeWithKnife (<https://typewithknife.com>)
- **Usage restriction:** The font files referenced in this repository are
  served under Flexcompute's commercial license and **may not** be downloaded,
  redistributed, self-hosted, or used for any purpose — commercial or otherwise
  — without a separate license from TypeWithKnife. If you fork or deploy this
  project independently you **must** either purchase your own TWK Everett
  license or replace the font with an alternative.
- **Fallback:** The UI will fall back to system sans-serif fonts if the
  typeface fails to load.

---

## Open-Source Components

### mfoil

- **Path:** `mfoil/mfoil.py`
- **Author:** Krzysztof J. Fidkowski
- **License:** MIT
- **Source:** <https://github.com/kjfid/mfoil>
- **Description:** Python-based airfoil analysis code used as a reference
  implementation during development.

---

## Algorithmic Heritage

FlexFoil's solver (RustFoil) implements algorithms originally described in Mark
Drela's XFOIL. RustFoil is an independent implementation in Rust — it does not
contain any XFOIL source code. The mathematical methods are documented in:

- Drela, M. "XFOIL: An Analysis and Design System for Low Reynolds Number
  Airfoils" (1989)
- Drela, M. & Giles, M. B. "Viscous-Inviscid Analysis of Transonic and Low
  Reynolds Number Airfoils" (1987)
