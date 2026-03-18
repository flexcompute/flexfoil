"""Type stubs for the native Rust extension module."""

def analyze_faithful(
    coords: list[float],
    alpha_deg: float,
    reynolds: float = 1e6,
    mach: float = 0.0,
    ncrit: float = 9.0,
    max_iterations: int = 100,
) -> dict[str, object]:
    """Viscous (XFOIL-faithful) analysis at a single operating point."""
    ...

def analyze_inviscid(
    coords: list[float],
    alpha_deg: float,
) -> dict[str, object]:
    """Inviscid panel-method analysis at a single angle of attack."""
    ...

def generate_naca4(
    designation: int,
    n_points_per_side: int | None = None,
) -> list[tuple[float, float]]:
    """Generate NACA 4-series airfoil using XFOIL's exact algorithm."""
    ...

def repanel_xfoil(
    coords: list[float],
    n_panels: int = 160,
    curv_param: float = 1.0,
    te_le_ratio: float = 0.15,
    te_spacing_ratio: float = 0.667,
) -> list[tuple[float, float]]:
    """Repanel airfoil using XFOIL's curvature-based algorithm."""
    ...

def deflect_flap(
    coords: list[float],
    hinge_x_frac: float,
    deflection_deg: float,
    hinge_y_frac: float = 0.5,
) -> list[tuple[float, float]]:
    """Deflect a flap by rotating points aft of the hinge."""
    ...

def analyze_faithful_batch(
    coords: list[float],
    alphas: list[float],
    reynolds: float = 1e6,
    mach: float = 0.0,
    ncrit: float = 9.0,
    max_iterations: int = 100,
) -> list[dict[str, object]]:
    """Batch viscous analysis: solve multiple alphas in parallel via rayon."""
    ...

def analyze_inviscid_batch(
    coords: list[float],
    alphas: list[float],
) -> list[dict[str, object]]:
    """Batch inviscid analysis: solve multiple alphas in parallel via rayon."""
    ...

def parse_dat_file(path: str) -> list[tuple[float, float]]:
    """Parse a Selig/Lednicer .dat file and return coordinate tuples."""
    ...
