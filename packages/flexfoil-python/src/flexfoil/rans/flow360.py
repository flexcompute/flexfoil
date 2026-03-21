"""Flow360 client wrapper for pseudo-2D airfoil RANS.

Uses the modern ``flow360`` SDK (v25+) for mesh upload and case submission.
Cases appear in the main Flow360 workspace under Project view.
Falls back to ``flow360client`` (v23) if the modern SDK is not installed.
"""

from __future__ import annotations

import json
import tempfile
import time
from pathlib import Path
from typing import Callable

from flexfoil.rans import RANSResult
from flexfoil.rans.config import build_case_config
from flexfoil.rans.mesh import generate_and_write_csm, generate_and_write_mesh


# ---------------------------------------------------------------------------
# SDK detection
# ---------------------------------------------------------------------------

def _has_modern_sdk() -> bool:
    """Check if the modern flow360 SDK (v25+) is available."""
    try:
        import flow360
        return hasattr(flow360, "VolumeMesh")
    except ImportError:
        return False


def _has_legacy_sdk() -> bool:
    """Check if the legacy flow360client SDK is available."""
    try:
        import flow360client
        return True
    except ImportError:
        return False


def check_auth() -> bool:
    """Verify that Flow360 credentials are configured."""
    try:
        import flow360client
        auth = flow360client.Config.auth
        return bool(auth and auth.get("accessToken"))
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Modern SDK (flow360 v25+)
# ---------------------------------------------------------------------------

def _submit_modern(
    ugrid_path: Path,
    case_config: dict,
    case_label: str,
    *,
    timeout: int = 3600,
    on_progress: Callable[[str, float], None] | None = None,
) -> tuple[str, str]:
    """Upload mesh and submit case using the modern flow360 SDK.

    Returns (case_id, mesh_id).
    """
    import flow360 as fl
    from flow360.component.v1.flow360_params import Flow360Params

    # Upload mesh
    if on_progress:
        on_progress("Uploading mesh", 0.1)

    draft = fl.VolumeMesh.from_file(
        str(ugrid_path),
        project_name=f"FlexFoil: {case_label}",
        solver_version="release-25.2",
    )
    vm = draft.submit()
    vm.wait()
    mesh_id = vm.id

    # Build Flow360Params from config dict
    if on_progress:
        on_progress("Submitting case", 0.15)

    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False)
    json.dump(case_config, tmp)
    tmp.close()
    params = Flow360Params.from_file(tmp.name)
    Path(tmp.name).unlink(missing_ok=True)

    # Submit case
    case_draft = fl.Case.create(
        name=case_label,
        params=params,
        volume_mesh_id=mesh_id,
    )
    case = case_draft.submit()
    case_id = case.id

    # Wait for completion
    if on_progress:
        on_progress("Running RANS solver", 0.2)

    # Poll until case completes (case.wait() can return prematurely)
    if on_progress:
        on_progress("Running RANS solver", 0.2)

    start = time.time()
    while True:
        info = case.get_info()
        status = info.caseStatus

        if on_progress:
            frac = {"preprocessing": 0.25, "queued": 0.3, "running": 0.5,
                     "postprocessing": 0.9, "completed": 1.0,
                     "error": 1.0, "diverged": 1.0}.get(status, 0.2)
            on_progress(status, frac)

        if status in ("completed", "error", "diverged"):
            break

        if time.time() - start > timeout:
            break

        time.sleep(10)

    return case_id, mesh_id


# ---------------------------------------------------------------------------
# Legacy SDK (flow360client v23)
# ---------------------------------------------------------------------------

def _submit_legacy(
    ugrid_path: Path,
    case_config: dict,
    case_label: str,
    *,
    timeout: int = 3600,
    on_progress: Callable[[str, float], None] | None = None,
) -> tuple[str, str]:
    """Upload mesh and submit case using the legacy flow360client SDK.

    Returns (case_id, mesh_id).
    """
    import flow360client

    if on_progress:
        on_progress("Uploading mesh", 0.1)

    mesh_id = flow360client.NewMesh(
        str(ugrid_path),
        meshName=case_label,
        meshJson={"boundaries": {"noSlipWalls": ["1"]}},
        fmat="aflr3",
        endianness="big",
    )

    if on_progress:
        on_progress("Submitting case", 0.15)

    # Legacy SDK uses integer boundary tags
    legacy_config = _config_with_integer_boundaries(case_config)
    case_id = flow360client.NewCase(
        meshId=mesh_id,
        config=legacy_config,
        caseName=case_label,
    )

    if on_progress:
        on_progress("Running RANS solver", 0.2)

    import flow360client.case as case_api
    start = time.time()
    while True:
        info = case_api.GetCaseInfo(case_id)
        status = info.get("status", "unknown")

        if on_progress:
            frac = {"preprocessing": 0.1, "queued": 0.15, "running": 0.5,
                     "postprocessing": 0.9, "completed": 1.0,
                     "error": 1.0, "diverged": 1.0}.get(status, 0.2)
            on_progress(status, frac)

        if status in ("completed", "error", "diverged"):
            break

        if time.time() - start > timeout:
            break

        time.sleep(10)

    return case_id, mesh_id


def _config_with_integer_boundaries(config: dict) -> dict:
    """Convert named boundaries to integer tags for legacy SDK."""
    config = dict(config)
    name_to_tag = {
        "wall": "1", "farfield": "2",
        "symmetry_y0": "3", "symmetry_y1": "4",
    }

    if "boundaries" in config:
        new_boundaries = {}
        for name, bc in config["boundaries"].items():
            tag = name_to_tag.get(name, name)
            new_boundaries[tag] = bc
        config["boundaries"] = new_boundaries

    if "surfaceOutput" in config and "surfaces" in config["surfaceOutput"]:
        new_surfaces = {}
        for name, so in config["surfaceOutput"]["surfaces"].items():
            tag = name_to_tag.get(name, name)
            new_surfaces[tag] = so
        config["surfaceOutput"]["surfaces"] = new_surfaces

    return config


# ---------------------------------------------------------------------------
# Result fetching (works with both SDKs via flow360client)
# ---------------------------------------------------------------------------

def fetch_results(case_id: str, *, alpha: float, Re: float, mach: float) -> RANSResult:
    """Download results from a completed Flow360 case.

    Uses the modern SDK to check status, falls back to legacy for force data.
    """
    # Check status via modern SDK if available
    status = "unknown"
    if _has_modern_sdk():
        try:
            import flow360 as fl
            case = fl.Case.from_cloud(case_id)
            info = case.get_info()
            status = info.caseStatus or str(info.status)
        except Exception:
            pass

    # Fall back to legacy SDK for status
    if status == "unknown" and _has_legacy_sdk():
        try:
            import flow360client.case as case_api
            info = case_api.GetCaseInfo(case_id)
            status = info.get("status", "unknown")
        except Exception:
            pass

    # Normalize status string (modern SDK may return enum or string)
    status = str(status).lower().replace("flow360status.", "")
    if "completed" not in status:
        return RANSResult(
            cl=0.0, cd=0.0, cm=0.0,
            alpha=alpha, reynolds=Re, mach=mach,
            converged=False, success=False,
            error=f"Case status: {status}",
            case_id=case_id,
        )

    # Fetch forces — try modern SDK first, then legacy
    forces = None

    if _has_modern_sdk():
        try:
            import flow360 as fl
            case = fl.Case.from_cloud(case_id)
            tf = case.results.total_forces
            df = tf.as_dataframe()
            forces = {col: df[col].tolist() for col in df.columns}
        except Exception:
            pass

    if forces is None and _has_legacy_sdk():
        try:
            import flow360client.case as case_api
            forces = case_api.GetCaseTotalForces(case_id)
        except Exception as e:
            return RANSResult(
                cl=0.0, cd=0.0, cm=0.0,
                alpha=alpha, reynolds=Re, mach=mach,
                converged=False, success=False,
                error=f"Failed to fetch forces: {e}",
                case_id=case_id,
            )

    if forces is None:
        return RANSResult(
            cl=0.0, cd=0.0, cm=0.0,
            alpha=alpha, reynolds=Re, mach=mach,
            converged=False, success=False,
            error="Could not fetch force data from either SDK",
            case_id=case_id,
        )

    cl = _get_last_value(forces, "CL")
    cd = _get_last_value(forces, "CD")
    cm = _get_last_value(forces, "CMz")
    cd_pressure = _get_last_value(forces, "CDPressure")
    cd_friction = _get_last_value(forces, "CDSkinFriction")

    return RANSResult(
        cl=cl, cd=cd, cm=cm,
        alpha=alpha, reynolds=Re, mach=mach,
        converged=True, success=True,
        cd_pressure=cd_pressure, cd_friction=cd_friction,
        case_id=case_id,
    )


def _get_last_value(forces: dict, key: str) -> float:
    """Extract the last (converged) value from a forces time-history dict."""
    if key in forces:
        val = forces[key]
        if isinstance(val, list):
            return float(val[-1]) if val else 0.0
        return float(val)
    return 0.0


# ---------------------------------------------------------------------------
# CSM-based meshing (uses Flow360's automated mesher)
# ---------------------------------------------------------------------------

def _submit_csm(
    csm_path: Path,
    case_config: dict,
    case_label: str,
    *,
    Re: float = 1e6,
    timeout: int = 3600,
    on_progress: Callable[[str, float], None] | None = None,
) -> tuple[str, str]:
    """Upload CSM geometry, auto-mesh, and submit case.

    Uses Flow360's automated meshing pipeline:
    CSM geometry → surface mesh → volume mesh → RANS case.

    Returns (case_id, mesh_id).
    """
    from flexfoil.rans.mesh import estimate_first_cell_height

    import flow360client

    # Step 1: Surface mesh from geometry
    if on_progress:
        on_progress("Generating surface mesh", 0.1)

    surface_config = {
        "maxEdgeLength": 0.05,
        "curvatureResolutionAngle": 15.0,
        "growthRate": 1.2,
    }

    import json as _json
    import tempfile as _tempfile

    surf_cfg_path = Path(_tempfile.mktemp(suffix=".json"))
    surf_cfg_path.write_text(_json.dumps(surface_config))

    surface_mesh_id = flow360client.NewSurfaceMeshFromGeometry(
        str(csm_path),
        str(surf_cfg_path),
    )
    surf_cfg_path.unlink(missing_ok=True)

    # Step 2: Volume mesh from surface
    if on_progress:
        on_progress("Generating volume mesh", 0.2)

    first_cell = estimate_first_cell_height(Re)
    volume_config = {
        "firstLayerThickness": first_cell,
        "growthRate": 1.15,
        "volume": {
            "firstLayerThickness": first_cell,
            "growthRate": 1.15,
        },
    }

    vol_cfg_path = Path(_tempfile.mktemp(suffix=".json"))
    vol_cfg_path.write_text(_json.dumps(volume_config))

    mesh_id = flow360client.NewMeshFromSurface(
        surfaceMeshId=surface_mesh_id,
        config=str(vol_cfg_path),
        meshName=case_label,
    )
    vol_cfg_path.unlink(missing_ok=True)

    # Step 3: Submit case
    if on_progress:
        on_progress("Submitting case", 0.3)

    case_id = flow360client.NewCase(
        meshId=mesh_id,
        config=case_config,
        caseName=case_label,
    )

    # Step 4: Wait
    if on_progress:
        on_progress("Running RANS solver", 0.35)

    import flow360client.case as case_api
    start = time.time()
    while True:
        info = case_api.GetCaseInfo(case_id)
        status = info.get("status", "unknown")
        if on_progress:
            frac = {"preprocessing": 0.4, "queued": 0.45, "running": 0.6,
                     "postprocessing": 0.9, "completed": 1.0,
                     "error": 1.0, "diverged": 1.0}.get(status, 0.35)
            on_progress(status, frac)
        if status in ("completed", "error", "diverged"):
            break
        if time.time() - start > timeout:
            break
        time.sleep(10)

    return case_id, mesh_id


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_rans(
    coords: list[tuple[float, float]],
    *,
    alpha: float = 0.0,
    Re: float = 1e6,
    mach: float = 0.2,
    airfoil_name: str = "airfoil",
    n_normal: int = 64,
    growth_rate: float = 1.15,
    farfield_radius: float = 100.0,
    span: float = 0.01,
    max_steps: int = 5000,
    turbulence_model: str = "SpalartAllmaras",
    timeout: int = 3600,
    on_progress: Callable[[str, float], None] | None = None,
    cleanup: bool = True,
    use_auto_mesh: bool = True,
) -> RANSResult:
    """Full RANS pipeline: coords → mesh → upload → solve → results.

    Parameters
    ----------
    use_auto_mesh : bool
        If True (default), use Flow360's automated meshing via CSM geometry.
        This produces higher-quality meshes. If False, generate a C-grid
        mesh locally and upload as UGRID.
    """
    if not check_auth():
        return RANSResult(
            cl=0.0, cd=0.0, cm=0.0,
            alpha=alpha, reynolds=Re, mach=mach,
            converged=False, success=False,
            error="Flow360 credentials not configured.",
        )

    use_modern = _has_modern_sdk()

    if not use_modern and not _has_legacy_sdk():
        return RANSResult(
            cl=0.0, cd=0.0, cm=0.0,
            alpha=alpha, reynolds=Re, mach=mach,
            converged=False, success=False,
            error="No Flow360 SDK installed. Run: pip install flow360",
        )

    tmpdir = tempfile.mkdtemp(prefix="flexfoil_rans_")

    try:
        case_label = f"{airfoil_name}_a{alpha:.1f}_Re{Re:.0e}_M{mach:.2f}"
        case_config = build_case_config(
            alpha=alpha, Re=Re, mach=mach, span=span,
            max_steps=max_steps, turbulence_model=turbulence_model,
        )

        if use_auto_mesh and _has_legacy_sdk():
            # Primary path: CSM geometry → Flow360 automated meshing
            if on_progress:
                on_progress("Generating geometry", 0.05)

            csm_path = generate_and_write_csm(
                coords, tmpdir, span=span, mesh_name=airfoil_name,
            )

            case_id, mesh_id = _submit_csm(
                csm_path, case_config, case_label,
                Re=Re, timeout=timeout, on_progress=on_progress,
            )

        else:
            # Fallback: generate C-grid mesh locally
            if on_progress:
                on_progress("Generating mesh", 0.05)

            ugrid_path, mapbc_path = generate_and_write_mesh(
                coords, tmpdir,
                Re=Re, n_normal=n_normal,
                growth_rate=growth_rate, farfield_radius=farfield_radius,
                span=span, mesh_name=airfoil_name,
            )

            if use_modern:
                case_id, mesh_id = _submit_modern(
                    ugrid_path, case_config, case_label,
                    timeout=timeout, on_progress=on_progress,
                )
            else:
                case_id, mesh_id = _submit_legacy(
                    ugrid_path, case_config, case_label,
                    timeout=timeout, on_progress=on_progress,
                )

        # Step 5: Fetch results
        if on_progress:
            on_progress("Fetching results", 0.95)

        result = fetch_results(case_id, alpha=alpha, Re=Re, mach=mach)
        result.mesh_id = mesh_id
        return result

    except Exception as e:
        return RANSResult(
            cl=0.0, cd=0.0, cm=0.0,
            alpha=alpha, reynolds=Re, mach=mach,
            converged=False, success=False,
            error=str(e),
        )

    finally:
        if cleanup:
            import shutil
            shutil.rmtree(tmpdir, ignore_errors=True)
