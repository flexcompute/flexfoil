"""Flow360 client wrapper for pseudo-2D airfoil RANS.

Handles mesh upload, case submission, polling, and result retrieval.
"""

from __future__ import annotations

import tempfile
import time
from pathlib import Path
from typing import Callable

from flexfoil.rans import RANSResult
from flexfoil.rans.config import build_case_config
from flexfoil.rans.mesh import generate_and_write_mesh


def _import_flow360():
    """Lazy-import flow360client, raising a helpful error if missing."""
    try:
        import flow360client
        import flow360client.case as case_api
        return flow360client, case_api
    except ImportError:
        raise ImportError(
            "flow360client is required for RANS analysis. "
            "Install it with: pip install flexfoil[rans]"
        ) from None


def check_auth() -> bool:
    """Verify that flow360client credentials are configured."""
    flow360client, _ = _import_flow360()
    try:
        auth = flow360client.Config.auth
        return bool(auth and auth.get("accessToken"))
    except Exception:
        return False


def submit_mesh(
    ugrid_path: str | Path,
    mapbc_path: str | Path,
    *,
    mesh_name: str = "flexfoil-airfoil",
) -> str:
    """Upload a UGRID mesh to Flow360.

    Returns the mesh ID.
    """
    flow360client, _ = _import_flow360()

    mesh_json = {
        "boundaries": {
            "noSlipWalls": ["1"],
        }
    }

    mesh_id = flow360client.NewMesh(
        str(ugrid_path),
        meshName=mesh_name,
        meshJson=mesh_json,
        fmat="aflr3",
        endianness="big",
    )

    return mesh_id


def submit_case(
    mesh_id: str,
    case_config: dict,
    *,
    case_name: str = "flexfoil-rans",
) -> str:
    """Submit a RANS case to Flow360.

    Returns the case ID.
    """
    flow360client, _ = _import_flow360()
    case_id = flow360client.NewCase(
        meshId=mesh_id,
        config=case_config,
        caseName=case_name,
    )
    return case_id


def wait_for_case(
    case_id: str,
    *,
    timeout: int = 3600,
    poll_interval: int = 10,
    on_progress: Callable[[str, float], None] | None = None,
) -> dict:
    """Wait for a Flow360 case to complete.

    Parameters
    ----------
    case_id : str
        Flow360 case ID.
    timeout : int
        Maximum wait time in seconds (default 1 hour).
    poll_interval : int
        Time between status checks in seconds.
    on_progress : callable or None
        Called with (status_string, progress_fraction) on each poll.

    Returns
    -------
    dict
        Case info from GetCaseInfo.
    """
    _, case_api = _import_flow360()

    start = time.time()
    while True:
        info = case_api.GetCaseInfo(case_id)
        status = info.get("status", "unknown")

        if on_progress:
            # Estimate progress from status
            progress_map = {
                "preprocessing": 0.1,
                "queued": 0.15,
                "running": 0.5,
                "postprocessing": 0.9,
                "completed": 1.0,
                "error": 1.0,
                "diverged": 1.0,
            }
            frac = progress_map.get(status, 0.2)
            on_progress(status, frac)

        if status in ("completed", "error", "diverged"):
            return info

        elapsed = time.time() - start
        if elapsed > timeout:
            return {"status": "timeout", "elapsed": elapsed}

        time.sleep(poll_interval)


def fetch_results(case_id: str, *, alpha: float, Re: float, mach: float) -> RANSResult:
    """Download results from a completed Flow360 case.

    Returns a RANSResult with integrated forces.
    """
    _, case_api = _import_flow360()

    info = case_api.GetCaseInfo(case_id)
    status = info.get("status", "unknown")

    if status != "completed":
        return RANSResult(
            cl=0.0, cd=0.0, cm=0.0,
            alpha=alpha, reynolds=Re, mach=mach,
            converged=False, success=False,
            error=f"Case {status}: {info.get('statusMessage', 'unknown error')}",
            case_id=case_id,
        )

    try:
        forces = case_api.GetCaseTotalForces(case_id)
    except Exception as e:
        return RANSResult(
            cl=0.0, cd=0.0, cm=0.0,
            alpha=alpha, reynolds=Re, mach=mach,
            converged=False, success=False,
            error=f"Failed to fetch forces: {e}",
            case_id=case_id,
        )

    # forces is a dict with time-history lists; take the last (converged) value
    cl = _get_last_value(forces, "CL")
    cd = _get_last_value(forces, "CD")
    cm = _get_last_value(forces, "CMz")
    cd_pressure = _get_last_value(forces, "CDPressure")
    cd_friction = _get_last_value(forces, "CDSkinFriction")

    return RANSResult(
        cl=cl,
        cd=cd,
        cm=cm,
        alpha=alpha,
        reynolds=Re,
        mach=mach,
        converged=True,
        success=True,
        cd_pressure=cd_pressure,
        cd_friction=cd_friction,
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
) -> RANSResult:
    """Full RANS pipeline: coords → mesh → upload → solve → results.

    Parameters
    ----------
    coords : list of (x, y)
        Airfoil coordinates (Selig ordering).
    alpha : float
        Angle of attack in degrees.
    Re : float
        Reynolds number.
    mach : float
        Freestream Mach number.
    airfoil_name : str
        Name for the mesh/case in Flow360.
    n_normal : int
        Mesh cells in wall-normal direction.
    growth_rate : float
        BL mesh growth rate.
    farfield_radius : float
        Farfield distance in chord lengths.
    span : float
        Pseudo-3D span.
    max_steps : int
        Max pseudo-time steps.
    turbulence_model : str
        'SpalartAllmaras' or 'kOmegaSST'.
    timeout : int
        Max wait time in seconds.
    on_progress : callable or None
        Progress callback: (status, fraction).
    cleanup : bool
        Remove temporary mesh files after upload.

    Returns
    -------
    RANSResult
    """
    if not check_auth():
        return RANSResult(
            cl=0.0, cd=0.0, cm=0.0,
            alpha=alpha, reynolds=Re, mach=mach,
            converged=False, success=False,
            error="Flow360 credentials not configured. Run flow360client.ChooseAccount() first.",
        )

    tmpdir = tempfile.mkdtemp(prefix="flexfoil_rans_")

    try:
        # Step 1: Generate mesh
        if on_progress:
            on_progress("Generating mesh", 0.05)

        ugrid_path, mapbc_path = generate_and_write_mesh(
            coords,
            tmpdir,
            Re=Re,
            n_normal=n_normal,
            growth_rate=growth_rate,
            farfield_radius=farfield_radius,
            span=span,
            mesh_name=airfoil_name,
        )

        # Step 2: Upload mesh
        if on_progress:
            on_progress("Uploading mesh", 0.1)

        case_label = f"{airfoil_name}_a{alpha:.1f}_Re{Re:.0e}_M{mach:.2f}"
        mesh_id = submit_mesh(ugrid_path, mapbc_path, mesh_name=case_label)

        # Step 3: Submit case
        if on_progress:
            on_progress("Submitting case", 0.15)

        case_config = build_case_config(
            alpha=alpha,
            Re=Re,
            mach=mach,
            span=span,
            max_steps=max_steps,
            turbulence_model=turbulence_model,
        )

        case_id = submit_case(mesh_id, case_config, case_name=case_label)

        # Step 4: Wait for completion
        if on_progress:
            on_progress("Running RANS solver", 0.2)

        info = wait_for_case(
            case_id,
            timeout=timeout,
            on_progress=on_progress,
        )

        if info.get("status") == "timeout":
            return RANSResult(
                cl=0.0, cd=0.0, cm=0.0,
                alpha=alpha, reynolds=Re, mach=mach,
                converged=False, success=False,
                error=f"Case timed out after {timeout}s",
                case_id=case_id, mesh_id=mesh_id,
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
