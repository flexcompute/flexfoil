"""Flow360 case construction: simulation.json (via the build's own flow360 API,
so it is version-matched) + the solver preprocessing chain that yields Flow360.json.

Wall patches use SlipWall on the two span (y) faces because their outward normals
are anti-parallel (a single SymmetryPlane group is rejected by MeshProcessor).
The reference area is set to chord*span so CL/CD are span-normalized.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
from pathlib import Path

from .config import CaseConfig

# Quiet the flow360 SDK: suppress the beta banner (checked at import) and the
# INFO chatter (unit system / default temperature / thermal state).
os.environ.setdefault("FLOW360_SUPPRESS_BETA_WARNING", "1")


def _write_params(cfg: CaseConfig, wall_names: list[str], boundary_names: list[str],
                  out_path: str | Path, *, make_time_stepping, udds=None,
                  output_frequency: int | None = None) -> None:
    """Build + write simulation.json (and inject the VolumeMesh asset-cache root item).

    Shared between the steady builder and the unsteady α-sweep builder. ``make_time_stepping``
    is ``(fl, u) -> TimeStepping``; ``udds`` an optional list of UserDefinedDynamic; and
    ``output_frequency`` (unsteady only) writes each output every N physical steps."""
    import flow360 as fl
    from flow360 import u
    from flow360.log import set_logging_level
    set_logging_level("ERROR")
    from flow360_schema.models.entities.surface_entities import Surface
    from flow360_schema.models.entities.volume_entities import GenericVolume
    from flow360_schema.models.entity_info import VolumeMeshEntityInfo

    walls = [Surface(name=n) for n in wall_names]
    farfield = Surface(name="fluid/farfield")
    sym1 = Surface(name="fluid/symmetry1")
    sym2 = Surface(name="fluid/symmetry2")
    span = cfg.mesh.span
    # `frequency` is rejected in steady runs, so only pass it when unsteady.
    freq = {} if output_frequency is None else {"frequency": output_frequency}

    with fl.SI_unit_system:
        params = fl.SimulationParams(
            reference_geometry=fl.ReferenceGeometry(
                area=span * u.m ** 2, moment_length=1 * u.m, moment_center=(0, 0, 0) * u.m),
            operating_condition=fl.AerospaceCondition.from_mach_reynolds(
                mach=cfg.flow.mach, reynolds_mesh_unit=cfg.flow.reynolds,
                project_length_unit=1 * u.m,
                alpha=cfg.flow.alpha_deg * u.deg, temperature=cfg.flow.temperature * u.K),
            time_stepping=make_time_stepping(fl, u),
            user_defined_dynamics=udds,
            models=[
                fl.Fluid(navier_stokes_solver=fl.NavierStokesSolver(absolute_tolerance=1e-9),
                         turbulence_model_solver=fl.SpalartAllmaras(absolute_tolerance=1e-8)),
                fl.Wall(entities=walls),
                fl.Freestream(entities=[farfield]),
                fl.SlipWall(entities=[sym1, sym2]),
            ],
            outputs=[
                fl.SurfaceOutput(entities=walls, output_fields=["Cp", "Cf", "yPlus"], **freq),
                fl.VolumeOutput(output_fields=["Mach", "primitiveVars"], **freq),
                fl.SliceOutput(
                    entities=[fl.Slice(name="centerSpan", origin=(0, -span / 2, 0) * u.m,
                                       normal=(0, 1, 0))],
                    output_fields=["Cp", "Mach", "primitiveVars"], **freq),
            ],
        )
    params.to_file(str(out_path))

    info = VolumeMeshEntityInfo(
        zones=[GenericVolume(name="fluid")],
        boundaries=[Surface(name=n, private_attribute_is_interface=False) for n in boundary_names])
    sim = json.loads(Path(out_path).read_text())
    ac = sim["private_attribute_asset_cache"]
    ac["project_entity_info"] = info.model_dump(by_alias=False, mode="json")
    ac["project_length_unit"] = {"value": 1.0, "units": "m"}
    Path(out_path).write_text(json.dumps(sim, indent=2))


def build_simulation_json(cfg: CaseConfig, wall_names: list[str],
                          boundary_names: list[str], out_path: str | Path) -> None:
    """Steady simulation.json. ``wall_names``/``boundary_names`` are CGNS patch names."""
    _write_params(cfg, wall_names, boundary_names, out_path,
                  make_time_stepping=lambda fl, u: fl.Steady(
                      max_steps=cfg.solver.max_steps, CFL=fl.AdaptiveCFL()))


def build_alpha_sweep_simulation_json(cfg: CaseConfig, wall_names: list[str],
                                      boundary_names: list[str], out_path: str | Path, *,
                                      alphas_deg: list[float], step_size_s: float = 1.0e6,
                                      max_pseudo_steps: int = 2000) -> None:
    """Unsteady-as-steady α sweep: one physical step per α, a huge time step (so the
    dual-time term → 0 ⇒ each step solves the steady equations via up to
    ``max_pseudo_steps`` pseudo-iterations), and a UserDefinedDynamic that drives the
    freestream ``alphaAngle`` as a function of the physical-step index. With
    ``output_frequency=1`` the solver writes forces + the center-span slice at every α,
    warm-started from the previous converged field. Requires evenly-spaced ``alphas_deg``."""
    import flow360 as fl

    a0 = alphas_deg[0]
    n = len(alphas_deg)
    da = (alphas_deg[-1] - a0) / (n - 1) if n > 1 else 0.0

    def make_ts(fl, u):
        # Each physical step is effectively a steady solve (huge Δt ⇒ no dual-time
        # damping), so use the conservative STEADY adaptive-CFL controller, NOT the
        # aggressive unsteady default (max 1e6 / maxRelChange 50 / convLimitFactor 1.0)
        # that Unsteady would otherwise fill in — that aggressive ramp, with no Δt
        # damping, blows up at stiff high-α steps. Steady defaults rarely diverge.
        return fl.Unsteady(steps=n, step_size=step_size_s * u.s,
                           max_pseudo_steps=max_pseudo_steps, order_of_accuracy=1,
                           CFL=fl.AdaptiveCFL.default_steady())

    # At the first pseudo-step of each physical step, set α = α0 + step·Δα; hold it
    # through the inner iterations. state[0] carries the current α. The solver's
    # `alphaAngle` control is in DEGREES (same convention as the freestream field),
    # so the constants are degrees — NOT radians (validated against steady runs).
    udd = fl.UserDefinedDynamic(
        name="alphaSweep",
        input_vars=["CL"],                      # required; unused here
        constants={"alpha0": a0, "dAlpha": da},
        output_vars={"alphaAngle": "state[0];"},
        state_vars_initial_value=["alpha0"],
        update_law=["if (pseudoStep == 0) alpha0 + physicalStep * dAlpha; else state[0];"],
    )
    _write_params(cfg, wall_names, boundary_names, out_path,
                  make_time_stepping=make_ts, udds=[udd], output_frequency=1)


def _run(cmd: list[str], cwd: Path, env: dict) -> None:
    subprocess.run(cmd, cwd=str(cwd), env=env, check=True, capture_output=True, text=True)


# The SDK case files that are mesh-independent (depend only on flow conditions +
# boundary names + solver settings) and therefore cacheable across geometry iterations.
_SDK_CACHE_FILES = ("simulation.json", "Flow360_processed.json", "columnar.json")


def preprocess(workdir: str | Path, mesh_name: str, find, env: dict, *,
               cfg, wall_names: list[str], boundary_names: list[str],
               timings: dict | None = None,
               sdk_cache_dir: str | Path | None = None,
               sim_builder=build_simulation_json) -> str:
    """Run the version-matched preprocessing chain in ``workdir`` and return the
    path to the solver-ready Flow360.json. Builds the (cacheable) SDK case JSONs
    itself. ``sdk_cache_dir`` (if set) caches/reuses the mesh-independent SDK outputs
    to skip the ~10 s flow360 SDK imports on a cache hit.
    (Note: preprocessCaseJson always populates auto-vis in this build, so the
    bounding-box step can't be skipped.)"""
    workdir = Path(workdir)
    b = find

    def step(name: str, cmd: list[str]) -> None:
        s = time.perf_counter()
        _run(cmd, workdir, env)
        if timings is not None:
            timings[name] = round(time.perf_counter() - s, 3)

    step("printMeshMetaData", [b("PrintMeshMetaData"), mesh_name])

    # SDK case JSONs (simulation.json + Flow360_processed.json + columnar.json) —
    # cacheable (mesh-independent). On a hit, copy them and skip the SDK imports.
    s = time.perf_counter()
    cache = Path(sdk_cache_dir) if sdk_cache_dir else None
    hit = cache is not None and all((cache / f).exists() for f in _SDK_CACHE_FILES)
    if hit:
        for f in _SDK_CACHE_FILES:
            shutil.copy(cache / f, workdir / f)
    else:
        sim_builder(cfg, wall_names, boundary_names, workdir / "simulation.json")
        _run([b("convertSimulationToSolverJSON.py"),
              "--inputSimulationJson", "simulation.json",
              "--inputMeshMetadataJson", "meshMetaData.json",
              "--outputConvertedJson", "Flow360_processed.json",
              "--columnarDataProcessorJson", "columnar.json"], workdir, env)
        if cache is not None:
            cache.mkdir(parents=True, exist_ok=True)
            for f in _SDK_CACHE_FILES:
                shutil.copy(workdir / f, cache / f)
    if timings is not None:
        timings["sdkCase"] = {"cache": "hit" if hit else "miss",
                              "t": round(time.perf_counter() - s, 3)}

    step("generateMeshJson", [b("generateMeshJson.py"),
          "--inputCaseJson", "Flow360_processed.json",
          "--inputMeshMetaDataJson", "meshMetaData.json",
          "--outputMeshJson", "Flow360Mesh.json"])
    step("meshPartitioner", [b("MeshPartitioner"), "--meshfile", mesh_name, "--partitions", "1"])
    step("meshProcessor", [b("MeshProcessor"), "--threads", "1", mesh_name])
    step("boundingBox", [b("PrintMeshBoundaryBoundingBox"), mesh_name])
    step("preprocessCaseJson", [b("preprocessCaseJson.py"),
          "-i", "Flow360_processed.json",
          "--inputMeshProcessedJson", f"{mesh_name}.json",
          "-o", "Flow360.json", "--meshName", mesh_name,
          "--simulationBasedJson", "TRUE", "--preprocessAutoVis", "TRUE"])
    # tell the solver this is a case run
    fj = workdir / "Flow360.json"
    d = json.loads(fj.read_text())
    d.setdefault("runControl", {})["caseType"] = 0
    fj.write_text(json.dumps(d, indent=2))
    return str(fj)
