"""Flow360 case construction: simulation.json (via the build's own flow360 API,
so it is version-matched) + the solver preprocessing chain that yields Flow360.json.

Wall patches use SlipWall on the two span (y) faces because their outward normals
are anti-parallel (a single SymmetryPlane group is rejected by MeshProcessor).
The reference area is set to chord*span so CL/CD are span-normalized.
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

from .config import CaseConfig

# Quiet the flow360 SDK: suppress the beta banner (checked at import) and the
# INFO chatter (unit system / default temperature / thermal state).
os.environ.setdefault("FLOW360_SUPPRESS_BETA_WARNING", "1")


def build_simulation_json(cfg: CaseConfig, wall_names: list[str],
                          boundary_names: list[str], out_path: str | Path) -> None:
    """Write simulation.json + inject the VolumeMesh asset-cache root item.
    ``wall_names``/``boundary_names`` are CGNS patch names (e.g. 'fluid/main')."""
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

    with fl.SI_unit_system:
        params = fl.SimulationParams(
            reference_geometry=fl.ReferenceGeometry(
                area=span * u.m ** 2, moment_length=1 * u.m, moment_center=(0, 0, 0) * u.m),
            operating_condition=fl.AerospaceCondition.from_mach_reynolds(
                mach=cfg.flow.mach, reynolds_mesh_unit=cfg.flow.reynolds,
                project_length_unit=1 * u.m,
                alpha=cfg.flow.alpha_deg * u.deg, temperature=cfg.flow.temperature * u.K),
            time_stepping=fl.Steady(max_steps=cfg.solver.max_steps, CFL=fl.AdaptiveCFL()),
            models=[
                fl.Fluid(navier_stokes_solver=fl.NavierStokesSolver(absolute_tolerance=1e-9),
                         turbulence_model_solver=fl.SpalartAllmaras(absolute_tolerance=1e-8)),
                fl.Wall(entities=walls),
                fl.Freestream(entities=[farfield]),
                fl.SlipWall(entities=[sym1, sym2]),
            ],
            outputs=[
                fl.SurfaceOutput(entities=walls, output_fields=["Cp", "Cf", "yPlus"]),
                fl.VolumeOutput(output_fields=["Mach", "primitiveVars"]),
                fl.SliceOutput(
                    entities=[fl.Slice(name="centerSpan", origin=(0, -span / 2, 0) * u.m,
                                       normal=(0, 1, 0))],
                    output_fields=["Cp", "Mach", "primitiveVars"]),
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


def _run(cmd: list[str], cwd: Path, env: dict) -> None:
    subprocess.run(cmd, cwd=str(cwd), env=env, check=True, capture_output=True, text=True)


def preprocess(workdir: str | Path, mesh_name: str, find, env: dict) -> str:
    """Run the version-matched preprocessing chain in ``workdir`` and return the
    path to the solver-ready Flow360.json. Assumes simulation.json + the CGNS mesh
    are already present. All of this runs fine in-session (only the solver does not).
    ``find`` resolves a tool name to its absolute path (see rans.env.make_env)."""
    workdir = Path(workdir)
    b = find
    _run([b("PrintMeshMetaData"), mesh_name], workdir, env)
    _run([b("convertSimulationToSolverJSON.py"),
          "--inputSimulationJson", "simulation.json",
          "--inputMeshMetadataJson", "meshMetaData.json",
          "--outputConvertedJson", "Flow360_processed.json",
          "--columnarDataProcessorJson", "columnar.json"], workdir, env)
    _run([b("generateMeshJson.py"),
          "--inputCaseJson", "Flow360_processed.json",
          "--inputMeshMetaDataJson", "meshMetaData.json",
          "--outputMeshJson", "Flow360Mesh.json"], workdir, env)
    _run([b("MeshPartitioner"), "--meshfile", mesh_name, "--partitions", "1"], workdir, env)
    _run([b("MeshProcessor"), "--threads", "1", mesh_name], workdir, env)
    _run([b("PrintMeshBoundaryBoundingBox"), mesh_name], workdir, env)
    _run([b("preprocessCaseJson.py"),
          "-i", "Flow360_processed.json",
          "--inputMeshProcessedJson", f"{mesh_name}.json",
          "-o", "Flow360.json", "--meshName", mesh_name,
          "--simulationBasedJson", "TRUE", "--preprocessAutoVis", "TRUE"], workdir, env)
    # tell the solver this is a case run
    fj = workdir / "Flow360.json"
    d = json.loads(fj.read_text())
    d.setdefault("runControl", {})["caseType"] = 0
    fj.write_text(json.dumps(d, indent=2))
    return str(fj)
