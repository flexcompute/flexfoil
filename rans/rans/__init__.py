"""Quasi-2D multi-element RANS pipeline for FlexFoil.

contours JSON -> in-house anisotropic mesh -> Flow360 case -> GPU solve -> CL/CD/Cp.
See README.md. Public entry point: ``rans.pipeline.run``.
"""
from .config import CaseConfig
from .pipeline import run

__all__ = ["CaseConfig", "run"]
