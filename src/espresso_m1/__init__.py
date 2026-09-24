"""Milestone-1 package: 2D axisymmetric single-phase Darcy-flow solver."""

from .driver import analytical_cylindrical_Q, run_case, run_sweep
from .params import SimParams, from_dict, from_json_file, load_catalog
from .tparams import TransportParams
from .transport import run_transport

__all__ = [
    "SimParams",
    "TransportParams",
    "analytical_cylindrical_Q",
    "from_dict",
    "from_json_file",
    "load_catalog",
    "run_case",
    "run_sweep",
    "run_transport",
]

__version__ = "0.1.0"
