"""Simulation parameters, validation, and machine-readable catalog access.

All physical quantities are SI. Gauge pressure (Pa) is used consistently:
p = 0 corresponds to atmospheric pressure at the outlet.
"""

from __future__ import annotations

import json
import math
import numbers
from dataclasses import asdict, dataclass, fields
from pathlib import Path

BAR_PA = 1.0e5  # 1 bar in Pa

CATALOG_FILENAME = "catalog.json"

PROVENANCE_LABELS = ("measured", "literature", "illustrative assumption", "fitted")

PERMEABILITY_MODES = ("uniform", "compaction")

INT_FIELDS = ("Nr", "Nz", "max_iter")
FLOAT_FIELDS = (
    "R", "L", "R_perf", "mu", "k0", "phi0", "p_in", "p_out",
    "sigma_eff_top", "sigma_c", "phi_min",
    "tol_p_rel", "tol_p_atol_Pa", "tol_q_rel", "tol_res_rel", "relaxation",
)


@dataclass
class SimParams:
    """Milestone-1 parameter set (SI units throughout)."""

    # Basket / bed geometry
    R: float = 0.029  # basket inner radius [m]
    L: float = 0.02  # bed thickness [m]
    Nr: int = 20  # radial cells [-]
    Nz: int = 40  # axial cells [-]
    R_perf: float = 0.029  # perforated radius of bottom plate [m]

    # Fluid
    mu: float = 3.0e-4  # liquid dynamic viscosity [Pa s]

    # Bed (reference state)
    k0: float = 5.0e-15  # reference permeability [m^2]
    phi0: float = 0.35  # reference (uncompacted) porosity [-]

    # Boundary pressures (gauge, Pa)
    p_in: float = 9.0e5  # inlet (top) pressure [Pa]
    p_out: float = 0.0  # outlet (perforation) pressure [Pa]

    # Permeability closure
    permeability_mode: str = "uniform"  # "uniform" | "compaction"

    # Compaction closure (used only when permeability_mode == "compaction")
    sigma_eff_top: float = 1.0e4  # effective axial stress at the top [Pa]
    sigma_c: float = 2.0e5  # characteristic compaction stress [Pa]
    phi_min: float = 0.20  # minimum (fully compacted) porosity [-]

    # Nonlinear (Picard) solver
    tol_p_rel: float = 1.0e-8  # pressure-update tolerance (relative, inf-norm)
    tol_p_atol_Pa: float = 1.0e-6  # pressure-update floor [Pa]
    tol_q_rel: float = 1.0e-8  # outlet-flux-update tolerance (relative)
    tol_res_rel: float = 1.0e-6  # mass-residual tolerance (relative)
    max_iter: int = 100  # maximum Picard iterations [-]
    relaxation: float = 1.0  # Picard relaxation factor on k, (0, 1] [-]

    # Linear solver (fixed choice for milestone 1, see assembly/solver notes)
    linear_solver: str = "spsolve"


def catalog_path() -> Path:
    return Path(__file__).with_name(CATALOG_FILENAME)


def load_catalog() -> dict:
    """Load the machine-readable parameter catalog (meaning, units, defaults,
    enforced admissible ranges, unenforced recommended ranges, provenance)."""
    with open(catalog_path(), "r", encoding="utf-8") as fh:
        return json.load(fh)


def to_dict(params: SimParams) -> dict:
    return asdict(params)


def _coerce_int(key: str, value) -> int:
    """Coerce grid/iteration counts; reject booleans and fractions first."""
    if isinstance(value, bool):
        raise TypeError(f"{key} must be an integer, got boolean {value!r}")
    if isinstance(value, numbers.Integral):
        return int(value)
    if isinstance(value, float) and value.is_integer():
        return int(value)
    raise ValueError(f"{key} must be an integer, got {value!r}")


def from_dict(data: dict) -> SimParams:
    """Build validated SimParams from a plain dict (e.g. parsed JSON).

    Raises KeyError on unknown keys, TypeError on booleans/non-numeric
    values where numbers are required, and ValueError on invalid values.
    """
    known = {f.name for f in fields(SimParams)}
    unknown = set(data) - known
    if unknown:
        raise KeyError(f"Unknown parameter(s): {sorted(unknown)}")
    kwargs = dict(data)
    for key in INT_FIELDS:
        if key in kwargs:
            kwargs[key] = _coerce_int(key, kwargs[key])
    for key in FLOAT_FIELDS:
        if key in kwargs and isinstance(kwargs[key], bool):
            raise TypeError(f"{key} must be a number, got boolean")
    params = SimParams(**kwargs)
    validate(params)
    return params


def from_json_file(path: str | Path) -> SimParams:
    with open(path, "r", encoding="utf-8") as fh:
        return from_dict(json.load(fh))


def apply_overrides(params: SimParams, items: list[str]) -> SimParams:
    """Apply CLI-style ``key=value`` overrides with type coercion.

    Example: ["Nr=40", "p_in=9e5", "permeability_mode=compaction"].
    Integer fields use int() parsing (fractional strings are rejected);
    non-finite float strings are rejected by validation.
    """
    data = to_dict(params)
    types = {f.name: f.type for f in fields(SimParams)}
    for item in items:
        if "=" not in item:
            raise ValueError(f"Override {item!r} is not of the form key=value")
        key, _, value = item.partition("=")
        key, value = key.strip(), value.strip()
        if key not in data:
            raise KeyError(f"Unknown parameter {key!r}")
        target = types[key]
        if target is int or target == "int":
            try:
                data[key] = int(value)
            except ValueError:
                raise ValueError(f"{key} must be an integer, got {value!r}")
        elif target is float or target == "float":
            data[key] = float(value)
        else:
            data[key] = value
    return from_dict(data)


def _require_finite(name: str, value, errors: list[str]) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        errors.append(f"{name} must be a number, got {value!r}")
    elif not math.isfinite(value):
        errors.append(f"{name} must be finite, got {value!r}")


def validate(params: SimParams) -> None:
    """Check parameter validity; raise ValueError describing all violations."""
    errors: list[str] = []
    for name in FLOAT_FIELDS:
        _require_finite(name, getattr(params, name), errors)
    for name in INT_FIELDS:
        value = getattr(params, name)
        if isinstance(value, bool) or not isinstance(value, numbers.Integral):
            errors.append(f"{name} must be an integer, got {value!r}")
    if errors:
        raise ValueError("Invalid parameters: " + "; ".join(errors))
    if not params.R > 0:
        errors.append("R must be positive")
    if not params.L > 0:
        errors.append("L must be positive")
    if not params.Nr >= 2:
        errors.append("Nr must be an integer >= 2")
    if not params.Nz >= 2:
        errors.append("Nz must be an integer >= 2")
    if not 0.0 <= params.R_perf <= params.R:
        errors.append("R_perf must satisfy 0 <= R_perf <= R")
    if not params.mu > 0:
        errors.append("mu must be positive")
    if not params.k0 > 0:
        errors.append("k0 must be positive")
    if not 0.0 < params.phi0 < 1.0:
        errors.append("phi0 must be in (0, 1)")
    if not 0.0 < params.phi_min < params.phi0:
        errors.append("phi_min must satisfy 0 < phi_min < phi0")
    if params.permeability_mode not in PERMEABILITY_MODES:
        errors.append(f"permeability_mode must be one of {PERMEABILITY_MODES}")
    if not params.sigma_eff_top >= 0:
        errors.append("sigma_eff_top must be non-negative")
    if not params.sigma_c > 0:
        errors.append("sigma_c must be positive")
    if not params.p_in >= 0:
        errors.append("p_in (gauge) must be non-negative")
    if not params.p_out >= 0:
        errors.append("p_out (gauge) must be non-negative")
    if not params.tol_p_rel > 0:
        errors.append("tol_p_rel must be positive")
    if not params.tol_p_atol_Pa > 0:
        errors.append("tol_p_atol_Pa must be positive")
    if not params.tol_q_rel > 0:
        errors.append("tol_q_rel must be positive")
    if not params.tol_res_rel > 0:
        errors.append("tol_res_rel must be positive")
    if not params.max_iter >= 1:
        errors.append("max_iter must be an integer >= 1")
    if not 0.0 < params.relaxation <= 1.0:
        errors.append("relaxation must be in (0, 1]")
    if params.linear_solver != "spsolve":
        errors.append("linear_solver must be 'spsolve' in milestone 1")
    if errors:
        raise ValueError("Invalid parameters: " + "; ".join(errors))
