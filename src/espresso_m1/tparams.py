"""Transport/exchange parameters, validation, and catalog access.

SI units throughout. Time zero means an already wetted bed; all values are
illustrative assumptions until calibrated (see catalog provenance labels and
REPORT_TRANSPORT.md). No measured shot data is used.
"""

from __future__ import annotations

import math
import numbers
from dataclasses import asdict, dataclass, fields

from .params import load_catalog

MODES = ("extraction", "tracer")

FLOAT_FIELDS = (
    "M_dose",
    "f_sol",
    "exchange_rate",
    "K",
    "D",
    "rho_solution",
    "t_end",
    "dt",
    "cfl_advective",
    "C_init",
    "C_inlet",
)


@dataclass
class TransportParams:
    """Milestone-2 parameter set (SI units throughout)."""

    mode: str = "extraction"  # "extraction" | "tracer"
    M_dose: float = 0.018  # dry coffee dose [kg]
    f_sol: float = 0.25  # accessible soluble mass fraction [-]
    exchange_rate: float = 0.1  # solid/liquid exchange rate lambda [1/s]
    K: float = 0.005  # equilibrium solid fraction / liquid conc [m^3/kg]
    D: float = 1e-9  # pore-liquid dispersion coefficient [m^2/s]
    rho_solution: float = 1000.0  # cup solution density [kg/m^3]
    t_end: float = 60.0  # final time [s]
    dt: float = 0.5  # requested time step [s]
    cfl_advective: float = 0.5  # accuracy limit factor; 0 disables [-]
    C_init: float = 0.0  # uniform initial liquid concentration [kg/m^3]
    C_inlet: float = 0.0  # inlet concentration for t < first breakpoint
    # Inlet schedule: [[t_seconds, C_kgm3], ...] sorted by t, piecewise const.
    C_in_schedule: list = None  # type: ignore[assignment]

    def __post_init__(self):
        if self.C_in_schedule is None:
            self.C_in_schedule = []


def to_dict(params: TransportParams) -> dict:
    return asdict(params)


def from_dict(data: dict) -> TransportParams:
    """Build validated TransportParams from a plain dict.

    Raises KeyError on unknown keys, TypeError on booleans/non-numeric
    values where numbers are required, and ValueError on invalid values.
    """
    known = {f.name for f in fields(TransportParams)}
    unknown = set(data) - known
    if unknown:
        raise KeyError(f"Unknown parameter(s): {sorted(unknown)}")
    kwargs = dict(data)
    for key in FLOAT_FIELDS:
        if key in kwargs and isinstance(kwargs[key], bool):
            raise TypeError(f"{key} must be a number, got boolean")
    schedule = kwargs.get("C_in_schedule", None)
    if schedule is not None:
        checked = []
        for item in schedule:
            try:
                t, c = item
            except (TypeError, ValueError):
                raise ValueError(
                    f"C_in_schedule entries must be [t, C] pairs, got {item!r}"
                )
            if isinstance(t, bool) or isinstance(c, bool):
                raise TypeError("C_in_schedule entries must be numbers")
            checked.append([float(t), float(c)])
        kwargs["C_in_schedule"] = checked
    params = TransportParams(**kwargs)
    validate(params)
    return params


def apply_overrides(params: TransportParams, items: list[str]) -> TransportParams:
    """Apply CLI-style ``key=value`` overrides (C_in_schedule is file-only)."""
    data = to_dict(params)
    types = {f.name: f.type for f in fields(TransportParams)}
    for item in items:
        if "=" not in item:
            raise ValueError(f"Override {item!r} is not of the form key=value")
        key, _, value = item.partition("=")
        key, value = key.strip(), value.strip()
        if key not in data:
            raise KeyError(f"Unknown parameter {key!r}")
        if key == "C_in_schedule":
            raise ValueError("C_in_schedule can only be set via the JSON file")
        target = types[key]
        if target is float or target == "float":
            data[key] = float(value)
        else:
            data[key] = value
    return from_dict(data)


def _require_finite(name: str, value, errors: list[str]) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        errors.append(f"{name} must be a number, got {value!r}")
    elif not math.isfinite(value):
        errors.append(f"{name} must be finite, got {value!r}")


def validate(params: TransportParams) -> None:
    """Check parameter validity; raise ValueError describing all violations."""
    errors: list[str] = []
    for name in FLOAT_FIELDS:
        _require_finite(name, getattr(params, name), errors)
    if errors:
        raise ValueError("Invalid transport parameters: " + "; ".join(errors))
    if params.mode not in MODES:
        errors.append(f"mode must be one of {MODES}")
    if not params.M_dose > 0:
        errors.append("M_dose must be positive")
    if not 0.0 <= params.f_sol <= 1.0:
        errors.append("f_sol must satisfy 0 <= f_sol <= 1")
    if not params.exchange_rate >= 0:
        errors.append("exchange_rate must be non-negative")
    if not params.K >= 0:
        errors.append("K must be non-negative")
    if not params.D >= 0:
        errors.append("D must be non-negative")
    if not params.rho_solution > 0:
        errors.append("rho_solution must be positive")
    if not params.t_end > 0:
        errors.append("t_end must be positive")
    if not params.dt > 0:
        errors.append("dt must be positive")
    if not params.cfl_advective >= 0:
        errors.append("cfl_advective must be non-negative")
    if not params.C_init >= 0:
        errors.append("C_init must be non-negative")
    if not params.C_inlet >= 0:
        errors.append("C_inlet must be non-negative")
    sched = params.C_in_schedule or []
    last_t = None
    for item in sched:
        t, c = item
        if not math.isfinite(t) or not math.isfinite(c):
            errors.append(f"C_in_schedule entry {item!r} must be finite")
        if t < 0 or c < 0:
            errors.append(f"C_in_schedule entry {item!r} must be non-negative")
        if last_t is not None and t < last_t:
            errors.append("C_in_schedule times must be non-decreasing")
        last_t = t
    if params.mode == "extraction":
        # Cup solute must be attributable to the coffee alone.
        if params.C_init != 0.0:
            errors.append("extraction mode requires C_init = 0")
        if params.C_inlet != 0.0 or any(c != 0.0 for _, c in sched):
            errors.append("extraction mode requires zero inlet solute")
    else:  # tracer: strictly passive, no solid pool, no exchange.
        if params.f_sol != 0.0:
            errors.append("tracer mode requires f_sol = 0")
        if params.exchange_rate != 0.0:
            errors.append("tracer mode requires exchange_rate = 0")
    if errors:
        raise ValueError("Invalid transport parameters: " + "; ".join(errors))


def inlet_concentration(params: TransportParams, t: float) -> float:
    """Piecewise-constant inlet history value at time t."""
    value = params.C_inlet
    for st, sc in params.C_in_schedule or []:
        if t >= st:
            value = sc
        else:
            break
    return float(value)


def schedule_breakpoints(params: TransportParams) -> list[float]:
    """Sorted positive schedule times (dt shortening targets)."""
    return sorted({float(st) for st, _ in (params.C_in_schedule or []) if st > 0.0})


def catalog_transport_entries() -> dict:
    """Transport entries of the shared machine-readable catalog."""
    return load_catalog()["parameters"]
