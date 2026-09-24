"""Shared test helpers."""

from espresso_m1.params import SimParams


def uniform_params(**kwargs) -> SimParams:
    base = dict(permeability_mode="uniform", R_perf=0.029, p_in=9.0e5, p_out=0.0)
    base.update(kwargs)
    return SimParams(**base)


def compaction_params(**kwargs) -> SimParams:
    base = dict(permeability_mode="compaction", R_perf=0.029, p_in=9.0e5, p_out=0.0)
    base.update(kwargs)
    return SimParams(**base)
