"""Parameter catalog consistency and validation (base coverage)."""

import pytest

from espresso_m1.params import (
    SimParams,
    apply_overrides,
    from_dict,
    to_dict,
)


def test_validation_rejects_bad_values():
    with pytest.raises(ValueError):
        from_dict({**to_dict(SimParams()), "R_perf": 1.0})  # > R
    with pytest.raises(ValueError):
        from_dict({**to_dict(SimParams()), "phi_min": 0.9})  # > phi0
    with pytest.raises(ValueError):
        from_dict({**to_dict(SimParams()), "permeability_mode": "darcy2"})
    with pytest.raises(KeyError):
        from_dict({**to_dict(SimParams()), "espresso_strength": 11})


def test_overrides_coerce_types():
    params = apply_overrides(SimParams(), ["Nr=40", "p_in=9e5",
                                           "permeability_mode=compaction"])
    assert params.Nr == 40 and isinstance(params.Nr, int)
    assert params.p_in == pytest.approx(9e5)
    assert params.permeability_mode == "compaction"
