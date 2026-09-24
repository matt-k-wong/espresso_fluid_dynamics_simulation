"""Task B regression tests: input validation and geometric intent."""

import pytest

from espresso_m1.cli import parse_bar_list
from espresso_m1.driver import run_case, run_sweep
from espresso_m1.grid import make_grid
from espresso_m1.params import (
    PROVENANCE_LABELS,
    SimParams,
    apply_overrides,
    from_dict,
    load_catalog,
    to_dict,
)


@pytest.mark.parametrize("values", [
    {"Nr": 3.9},
    {"Nr": 3.5},
    {"Nz": 8.5},
    {"max_iter": 2.5},
    {"mu": float("inf")},
    {"mu": float("-inf")},
    {"mu": float("nan")},
    {"p_in": float("inf")},
    {"p_in": float("nan")},
    {"k0": float("nan")},
    {"max_iter": True},
    {"max_iter": False},
    {"Nr": True},
    {"mu": True},
    {"p_in": True},
    {"R": "0.029"},
])
def test_from_dict_rejects_bad_numerics(values):
    with pytest.raises((ValueError, TypeError)):
        from_dict(values)


def test_integral_floats_still_coerce():
    params = from_dict({"Nr": 40.0, "Nz": 80.0, "max_iter": 10.0})
    assert (params.Nr, params.Nz, params.max_iter) == (40, 80, 10)


def test_run_case_rejects_nonfinite_direct_construction():
    import dataclasses
    base = SimParams()
    with pytest.raises(ValueError):
        run_case(dataclasses.replace(base, mu=float("inf")))
    with pytest.raises(ValueError):
        run_case(dataclasses.replace(base, p_in=float("nan")))


def test_cli_integer_parsing_retained():
    params = apply_overrides(SimParams(), ["Nr=40", "max_iter=10"])
    assert (params.Nr, params.max_iter) == (40, 10)
    with pytest.raises(ValueError):
        apply_overrides(SimParams(), ["Nr=3.9"])


def test_unresolved_opening_rejected_with_resolution():
    with pytest.raises(ValueError, match="(?i)under-resolved|increase Nr|R_perf"):
        make_grid(8, 16, 0.029, 0.02, 0.029 * 0.007)
    with pytest.raises(ValueError):
        run_case(SimParams(Nr=8, Nz=16, R_perf=0.029 * 0.007))


def test_sealed_bed_is_valid_and_distinct():
    grid = make_grid(8, 16, 0.029, 0.02, 0.0)
    assert grid.open_area == 0.0
    case = run_case(SimParams(Nr=8, Nz=16, R_perf=0.0))
    assert case.converged


def test_catalog_constraints_are_honest():
    catalog = load_catalog()
    assert catalog["parameters"]["R"]["recommended_is_enforced"] is False
    for name in to_dict(SimParams()):
        entry = catalog["parameters"][name]
        assert entry["units"], name
        assert entry["description"], name
        assert entry["provenance"] in PROVENANCE_LABELS, name
        assert "admissible_range" in entry, name
        assert "recommended_range" in entry, name
        assert "provenance_detail" in entry, name
        default = to_dict(SimParams())[name]
        if isinstance(default, float):
            assert entry["default"] == pytest.approx(default), name
        else:
            assert entry["default"] == default, name
    # mu keeps a literature provenance with an actual source cited.
    mu = catalog["parameters"]["mu"]
    assert mu["provenance"] == "literature"
    assert "IAPWS" in mu["provenance_detail"]
    # Uncalibrated closure coefficients stay illustrative.
    for name in ("k0", "phi0", "sigma_c", "sigma_eff_top", "phi_min"):
        assert catalog["parameters"][name]["provenance"] == "illustrative assumption"


def test_range_parsing_never_overshoots():
    assert parse_bar_list("1:2:0.6") == pytest.approx([1e5, 1.6e5])
    assert parse_bar_list("3:1:-1") == pytest.approx([3e5, 2e5, 1e5])
    assert parse_bar_list("2:2:1") == pytest.approx([2e5])
    with pytest.raises(ValueError):
        parse_bar_list("3:1:1")
    with pytest.raises(ValueError):
        parse_bar_list("1:2:0")
    with pytest.raises(ValueError):
        parse_bar_list("1:2:nan")
    with pytest.raises(ValueError):
        parse_bar_list("")


def test_empty_and_negative_sweeps_rejected():
    with pytest.raises(ValueError):
        run_sweep(SimParams(), [])
    with pytest.raises(ValueError):
        run_sweep(SimParams(), [1e5, -1e5])
