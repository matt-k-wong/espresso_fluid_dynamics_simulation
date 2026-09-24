"""Task E8: accuracy evidence, separated by axis.

dt refinement at fixed grid (explicit steps, advective limit off) and
spatial refinement at a controlled small dt, three resolutions each. Cup
solute mass/EY plus a volume-weighted field norm are reported; K and lambda
are identical across grids. The verdict states whether the examples are
adequately resolved; budget acceptance alone is not spatial accuracy.
"""

import numpy as np
import pytest

from espresso_m1.params import SimParams
from espresso_m1.tparams import TransportParams
from espresso_m1.transport import run_transport

T_END = 20.0


def extraction_run(Nr, dt, cfl=0.0):
    h = SimParams(permeability_mode="compaction", Nr=Nr, Nz=2 * Nr,
                  R_perf=0.0145)
    t = TransportParams(mode="extraction", t_end=T_END, dt=dt,
                        cfl_advective=cfl)
    case = run_transport(h, t)
    assert case.converged, case.error
    return case


def test_dt_refinement_first_order_cup_mass():
    masses = []
    for dt in (0.25, 0.125, 0.0625):
        case = extraction_run(20, dt)
        masses.append(case.result.M_cup)
        print(f"\ndt={dt}: M_cup={case.result.M_cup:.8e}")
    d1, d2 = abs(masses[1] - masses[0]), abs(masses[2] - masses[1])
    assert 1.5 < d1 / d2 < 2.6  # backward Euler is first-order in time


def restrict_to_coarse(fine_C, fine_volumes, factor=2):
    """Volume-average a fine field onto the once-coarser grid."""
    Nr2, Nz2 = fine_C.shape
    coarse = (
        fine_C[0::2, 0::2] * fine_volumes[0::2, 0::2]
        + fine_C[1::2, 0::2] * fine_volumes[1::2, 0::2]
        + fine_C[0::2, 1::2] * fine_volumes[0::2, 1::2]
        + fine_C[1::2, 1::2] * fine_volumes[1::2, 1::2]
    ) / (
        fine_volumes[0::2, 0::2] + fine_volumes[1::2, 0::2]
        + fine_volumes[0::2, 1::2] + fine_volumes[1::2, 1::2]
    )
    assert coarse.shape == (Nr2 // factor, Nz2 // factor)
    return coarse


def test_grid_refinement_cup_and_field():
    cases = [extraction_run(Nr, 0.0625) for Nr in (10, 20, 40)]
    masses = [c.result.M_cup for c in cases]
    for Nr, c in zip((10, 20, 40), cases):
        print(f"\n{Nr}x{2*Nr}: M_cup={c.result.M_cup:.8e} "
              f"EY={c.result.EY_pct:.5f} e_var={c.result.e_var:.4e}")
    d1, d2 = abs(masses[1] - masses[0]), abs(masses[2] - masses[1])
    assert d1 / d2 > 1.3  # first-order upwind spatial error shrinks
    rel_change = d2 / masses[2]
    print(f"fine relative cup-mass change: {rel_change:.3e}")
    # Field norm between successive grids (fine restricted to coarse).
    for coarse_c, fine_c in zip(cases, cases[1:]):
        coarse_C = coarse_c.result.C
        restricted = restrict_to_coarse(
            fine_c.result.C, fine_c.result.W, factor=2)
        vol = coarse_c.result.W / np.sum(coarse_c.result.W)
        norm = float(np.sqrt(np.sum(vol * (restricted - coarse_C) ** 2)))
        scale = float(np.sqrt(np.sum(vol * coarse_C ** 2)))
        print(f"field change/coarse norm: {norm / scale:.4e}")
    # Verdict: cup mass changes 1.5% on the last refinement (Richardson
    # remaining-error estimate ~1-2% at 40x80); the e_var field metric is
    # still grid-sensitive (see REPORT_TRANSPORT.md).
    assert rel_change < 0.03
