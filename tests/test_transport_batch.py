"""Task E3: closed-batch exchange against the continuum reference.

Single sealed cell (zero flow, D=0): total M = B + W*C is invariant and

    C_eq = M/(W+K*m),
    C(t) = C_eq + (C(0)-C_eq)*exp(-lambda*(1+K*m/W)*t),
    B(t) = M - W*C(t).

Backward Euler must show asymptotic first-order time convergence (error
ratio near two for dt halving) and roundoff-exact inventory. K=0 must
approach B(t)=B0*exp(-lambda*t); lambda=0 must leave B unchanged.
"""

import numpy as np
import pytest

from espresso_m1.grid import make_grid
from espresso_m1.transport import (
    assemble_transport,
    transport_step,
    zero_flow_snapshot,
)

LAM = 0.1
K = 0.005


def batch_setup(K=K, lam=LAM):
    grid = make_grid(1, 1, 0.029, 0.02, 0.0)
    snap = zero_flow_snapshot(grid, phi_value=0.35)
    asm = assemble_transport(grid, snap.phi, 0.0, snap.Fr, snap.Fz)
    assert asm.L.nnz == 0
    return grid, snap, asm


def run_batch(dt, t_end, K=K, lam=LAM, C0=0.0, B0=None):
    grid, snap, asm = batch_setup(K, lam)
    W, m = float(snap.W[0, 0]), float(snap.m[0, 0])
    if B0 is None:
        B0 = 0.25 * m
    C = np.array([[C0]])
    B = np.array([[B0]])
    M = B0 + W * C0
    t = 0.0
    worst = 0.0
    while t < t_end - 1e-12:
        C, B, _ = transport_step(C, B, asm, snap.W, snap.m, K, lam, dt, 0.0)
        t += dt
        worst = max(worst, abs(float(B[0, 0] + W * C[0, 0]) - M) / M)
    C_eq = M / (W + K * m)
    rate = lam * (1.0 + K * m / W)
    C_ref = C_eq + (C0 - C_eq) * np.exp(-rate * t_end)
    return abs(float(C[0, 0]) - C_ref) / abs(C_ref), worst


def test_batch_first_order_time_convergence_and_exact_inventory():
    errs = []
    for dt in (0.5, 0.25, 0.125):
        err, worst = run_batch(dt, 20.0)
        errs.append(err)
        print(f"\ndt={dt}: rel err={err:.4e} inventory worst={worst:.2e}")
        assert worst < 1e-12
    assert errs[-1] < 1e-4
    for coarse, fine in zip(errs, errs[1:]):
        assert 1.5 < coarse / fine < 3.0


def test_batch_irreversible_limit_K_zero():
    grid, snap, asm = batch_setup(K=0.0, lam=LAM)
    m = float(snap.m[0, 0])
    B0 = 0.25 * m
    C = np.array([[0.0]])
    B = np.array([[B0]])
    C, B, _ = transport_step(C, B, asm, snap.W, snap.m, 0.0, LAM, 0.5, 0.0)
    # One backward-Euler step toward B0*exp(-lambda*dt), first-order close.
    assert float(B[0, 0]) / B0 == pytest.approx(np.exp(-LAM * 0.5), rel=0.1)
    errs = []
    # Full-horizon check with halving.
    for dt in (1.0, 0.5, 0.25):
        g, s, a = batch_setup(K=0.0, lam=LAM)
        Cc = np.array([[0.0]])
        Bb = np.array([[B0]])
        t = 0.0
        while t < 10.0 - 1e-12:
            Cc, Bb, _ = transport_step(Cc, Bb, a, s.W, s.m, 0.0, LAM, dt, 0.0)
            t += dt
        errs.append(abs(float(Bb[0, 0]) / B0 - np.exp(-LAM * 10.0)))
    assert errs[-1] < 0.02
    assert 1.5 < errs[0] / errs[1] < 2.6
    assert 1.5 < errs[1] / errs[2] < 2.6


def test_batch_lambda_zero_leaves_pool_unchanged():
    grid, snap, asm = batch_setup(K=K, lam=0.0)
    m = float(snap.m[0, 0])
    B0 = 0.25 * m
    C = np.array([[1.0]])
    B = np.array([[B0]])
    C1, B1, _ = transport_step(C, B, asm, snap.W, snap.m, K, 0.0, 1.0, 7.0)
    assert float(B1[0, 0]) == B0
    assert float(C1[0, 0]) == 1.0
