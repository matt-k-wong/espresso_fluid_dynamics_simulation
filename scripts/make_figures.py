#!/usr/bin/env python3
"""Regenerate README figures from synthetic fixtures (deterministic).

Outputs docs/figures/*.svg via src/espresso_m1/shots.py plotting.
All figures are synthetic/illustrative and labeled as such; rerun with:
    PYTHONPATH=src ./.venv/bin/python scripts/make_figures.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from espresso_m1.shots import (compare_shots, constant_flow_baseline,
                               derive_flow, load_shot_file, write_svg_plot)

FIG = ROOT / "docs" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

LIN = str(ROOT / "examples/shots/synthetic_linear.json")
IRR = str(ROOT / "examples/shots/synthetic_irregular.json")
EY = str(ROOT / "examples/shots/synthetic_ey.json")

linear = load_shot_file(LIN)
irregular = load_shot_file(IRR)
ey_shot = load_shot_file(EY)

# 1. Mass vs time: two observed traces + simulated constant-flow baseline.
comp = compare_shots([linear], [irregular], 0.0)
q = comp["fit"]["q_g_s"]
mass_series = []
for s, tag in ((linear, "linear 1.2 g/s"), (irregular, "irregular 1.2 g/s")):
    times = [p.time_s for p in s.samples]
    masses = [p.beverage_mass_g for p in s.samples]
    mass_series.append((f"{tag} (observed)", times, masses, "observed"))
    mass_series.append(("constant-flow baseline (simulated)", times,
                        constant_flow_baseline(times, q, 0.0, 0.0), "simulated"))
write_svg_plot(str(FIG / "shot_mass.svg"), mass_series,
               xlabel="time_s (s)", ylabel="beverage_mass_g (g)",
               title="Synthetic shots: mass vs time (illustrative)")

# 2. Derived flow for both traces.
flow_series = []
for s, tag in ((linear, "linear"), (irregular, "irregular")):
    times = [p.time_s for p in s.samples]
    flows, _ = derive_flow(s.samples, 2.0)
    flow_series.append((f"{tag} derived flow (observed)", times, flows,
                        "observed"))
write_svg_plot(str(FIG / "shot_flow.svg"), flow_series,
               xlabel="time_s (s)", ylabel="derived_flow_g_s (g/s)",
               title="Synthetic shots: derived flow, 2 s window (illustrative)")

# 3. Held-out residuals of the irregular trace vs the trained baseline.
resid = comp["heldout_residuals"][0]
times = [p.time_s for p in irregular.samples]
write_svg_plot(str(FIG / "heldout_residuals.svg"),
               [("irregular residuals (observed - simulated)", times,
                 resid["residuals_g"], "observed")],
               xlabel="time_s (s)", ylabel="residual_g (g)",
               title="Held-out residuals vs constant-flow baseline (synthetic)")

print(f"q={q:.4f} g/s rmse={resid['rmse_g']:.4f} g "
      f"ey_demo={100.0 * 36.0 * 0.085 / 18.0:.1f} %")
print(f"wrote {sorted(p.name for p in FIG.glob('*.svg'))}")
