"""Static experiment viewer: one self-contained HTML page, no backend.

``build_site`` reads one or more saved experiment bundles (see
``experiment.py``) and writes a single ``index.html`` with an experiment
dropdown, inline SVG plots, metrics tables, per-input provenance, and the
limitations/calibration text. No network, no live fitting, no machine I/O:
open the file directly in a browser. All text is HTML-escaped; synthetic
fixtures stay labeled synthetic.
"""

from __future__ import annotations

import csv
import html as _html
import json
from pathlib import Path

from .shots import ShotError


def _load_bundle(bundle_dir: str) -> dict:
    root = Path(bundle_dir)
    try:
        manifest = json.loads((root / "manifest.json").read_text())
        comparison = json.loads((root / "comparison.json").read_text())
    except (FileNotFoundError, json.JSONDecodeError) as exc:
        raise ShotError(f"bundle unreadable: {exc}") from exc
    mass_svg = (root / "mass.svg").read_text() if (root / "mass.svg").exists() else ""
    flow_svg = (root / "flow.svg").read_text() if (root / "flow.svg").exists() else ""
    metrics: list[dict] = []
    mp = root / "export" / "metrics.csv"
    if mp.exists():
        with open(mp, newline="") as fh:
            metrics = list(csv.DictReader(fh))
    return {"dir": str(bundle_dir), "manifest": manifest,
            "comparison": comparison, "mass_svg": mass_svg,
            "flow_svg": flow_svg, "metrics": metrics}


def _fmt(x, nd=4) -> str:
    if x is None or x == "":
        return "—"
    try:
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return _html.escape(str(x))


def build_site(bundle_dirs: list[str], out_path: str, title: str = "Espresso experiments") -> dict:
    """Build a static viewer page from saved bundles. Returns summary."""
    if not bundle_dirs:
        raise ShotError("need >= 1 bundle directory")
    bundles = [_load_bundle(b) for b in bundle_dirs]
    esc_title = _html.escape(title)
    options, sections = [], []
    for i, b in enumerate(bundles):
        m, c = b["manifest"], b["comparison"]
        label = _html.escape(str(m.get("label") or f"experiment-{i}"))
        options.append(f'<option value="exp{i}">{label}</option>')
        prov_rows = "".join(
            f"<tr><td>{_html.escape(str(e.get('shot_id', '?')))}</td>"
            f"<td>{_html.escape(str(e.get('provenance', '?')))}</td>"
            f"<td><code>{_html.escape(str(e.get('sha256', ''))[:12])}…</code></td></tr>"
            for e in m.get("inputs", []))
        metric_rows = "".join(
            f"<tr><td>{_html.escape(str(r.get('shot_id', '?')))}</td>"
            f"<td>{_html.escape(str(r.get('provenance', '?')))}</td>"
            f"<td>{_fmt(r.get('final_mass_g'), 3)}</td>"
            f"<td>{_fmt(r.get('brew_ratio'))}</td>"
            f"<td>{_fmt(r.get('cup_EY_pct'), 3)}</td></tr>"
            for r in b["metrics"])
        resid_rows = "".join(
            f"<tr><td>{_html.escape(str(r.get('shot_id', '?')))}</td>"
            f"<td>{_fmt(r.get('rmse_g'))}</td>"
            f"<td>{_fmt(r.get('max_abs_resid_g'))}</td>"
            f"<td>{_fmt(r.get('mean_resid_g'))}</td></tr>"
            for r in c.get("heldout_residuals", []))
        tr = c.get("train_repeatability", {}).get("final_mass_g", {})
        fit = c.get("fit", {})
        hidden = "" if i == 0 else ' hidden="hidden"'
        sections.append(
            f"<section id=\"exp{i}\"{hidden}>"
            f"<h2>{label}</h2>"
            f"<p>Train: <code>{_html.escape(', '.join(m.get('train_ids', [])))}</code> "
            f"| Holdout (whole shots): <code>{_html.escape(', '.join(m.get('holdout_ids', [])))}</code> "
            f"| t0={_html.escape(str(m.get('already_wet_t0_s')))} s "
            f"| fitted q={_fmt(fit.get('q_g_s'))} g/s (empirical baseline, train only; "
            f"permeability NOT varied)</p>"
            f"<p>Train repeatability: n={_html.escape(str(c.get('train_repeatability', {}).get('n_shots')))} "
            f"final-mass mean={_fmt(tr.get('mean'), 3)} g sd={_fmt(tr.get('stdev'), 3)} g</p>"
            f"<h3>Mass vs time (observed solid, simulated dashed)</h3>{b['mass_svg']}"
            f"<h3>Derived flow, window {_html.escape(str(m.get('flow_window_s')))} s (derived, not measured)</h3>{b['flow_svg']}"
            f"<h3>Per-shot metrics</h3><table><tr><th>shot</th><th>provenance</th>"
            f"<th>final g</th><th>brew ratio</th><th>EY %</th></tr>{metric_rows}</table>"
            f"<h3>Held-out residuals vs constant-flow baseline</h3><table><tr><th>shot</th>"
            f"<th>RMSE g</th><th>max|res| g</th><th>mean g</th></tr>{resid_rows}</table>"
            f"<h3>Input provenance</h3><table><tr><th>shot</th><th>class</th><th>sha256</th></tr>"
            f"{prov_rows}</table>"
            f"<p><b>Limitations:</b> {_html.escape(str(c.get('limitations', '')))}</p>"
            f"<p><b>Calibration:</b> {_html.escape(str(c.get('calibration_status', '')))}</p>"
            f"</section>")
    page = ("<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"utf-8\">"
            f"<title>{esc_title}</title>"
            "<style>body{font-family:sans-serif;max-width:900px;margin:2em auto;padding:0 1em}"
            "table{border-collapse:collapse;margin:1em 0}th,td{border:1px solid #999;padding:4px 8px}"
            "svg{max-width:100%;height:auto;border:1px solid #ddd}</style></head><body>"
            f"<h1>{esc_title}</h1>"
            f"<label>Experiment: <select id=\"picker\" onchange=\""
            "document.querySelectorAll('section').forEach(s=>s.hidden=true);"
            "document.getElementById(this.value).hidden=false;\">"
            + "".join(options) + "</select></label>"
            + "".join(sections) +
            "<hr><p>Static bundle viewer: no backend, no live fitting, no machine connection. "
            "Regenerate with <code>experiment-site --bundles … --out index.html</code>.</p>"
            "</body></html>")
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page, encoding="utf-8")
    return {"out": str(out), "n_experiments": len(bundles)}
