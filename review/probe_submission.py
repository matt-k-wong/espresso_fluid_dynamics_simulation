"""Independent review probes. Run from repository root with PYTHONPATH=src.

This is a review artifact, not part of the numerical engine. It writes evidence
to a caller-specified path and exits nonzero if any correctness check fails.
Keep submission_results.json as the original snapshot; use a new output path
after repairs. --wheel checks a wheel built from a temporary project copy.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import pathlib
import shutil
import subprocess
import sys
import tempfile
import zipfile
from dataclasses import replace

import numpy as np
import scipy
from scipy.integrate import quad

from espresso_m1.assembly import assemble
from espresso_m1.cli import main as cli_main, parse_bar_list
from espresso_m1.driver import run_case
from espresso_m1.grid import make_grid
from espresso_m1.params import SimParams, from_dict
from espresso_m1.solver import solve_linear


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--wheel", action="store_true")
    args = parser.parse_args()
    root = pathlib.Path(__file__).resolve().parents[1]
    if pathlib.Path.cwd() != root:
        raise SystemExit("Run from the repository root.")
    checks: list[dict] = []
    observations: dict = {}

    def check(name, ok, **details):
        checks.append({"name": name, "passed": bool(ok), **details})

    def cli(argv):
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                return cli_main(argv)
        except (ValueError, KeyError, SystemExit) as exc:
            return {"rejected": str(exc)}

    base = SimParams(permeability_mode="compaction", Nr=8, Nz=16, R_perf=0.0145)
    for relaxation in (1.0, 1e-8, 1e-12):
        c = run_case(replace(base, relaxation=relaxation))
        details = {"converged": c.converged, "iters": c.iters}
        ok = not c.converged  # honest nonconvergence is allowed
        if c.converged:
            A, b = assemble(c.params, c.grid, c.sol.k)
            raw_res = float(np.max(np.abs(A @ c.sol.p.T.ravel() - b)))
            details.update(
                Q_mLs=c.diag.Q_out_mLs,
                reported_res_rel=c.sol.res_rel,
                recomputed_Ap_minus_b_m3s=raw_res,
                mismatch_rel=c.diag.mismatch_rel,
                max_cell_residual_rel=c.diag.max_residual_rel,
            )
            ok = (c.diag.mismatch_rel < 1e-6 and c.diag.max_residual_rel < 1e-6)
        check(f"honest_convergence_relaxation_{relaxation:g}", ok, **details)

    c = run_case(replace(base, p_in=1e5, p_out=1e5))
    check("equal_pressure_equilibrium", c.converged and abs(c.diag.Q_out_m3s) < 1e-18,
          converged=c.converged, error=c.error)

    try:
        c = run_case(replace(base, R_perf=base.R * 0.007))
        check("reject_unresolved_nonzero_opening", not c.converged,
              converged=c.converged, ideal_open_area=c.grid.open_area_ideal,
              discrete_open_area=c.grid.open_area)
    except ValueError as exc:
        check("reject_unresolved_nonzero_opening", True, rejection=str(exc))

    for name, values in [
        ("fractional_grid_size", {"Nr": 3.9}),
        ("infinite_viscosity", {"mu": float("inf")}),
        ("infinite_pressure", {"p_in": float("inf")}),
        ("boolean_iteration_count", {"max_iter": True}),
    ]:
        try:
            from_dict(values)
            check(f"reject_{name}", False, result="accepted")
        except (ValueError, TypeError):
            check(f"reject_{name}", True)

    vals = parse_bar_list("1:2:0.6")
    check("pressure_range_does_not_overshoot", bool(vals) and max(vals) <= 2e5 * (1 + 1e-14),
          values_Pa=vals)

    with tempfile.TemporaryDirectory(prefix="espresso-review-") as td:
        temp = pathlib.Path(td)
        out = temp / "case"
        common = ["run", "--config", "configs/default.json", "--out-dir", str(out)]
        success = cli(common)
        failure = cli(common + ["--set", "max_iter=1"])
        stale_success = [p.name for p in out.iterdir() if p.name in ("summary.json", "fields.npz")]
        check("failure_removes_previous_success_artifacts", success == 0 and failure != 0 and not stale_success,
              files_after_failure=sorted(p.name for p in out.iterdir()))
        success_again = cli(common)
        check("success_removes_previous_failure_artifact", success_again == 0 and not (out / "failure.json").exists())
        empty = cli(["sweep", "--config", "configs/default.json", "--p-in-bar", "3:1:1",
                     "--out", str(temp / "empty.csv")])
        check("empty_sweep_is_not_success", empty != 0, result=empty)
        zero = temp / "zero"
        zero_exit = cli(["run", "--config", "configs/default.json", "--set", "p_in=0",
                         "--out-dir", str(zero)])
        try:
            def reject_constant(value):
                raise ValueError(value)
            json.loads((zero / "summary.json").read_text(), parse_constant=reject_constant)
            check("zero_flow_summary_is_strict_json", zero_exit == 0)
        except (ValueError, FileNotFoundError) as exc:
            check("zero_flow_summary_is_strict_json", False, error=str(exc))

    # Independent quadrature of the stated closure, without calling closures.py.
    p = SimParams(permeability_mode="compaction")
    f0 = p.phi0**3 / (1 - p.phi0)**2
    def normalized_k(t):
        phi = p.phi_min + (p.phi0-p.phi_min) * np.exp(
            -(p.sigma_eff_top+t*(p.p_in-p.p_out))/p.sigma_c)
        return phi**3 / (1-phi)**2 / f0
    ratio = quad(normalized_k, 0, 1, epsabs=1e-12, epsrel=1e-12)[0]
    q_exact = np.pi*p.R**2*p.k0*(p.p_in-p.p_out)/(p.mu*p.L)*ratio
    ref_rows = []
    for n in (10, 20, 40, 80):
        for radius in (p.R, p.R/2):
            params = replace(p, Nr=n, Nz=2*n, R_perf=radius)
            c = run_case(params)
            u = run_case(replace(params, permeability_mode="uniform"))
            if not (c.converged and u.converged):
                raise RuntimeError("Reference runs failed: " + c.error + u.error)
            ref_rows.append(dict(Nr=n, Nz=2*n, R_perf_over_R=radius/p.R,
                                 Q_mLs=c.diag.Q_out_mLs,
                                 relative_open_area_error=c.grid.open_area/c.grid.open_area_ideal-1,
                                 normalized_flow_ratio_error=c.diag.Q_out_m3s/u.diag.Q_out_m3s/ratio-1))
    observations["nonlinear_reference"] = dict(Q_exact_full_mLs=q_exact*1e6,
                                               compaction_to_uniform_ratio=ratio,
                                               grid_runs=ref_rows)
    full_errors = [abs(row["normalized_flow_ratio_error"]) for row in ref_rows if row["R_perf_over_R"] == 1]
    check("nonlinear_full_open_matches_independent_integral", full_errors[-1] < 3e-5 and
          all(a/b > 3.5 for a,b in zip(full_errors, full_errors[1:])), errors=full_errors)
    perf_errors = [abs(row["normalized_flow_ratio_error"]) for row in ref_rows if row["R_perf_over_R"] == .5]
    check("perforated_solution_approaches_transform_reference", perf_errors[-1] < 2e-4 and
          all(a > b for a,b in zip(perf_errors, perf_errors[1:])), errors=perf_errors)

    # Nontrivial radial and axial manufactured solution. Source integrates the
    # continuum Laplacian exactly over each cylindrical finite volume; it does
    # not obtain a 'reference' by applying the discrete matrix to sampled p.
    mms = []
    for n in (8, 16, 32, 64):
        par = replace(p, Nr=n, Nz=2*n, permeability_mode="uniform")
        g = make_grid(par.Nr, par.Nz, par.R, par.L, par.R_perf)
        a = np.pi/par.L
        amplitude = par.p_in/5
        rf = g.r_faces
        zf = np.linspace(0, par.L, par.Nz+1)
        int_f = np.pi*np.diff(rf**2-rf**4/par.R**2+rf**6/(3*par.R**4))
        int_lap_f = 2*np.pi*np.diff(-4*rf**2/par.R**2+4*rf**4/par.R**4)
        int_sin = (np.cos(a*zf[:-1])-np.cos(a*zf[1:]))/a
        source = -par.k0/par.mu*amplitude*(int_lap_f-a*a*int_f)[:,None]*int_sin[None,:]
        A,b = assemble(par,g,np.full((par.Nr,par.Nz),par.k0))
        numeric = solve_linear(A,b+source.T.ravel()).reshape(par.Nz,par.Nr).T
        exact = par.p_in*(1-g.z_centers[None,:]/par.L) + amplitude*(1-(g.r_centers[:,None]/par.R)**2)**2*np.sin(a*g.z_centers[None,:])
        error = float(np.sqrt(np.sum(g.volumes*(numeric-exact)**2)/np.sum(g.volumes))/amplitude)
        mms.append(dict(Nr=n,Nz=2*n,volume_weighted_L2_over_amplitude=error))
    errors = [row["volume_weighted_L2_over_amplitude"] for row in mms]
    check("manufactured_2d_second_order_accuracy", all(a/b > 3.8 for a,b in zip(errors,errors[1:])), grids=mms)

    if args.wheel:
        with tempfile.TemporaryDirectory(prefix="espresso-wheel-review-") as td:
            temp = pathlib.Path(td)
            project = temp / "project"
            project.mkdir()
            shutil.copy("pyproject.toml", project / "pyproject.toml")
            shutil.copytree("src", project / "src", ignore=shutil.ignore_patterns("__pycache__", "*.egg-info"))
            build = subprocess.run([sys.executable,"-m","pip","wheel",str(project),"--no-deps","-w",str(temp/"dist")],capture_output=True,text=True)
            if build.returncode:
                check("installed_package_contains_catalog", False, build_error=build.stderr[-2000:])
            else:
                wheel = next((temp/"dist").glob("*.whl"))
                with zipfile.ZipFile(wheel) as z:
                    names = z.namelist()
                check("installed_package_contains_catalog", "espresso_m1/catalog.json" in names, wheel_members=names)

    files = [path for folder in ("src", "tests", "configs") for path in (root/folder).rglob("*")
             if path.is_file() and "__pycache__" not in path.parts and path.suffix in (".py", ".json")]
    files += [root/name for name in ("README.md", "ARCHITECTURE.md", "REPORT.md", "pyproject.toml")]
    evidence = dict(environment=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
                    checks=checks, observations=observations,
                    submission_sha256={str(path.relative_to(root)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(files)})
    outpath = pathlib.Path(args.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    outpath.write_text(json.dumps(evidence,indent=2,allow_nan=False)+"\n")
    failures = [c["name"] for c in checks if not c["passed"]]
    print(f"{len(checks)-len(failures)}/{len(checks)} review checks passed; evidence: {outpath}")
    for name in failures:
        print("FAIL:",name)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
