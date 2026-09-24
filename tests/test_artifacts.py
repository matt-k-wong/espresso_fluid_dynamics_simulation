"""Task C regression tests: trustworthy artifacts and packaging."""

import csv
import json
import os
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

import pytest

from espresso_m1.cli import main as cli_main

ROOT = Path(__file__).resolve().parents[1]


def run_cli(argv):
    return cli_main(argv)


def test_success_failure_success_artifact_hygiene(tmp_path):
    out = tmp_path / "case"
    common = ["run", "--config", "configs/uniform.json", "--out-dir", str(out)]
    assert run_cli(common) == 0
    assert (out / "summary.json").exists()
    assert (out / "fields.npz").exists()
    assert not (out / "failure.json").exists()

    assert run_cli(common + ["--set", "max_iter=1"]) == 1
    left = sorted(p.name for p in out.iterdir())
    assert "summary.json" not in left and "fields.npz" not in left
    assert "failure.json" in left
    failure = json.loads((out / "failure.json").read_text())
    assert failure["converged"] is False and "params" in failure

    assert run_cli(common) == 0
    assert (out / "summary.json").exists()
    assert (out / "fields.npz").exists()
    assert not (out / "failure.json").exists()


def test_run_leaves_no_foreign_files_alone(tmp_path):
    out = tmp_path / "case"
    out.mkdir()
    sentinel = out / "notes.txt"
    sentinel.write_text("user file")
    assert run_cli(["run", "--config", "configs/uniform.json",
                    "--out-dir", str(out)]) == 0
    assert sentinel.read_text() == "user file"


def test_zero_flow_summary_is_strict_json(tmp_path):
    out = tmp_path / "zero"
    rc = run_cli(["run", "--config", "configs/default.json",
                  "--set", "p_in=0", "--out-dir", str(out)])
    assert rc == 0

    def reject_constant(value):
        raise ValueError(f"non-standard constant {value}")

    summary = json.loads((out / "summary.json").read_text(),
                         parse_constant=reject_constant)
    assert summary["R_hyd"] is None
    assert summary["R_hyd_status"] == "undefined_zero_flow"
    assert summary["Q_m3s"] == 0.0


def test_empty_sweep_is_not_success(tmp_path):
    rc = run_cli(["sweep", "--config", "configs/default.json",
                  "--p-in-bar", "3:1:1", "--out", str(tmp_path / "e.csv")])
    assert rc != 0


def test_failed_sweep_points_recorded_with_nonzero_exit(tmp_path):
    csv_path = tmp_path / "sweep.csv"
    rc = run_cli(["sweep", "--config", "configs/default.json",
                  "--set", "max_iter=1", "--p-in-bar", "1,2",
                  "--out", str(csv_path)])
    assert rc != 0
    with open(csv_path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 2
    assert all(r["converged"] == "False" for r in rows)
    assert all(r["error"] for r in rows)


def _build_wheel(dist_dir: Path) -> Path:
    """Build a wheel the way the review probe does (isolated build)."""
    with tempfile.TemporaryDirectory(prefix="espresso-pkg-") as staging:
        project = Path(staging) / "project"
        project.mkdir()
        for name in ("pyproject.toml",):
            (project / name).write_bytes((ROOT / name).read_bytes())
        import shutil
        shutil.copytree(ROOT / "src", project / "src",
                        ignore=shutil.ignore_patterns("__pycache__",
                                                      "*.egg-info"))
        build = subprocess.run(
            [sys.executable, "-m", "pip", "wheel", str(project),
             "--no-deps", "-w", str(dist_dir)],
            capture_output=True, text=True, cwd=str(ROOT),
        )
        if build.returncode != 0:
            # Offline fallback: reuse the local toolchain without isolation.
            build = subprocess.run(
                [sys.executable, "-m", "pip", "wheel", str(project),
                 "--no-deps", "--no-build-isolation", "-w", str(dist_dir)],
                capture_output=True, text=True, cwd=str(ROOT),
            )
        assert build.returncode == 0, build.stderr[-2000:]
    wheels = list(dist_dir.glob("*.whl"))
    assert len(wheels) == 1
    return wheels[0]


def test_wheel_contains_catalog_and_imports_off_tree(tmp_path):
    dist = tmp_path / "dist"
    dist.mkdir()
    wheel = _build_wheel(dist)
    with zipfile.ZipFile(wheel) as archive:
        assert "espresso_m1/catalog.json" in archive.namelist()

    target = tmp_path / "install"
    install = subprocess.run(
        [sys.executable, "-m", "pip", "install", str(wheel), "--no-deps",
         "--no-build-isolation", "-t", str(target), "-q"],
        capture_output=True, text=True,
    )
    assert install.returncode == 0, install.stderr[-2000:]

    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    env["PYTHONPATH"] = str(target)
    # Run outside the source tree with the source tree unreachable.
    check = subprocess.run(
        [sys.executable, "-c",
         "from espresso_m1.params import load_catalog; "
         "c = load_catalog(); "
         "assert 'mu' in c['parameters']; "
         "print('catalog ok:', len(c['parameters']))"],
        capture_output=True, text=True, cwd=str(tmp_path), env=env,
    )
    assert check.returncode == 0, check.stderr[-2000:]
    assert "catalog ok" in check.stdout
