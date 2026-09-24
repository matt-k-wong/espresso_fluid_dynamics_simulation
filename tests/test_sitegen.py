"""Static viewer built from synthetic bundles only."""

from pathlib import Path

from espresso_m1.cli import main as cli_main

LIN = "examples/shots/synthetic_linear.json"
IRR = "examples/shots/synthetic_irregular.json"


def test_site_two_bundles(tmp_path):
    b1 = tmp_path / "a"
    b2 = tmp_path / "b"
    assert cli_main(["experiment-save", "--shots", LIN, IRR,
                     "--already-wet-t0", "0.0", "--label", "alpha",
                     "--out-dir", str(b1)]) == 0
    assert cli_main(["experiment-save", "--shots", LIN, IRR,
                     "--already-wet-t0", "0.0", "--holdout",
                     "synthetic-linear-01",
                     "--label", "beta", "--out-dir", str(b2)]) == 0
    page = tmp_path / "site" / "index.html"
    assert cli_main(["experiment-site", "--bundles", str(b1), str(b2),
                     "--out", str(page), "--title", "Demo"]) == 0
    text = page.read_text()
    assert "<html" in text and "alpha" in text and "beta" in text
    assert text.count("<option") == 2
    assert "1.2000" in text  # fitted q shown
    assert "observed" in text and "simulated" in text
    assert "Limitations" in text and "Calibration" in text
    assert "synthetic" in text
    import os as _os
    assert _os.path.expanduser("~") not in text  # no builder home leaked
    assert "C:" + chr(92) not in text  # no Windows absolute path
    # Note: the inline SVG xmlns ("http://www.w3.org/...") is a static
    # namespace identifier, not a network fetch; no backend exists.
    assert cli_main(["experiment-site", "--bundles", str(tmp_path / "nope"),
                     "--out", str(tmp_path / "x.html")]) == 2
    import pytest as _pytest
    with _pytest.raises(SystemExit):
        cli_main(["experiment-site", "--bundles",
                  "--out", str(tmp_path / "y.html")])
