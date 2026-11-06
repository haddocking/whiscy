import filecmp
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from . import GOLDEN_DATA_PATH

HADDOCK_BIN = Path(Path(__file__).parent.parent, "src", "whiscy", "cli_haddock.py")

env = os.environ.copy()
env["PYTHONPATH"] = str(
    Path(Path(__file__).parent.parent, "src"),
)

DATA_PATH = Path(GOLDEN_DATA_PATH, "regression_haddock")


def read_list(file_name):
    """Residue numbers of a comma-separated active.txt or passive.txt file."""
    return {int(x) for x in Path(file_name).read_text().split(",") if x}


def run_haddock(scratch_path, *extra):
    """Run whiscy_haddock on 3qi0 chain C in scratch_path."""
    cmd_line = [
        sys.executable,
        str(HADDOCK_BIN),
        str(Path(DATA_PATH, "3qi0_C.pdb")),
        str(Path(DATA_PATH, "3qi0_C.cons")),
        *extra,
    ]
    return subprocess.run(cmd_line, capture_output=True, env=env, cwd=scratch_path)


@pytest.mark.regression
@pytest.mark.parametrize(
    "scratch_path", ["scratch_regression_haddock_surface"], indirect=True
)
def test_regression_3qi0_C_surface_list(scratch_path):
    result = run_haddock(
        scratch_path, "--surface_list", str(Path(DATA_PATH, "3qi0_C.sur"))
    )

    assert result.returncode == 0
    assert filecmp.cmp(
        Path(scratch_path, "haddock-cutoff.pdb"), Path(DATA_PATH, "3qi0_C_haddock.pdb")
    )
    assert filecmp.cmp(
        Path(scratch_path, "active.txt"), Path(DATA_PATH, "3qi0_C_active.txt")
    )
    assert filecmp.cmp(
        Path(scratch_path, "passive.txt"), Path(DATA_PATH, "3qi0_C_passive.txt")
    )


@pytest.mark.regression
@pytest.mark.parametrize(
    "scratch_path", ["scratch_regression_haddock_scores"], indirect=True
)
def test_regression_3qi0_C_scored_residues(scratch_path):
    # Without a surface list, the scored residues are the surface; for 3qi0
    # chain C these are the same residues as in the surface list
    result = run_haddock(scratch_path)

    assert result.returncode == 0
    assert filecmp.cmp(
        Path(scratch_path, "passive.txt"), Path(DATA_PATH, "3qi0_C_passive.txt")
    )


@pytest.mark.regression
@pytest.mark.parametrize(
    "scratch_path", ["scratch_regression_haddock_buried"], indirect=True
)
def test_passive_residues_are_on_the_surface(scratch_path):
    surface = {
        int(line.split()[0])
        for line in Path(DATA_PATH, "3qi0_C.sur").read_text().splitlines()
        if line.strip()
    }
    result = run_haddock(
        scratch_path, "--surface_list", str(Path(DATA_PATH, "3qi0_C.sur"))
    )

    assert result.returncode == 0
    active = read_list(Path(scratch_path, "active.txt"))
    passive = read_list(Path(scratch_path, "passive.txt"))
    assert passive
    assert passive <= surface
    assert not passive & active


@pytest.mark.regression
@pytest.mark.parametrize("scratch_path", ["scratch_haddock_report"], indirect=True)
def test_report(scratch_path):
    result = run_haddock(scratch_path)

    assert result.returncode == 0
    state = json.loads(Path(scratch_path, "state.mvsj").read_text())
    snapshots = state["snapshots"]
    assert [s["metadata"]["key"] for s in snapshots] == [
        "scores",
        "haddock1",
        "haddock2",
    ]
    # Mol* loads the structure relative to the report, in the same directory
    for s in snapshots:
        assert Path(scratch_path, s["root"]["children"][0]["params"]["url"]).exists()
    html = Path(scratch_path, "report.html").read_text()
    assert Path(scratch_path, "active.txt").read_text() in html
    # The input structure is copied next to the report and linked
    assert filecmp.cmp(Path(scratch_path, "3qi0_C.pdb"), Path(DATA_PATH, "3qi0_C.pdb"))
    assert 'href="3qi0_C.pdb"' in html
    # The terminal says how to view the report
    assert b"python -m http.server" in result.stdout


@pytest.mark.regression
@pytest.mark.parametrize("scratch_path", ["scratch_haddock_log"], indirect=True)
def test_log_file(scratch_path):
    # Each command appends to whiscy.log, starting with the version and command
    run_haddock(scratch_path)
    run_haddock(scratch_path, "--no-report")

    log = Path(scratch_path, "whiscy.log").read_text()
    assert log.count("whiscy_haddock [INFO] WHISCY ") == 2
    assert "Report written to" in log
    assert 'href="whiscy.log"' in Path(scratch_path, "report.html").read_text()


@pytest.mark.regression
@pytest.mark.parametrize("scratch_path", ["scratch_haddock_top10"], indirect=True)
def test_top10(scratch_path):
    # The 10 best-scored residues are active, and none of them is passive
    sur = Path(DATA_PATH, "3qi0_C.sur")
    result = run_haddock(scratch_path, "--no-report", "--surface_list", str(sur))

    assert result.returncode == 0
    lines = Path(DATA_PATH, "3qi0_C.cons").read_text().split("\n")
    scores = [(float(x.split()[0]), int(x.split()[1][1:])) for x in lines if x]
    top10 = {resi for _, resi in sorted(scores, reverse=True)[:10]}
    assert read_list(Path(scratch_path, "active_top10.txt")) == top10
    passive = read_list(Path(scratch_path, "passive_top10.txt"))
    assert passive and not passive & top10
    lines = Path(scratch_path, "haddock-top10.pdb").read_text().splitlines()
    atoms = [x for x in lines if x.startswith("ATOM  ")]
    assert {int(x[22:26]) for x in atoms if float(x[60:66]) == 100} == top10
