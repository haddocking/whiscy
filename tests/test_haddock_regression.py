import filecmp
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
        "test_haddock.pdb",
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
        Path(scratch_path, "test_haddock.pdb"), Path(DATA_PATH, "3qi0_C_haddock.pdb")
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
