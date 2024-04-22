import os
import subprocess
import sys
from pathlib import Path

import pytest

from . import GOLDEN_DATA_PATH

env = os.environ.copy()
env["PYTHONPATH"] = str(
    Path(Path(__file__).parent.parent, "src"),
)

RESDIST_BIN = Path(Path(__file__).parent.parent, "src", "whiscy", "cli_resdist.py")


@pytest.fixture
def pdb_file():
    return Path(GOLDEN_DATA_PATH, "regression_residue_distance", "1ppe_I.pdb")


@pytest.fixture
def conversion_file():
    return Path(GOLDEN_DATA_PATH, "regression_residue_distance", "1ppe_I.conv")


@pytest.fixture
def rd_file():
    return Path(GOLDEN_DATA_PATH, "regression_residue_distance", "1ppe_I.rd")


@pytest.mark.regression
@pytest.mark.parametrize(
    "scratch_path", ["scratch_regression_residue_distance"], indirect=True
)
def test_rd_regression_1PPEI(scratch_path, pdb_file, conversion_file, rd_file):

    test_prediction_output = Path(scratch_path, "test.prediction")
    cmd_line = f"{sys.executable} {RESDIST_BIN} {pdb_file} {conversion_file} {test_prediction_output}"

    result = subprocess.run(
        cmd_line.split(),
        capture_output=True,
        env=env,
    )

    assert result.returncode == 0

    with open(test_prediction_output, "r") as f1, open(rd_file, "r") as f2:
        f1_lines = f1.readlines()
        f2_lines = f2.readlines()

        assert len(f1_lines) == len(f2_lines)

        # Check if the lines are the same
        for l1, l2 in zip(f1_lines, f2_lines):
            r1_1, r1_2, v1 = l1.split()
            r2_1, r2_2, v2 = l2.split()

            assert r1_1 == r2_1
            assert r1_2 == r2_2
            assert float(v1) == pytest.approx(float(v2), abs=0.001)
