import filecmp
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

CONSADJUST_BIN = Path(
    Path(__file__).parent.parent, "src", "whiscy", "cli_consadjust.py"
)


@pytest.fixture
def cons_file():
    return Path(GOLDEN_DATA_PATH, "regression_consadjust", "2SNIE.cons")


@pytest.fixture
def acons_file():
    return Path(GOLDEN_DATA_PATH, "regression_consadjust", "2SNIE.acons")


# @pytest.fixture
# def consadjust_bin():
#     return Path(WHISCY_PATH, "bin", "consadjust.py")


@pytest.mark.regression
@pytest.mark.parametrize(
    "scratch_path", ["scratch_regression_consadjust"], indirect=True
)
def test_regression_consadjust2SNIE(scratch_path, cons_file, acons_file):

    test_prediction_output = Path(scratch_path, "test.prediction")
    cmd_line = (
        f"{sys.executable} {CONSADJUST_BIN} {cons_file} -o {test_prediction_output}"
    )

    result = subprocess.run(
        cmd_line.split(),
        capture_output=True,
        env=env,
    )

    assert result.returncode == 0

    assert filecmp.cmp(test_prediction_output, acons_file)
