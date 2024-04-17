import filecmp
import subprocess
from pathlib import Path
import sys

import pytest

from . import GOLDEN_DATA_PATH, WHISCY_PATH


@pytest.fixture
def cons_file():
    return Path(GOLDEN_DATA_PATH, "regression_consadjust", "2SNIE.cons")


@pytest.fixture
def weight_ma_file():
    return Path(WHISCY_PATH, "param", "weight_ma.txt")


@pytest.fixture
def ztable_file():
    return Path(WHISCY_PATH, "param", "ztable.txt")


@pytest.fixture
def acons_file():
    return Path(GOLDEN_DATA_PATH, "regression_consadjust", "2SNIE.acons")


@pytest.fixture
def consadjust_bin():
    return Path(WHISCY_PATH, "bin", "consadjust.py")


@pytest.mark.regression
@pytest.mark.parametrize(
    "scratch_path", ["scratch_regression_consadjust"], indirect=True
)
def test_regression_consadjust2SNIE(
    scratch_path, consadjust_bin, cons_file, weight_ma_file, ztable_file, acons_file
):

    test_prediction_output = Path(scratch_path, "test.prediction")
    cmd_line = f"{sys.executable} {consadjust_bin} {cons_file} {weight_ma_file} {ztable_file} -o {test_prediction_output}"

    result = subprocess.run(
        cmd_line.split(),
        capture_output=True,
        env={"PYTHONPATH": WHISCY_PATH},
    )

    assert result.returncode == 0

    assert filecmp.cmp(test_prediction_output, acons_file)
