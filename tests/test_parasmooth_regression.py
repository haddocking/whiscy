import filecmp
import sys
import subprocess
from pathlib import Path

import pytest

from . import GOLDEN_DATA_PATH, WHISCY_PATH


@pytest.fixture
def acons_file():
    return Path(GOLDEN_DATA_PATH, "regression_parasmooth", "2SNIE.acons")


@pytest.fixture
def lcons_file():
    return Path(GOLDEN_DATA_PATH, "regression_parasmooth", "2SNIE.lcons")


@pytest.fixture
def rd_file():
    return Path(GOLDEN_DATA_PATH, "regression_parasmooth", "2SNIE.rd")


@pytest.fixture
def par_file():
    return Path(WHISCY_PATH, "param", "parasmooth.par")


@pytest.fixture
def prediction_file():
    return Path(GOLDEN_DATA_PATH, "regression_parasmooth", "2SNIE.parasmooth")


@pytest.fixture
def parasmooth_bin():
    return Path(WHISCY_PATH, "bin", "parasmooth.py")


@pytest.mark.parametrize(
    "scratch_path", ["scratch_regression_parasmooth"], indirect=True
)
def test_regression_2SNIE(
    scratch_path,
    acons_file,
    lcons_file,
    rd_file,
    par_file,
    prediction_file,
    parasmooth_bin,
):

    test_prediction_output = Path(scratch_path, "test.prediction")
    cmd_line = f"{sys.executable} {parasmooth_bin} {acons_file} {lcons_file} {rd_file} {par_file} -o {test_prediction_output}"

    result = subprocess.run(
        cmd_line.split(),
        capture_output=True,
        env={"PYTHONPATH": WHISCY_PATH},
    )

    assert result.returncode == 0

    assert filecmp.cmp(test_prediction_output, prediction_file)
