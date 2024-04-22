import filecmp
import os
import subprocess
import sys
from pathlib import Path

import pytest

from . import GOLDEN_DATA_PATH

PARASMOOTH_BIN = Path(
    Path(__file__).parent.parent, "src", "whiscy", "cli_parasmooth.py"
)

env = os.environ.copy()
env["PYTHONPATH"] = str(
    Path(Path(__file__).parent.parent, "src"),
)


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
def prediction_file():
    return Path(GOLDEN_DATA_PATH, "regression_parasmooth", "2SNIE.parasmooth")


@pytest.mark.regression
@pytest.mark.parametrize(
    "scratch_path", ["scratch_regression_parasmooth"], indirect=True
)
def test_regression_2SNIE(
    scratch_path,
    acons_file,
    lcons_file,
    rd_file,
    prediction_file,
):

    test_prediction_output = Path(scratch_path, "test.prediction")
    cmd_line = f"{sys.executable} {PARASMOOTH_BIN} {acons_file} {lcons_file} {rd_file} -o {test_prediction_output}"

    result = subprocess.run(
        cmd_line.split(),
        capture_output=True,
        env=env,
    )

    assert result.returncode == 0

    assert filecmp.cmp(test_prediction_output, prediction_file)
