import filecmp
import os
import shutil
import subprocess
import os
import sys
from pathlib import Path

import pytest

from . import GOLDEN_DATA_PATH

env = os.environ.copy()
env["PYTHONPATH"] = str(
    Path(Path(__file__).parent.parent, "src"),
)

WHISCY_BIN = Path(Path(__file__).parent.parent, "src", "whiscy", "cli.py")


@pytest.fixture
def surface_file_2snie():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy", "2SNIE.sur")


@pytest.fixture
def conversion_file_2snie():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy", "2SNIE.conv")


@pytest.fixture
def alignment_file_2snie():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy", "2SNIE.phylseq")


@pytest.fixture
def distance_file_2snie():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy", "2SNIE.out")


@pytest.fixture
def prediction_file_2snie():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy", "2SNIE.whiscy")


@pytest.fixture
def surface_file_3qi0_c():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy", "3qi0_C.sur")


@pytest.fixture
def conversion_file_3qi0_c():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy", "3qi0_C.conv")


@pytest.fixture
def alignment_file_3qi0_c():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy", "3qi0_C.phylseq")


@pytest.fixture
def distance_file_3qi0_c():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy", "3qi0_C.out")


@pytest.fixture
def prediction_file_3qi0_c():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy", "3qi0_C.cons")


@pytest.mark.regression
@pytest.mark.parametrize(
    "scratch_path", ["scratch_regression_whiscy_2snie"], indirect=True
)
def test_regression_2SNIE(
    scratch_path,
    surface_file_2snie,
    conversion_file_2snie,
    alignment_file_2snie,
    distance_file_2snie,
    prediction_file_2snie,
):

    test_prediction_output = Path(scratch_path, "test.prediction")
    cmd_line = f"{sys.executable} {WHISCY_BIN} {surface_file_2snie} {conversion_file_2snie} {alignment_file_2snie} {distance_file_2snie} -o {test_prediction_output}"

    result = subprocess.run(
        cmd_line.split(),
        capture_output=True,
        env=env,
    )

    assert result.returncode == 0

    assert filecmp.cmp(test_prediction_output, prediction_file_2snie)


@pytest.mark.regression
@pytest.mark.parametrize(
    "scratch_path", ["scratch_regression_whiscy_3qio_c"], indirect=True
)
def test_regression_3QI0C(
    scratch_path,
    surface_file_3qi0_c,
    conversion_file_3qi0_c,
    alignment_file_3qi0_c,
    distance_file_3qi0_c,
    prediction_file_3qi0_c,
):

    test_prediction_output = Path(scratch_path, "test.prediction")
    cmd_line = f"{sys.executable} {WHISCY_BIN} {surface_file_3qi0_c} {conversion_file_3qi0_c} {alignment_file_3qi0_c} {distance_file_3qi0_c} -o {test_prediction_output}"

    result = subprocess.run(
        cmd_line.split(),
        capture_output=True,
        env=env,
    )

    assert result.returncode == 0

    assert filecmp.cmp(test_prediction_output, prediction_file_3qi0_c)
