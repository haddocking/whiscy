import filecmp
import os
import shutil
from pathlib import Path

import pytest

from whiscy.modules import access

from . import GOLDEN_DATA_PATH


@pytest.fixture
def rsa_file():
    return Path(GOLDEN_DATA_PATH, "access", "2sni_E.rsa")


@pytest.fixture
def sur_file():
    return Path(GOLDEN_DATA_PATH, "access", "2sni_E.sur")


@pytest.fixture
def suract_file():
    return Path(GOLDEN_DATA_PATH, "access", "2sni_E.suract")


@pytest.fixture
def lac_file():
    return Path(GOLDEN_DATA_PATH, "access", "2sni_E.lac")


@pytest.fixture
def pdb_file():
    return Path(GOLDEN_DATA_PATH, "access", "2sni_E.pdb")


@pytest.mark.unit
@pytest.mark.parametrize("scratch_path", ["scratch_access"], indirect=True)
def test_create_cutoff_files(scratch_path, sur_file, suract_file, lac_file, rsa_file):

    access.create_cutoff_files(
        rsa_file,
        "2sni",
        "E",
        cutoffs={"sa_pred_cutoff": 15.0, "sa_act_cutoff": 40.0},
        path=scratch_path,
    )

    expected_sur = Path(scratch_path, "2sni_E.sur")
    expected_suract = Path(scratch_path, "2sni_E.suract")
    expected_lac = Path(scratch_path, "2sni_E.lac")

    assert filecmp.cmp(sur_file, expected_sur)
    assert filecmp.cmp(suract_file, expected_suract)
    assert filecmp.cmp(lac_file, expected_lac)


@pytest.mark.unit
@pytest.mark.parametrize("scratch_path", ["scratch_access"], indirect=True)
def test_calculate_accessibility(scratch_path, pdb_file, rsa_file):

    # Move to this directory in order to the output of freesasa be the same
    ori_dir = os.getcwd()
    os.chdir(scratch_path)

    output_rsa = Path(scratch_path, "2sni_E.rsa")

    # Copy the pdb to the scratch directory
    src = pdb_file
    dst = Path(scratch_path, pdb_file.name)
    shutil.copy(src, dst)

    access.calculate_accessibility(
        pdb_file_name=pdb_file.name, output_file_name=output_rsa
    )

    os.chdir(ori_dir)

    assert filecmp.cmp(rsa_file, output_rsa)
