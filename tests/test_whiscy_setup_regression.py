import filecmp
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from . import FREESASA_PATH, GOLDEN_DATA_PATH, WHISCY_PATH

env = os.environ.copy()
env["PATH"] += os.pathsep + str(FREESASA_PATH)
env["WHISCY_PATH"] = str(WHISCY_PATH)


@pytest.fixture
def hssp_file_1ppe():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "1ppe.hssp")


@pytest.fixture
def pdb_file_1ppe():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "1ppe.pdb")


@pytest.fixture
def conversion_file_1ppe():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "1ppe_I.conv")


@pytest.fixture
def fasta_file_1ppe():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "1ppe_I.fasta")


@pytest.fixture
def lac_file_1ppe():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "1ppe_I.lac")


@pytest.fixture
def protdist_file_1ppe():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "1ppe_I.out")


@pytest.fixture
def pdb_parsed_file_1ppe():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "1ppe_I.pdb")


@pytest.fixture
def alignment_file_1ppe():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "1ppe_I.phylseq")


@pytest.fixture
def rsa_file_1ppe():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "1ppe_I.rsa")


@pytest.fixture
def surface_file_1ppe():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "1ppe_I.sur")


@pytest.fixture
def suract_file_1ppe():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "1ppe_I.suract")


@pytest.fixture
def hssp_file_3qi0():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "3qi0.hssp")


@pytest.fixture
def pdb_file_3qi0():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "3qi0.pdb")


@pytest.fixture
def conversion_file_3qi0():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "3qi0_C.conv")


@pytest.fixture
def fasta_file_3qi0():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "3qi0_C.fasta")


@pytest.fixture
def lac_file_3qi0():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "3qi0_C.lac")


@pytest.fixture
def protdist_file_3qi0():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "3qi0_C.out")


@pytest.fixture
def pdb_parsed_file_3qi0():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "3qi0_C.pdb")


@pytest.fixture
def alignment_file_3qi0():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "3qi0_C.phylseq")


@pytest.fixture
def rsa_file_3qi0():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "3qi0_C.rsa")


@pytest.fixture
def surface_file_3qi0():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "3qi0_C.sur")


@pytest.fixture
def suract_file_3qi0():
    return Path(GOLDEN_DATA_PATH, "regression_whiscy_setup", "3qi0_C.suract")


@pytest.fixture
def whiscy_setup_bin():
    return Path(WHISCY_PATH, "whiscy_setup.py")


@pytest.mark.skip(reason="HSSP BROKEN")
@pytest.mark.parametrize(
    "scratch_path", ["scratch_regression_whiscy_setup_1ppe"], indirect=True
)
def test_regression_1PPEI(
    scratch_path,
    hssp_file_1ppe,
    pdb_file_1ppe,
    conversion_file_1ppe,
    fasta_file_1ppe,
    lac_file_1ppe,
    protdist_file_1ppe,
    pdb_parsed_file_1ppe,
    alignment_file_1ppe,
    rsa_file_1ppe,
    surface_file_1ppe,
    suract_file_1ppe,
    whiscy_setup_bin,
):
    # All the files produced by the whiscy_setup protocol:
    # 1ppe.hssp      1ppe.pdb       1ppe_I.conv    1ppe_I.fasta   1ppe_I.lac
    # 1ppe_I.out     1ppe_I.pdb     1ppe_I.phylseq 1ppe_I.rsa     1ppe_I.sur     1ppe_I.suract

    pred_hssp_file = Path(scratch_path, "1ppe.hssp")
    pred_pdb_file = Path(scratch_path, "1ppe.pdb")
    pred_conversion_file = Path(scratch_path, "1ppe_I.conv")
    pred_fasta_file = Path(scratch_path, "1ppe_I.fasta")
    pred_lac_file = Path(scratch_path, "1ppe_I.lac")
    pred_protdist_file = Path(scratch_path, "1ppe_I.out")
    pred_pdb_parsed_file = Path(scratch_path, "1ppe_I.pdb")
    pred_alignment_file = Path(scratch_path, "1ppe_I.phylseq")
    pred_rsa_file = Path(scratch_path, "1ppe_I.rsa")
    pred_surface_file = Path(scratch_path, "1ppe_I.sur")
    pred_suract_file = Path(scratch_path, "1ppe_I.suract")

    cmd_line = f"{sys.executable} {whiscy_setup_bin} 1ppe i"
    result = subprocess.run(
        cmd_line.split(),
        capture_output=True,
        env=env,
    )

    assert result.returncode == 0, "whiscy_setup failed."

    assert filecmp.cmp(hssp_file_1ppe, pred_hssp_file)
    assert filecmp.cmp(pdb_file_1ppe, pred_pdb_file)
    assert filecmp.cmp(conversion_file_1ppe, pred_conversion_file)
    assert filecmp.cmp(fasta_file_1ppe, pred_fasta_file)
    assert filecmp.cmp(lac_file_1ppe, pred_lac_file)
    assert filecmp.cmp(protdist_file_1ppe, pred_protdist_file)
    assert filecmp.cmp(pdb_parsed_file_1ppe, pred_pdb_parsed_file)
    assert filecmp.cmp(alignment_file_1ppe, pred_alignment_file)
    assert filecmp.cmp(rsa_file_1ppe, pred_rsa_file)
    assert filecmp.cmp(surface_file_1ppe, pred_surface_file)
    assert filecmp.cmp(suract_file_1ppe, pred_suract_file)


@pytest.mark.skip(reason="HSSP BROKEN")
@pytest.mark.parametrize(
    "scratch_path", ["scratch_regression_whiscy_setup_3qi0"], indirect=True
)
def test_regression_3QI0(
    scratch_path,
    hssp_file_3qi0,
    pdb_file_3qi0,
    conversion_file_3qi0,
    fasta_file_3qi0,
    lac_file_3qi0,
    protdist_file_3qi0,
    pdb_parsed_file_3qi0,
    alignment_file_3qi0,
    rsa_file_3qi0,
    surface_file_3qi0,
    suract_file_3qi0,
    whiscy_setup_bin,
):

    pred_hssp_file = Path(scratch_path, "3qi0.hssp")
    pred_pdb_file = Path(scratch_path, "3qi0.pdb")
    pred_conversion_file = Path(scratch_path, "3qi0_C.conv")
    pred_fasta_file = Path(scratch_path, "3qi0_C.fasta")
    pred_lac_file = Path(scratch_path, "3qi0_C.lac")
    pred_protdist_file = Path(scratch_path, "3qi0_C.out")
    pred_pdb_parsed_file = Path(scratch_path, "3qi0_C.pdb")
    pred_alignment_file = Path(scratch_path, "3qi0_C.phylseq")
    pred_rsa_file = Path(scratch_path, "3qi0_C.rsa")
    pred_surface_file = Path(scratch_path, "3qi0_C.sur")
    pred_suract_file = Path(scratch_path, "3qi0_C.suract")

    cmd_line = f"{sys.executable} {whiscy_setup_bin} 3qi0 C"
    result = subprocess.run(
        cmd_line.split(),
        capture_output=True,
        env=env,
    )

    assert result.returncode == 0, "whiscy_setup failed."

    assert filecmp.cmp(hssp_file_3qi0, pred_hssp_file)
    assert filecmp.cmp(pdb_file_3qi0, pred_pdb_file)
    assert filecmp.cmp(conversion_file_3qi0, pred_conversion_file)
    assert filecmp.cmp(fasta_file_3qi0, pred_fasta_file)
    assert filecmp.cmp(lac_file_3qi0, pred_lac_file)
    assert filecmp.cmp(protdist_file_3qi0, pred_protdist_file)
    assert filecmp.cmp(pdb_parsed_file_3qi0, pred_pdb_parsed_file)
    assert filecmp.cmp(alignment_file_3qi0, pred_alignment_file)
    assert filecmp.cmp(rsa_file_3qi0, pred_rsa_file)
    assert filecmp.cmp(surface_file_3qi0, pred_surface_file)
    assert filecmp.cmp(suract_file_3qi0, pred_suract_file)
