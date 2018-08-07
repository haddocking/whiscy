import os
import shutil
import filecmp
from libwhiscy import access


golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
scratch_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/scratch_access'


def test_create_cutoff_files():
    rsa_file = os.path.join(golden_data_path, 'access', '2sni_E.rsa')
    sur_file = os.path.join(golden_data_path, 'access', '2sni_E.sur')
    suract_file = os.path.join(golden_data_path, 'access', '2sni_E.suract')
    lac_file = os.path.join(golden_data_path, 'access', '2sni_E.lac')
    cutoffs = {"sa_pred_cutoff": 15.0, "sa_act_cutoff": 40.0}

    if os.path.exists(scratch_path):
        shutil.rmtree(scratch_path)
    os.mkdir(scratch_path)

    access.create_cutoff_files(rsa_file, '2sni', 'E', cutoffs, path=scratch_path)

    assert filecmp.cmp(sur_file, os.path.join(scratch_path, '2sni_E.sur'))
    assert filecmp.cmp(suract_file, os.path.join(scratch_path, '2sni_E.suract'))
    assert filecmp.cmp(lac_file, os.path.join(scratch_path, '2sni_E.lac'))

    shutil.rmtree(scratch_path)


def test_calculate_accessibility():
    pdb_file = '2sni_E.pdb'
    rsa_file = os.path.join(golden_data_path, 'access', '2sni_E.rsa')

    if os.path.exists(scratch_path):
        shutil.rmtree(scratch_path)
    os.mkdir(scratch_path)

    # Move to this directory in order to the output of freesasa be the same
    os.chdir(os.path.join(golden_data_path, 'access'))

    access.calculate_accessibility(pdb_file, os.path.join(scratch_path, '2sni_E.rsa'))

    assert filecmp.cmp(rsa_file, os.path.join(scratch_path, '2sni_E.rsa'))

    shutil.rmtree(scratch_path)
