import os
import shutil
import filecmp


golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
scratch_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/scratch_regression_whiscy_setup'


def test_regression_1PPEI():
    # All the files produced by the whiscy_setup protocol:
    # 1ppe.hssp      1ppe.pdb       1ppe_I.conv    1ppe_I.fasta   1ppe_I.lac     
    # 1ppe_I.out     1ppe_I.pdb     1ppe_I.phylseq 1ppe_I.rsa     1ppe_I.sur     1ppe_I.suract
    hssp_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '1ppe.hssp')
    pdb_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '1ppe.pdb')
    conversion_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '1ppe_I.conv')
    fasta_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '1ppe_I.fasta')
    lac_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '1ppe_I.lac')
    protdist_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '1ppe_I.out')
    pdb_parsed_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '1ppe_I.pdb')
    alignment_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '1ppe_I.phylseq')
    rsa_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '1ppe_I.rsa')
    surface_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '1ppe_I.sur')
    suract_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '1ppe_I.suract')

    pred_hssp_file = os.path.join(scratch_path, '1ppe.hssp')
    pred_pdb_file = os.path.join(scratch_path, '1ppe.pdb')
    pred_conversion_file = os.path.join(scratch_path, '1ppe_I.conv')
    pred_fasta_file = os.path.join(scratch_path, '1ppe_I.fasta')
    pred_lac_file = os.path.join(scratch_path, '1ppe_I.lac')
    pred_protdist_file = os.path.join(scratch_path, '1ppe_I.out')
    pred_pdb_parsed_file = os.path.join(scratch_path, '1ppe_I.pdb')
    pred_alignment_file = os.path.join(scratch_path, '1ppe_I.phylseq')
    pred_rsa_file = os.path.join(scratch_path, '1ppe_I.rsa')
    pred_surface_file = os.path.join(scratch_path, '1ppe_I.sur')
    pred_suract_file = os.path.join(scratch_path, '1ppe_I.suract')
    
    whiscy_setup_bin =  os.path.join(os.environ['WHISCY_PATH'], 'whiscy_setup.py')
    if os.path.exists(scratch_path):
        shutil.rmtree(scratch_path)
    os.mkdir(scratch_path)
    os.chdir(scratch_path)

    cmd_line = "{0} 1ppe i > /dev/null 2>&1".format(whiscy_setup_bin)
    os.system(cmd_line)

    assert filecmp.cmp(hssp_file, pred_hssp_file)
    assert filecmp.cmp(pdb_file, pred_pdb_file)
    assert filecmp.cmp(conversion_file, pred_conversion_file)
    assert filecmp.cmp(fasta_file, pred_fasta_file)
    assert filecmp.cmp(lac_file, pred_lac_file)
    assert filecmp.cmp(protdist_file, pred_protdist_file)
    assert filecmp.cmp(pdb_parsed_file, pred_pdb_parsed_file)
    assert filecmp.cmp(alignment_file, pred_alignment_file)
    assert filecmp.cmp(rsa_file, pred_rsa_file)
    assert filecmp.cmp(surface_file, pred_surface_file)
    assert filecmp.cmp(suract_file, pred_suract_file)

    shutil.rmtree(scratch_path)


def test_regression_3QIC():
    # All the files produced by the whiscy_setup protocol
    hssp_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '3qi0.hssp')
    pdb_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '3qi0.pdb')
    conversion_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '3qi0_C.conv')
    fasta_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '3qi0_C.fasta')
    lac_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '3qi0_C.lac')
    protdist_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '3qi0_C.out')
    pdb_parsed_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '3qi0_C.pdb')
    alignment_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '3qi0_C.phylseq')
    rsa_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '3qi0_C.rsa')
    surface_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '3qi0_C.sur')
    suract_file = os.path.join(golden_data_path, 'regression_whiscy_setup', '3qi0_C.suract')

    pred_hssp_file = os.path.join(scratch_path, '3qi0.hssp')
    pred_pdb_file = os.path.join(scratch_path, '3qi0.pdb')
    pred_conversion_file = os.path.join(scratch_path, '3qi0_C.conv')
    pred_fasta_file = os.path.join(scratch_path, '3qi0_C.fasta')
    pred_lac_file = os.path.join(scratch_path, '3qi0_C.lac')
    pred_protdist_file = os.path.join(scratch_path, '3qi0_C.out')
    pred_pdb_parsed_file = os.path.join(scratch_path, '3qi0_C.pdb')
    pred_alignment_file = os.path.join(scratch_path, '3qi0_C.phylseq')
    pred_rsa_file = os.path.join(scratch_path, '3qi0_C.rsa')
    pred_surface_file = os.path.join(scratch_path, '3qi0_C.sur')
    pred_suract_file = os.path.join(scratch_path, '3qi0_C.suract')
    
    whiscy_setup_bin =  os.path.join(os.environ['WHISCY_PATH'], 'whiscy_setup.py')
    if os.path.exists(scratch_path):
        shutil.rmtree(scratch_path)
    os.mkdir(scratch_path)
    os.chdir(scratch_path)

    cmd_line = "{0} 3qi0 C > /dev/null 2>&1".format(whiscy_setup_bin)
    os.system(cmd_line)

    assert filecmp.cmp(hssp_file, pred_hssp_file)
    assert filecmp.cmp(pdb_file, pred_pdb_file)
    assert filecmp.cmp(conversion_file, pred_conversion_file)
    assert filecmp.cmp(fasta_file, pred_fasta_file)
    assert filecmp.cmp(lac_file, pred_lac_file)
    assert filecmp.cmp(protdist_file, pred_protdist_file)
    assert filecmp.cmp(pdb_parsed_file, pred_pdb_parsed_file)
    assert filecmp.cmp(alignment_file, pred_alignment_file)
    assert filecmp.cmp(rsa_file, pred_rsa_file)
    assert filecmp.cmp(surface_file, pred_surface_file)
    assert filecmp.cmp(suract_file, pred_suract_file)

    shutil.rmtree(scratch_path)
