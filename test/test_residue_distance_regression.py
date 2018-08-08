import os
import shutil
import filecmp


golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
scratch_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/scratch_regression_residue_distance'


def test_rd_regression_1PPEI():
    pdb_file = os.path.join(golden_data_path, 'regression_residue_distance', '1ppe_I.pdb')
    conversion_file = os.path.join(golden_data_path, 'regression_residue_distance', '1ppe_I.conv')
    rd_file = os.path.join(golden_data_path, 'regression_residue_distance', '1ppe_I.rd')

    rd_bin =  os.path.join(os.environ['WHISCY_PATH'], 'bin', 'residue_distance.py')
    if os.path.exists(scratch_path):
        shutil.rmtree(scratch_path)
    os.mkdir(scratch_path)

    test_prediction_output = os.path.join(scratch_path, 'test.prediction')
    cmd_line = "{0} {1} {2} {3} > /dev/null 2>&1".format(rd_bin,
                                                         pdb_file,
                                                         conversion_file,
                                                         test_prediction_output)
    os.system(cmd_line)

    assert filecmp.cmp(test_prediction_output, rd_file)

    shutil.rmtree(scratch_path)
