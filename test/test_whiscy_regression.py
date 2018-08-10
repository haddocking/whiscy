import os
import shutil
import filecmp


golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
scratch_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/scratch_regression_whiscy'


def test_regression_2SNIE():
    surface_file = os.path.join(golden_data_path, 'regression_whiscy', '2SNIE.sur')
    conversion_file = os.path.join(golden_data_path, 'regression_whiscy', '2SNIE.conv')
    alignment_file = os.path.join(golden_data_path, 'regression_whiscy', '2SNIE.phylseq')
    distance_file = os.path.join(golden_data_path, 'regression_whiscy', '2SNIE.out')
    prediction_file = os.path.join(golden_data_path, 'regression_whiscy', '2SNIE.whiscy')

    whiscy_bin =  os.environ['WHISCY_BIN']
    if os.path.exists(scratch_path):
        shutil.rmtree(scratch_path)
    os.mkdir(scratch_path)

    test_prediction_output = os.path.join(scratch_path, 'test.prediction')
    cmd_line = "{0} {1} {2} {3} {4} -o {5} > /dev/null 2>&1".format(whiscy_bin,
                                                                    surface_file,
                                                                    conversion_file,
                                                                    alignment_file,
                                                                    distance_file,
                                                                    test_prediction_output)
    os.system(cmd_line)

    assert filecmp.cmp(test_prediction_output, prediction_file)

    shutil.rmtree(scratch_path)
