import os
import shutil
import filecmp


golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
scratch_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/scratch_regression_parasmooth'

def test_regression_2SNIE():
    acons_file = os.path.join(golden_data_path, 'regression_parasmooth', '2SNIE.acons')
    lcons_file = os.path.join(golden_data_path, 'regression_parasmooth', '2SNIE.lcons')
    rd_file = os.path.join(golden_data_path, 'regression_parasmooth', '2SNIE.rd')
    par_file = os.path.join(os.environ['WHISCY_PATH'], 'param', 'parasmooth.par')
    prediction_file = os.path.join(golden_data_path, 'regression_parasmooth', '2SNIE.parasmooth')

    parasmooth_bin =  os.path.join(os.environ['WHISCY_PATH'], 'bin', 'parasmooth.py')
    if os.path.exists(scratch_path):
        shutil.rmtree(scratch_path)
    os.mkdir(scratch_path)

    test_prediction_output = os.path.join(scratch_path, 'test.prediction')
    cmd_line = "{0} {1} {2} {3} {4} -o {5} > /dev/null 2>&1".format(parasmooth_bin,
                                                                    acons_file,
                                                                    lcons_file,
                                                                    rd_file,
                                                                    par_file,
                                                                    test_prediction_output)
    os.system(cmd_line)

    assert filecmp.cmp(test_prediction_output, prediction_file)

    shutil.rmtree(scratch_path)
