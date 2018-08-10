import os
import shutil
import filecmp


golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
scratch_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/scratch_regression_consadjust'


def test_regression_consadjust2SNIE():
    cons_file = os.path.join(golden_data_path, 'regression_consadjust', '2SNIE.cons')
    weight_ma_file = os.path.join(os.environ['WHISCY_PATH'], 'param', 'weight_ma.txt')
    ztable_file = os.path.join(os.environ['WHISCY_PATH'], 'param', 'ztable.txt')
    acons_file = os.path.join(golden_data_path, 'regression_consadjust', '2SNIE.acons')

    consadjust_bin =  os.path.join(os.environ['WHISCY_PATH'], 'bin', 'consadjust.py')
    if os.path.exists(scratch_path):
        shutil.rmtree(scratch_path)
    os.mkdir(scratch_path)

    test_prediction_output = os.path.join(scratch_path, 'test.prediction')
    cmd_line = "{0} {1} {2} {3} -o {4} > /dev/null 2>&1".format(consadjust_bin,
                                                                cons_file,
                                                                weight_ma_file,
                                                                ztable_file,
                                                                test_prediction_output)
    os.system(cmd_line)

    assert filecmp.cmp(test_prediction_output, acons_file)

    shutil.rmtree(scratch_path)
