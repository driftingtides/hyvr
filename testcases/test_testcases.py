"""
This is a simple testcase for HyVR. It runs the ``made.ini`` testcase with a
seeded random number generator and compares the output files in ``small`` to
the reference output in ``small_reference``.

This can simply be adapted to other test cases by replacing the respective filenames.
"""

import hyvr
import os
import shutil
import numpy as np
import scipy.io as sio
import scipy.io as sio
import filecmp

def run_testcase(filename, runname):

    testcasedir = os.path.relpath(os.path.join(__file__, '../'))
    testfile = os.path.join(testcasedir, filename)

    # remove old output
    if os.path.exists(os.path.join(testcasedir, runname)):
        shutil.rmtree(os.path.join(testcasedir, runname))

    # run seeded
    np.random.seed(0)
    hyvr.run(testfile)

    # check output
    # ============
    outputdir = os.path.join(testcasedir, runname)
    refdir = os.path.join(testcasedir, runname + '_reference')

    # parameter file
    # --------------
    assert filecmp.cmp(os.path.join(outputdir, runname + '_parameters.ini'),
                       os.path.join(refdir, runname + '_parameters.ini'))

    # numpy
    # ------
    numpy_output = np.load(os.path.join(outputdir, runname + '.npz'))
    numpy_ref_output = np.load(os.path.join(refdir, runname + '.npz'))

    assert numpy_output.files == numpy_ref_output.files
    for name in numpy_output.files:
        assert np.all(numpy_output[name] == numpy_ref_output[name])

    numpy_output.close()
    numpy_ref_output.close()

    # matlab
    # ------
    matlab_output = sio.loadmat(os.path.join(outputdir, runname + '.mat'))
    matlab_ref_output = sio.loadmat(os.path.join(refdir, runname + '.mat'))

    del matlab_output['__header__']
    del matlab_ref_output['__header__']

    assert matlab_output.keys() == matlab_ref_output.keys()
    for name in matlab_output.keys():
        assert np.all(matlab_output[name] == matlab_ref_output[name])

    # vtr
    # ---
    # this seems to work
    assert filecmp.cmp(os.path.join(outputdir, runname + '.vtr'),
                       os.path.join(refdir, runname + '.vtr'))

    # TODO: modflow output

    print("Everything okay!")


def test():
    run_testcase('made.ini', 'small')
    # run_testcase('test_lu.ini', 'test_lu')
