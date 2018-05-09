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

def test_made():

    testcasedir = os.path.relpath(os.path.join(__file__, '../'))
    testfile = os.path.join(testcasedir, 'made.ini')

    # remove old output
    shutil.rmtree(os.path.join(testcasedir, 'small'))

    # run seeded
    np.random.seed(0)
    hyvr.run(testfile)

    # check output
    # ============
    outputdir = os.path.join(testcasedir, 'small')
    refdir = os.path.join(testcasedir, 'small_reference')

    # parameter file
    # --------------
    assert filecmp.cmp(os.path.join(outputdir, 'small_parameters.ini'),
                       os.path.join(refdir, 'small_parameters.ini'))

    # numpy
    # ------
    numpy_output = np.load(os.path.join(outputdir, 'small.npz'))
    numpy_ref_output = np.load(os.path.join(refdir, 'small.npz'))

    assert numpy_output.files == numpy_ref_output.files
    for name in numpy_output.files:
        assert np.all(numpy_output[name] == numpy_ref_output[name])

    numpy_output.close()
    numpy_ref_output.close()

    # matlab
    # ------
    matlab_output = sio.loadmat(os.path.join(outputdir, 'small.mat'))
    matlab_ref_output = sio.loadmat(os.path.join(refdir, 'small.mat'))

    del matlab_output['__header__']
    del matlab_ref_output['__header__']

    assert matlab_output.keys() == matlab_ref_output.keys()
    for name in matlab_output.keys():
        assert np.all(matlab_output[name] == matlab_ref_output[name])

    # vtr
    # ---
    # this seems to work
    assert filecmp.cmp(os.path.join(outputdir, 'small.vtr'),
                       os.path.join(refdir, 'small.vtr'))

    # TODO: modflow output

    print("Everything okay!")
