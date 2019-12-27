"""
This runs several test cases for integration testing. It sets a fixed seed for
the RNG so that its possible to directly compare the new output to a reference
output.

This can simply be adapted to other test cases by replacing the respective filenames.
"""

import pytest
import hyvr
import os
import shutil
import numpy as np
import scipy.io as sio
import scipy.io as sio
import filecmp

from tests.run_testcase import run_testcase

if __name__ == '__main__':
    run_testcase('made.ini', 'made', 'made_reference')
    # run_testcase('test_lu.ini', 'test_lu')
