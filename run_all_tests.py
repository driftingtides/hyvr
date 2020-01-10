"""
This is intended to be run with coverage:

$ coverage run run_all_tests.py

It runs both the pytest-tests and the MADE test case.
"""

import pytest
import hyvr

pytest.main()
hyvr.run(0)
