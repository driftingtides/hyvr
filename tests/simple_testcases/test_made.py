import os
import numpy as np
import pytest
from tests.run_testcase import run_testcase

@pytest.mark.default
def test_made_seeded():
    """
    Run the simple made testcase with a seeded RNG
    """

    inifile = os.path.join(os.path.dirname(__file__), 'made_simple.ini')
    runname = 'made_simple'

    np.random.seed(42)
    run_testcase(inifile, runname, runname + '_reference')


@pytest.fixture
def mock_rng(monkeypatch):

    def mock_randn(*args):
        return np.zeros(args)

    def mock_rand(*args):
        return 0.5*np.ones(args)

    # randint is only called for choosing a trough structure
    def mock_randint(high):
        return 0 # only flat troughs

    def mock_choice(values, **kwargs):
        return values[0] # use first value

    def mock_normal(loc, scale):
        return loc

    def mock_poisson(lam=1):
        return int(np.round(lam))

    def mock_uniform(low=0.0, high=1.0, size=None):
        if size is None:
            # single value
            return 0.5*(low + high)
        else:
            return np.ones(size)*0.5*(low+high)

    monkeypatch.setattr(np.random, "randn", mock_randn)
    monkeypatch.setattr(np.random, "rand", mock_rand)
    monkeypatch.setattr(np.random, "randint", mock_randint)
    monkeypatch.setattr(np.random, "choice", mock_choice)
    monkeypatch.setattr(np.random, "normal", mock_normal)
    monkeypatch.setattr(np.random, "poisson", mock_poisson)
    monkeypatch.setattr(np.random, "uniform", mock_uniform)
        

@pytest.mark.no_random
def test_made_no_random(mock_rng):
    """
    Run the simple made testcase without any random numbers
    """

    inifile = os.path.join(os.path.dirname(__file__), 'made_simple.ini')
    runname = 'made_simple'

    run_testcase(inifile, runname, runname + '_no_random_reference')

    
