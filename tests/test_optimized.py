import pickle
import os
import hyvr
import numpy as np

def test_reindex():
    with open(os.path.join(os.path.dirname(__file__), 'reindex.inray'), 'rb') as f:
        inray = pickle.load(f)
    with open(os.path.join(os.path.dirname(__file__),'reindex.outray'), 'rb') as f:
        outray = pickle.load(f)
    outray_new = hyvr.optimized.reindex(inray)
    assert np.all(outray == outray_new)


