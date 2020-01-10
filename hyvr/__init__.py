import pkg_resources

__version__ = pkg_resources.get_distribution('hyvr').version

import hyvr.sim
from hyvr.sim import run
import hyvr.postprocess
