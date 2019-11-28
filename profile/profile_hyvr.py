import cProfile
import hyvr

cProfile.run("hyvr.run('made.ini', True)", "hyvr.profile")
