import sys
import hyvr

if len(sys.argv) <= 1:
    filename = 0
else:
    filename = sys.argv[1]
try:
    hyvr.run(filename)
except Exception as e:
    raise e
    print(e, file=sys.stderr)
    sys.exit(1)
