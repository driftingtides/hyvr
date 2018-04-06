import sys
import hyvr

if len(sys.argv) <= 1:
    print("Please pass a *.ini file as argument.")
    sys.exit(1)

try:
    hyvr.run(sys.argv[1])
except Exception as e:
    print(e, file=sys.stderr)
    sys.exit(1)
