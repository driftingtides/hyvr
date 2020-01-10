import sys
import hyvr
import argparse

parser = argparse.ArgumentParser(description="HyVR: Turning your geofantasy into reality!")
parser.add_argument('--overwrite', action='store_true',
        help='whether to overwrite the old run directory')
parser.add_argument('filename', default=0, action="store", nargs='?',
        help='HyVR config file. Defaults to a testcase.')
args = parser.parse_args()


if args.overwrite == False:
    args.overwrite = None
try:
    hyvr.run(args.filename, args.overwrite)
except Exception as e:
    print("Something went wrong:", e, file=sys.stderr)
    if args.filename == 0:
        raise
    else:
        sys.exit(1)
