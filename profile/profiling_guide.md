# Guide to profiling HyVR

The script `profile_hyvr.py` can be run to profile HyVR. It creates an output
file called `hyvr.profile`, which can be opened with `pstats`:
```
python -m pstats hyvr.profile
```
