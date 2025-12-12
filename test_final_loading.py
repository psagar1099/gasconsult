#!/usr/bin/env python3
"""
Final VitalDB data loading test
Testing the methods that actually work
"""

import vitaldb
import numpy as np

print("=" * 70)
print("FINAL VITALDB DATA LOADING TEST")
print("=" * 70)
print()

case_id = 1

# Test 1: vitaldb.load_case with tracks specified
print("Test 1: vitaldb.load_case() with track names")
try:
    data = vitaldb.load_case(case_id, ['Solar8000/ART_MBP', 'Solar8000/HR'])
    print(f"  Type: {type(data)}")
    print(f"  Shape: {data.shape}")
    print(f"  First 10 rows:")
    print(data[:10])
    print()
except Exception as e:
    print(f"  Error: {e}")
    print()

# Test 2: VitalFile.load_opendata()
print("Test 2: VitalFile.load_opendata()")
try:
    vf = vitaldb.VitalFile()
    vf.load_opendata(case_id)

    tracks = vf.get_track_names()
    print(f"  Loaded {len(tracks)} tracks")

    # Find MAP track
    map_tracks = [t for t in tracks if 'MBP' in t or 'MAP' in t]
    print(f"  MAP tracks: {map_tracks}")

    if map_tracks:
        # Try to get samples with interval
        track = map_tracks[0]
        print(f"\n  Testing track: {track}")
        samples = vf.get_samples(track, interval=10)
        print(f"  Samples type: {type(samples)}")
        if hasattr(samples, '__len__'):
            print(f"  Number of samples: {len(samples)}")
            if len(samples) > 0:
                print(f"  First 10 samples: {samples[:10]}")
except Exception as e:
    print(f"  Error: {e}")
    print()

# Test 3: VitalFile.to_pandas()
print("\nTest 3: VitalFile.to_pandas()")
try:
    vf = vitaldb.VitalFile()
    vf.load_opendata(case_id)

    df = vf.to_pandas(['Solar8000/ART_MBP', 'Solar8000/HR'], interval=10)
    print(f"  Type: {type(df)}")
    print(f"  Shape: {df.shape}")
    print(f"  Columns: {list(df.columns)}")
    print(f"\n  First 10 rows:")
    print(df.head(10))
except Exception as e:
    print(f"  Error: {e}")

print()
print("=" * 70)
