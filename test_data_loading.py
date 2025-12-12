#!/usr/bin/env python3
"""
VitalDB Data Loading Test
Tests different methods of loading MAP data
"""

import vitaldb
import numpy as np

print("=" * 70)
print("VITALDB DATA LOADING TEST")
print("=" * 70)
print()

case_id = 1
print(f"Testing Case {case_id}:")
print()

vf = vitaldb.VitalFile(case_id)

# Find MAP track
track = 'Solar8000/ART_MBP'
print(f"Track: {track}")
print()

# Method 1: Load with interval parameter
print("Method 1: to_numpy with interval=10")
try:
    data = vf.to_numpy([track], interval=10)
    if data is not None and len(data) > 0:
        print(f"  Result: {len(data[0])} samples")
        if len(data[0]) > 0:
            print(f"  First 10 values: {data[0][:10]}")
    else:
        print("  Result: No data")
except Exception as e:
    print(f"  Error: {e}")

print()

# Method 2: Load raw samples
print("Method 2: to_numpy without interval")
try:
    data = vf.to_numpy([track])
    if data is not None and len(data) > 0:
        print(f"  Result: {len(data[0])} samples")
        if len(data[0]) > 0:
            print(f"  First 10 values: {data[0][:10]}")
            # Manually downsample to 10-second intervals
            downsampled = data[0][::10]  # Every 10th sample
            print(f"  Downsampled (every 10th): {len(downsampled)} samples")
    else:
        print("  Result: No data")
except Exception as e:
    print(f"  Error: {e}")

print()

# Method 3: Using vitaldb.load_case (if available)
print("Method 3: Check vitaldb module functions")
print(f"  Available functions: {[x for x in dir(vitaldb) if not x.startswith('_')]}")

print()
print("=" * 70)
