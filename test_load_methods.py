#!/usr/bin/env python3
"""
Test VitalDB's load_case and vital_recs functions
"""

import vitaldb

print("=" * 70)
print("TESTING VITALDB load_case AND vital_recs")
print("=" * 70)
print()

case_id = 1

# Test load_case
print("Method 1: vitaldb.load_case()")
try:
    data = vitaldb.load_case(case_id, ['Solar8000/ART_MBP', 'Solar8000/HR'])
    print(f"  Type: {type(data)}")
    if isinstance(data, dict):
        print(f"  Keys: {list(data.keys())}")
        for key, val in data.items():
            if hasattr(val, '__len__'):
                print(f"  {key}: {len(val)} samples")
                if len(val) > 0:
                    print(f"    First 10: {val[:10]}")
    elif isinstance(data, list):
        print(f"  Length: {len(data)}")
        if len(data) > 0:
            print(f"  First item: {data[0]}")
except Exception as e:
    print(f"  Error: {e}")

print()

# Test vital_recs
print("Method 2: vitaldb.vital_recs()")
try:
    # vital_recs might take case_id and track name
    data = vitaldb.vital_recs(case_id, 'Solar8000/ART_MBP')
    print(f"  Type: {type(data)}")
    if hasattr(data, '__len__'):
        print(f"  Length: {len(data)}")
        if len(data) > 0:
            print(f"  First 10: {data[:10]}")
    else:
        print(f"  Value: {data}")
except Exception as e:
    print(f"  Error: {e}")

print()

# Test VitalFile with different approach
print("Method 3: VitalFile - check attributes")
try:
    vf = vitaldb.VitalFile(case_id)
    print(f"  VitalFile attributes: {[x for x in dir(vf) if not x.startswith('_')]}")

    # Try to get duration
    if hasattr(vf, 'get_duration'):
        duration = vf.get_duration()
        print(f"  Duration: {duration} seconds")

    # Try to get samples differently
    if hasattr(vf, 'get_samples'):
        print("\n  Testing get_samples():")
        samples = vf.get_samples('Solar8000/ART_MBP')
        if samples:
            print(f"    Type: {type(samples)}")
            if hasattr(samples, '__len__'):
                print(f"    Length: {len(samples)}")
                if len(samples) > 0:
                    print(f"    First 10: {samples[:10]}")

except Exception as e:
    print(f"  Error: {e}")

print()
print("=" * 70)
