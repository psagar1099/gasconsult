#!/usr/bin/env python3
"""
VitalDB Diagnostic Script
Tests VitalDB library and shows what data is available
"""

import vitaldb

print("=" * 70)
print("VITALDB DIAGNOSTIC TEST")
print("=" * 70)
print()

# Test with a few known case IDs
test_cases = [1, 2, 3, 100, 200]

for case_id in test_cases:
    print(f"\nTesting Case {case_id}:")
    print("-" * 50)

    try:
        # Load without specifying tracks
        vf = vitaldb.VitalFile(case_id)

        # Get all available tracks
        tracks = vf.get_track_names()

        if not tracks:
            print("  ✗ No tracks available")
            continue

        print(f"  ✓ {len(tracks)} tracks available")

        # Show first 10 tracks
        print("  Sample tracks:")
        for track in tracks[:10]:
            print(f"    - {track}")

        # Look for MAP-related tracks
        map_tracks = [t for t in tracks if 'MBP' in t or 'MAP' in t or 'ABP' in t]
        if map_tracks:
            print(f"  ✓ Found {len(map_tracks)} MAP tracks:")
            for track in map_tracks:
                print(f"    - {track}")

                # Try to load a small sample
                try:
                    data = vf.to_numpy([track], interval=10)
                    if data is not None and len(data) > 0:
                        print(f"      → {len(data[0])} samples available")
                    else:
                        print(f"      → No data")
                except Exception as e:
                    print(f"      → Error loading: {str(e)[:50]}")
        else:
            print("  ✗ No MAP tracks found")

    except Exception as e:
        print(f"  ✗ Error: {str(e)[:100]}")

print()
print("=" * 70)
print("TEST COMPLETE")
print("=" * 70)
