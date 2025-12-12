#!/usr/bin/env python3
"""
VitalDB IOH Downloader - Working Version
Uses vitaldb.load_case() which actually returns data!
"""

import vitaldb
import numpy as np
import pickle
from typing import List, Dict

# IOH threshold
IOH_THRESHOLD = 65
PREDICTION_WINDOW_SEC = 300

def download_vitaldb_data(max_cases: int = 100) -> Dict:
    """Download VitalDB data using load_case()"""

    print("=" * 70)
    print("VITALDB IOH DATASET CREATION")
    print("=" * 70)
    print()

    all_features = []
    ioh_count = 0
    processed_count = 0

    # Try cases 1-6388
    print(f"Processing up to {max_cases} cases with valid MAP data...")
    print()

    cases_attempted = 0
    case_id = 1

    while processed_count < max_cases and case_id <= 6388:
        cases_attempted += 1

        if cases_attempted % 10 == 0:
            print(f"[Attempted {cases_attempted}, Processed {processed_count}]")

        try:
            # Load case with MAP and HR tracks
            data = vitaldb.load_case(
                case_id,
                ['Solar8000/ART_MBP', 'Solar8000/NIBP_MBP', 'Solar8000/HR']
            )

            if data is None or len(data) == 0:
                case_id += 1
                continue

            # Check which columns have data
            # data shape: (n_samples, n_tracks)
            # Tracks: [ART_MBP, NIBP_MBP, HR]

            # Find MAP column (prefer ART, fallback to NIBP)
            map_col = None
            map_name = None

            # Check ART_MBP (column 0)
            if data.shape[1] > 0:
                art_map = data[:, 0]
                valid_art = art_map[(art_map > 0) & (art_map < 200) & (~np.isnan(art_map))]
                if len(valid_art) > 100:  # Need at least 100 valid samples
                    map_col = 0
                    map_name = 'ART_MBP'

            # Check NIBP_MBP (column 1) if ART not available
            if map_col is None and data.shape[1] > 1:
                nibp_map = data[:, 1]
                valid_nibp = nibp_map[(nibp_map > 0) & (nibp_map < 200) & (~np.isnan(nibp_map))]
                if len(valid_nibp) > 100:
                    map_col = 1
                    map_name = 'NIBP_MBP'

            if map_col is None:
                case_id += 1
                continue

            # Get HR (column 2)
            hr_col = 2 if data.shape[1] > 2 else None

            # Extract MAP values
            map_values = data[:, map_col]
            hr_values = data[:, hr_col] if hr_col is not None else None

            # Clean data (replace VitalDB's -8, -9 codes with None)
            map_values = np.where((map_values < 0) | (map_values > 200) | np.isnan(map_values),
                                  np.nan, map_values)

            if hr_values is not None:
                hr_values = np.where((hr_values < 20) | (hr_values > 200) | np.isnan(hr_values),
                                     np.nan, hr_values)

            # Extract features
            case_features = extract_features(
                case_id, map_values, hr_values
            )

            if len(case_features) > 0:
                all_features.extend(case_features)
                case_ioh = sum(1 for f in case_features if f['ioh_label'] == 1)
                ioh_count += case_ioh
                processed_count += 1
                print(f"[{processed_count}/{max_cases}] Case {case_id} ({map_name}): "
                      f"{len(case_features)} samples ({case_ioh} IOH)")

        except Exception as e:
            if 'attempted' not in str(e).lower():
                pass  # Silent skip

        case_id += 1

        if cases_attempted >= 1000:  # Safety limit
            print(f"\nReached attempt limit. Processed {processed_count} valid cases.")
            break

    print()
    print("=" * 70)
    print("DOWNLOAD COMPLETE")
    print("=" * 70)
    print(f"Cases attempted: {cases_attempted}")
    print(f"Cases processed: {processed_count}")
    print(f"Total samples: {len(all_features)}")
    print(f"IOH events: {ioh_count} ({ioh_count/len(all_features)*100:.1f}% if all_features else 0)")
    print()

    # Save dataset
    dataset = {
        'features': all_features,
        'metadata': {
            'total_cases': processed_count,
            'total_samples': len(all_features),
            'ioh_events': ioh_count,
            'ioh_rate': ioh_count / len(all_features) if all_features else 0,
        }
    }

    with open('vitaldb_ioh_dataset.pkl', 'wb') as f:
        pickle.dump(dataset, f)

    print("✓ Dataset saved to: vitaldb_ioh_dataset.pkl")
    print()

    return dataset

def extract_features(case_id, map_values, hr_values):
    """Extract IOH features from MAP time series"""

    features = []

    # Downsample to 10-second intervals (VitalDB is 1/500 Hz)
    # Assume data is already at reasonable sampling rate
    interval = 10  # samples (assuming ~1 Hz data)

    map_ds = map_values[::interval]
    hr_ds = hr_values[::interval] if hr_values is not None else None

    # Need at least 15 minutes of data
    if len(map_ds) < 90:
        return []

    # Calculate baseline (first 60 seconds)
    baseline_map = np.nanmean(map_ds[:6])
    baseline_hr = np.nanmean(hr_ds[:6]) if hr_ds is not None else 75

    if np.isnan(baseline_map):
        return []

    # Extract features for each time point
    prediction_steps = PREDICTION_WINDOW_SEC // 10  # 30 steps

    for i in range(60, len(map_ds) - prediction_steps):
        current_map = map_ds[i]

        if np.isnan(current_map):
            continue

        # MAP trends
        map_5min = map_ds[max(0, i-30)] if i >= 30 else current_map
        map_10min = map_ds[max(0, i-60)] if i >= 60 else current_map

        if np.isnan(map_5min):
            map_5min = current_map
        if np.isnan(map_10min):
            map_10min = current_map

        # Future MAP (for label)
        future_maps = map_ds[i+1:i+1+prediction_steps]
        future_maps_valid = future_maps[~np.isnan(future_maps)]

        if len(future_maps_valid) < 10:  # Need some future data
            continue

        ioh_label = 1 if np.min(future_maps_valid) < IOH_THRESHOLD else 0

        # Create feature dict
        features.append({
            'caseid': case_id,
            'time_idx': i,
            'age': 60,  # Default (not available)
            'sex': 1,
            'bmi': 25,
            'asa': 2,
            'baseline_map': baseline_map,
            'baseline_hr': baseline_hr,
            'current_map': current_map,
            'map_5min': map_5min,
            'map_10min': map_10min,
            'surgery_duration': i,
            'vasopressor': 0,
            'surgery_type': 1,
            'induction_agent': 0,
            'emergency': 0,
            'ioh_label': ioh_label
        })

    return features

if __name__ == "__main__":
    dataset = download_vitaldb_data(max_cases=100)
