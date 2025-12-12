#!/usr/bin/env python3
"""
VitalDB Data Downloader - Using Official VitalDB Library

This version uses the official vitaldb Python package which handles
authentication and access properly.

Installation: pip install vitaldb
"""

import pickle
import numpy as np
from typing import List, Dict

# Check if vitaldb is available
try:
    import vitaldb
    VITALDB_AVAILABLE = True
except ImportError:
    VITALDB_AVAILABLE = False
    print("=" * 70)
    print("ERROR: VitalDB library not installed")
    print("=" * 70)
    print()
    print("Please install the VitalDB library first:")
    print("  pip install vitaldb")
    print()
    print("Then run this script again.")
    print()
    exit(1)

# IOH definition: MAP < 65 mmHg
IOH_THRESHOLD = 65

# Prediction window: predict IOH in next 5 minutes
PREDICTION_WINDOW_SEC = 300

def download_vitaldb_data(max_cases: int = 100) -> Dict:
    """
    Download VitalDB data using official library

    Args:
        max_cases: Maximum number of cases to process

    Returns:
        Dictionary with features and metadata
    """
    print("=" * 70)
    print("VITALDB IOH DATASET CREATION (Official Library)")
    print("=" * 70)
    print()

    # Get list of all cases
    print("Loading VitalDB case list...")
    vf = vitaldb.VitalFile()
    all_cases = vf.get_cases()

    print(f"✓ Found {len(all_cases)} total cases")
    print(f"Processing first {max_cases} cases...")
    print()

    all_features = []
    ioh_count = 0
    processed_count = 0

    for idx, case_id in enumerate(all_cases[:max_cases]):
        print(f"[{idx+1}/{max_cases}] Processing case {case_id}...", end=' ')

        try:
            # Load case data
            vf = vitaldb.VitalFile(case_id)

            # Get clinical info
            clinical_info = vf.get_clinical_info()

            # Extract demographics
            age = clinical_info.get('age', None)
            sex = 1 if clinical_info.get('sex') == 'M' else 0
            height = clinical_info.get('height', None)
            weight = clinical_info.get('weight', None)
            bmi = weight / ((height / 100) ** 2) if weight and height else None
            asa = clinical_info.get('asa', None)

            # Get MAP data (try different track names)
            map_track = None
            for track_name in ['SNUADC/ART_MBP', 'Solar8000/ART_MBP', 'SNUADC/NIBP_MBP']:
                if track_name in vf.get_track_names():
                    map_track = track_name
                    break

            if not map_track:
                print("(no MAP data)")
                continue

            # Load MAP data at 10-second intervals
            map_data = vf.get_samples(map_track, interval=10)

            if map_data is None or len(map_data) < 90:
                print("(insufficient MAP data)")
                continue

            # Try to get HR data
            hr_track = None
            for track_name in ['SNUADC/ECG_HR', 'Solar8000/HR']:
                if track_name in vf.get_track_names():
                    hr_track = track_name
                    break

            hr_data = vf.get_samples(hr_track, interval=10) if hr_track else None

            # Process this case
            case_features = extract_features_from_case(
                case_id, clinical_info, map_data, hr_data
            )

            if case_features:
                all_features.extend(case_features)
                case_ioh = sum(1 for f in case_features if f['ioh_label'] == 1)
                ioh_count += case_ioh
                processed_count += 1
                print(f"✓ {len(case_features)} samples ({case_ioh} IOH)")
            else:
                print("(processing failed)")

        except Exception as e:
            print(f"✗ Error: {str(e)[:50]}")
            continue

        # Progress update every 10 cases
        if (idx + 1) % 10 == 0:
            print()
            print(f"Progress: {processed_count} cases, {len(all_features)} samples, {ioh_count} IOH events")
            print()

    print()
    print("=" * 70)
    print("DOWNLOAD COMPLETE")
    print("=" * 70)
    print(f"Processed cases: {processed_count}")
    print(f"Total samples: {len(all_features)}")
    print(f"IOH events: {ioh_count} ({ioh_count/len(all_features)*100:.1f}%)")
    print()

    # Create dataset
    dataset = {
        'features': all_features,
        'metadata': {
            'total_cases': processed_count,
            'total_samples': len(all_features),
            'ioh_events': ioh_count,
            'ioh_rate': ioh_count / len(all_features) if all_features else 0,
            'prediction_window_sec': PREDICTION_WINDOW_SEC,
            'ioh_threshold': IOH_THRESHOLD
        }
    }

    # Save to pickle
    output_file = 'vitaldb_ioh_dataset.pkl'
    with open(output_file, 'wb') as f:
        pickle.dump(dataset, f)

    print(f"✓ Dataset saved to: {output_file}")
    print()

    return dataset

def extract_features_from_case(case_id, clinical_info, map_data, hr_data):
    """Extract IOH prediction features from a single case"""

    features_list = []

    # Parse clinical info
    age = int(clinical_info.get('age', 0)) if clinical_info.get('age') else None
    sex = 1 if clinical_info.get('sex') == 'M' else 0
    height = clinical_info.get('height', None)
    weight = clinical_info.get('weight', None)
    bmi = weight / ((height / 100) ** 2) if weight and height else None
    asa = int(clinical_info.get('asa', 0)) if clinical_info.get('asa') else None

    # Surgery type (simplified)
    optype = clinical_info.get('optype', '').lower()
    if 'cardiac' in optype or 'cabg' in optype:
        surgery_type = 3
    elif 'vascular' in optype:
        surgery_type = 4
    elif 'abdom' in optype or 'laparotomy' in optype:
        surgery_type = 2
    else:
        surgery_type = 1

    emergency = 1 if clinical_info.get('emop') == 'Y' else 0

    # Calculate baseline MAP (first 60 seconds)
    baseline_maps = [m for m in map_data[:6] if m is not None]
    if not baseline_maps:
        return []
    baseline_map = sum(baseline_maps) / len(baseline_maps)

    # Baseline HR
    if hr_data:
        baseline_hrs = [h for h in hr_data[:6] if h is not None]
        baseline_hr = sum(baseline_hrs) / len(baseline_hrs) if baseline_hrs else 75
    else:
        baseline_hr = 75

    # Process each time point
    prediction_window_steps = PREDICTION_WINDOW_SEC // 10  # 30 steps for 5 min

    for i in range(60, len(map_data) - prediction_window_steps):
        if map_data[i] is None:
            continue

        current_map = map_data[i]

        # MAP trends
        idx_5min = max(0, i - 30)
        idx_10min = max(0, i - 60)

        map_5min = map_data[idx_5min] if map_data[idx_5min] is not None else current_map
        map_10min = map_data[idx_10min] if map_data[idx_10min] is not None else current_map

        # Surgery duration
        surgery_duration = (i * 10) // 60  # in minutes

        # Vasopressor detection (simplified)
        vasopressor = 0
        if i > 10:
            recent_maps = [m for m in map_data[i-10:i] if m is not None]
            if recent_maps and min(recent_maps) < 70 and current_map > min(recent_maps) + 5:
                vasopressor = 1

        # IOH label (MAP < 65 in next 5 minutes)
        future_maps = [m for m in map_data[i+1:i+1+prediction_window_steps] if m is not None]
        if not future_maps:
            continue

        ioh_label = 1 if min(future_maps) < IOH_THRESHOLD else 0

        # Create feature dict
        features = {
            'caseid': case_id,
            'time_idx': i,
            'age': age,
            'sex': sex,
            'bmi': bmi,
            'asa': asa,
            'baseline_map': baseline_map,
            'baseline_hr': baseline_hr,
            'current_map': current_map,
            'map_5min': map_5min,
            'map_10min': map_10min,
            'surgery_duration': surgery_duration,
            'vasopressor': vasopressor,
            'surgery_type': surgery_type,
            'induction_agent': 0,  # Not available
            'emergency': emergency,
            'ioh_label': ioh_label
        }

        # Only include if all critical features available
        if all([age, bmi, asa, baseline_map, current_map]):
            features_list.append(features)

    return features_list

if __name__ == "__main__":
    # Download 100 cases for testing
    dataset = download_vitaldb_data(max_cases=100)
