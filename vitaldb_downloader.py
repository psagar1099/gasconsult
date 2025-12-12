#!/usr/bin/env python3
"""
VitalDB Data Downloader for IOH Prediction

Downloads VitalDB cases and prepares them for IOH prediction model training.
Uses VitalDB REST API (no package installation required).

VitalDB API Documentation: https://vitaldb.net/docs/
Dataset: 6,388 surgical cases with high-resolution vital signs
"""

import requests
import json
import gzip
import io
import csv
import os
import time
from typing import List, Dict, Optional
import pickle

# VitalDB API endpoints
VITALDB_API_BASE = "https://api.vitaldb.net"
VITALDB_CASES_ENDPOINT = f"{VITALDB_API_BASE}/cases"
VITALDB_TRKS_ENDPOINT = f"{VITALDB_API_BASE}/trks"
VITALDB_VALS_ENDPOINT = f"{VITALDB_API_BASE}/vals"

# IOH definition: MAP < 65 mmHg
IOH_THRESHOLD = 65

# Prediction window: predict IOH in next 5 minutes
PREDICTION_WINDOW_SEC = 300

def download_case_list() -> List[Dict]:
    """
    Download list of all VitalDB cases with clinical information

    Returns:
        List of dictionaries containing case metadata
    """
    print("Downloading VitalDB case list...")

    try:
        response = requests.get(VITALDB_CASES_ENDPOINT, timeout=30)
        response.raise_for_status()

        # Response is gzip-compressed CSV
        with gzip.GzipFile(fileobj=io.BytesIO(response.content)) as f:
            csv_data = f.read().decode('utf-8')

        # Parse CSV
        reader = csv.DictReader(io.StringIO(csv_data))
        cases = list(reader)

        print(f"✓ Downloaded {len(cases)} cases")
        return cases

    except Exception as e:
        print(f"✗ Error downloading case list: {e}")
        return []

def get_available_tracks(caseid: str) -> List[str]:
    """
    Get list of available tracks (vital sign parameters) for a case

    Args:
        caseid: VitalDB case ID

    Returns:
        List of track names
    """
    try:
        url = f"{VITALDB_TRKS_ENDPOINT}/{caseid}"
        response = requests.get(url, timeout=10)
        response.raise_for_status()

        # Parse gzip CSV
        with gzip.GzipFile(fileobj=io.BytesIO(response.content)) as f:
            csv_data = f.read().decode('utf-8')

        reader = csv.DictReader(io.StringIO(csv_data))
        tracks = [row['name'] for row in reader]

        return tracks

    except Exception as e:
        print(f"  Warning: Could not get tracks for case {caseid}: {e}")
        return []

def download_track_data(caseid: str, track: str, interval: int = 10) -> Optional[List[float]]:
    """
    Download time-series data for a specific track

    Args:
        caseid: VitalDB case ID
        track: Track name (e.g., "Solar8000/ART_MBP")
        interval: Sampling interval in seconds (default: 10s)

    Returns:
        List of values, or None if error
    """
    try:
        url = f"{VITALDB_VALS_ENDPOINT}/{caseid}/{track}/{interval}"
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        # Parse gzip CSV
        with gzip.GzipFile(fileobj=io.BytesIO(response.content)) as f:
            csv_data = f.read().decode('utf-8')

        # Values are comma-separated
        values = [float(v) if v and v != '' else None
                 for v in csv_data.strip().split(',')]

        return values

    except Exception as e:
        print(f"  Warning: Could not download {track} for case {caseid}: {e}")
        return None

def extract_ioh_features(caseid: str, case_info: Dict,
                         map_data: List[float],
                         hr_data: List[float],
                         interval: int = 10) -> List[Dict]:
    """
    Extract features and labels for IOH prediction from time-series data

    For each time point, we create a feature vector and label (IOH in next 5 min)

    Args:
        caseid: Case ID
        case_info: Clinical information dictionary
        map_data: Mean arterial pressure time series (every 10s)
        hr_data: Heart rate time series (every 10s)
        interval: Sampling interval in seconds

    Returns:
        List of feature dictionaries with labels
    """
    features_list = []

    # Parse clinical info
    age = int(case_info.get('age', 0)) if case_info.get('age') else None
    sex = 1 if case_info.get('sex') == 'M' else 0
    height = float(case_info.get('height', 0)) if case_info.get('height') else None
    weight = float(case_info.get('weight', 0)) if case_info.get('weight') else None
    bmi = weight / ((height / 100) ** 2) if weight and height else None
    asa = int(case_info.get('asa', 0)) if case_info.get('asa') else None

    # Surgical type mapping (simplified)
    optype = case_info.get('optype', '').lower()
    if 'cardiac' in optype or 'cabg' in optype:
        surgery_type = 3  # cardiac
    elif 'vascular' in optype or 'aorta' in optype:
        surgery_type = 4  # vascular
    elif 'abdom' in optype or 'laparotomy' in optype:
        surgery_type = 2  # major abdominal
    elif 'minor' in optype or 'superficial' in optype:
        surgery_type = 0  # minor
    else:
        surgery_type = 1  # moderate

    emergency = 1 if case_info.get('emop') == 'Y' else 0

    # Need at least MAP data
    if not map_data or len(map_data) < 90:  # At least 15 minutes of data
        return features_list

    # Iterate through time points (start at 10 minutes to have baseline)
    prediction_window_steps = PREDICTION_WINDOW_SEC // interval  # 30 steps for 5 min
    baseline_steps = 60 // interval  # First 60 seconds for baseline

    for i in range(60, len(map_data) - prediction_window_steps):
        # Skip if MAP is missing at current time or in future window
        if map_data[i] is None:
            continue

        current_map = map_data[i]

        # Baseline MAP (first 60 seconds)
        baseline_maps = [m for m in map_data[:baseline_steps] if m is not None]
        if not baseline_maps:
            continue
        baseline_map = sum(baseline_maps) / len(baseline_maps)

        # MAP 5 and 10 minutes ago
        idx_5min = max(0, i - (300 // interval))
        idx_10min = max(0, i - (600 // interval))

        map_5min = map_data[idx_5min] if map_data[idx_5min] is not None else current_map
        map_10min = map_data[idx_10min] if map_data[idx_10min] is not None else current_map

        # Heart rate (current and baseline)
        if hr_data and i < len(hr_data) and hr_data[i] is not None:
            current_hr = hr_data[i]
            baseline_hrs = [h for h in hr_data[:baseline_steps] if h is not None]
            baseline_hr = sum(baseline_hrs) / len(baseline_hrs) if baseline_hrs else current_hr
        else:
            current_hr = 75  # Default
            baseline_hr = 75

        # Surgery duration (in minutes)
        surgery_duration = (i * interval) // 60

        # Induction agent (not available in VitalDB, use propofol as default)
        induction_agent = 0  # propofol

        # Vasopressor use (simplified: if MAP recovering after drop, assume vasopressor)
        vasopressor = 0
        if i > 10:
            recent_maps = [m for m in map_data[i-10:i] if m is not None]
            if recent_maps:
                min_recent_map = min(recent_maps)
                if min_recent_map < 70 and current_map > min_recent_map + 5:
                    vasopressor = 1  # Likely on vasopressor

        # Label: IOH in next 5 minutes?
        future_maps = [m for m in map_data[i+1:i+1+prediction_window_steps] if m is not None]
        if not future_maps:
            continue

        ioh_label = 1 if min(future_maps) < IOH_THRESHOLD else 0

        # Create feature vector
        features = {
            'caseid': caseid,
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
            'induction_agent': induction_agent,
            'emergency': emergency,
            'ioh_label': ioh_label
        }

        # Only include if all critical features are available
        if all([age, bmi, asa, baseline_map, current_map]):
            features_list.append(features)

    return features_list

def download_and_process_vitaldb(max_cases: int = 1000,
                                  output_file: str = 'vitaldb_ioh_dataset.pkl') -> Dict:
    """
    Download VitalDB data and process for IOH prediction

    Args:
        max_cases: Maximum number of cases to download
        output_file: Pickle file to save processed data

    Returns:
        Dictionary with processed features and metadata
    """
    print("=" * 70)
    print("VITALDB IOH DATASET CREATION")
    print("=" * 70)
    print()

    # Download case list
    all_cases = download_case_list()
    if not all_cases:
        print("✗ Failed to download case list")
        return None

    # Limit cases
    cases_to_process = all_cases[:max_cases]
    print(f"Processing {len(cases_to_process)} cases...")
    print()

    all_features = []
    ioh_count = 0
    processed_count = 0

    for idx, case_info in enumerate(cases_to_process):
        caseid = case_info.get('caseid')
        if not caseid:
            continue

        print(f"[{idx+1}/{len(cases_to_process)}] Processing case {caseid}...", end=' ')

        # Get available tracks
        tracks = get_available_tracks(caseid)

        # Look for MAP track (various names in VitalDB)
        map_track = None
        for track in tracks:
            if 'ART_MBP' in track or 'NIBP_MBP' in track or 'ART' in track.upper():
                map_track = track
                break

        if not map_track:
            print("(no MAP data)")
            continue

        # Look for HR track
        hr_track = None
        for track in tracks:
            if 'HR' in track.upper() and 'ECG' in track.upper():
                hr_track = track
                break

        # Download MAP data (10 second intervals)
        map_data = download_track_data(caseid, map_track, interval=10)
        if not map_data:
            print("(MAP download failed)")
            continue

        # Download HR data
        hr_data = download_track_data(caseid, hr_track, interval=10) if hr_track else None

        # Extract features
        case_features = extract_ioh_features(caseid, case_info, map_data, hr_data)

        if case_features:
            all_features.extend(case_features)
            case_ioh = sum(1 for f in case_features if f['ioh_label'] == 1)
            ioh_count += case_ioh
            processed_count += 1
            print(f"✓ {len(case_features)} samples ({case_ioh} IOH)")
        else:
            print("(insufficient data)")

        # Rate limiting
        time.sleep(0.5)

        # Save progress every 50 cases
        if (idx + 1) % 50 == 0:
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

    # Save to pickle
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

    with open(output_file, 'wb') as f:
        pickle.dump(dataset, f)

    print(f"✓ Dataset saved to: {output_file}")
    print()

    return dataset

if __name__ == "__main__":
    # Download and process VitalDB data
    # Start with 100 cases for testing (full dataset = 6,388 cases)
    dataset = download_and_process_vitaldb(max_cases=100,
                                          output_file='vitaldb_ioh_dataset.pkl')
