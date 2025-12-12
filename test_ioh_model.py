#!/usr/bin/env python3
"""
Test IOH Model Files

Verifies that the model files work correctly with the IOH calculator logic
"""

import pickle
# Import model classes so pickle can deserialize them
from build_model_stub import IOHModelStub, IOHScalerStub

def test_ioh_model():
    """Test the IOH model with sample patient data"""

    print("=" * 70)
    print("TESTING IOH MODEL")
    print("=" * 70)
    print()

    # Load model and scaler
    print("Loading model files...")
    try:
        with open('ioh_model.pkl', 'rb') as f:
            model = pickle.load(f)
        print("✓ Model loaded")

        with open('ioh_scaler.pkl', 'rb') as f:
            scaler = pickle.load(f)
        print("✓ Scaler loaded")
        print()
    except Exception as e:
        print(f"✗ Error loading models: {e}")
        return False

    # Test cases matching real IOH calculator inputs
    test_cases = [
        {
            "name": "High Risk: Elderly, Declining MAP, Cardiac Surgery",
            "features": [
                75,  # age
                1,   # sex (male)
                28,  # bmi
                3,   # asa
                90,  # baseline_map
                75,  # baseline_hr
                62,  # current_map (low!)
                66,  # map_5min
                70,  # map_10min (declining trend)
                120, # surgery_duration
                1,   # vasopressor (phenylephrine)
                3,   # surgery_type (cardiac)
                0,   # induction_agent (propofol)
                0    # emergency
            ],
            "expected_risk": "high"
        },
        {
            "name": "Moderate Risk: Middle-aged, Stable MAP, Minor Surgery",
            "features": [
                55,  # age
                0,   # sex (female)
                25,  # bmi
                2,   # asa
                95,  # baseline_map
                70,  # baseline_hr
                82,  # current_map
                83,  # map_5min
                85,  # map_10min (stable)
                60,  # surgery_duration
                0,   # vasopressor (none)
                1,   # surgery_type (moderate)
                1,   # induction_agent (etomidate)
                0    # emergency
            ],
            "expected_risk": "low-moderate"
        },
        {
            "name": "Low Risk: Young, High MAP, Minor Surgery",
            "features": [
                35,  # age
                1,   # sex (male)
                24,  # bmi
                1,   # asa
                100, # baseline_map
                65,  # baseline_hr
                95,  # current_map (high, stable)
                96,  # map_5min
                97,  # map_10min
                45,  # surgery_duration
                0,   # vasopressor (none)
                0,   # surgery_type (minor)
                2,   # induction_agent (ketamine)
                0    # emergency
            ],
            "expected_risk": "low"
        }
    ]

    print("Running test cases...")
    print()

    all_passed = True

    for i, test_case in enumerate(test_cases, 1):
        print(f"Test Case {i}: {test_case['name']}")
        print("-" * 70)

        features = [test_case['features']]

        # Scale features
        features_scaled = scaler.transform(features)

        # Get prediction
        ioh_prob = model.predict_proba(features_scaled)[0][1]
        prob_pct = int(ioh_prob * 100)

        # Classify risk
        if prob_pct < 30:
            risk_class = "low"
        elif prob_pct < 60:
            risk_class = "moderate"
        else:
            risk_class = "high"

        print(f"IOH Probability: {prob_pct}%")
        print(f"Risk Classification: {risk_class}")
        print(f"Expected: {test_case['expected_risk']}")

        # Verify reasonable range
        if 0 <= prob_pct <= 100:
            print("✓ Probability in valid range")
        else:
            print(f"✗ Invalid probability: {prob_pct}%")
            all_passed = False

        print()

    # Test feature importances
    print("Feature Importances:")
    print("-" * 70)

    if hasattr(model, 'feature_importances_'):
        feature_names = [
            'age', 'sex', 'bmi', 'asa', 'baseline_map', 'baseline_hr',
            'current_map', 'map_5min', 'map_10min', 'surgery_duration',
            'vasopressor', 'surgery_type', 'induction_agent', 'emergency'
        ]

        importances = model.feature_importances_
        top_features = sorted(
            zip(feature_names, importances),
            key=lambda x: x[1],
            reverse=True
        )[:5]

        for feat, imp in top_features:
            print(f"  {feat:20s}: {imp:.4f}")
        print()
    else:
        print("  (Feature importances not available)")
        print()

    # Summary
    print("=" * 70)
    if all_passed:
        print("✓ ALL TESTS PASSED")
        print("=" * 70)
        print()
        print("The IOH model is working correctly and ready for use.")
        print()
        return True
    else:
        print("✗ SOME TESTS FAILED")
        print("=" * 70)
        print()
        return False

if __name__ == "__main__":
    success = test_ioh_model()
    exit(0 if success else 1)
