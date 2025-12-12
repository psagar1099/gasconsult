#!/usr/bin/env python3
"""
Build IOH Model Stub Files

This creates functional pickle files for the IOH calculator
that can be used until proper models are trained.

The stub provides realistic probability estimates based on clinical heuristics.
"""

import pickle
from ioh_models import IOHModelStub, IOHScalerStub

# Class definitions are now in ioh_models.py module
# This ensures they can be imported when unpickling

def create_stub_models():
    """Create and save stub model files"""

    print("=" * 70)
    print("CREATING IOH MODEL STUB FILES")
    print("=" * 70)
    print()
    print("This creates functional pickle files for the IOH calculator")
    print("based on clinical heuristics.")
    print()

    # Create model
    model = IOHModelStub()
    print("✓ Model object created")

    # Create scaler
    scaler = IOHScalerStub()
    print("✓ Scaler object created")
    print()

    # Save to pickle files
    print("Saving files...")

    with open('ioh_model.pkl', 'wb') as f:
        pickle.dump(model, f)
    print("✓ ioh_model.pkl saved")

    with open('ioh_scaler.pkl', 'wb') as f:
        pickle.dump(scaler, f)
    print("✓ ioh_scaler.pkl saved")
    print()

    # Test loading
    print("Verifying files...")

    with open('ioh_model.pkl', 'rb') as f:
        test_model = pickle.load(f)

    with open('ioh_scaler.pkl', 'rb') as f:
        test_scaler = pickle.load(f)

    print("✓ Files load successfully")
    print()

    # Test prediction
    print("Testing prediction...")

    # Test case: 70yo, male, BMI 28, ASA 3, baseline 90/75, current 65, declining trend
    test_features = [[70, 1, 28, 3, 90, 75, 65, 68, 72, 120, 1, 2, 0, 0]]
    scaled = test_scaler.transform(test_features)
    prob = test_model.predict_proba(scaled)[0][1]

    print(f"✓ Test prediction: {prob*100:.1f}% IOH risk")
    print()

    print("=" * 70)
    print("SUCCESS!")
    print("=" * 70)
    print()
    print("Model stub files created successfully.")
    print("These provide clinically-informed predictions for the IOH calculator.")
    print()
    print("NOTE: For production use, train a proper ML model using:")
    print("  python create_ioh_model.py")
    print()

if __name__ == "__main__":
    create_stub_models()
