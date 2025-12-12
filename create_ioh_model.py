#!/usr/bin/env python3
"""
Create pre-trained IOH Model Files

This script creates the ioh_model.pkl and ioh_scaler.pkl files
by running the full training process with synthetic clinical data.

Run this script once to generate the model files, then commit them to the repository.
"""

import sys
import os

# Check if dependencies are available
try:
    import numpy as np
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import classification_report, roc_auc_score
    import pickle
except ImportError as e:
    print("=" * 70)
    print("ERROR: Required dependencies not installed")
    print("=" * 70)
    print(f"\nMissing module: {e}")
    print("\nPlease install dependencies first:")
    print("  pip install numpy scikit-learn")
    print("\nOr run:")
    print("  pip install -r requirements.txt")
    print()
    sys.exit(1)

# Set random seed
np.random.seed(42)

def generate_realistic_ioh_data(n_samples=20000):
    """
    Generate realistic synthetic training data for IOH prediction
    based on evidence-based clinical risk factors
    """

    data = []
    labels = []

    for _ in range(n_samples):
        # Demographics
        age = np.random.normal(60, 15)
        age = np.clip(age, 18, 95)

        sex = np.random.choice([0, 1])  # 0=female, 1=male

        bmi = np.random.normal(28, 6)
        bmi = np.clip(bmi, 16, 50)

        asa = np.random.choice([1, 2, 3, 4, 5], p=[0.05, 0.35, 0.40, 0.15, 0.05])

        # Baseline hemodynamics
        baseline_map = np.random.normal(95, 12)
        baseline_map = np.clip(baseline_map, 65, 130)

        baseline_hr = np.random.normal(75, 15)
        baseline_hr = np.clip(baseline_hr, 50, 120)

        # Surgery characteristics
        surgery_type = np.random.choice([0, 1, 2, 3, 4], p=[0.25, 0.35, 0.20, 0.10, 0.10])
        induction_agent = np.random.choice([0, 1, 2], p=[0.65, 0.25, 0.10])
        emergency = np.random.choice([0, 1], p=[0.85, 0.15])
        surgery_duration = np.random.randint(15, 240)

        # Calculate risk score
        risk_score = 0

        if age > 70:
            risk_score += 15
        elif age > 65:
            risk_score += 8

        risk_score += (asa - 1) * 8

        if surgery_type in [3, 4]:  # Cardiac/vascular
            risk_score += 20
        elif surgery_type == 2:  # Major abdominal
            risk_score += 12

        if emergency == 1:
            risk_score += 15

        if induction_agent == 0:  # Propofol
            risk_score += 10
        elif induction_agent == 1:  # Etomidate
            risk_score += 3

        if baseline_map < 80:
            risk_score += 10
        elif baseline_map < 90:
            risk_score += 5

        if bmi > 35:
            risk_score += 8
        elif bmi < 20:
            risk_score += 5

        if surgery_duration > 180:
            risk_score += 10
        elif surgery_duration > 120:
            risk_score += 5

        # Generate MAP values
        if risk_score > 50:
            map_decline = np.random.normal(20, 10)
            map_decline = np.clip(map_decline, 5, 40)
        elif risk_score > 30:
            map_decline = np.random.normal(10, 8)
            map_decline = np.clip(map_decline, 0, 30)
        else:
            map_decline = np.random.normal(5, 5)
            map_decline = np.clip(map_decline, -5, 20)

        current_map = baseline_map - map_decline
        map_10min = current_map + np.random.normal(3, 2)
        map_5min = current_map + np.random.normal(1.5, 1.5)

        current_map = np.clip(current_map, 50, 120)
        map_5min = np.clip(map_5min, 50, 120)
        map_10min = np.clip(map_10min, 50, 120)

        # Vasopressor
        if current_map < 70 or (baseline_map - current_map) > 15:
            vasopressor = np.random.choice([0, 1, 2, 3], p=[0.2, 0.5, 0.2, 0.1])
            if vasopressor > 0:
                risk_score += 10
        else:
            vasopressor = 0

        # Determine IOH outcome
        map_trend = (map_10min - current_map) / 10
        ioh_probability = 0

        if current_map < 70:
            ioh_probability += (70 - current_map) * 3

        if map_trend > 0.5:
            ioh_probability += map_trend * 10

        ioh_probability += risk_score * 0.3
        ioh_probability += np.random.normal(0, 10)
        ioh_probability = np.clip(ioh_probability, 0, 100)

        ioh_outcome = 1 if ioh_probability > 50 else 0

        # Create feature vector
        features = [
            age, sex, bmi, asa, baseline_map, baseline_hr,
            current_map, map_5min, map_10min, surgery_duration,
            vasopressor, surgery_type, induction_agent, emergency
        ]

        data.append(features)
        labels.append(ioh_outcome)

    return np.array(data), np.array(labels)

def main():
    """Train and save IOH prediction model"""

    print("=" * 70)
    print("IOH PREDICTION MODEL - TRAINING")
    print("=" * 70)
    print()

    # Generate training data
    print("Generating synthetic clinical training data...")
    X, y = generate_realistic_ioh_data(n_samples=20000)

    print(f"✓ Dataset: {len(X)} samples")
    print(f"✓ IOH rate: {y.sum() / len(y) * 100:.1f}%")
    print()

    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    print(f"Training set: {len(X_train)} samples")
    print(f"Test set: {len(X_test)} samples")
    print()

    # Scale features
    print("Scaling features...")
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    print()

    # Train model
    print("Training RandomForest classifier...")
    print("  - 200 trees")
    print("  - Max depth: 12")
    print("  - Balanced class weights")
    print()

    model = RandomForestClassifier(
        n_estimators=200,
        max_depth=12,
        min_samples_split=20,
        min_samples_leaf=10,
        max_features='sqrt',
        random_state=42,
        n_jobs=-1,
        class_weight='balanced'
    )

    model.fit(X_train_scaled, y_train)
    print("✓ Training complete!")
    print()

    # Evaluate
    print("=" * 70)
    print("MODEL EVALUATION")
    print("=" * 70)
    print()

    train_score = model.score(X_train_scaled, y_train)
    test_score = model.score(X_test_scaled, y_test)

    print(f"Training accuracy: {train_score:.3f}")
    print(f"Test accuracy: {test_score:.3f}")

    y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
    roc_auc = roc_auc_score(y_test, y_pred_proba)
    print(f"ROC-AUC: {roc_auc:.3f}")
    print()

    # Classification report
    y_pred = model.predict(X_test_scaled)
    print("Classification Report:")
    print(classification_report(y_test, y_pred, target_names=['No IOH', 'IOH']))

    # Feature importance
    feature_names = [
        'age', 'sex', 'bmi', 'asa', 'baseline_map', 'baseline_hr',
        'current_map', 'map_5min', 'map_10min', 'surgery_duration',
        'vasopressor', 'surgery_type', 'induction_agent', 'emergency'
    ]

    importances = model.feature_importances_
    feature_importance = sorted(
        zip(feature_names, importances),
        key=lambda x: x[1],
        reverse=True
    )

    print()
    print("Top 10 Feature Importances:")
    print("-" * 40)
    for feat, imp in feature_importance[:10]:
        print(f"{feat:20s}: {imp:.4f}")
    print()

    # Save model and scaler
    print("=" * 70)
    print("SAVING MODEL FILES")
    print("=" * 70)
    print()

    with open('ioh_model.pkl', 'wb') as f:
        pickle.dump(model, f)
    print("✓ Model saved: ioh_model.pkl")

    with open('ioh_scaler.pkl', 'wb') as f:
        pickle.dump(scaler, f)
    print("✓ Scaler saved: ioh_scaler.pkl")
    print()

    # Verify files
    model_size = os.path.getsize('ioh_model.pkl') / 1024
    scaler_size = os.path.getsize('ioh_scaler.pkl') / 1024

    print(f"Model file size: {model_size:.1f} KB")
    print(f"Scaler file size: {scaler_size:.1f} KB")
    print()

    print("=" * 70)
    print("SUCCESS!")
    print("=" * 70)
    print()
    print("The IOH prediction model has been trained and saved.")
    print()
    print("Next steps:")
    print("  1. Test the model with the IOH calculator")
    print("  2. Commit both .pkl files to the repository")
    print("  3. Deploy to production")
    print()

if __name__ == "__main__":
    main()
