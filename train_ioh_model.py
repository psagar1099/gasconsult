"""
Train Machine Learning Model for Intraoperative Hypotension Prediction
Generates synthetic training data based on clinical knowledge and trains a RandomForestClassifier
"""

import numpy as np
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

# Set random seed for reproducibility
np.random.seed(42)

def generate_synthetic_ioh_data(n_samples=10000):
    """
    Generate synthetic training data based on clinical IOH risk factors

    Features:
    - age (20-90 years)
    - sex (0=female, 1=male)
    - bmi (18-45)
    - asa (1-4)
    - baseline_map (65-110 mmHg)
    - baseline_hr (50-120 bpm)
    - current_map (45-100 mmHg)
    - map_5min (45-105 mmHg)
    - map_10min (50-110 mmHg)
    - surgery_duration (30-480 minutes)
    - vasopressor (0=none, 1=phenylephrine, 2=ephedrine, 3=norepinephrine)
    - surgery_type (0=minor, 1=moderate, 2=major_abdominal, 3=cardiac, 4=vascular)
    - induction_agent (0=propofol, 1=etomidate, 2=ketamine)
    - emergency (0=no, 1=yes)

    Target: IOH in next 5 minutes (0=no, 1=yes)
    """

    data = []
    labels = []

    for _ in range(n_samples):
        # Generate patient demographics
        age = np.random.randint(20, 91)
        sex = np.random.randint(0, 2)
        bmi = np.random.normal(28, 6)
        bmi = max(18, min(45, bmi))
        asa = np.random.choice([1, 2, 3, 4], p=[0.1, 0.4, 0.4, 0.1])

        # Generate baseline vitals
        baseline_map = np.random.normal(85, 12)
        baseline_map = max(65, min(110, baseline_map))
        baseline_hr = np.random.normal(75, 15)
        baseline_hr = max(50, min(120, baseline_hr))

        # Generate MAP trend (key predictor)
        map_10min = np.random.normal(baseline_map - 5, 10)
        map_5min = map_10min + np.random.normal(-3, 5)
        current_map = map_5min + np.random.normal(-3, 5)

        # Ensure realistic ranges
        map_10min = max(50, min(110, map_10min))
        map_5min = max(45, min(105, map_5min))
        current_map = max(45, min(100, current_map))

        # Surgical factors
        surgery_duration = np.random.randint(30, 481)
        vasopressor = np.random.choice([0, 1, 2, 3], p=[0.5, 0.3, 0.15, 0.05])
        surgery_type = np.random.choice([0, 1, 2, 3, 4], p=[0.3, 0.3, 0.2, 0.1, 0.1])
        induction_agent = np.random.choice([0, 1, 2], p=[0.6, 0.3, 0.1])
        emergency = np.random.choice([0, 1], p=[0.8, 0.2])

        # Calculate IOH risk based on clinical evidence
        risk_score = 0

        # MAP trend (most important)
        map_trend = (current_map - map_5min) + (current_map - map_10min) / 2
        if map_trend < -10:
            risk_score += 35
        elif map_trend < -5:
            risk_score += 20
        elif map_trend < 0:
            risk_score += 10

        # Current MAP level
        if current_map < 65:
            risk_score += 30
        elif current_map < 70:
            risk_score += 20
        elif current_map < 75:
            risk_score += 10

        # MAP drop from baseline
        map_drop_pct = ((baseline_map - current_map) / baseline_map) * 100
        if map_drop_pct > 30:
            risk_score += 25
        elif map_drop_pct > 20:
            risk_score += 15
        elif map_drop_pct > 10:
            risk_score += 8

        # Age risk
        if age > 75:
            risk_score += 15
        elif age > 65:
            risk_score += 10
        elif age > 55:
            risk_score += 5

        # ASA class
        if asa >= 4:
            risk_score += 20
        elif asa == 3:
            risk_score += 12
        elif asa == 2:
            risk_score += 5

        # Surgery type
        if surgery_type in [3, 4]:  # cardiac, vascular
            risk_score += 12
        elif surgery_type == 2:  # major abdominal
            risk_score += 8

        # Vasopressor use (indicates instability)
        if vasopressor > 0:
            risk_score += 10

        # Induction agent (propofol)
        if induction_agent == 0:
            risk_score += 8

        # Emergency surgery
        if emergency == 1:
            risk_score += 10

        # BMI extremes
        if bmi < 20 or bmi > 35:
            risk_score += 5

        # Surgery duration
        if surgery_duration > 240:
            risk_score += 8
        elif surgery_duration > 180:
            risk_score += 5

        # Convert risk score to probability (with some randomness)
        base_prob = min(0.95, risk_score / 120)
        random_factor = np.random.normal(0, 0.1)
        ioh_prob = max(0, min(1, base_prob + random_factor))

        # Generate binary label
        ioh_label = 1 if np.random.random() < ioh_prob else 0

        # Store features
        features = [
            age, sex, bmi, asa, baseline_map, baseline_hr,
            current_map, map_5min, map_10min, surgery_duration,
            vasopressor, surgery_type, induction_agent, emergency
        ]

        data.append(features)
        labels.append(ioh_label)

    return np.array(data), np.array(labels)

def train_model():
    """Train and save the IOH prediction model"""

    print("Generating synthetic training data...")
    X, y = generate_synthetic_ioh_data(n_samples=10000)

    print(f"Dataset: {X.shape[0]} samples, {X.shape[1]} features")
    print(f"IOH prevalence: {y.mean()*100:.1f}%")

    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    # Scale features
    print("\nScaling features...")
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Train Random Forest model
    print("\nTraining RandomForestClassifier...")
    model = RandomForestClassifier(
        n_estimators=200,
        max_depth=15,
        min_samples_split=10,
        min_samples_leaf=5,
        random_state=42,
        n_jobs=-1,
        class_weight='balanced'
    )

    model.fit(X_train_scaled, y_train)

    # Evaluate
    train_score = model.score(X_train_scaled, y_train)
    test_score = model.score(X_test_scaled, y_test)

    print(f"\nTraining accuracy: {train_score*100:.2f}%")
    print(f"Test accuracy: {test_score*100:.2f}%")

    # Feature importance
    feature_names = [
        'age', 'sex', 'bmi', 'asa', 'baseline_map', 'baseline_hr',
        'current_map', 'map_5min', 'map_10min', 'surgery_duration',
        'vasopressor', 'surgery_type', 'induction_agent', 'emergency'
    ]

    importances = model.feature_importances_
    indices = np.argsort(importances)[::-1]

    print("\nTop 5 Feature Importances:")
    for i in range(5):
        idx = indices[i]
        print(f"  {feature_names[idx]}: {importances[idx]:.4f}")

    # Save model and scaler
    print("\nSaving model and scaler...")
    with open('ioh_model.pkl', 'wb') as f:
        pickle.dump(model, f)

    with open('ioh_scaler.pkl', 'wb') as f:
        pickle.dump(scaler, f)

    print("\n✓ Model saved as 'ioh_model.pkl'")
    print("✓ Scaler saved as 'ioh_scaler.pkl'")
    print("\nModel training complete!")

if __name__ == "__main__":
    train_model()
