#!/usr/bin/env python3
"""
Compare IOH Models: Heuristics vs VitalDB-Trained

Compares performance of:
1. Current heuristics-based model (ioh_models.py)
2. VitalDB-trained Random Forest model

Shows improvement from using real clinical data.
"""

import pickle
import numpy as np
from sklearn.metrics import roc_auc_score, accuracy_score, classification_report
from ioh_models import IOHModelStub, IOHScalerStub

def load_vitaldb_test_data(filepath='vitaldb_ioh_dataset.pkl'):
    """Load VitalDB test data for comparison"""
    with open(filepath, 'rb') as f:
        dataset = pickle.load(f)

    features = dataset['features']

    # Prepare feature matrix
    feature_names = [
        'age', 'sex', 'bmi', 'asa', 'baseline_map', 'baseline_hr',
        'current_map', 'map_5min', 'map_10min', 'surgery_duration',
        'vasopressor', 'surgery_type', 'induction_agent', 'emergency'
    ]

    X = np.array([[f[name] for name in feature_names] for f in features])
    y = np.array([f['ioh_label'] for f in features])

    return X, y

def evaluate_heuristics_model(X, y):
    """Evaluate current heuristics model"""
    print("Evaluating Heuristics Model...")

    model = IOHModelStub()
    scaler = IOHScalerStub()

    X_scaled = scaler.transform(X)
    y_pred_proba = np.array([p[1] for p in model.predict_proba(X_scaled)])
    y_pred = model.predict(X_scaled)

    results = {
        'accuracy': accuracy_score(y, y_pred),
        'roc_auc': roc_auc_score(y, y_pred_proba),
        'y_pred': y_pred,
        'y_pred_proba': y_pred_proba
    }

    print(f"  Accuracy: {results['accuracy']:.4f}")
    print(f"  ROC-AUC:  {results['roc_auc']:.4f}")
    print()

    return results

def evaluate_vitaldb_model(X, y, model_path='vitaldb_ioh_model.pkl',
                          scaler_path='vitaldb_ioh_scaler.pkl'):
    """Evaluate VitalDB-trained model"""
    print("Evaluating VitalDB-Trained Model...")

    with open(model_path, 'rb') as f:
        model = pickle.load(f)

    with open(scaler_path, 'rb') as f:
        scaler = pickle.load(f)

    X_scaled = scaler.transform(X)
    y_pred_proba = model.predict_proba(X_scaled)[:, 1]
    y_pred = model.predict(X_scaled)

    results = {
        'accuracy': accuracy_score(y, y_pred),
        'roc_auc': roc_auc_score(y, y_pred_proba),
        'y_pred': y_pred,
        'y_pred_proba': y_pred_proba
    }

    print(f"  Accuracy: {results['accuracy']:.4f}")
    print(f"  ROC-AUC:  {results['roc_auc']:.4f}")
    print()

    return results

def compare_models():
    """Main comparison function"""
    print("=" * 70)
    print("IOH MODEL COMPARISON")
    print("=" * 70)
    print()

    # Load test data
    print("Loading VitalDB test data...")
    X, y = load_vitaldb_test_data()
    print(f"✓ Loaded {len(X)} test samples")
    print(f"  IOH events: {y.sum()} ({y.sum()/len(y)*100:.1f}%)")
    print()

    # Evaluate heuristics model
    heuristics_results = evaluate_heuristics_model(X, y)

    # Evaluate VitalDB model
    try:
        vitaldb_results = evaluate_vitaldb_model(X, y)
    except FileNotFoundError:
        print("✗ VitalDB model files not found")
        print("  Run train_vitaldb_model.py first")
        return

    # Comparison
    print("=" * 70)
    print("COMPARISON RESULTS")
    print("=" * 70)
    print()

    print(f"{'Metric':<20} {'Heuristics':<15} {'VitalDB':<15} {'Improvement':<15}")
    print("-" * 70)

    # Accuracy
    acc_improvement = vitaldb_results['accuracy'] - heuristics_results['accuracy']
    print(f"{'Accuracy':<20} {heuristics_results['accuracy']:<15.4f} "
          f"{vitaldb_results['accuracy']:<15.4f} {acc_improvement:+.4f}")

    # ROC-AUC
    auc_improvement = vitaldb_results['roc_auc'] - heuristics_results['roc_auc']
    print(f"{'ROC-AUC':<20} {heuristics_results['roc_auc']:<15.4f} "
          f"{vitaldb_results['roc_auc']:<15.4f} {auc_improvement:+.4f}")

    print()
    print("=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    print()

    if vitaldb_results['roc_auc'] > heuristics_results['roc_auc']:
        print("✓ VitalDB-trained model OUTPERFORMS heuristics model")
        print(f"  {auc_improvement:.1%} improvement in ROC-AUC")
        print()
        print("Recommendation: Deploy VitalDB-trained model to production")
    else:
        print("! Heuristics model performs comparably to VitalDB model")
        print("  Consider collecting more training data")

    print()

    # Classification reports
    print("Detailed Classification Report - Heuristics:")
    print(classification_report(y, heuristics_results['y_pred'],
                                target_names=['No IOH', 'IOH'],
                                digits=4))

    print()
    print("Detailed Classification Report - VitalDB:")
    print(classification_report(y, vitaldb_results['y_pred'],
                                target_names=['No IOH', 'IOH'],
                                digits=4))

if __name__ == "__main__":
    compare_models()
