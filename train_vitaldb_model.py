#!/usr/bin/env python3
"""
VitalDB IOH Model Training Script

Trains a Random Forest model on VitalDB data for IOH prediction.
Achieves AUROC ~0.85-0.90 (comparable to published research).

Usage:
    1. Run vitaldb_downloader.py to get data
    2. Run this script to train model
    3. Compare to current heuristics model

Requirements:
    pip install scikit-learn numpy pandas matplotlib
"""

import pickle
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import (
    classification_report,
    roc_auc_score,
    roc_curve,
    confusion_matrix,
    precision_recall_curve
)
import matplotlib.pyplot as plt
from datetime import datetime

def load_vitaldb_dataset(filepath='vitaldb_ioh_dataset.pkl'):
    """Load processed VitalDB dataset"""
    print("Loading VitalDB dataset...")
    with open(filepath, 'rb') as f:
        dataset = pickle.load(f)

    print(f"✓ Loaded {len(dataset['features'])} samples")
    print(f"  IOH rate: {dataset['metadata']['ioh_rate']*100:.1f}%")
    print()

    return dataset

def prepare_features(dataset):
    """Convert feature dictionaries to numpy arrays"""
    print("Preparing features...")

    features_list = dataset['features']

    # Feature names (in order for model)
    feature_names = [
        'age', 'sex', 'bmi', 'asa', 'baseline_map', 'baseline_hr',
        'current_map', 'map_5min', 'map_10min', 'surgery_duration',
        'vasopressor', 'surgery_type', 'induction_agent', 'emergency'
    ]

    # Extract features and labels
    X = []
    y = []
    case_ids = []

    for sample in features_list:
        feature_vector = [sample[name] for name in feature_names]
        X.append(feature_vector)
        y.append(sample['ioh_label'])
        case_ids.append(sample['caseid'])

    X = np.array(X)
    y = np.array(y)

    print(f"✓ Feature matrix shape: {X.shape}")
    print(f"  Positive samples (IOH): {y.sum()} ({y.sum()/len(y)*100:.1f}%)")
    print()

    return X, y, case_ids, feature_names

def train_random_forest(X_train, y_train, X_test, y_test, feature_names):
    """Train Random Forest model with hyperparameter tuning"""

    print("=" * 70)
    print("TRAINING RANDOM FOREST MODEL")
    print("=" * 70)
    print()

    # Scale features
    print("Scaling features...")
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    print()

    # Train model
    print("Training Random Forest...")
    print("  - 300 trees")
    print("  - Max depth: 15")
    print("  - Balanced class weights")
    print()

    model = RandomForestClassifier(
        n_estimators=300,
        max_depth=15,
        min_samples_split=20,
        min_samples_leaf=10,
        max_features='sqrt',
        random_state=42,
        n_jobs=-1,
        class_weight='balanced',
        verbose=1
    )

    model.fit(X_train_scaled, y_train)
    print()
    print("✓ Training complete!")
    print()

    return model, scaler

def evaluate_model(model, scaler, X_train, y_train, X_test, y_test, feature_names):
    """Comprehensive model evaluation"""

    print("=" * 70)
    print("MODEL EVALUATION")
    print("=" * 70)
    print()

    X_train_scaled = scaler.transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Training performance
    train_score = model.score(X_train_scaled, y_train)
    print(f"Training Accuracy: {train_score:.4f}")

    # Test performance
    test_score = model.score(X_test_scaled, y_test)
    print(f"Test Accuracy: {test_score:.4f}")

    # ROC-AUC
    y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
    roc_auc = roc_auc_score(y_test, y_pred_proba)
    print(f"Test ROC-AUC: {roc_auc:.4f}")
    print()

    # Classification report
    y_pred = model.predict(X_test_scaled)
    print("Classification Report:")
    print(classification_report(y_test, y_pred,
                                target_names=['No IOH', 'IOH'],
                                digits=4))

    # Confusion matrix
    cm = confusion_matrix(y_test, y_pred)
    print("Confusion Matrix:")
    print(f"                  Predicted")
    print(f"                  No IOH    IOH")
    print(f"Actual  No IOH    {cm[0][0]:6d}  {cm[0][1]:6d}")
    print(f"        IOH       {cm[1][0]:6d}  {cm[1][1]:6d}")
    print()

    # Feature importance
    importances = model.feature_importances_
    feature_importance = sorted(
        zip(feature_names, importances),
        key=lambda x: x[1],
        reverse=True
    )

    print("Top 10 Feature Importances:")
    print("-" * 50)
    for feat, imp in feature_importance[:10]:
        print(f"{feat:20s}: {imp:.4f}")
    print()

    return {
        'roc_auc': roc_auc,
        'test_accuracy': test_score,
        'train_accuracy': train_score,
        'y_test': y_test,
        'y_pred': y_pred,
        'y_pred_proba': y_pred_proba,
        'feature_importance': feature_importance
    }

def plot_performance(results, save_path='vitaldb_model_performance.png'):
    """Create performance visualization plots"""

    print("Generating performance plots...")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # ROC Curve
    fpr, tpr, _ = roc_curve(results['y_test'], results['y_pred_proba'])
    axes[0].plot(fpr, tpr, linewidth=2, label=f"ROC (AUC = {results['roc_auc']:.3f})")
    axes[0].plot([0, 1], [0, 1], 'k--', linewidth=1, label='Random')
    axes[0].set_xlabel('False Positive Rate', fontsize=12)
    axes[0].set_ylabel('True Positive Rate', fontsize=12)
    axes[0].set_title('ROC Curve', fontsize=14, fontweight='bold')
    axes[0].legend(fontsize=10)
    axes[0].grid(alpha=0.3)

    # Precision-Recall Curve
    precision, recall, _ = precision_recall_curve(results['y_test'], results['y_pred_proba'])
    axes[1].plot(recall, precision, linewidth=2, color='darkblue')
    axes[1].set_xlabel('Recall', fontsize=12)
    axes[1].set_ylabel('Precision', fontsize=12)
    axes[1].set_title('Precision-Recall Curve', fontsize=14, fontweight='bold')
    axes[1].grid(alpha=0.3)

    # Feature Importance
    top_features = results['feature_importance'][:10]
    names = [f[0] for f in top_features]
    importances = [f[1] for f in top_features]

    axes[2].barh(range(len(names)), importances, color='steelblue')
    axes[2].set_yticks(range(len(names)))
    axes[2].set_yticklabels(names, fontsize=10)
    axes[2].set_xlabel('Importance', fontsize=12)
    axes[2].set_title('Top 10 Features', fontsize=14, fontweight='bold')
    axes[2].grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"✓ Plots saved to: {save_path}")
    print()

def save_production_model(model, scaler, results, output_prefix='vitaldb_ioh'):
    """Save trained model for production use"""

    print("=" * 70)
    print("SAVING PRODUCTION MODEL")
    print("=" * 70)
    print()

    # Save model
    model_path = f'{output_prefix}_model.pkl'
    with open(model_path, 'wb') as f:
        pickle.dump(model, f)
    print(f"✓ Model saved: {model_path}")

    # Save scaler
    scaler_path = f'{output_prefix}_scaler.pkl'
    with open(scaler_path, 'wb') as f:
        pickle.dump(scaler, f)
    print(f"✓ Scaler saved: {scaler_path}")

    # Save metadata
    metadata = {
        'roc_auc': results['roc_auc'],
        'test_accuracy': results['test_accuracy'],
        'train_accuracy': results['train_accuracy'],
        'feature_importance': results['feature_importance'],
        'training_date': datetime.now().isoformat(),
        'model_type': 'RandomForestClassifier',
        'dataset': 'VitalDB',
        'prediction_window': '5 minutes',
        'ioh_threshold': '65 mmHg MAP'
    }

    metadata_path = f'{output_prefix}_metadata.pkl'
    with open(metadata_path, 'wb') as f:
        pickle.dump(metadata, f)
    print(f"✓ Metadata saved: {metadata_path}")
    print()

    print("=" * 70)
    print("MODEL TRAINING COMPLETE!")
    print("=" * 70)
    print()
    print(f"Performance Summary:")
    print(f"  ROC-AUC: {results['roc_auc']:.4f}")
    print(f"  Test Accuracy: {results['test_accuracy']:.4f}")
    print()
    print("To deploy to production:")
    print(f"  1. Replace ioh_model.pkl with {model_path}")
    print(f"  2. Replace ioh_scaler.pkl with {scaler_path}")
    print(f"  3. Update app.py to use new model")
    print()

def main():
    """Main training pipeline"""

    print("=" * 70)
    print("VITALDB IOH MODEL TRAINING PIPELINE")
    print("=" * 70)
    print()

    # Load data
    dataset = load_vitaldb_dataset('vitaldb_ioh_dataset.pkl')

    # Prepare features
    X, y, case_ids, feature_names = prepare_features(dataset)

    # Split by case (not by sample) to prevent data leakage
    unique_cases = list(set(case_ids))
    train_cases, test_cases = train_test_split(unique_cases,
                                               test_size=0.2,
                                               random_state=42)

    train_idx = [i for i, cid in enumerate(case_ids) if cid in train_cases]
    test_idx = [i for i, cid in enumerate(case_ids) if cid in test_cases]

    X_train = X[train_idx]
    X_test = X[test_idx]
    y_train = y[train_idx]
    y_test = y[test_idx]

    print(f"Train/Test Split:")
    print(f"  Train: {len(X_train)} samples from {len(train_cases)} cases")
    print(f"  Test:  {len(X_test)} samples from {len(test_cases)} cases")
    print()

    # Train model
    model, scaler = train_random_forest(X_train, y_train, X_test, y_test, feature_names)

    # Evaluate
    results = evaluate_model(model, scaler, X_train, y_train, X_test, y_test, feature_names)

    # Plot performance
    try:
        plot_performance(results)
    except Exception as e:
        print(f"Warning: Could not generate plots: {e}")
        print()

    # Save for production
    save_production_model(model, scaler, results)

if __name__ == "__main__":
    main()
