"""
IOH Model Classes

This module contains the model classes used for IOH prediction.
It must be importable for pickle to deserialize the model files.
"""

class IOHModelStub:
    """
    RandomForest-compatible model for IOH prediction
    Provides realistic predictions based on clinical heuristics
    """

    def __init__(self):
        self.feature_importances_ = [
            0.08,  # age
            0.02,  # sex
            0.03,  # bmi
            0.07,  # asa
            0.10,  # baseline_map
            0.04,  # baseline_hr
            0.18,  # current_map (most important)
            0.15,  # map_5min
            0.12,  # map_10min
            0.03,  # surgery_duration
            0.05,  # vasopressor
            0.06,  # surgery_type
            0.05,  # induction_agent
            0.02,  # emergency
        ]

    def predict_proba(self, X):
        """
        Predict IOH probability based on clinical features

        Features (in order):
        0: age, 1: sex, 2: bmi, 3: asa, 4: baseline_map, 5: baseline_hr,
        6: current_map, 7: map_5min, 8: map_10min, 9: surgery_duration,
        10: vasopressor, 11: surgery_type, 12: induction_agent, 13: emergency
        """
        results = []

        for features in X:
            risk_score = 0

            age = features[0]
            asa = features[3]
            baseline_map = features[4]
            current_map = features[6]
            map_5min = features[7]
            map_10min = features[8]
            vasopressor = features[10]
            surgery_type = features[11]
            induction_agent = features[12]
            emergency = features[13]
            bmi = features[2]

            # MAP trend (most important predictor)
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
            if baseline_map > 0:
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

            # Vasopressor use (indicates existing instability)
            if vasopressor > 0:
                risk_score += 10

            # Induction agent (propofol causes more hypotension)
            if induction_agent == 0:  # propofol
                risk_score += 8

            # Emergency surgery
            if emergency == 1:
                risk_score += 10

            # BMI extremes
            if bmi < 20 or bmi > 35:
                risk_score += 5

            # Convert risk score to probability (0-1)
            # Max risk score ~200, but cap at 120 for 95% prob
            ioh_prob = min(0.95, max(0.05, risk_score / 120))

            # Return [prob_no_ioh, prob_ioh]
            results.append([1 - ioh_prob, ioh_prob])

        return results

    def predict(self, X):
        """Binary prediction"""
        proba = self.predict_proba(X)
        return [1 if p[1] > 0.5 else 0 for p in proba]


class IOHScalerStub:
    """
    StandardScaler-compatible scaler for IOH model
    Pass-through scaler (features are already in appropriate ranges)
    """

    def __init__(self):
        # Mock scaler parameters (14 features)
        self.mean_ = [0] * 14
        self.scale_ = [1] * 14

    def transform(self, X):
        """Pass through without scaling"""
        return X

    def fit_transform(self, X):
        """Pass through without scaling"""
        return X
