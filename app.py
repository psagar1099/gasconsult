from flask import Flask, request, render_template_string, session, redirect, url_for
from Bio import Entrez
import openai
import os
import re
import secrets

app = Flask(__name__)
app.secret_key = os.getenv('FLASK_SECRET_KEY', secrets.token_hex(32))

Entrez.email = "psagar1099@gmail.com"
Entrez.api_key = "f239ceccf2b12431846e6c03ffe29691ac08"

openai_client = openai.OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

def clean_query(query):
    """Strip conversational filler from user queries to extract the core medical topic."""
    q = query.lower().strip()

    # Remove common conversational phrases
    filler_phrases = [
        r"^(can you |could you |please |pls )",
        r"^(tell me about |tell me |explain |describe )",
        r"^(what is |what's |what are |whats )",
        r"^(what do you know about |what does the evidence say about )",
        r"^(what's the evidence for |what is the evidence for |evidence for )",
        r"^(i want to know about |i'd like to know about |i need to know about )",
        r"^(give me info on |give me information on |info on )",
        r"^(how does |how do |how is |how are )",
        r"^(why does |why do |why is |why are )",
        r"^(search for |look up |find |show me )",
        r"^(help me understand |help me with )",
        r"(indications to use |indications for |indication for |when to use )",
    ]

    # Apply each pattern repeatedly until no more matches
    for pattern in filler_phrases:
        q = re.sub(pattern, "", q, flags=re.IGNORECASE)

    # Remove trailing question marks and extra whitespace
    q = re.sub(r"\?+$", "", q).strip()
    q = re.sub(r"\s+", " ", q)

    return q if q else query  # Return original if cleaning removed everything

def detect_and_calculate(query, context_hint=None):
    """Detect calculation requests and perform medical calculations."""
    q = query.lower()

    # If context_hint provided, include it in detection
    if context_hint:
        q = context_hint + " " + q

    # Extract all numbers from the query
    numbers = re.findall(r'\d+\.?\d*', query)

    # Maximum Allowable Blood Loss (MABL)
    if any(term in q for term in ['mabl', 'maximum allowable blood loss', 'max blood loss']):
        if len(numbers) >= 3:
            try:
                ebv = float(numbers[0])  # Estimated blood volume (mL)
                hi = float(numbers[1])   # Initial Hct/Hgb
                hf = float(numbers[2])   # Final Hct/Hgb
                mabl = ebv * (hi - hf) / hi
                return f"""
                <h3>Maximum Allowable Blood Loss (MABL) Calculation</h3>
                <p><strong>Formula:</strong> MABL = EBV × (Hi - Hf) / Hi</p>
                <p><strong>Given:</strong></p>
                <ul>
                    <li>Estimated Blood Volume (EBV): {ebv:.0f} mL</li>
                    <li>Initial Hematocrit/Hemoglobin: {hi:.1f}</li>
                    <li>Final (acceptable) Hematocrit/Hemoglobin: {hf:.1f}</li>
                </ul>
                <p><strong>Result: MABL = {mabl:.0f} mL</strong></p>
                <p><em>Note: This is an estimate. Clinical judgment and patient condition should guide transfusion decisions.</em></p>
                """
            except:
                pass
        else:
            # Not enough numbers provided
            return f"""
            <h3>Maximum Allowable Blood Loss (MABL) Calculator</h3>
            <p>To calculate MABL, I need three values:</p>
            <ol>
                <li><strong>Estimated Blood Volume (EBV)</strong> in mL</li>
                <li><strong>Initial Hematocrit/Hemoglobin</strong> (e.g., 42 for Hct or 14 for Hgb)</li>
                <li><strong>Final/Acceptable Hematocrit/Hemoglobin</strong> (e.g., 30 for Hct or 10 for Hgb)</li>
            </ol>
            <p><strong>Example query:</strong> "Calculate MABL for 5000 mL blood volume, 42 initial Hct, 30 final Hct"</p>
            <p><em>You provided {len(numbers)} number(s). Please provide all three values.</em></p>
            """

    # Ideal Body Weight (IBW)
    if any(term in q for term in ['ibw', 'ideal body weight', 'ideal weight']):
        if len(numbers) >= 1:
            try:
                height_cm = float(numbers[0])
                is_male = any(word in q for word in ['male', 'man', 'm,'])
                is_female = any(word in q for word in ['female', 'woman', 'f,'])

                if is_male:
                    ibw = 50 + 0.91 * (height_cm - 152.4)
                    sex = "Male"
                elif is_female:
                    ibw = 45.5 + 0.91 * (height_cm - 152.4)
                    sex = "Female"
                else:
                    ibw_m = 50 + 0.91 * (height_cm - 152.4)
                    ibw_f = 45.5 + 0.91 * (height_cm - 152.4)
                    return f"""
                    <h3>Ideal Body Weight (IBW) Calculation</h3>
                    <p><strong>Height:</strong> {height_cm:.1f} cm</p>
                    <p><strong>Results:</strong></p>
                    <ul>
                        <li>Male IBW: {ibw_m:.1f} kg</li>
                        <li>Female IBW: {ibw_f:.1f} kg</li>
                    </ul>
                    <p><em>Tip: Specify sex (male/female) for a more specific result.</em></p>
                    """

                return f"""
                <h3>Ideal Body Weight (IBW) Calculation</h3>
                <p><strong>Formula ({sex}):</strong> {sex} formula using Devine equation</p>
                <p><strong>Height:</strong> {height_cm:.1f} cm</p>
                <p><strong>Result: IBW = {ibw:.1f} kg</strong></p>
                """
            except:
                pass
        else:
            return f"""
            <h3>Ideal Body Weight (IBW) Calculator</h3>
            <p>To calculate IBW, I need:</p>
            <ol>
                <li><strong>Height</strong> in cm</li>
                <li><strong>Sex</strong> (male or female) - optional, will show both if not specified</li>
            </ol>
            <p><strong>Example queries:</strong></p>
            <ul>
                <li>"Calculate IBW for 175 cm male"</li>
                <li>"Ideal body weight 165 cm female"</li>
            </ul>
            """

    # Body Surface Area (BSA)
    if any(term in q for term in ['bsa', 'body surface area', 'surface area']):
        if len(numbers) >= 2:
            try:
                weight = float(numbers[0])
                height = float(numbers[1])
                # Mosteller formula
                bsa = ((weight * height) / 3600) ** 0.5
                return f"""
                <h3>Body Surface Area (BSA) Calculation</h3>
                <p><strong>Formula:</strong> Mosteller formula: √((weight × height) / 3600)</p>
                <p><strong>Given:</strong></p>
                <ul>
                    <li>Weight: {weight:.1f} kg</li>
                    <li>Height: {height:.1f} cm</li>
                </ul>
                <p><strong>Result: BSA = {bsa:.2f} m²</strong></p>
                """
            except:
                pass
        else:
            return f"""
            <h3>Body Surface Area (BSA) Calculator</h3>
            <p>To calculate BSA, I need two values:</p>
            <ol>
                <li><strong>Weight</strong> in kg</li>
                <li><strong>Height</strong> in cm</li>
            </ol>
            <p><strong>Example query:</strong> "Calculate BSA for 70 kg and 175 cm"</p>
            <p><em>You provided {len(numbers)} number(s). Please provide both weight and height.</em></p>
            """

    # Maintenance Fluids (4-2-1 rule)
    if any(term in q for term in ['maintenance fluid', 'fluid requirement', '4-2-1', 'hourly fluid']):
        if len(numbers) >= 1:
            try:
                weight = float(numbers[0])
                if weight <= 10:
                    rate = weight * 4
                elif weight <= 20:
                    rate = 40 + (weight - 10) * 2
                else:
                    rate = 60 + (weight - 20) * 1

                return f"""
                <h3>Maintenance Fluid Requirement (4-2-1 Rule)</h3>
                <p><strong>Weight:</strong> {weight:.1f} kg</p>
                <p><strong>Calculation:</strong></p>
                <ul>
                    <li>First 10 kg: 4 mL/kg/hr</li>
                    <li>Second 10 kg: 2 mL/kg/hr</li>
                    <li>Each kg above 20: 1 mL/kg/hr</li>
                </ul>
                <p><strong>Result: {rate:.0f} mL/hr</strong></p>
                <p><strong>Daily requirement: {rate * 24:.0f} mL/day</strong></p>
                """
            except:
                pass
        else:
            return f"""
            <h3>Maintenance Fluid Requirement Calculator</h3>
            <p>To calculate maintenance fluids using the 4-2-1 rule, I need:</p>
            <ol>
                <li><strong>Patient weight</strong> in kg</li>
            </ol>
            <p><strong>Example queries:</strong></p>
            <ul>
                <li>"Maintenance fluids for 25 kg"</li>
                <li>"Calculate hourly fluid requirement for 70 kg patient"</li>
            </ul>
            """

    # QTc (Corrected QT interval)
    if any(term in q for term in ['qtc', 'corrected qt', 'qt interval']):
        if len(numbers) >= 2:
            try:
                qt = float(numbers[0])
                rr = float(numbers[1])
                # Bazett's formula
                qtc = qt / (rr ** 0.5)
                interpretation = "Normal" if qtc < 450 else "Prolonged (>450ms - risk of arrhythmia)"

                return f"""
                <h3>QTc (Corrected QT Interval) Calculation</h3>
                <p><strong>Formula:</strong> Bazett's formula: QTc = QT / √RR</p>
                <p><strong>Given:</strong></p>
                <ul>
                    <li>QT interval: {qt:.0f} ms</li>
                    <li>RR interval: {rr:.0f} ms</li>
                </ul>
                <p><strong>Result: QTc = {qtc:.0f} ms</strong></p>
                <p><strong>Interpretation:</strong> {interpretation}</p>
                """
            except:
                pass
        else:
            return f"""
            <h3>QTc (Corrected QT Interval) Calculator</h3>
            <p>To calculate QTc using Bazett's formula, I need:</p>
            <ol>
                <li><strong>QT interval</strong> in milliseconds</li>
                <li><strong>RR interval</strong> in milliseconds</li>
            </ol>
            <p><strong>Example query:</strong> "Calculate QTc for QT 400 and RR 800"</p>
            <p><em>You provided {len(numbers)} number(s). Please provide both QT and RR intervals.</em></p>
            """

    return None  # No calculation detected

PREOP_HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Pre-Op Assessment — gasconsult.ai</title>
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=3">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=3">
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: -apple-system, BlinkMacSystemFont, 'SF Pro Display', 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background: #ffffff;
            color: #0A3D62;
            line-height: 1.6;
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
        }

        /* Navigation */
        nav {
            background: #ffffff;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.08);
            padding: 18px 0;
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            z-index: 1000;
            border-bottom: 1px solid #f0f0f0;
        }

        nav .container {
            max-width: 1200px;
            margin: 0 auto;
            padding: 0 40px;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }

        .logo-text {
            font-size: 1.5rem;
            font-weight: 700;
            color: #0A3D62;
            letter-spacing: -0.5px;
            text-decoration: none;
        }

        .logo-symbol {
            display: inline-block;
            margin-right: 10px;
            font-size: 1.6rem;
            color: #FF6B35;
        }

        .nav-actions {
            display: flex;
            gap: 12px;
            align-items: center;
        }

        .nav-link {
            color: #0A3D62;
            text-decoration: none;
            font-size: 0.95rem;
            font-weight: 600;
            padding: 8px 16px;
            border-radius: 8px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            background: #E8F4FD;
        }

        /* Main Content */
        .preop-container {
            max-width: 900px;
            margin: 0 auto;
            padding: 100px 20px 60px;
        }

        .preop-header {
            text-align: center;
            margin-bottom: 40px;
        }

        .preop-header h1 {
            font-size: 2.5rem;
            color: #0A3D62;
            margin-bottom: 12px;
        }

        .preop-header p {
            font-size: 1.1rem;
            color: #666;
        }

        /* Form Sections */
        .form-section {
            background: #E8F4FD;
            border-radius: 12px;
            padding: 28px;
            margin-bottom: 24px;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.04);
        }

        .form-section h2 {
            color: #0A3D62;
            font-size: 1.4rem;
            margin-bottom: 20px;
            font-weight: 700;
        }

        .form-row {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 16px;
            margin-bottom: 16px;
        }

        .form-group {
            display: flex;
            flex-direction: column;
        }

        label {
            font-size: 0.9rem;
            font-weight: 600;
            color: #0A3D62;
            margin-bottom: 6px;
        }

        input[type="text"],
        input[type="number"],
        textarea,
        select {
            padding: 12px;
            border: 2px solid #d0e8f7;
            border-radius: 8px;
            font-size: 1rem;
            font-family: inherit;
            background: white;
            color: #0A3D62;
            transition: all 0.2s ease;
        }

        input:focus,
        textarea:focus,
        select:focus {
            outline: none;
            border-color: #FF6B35;
            box-shadow: 0 0 0 3px rgba(255, 107, 53, 0.1);
        }

        textarea {
            resize: vertical;
            min-height: 80px;
        }

        .checkbox-group {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(180px, 1fr));
            gap: 12px;
            margin-top: 8px;
        }

        .checkbox-item {
            display: flex;
            align-items: center;
            gap: 8px;
        }

        .checkbox-item input[type="checkbox"] {
            width: 18px;
            height: 18px;
            cursor: pointer;
        }

        .checkbox-item label {
            margin: 0;
            cursor: pointer;
            font-weight: 500;
        }

        .submit-btn {
            background: #FF6B35;
            color: white;
            padding: 14px 32px;
            border-radius: 24px;
            font-size: 1.05rem;
            font-weight: 600;
            border: none;
            cursor: pointer;
            transition: all 0.2s ease;
            width: 100%;
            margin-top: 20px;
        }

        .submit-btn:hover {
            background: #ff5722;
            transform: translateY(-2px);
            box-shadow: 0 6px 20px rgba(255, 107, 53, 0.3);
        }

        /* Summary Display */
        .summary-container {
            background: #E8F4FD;
            border-radius: 12px;
            padding: 32px;
            margin-top: 40px;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.04);
        }

        .summary-container h2 {
            color: #0A3D62;
            font-size: 1.8rem;
            margin-bottom: 24px;
        }

        .summary-content {
            color: #0A3D62;
            line-height: 1.8;
        }

        .summary-content h3 {
            color: #FF6B35;
            font-size: 1.3rem;
            margin-top: 20px;
            margin-bottom: 12px;
        }

        .summary-content strong {
            color: #0A3D62;
        }

        .ref-item {
            padding: 8px 0;
            transition: padding-left 0.2s ease;
        }

        .ref-item:hover {
            padding-left: 6px;
        }

        .ref-item a {
            color: #0A3D62;
            text-decoration: none;
            font-weight: 500;
            transition: color 0.2s ease;
        }

        .ref-item a:hover {
            color: #FF6B35;
            text-decoration: underline;
        }

        .auto-calc {
            background: white;
            padding: 12px;
            border-radius: 8px;
            margin-top: 12px;
            border: 2px solid #d0e8f7;
        }

        .auto-calc strong {
            color: #FF6B35;
        }
    </style>
    <script>
        // Auto-calculate BMI and IBW
        function calculateMetrics() {
            const weight = parseFloat(document.getElementById('weight').value);
            const height = parseFloat(document.getElementById('height').value);
            const sex = document.querySelector('input[name="sex"]:checked')?.value;

            let results = '';

            if (weight && height) {
                // BMI
                const bmi = (weight / ((height / 100) ** 2)).toFixed(1);
                results += `<strong>BMI:</strong> ${bmi} kg/m²<br>`;

                // IBW
                if (sex === 'male') {
                    const ibw = (50 + 0.91 * (height - 152.4)).toFixed(1);
                    results += `<strong>IBW (Male):</strong> ${ibw} kg`;
                } else if (sex === 'female') {
                    const ibw = (45.5 + 0.91 * (height - 152.4)).toFixed(1);
                    results += `<strong>IBW (Female):</strong> ${ibw} kg`;
                }
            }

            document.getElementById('autoCalc').innerHTML = results || 'Enter weight, height, and sex';
        }
    </script>
</head>
<body>
    <nav>
        <div class="container">
            <a href="/" class="logo-text">
                <span class="logo-symbol">⚕</span>gasconsult.ai
            </a>
            <div class="nav-actions">
                <a href="/" class="nav-link">Chat</a>
                <a href="/preop" class="nav-link">Pre-Op</a>
            </div>
        </div>
    </nav>

    <div class="preop-container">
        <div class="preop-header">
            <h1>Pre-Operative Assessment</h1>
            <p>Evidence-based risk stratification and recommendations</p>
        </div>

        {% if not summary %}
        <form method="post" action="/preop">
            <!-- Demographics -->
            <div class="form-section">
                <h2>1. Patient Demographics</h2>
                <div class="form-row">
                    <div class="form-group">
                        <label for="age">Age (years)</label>
                        <input type="number" id="age" name="age" required>
                    </div>
                    <div class="form-group">
                        <label for="weight">Weight (kg)</label>
                        <input type="number" id="weight" name="weight" step="0.1" required onchange="calculateMetrics()">
                    </div>
                    <div class="form-group">
                        <label for="height">Height (cm)</label>
                        <input type="number" id="height" name="height" step="0.1" required onchange="calculateMetrics()">
                    </div>
                </div>
                <div class="form-row">
                    <div class="form-group">
                        <label>Sex</label>
                        <div style="display: flex; gap: 20px; margin-top: 8px;">
                            <div class="checkbox-item">
                                <input type="radio" id="male" name="sex" value="male" onchange="calculateMetrics()" required>
                                <label for="male">Male</label>
                            </div>
                            <div class="checkbox-item">
                                <input type="radio" id="female" name="sex" value="female" onchange="calculateMetrics()" required>
                                <label for="female">Female</label>
                            </div>
                        </div>
                    </div>
                </div>
                <div class="auto-calc" id="autoCalc">
                    Enter weight, height, and sex to calculate BMI and IBW
                </div>
            </div>

            <!-- Comorbidities -->
            <div class="form-section">
                <h2>2. Comorbidities</h2>
                <div class="checkbox-group">
                    <div class="checkbox-item">
                        <input type="checkbox" id="dm" name="comorbidities" value="Diabetes Mellitus">
                        <label for="dm">Diabetes Mellitus</label>
                    </div>
                    <div class="checkbox-item">
                        <input type="checkbox" id="htn" name="comorbidities" value="Hypertension">
                        <label for="htn">Hypertension</label>
                    </div>
                    <div class="checkbox-item">
                        <input type="checkbox" id="cad" name="comorbidities" value="Coronary Artery Disease">
                        <label for="cad">CAD</label>
                    </div>
                    <div class="checkbox-item">
                        <input type="checkbox" id="chf" name="comorbidities" value="Heart Failure">
                        <label for="chf">Heart Failure</label>
                    </div>
                    <div class="checkbox-item">
                        <input type="checkbox" id="copd" name="comorbidities" value="COPD">
                        <label for="copd">COPD</label>
                    </div>
                    <div class="checkbox-item">
                        <input type="checkbox" id="asthma" name="comorbidities" value="Asthma">
                        <label for="asthma">Asthma</label>
                    </div>
                    <div class="checkbox-item">
                        <input type="checkbox" id="osa" name="comorbidities" value="Obstructive Sleep Apnea">
                        <label for="osa">OSA</label>
                    </div>
                    <div class="checkbox-item">
                        <input type="checkbox" id="ckd" name="comorbidities" value="Chronic Kidney Disease">
                        <label for="ckd">CKD</label>
                    </div>
                    <div class="checkbox-item">
                        <input type="checkbox" id="stroke" name="comorbidities" value="Prior Stroke">
                        <label for="stroke">Prior Stroke</label>
                    </div>
                    <div class="checkbox-item">
                        <input type="checkbox" id="afib" name="comorbidities" value="Atrial Fibrillation">
                        <label for="afib">Atrial Fibrillation</label>
                    </div>
                </div>
            </div>

            <!-- Medications -->
            <div class="form-section">
                <h2>3. Current Medications</h2>
                <div class="form-group">
                    <label for="medications">List all medications (include anticoagulants, antiplatelets, insulin, etc.)</label>
                    <textarea id="medications" name="medications" placeholder="e.g., Aspirin 81mg daily, Metoprolol 50mg BID, Apixaban 5mg BID..."></textarea>
                </div>
            </div>

            <!-- Labs -->
            <div class="form-section">
                <h2>4. Laboratory Values</h2>
                <div class="form-row">
                    <div class="form-group">
                        <label for="hgb">Hemoglobin (g/dL)</label>
                        <input type="number" id="hgb" name="hgb" step="0.1">
                    </div>
                    <div class="form-group">
                        <label for="plt">Platelets (×10³/μL)</label>
                        <input type="number" id="plt" name="plt">
                    </div>
                    <div class="form-group">
                        <label for="cr">Creatinine (mg/dL)</label>
                        <input type="number" id="cr" name="cr" step="0.01">
                    </div>
                    <div class="form-group">
                        <label for="inr">INR</label>
                        <input type="number" id="inr" name="inr" step="0.1">
                    </div>
                </div>
            </div>

            <!-- Procedure -->
            <div class="form-section">
                <h2>5. Surgical Procedure</h2>
                <div class="form-group">
                    <label for="procedure">Procedure Type</label>
                    <input type="text" id="procedure" name="procedure" placeholder="e.g., Total Knee Arthroplasty, CABG, Laparoscopic Cholecystectomy..." required>
                </div>
                <div class="form-row">
                    <div class="form-group">
                        <label for="surgery_risk">Surgery Risk Category</label>
                        <select id="surgery_risk" name="surgery_risk" required>
                            <option value="">Select...</option>
                            <option value="Low">Low Risk (&lt;1% cardiac risk)</option>
                            <option value="Intermediate">Intermediate Risk (1-5% cardiac risk)</option>
                            <option value="High">High Risk (&gt;5% cardiac risk)</option>
                        </select>
                    </div>
                </div>
            </div>

            <!-- Additional Info -->
            <div class="form-section">
                <h2>6. Additional Information</h2>
                <div class="form-row">
                    <div class="form-group">
                        <label for="npo">NPO Status</label>
                        <input type="text" id="npo" name="npo" placeholder="e.g., NPO since midnight, last PO intake 6 hours ago...">
                    </div>
                    <div class="form-group">
                        <label for="allergies">Allergies</label>
                        <input type="text" id="allergies" name="allergies" placeholder="e.g., PCN (rash), Morphine (nausea), NKDA...">
                    </div>
                </div>
            </div>

            <button type="submit" class="submit-btn">Generate Evidence-Based Assessment</button>
        </form>
        {% else %}
        <div class="summary-container">
            <h2>Pre-Operative Assessment Summary</h2>
            <div class="summary-content">
                {{ summary|safe }}
            </div>
            {% if references %}
            <div style="margin-top: 30px; padding-top: 20px; border-top: 2px solid rgba(10, 61, 98, 0.15);">
                <h3 style="color: #FF6B35; margin-bottom: 16px;">References:</h3>
                {% for ref in references %}
                <div class="ref-item">
                    <a href="https://pubmed.ncbi.nlm.nih.gov/{{ ref.pmid }}/" target="_blank">
                        [{{ loop.index }}] {{ ref.title }} ({{ ref.year }})
                    </a>
                </div>
                {% endfor %}
            </div>
            {% endif %}
        </div>
        <div style="text-align: center; margin-top: 30px;">
            <a href="/preop" class="submit-btn" style="display: inline-block; width: auto; text-decoration: none;">New Assessment</a>
        </div>
        {% endif %}
    </div>
</body>
</html>
"""

HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>gasconsult.ai — Evidence-Based Anesthesiology</title>
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=3">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=3">
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: -apple-system, BlinkMacSystemFont, 'SF Pro Display', 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background: #ffffff;
            color: #0A3D62;
            line-height: 1.6;
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
        }

        /* Navigation */
        nav {
            background: #ffffff;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.08);
            padding: 18px 0;
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            z-index: 1000;
            border-bottom: 1px solid #f0f0f0;
        }

        nav .container {
            max-width: 1200px;
            margin: 0 auto;
            padding: 0 40px;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }

        .logo-text {
            font-size: 1.5rem;
            font-weight: 700;
            color: #0A3D62;
            letter-spacing: -0.5px;
            transition: opacity 0.2s ease;
        }

        .logo-text:hover {
            opacity: 0.8;
        }

        .logo-symbol {
            display: inline-block;
            margin-right: 10px;
            font-size: 1.6rem;
            color: #FF6B35;
        }

        .nav-actions {
            display: flex;
            gap: 12px;
            align-items: center;
        }

        .nav-link {
            color: #0A3D62;
            text-decoration: none;
            font-size: 0.95rem;
            font-weight: 600;
            padding: 8px 16px;
            border-radius: 8px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            background: #E8F4FD;
        }

        .new-chat-btn {
            background: #FF6B35;
            color: white;
            padding: 10px 24px;
            border-radius: 20px;
            font-size: 0.95rem;
            font-weight: 600;
            text-decoration: none;
            transition: all 0.2s ease;
            border: none;
            cursor: pointer;
        }

        .new-chat-btn:hover {
            background: #ff5722;
            transform: translateY(-1px);
            box-shadow: 0 4px 12px rgba(255, 107, 53, 0.3);
        }

        /* Welcome Screen */
        .welcome-screen {
            padding: 100px 40px 60px;
            text-align: center;
            max-width: 1000px;
            margin: 0 auto;
            animation: fadeIn 0.5s ease;
        }

        @keyframes fadeIn {
            from { opacity: 0; transform: translateY(20px); }
            to { opacity: 1; transform: translateY(0); }
        }

        .welcome-screen .hero-logo {
            font-size: 5rem;
            margin-bottom: 24px;
            color: #FF6B35;
        }

        .welcome-screen .tagline {
            font-size: 2rem;
            color: #0A3D62;
            margin-bottom: 50px;
            font-weight: 700;
            letter-spacing: -0.5px;
        }

        .welcome-screen .feature-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 24px;
            margin-bottom: 40px;
        }

        .feature-card {
            background: #E8F4FD;
            border-radius: 12px;
            padding: 32px 28px;
            transition: all 0.3s ease;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.04);
        }

        .feature-card:hover {
            transform: translateY(-4px);
            box-shadow: 0 8px 20px rgba(0, 0, 0, 0.1);
        }

        .feature-card h3 {
            color: #0A3D62;
            font-size: 1.15rem;
            margin-bottom: 12px;
            font-weight: 700;
        }

        .feature-card p {
            color: #555;
            font-size: 0.95rem;
            line-height: 1.6;
        }

        /* Main Content Area */
        .main-content {
            padding-top: 72px;
            min-height: 100vh;
            display: flex;
            flex-direction: column;
        }

        /* Chat Container */
        .chat-container {
            max-width: 900px;
            margin: 0 auto;
            width: 100%;
            flex: 1;
            display: flex;
            flex-direction: column;
            padding: 0 20px;
        }

        .chat-messages {
            flex: 1;
            overflow-y: auto;
            padding: 30px 20px;
            scroll-behavior: smooth;
        }

        .chat-messages::-webkit-scrollbar {
            width: 8px;
        }

        .chat-messages::-webkit-scrollbar-track {
            background: #f5f5f5;
        }

        .chat-messages::-webkit-scrollbar-thumb {
            background: #ddd;
            border-radius: 10px;
        }

        .chat-messages::-webkit-scrollbar-thumb:hover {
            background: #bbb;
        }

        /* Chat Messages */
        .message {
            margin-bottom: 28px;
            display: flex;
            animation: slideIn 0.3s ease;
        }

        @keyframes slideIn {
            from { opacity: 0; transform: translateY(10px); }
            to { opacity: 1; transform: translateY(0); }
        }

        .message.user {
            justify-content: flex-end;
        }

        .message.assistant {
            justify-content: flex-start;
        }

        .message-content {
            max-width: 85%;
            padding: 18px 24px;
            border-radius: 24px;
            font-size: 1rem;
            line-height: 1.7;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.08);
            transition: all 0.2s ease;
        }

        .message-content:hover {
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.12);
        }

        .message.user .message-content {
            background: #0A3D62;
            color: white;
        }

        .message.assistant .message-content {
            background: #E8F4FD;
            color: #0A3D62;
            border: 1px solid #d0e8f7;
        }

        .message-text h3 {
            font-size: 1.2rem;
            margin-top: 16px;
            margin-bottom: 12px;
            color: #0A3D62;
            font-weight: 700;
        }

        .message-text ul, .message-text ol {
            margin-left: 20px;
            margin-bottom: 12px;
        }

        .message-text p {
            margin-bottom: 12px;
        }

        .message-refs {
            margin-top: 20px;
            padding-top: 16px;
            border-top: 2px solid rgba(10, 61, 98, 0.15);
            font-size: 0.92rem;
        }

        .message-refs strong {
            display: block;
            margin-bottom: 12px;
            color: #FF6B35;
            font-weight: 700;
        }

        .ref-item {
            padding: 8px 0;
            transition: padding-left 0.2s ease;
        }

        .ref-item:hover {
            padding-left: 6px;
        }

        .ref-item a {
            color: #0A3D62;
            text-decoration: none;
            font-weight: 500;
            transition: color 0.2s ease;
        }

        .ref-item a:hover {
            color: #FF6B35;
            text-decoration: underline;
        }

        .message-meta {
            margin-top: 12px;
            font-size: 0.85rem;
            color: #0A3D62;
            opacity: 0.7;
            font-weight: 600;
        }

        /* Chat Input */
        .chat-input-container {
            background: #ffffff;
            padding: 20px;
            border-top: 1px solid #f0f0f0;
            position: sticky;
            bottom: 0;
        }

        .chat-form {
            position: relative;
            max-width: 900px;
            margin: 0 auto;
        }

        .chat-form textarea {
            width: 100%;
            padding: 15px 100px 15px 20px;
            font-size: 1rem;
            font-family: inherit;
            border: 2px solid #e0e0e0;
            border-radius: 24px;
            resize: none;
            transition: all 0.2s ease;
            background: #ffffff;
            color: #0A3D62;
            height: 52px;
            line-height: 1.3;
        }

        .chat-form textarea:focus {
            outline: none;
            border-color: #FF6B35;
            box-shadow: 0 0 0 3px rgba(255, 107, 53, 0.1);
        }

        .chat-form textarea::placeholder {
            color: #999;
        }

        .send-btn {
            position: absolute;
            right: 6px;
            top: 50%;
            transform: translateY(-50%);
            background: #FF6B35;
            color: white;
            border: none;
            width: 40px;
            height: 40px;
            border-radius: 50%;
            font-size: 1.2rem;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s ease;
            display: flex;
            align-items: center;
            justify-content: center;
            padding: 0;
        }

        .send-btn:hover {
            background: #ff5722;
            transform: translateY(-50%) scale(1.02);
            box-shadow: 0 4px 12px rgba(255, 107, 53, 0.3);
        }

        .send-btn:active {
            transform: translateY(-50%) scale(0.98);
        }

        /* Loading Animation */
        .loading-container {
            display: none;
            padding: 20px;
            text-align: center;
        }

        .loading-container.active {
            display: block;
            animation: fadeIn 0.3s ease;
        }

        .loading-dots {
            display: inline-flex;
            gap: 8px;
            align-items: center;
        }

        .loading-dots span {
            width: 12px;
            height: 12px;
            background: #FF6B35;
            border-radius: 50%;
            animation: bounce 1.4s infinite ease-in-out both;
        }

        .loading-dots span:nth-child(1) {
            animation-delay: -0.32s;
        }

        .loading-dots span:nth-child(2) {
            animation-delay: -0.16s;
        }

        @keyframes bounce {
            0%, 80%, 100% {
                transform: scale(0.8);
                opacity: 0.5;
            }
            40% {
                transform: scale(1.2);
                opacity: 1;
            }
        }

        /* Responsive Design */
        @media (max-width: 768px) {
            .welcome-screen .hero-logo {
                font-size: 4rem;
            }

            .welcome-screen .tagline {
                font-size: 1.5rem;
            }

            .welcome-screen .feature-grid {
                grid-template-columns: 1fr;
            }

            .message-content {
                max-width: 92%;
            }

            .chat-messages {
                padding: 20px 10px;
            }

            .chat-input-container {
                padding: 15px;
            }

            nav .container {
                padding: 0 20px;
            }

            .send-btn {
                position: relative;
                transform: none;
                width: 48px;
                height: 48px;
                margin-top: 10px;
                margin-left: auto;
                display: flex;
            }

            .send-btn:hover {
                transform: scale(1.05);
            }

            .send-btn:active {
                transform: scale(0.95);
            }

            .chat-form textarea {
                padding: 14px 20px;
                height: 50px;
            }

            .logo-text {
                font-size: 1.2rem;
            }
        }
    </style>
    <script>
        // Auto-scroll to bottom of chat on page load
        window.onload = function() {
            const chatMessages = document.getElementById('chatMessages');
            if (chatMessages) {
                chatMessages.scrollTop = chatMessages.scrollHeight;
            }
        };

        // Show loading animation when form submits
        document.addEventListener('DOMContentLoaded', function() {
            const chatForm = document.querySelector('.chat-form');
            const loadingContainer = document.getElementById('loadingIndicator');

            if (chatForm) {
                chatForm.addEventListener('submit', function(e) {
                    const textarea = chatForm.querySelector('textarea');
                    if (textarea.value.trim()) {
                        if (loadingContainer) {
                            loadingContainer.classList.add('active');
                        }
                    }
                });
            }
        });
    </script>
</head>
<body>
    <nav>
        <div class="container">
            <div class="logo-text">
                <span class="logo-symbol">⚕</span>gasconsult.ai
            </div>
            <div class="nav-actions">
                <a href="/" class="nav-link">Chat</a>
                <a href="/preop" class="nav-link">Pre-Op</a>
                {% if messages %}
                <a href="/clear" class="new-chat-btn">+ New Chat</a>
                {% endif %}
            </div>
        </div>
    </nav>

    <div class="main-content">
        {% if not messages %}
        <!-- Welcome Screen -->
        <div class="welcome-screen">
            <div class="hero-logo">⚕</div>
            <p class="tagline">Strictly Evidence-Based Anesthesiology</p>

            <div class="feature-grid">
                <div class="feature-card">
                    <h3>📚 PubMed-Backed</h3>
                    <p>Every answer sourced from peer-reviewed research and clinical guidelines</p>
                </div>
                <div class="feature-card">
                    <h3>🧮 Medical Calculators</h3>
                    <p>MABL, IBW, BSA, QTc, maintenance fluids, and more</p>
                </div>
                <div class="feature-card">
                    <h3>💬 Conversational AI</h3>
                    <p>Ask follow-ups, refine your questions, and explore topics naturally</p>
                </div>
            </div>
        </div>
        {% else %}
        <!-- Chat Container -->
        <div class="chat-container">
            <div class="chat-messages" id="chatMessages">
                {% for msg in messages %}
                    <div class="message {{ msg.role }}">
                        <div class="message-content">
                            {% if msg.role == 'user' %}
                                <div class="message-text">{{ msg.content }}</div>
                            {% else %}
                                <div class="message-text">{{ msg.content|safe }}</div>
                                {% if msg.references %}
                                <div class="message-refs">
                                    <strong>References:</strong>
                                    {% for ref in msg.references %}
                                    <div class="ref-item">
                                        <a href="https://pubmed.ncbi.nlm.nih.gov/{{ ref.pmid }}/" target="_blank">
                                            [{{ loop.index }}] {{ ref.title }} ({{ ref.year }})
                                        </a>
                                    </div>
                                    {% endfor %}
                                </div>
                                {% endif %}
                                {% if msg.num_papers > 0 %}
                                <div class="message-meta">📊 {{ msg.num_papers }} papers from PubMed</div>
                                {% endif %}
                            {% endif %}
                        </div>
                    </div>
                {% endfor %}
            </div>

            <!-- Loading Indicator -->
            <div id="loadingIndicator" class="loading-container">
                <div class="loading-dots">
                    <span></span>
                    <span></span>
                    <span></span>
                </div>
            </div>
        </div>
        {% endif %}

        <!-- Chat Input - Always Visible -->
        <div class="chat-input-container">
            <form method="post" action="/" class="chat-form">
                <textarea name="query" id="chatInput" placeholder="Ask anything about anesthesiology..." required rows="2"></textarea>
                <button type="submit" class="send-btn">↑</button>
            </form>
        </div>
    </div>
</body>
</html>
"""

@app.route("/", methods=["GET", "POST"])
def index():
    # Initialize conversation history in session
    if 'messages' not in session:
        session['messages'] = []

    if request.method == "POST":
        raw_query = request.form["query"].strip()

        # Add user message to conversation
        session['messages'].append({"role": "user", "content": raw_query})

        # Check if this is a calculation request
        # Look at previous messages for context
        context_hint = None
        if len(session['messages']) >= 3:
            # Check last few messages for calculator keywords
            last_msgs = session['messages'][-4:]
            for msg in last_msgs:
                content = msg.get('content', '').lower()
                if any(term in content for term in ['mabl', 'ibw', 'bsa', 'qtc', 'maintenance fluid', 'ideal body weight', 'body surface']):
                    # Extract the calculator type
                    for term in ['mabl', 'ibw', 'bsa', 'qtc', 'maintenance fluid']:
                        if term in content:
                            context_hint = term
                            break
                    break

        calc_result = detect_and_calculate(raw_query, context_hint=context_hint)

        if calc_result:
            # Add calculation result to conversation
            session['messages'].append({
                "role": "assistant",
                "content": calc_result,
                "references": [],
                "num_papers": 0
            })
            session.modified = True
            return redirect(url_for('index'))

        query = clean_query(raw_query)  # Strip conversational filler

        # Expand synonyms and medical abbreviations
        q = query.lower()
        q = q.replace(" versus ", " OR ")
        q = q.replace(" vs ", " OR ")
        q = q.replace(" vs. ", " OR ")
        q = q.replace(" over ", " OR ")
        q = q.replace(" compared to ", " OR ")
        q = q.replace(" compared with ", " OR ")
        q = q.replace("txa", '"tranexamic acid" OR TXA')
        q = q.replace("blood loss", '"blood loss" OR hemorrhage OR transfusion')
        q = q.replace("spine surgery", '"spine surgery" OR "spinal fusion" OR scoliosis')
        q = q.replace("peds", 'pediatric OR children OR peds')
        q = q.replace("pediatric", 'pediatric OR children OR peds')
        q = q.replace("ponv", 'PONV OR "postoperative nausea"')
        q = q.replace("propofol", '"propofol"[MeSH Terms] OR propofol')
        q = q.replace("etomidate", '"etomidate"[MeSH Terms] OR etomidate')

        search_term = (
            f'({q}) AND '
            f'(systematic review[pt] OR meta-analysis[pt] OR "randomized controlled trial"[pt] OR '
            f'"Cochrane Database Syst Rev"[ta] OR guideline[pt]) AND '
            f'("2015/01/01"[PDAT] : "3000"[PDAT])'
        )

        # Try anesthesiology-specific high-quality evidence first
        handle = Entrez.esearch(db="pubmed", term=f'anesthesiology[MeSH Terms] AND {search_term}', retmax=15, sort="relevance", api_key=Entrez.api_key)
        result = Entrez.read(handle)
        ids = result["IdList"]

        # Fallback 1: Try without anesthesiology restriction
        if not ids:
            handle = Entrez.esearch(db="pubmed", term=search_term, retmax=15, sort="relevance", api_key=Entrez.api_key)
            result = Entrez.read(handle)
            ids = result["IdList"]

        # Fallback 2: Drop publication type restrictions, keep recent papers
        if not ids:
            broader_search = f'({q}) AND ("2015/01/01"[PDAT] : "3000"[PDAT])'
            handle = Entrez.esearch(db="pubmed", term=broader_search, retmax=15, sort="relevance", api_key=Entrez.api_key)
            result = Entrez.read(handle)
            ids = result["IdList"]

        # Fallback 3: Even broader - any recent paper in anesthesiology
        if not ids:
            broadest_search = f'({q}) AND anesthesiology AND ("2010/01/01"[PDAT] : "3000"[PDAT])'
            handle = Entrez.esearch(db="pubmed", term=broadest_search, retmax=15, sort="relevance", api_key=Entrez.api_key)
            result = Entrez.read(handle)
            ids = result["IdList"]

        if not ids:
            error_msg = "<p>No relevant evidence found. Try rephrasing your question or using different medical terms.</p>"
            session['messages'].append({
                "role": "assistant",
                "content": error_msg,
                "references": [],
                "num_papers": 0
            })
            session.modified = True
            return redirect(url_for('index'))

        handle = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml", api_key=Entrez.api_key)
        papers = Entrez.read(handle)["PubmedArticle"]

        refs = []
        context = ""
        for p in papers[:12]:
            try:
                art = p["MedlineCitation"]["Article"]
                title = art.get("ArticleTitle", "No title")
                abstract = " ".join(str(t) for t in art.get("Abstract", {}).get("AbstractText", [])) if art.get("Abstract") else ""
                authors = ", ".join([a.get("LastName","") + " " + (a.get("ForeName","")[:1]+"." if a.get("ForeName") else "") for a in art.get("AuthorList",[])[:5]])
                journal = art["Journal"].get("Title", "Unknown")
                year = art["Journal"]["JournalIssue"]["PubDate"].get("Year", "N/A")
                pmid = p["MedlineCitation"]["PMID"]

                refs.append({"title": title, "authors": authors, "journal": journal, "year": year, "pmid": pmid})
                context += f"Title: {title}\nAbstract: {abstract}\nAuthors: {authors}\nJournal: {journal} ({year})\nPMID: {pmid}\n\n"
            except:
                continue

        num_papers = len(refs)

        # Build conversation context for GPT (last 10 messages)
        conversation_context = ""
        recent_messages = session['messages'][-10:]
        for msg in recent_messages:
            if msg['role'] == 'user':
                conversation_context += f"User: {msg['content']}\n"
            else:
                # Strip HTML for cleaner context
                content_text = re.sub('<[^<]+?>', '', msg.get('content', ''))
                conversation_context += f"Assistant: {content_text[:200]}...\n"

        # Create numbered reference list for citation
        ref_list = ""
        for i, ref in enumerate(refs, 1):
            ref_list += f"[{i}] {ref['title']} - {ref['authors']} ({ref['year']}) PMID: {ref['pmid']}\n"

        prompt = f"""You are a clinical expert anesthesiologist AI assistant designed for real-time decision support. Your responses must be immediately actionable and clinically complete.

Previous conversation:
{conversation_context if len(session['messages']) > 1 else "This is the start of the conversation."}

Current question: {raw_query}

Available research papers (use ONLY numbered citations like [1], [2], etc.):
{ref_list}

Paper details for context:
{context}

CRITICAL INSTRUCTIONS:
1. **Be maximally informative**: Include specific dosages (mg/kg, total doses, infusion rates), contraindications, side effects, patient risk factors, procedural steps, and monitoring parameters whenever relevant
2. **Emergency-ready responses**: If asked about acute situations (bronchospasm, hypotension, anaphylaxis, etc.), provide step-by-step management protocols with specific drugs and doses
3. **Never defer to external literature**: NEVER say "you might want to look into more literature" or "consider consulting guidelines" - extract ALL relevant information from the provided papers and present it
4. **Use numbered citations [1], [2], [3]**: Match the reference list exactly - NO author names in text
5. **Include clinical context**: Mention patient populations studied, contraindications, relative vs absolute benefits, NNT/NNH when available
6. **Be conversational but complete**: Natural tone like talking to a colleague, but don't sacrifice clinical detail for brevity
7. **Risk stratification**: When relevant, discuss patient-specific risk factors, ASA classification implications, comorbidity considerations
8. **HTML FORMATTING REQUIRED**: Use HTML tags for formatting:
   - Use <h3>Section Header</h3> for major sections
   - Use <strong>bold text</strong> for emphasis
   - Use <br><br> for paragraph breaks (double line breaks between sections)
   - Use <ul><li>item</li></ul> for bullet lists
   - Use <p>paragraph text</p> for paragraphs
   - Keep it well-structured and easy to scan

Example response for "What should I do for bronchospasm?":
"<h3>Acute Bronchospasm Management</h3>

<p>Here's the step-by-step approach [1][2]:</p>

<p><strong>1. Immediate Actions:</strong><br>
Deepen anesthesia with propofol 0.5-1 mg/kg bolus, or increase volatile to 2+ MAC [1]</p>

<p><strong>2. Bronchodilators:</strong><br>
Albuterol 4-8 puffs via ETT or 2.5mg nebulized [2]. Add ipratropium if severe [2]</p>

<p><strong>3. Steroids:</strong><br>
Methylprednisolone 1-2 mg/kg IV or hydrocortisone 100mg IV [1]</p>

<p><strong>4. Epinephrine:</strong><br>
If refractory - start with 10-20mcg IV boluses, titrate to effect. Consider infusion 0.01-0.05 mcg/kg/min [2][3]</p>

<p><strong>5. Ketamine:</strong><br>
0.5-1 mg/kg if resistant to above [3]</p>

<h3>Risk Factors</h3>
<p>Assess for: asthma history, recent URI, smoking, COPD [1]. Ensure adequate depth before any airway manipulation [2]</p>"

Respond with maximum clinical utility using HTML formatting:"""

        try:
            response = openai_client.chat.completions.create(
                model="gpt-4o",
                messages=[{"role": "user", "content": prompt}],
                temperature=0.1
            ).choices[0].message.content
        except Exception as e:
            response = f"AI error: {str(e)}"

        # Add assistant response to conversation
        session['messages'].append({
            "role": "assistant",
            "content": response,
            "references": refs,
            "num_papers": num_papers
        })
        session.modified = True

        return redirect(url_for('index'))

    return render_template_string(HTML, messages=session.get('messages', []))

@app.route("/clear")
def clear():
    """Clear conversation history"""
    session.pop('messages', None)
    return redirect(url_for('index'))

@app.route("/preop", methods=["GET", "POST"])
def preop_assessment():
    """Pre-operative assessment with evidence-based risk stratification"""
    if request.method == "GET":
        return render_template_string(PREOP_HTML, summary=None, references=None)

    # Collect form data
    age = int(request.form.get("age", 0))
    weight = float(request.form.get("weight", 0))
    height = float(request.form.get("height", 0))
    sex = request.form.get("sex", "")
    comorbidities = request.form.getlist("comorbidities")
    medications = request.form.get("medications", "")
    hgb = request.form.get("hgb", "")
    plt = request.form.get("plt", "")
    cr = request.form.get("cr", "")
    inr = request.form.get("inr", "")
    procedure = request.form.get("procedure", "")
    surgery_risk = request.form.get("surgery_risk", "")
    npo = request.form.get("npo", "")
    allergies = request.form.get("allergies", "")

    # Calculate BMI and IBW
    bmi = round(weight / ((height / 100) ** 2), 1) if weight and height else None
    if sex == 'male':
        ibw = round(50 + 0.91 * (height - 152.4), 1)
    else:
        ibw = round(45.5 + 0.91 * (height - 152.4), 1)

    # Build targeted PubMed searches based on patient risk factors
    search_queries = []

    # Anticoagulation management
    if any(drug in medications.lower() for drug in ['apixaban', 'rivaroxaban', 'warfarin', 'dabigatran', 'edoxaban', 'eliquis', 'xarelto', 'coumadin']):
        search_queries.append(f"perioperative anticoagulation management {procedure}")

    # Antiplatelet management
    if any(drug in medications.lower() for drug in ['aspirin', 'plavix', 'clopidogrel', 'ticagrelor', 'brilinta']):
        search_queries.append(f"perioperative antiplatelet management {procedure}")

    # Obesity + OSA
    if bmi and bmi >= 30 and "Obstructive Sleep Apnea" in comorbidities:
        search_queries.append("obese patient OSA perioperative anesthesia management")

    # Diabetes management
    if "Diabetes Mellitus" in comorbidities:
        search_queries.append(f"perioperative diabetes insulin management {procedure}")

    # Cardiac risk
    if any(c in comorbidities for c in ["Coronary Artery Disease", "Heart Failure", "Prior Stroke"]):
        search_queries.append(f"perioperative cardiac risk {procedure} guidelines")

    # CKD
    if "Chronic Kidney Disease" in comorbidities:
        search_queries.append(f"chronic kidney disease perioperative management {procedure}")

    # General perioperative guidelines for procedure
    search_queries.append(f"{procedure} anesthesia perioperative management guidelines")

    # Search PubMed for all queries and collect papers
    all_refs = []
    all_context = ""

    for query in search_queries[:3]:  # Limit to 3 searches to avoid overwhelming
        try:
            q_expanded = query.replace(" ", " AND ")
            search_term = (
                f'({q_expanded}) AND '
                f'(systematic review[pt] OR meta-analysis[pt] OR guideline[pt] OR '
                f'"randomized controlled trial"[pt]) AND '
                f'("2015/01/01"[PDAT] : "3000"[PDAT])'
            )

            handle = Entrez.esearch(db="pubmed", term=search_term, retmax=5, sort="relevance", api_key=Entrez.api_key)
            result = Entrez.read(handle)
            ids = result["IdList"]

            if ids:
                handle = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml", api_key=Entrez.api_key)
                papers = Entrez.read(handle)["PubmedArticle"]

                for p in papers:
                    try:
                        art = p["MedlineCitation"]["Article"]
                        title = art.get("ArticleTitle", "No title")
                        abstract = " ".join(str(t) for t in art.get("Abstract", {}).get("AbstractText", [])) if art.get("Abstract") else ""
                        authors = ", ".join([a.get("LastName","") + " " + (a.get("ForeName","")[:1]+"." if a.get("ForeName") else "") for a in art.get("AuthorList",[])[:3]])
                        journal = art["Journal"].get("Title", "Unknown")
                        year = art["Journal"]["JournalIssue"]["PubDate"].get("Year", "N/A")
                        pmid = p["MedlineCitation"]["PMID"]

                        all_refs.append({"title": title, "authors": authors, "journal": journal, "year": year, "pmid": pmid})
                        all_context += f"Title: {title}\nAbstract: {abstract}\nAuthors: {authors}\nJournal: {journal} ({year})\nPMID: {pmid}\n\n"
                    except:
                        continue
        except:
            continue

    # Remove duplicate references by PMID
    seen_pmids = set()
    unique_refs = []
    for ref in all_refs:
        if ref['pmid'] not in seen_pmids:
            seen_pmids.add(ref['pmid'])
            unique_refs.append(ref)

    # Create numbered reference list for GPT
    ref_list = ""
    for i, ref in enumerate(unique_refs, 1):
        ref_list += f"[{i}] {ref['title']} - {ref['authors']} ({ref['year']}) PMID: {ref['pmid']}\n"

    # Build patient summary for GPT
    patient_data = f"""
Patient Demographics:
- Age: {age} years
- Weight: {weight} kg
- Height: {height} cm
- Sex: {sex}
- BMI: {bmi} kg/m²
- IBW: {ibw} kg

Comorbidities: {', '.join(comorbidities) if comorbidities else 'None reported'}

Medications: {medications if medications else 'None reported'}

Laboratory Values:
- Hemoglobin: {hgb} g/dL
- Platelets: {plt} ×10³/μL
- Creatinine: {cr} mg/dL
- INR: {inr}

Procedure: {procedure}
Surgery Risk Category: {surgery_risk}

NPO Status: {npo if npo else 'Not specified'}
Allergies: {allergies if allergies else 'NKDA'}
"""

    # Generate GPT summary
    prompt = f"""You are an expert anesthesiologist performing a comprehensive pre-operative assessment. Based on the patient data and evidence from recent literature, provide a detailed evidence-based assessment.

Patient Information:
{patient_data}

Available Evidence (use numbered citations [1], [2], etc.):
{ref_list}

Paper Details:
{all_context}

Generate a comprehensive pre-operative assessment including:

1. **ASA Physical Status Classification**: Assign ASA class (I-V) with detailed justification based on comorbidities and functional status

2. **Cardiac Risk Stratification**: Calculate RCRI score if applicable (high-risk surgery + cardiac disease). Discuss perioperative cardiac risk.

3. **Perioperative Recommendations**:
   - Medication management (which to continue, hold, or adjust with specific timing)
   - Airway considerations (OSA, obesity, difficult airway predictors)
   - Hemodynamic management strategies
   - VTE prophylaxis recommendations
   - Glycemic control if diabetic
   - Renal protection if CKD

4. **Anesthetic Considerations**:
   - Preferred anesthetic technique with rationale
   - Drug selection and dosing adjustments
   - Monitoring requirements (standard vs advanced)
   - Postoperative disposition (PACU vs ICU)

5. **Evidence-Based Citations**: Use [1], [2], [3] format referencing the papers provided above

Use HTML formatting:
- <h3>Section Headers</h3>
- <p>Paragraphs</p>
- <strong>Bold for emphasis</strong>
- <br><br> for spacing

Provide maximum clinical utility with specific, actionable recommendations backed by evidence."""

    try:
        response = openai_client.chat.completions.create(
            model="gpt-4o",
            messages=[{"role": "user", "content": prompt}],
            temperature=0.1
        ).choices[0].message.content
    except Exception as e:
        response = f"<p>Error generating assessment: {str(e)}</p>"

    return render_template_string(PREOP_HTML, summary=response, references=unique_refs)

if __name__ == "__main__":
    app.run(debug=True)