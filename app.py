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

    # Remove trailing question marks FIRST so other patterns can match
    q = re.sub(r"[?!.]+$", "", q).strip()

    # Remove common conversational phrases
    filler_phrases = [
        r"^(can you |could you |please |pls )",
        r"^(tell me about |tell me |explain |describe )",
        r"^(what is |what's |what are |whats |what about )",
        r"^(what do you know about |what does the evidence say about )",
        r"^(what's the evidence for |what is the evidence for |evidence for )",
        r"^(i want to know about |i'd like to know about |i need to know about )",
        r"^(give me info on |give me information on |info on )",
        r"^(how does |how do |how is |how are |how about )",
        r"^(why does |why do |why is |why are )",
        r"^(search for |look up |find |show me )",
        r"^(help me understand |help me with )",
        r"(indications to use |indications for |indication for |when to use )",
        r"\s+(instead|as well|also|too)\s*$",  # Match with whitespace
    ]

    # Apply each pattern repeatedly until no more matches
    for pattern in filler_phrases:
        q = re.sub(pattern, "", q, flags=re.IGNORECASE).strip()

    # Remove extra whitespace
    q = re.sub(r"\s+", " ", q).strip()

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
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=3">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=3">
    <style>
        :root {
            --primary: #F97316;
            --primary-dark: #EA580C;
            --secondary: #3B82F6;
            --text-primary: #111827;
            --text-secondary: #4B5563;
            --text-muted: #9CA3AF;
            --background: #FFFFFF;
            --background-secondary: #F9FAFB;
            --border: #E5E7EB;
        }

        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'SF Pro Display', 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background: var(--background);
            color: #0A3D62;
            line-height: 1.6;
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
        }

        /* Navigation */
        nav {
            background: #FFFFFF;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.08);
            padding: 18px 40px;
            position: sticky;
            top: 0;
            left: 0;
            right: 0;
            z-index: 100;
            border-bottom: 1px solid #E5E7EB;
        }

        nav .container {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }

        .logo-text {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, system-ui, sans-serif;
            font-size: 24px;
            font-weight: 600;
            background: linear-gradient(135deg, #FF6B35 0%, #F97316 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            letter-spacing: -0.5px;
            transition: transform 0.2s ease;
            text-decoration: none;
            display: flex;
            align-items: center;
            gap: 10px;
            cursor: pointer;
        }

        .logo-text:hover {
            transform: scale(1.02);
        }

        .logo-symbol {
            font-size: 32px;
            background: linear-gradient(135deg, #FF6B35 0%, #F97316 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
        }

        .nav-actions {
            display: flex;
            gap: 12px;
            align-items: center;
        }

        .nav-link {
            color: #64748B;
            text-decoration: none;
            font-size: 15px;
            font-weight: 500;
            padding: 10px 20px;
            border-radius: 10px;
            transition: all 0.2s ease;
            background: transparent;
            cursor: pointer;
        }

        .nav-link:hover {
            background: #FFF5F0;
            color: #FF6B35;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.08);
            transform: translateY(-1px);
        }

        .nav-link.active {
            background: #FFF5F0;
            color: #FF6B35;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.08);
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
            color: #0066CC;
            margin-bottom: 12px;
            font-weight: 700;
        }

        .preop-header p {
            font-size: 1.1rem;
            color: #4B5563;
        }

        /* Form Sections */
        .form-section {
            background: rgba(255, 255, 255, 0.7);
            backdrop-filter: blur(20px);
            border-radius: 16px;
            padding: 28px;
            margin-bottom: 24px;
            box-shadow: 0 2px 12px rgba(0, 102, 204, 0.08);
            border: 1px solid rgba(0, 102, 204, 0.12);
        }

        .form-section h2 {
            color: #0066CC;
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
            color: #0066CC;
            margin-bottom: 6px;
        }

        input[type="text"],
        input[type="number"],
        textarea,
        select {
            padding: 12px;
            border: 2px solid rgba(0, 102, 204, 0.2);
            border-radius: 10px;
            font-size: 1rem;
            font-family: inherit;
            background: rgba(255, 255, 255, 0.9);
            color: #1F2937;
            transition: all 0.2s ease;
        }

        input:focus,
        textarea:focus,
        select:focus {
            outline: none;
            border-color: #0066CC;
            box-shadow: 0 0 0 3px rgba(0, 102, 204, 0.1);
            background: white;
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
            background: linear-gradient(135deg, #0066CC 0%, #0052A3 100%);
            color: white;
            padding: 14px 32px;
            border-radius: 12px;
            font-size: 1.05rem;
            font-weight: 600;
            border: none;
            cursor: pointer;
            transition: all 0.2s ease;
            width: 100%;
            margin-top: 20px;
            box-shadow: 0 4px 12px rgba(0, 102, 204, 0.25);
        }

        .submit-btn:hover {
            background: linear-gradient(135deg, #0052A3 0%, #003D7A 100%);
            transform: translateY(-2px);
            box-shadow: 0 6px 20px rgba(0, 102, 204, 0.35);
        }

        /* Summary Display */
        .summary-container {
            background: rgba(255, 255, 255, 0.8);
            backdrop-filter: blur(20px);
            border-radius: 16px;
            padding: 32px;
            margin-top: 40px;
            box-shadow: 0 4px 16px rgba(0, 102, 204, 0.1);
            border: 1px solid rgba(0, 102, 204, 0.15);
        }

        .summary-container h2 {
            color: #0066CC;
            font-size: 1.8rem;
            margin-bottom: 24px;
            font-weight: 700;
        }

        .summary-content {
            color: #1F2937;
            line-height: 1.8;
        }

        .summary-content h3 {
            color: #0066CC;
            font-size: 1.3rem;
            margin-top: 20px;
            margin-bottom: 12px;
            font-weight: 700;
        }

        .summary-content strong {
            color: #0052A3;
        }

        .ref-item {
            padding: 8px 0;
            transition: padding-left 0.2s ease;
        }

        .ref-item:hover {
            padding-left: 6px;
        }

        .ref-item a {
            color: #0066CC;
            text-decoration: none;
            font-weight: 500;
            transition: color 0.2s ease;
        }

        .ref-item a:hover {
            color: #0052A3;
            text-decoration: underline;
        }

        .auto-calc {
            background: rgba(0, 102, 204, 0.05);
            backdrop-filter: blur(10px);
            padding: 12px;
            border-radius: 10px;
            margin-top: 12px;
            border: 2px solid rgba(0, 102, 204, 0.15);
        }

        .auto-calc strong {
            color: #0066CC;
        }

        /* Footer */
        footer {
            text-align: center;
            padding: 50px 40px 40px;
            color: var(--text-muted);
            font-size: 0.875rem;
            border-top: 1px solid var(--border);
            margin-top: 80px;
            background: #FAFBFC;
        }

        footer p {
            margin-bottom: 8px;
            line-height: 1.6;
        }

        footer .disclaimer {
            max-width: 800px;
            margin: 16px auto 0;
            font-size: 0.8rem;
            color: #9CA3AF;
            line-height: 1.7;
        }

        footer a {
            color: var(--primary);
            text-decoration: none;
            font-weight: 500;
        }

        footer a:hover {
            text-decoration: underline;
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
                <span class="logo-symbol">⚕</span>
                <span>gasconsult.ai</span>
            </a>
            <div class="nav-actions">
                <a href="/" class="nav-link">Ask</a>
                <a href="/preop" class="nav-link active">Pre-Op Assessment</a>
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
                <div class="form-group" style="margin-top: 16px;">
                    <label for="other_comorbidities">Other Comorbidities (if not listed above)</label>
                    <textarea id="other_comorbidities" name="other_comorbidities" placeholder="e.g., GERD, Hypothyroidism, Chronic Pain..." rows="2"></textarea>
                </div>
            </div>

            <!-- Functional Status -->
            <div class="form-section">
                <h2>2b. Functional Status</h2>
                <div class="form-row">
                    <div class="form-group">
                        <label for="mets">Metabolic Equivalents (METs)</label>
                        <select id="mets" name="mets" required>
                            <option value="">Select...</option>
                            <option value="Unknown">Unknown / Not documented</option>
                            <option value="<4 METs">&lt;4 METs (Cannot climb 2 flights of stairs or walk 2 blocks)</option>
                            <option value="4-10 METs">4-10 METs (Can climb 2 flights of stairs)</option>
                            <option value=">10 METs">&gt;10 METs (Very active, can run or do strenuous sports)</option>
                        </select>
                    </div>
                </div>
            </div>

            <!-- Anesthesia History -->
            <div class="form-section">
                <h2>2c. Previous Anesthesia History</h2>
                <div class="form-group">
                    <label for="previous_anesthesia">Previous Anesthetics & Complications</label>
                    <textarea id="previous_anesthesia" name="previous_anesthesia" placeholder="e.g., General anesthesia for appendectomy 2015 - no complications. Family history of malignant hyperthermia..." rows="3"></textarea>
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

    <footer>
        <p>&copy; 2025 gasconsult.ai. All rights reserved.</p>
        <p class="disclaimer">For informational and educational purposes only. This tool is not intended to replace professional medical judgment and should only be used by qualified healthcare providers as a clinical decision support aid. All treatment decisions should be made in consultation with appropriate medical professionals.</p>
    </footer>

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
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=3">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=3">
    <style>
        :root {
            --primary: #F97316;
            --primary-dark: #EA580C;
            --secondary: #3B82F6;
            --text-primary: #111827;
            --text-secondary: #4B5563;
            --text-muted: #9CA3AF;
            --background: #FFFFFF;
            --background-secondary: #F9FAFB;
            --border: #E5E7EB;
        }

        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'SF Pro Display', 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background: var(--background);
            color: #0A3D62;
            line-height: 1.6;
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
        }

        /* Navigation */
        nav {
            background: #FFFFFF;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.08);
            padding: 18px 40px;
            position: sticky;
            top: 0;
            left: 0;
            right: 0;
            z-index: 100;
            border-bottom: 1px solid #E5E7EB;
        }

        nav .container {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }

        .logo-text {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, system-ui, sans-serif;
            font-size: 24px;
            font-weight: 600;
            background: linear-gradient(135deg, #FF6B35 0%, #F97316 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            letter-spacing: -0.5px;
            transition: transform 0.2s ease;
            text-decoration: none;
            display: flex;
            align-items: center;
            gap: 10px;
            cursor: pointer;
        }

        .logo-text:hover {
            transform: scale(1.02);
        }

        .logo-symbol {
            font-size: 32px;
            background: linear-gradient(135deg, #FF6B35 0%, #F97316 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
        }

        .nav-actions {
            display: flex;
            gap: 12px;
            align-items: center;
        }

        .nav-link {
            color: #64748B;
            text-decoration: none;
            font-size: 15px;
            font-weight: 500;
            padding: 10px 20px;
            border-radius: 10px;
            transition: all 0.2s ease;
            background: transparent;
            cursor: pointer;
        }

        .nav-link:hover {
            background: #FFF5F0;
            color: #FF6B35;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.08);
            transform: translateY(-1px);
        }

        .nav-link.active {
            background: #FFF5F0;
            color: #FF6B35;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.08);
        }

        .new-chat-btn {
            background: linear-gradient(135deg, #FF6B35 0%, #F97316 100%);
            color: white;
            padding: 10px 24px;
            border-radius: 8px;
            font-size: 15px;
            font-weight: 500;
            text-decoration: none;
            transition: all 0.2s ease;
            border: none;
            cursor: pointer;
            box-shadow: 0 2px 8px rgba(255, 107, 53, 0.2);
        }

        .new-chat-btn:hover {
            background: linear-gradient(135deg, #ff5722 0%, #EA580C 100%);
            transform: translateY(-1px);
            box-shadow: 0 4px 12px rgba(255, 107, 53, 0.35);
        }

        /* Welcome Screen */
        .welcome-screen {
            padding: 80px 40px 60px;
            margin: 0 auto;
            max-width: 1100px;
            text-align: center;
        }

        .hero-headline {
            font-size: 3.5rem;
            font-weight: 700;
            background: linear-gradient(135deg, #FF6B35 0%, #F97316 50%, #3B82F6 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            margin-bottom: 24px;
            letter-spacing: -1.5px;
            line-height: 1.1;
        }

        .hero-subtitle {
            font-size: 1.25rem;
            color: var(--text-secondary);
            margin-bottom: 48px;
            font-weight: 400;
            line-height: 1.6;
            max-width: 700px;
            margin-left: auto;
            margin-right: auto;
        }

        .preop-cta {
            display: inline-block;
            margin-bottom: 60px;
            padding: 14px 32px;
            background: linear-gradient(135deg, #FF6B35 0%, #F97316 100%);
            color: white;
            text-decoration: none;
            border-radius: 10px;
            font-weight: 600;
            font-size: 1rem;
            transition: all 0.3s ease;
            box-shadow: 0 4px 14px rgba(255, 107, 53, 0.25);
        }

        .preop-cta:hover {
            transform: translateY(-2px);
            box-shadow: 0 6px 20px rgba(255, 107, 53, 0.35);
        }

        .welcome-screen .feature-grid {
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 20px;
            margin: 0 auto;
        }

        .feature-card {
            background: rgba(255, 255, 255, 0.6);
            backdrop-filter: blur(10px);
            border-radius: 16px;
            padding: 32px 24px;
            transition: all 0.3s ease;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.04);
            border: 1px solid rgba(255, 107, 53, 0.1);
            text-align: center;
        }

        .feature-card:hover {
            transform: translateY(-4px);
            box-shadow: 0 12px 24px rgba(0, 0, 0, 0.08);
            border-color: rgba(255, 107, 53, 0.3);
        }

        .feature-icon {
            width: 56px;
            height: 56px;
            margin: 0 auto 20px;
            display: flex;
            align-items: center;
            justify-content: center;
            border-radius: 14px;
            position: relative;
        }

        .feature-icon.orange {
            background: linear-gradient(135deg, rgba(255, 107, 53, 0.1) 0%, rgba(249, 115, 22, 0.15) 100%);
            border: 2px solid rgba(255, 107, 53, 0.2);
        }

        .feature-icon.blue {
            background: linear-gradient(135deg, rgba(59, 130, 246, 0.1) 0%, rgba(37, 99, 235, 0.15) 100%);
            border: 2px solid rgba(59, 130, 246, 0.2);
        }

        .feature-icon svg {
            width: 28px;
            height: 28px;
        }

        .feature-card h3 {
            color: var(--text-primary);
            font-size: 1.1rem;
            margin-bottom: 10px;
            font-weight: 600;
        }

        .feature-card p {
            color: var(--text-secondary);
            font-size: 0.9rem;
            line-height: 1.6;
        }

        /* Main Content Area */
        .main-content {
            padding-top: 72px;
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
            max-width: 75%;
            padding: 12px 18px;
            border-radius: 18px;
            font-size: 0.95rem;
            line-height: 1.6;
            box-shadow: 0 1px 4px rgba(0, 0, 0, 0.08);
            transition: all 0.2s ease;
        }

        .message-content:hover {
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.12);
        }

        .message.user .message-content {
            background: linear-gradient(135deg, #0066CC 0%, #0052A3 100%);
            color: white;
            backdrop-filter: blur(10px);
            box-shadow: 0 2px 8px rgba(0, 102, 204, 0.25);
        }

        .message.assistant .message-content {
            background: rgba(255, 255, 255, 0.85);
            backdrop-filter: blur(20px);
            color: #1F2937;
            border: 1px solid rgba(0, 102, 204, 0.1);
            box-shadow: 0 1px 4px rgba(0, 0, 0, 0.06);
        }

        .message-text h3 {
            font-size: 1.15rem;
            margin-top: 16px;
            margin-bottom: 12px;
            color: #0066CC;
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
            border-top: 2px solid rgba(0, 102, 204, 0.15);
            font-size: 0.92rem;
        }

        .message-refs strong {
            display: block;
            margin-bottom: 12px;
            color: #0066CC;
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
            color: #0066CC;
            text-decoration: none;
            font-weight: 500;
            transition: color 0.2s ease;
        }

        .ref-item a:hover {
            color: #0052A3;
            text-decoration: underline;
        }

        .message-meta {
            margin-top: 12px;
            font-size: 0.85rem;
            color: #4B5563;
            opacity: 0.8;
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
            padding: 14px 56px 14px 20px;
            font-size: 1rem;
            font-family: inherit;
            border: 2px solid #e0e0e0;
            border-radius: 24px;
            resize: none;
            overflow: hidden;
            transition: all 0.2s ease;
            background: #ffffff;
            color: #0A3D62;
            height: 52px;
            line-height: 1.5;
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
            right: 8px;
            top: 6px;
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
            line-height: 1;
        }

        .send-btn:hover {
            background: #ff5722;
            transform: scale(1.05);
            box-shadow: 0 4px 12px rgba(255, 107, 53, 0.3);
        }

        .send-btn:active {
            transform: scale(0.95);
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

        /* Footer */
        footer {
            text-align: center;
            padding: 50px 40px 40px;
            color: var(--text-muted);
            font-size: 0.875rem;
            border-top: 1px solid var(--border);
            margin-top: 80px;
            background: #FAFBFC;
        }

        footer p {
            margin-bottom: 8px;
            line-height: 1.6;
        }

        footer .disclaimer {
            max-width: 800px;
            margin: 16px auto 0;
            font-size: 0.8rem;
            color: #9CA3AF;
            line-height: 1.7;
        }

        footer a {
            color: var(--primary);
            text-decoration: none;
            font-weight: 500;
        }

        footer a:hover {
            text-decoration: underline;
        }

        /* Responsive Design */
        @media (max-width: 768px) {
            .hero-section {
                padding: 80px 24px 60px;
            }

            .welcome-screen .hero-logo {
                font-size: 4rem;
            }

            .welcome-screen .tagline {
                font-size: 2rem;
            }

            .hero-subtitle {
                font-size: 1.1rem;
            }

            .features-section {
                padding: 60px 24px;
            }

            .welcome-screen .feature-grid {
                grid-template-columns: 1fr;
                gap: 24px;
            }

            .feature-card {
                padding: 32px 24px;
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
            <a href="/" class="logo-text">
                <span class="logo-symbol">⚕</span>
                <span>gasconsult.ai</span>
            </a>
            <div class="nav-actions">
                <a href="/" class="nav-link active">Ask</a>
                <a href="/preop" class="nav-link">Pre-Op Assessment</a>
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
            <!-- Hero Section -->
            <h1 class="hero-headline">Evidence-Based Anesthesiology Consultation</h1>
            <p class="hero-subtitle">Get instant, citation-backed clinical answers powered by PubMed research. Real evidence. Real citations. Zero hallucinations.</p>
            <a href="/preop" class="preop-cta">Pre-Operative Assessment Tool →</a>

            <!-- Features Section -->
            <div class="features-section">
                <div class="feature-grid">
                    <div class="feature-card">
                        <div class="feature-icon orange">
                            <svg viewBox="0 0 24 24" fill="none" stroke="#FF6B35" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path>
                                <path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path>
                            </svg>
                        </div>
                        <h3>PubMed-Backed Answers</h3>
                        <p>Every answer sourced from peer-reviewed research, systematic reviews, and clinical guidelines—with full citations you can verify.</p>
                    </div>
                    <div class="feature-card">
                        <div class="feature-icon blue">
                            <svg viewBox="0 0 24 24" fill="none" stroke="#3B82F6" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <rect x="4" y="2" width="16" height="20" rx="2"></rect>
                                <line x1="8" y1="6" x2="16" y2="6"></line>
                                <line x1="8" y1="10" x2="16" y2="10"></line>
                                <line x1="8" y1="14" x2="12" y2="14"></line>
                                <line x1="8" y1="18" x2="12" y2="18"></line>
                            </svg>
                        </div>
                        <h3>Medical Calculators</h3>
                        <p>Built-in calculators for MABL, IBW, BSA, QTc, maintenance fluids, and more. Just type your values and get instant results.</p>
                    </div>
                    <div class="feature-card">
                        <div class="feature-icon orange">
                            <svg viewBox="0 0 24 24" fill="none" stroke="#FF6B35" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <path d="M21 15a2 2 0 0 1-2 2H7l-4 4V5a2 2 0 0 1 2-2h14a2 2 0 0 1 2 2z"></path>
                                <line x1="9" y1="10" x2="15" y2="10"></line>
                                <line x1="12" y1="7" x2="12" y2="13"></line>
                            </svg>
                        </div>
                        <h3>Conversational AI</h3>
                        <p>Ask follow-up questions, refine your queries, and explore topics naturally—like talking to a colleague who knows the literature.</p>
                    </div>
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

    <footer>
        <p>&copy; 2025 gasconsult.ai. All rights reserved.</p>
        <p class="disclaimer">For informational and educational purposes only. This tool is not intended to replace professional medical judgment and should only be used by qualified healthcare providers as a clinical decision support aid. All treatment decisions should be made in consultation with appropriate medical professionals.</p>
    </footer>

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
        session.modified = True

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

        # IMPROVED FOLLOW-UP HANDLING: Detect follow-up questions
        # A follow-up is any question after the first exchange
        is_followup = len(session['messages']) >= 3  # At least one prior Q&A pair exists

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
        q = q.replace("barbiturates", '"barbiturates"[MeSH Terms] OR barbiturates OR thiopental')
        q = q.replace("barbs", '"barbiturates"[MeSH Terms] OR barbiturates OR thiopental')
        q = q.replace("ketamine", '"ketamine"[MeSH Terms] OR ketamine')
        q = q.replace("induction", '"anesthesia induction" OR induction OR "induction agent"')

        # Use broader search for follow-up questions - more likely to find results
        if is_followup:
            # For follow-ups, be very broad - include clinical trials and go back to 2005
            search_term = f'({q}) AND ("2005/01/01"[PDAT] : "3000"[PDAT])'
        else:
            # For initial questions, use high-quality evidence filters
            search_term = (
                f'({q}) AND '
                f'(systematic review[pt] OR meta-analysis[pt] OR "randomized controlled trial"[pt] OR '
                f'"Cochrane Database Syst Rev"[ta] OR guideline[pt]) AND '
                f'("2015/01/01"[PDAT] : "3000"[PDAT])'
            )

        # Try anesthesiology-specific search first
        try:
            handle = Entrez.esearch(db="pubmed", term=f'anesthesiology[MeSH Terms] AND {search_term}', retmax=15, sort="relevance", api_key=Entrez.api_key)
            result = Entrez.read(handle)
            ids = result["IdList"]
        except Exception as e:
            ids = []

        # Fallback 1: Try without anesthesiology restriction
        if not ids:
            try:
                handle = Entrez.esearch(db="pubmed", term=search_term, retmax=15, sort="relevance", api_key=Entrez.api_key)
                result = Entrez.read(handle)
                ids = result["IdList"]
            except:
                ids = []

        # Fallback 2: For initial queries, drop publication type restrictions
        if not ids and not is_followup:
            try:
                broader_search = f'({q}) AND ("2015/01/01"[PDAT] : "3000"[PDAT])'
                handle = Entrez.esearch(db="pubmed", term=broader_search, retmax=15, sort="relevance", api_key=Entrez.api_key)
                result = Entrez.read(handle)
                ids = result["IdList"]
            except:
                ids = []

        # If no papers found, handle gracefully
        if not ids:
            if is_followup:
                # For follow-ups with no new papers, use GPT with conversation context only
                conversation_context = ""
                recent_messages = session['messages'][-8:]  # Last 4 Q&A pairs
                for msg in recent_messages:
                    if msg['role'] == 'user':
                        conversation_context += f"User: {msg['content']}\n"
                    else:
                        content_text = re.sub('<[^<]+?>', '', msg.get('content', ''))
                        conversation_context += f"Assistant: {content_text[:400]}\n"

                prompt = f"""You are a clinical expert anesthesiologist AI assistant. The user is asking a follow-up question based on the conversation below.

Previous conversation:
{conversation_context}

Current follow-up question: {raw_query}

Provide a comprehensive, evidence-based answer that:
1. Builds naturally on the previous discussion
2. Includes specific clinical details (dosages, indications, contraindications, side effects)
3. Uses HTML formatting (<h3> for sections, <p> for paragraphs, <strong> for emphasis, <ul><li> for lists)
4. Is conversational but clinically complete
5. Notes that this draws from general anesthesiology knowledge and the previous discussion

Answer as if you're a colleague continuing the conversation:"""

                try:
                    response = openai_client.chat.completions.create(
                        model="gpt-4o",
                        messages=[{"role": "user", "content": prompt}],
                        temperature=0.2
                    ).choices[0].message.content

                    if not response or len(response.strip()) < 10:
                        response = "<p>I apologize, but I'm having trouble generating a response. Could you rephrase your question?</p>"

                except Exception as e:
                    response = f"<p>Error: Unable to generate response. Please try rephrasing your question.</p>"

                session['messages'].append({
                    "role": "assistant",
                    "content": response,
                    "references": [],
                    "num_papers": 0
                })
                session.modified = True
                return redirect(url_for('index'))
            else:
                # For initial queries with no results
                error_msg = "<p>No relevant evidence found in recent literature. Try rephrasing your question or using different medical terms.</p>"
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
{'NOTE: This is a FOLLOW-UP question - reference the previous discussion naturally and build on it.' if is_followup else ''}

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
    other_comorbidities = request.form.get("other_comorbidities", "")
    mets = request.form.get("mets", "")
    previous_anesthesia = request.form.get("previous_anesthesia", "")
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
    all_comorbidities = ', '.join(comorbidities) if comorbidities else 'None'
    if other_comorbidities:
        all_comorbidities += f"; {other_comorbidities}"

    patient_data = f"""
Patient Demographics:
- Age: {age} years
- Weight: {weight} kg
- Height: {height} cm
- Sex: {sex}
- BMI: {bmi} kg/m²
- IBW: {ibw} kg

Comorbidities: {all_comorbidities}

Functional Status:
- METs: {mets}

Previous Anesthesia History: {previous_anesthesia if previous_anesthesia else 'None reported'}

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

1. **ASA Physical Status Classification**: Assign ASA class (I-V) with detailed justification based on comorbidities and functional status (METs)

2. **Cardiac Risk Stratification**:
   - Calculate RCRI score if applicable (high-risk surgery + cardiac disease)
   - Reference ACS NSQIP Surgical Risk Calculator considerations for this patient's specific risk profile (age, comorbidities, functional status, procedure type)
   - Discuss perioperative cardiac risk with specific percentages when possible

3. **Perioperative Recommendations**:
   - Medication management (which to continue, hold, or adjust with specific timing)
   - Airway considerations (OSA, obesity, difficult airway predictors, previous anesthesia complications)
   - Hemodynamic management strategies
   - VTE prophylaxis recommendations
   - Glycemic control if diabetic
   - Renal protection if CKD
   - Special considerations based on previous anesthesia history

4. **Anesthetic Considerations**:
   - Preferred anesthetic technique with rationale
   - Drug selection and dosing adjustments
   - Monitoring requirements (standard vs advanced)
   - Postoperative disposition (PACU vs ICU)
   - Risk mitigation strategies

5. **Evidence-Based Citations**: Use [1], [2], [3] format referencing the papers provided above

Use HTML formatting:
- <h3>Section Headers</h3>
- <p>Paragraphs</p>
- <strong>Bold for emphasis</strong>
- <br><br> for spacing

Provide maximum clinical utility with specific, actionable recommendations backed by evidence. When discussing risk, reference the ACS NSQIP risk calculator framework and provide estimated risk percentages for major complications when relevant based on the patient's profile."""

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