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

HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>gasconsult.ai — Evidence-Based Anesthesiology</title>
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=2">
    <link rel="shortcut icon" type="image/x-icon" href="/static/favicon.ico?v=2">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=2">
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: -apple-system, BlinkMacSystemFont, 'SF Pro Display', 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            background-attachment: fixed;
            color: #1d1d1f;
            line-height: 1.6;
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
            overflow: hidden;
        }

        /* Navigation */
        nav {
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: saturate(180%) blur(20px);
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
            padding: 16px 0;
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            z-index: 1000;
            transition: all 0.3s ease;
        }

        nav .container {
            max-width: 1400px;
            margin: 0 auto;
            padding: 0 40px;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }

        .logo-text {
            font-size: 1.4rem;
            font-weight: 700;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            letter-spacing: -0.5px;
            transition: transform 0.2s ease;
        }

        .logo-text:hover {
            transform: scale(1.02);
        }

        .logo-symbol {
            display: inline-block;
            margin-right: 8px;
            font-size: 1.6rem;
            filter: grayscale(0);
        }

        .nav-actions {
            display: flex;
            gap: 12px;
            align-items: center;
        }

        .new-chat-btn {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 10px 24px;
            border-radius: 24px;
            font-size: 0.95rem;
            font-weight: 600;
            text-decoration: none;
            transition: all 0.3s ease;
            box-shadow: 0 4px 15px rgba(102, 126, 234, 0.3);
        }

        .new-chat-btn:hover {
            transform: translateY(-2px);
            box-shadow: 0 6px 20px rgba(102, 126, 234, 0.4);
        }

        /* Welcome Screen */
        .welcome-screen {
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            min-height: 100vh;
            padding: 80px 40px 40px;
            text-align: center;
            animation: fadeIn 0.6s ease;
        }

        @keyframes fadeIn {
            from { opacity: 0; transform: translateY(20px); }
            to { opacity: 1; transform: translateY(0); }
        }

        .welcome-screen .hero-logo {
            font-size: 4rem;
            margin-bottom: 24px;
            filter: drop-shadow(0 4px 20px rgba(255, 255, 255, 0.3));
            animation: float 3s ease-in-out infinite;
        }

        @keyframes float {
            0%, 100% { transform: translateY(0); }
            50% { transform: translateY(-10px); }
        }

        .welcome-screen h1 {
            font-size: 3.5rem;
            font-weight: 800;
            color: white;
            margin-bottom: 20px;
            letter-spacing: -1px;
            text-shadow: 0 2px 20px rgba(0, 0, 0, 0.2);
        }

        .welcome-screen .tagline {
            font-size: 1.4rem;
            color: rgba(255, 255, 255, 0.95);
            margin-bottom: 50px;
            font-weight: 500;
            max-width: 600px;
        }

        .welcome-screen .feature-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            max-width: 900px;
            margin-bottom: 40px;
        }

        .feature-card {
            background: rgba(255, 255, 255, 0.15);
            backdrop-filter: blur(10px);
            border: 1px solid rgba(255, 255, 255, 0.2);
            border-radius: 16px;
            padding: 24px;
            transition: all 0.3s ease;
        }

        .feature-card:hover {
            background: rgba(255, 255, 255, 0.25);
            transform: translateY(-4px);
        }

        .feature-card h3 {
            color: white;
            font-size: 1.1rem;
            margin-bottom: 8px;
            font-weight: 600;
        }

        .feature-card p {
            color: rgba(255, 255, 255, 0.9);
            font-size: 0.95rem;
        }

        /* Chat Container - Full Page */
        .chat-container {
            display: flex;
            flex-direction: column;
            height: 100vh;
            padding-top: 64px;
            max-width: 1000px;
            margin: 0 auto;
        }

        .chat-messages {
            flex: 1;
            overflow-y: auto;
            padding: 40px 40px 20px;
            scroll-behavior: smooth;
        }

        .chat-messages::-webkit-scrollbar {
            width: 8px;
        }

        .chat-messages::-webkit-scrollbar-track {
            background: transparent;
        }

        .chat-messages::-webkit-scrollbar-thumb {
            background: rgba(255, 255, 255, 0.3);
            border-radius: 10px;
        }

        .chat-messages::-webkit-scrollbar-thumb:hover {
            background: rgba(255, 255, 255, 0.5);
        }

        /* Chat Messages */
        .message {
            margin-bottom: 24px;
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
            max-width: 80%;
            padding: 18px 22px;
            border-radius: 20px;
            font-size: 1.05rem;
            line-height: 1.7;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
            transition: all 0.2s ease;
        }

        .message-content:hover {
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.15);
        }

        .message.user .message-content {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-bottom-right-radius: 6px;
        }

        .message.assistant .message-content {
            background: white;
            color: #1d1d1f;
            border-bottom-left-radius: 6px;
        }

        .message-text h3 {
            font-size: 1.25rem;
            margin-top: 16px;
            margin-bottom: 12px;
            color: #1d1d1f;
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
            margin-top: 16px;
            padding-top: 16px;
            border-top: 2px solid #f0f0f0;
            font-size: 0.95rem;
        }

        .message-refs strong {
            display: block;
            margin-bottom: 12px;
            color: #667eea;
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
            color: #667eea;
            text-decoration: none;
            font-weight: 500;
            transition: color 0.2s ease;
        }

        .ref-item a:hover {
            color: #764ba2;
            text-decoration: underline;
        }

        .message-meta {
            margin-top: 12px;
            font-size: 0.85rem;
            color: #888;
            font-weight: 600;
        }

        /* Chat Input */
        .chat-input-container {
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(20px);
            padding: 20px 40px 30px;
            border-top: 1px solid rgba(255, 255, 255, 0.3);
        }

        .chat-form {
            position: relative;
        }

        .chat-form textarea {
            width: 100%;
            padding: 18px 60px 18px 20px;
            font-size: 1.05rem;
            font-family: inherit;
            border: 2px solid rgba(102, 126, 234, 0.3);
            border-radius: 24px;
            resize: none;
            transition: all 0.3s ease;
            background: white;
            color: #1d1d1f;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.05);
        }

        .chat-form textarea:focus {
            outline: none;
            border-color: #667eea;
            box-shadow: 0 4px 20px rgba(102, 126, 234, 0.2);
        }

        .chat-form textarea::placeholder {
            color: #aaa;
        }

        .send-btn {
            position: absolute;
            right: 8px;
            bottom: 8px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border: none;
            padding: 12px 24px;
            border-radius: 18px;
            font-size: 0.95rem;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.3s ease;
            box-shadow: 0 2px 10px rgba(102, 126, 234, 0.3);
        }

        .send-btn:hover {
            transform: scale(1.05);
            box-shadow: 0 4px 15px rgba(102, 126, 234, 0.4);
        }

        .send-btn:active {
            transform: scale(0.98);
        }

        /* Responsive Design */
        @media (max-width: 768px) {
            .welcome-screen h1 {
                font-size: 2.5rem;
            }

            .welcome-screen .tagline {
                font-size: 1.1rem;
            }

            .welcome-screen .feature-grid {
                grid-template-columns: 1fr;
            }

            .message-content {
                max-width: 90%;
            }

            .chat-messages {
                padding: 20px;
            }

            .chat-input-container {
                padding: 15px 20px 20px;
            }

            nav .container {
                padding: 0 20px;
            }

            .send-btn {
                position: static;
                width: 100%;
                margin-top: 12px;
            }

            .chat-form textarea {
                padding: 16px 20px;
            }
        }

        /* Hide elements when chatting */
        .chatting .welcome-screen {
            display: none;
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
    </script>
</head>
<body class="{% if messages %}chatting{% endif %}">
    <nav>
        <div class="container">
            <div class="logo-text">
                <span class="logo-symbol">⚕</span>gasconsult.ai
            </div>
            <div class="nav-actions">
                {% if messages %}
                <a href="/clear" class="new-chat-btn">+ New Chat</a>
                {% endif %}
            </div>
        </div>
    </nav>

    {% if not messages %}
    <!-- Welcome Screen -->
    <div class="welcome-screen">
        <div class="hero-logo">⚕</div>
        <h1>gasconsult.ai</h1>
        <p class="tagline">Strictly Evidence-Based Anesthesiology Consults</p>

        <div class="feature-grid">
            <div class="feature-card">
                <h3>📚 PubMed-Backed</h3>
                <p>Every answer sourced from peer-reviewed research</p>
            </div>
            <div class="feature-card">
                <h3>🧮 Medical Calculators</h3>
                <p>MABL, IBW, BSA, QTc, and more</p>
            </div>
            <div class="feature-card">
                <h3>💬 Conversational</h3>
                <p>Ask follow-ups, refine your questions naturally</p>
            </div>
        </div>
    </div>
    {% endif %}

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
                                        {{ ref.title }} ({{ ref.year }})
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

        <div class="chat-input-container">
            <form method="post" action="/" class="chat-form">
                <textarea name="query" id="chatInput" placeholder="Ask anything about anesthesiology..." required rows="2"></textarea>
                <button type="submit" class="send-btn">Send</button>
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

        prompt = f"""You are an expert anesthesiologist AI assistant. You provide evidence-based answers and can also perform medical calculations conversationally.

Previous conversation:
{conversation_context if len(session['messages']) > 1 else "This is the start of the conversation."}

Current question: {raw_query}

I found the following research papers that may be relevant:

{context}

Your task:
1. If the user is asking a follow-up question related to previous conversation, reference that context naturally
2. Answer conversationally and concisely using the evidence provided
3. Cite sources by first author + year or PMID
4. If the user asked a calculation question but the papers aren't relevant, explain you couldn't find specific evidence but can help with the calculation if they provide values

Respond naturally as if having a conversation with a colleague.

Answer:"""

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

if __name__ == "__main__":
    app.run(debug=True)