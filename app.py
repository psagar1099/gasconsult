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

        .new-chat-btn {
            background: #FF6B35;
            color: white;
            padding: 10px 24px;
            border-radius: 10px;
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
            padding: 12px 100px 12px 20px;
            font-size: 1rem;
            font-family: inherit;
            border: 2px solid #e0e0e0;
            border-radius: 24px;
            resize: none;
            transition: all 0.2s ease;
            background: #ffffff;
            color: #0A3D62;
            height: 48px;
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
            right: 5px;
            top: 50%;
            transform: translateY(-50%);
            background: #FF6B35;
            color: white;
            border: none;
            padding: 9px 22px;
            border-radius: 18px;
            font-size: 0.95rem;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s ease;
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
                width: 100%;
                margin-top: 10px;
            }

            .send-btn:hover {
                transform: translateY(-1px);
            }

            .send-btn:active {
                transform: translateY(0);
            }

            .chat-form textarea {
                padding: 16px 20px;
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

Example response for "What should I do for bronchospasm?":
"For acute bronchospasm, here's the step-by-step approach [1][2]:

1. **Immediate**: Deepen anesthesia (propofol 0.5-1 mg/kg bolus, or increase volatile to 2+ MAC) [1]
2. **Bronchodilators**: Albuterol 4-8 puffs via ETT or 2.5mg nebulized [2], plus ipratropium if severe [2]
3. **Steroids**: Methylprednisolone 1-2 mg/kg IV or hydrocortisone 100mg IV [1]
4. **Epinephrine**: If refractory - start with 10-20mcg IV boluses, titrate to effect. Consider infusion 0.01-0.05 mcg/kg/min [2][3]
5. **Ketamine**: 0.5-1 mg/kg if resistant to above [3]

Risk factors to assess: asthma history, recent URI, smoking, COPD [1]. Ensure adequate depth before any airway manipulation [2]."

Respond with maximum clinical utility:"""

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