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

def detect_and_calculate(query):
    """Detect calculation requests and perform medical calculations."""
    q = query.lower()

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
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg">
    <link rel="icon" type="image/x-icon" href="/static/favicon.ico">
    <link rel="apple-touch-icon" href="/static/favicon.svg">
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: -apple-system, BlinkMacSystemFont, 'SF Pro Display', 'SF Pro Text', 'Helvetica Neue', Arial, sans-serif;
            background: #fafafa;
            color: #1d1d1f;
            line-height: 1.6;
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
        }

        /* Navigation */
        nav {
            background: rgba(255, 255, 255, 0.8);
            backdrop-filter: saturate(180%) blur(20px);
            border-bottom: 1px solid rgba(0, 0, 0, 0.1);
            padding: 20px 0;
            position: sticky;
            top: 0;
            z-index: 1000;
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
            font-weight: 600;
            background: linear-gradient(135deg, #0071e3 0%, #00c3ff 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            letter-spacing: -0.5px;
        }

        .logo-symbol {
            display: inline-block;
            margin-right: 8px;
            font-size: 1.8rem;
        }

        /* Hero Section */
        .hero {
            max-width: 980px;
            margin: 80px auto 60px;
            padding: 0 40px;
            text-align: center;
        }

        h1 {
            font-size: 4rem;
            font-weight: 700;
            letter-spacing: -1.5px;
            margin-bottom: 20px;
            color: #1d1d1f;
            line-height: 1.1;
        }

        .subtitle {
            font-size: 1.5rem;
            color: #6e6e73;
            font-weight: 400;
            margin-bottom: 50px;
            line-height: 1.5;
        }

        /* Chat Container */
        .chat-container {
            max-width: 980px;
            margin: 0 auto 40px;
            padding: 0 40px;
            display: flex;
            flex-direction: column;
            height: calc(100vh - 500px);
            min-height: 400px;
        }

        .chat-messages {
            flex: 1;
            overflow-y: auto;
            padding: 20px;
            background: white;
            border-radius: 20px 20px 0 0;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
        }

        .welcome-message {
            text-align: center;
            padding: 60px 40px;
            color: #6e6e73;
        }

        .welcome-message h3 {
            font-size: 2rem;
            color: #1d1d1f;
            margin-bottom: 20px;
        }

        .welcome-message ul {
            text-align: left;
            max-width: 600px;
            margin: 30px auto;
            list-style: none;
            padding: 0;
        }

        .welcome-message li {
            padding: 12px 0;
            border-bottom: 1px solid #e5e5e7;
        }

        .welcome-message li:last-child {
            border-bottom: none;
        }

        /* Chat Messages */
        .message {
            margin-bottom: 24px;
            display: flex;
        }

        .message.user {
            justify-content: flex-end;
        }

        .message.assistant {
            justify-content: flex-start;
        }

        .message-content {
            max-width: 75%;
            padding: 16px 20px;
            border-radius: 18px;
            font-size: 1.05rem;
            line-height: 1.6;
        }

        .message.user .message-content {
            background: #0071e3;
            color: white;
            border-bottom-right-radius: 4px;
        }

        .message.assistant .message-content {
            background: #f5f5f7;
            color: #1d1d1f;
            border-bottom-left-radius: 4px;
        }

        .message-text {
            margin-bottom: 8px;
        }

        .message-text h3 {
            font-size: 1.3rem;
            margin-top: 16px;
            margin-bottom: 12px;
            color: #1d1d1f;
        }

        .message-text ul, .message-text ol {
            margin-left: 20px;
            margin-bottom: 12px;
        }

        .message-text p {
            margin-bottom: 12px;
        }

        .message-refs {
            margin-top: 12px;
            padding-top: 12px;
            border-top: 1px solid #e5e5e7;
            font-size: 0.95rem;
        }

        .message-refs strong {
            display: block;
            margin-bottom: 8px;
            color: #1d1d1f;
        }

        .ref-item {
            padding: 6px 0;
        }

        .ref-item a {
            color: #0071e3;
            text-decoration: none;
            font-weight: 500;
        }

        .ref-item a:hover {
            text-decoration: underline;
        }

        .message-meta {
            margin-top: 8px;
            font-size: 0.85rem;
            color: #6e6e73;
        }

        /* Chat Input */
        .chat-input-container {
            background: white;
            padding: 20px;
            border-radius: 0 0 20px 20px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
        }

        .chat-form textarea {
            width: 100%;
            padding: 12px 16px;
            font-size: 1.05rem;
            font-family: inherit;
            border: 2px solid #e5e5e7;
            border-radius: 12px;
            resize: none;
            transition: all 0.2s ease;
            background: #fafafa;
            color: #1d1d1f;
        }

        .chat-form textarea:focus {
            outline: none;
            border-color: #0071e3;
            background: white;
            box-shadow: 0 0 0 4px rgba(0, 113, 227, 0.1);
        }

        .chat-form textarea::placeholder {
            color: #86868b;
        }

        .chat-buttons {
            display: flex;
            gap: 12px;
            margin-top: 12px;
            justify-content: flex-end;
        }

        .send-btn, .clear-btn {
            padding: 12px 32px;
            border-radius: 980px;
            font-size: 1rem;
            font-weight: 500;
            cursor: pointer;
            transition: all 0.3s ease;
            text-decoration: none;
            display: inline-block;
        }

        .send-btn {
            background: #0071e3;
            color: white;
            border: none;
        }

        .send-btn:hover {
            background: #0077ed;
            transform: scale(1.02);
        }

        .clear-btn {
            background: #f5f5f7;
            color: #1d1d1f;
            border: 2px solid #e5e5e7;
        }

        .clear-btn:hover {
            background: #e8e8ed;
        }

        /* Footer */
        footer {
            max-width: 980px;
            margin: 80px auto 40px;
            padding: 40px 40px 0;
            text-align: center;
            border-top: 1px solid #d2d2d7;
            color: #86868b;
            font-size: 0.9rem;
        }

        /* Responsive Design */
        @media (max-width: 768px) {
            h1 {
                font-size: 2.5rem;
            }

            .subtitle {
                font-size: 1.2rem;
            }

            .hero {
                margin: 40px auto 40px;
            }

            .message-content {
                max-width: 90%;
            }

            .chat-container {
                padding: 0 20px;
                height: calc(100vh - 400px);
            }

            nav .container,
            .hero,
            footer {
                padding-left: 20px;
                padding-right: 20px;
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
    </script>
</head>
<body>
    <nav>
        <div class="container">
            <div class="logo-text">
                <span class="logo-symbol">⚕</span>gasconsult.ai
            </div>
        </div>
    </nav>

    <div class="hero">
        <h1>Evidence-Based<br>Anesthesiology Answers</h1>
        <p class="subtitle">Chat with an AI consultant backed by peer-reviewed research. No hallucinations, just evidence.</p>
    </div>

    <div class="chat-container">
        <div class="chat-messages" id="chatMessages">
            {% if messages %}
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
            {% else %}
                <div class="welcome-message">
                    <h3>👋 Welcome to gasconsult.ai</h3>
                    <p>Ask me anything about anesthesiology, or use my medical calculators:</p>
                    <ul>
                        <li><strong>Evidence questions:</strong> "TXA in spine surgery", "propofol vs etomidate in peds"</li>
                        <li><strong>Calculations:</strong> "MABL", "IBW", "BSA", "maintenance fluids", "QTc"</li>
                    </ul>
                    <p>I'll remember our conversation, so feel free to ask follow-up questions!</p>
                </div>
            {% endif %}
        </div>

        <div class="chat-input-container">
            <form method="post" action="/" class="chat-form">
                <textarea name="query" id="chatInput" placeholder="Ask a question or request a calculation..." required rows="2"></textarea>
                <div class="chat-buttons">
                    <button type="submit" class="send-btn">Send</button>
                    {% if messages %}
                    <a href="/clear" class="clear-btn">Clear Chat</a>
                    {% endif %}
                </div>
            </form>
        </div>
    </div>

    <footer>
        <p>Not medical advice • For educational use only<br>
        © 2025 gasconsult.ai</p>
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

        # Check if this is a calculation request first
        calc_result = detect_and_calculate(raw_query)
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

        # Build conversation context for GPT (last 5 exchanges to keep context manageable)
        conversation_context = ""
        recent_messages = session['messages'][-10:]  # Last 5 user + 5 assistant messages
        for msg in recent_messages:
            if msg['role'] == 'user':
                conversation_context += f"User: {msg['content']}\n"
            else:
                conversation_context += f"Assistant: {msg['content']}\n"

        prompt = f"""You are an expert anesthesiologist providing evidence-based consultations.

Previous conversation:
{conversation_context if len(session['messages']) > 1 else "This is the start of the conversation."}

Current question: {raw_query}

The references below are real, recent, high-quality papers (systematic reviews, meta-analyses, RCTs) that answer this clinical question:

References:
{context}

Answer concisely and directly using the evidence. If the question references previous conversation, use that context. Cite by first author + year or PMID.

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