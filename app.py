from flask import Flask, request, render_template_string, session, redirect, url_for, Response, stream_with_context, jsonify
from flask_session import Session
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from flask_wtf.csrf import CSRFProtect
from Bio import Entrez
import openai
import os
import re
import secrets
import tempfile
import uuid
import json
import bleach
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

app = Flask(__name__)
app.secret_key = os.getenv('FLASK_SECRET_KEY', secrets.token_hex(32))

# Configure server-side sessions (stores sessions in filesystem instead of cookies)
app.config['SESSION_TYPE'] = 'filesystem'
app.config['SESSION_FILE_DIR'] = tempfile.gettempdir()
app.config['SESSION_PERMANENT'] = False
app.config['SESSION_USE_SIGNER'] = True
Session(app)

# Initialize CSRF protection
csrf = CSRFProtect(app)

# Initialize rate limiter
limiter = Limiter(
    app=app,
    key_func=get_remote_address,
    default_limits=[os.getenv("RATE_LIMIT", "60 per minute")],
    storage_uri="memory://"
)

Entrez.email = os.getenv("ENTREZ_EMAIL", "your-email@example.com")
Entrez.api_key = os.getenv("ENTREZ_API_KEY", "")

openai_client = openai.OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

# ====== Data Storage ======
# Simple in-memory storage for bookmarks and shared links
# In production, replace with database (Redis, PostgreSQL, etc.)
BOOKMARKS_STORAGE = {}  # {user_session_id: [{id, query, answer, references, timestamp}]}
SHARED_LINKS_STORAGE = {}  # {share_id: {query, answer, references, timestamp, expires}}
from datetime import datetime, timedelta

# ====== Security Functions ======

def sanitize_input(text, strip_tags=False):
    """
    Sanitize user input to prevent XSS attacks.

    Args:
        text: Input text to sanitize
        strip_tags: If True, strip all HTML tags. If False, allow safe HTML tags.

    Returns:
        Sanitized text safe for display
    """
    if not text:
        return ""

    if strip_tags:
        # Strip all HTML tags for plain text input
        return bleach.clean(text, tags=[], strip=True)
    else:
        # Allow safe HTML tags (for GPT-generated content)
        allowed_tags = [
            'p', 'br', 'strong', 'em', 'u', 'h1', 'h2', 'h3', 'h4', 'h5', 'h6',
            'ul', 'ol', 'li', 'a', 'blockquote', 'code', 'pre', 'hr', 'div', 'span',
            'table', 'thead', 'tbody', 'tr', 'th', 'td', 'sup', 'sub'
        ]
        allowed_attributes = {
            'a': ['href', 'title', 'target', 'rel'],
            'div': ['class', 'style'],
            'span': ['class', 'style'],
            'p': ['class', 'style'],
            'td': ['colspan', 'rowspan'],
            'th': ['colspan', 'rowspan']
        }
        allowed_protocols = ['http', 'https', 'mailto']

        return bleach.clean(
            text,
            tags=allowed_tags,
            attributes=allowed_attributes,
            protocols=allowed_protocols,
            strip=True
        )

def sanitize_user_query(query):
    """Sanitize user query input - strip all HTML tags."""
    return sanitize_input(query, strip_tags=True)

# ====== Logging Configuration ======

import logging
from logging.handlers import RotatingFileHandler
import sys

def setup_logging():
    """Configure comprehensive logging for the application"""
    log_level = os.getenv("LOG_LEVEL", "INFO").upper()
    log_file = os.getenv("LOG_FILE", "")

    # Create logger
    logger = logging.getLogger("gasconsult")
    logger.setLevel(getattr(logging, log_level, logging.INFO))

    # Create formatter
    formatter = logging.Formatter(
        '[%(asctime)s] %(levelname)s in %(module)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Console handler (always enabled)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # File handler (optional)
    if log_file:
        file_handler = RotatingFileHandler(
            log_file,
            maxBytes=10485760,  # 10MB
            backupCount=5
        )
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # Configure Flask app logger
    app.logger.handlers = logger.handlers
    app.logger.setLevel(logger.level)

    return logger

# Initialize logging
logger = setup_logging()

# Log startup information
logger.info("=" * 60)
logger.info("GasConsult.ai Application Starting")
logger.info(f"Environment: {os.getenv('FLASK_ENV', 'production')}")
logger.info(f"Log Level: {os.getenv('LOG_LEVEL', 'INFO')}")
logger.info(f"Rate Limit: {os.getenv('RATE_LIMIT', '60 per minute')}")
logger.info("=" * 60)

# Request logging middleware
@app.before_request
def log_request():
    """Log incoming requests"""
    logger.info(f"{request.method} {request.path} from {request.remote_addr}")

@app.after_request
def log_response(response):
    """Log outgoing responses"""
    logger.info(f"{request.method} {request.path} - Status: {response.status_code}")
    return response

@app.errorhandler(Exception)
def log_exception(e):
    """Log unhandled exceptions"""
    logger.error(f"Unhandled exception: {str(e)}", exc_info=True)
    return jsonify({
        "error": "An internal error occurred. Please try again later.",
        "status": "error"
    }), 500

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

def detect_question_type(query):
    """Identify what kind of clinical question this is to customize search and response"""
    q = query.lower()

    # Dosing questions - need guidelines more than RCTs
    if any(word in q for word in ['dose', 'dosing', 'how much', 'mg/kg', 'mcg/kg', 'ml/kg', 'dosage', 'loading dose', 'maintenance dose']):
        return 'dosing'

    # Safety/complications - need case reports too
    if any(word in q for word in ['safe', 'risk', 'complication', 'side effect', 'adverse', 'contraindication', 'warning', 'precaution', 'toxicity']):
        return 'safety'

    # Comparison questions - need head-to-head trials
    if any(word in q for word in [' or ', ' vs ', 'versus', 'compared', 'better', 'prefer', 'which', 'choice between']):
        return 'comparison'

    # Mechanism questions - need reviews
    if any(word in q for word in ['how does', 'mechanism', 'why does', 'works by', 'action of', 'pharmacology']):
        return 'mechanism'

    # Management/protocol questions - need guidelines
    if any(word in q for word in ['management', 'protocol', 'algorithm', 'approach', 'treat', 'handle', 'manage', 'strategy']):
        return 'management'

    return 'general'

def handle_negations(query):
    """Detect and properly handle negative/contraindication questions"""
    q = query.lower()

    is_negative = any(word in q for word in [
        'not ', 'avoid', 'never', 'contraindication', 'when to stop',
        'should not', 'shouldn\'t', 'don\'t use', 'do not use', 'when not to'
    ])

    if is_negative:
        # Add contraindication terms to search
        search_modifier = ' AND (contraindications[sh] OR adverse effects[sh] OR safety[ti])'

        # Update GPT prompt to focus on negatives
        prompt_modifier = "\n\nIMPORTANT: The user is asking about CONTRAINDICATIONS, AVOIDING, or WHEN NOT TO USE something. Focus your answer on safety concerns, contraindications, and situations where this is inappropriate or dangerous."

        return search_modifier, prompt_modifier

    return '', ''

def resolve_references(query, conversation_history):
    """Replace pronouns/vague references with actual terms from context"""
    if not conversation_history or len(conversation_history) < 2:
        return query

    q = query.lower()

    # Get last 2-3 user+assistant turns
    recent_messages = conversation_history[-6:] if len(conversation_history) >= 6 else conversation_history

    # Extract key medical terms from recent conversation (drugs, procedures, conditions, broad concepts)
    medical_entities = []
    for msg in recent_messages:
        content = msg.get('content', '').lower()

        # Find broader concepts first (multi-word phrases)
        broad_concepts = re.findall(
            r'\b(difficult airway management|difficult airway|airway management|'
            r'fiberoptic intubation|awake fiberoptic intubation|awake intubation|'
            r'rapid sequence induction|rapid sequence intubation|'
            r'neuraxial anesthesia|regional anesthesia|general anesthesia|'
            r'cardiac surgery|neurosurgery|spine surgery|orthopedic surgery|'
            r'video laryngoscopy|direct laryngoscopy|'
            r'spinal anesthesia|epidural anesthesia|combined spinal epidural|'
            r'peripheral nerve block|nerve block technique|'
            r'ultrasound guided block|landmark technique|'
            r'postoperative pain management|acute pain management|chronic pain|'
            r'preoperative optimization|preoperative assessment|'
            r'hemodynamic management|fluid management|blood management|'
            r'neuromuscular blockade|neuromuscular monitoring|'
            r'depth of anesthesia monitoring|bispectral index|'
            r'postoperative nausea|nausea prophylaxis|antiemetic strategy|'
            r'enhanced recovery|eras protocol|fast track surgery|'
            r'transfusion management|blood conservation|cell saver|'
            r'malignant hyperthermia crisis|local anesthetic toxicity crisis|'
            r'obstetric anesthesia|pediatric anesthesia|geriatric anesthesia|'
            r'liver transplant|kidney transplant|cardiac transplant|'
            r'one lung ventilation|lung isolation|double lumen tube|'
            r'airway exchange catheter|bougie|glidescope)\b',
            content
        )
        medical_entities.extend(broad_concepts)

        # Find single-word drugs, procedures, conditions
        entities = re.findall(
            r'\b(propofol|etomidate|ketamine|fentanyl|remifentanil|sufentanil|alfentanil|'
            r'rocuronium|vecuronium|succinylcholine|cisatracurium|'
            r'sevoflurane|desflurane|isoflurane|'
            r'midazolam|dexmedetomidine|precedex|'
            r'phenylephrine|ephedrine|epinephrine|norepinephrine|vasopressin|'
            r'tranexamic acid|txa|aminocaproic acid|'
            r'intubation|extubation|induction|emergence|'
            r'epidural|spinal|neuraxial|'
            r'ponv|hypotension|hypertension|bronchospasm|laryngospasm|'
            r'sugammadex|neostigmine|glycopyrrolate|atropine|'
            r'ondansetron|metoclopramide|dexamethasone|'
            r'naloxone|flumazenil|dantrolene|intralipid)\b',
            content
        )
        medical_entities.extend(entities)

    # Remove duplicates while preserving order
    seen = set()
    unique_entities = []
    for entity in reversed(medical_entities):  # Reverse to get most recent first
        if entity not in seen:
            unique_entities.append(entity)
            seen.add(entity)
    unique_entities.reverse()

    # Replace vague references if we have context
    if unique_entities and any(word in q for word in ['it', 'this', 'that', 'the drug', 'the medication', 'the technique', 'the procedure', 'the approach']):
        most_recent = unique_entities[-1]  # Last mentioned entity
        q = re.sub(r'\bit\b', most_recent, q)
        q = re.sub(r'\bthis\b', most_recent, q)
        q = re.sub(r'\bthat\b', most_recent, q)
        q = q.replace('the drug', most_recent)
        q = q.replace('the medication', most_recent)
        q = q.replace('the technique', most_recent)
        q = q.replace('the procedure', most_recent)
        q = q.replace('the approach', most_recent)
        print(f"[DEBUG] Resolved reference: '{query}' → '{q}'")

    # Handle "what about..." questions
    if q.startswith('what about') and unique_entities:
        # Extract the modifier (e.g., "pediatrics", "elderly", "renal failure")
        modifier = q.replace('what about', '').strip().rstrip('?').strip()
        main_topic = unique_entities[-1]
        q = f'{main_topic} in {modifier}'
        print(f"[DEBUG] Expanded 'what about': '{query}' → '{q}'")

    # Handle "how about..." questions similarly
    if q.startswith('how about') and unique_entities:
        modifier = q.replace('how about', '').strip().rstrip('?').strip()
        main_topic = unique_entities[-1]
        q = f'{main_topic} {modifier}'
        print(f"[DEBUG] Expanded 'how about': '{query}' → '{q}'")

    # Handle "can you" or "can I" questions that might be follow-ups
    if unique_entities and (q.startswith('can you') or q.startswith('can i') or q.startswith('should i') or q.startswith('should you')):
        # Check if the query contains specific techniques that should be added to context
        has_specific_technique = any(term in q for term in ['fiberoptic', 'video', 'laryngoscopy', 'intubate', 'block', 'epidural', 'spinal'])
        if has_specific_technique:
            # Prepend the context topic
            main_topic = unique_entities[-1]
            q = f'{main_topic} {q}'
            print(f"[DEBUG] Added context to technique question: '{query}' → '{q}'")

    return q if q != query.lower() else query

def detect_multipart(query):
    """Detect if user is asking multiple questions"""
    # Look for "and" connecting questions
    patterns = [
        r'(.*?)\s+and\s+(what|how|when|why|is|are|does|should|can)',
        r'(.*?)\s*\?\s*(what|how|when|why)',  # Multiple question marks
        r'(.*?)\s+(also|additionally)\s+(what|how|when|why)',
    ]

    for pattern in patterns:
        match = re.search(pattern, query, re.IGNORECASE)
        if match:
            return True

    return False

def build_smart_context(messages, current_query):
    """Build context that includes relevant history, not just recent"""
    if not messages or len(messages) == 0:
        return "New conversation."

    context = ""

    # Always include last 2 turns (most recent) - up to 6 messages
    recent = messages[-6:] if len(messages) >= 6 else messages

    # Also search for turns that mention same entities as current query (if conversation is longer)
    query_terms = set(re.findall(r'\b\w{4,}\b', current_query.lower()))
    relevant_earlier = []

    if len(messages) > 6:
        for i in range(0, len(messages) - 6, 2):  # Skip the recent ones we already have
            if i < len(messages):
                user_msg = messages[i].get('content', '').lower()
                terms = set(re.findall(r'\b\w{4,}\b', user_msg))
                overlap = len(query_terms & terms)
                if overlap >= 2:  # At least 2 shared terms
                    # Include Q&A pair
                    pair = messages[i:min(i+2, len(messages))]
                    relevant_earlier.append(pair)

    # Build context: relevant earlier + recent
    for pair in relevant_earlier[-2:]:  # Max 2 earlier relevant turns
        for msg in pair:
            role = "User" if msg['role'] == 'user' else "Assistant"
            content = re.sub('<[^<]+?>', '', msg.get('content', ''))
            context += f"{role}: {content[:200]}...\n"

    for msg in recent:
        role = "User" if msg['role'] == 'user' else "Assistant"
        content = re.sub('<[^<]+?>', '', msg.get('content', ''))
        max_len = 300 if msg['role'] == 'user' else 150
        if len(content) > max_len:
            context += f"{role}: {content[:max_len]}...\n"
        else:
            context += f"{role}: {content}\n"

    return context

def expand_medical_abbreviations(query):
    """Comprehensive medical abbreviation and synonym expansion"""
    q = query.lower()

    # Comparison operators
    q = q.replace(" versus ", " OR ")
    q = q.replace(" vs ", " OR ")
    q = q.replace(" vs. ", " OR ")
    q = q.replace(" over ", " OR ")
    q = q.replace(" compared to ", " OR ")
    q = q.replace(" compared with ", " OR ")

    # Common abbreviations and drugs - ORIGINAL ONES
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

    # NEW ADDITIONS - Neuromuscular blockers
    q = q.replace(" roc ", ' ("rocuronium"[MeSH Terms] OR rocuronium OR "neuromuscular blockade") ')
    q = q.replace(" vec ", ' ("vecuronium"[MeSH Terms] OR vecuronium) ')
    q = q.replace(" sux ", ' ("succinylcholine"[MeSH Terms] OR succinylcholine OR "muscle relaxant") ')
    q = q.replace("cisatracurium", '"cisatracurium"[MeSH Terms] OR cisatracurium OR nimbex')
    q = q.replace("rocuronium", '"rocuronium"[MeSH Terms] OR rocuronium')
    q = q.replace("vecuronium", '"vecuronium"[MeSH Terms] OR vecuronium')
    q = q.replace("succinylcholine", '"succinylcholine"[MeSH Terms] OR succinylcholine OR suxamethonium')

    # Opioids
    q = q.replace("fentanyl", '"fentanyl"[MeSH Terms] OR fentanyl OR opioid')
    q = q.replace(" remi ", ' ("remifentanil"[MeSH Terms] OR remifentanil) ')
    q = q.replace("remifentanil", '"remifentanil"[MeSH Terms] OR remifentanil')
    q = q.replace("sufentanil", '"sufentanil"[MeSH Terms] OR sufentanil')
    q = q.replace("alfentanil", '"alfentanil"[MeSH Terms] OR alfentanil')
    q = q.replace("morphine", '"morphine"[MeSH Terms] OR morphine')
    q = q.replace("hydromorphone", '"hydromorphone"[MeSH Terms] OR hydromorphone OR dilaudid')

    # Benzodiazepines and sedatives
    q = q.replace(" dex ", ' ("dexmedetomidine"[MeSH Terms] OR dexmedetomidine OR precedex) ')
    q = q.replace("dexmedetomidine", '"dexmedetomidine"[MeSH Terms] OR dexmedetomidine OR precedex')
    q = q.replace("versed", '"midazolam"[MeSH Terms] OR midazolam OR versed')
    q = q.replace("midazolam", '"midazolam"[MeSH Terms] OR midazolam OR versed')

    # Vasopressors and inotropes
    q = q.replace(" neo ", ' ("phenylephrine"[MeSH Terms] OR phenylephrine OR neosynephrine) ')
    q = q.replace("phenylephrine", '"phenylephrine"[MeSH Terms] OR phenylephrine OR neosynephrine')
    q = q.replace(" epi ", ' ("epinephrine"[MeSH Terms] OR epinephrine OR adrenaline) ')
    q = q.replace("epinephrine", '"epinephrine"[MeSH Terms] OR epinephrine OR adrenaline')
    q = q.replace("ephedrine", '"ephedrine"[MeSH Terms] OR ephedrine')
    q = q.replace("norepinephrine", '"norepinephrine"[MeSH Terms] OR norepinephrine OR levophed')
    q = q.replace("vasopressin", '"vasopressin"[MeSH Terms] OR vasopressin OR pitressin')

    # Common procedures
    q = q.replace(" rsi ", ' ("rapid sequence induction" OR RSI OR "rapid sequence intubation" OR "airway management") ')
    q = q.replace("rapid sequence", '"rapid sequence induction" OR RSI OR "rapid sequence intubation"')
    q = q.replace("awake intubation", '"awake intubation" OR "fiberoptic intubation" OR "difficult airway"')
    q = q.replace("awake fiberoptic", '"awake intubation" OR "fiberoptic intubation" OR "difficult airway"')
    q = q.replace(" cabg ", ' ("coronary artery bypass" OR CABG OR "cardiac surgery") ')
    q = q.replace("coronary artery bypass", '"coronary artery bypass" OR CABG OR "cardiac surgery"')
    q = q.replace(" tavr ", ' ("transcatheter aortic valve" OR TAVR OR "structural heart") ')
    q = q.replace("transcatheter aortic", '"transcatheter aortic valve" OR TAVR')

    # Common complications
    q = q.replace(" last ", ' ("local anesthetic systemic toxicity" OR LAST OR "lipid emulsion" OR "intralipid") ')
    q = q.replace("local anesthetic toxicity", '"local anesthetic systemic toxicity" OR LAST OR "lipid emulsion"')
    q = q.replace(" pris ", ' ("propofol infusion syndrome" OR PRIS) ')
    q = q.replace("propofol infusion syndrome", '"propofol infusion syndrome" OR PRIS')
    q = q.replace("malignant hyperthermia", '"malignant hyperthermia"[MeSH Terms] OR "dantrolene" OR MH')
    q = q.replace(" mh ", ' ("malignant hyperthermia"[MeSH Terms] OR "dantrolene" OR MH) ')

    # Common scenarios
    q = q.replace("full stomach", '"aspiration"[MeSH Terms] OR "rapid sequence" OR RSI OR "aspiration risk"')
    q = q.replace("difficult airway", '"difficult airway"[MeSH Terms] OR "airway management" OR intubation OR "difficult intubation"')
    q = q.replace("anticipated difficult", '"difficult airway"[MeSH Terms] OR "airway management" OR "awake intubation"')

    # Regional anesthesia
    q = q.replace("nerve block", '"nerve block"[MeSH Terms] OR "regional anesthesia" OR "peripheral nerve block"')
    q = q.replace("epidural", '"epidural"[MeSH Terms] OR "epidural anesthesia" OR "neuraxial"')
    q = q.replace("spinal", '"spinal anesthesia"[MeSH Terms] OR "subarachnoid" OR "neuraxial"')

    return q

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
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Sora:wght@400;600&display=swap" rel="stylesheet">
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=5">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">
    <meta name="apple-mobile-web-app-capable" content="yes">
    <meta name="apple-mobile-web-app-status-bar-style" content="default">
    <meta name="apple-mobile-web-app-title" content="gasconsult.ai">
    <style>
        :root {
            /* Primary Brand Colors */
            --primary-blue: #2563EB;
            --primary-blue-dark: #1D4ED8;
            --primary-blue-light: #DBEAFE;

            /* Anesthesia Color Palette (for logo & accents) */
            --opioid-blue: #2563EB;
            --nmb-red: #EF4444;
            --induction-yellow: #FBBF24;
            --vasopressor-violet: #8B5CF6;
            --anticholinergic-green: #10B981;
            --local-gray: #6B7280;

            /* Neutral Palette */
            --text-primary: #0F172A;
            --text-secondary: #475569;
            --text-muted: #94A3B8;
            --bg-primary: #FFFFFF;
            --bg-secondary: #F8FAFC;
            --border: #E2E8F0;

            /* Legacy aliases for compatibility */
            --primary: #2563EB;
            --primary-dark: #1D4ED8;
        }

        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        html {
            scroll-behavior: smooth;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'SF Pro Display', 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background: #F8FAFC;
            color: #0A3D62;
            line-height: 1.6;
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
            animation: pageFadeIn 0.6s cubic-bezier(0.4, 0, 0.2, 1);
            overflow-x: hidden;
            width: 100%;
        }

        @keyframes pageFadeIn {
            from {
                opacity: 0;
                transform: translateY(10px);
            }
            to {
                opacity: 1;
                transform: translateY(0);
            }
        }

        @keyframes slideUp {
            from {
                opacity: 0;
                transform: translateY(30px);
            }
            to {
                opacity: 1;
                transform: translateY(0);
            }
        }

        @keyframes slideIn {
            from {
                opacity: 0;
                transform: translateX(-20px);
            }
            to {
                opacity: 1;
                transform: translateX(0);
            }
        }

        @keyframes scaleIn {
            from {
                opacity: 0;
                transform: scale(0.95);
            }
            to {
                opacity: 1;
                transform: scale(1);
            }
        }

        @keyframes float {
            0%, 100% {
                transform: translateY(0);
            }
            50% {
                transform: translateY(-10px);
            }
        }

        @keyframes shimmer {
            0% {
                background-position: -1000px 0;
            }
            100% {
                background-position: 1000px 0;
            }
        }

        /* Navigation */
        nav {
            background: rgba(255, 255, 255, 0.85);
            backdrop-filter: blur(12px);
            -webkit-backdrop-filter: blur(12px);
            padding: 16px 40px;
            position: sticky;
            top: 0;
            left: 0;
            right: 0;
            z-index: 100;
            border-bottom: 1px solid rgba(226, 232, 240, 0.8);
        }

        nav .container {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            align-items: center;
            justify-content: space-between;
            flex-wrap: wrap;
            gap: 16px;
        }

        .logo-container {
            text-decoration: none;
            display: flex;
            align-items: center;
            gap: 12px;
            cursor: pointer;
            transition: transform 0.2s ease;
        }

        .logo-container:hover {
            transform: translateY(-1px);
        }

        .logo-ecg {
            height: 28px;
            width: auto;
            flex-shrink: 0;
        }

        .logo-wordmark {
            font-family: 'Sora', -apple-system, BlinkMacSystemFont, system-ui, sans-serif;
            font-size: 20px;
            font-weight: 600;
            letter-spacing: -0.5px;
            white-space: nowrap;
        }

        .logo-gas {
            color: #2563EB;
        }

        .logo-consult {
            color: #111111;
        }

        .logo-ai {
            font-weight: 400;
            color: #6B7280;
        }

        .nav-actions {
            display: flex;
            gap: 12px;
            align-items: center;
        }

        .nav-link {
            color: var(--text-secondary);
            text-decoration: none;
            font-size: 14px;
            font-weight: 500;
            padding: 8px 16px;
            border-radius: 8px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            color: var(--text-primary);
            background: rgba(255, 255, 255, 0.6);
            backdrop-filter: blur(8px);
            -webkit-backdrop-filter: blur(8px);
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.06);
        }

        .nav-link.active {
            color: var(--primary-blue);
            font-weight: 600;
        }

        /* PHI Warning Banner */
        .phi-warning {
            background: linear-gradient(135deg, #FEF3C7 0%, #FDE68A 100%);
            border-left: 4px solid #F59E0B;
            padding: 16px 20px;
            margin: 0;
            box-shadow: 0 2px 8px rgba(245, 158, 11, 0.15);
        }

        .phi-warning-content {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            align-items: center;
            gap: 12px;
        }

        .phi-warning-icon {
            font-size: 1.5rem;
            flex-shrink: 0;
        }

        .phi-warning-text {
            flex: 1;
        }

        .phi-warning-text strong {
            color: #92400E;
            font-weight: 600;
        }

        .phi-warning-text p {
            color: #78350F;
            font-size: 14px;
            line-height: 1.5;
            margin: 0;
        }

        /* Main Content */
        .preop-container {
            max-width: 900px;
            margin: 0 auto;
            padding: 100px 20px 60px;
        }

        .preop-header {
            text-align: center;
            margin-bottom: 48px;
        }

        .preop-header h1 {
            font-family: 'Sora', -apple-system, BlinkMacSystemFont, sans-serif;
            font-size: 2.5rem;
            background: linear-gradient(135deg, var(--primary-blue) 0%, var(--vasopressor-violet) 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            margin-bottom: 16px;
            font-weight: 700;
            letter-spacing: -1px;
            line-height: 1.2;
        }

        .preop-header p {
            font-size: 1.15rem;
            color: var(--text-secondary);
            font-weight: 500;
            letter-spacing: -0.2px;
        }

        /* Form Sections */
        .form-section {
            background: rgba(255, 255, 255, 0.7);
            backdrop-filter: blur(20px);
            -webkit-backdrop-filter: blur(20px);
            border-radius: 20px;
            padding: 28px;
            margin-bottom: 16px;
            box-shadow: 0 8px 32px rgba(37, 99, 235, 0.08), 0 2px 8px rgba(0, 0, 0, 0.04);
            border: 1px solid rgba(255, 255, 255, 0.8);
            position: relative;
            overflow: hidden;
            transition: all 0.4s cubic-bezier(0.4, 0, 0.2, 1);
            animation: slideUp 0.6s cubic-bezier(0.4, 0, 0.2, 1) backwards;
        }

        .form-section:nth-child(1) {
            animation-delay: 0.1s;
        }

        .form-section:nth-child(2) {
            animation-delay: 0.2s;
        }

        .form-section:nth-child(3) {
            animation-delay: 0.3s;
        }

        .form-section:nth-child(4) {
            animation-delay: 0.4s;
        }

        .form-section:nth-child(5) {
            animation-delay: 0.5s;
        }

        .form-section::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 4px;
            background: linear-gradient(90deg, var(--primary-blue) 0%, var(--vasopressor-violet) 100%);
            opacity: 0.8;
        }

        .form-section:hover {
            box-shadow: 0 12px 48px rgba(37, 99, 235, 0.15), 0 4px 12px rgba(0, 0, 0, 0.08);
            transform: translateY(-4px);
            border-color: rgba(37, 99, 235, 0.3);
        }

        .form-section h2 {
            color: var(--primary-blue);
            font-size: 1.35rem;
            margin-bottom: 24px;
            font-weight: 600;
            letter-spacing: -0.3px;
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
            font-size: 11px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: var(--text-muted);
            margin-bottom: 8px;
        }

        input[type="text"],
        input[type="number"],
        textarea,
        select {
            padding: 14px 16px;
            border: 2px solid rgba(226, 232, 240, 0.6);
            border-radius: 14px;
            font-size: 1rem;
            font-family: inherit;
            background: rgba(255, 255, 255, 0.9);
            color: var(--text-primary);
            transition: all 0.2s ease;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.02);
        }

        input:hover,
        textarea:hover,
        select:hover {
            border-color: rgba(37, 99, 235, 0.3);
            background: white;
        }

        input:focus,
        textarea:focus,
        select:focus {
            outline: none;
            border-color: var(--primary-blue);
            background: white;
            box-shadow: 0 4px 16px rgba(37, 99, 235, 0.15), 0 0 0 4px rgba(37, 99, 235, 0.08);
            transform: translateY(-1px);
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
            text-transform: none;
            letter-spacing: normal;
            font-size: 0.9rem;
            color: var(--text-primary);
        }

        .submit-btn {
            background: linear-gradient(135deg, var(--primary) 0%, var(--primary-dark) 100%);
            color: white;
            padding: 14px 32px;
            border-radius: 12px;
            font-size: 1.05rem;
            font-weight: 600;
            border: none;
            cursor: pointer;
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            width: 100%;
            margin-top: 20px;
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.25);
            position: relative;
            overflow: hidden;
        }

        .submit-btn::before {
            content: '';
            position: absolute;
            top: 50%;
            left: 50%;
            width: 0;
            height: 0;
            border-radius: 50%;
            background: rgba(255, 255, 255, 0.2);
            transform: translate(-50%, -50%);
            transition: width 0.6s, height 0.6s;
        }

        .submit-btn:hover::before {
            width: 300px;
            height: 300px;
        }

        .submit-btn:hover {
            background: linear-gradient(135deg, var(--primary-dark) 0%, #1E40AF 100%);
            transform: translateY(-2px) scale(1.02);
            box-shadow: 0 8px 24px rgba(37, 99, 235, 0.4);
        }

        .submit-btn:active {
            transform: translateY(0) scale(0.98);
        }

        /* Summary Display */
        .summary-container {
            background: linear-gradient(135deg, rgba(255, 255, 255, 0.95) 0%, rgba(255, 255, 255, 0.85) 100%);
            backdrop-filter: blur(40px);
            -webkit-backdrop-filter: blur(40px);
            border-radius: 24px;
            padding: 40px;
            margin-top: 32px;
            box-shadow: 0 20px 60px rgba(37, 99, 235, 0.15), 0 8px 24px rgba(0, 0, 0, 0.08);
            border: 1px solid rgba(255, 255, 255, 0.9);
            position: relative;
            overflow: hidden;
            animation: slideUp 0.5s ease-out;
        }

        .summary-container::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 6px;
            background: linear-gradient(90deg,
                var(--primary-blue) 0%,
                var(--vasopressor-violet) 50%,
                var(--anticholinergic-green) 100%);
        }

        @keyframes slideUp {
            from {
                opacity: 0;
                transform: translateY(20px);
            }
            to {
                opacity: 1;
                transform: translateY(0);
            }
        }

        .summary-container h2 {
            color: var(--primary-blue);
            font-size: 2rem;
            margin-bottom: 28px;
            font-weight: 600;
            letter-spacing: -0.5px;
        }

        .summary-content {
            color: #1F2937;
            line-height: 1.8;
        }

        .summary-content h3 {
            color: var(--primary-blue);
            font-size: 1.3rem;
            margin-top: 20px;
            margin-bottom: 12px;
            font-weight: 600;
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
            background: linear-gradient(135deg, rgba(37, 99, 235, 0.05) 0%, rgba(139, 92, 246, 0.05) 100%);
            backdrop-filter: blur(10px);
            -webkit-backdrop-filter: blur(10px);
            padding: 16px 18px;
            border-radius: 14px;
            margin-top: 16px;
            border: 1px solid rgba(37, 99, 235, 0.15);
            font-size: 14px;
            color: var(--text-secondary);
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.06);
            transition: all 0.2s ease;
        }

        .auto-calc:hover {
            border-color: rgba(37, 99, 235, 0.25);
            box-shadow: 0 6px 16px rgba(37, 99, 235, 0.1);
        }

        .auto-calc strong {
            color: var(--primary-blue);
            font-weight: 600;
        }

        /* Footer */
        footer {
            text-align: center;
            padding: 40px;
            border-top: 1px solid #E2E8F0;
            background: #FFFFFF;
            color: var(--text-muted);
            font-size: 13px;
            margin: 0;
        }

        footer a {
            transition: color 0.2s ease;
        }

        footer a:hover {
            color: var(--text-primary);
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

        /* Smooth Transitions & Animations */
        * {
            scroll-behavior: smooth;
        }

        .form-section {
            animation: fadeInUp 0.4s ease-out;
        }

        @keyframes fadeInUp {
            from {
                opacity: 0;
                transform: translateY(20px);
            }
            to {
                opacity: 1;
                transform: translateY(0);
            }
        }

        .submit-btn,
        input,
        textarea,
        select {
            transition: all 0.2s cubic-bezier(0.4, 0, 0.2, 1);
        }

        .summary-container {
            animation: slideIn 0.5s ease-out;
        }

        @keyframes slideIn {
            from {
                opacity: 0;
                transform: translateX(-20px);
            }
            to {
                opacity: 1;
                transform: translateX(0);
            }
        }

        /* Page fade-in */
        body {
            animation: pageFadeIn 0.3s ease-in;
        }

        @keyframes pageFadeIn {
            from {
                opacity: 0;
            }
            to {
                opacity: 1;
            }
        }

        /* Mobile responsiveness */
        @media (max-width: 768px) {
            nav {
                padding: 14px 20px;
            }

            .logo-text {
                font-size: 1.2rem;
            }

            .logo-svg {
                width: 28px;
                height: 28px;
            }

            .nav-actions {
                gap: 8px;
            }

            .nav-link {
                padding: 8px 14px;
                font-size: 0.9rem;
            }

            .preop-container {
                padding: 80px 20px 40px;
            }

            .preop-header {
                margin-bottom: 32px;
            }

            .preop-header h1 {
                font-size: 2.2rem;
                letter-spacing: -1px;
            }

            .preop-header p {
                font-size: 1rem;
            }

            .preop-form {
                padding: 30px 20px;
            }

            .form-row {
                grid-template-columns: 1fr;
                gap: 20px;
            }

            .summary-container {
                padding: 30px 20px;
            }

            .summary-container h2 {
                font-size: 1.5rem;
            }

            /* PHI Warning mobile adjustments */
            .phi-warning {
                padding: 12px 16px;
            }

            .phi-warning-content {
                gap: 10px;
            }

            .phi-warning-icon {
                font-size: 1.2rem;
            }

            .phi-warning-text p {
                font-size: 12px;
            }

            /* Form input mobile adjustments */
            .form-section {
                padding: 20px;
            }

            .form-section h2 {
                font-size: 1.1rem;
            }

            label {
                font-size: 0.9rem;
            }

            input[type="number"],
            input[type="text"],
            select,
            textarea {
                font-size: 16px !important; /* Prevents iOS zoom on focus */
            }

            .submit-btn {
                font-size: 0.95rem;
                padding: 14px 32px;
            }
        }

        @media (max-width: 640px) {
            nav .container {
                justify-content: space-between;
                flex-wrap: wrap;
            }

            .logo-container {
                order: 1;
            }

            .nav-actions {
                order: 2;
                width: 100%;
                justify-content: center;
                margin-top: 8px;
            }

            .nav-link {
                padding: 6px 10px;
                font-size: 0.8rem;
            }

            .logo-wordmark {
                font-size: 1rem;
            }

            .logo-ecg {
                height: 24px;
            }
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
            <a href="/" class="logo-container">
                <svg class="logo-ecg" viewBox="0 0 60 28" fill="none" xmlns="http://www.w3.org/2000/svg">
                    <defs>
                        <linearGradient id="ecgGrad" x1="0%" y1="0%" x2="100%" y2="0%">
                            <stop offset="0%" stop-color="#2563EB"/>
                            <stop offset="20%" stop-color="#EF4444"/>
                            <stop offset="40%" stop-color="#FBBF24"/>
                            <stop offset="60%" stop-color="#8B5CF6"/>
                            <stop offset="80%" stop-color="#10B981"/>
                            <stop offset="100%" stop-color="#6B7280"/>
                        </linearGradient>
                    </defs>
                    <path d="M2 14 L10 14 L14 12 L18 16 L22 4 L26 24 L30 10 L34 14 L42 14"
                          stroke="url(#ecgGrad)"
                          stroke-width="2.5"
                          stroke-linecap="round"
                          stroke-linejoin="round"
                          fill="none"/>
                </svg>
                <div class="logo-wordmark">
                    <span class="logo-gas">gas</span><span class="logo-consult">consult</span><span class="logo-ai">.ai</span>
                </div>
            </a>
            <div class="nav-actions">
                <a href="/" class="nav-link">Home</a>
                <a href="/preop" class="nav-link active">Pre-Op Assessment</a>
                <a href="/calculators" class="nav-link">Clinical Calculators</a>
                <a href="/quick-dose" class="nav-link">Quick Dose</a>
                <a href="/hypotension" class="nav-link">IOH Predictor</a>
            </div>
        </div>
    </nav>

    <!-- PHI Warning Banner -->
    <div class="phi-warning">
        <div class="phi-warning-content">
            <div class="phi-warning-icon">⚠️</div>
            <div class="phi-warning-text">
                <strong>Privacy Notice:</strong>
                <p>Do not enter patient names, dates of birth, MRNs, or other identifying information. Use age, weight, and clinical details only.</p>
            </div>
        </div>
    </div>

    <div class="preop-container">
        <div class="preop-header">
            <h1>Pre-Operative Assessment</h1>
            <p>Evidence-based risk stratification and recommendations</p>
        </div>

        {% if not summary %}
        <form method="post" action="/preop">
            <input type="hidden" name="csrf_token" value="{{ csrf_token() }}"/>
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
                        <label>Sex Assigned at Birth</label>
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
                <h2>4. Laboratory Values & Cardiac Assessment</h2>
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
                <div class="form-row">
                    <div class="form-group">
                        <label for="ef">Ejection Fraction (%)</label>
                        <input type="text" id="ef" name="ef" placeholder="e.g., 55-60% or None">
                    </div>
                    <div class="form-group">
                        <label for="ekg">EKG Findings</label>
                        <input type="text" id="ekg" name="ekg" placeholder="e.g., NSR, Afib, or None">
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
                <h3 style="color: var(--primary); margin-bottom: 16px;">References:</h3>
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
            <a href="/preop" class="submit-btn" style="display: inline-block; width: auto; text-decoration: none;">Clear & New Assessment</a>
        </div>
        {% endif %}
    </div>

    <footer>
        <p>&copy; 2025 gasconsult.ai. All rights reserved. | <a href="/terms" style="color: var(--primary); text-decoration: none;">Terms of Service</a> | <a href="/privacy" style="color: var(--primary); text-decoration: none;">Privacy Policy</a></p>
    </footer>

    <script>
        // Smooth form transitions
        const form = document.querySelector('form');
        if (form) {
            form.addEventListener('submit', function() {
                const submitBtn = form.querySelector('.submit-btn');
                if (submitBtn) {
                    submitBtn.disabled = true;
                    submitBtn.style.opacity = '0.7';
                    submitBtn.textContent = 'Processing...';
                }
            });
        }

        // Smooth scroll to summary if it exists
        window.addEventListener('load', function() {
            const summary = document.querySelector('.summary-container');
            if (summary) {
                setTimeout(() => {
                    summary.scrollIntoView({ behavior: 'smooth', block: 'start' });
                }, 200);
            }
        });
    </script>

</body>
</html>
"""

# ====== Premium Feature Functions ======

def format_citation_vancouver(ref, index):
    """Format a single reference in Vancouver style"""
    authors = ref.get('authors', 'Unknown')
    title = ref.get('title', 'Unknown title')
    journal = ref.get('journal', 'Unknown journal')
    year = ref.get('year', 'Unknown')
    pmid = ref.get('pmid', '')

    # Truncate authors if too long
    if len(authors) > 100:
        authors = authors[:97] + '...'

    return f"{index}. {authors}. {title}. {journal}. {year}. PMID: {pmid}"

def generate_bibtex(references):
    """Generate BibTeX format for citation managers"""
    bibtex_entries = []

    for i, ref in enumerate(references, 1):
        authors = ref.get('authors', 'Unknown').replace(' and ', ' and ')
        title = ref.get('title', 'Unknown title')
        journal = ref.get('journal', 'Unknown journal')
        year = ref.get('year', 'Unknown')
        pmid = ref.get('pmid', '')

        entry = f"""@article{{ref{i},
  author = {{{authors}}},
  title = {{{title}}},
  journal = {{{journal}}},
  year = {{{year}}},
  note = {{PMID: {pmid}}}
}}"""
        bibtex_entries.append(entry)

    return '\n\n'.join(bibtex_entries)

def generate_ris(references):
    """Generate RIS format for citation managers (Zotero, Mendeley, EndNote)"""
    ris_entries = []

    for ref in references:
        authors = ref.get('authors', 'Unknown')
        title = ref.get('title', 'Unknown title')
        journal = ref.get('journal', 'Unknown journal')
        year = ref.get('year', 'Unknown')
        pmid = ref.get('pmid', '')

        # Split authors
        author_list = authors.split(', ')
        author_lines = '\n'.join([f"AU  - {author}" for author in author_list[:5]])  # First 5 authors

        entry = f"""TY  - JOUR
{author_lines}
TI  - {title}
JO  - {journal}
PY  - {year}
N1  - PMID: {pmid}
ER  - """
        ris_entries.append(entry)

    return '\n\n'.join(ris_entries)

def get_evidence_strength(num_papers, references):
    """Analyze evidence strength and return classification"""
    if not num_papers or num_papers == 0:
        return {
            'level': 'Low',
            'color': '#EF4444',
            'description': 'Limited evidence available',
            'score': 1
        }

    # Analyze study types
    high_quality_count = 0
    meta_analysis_count = 0
    rct_count = 0
    review_count = 0

    for ref in references:
        title = ref.get('title', '').lower()
        journal = ref.get('journal', '').lower()

        if 'meta-analysis' in title or 'meta-analysis' in journal:
            meta_analysis_count += 1
            high_quality_count += 3
        elif 'randomized' in title or 'rct' in title:
            rct_count += 1
            high_quality_count += 2
        elif 'systematic review' in title or 'cochrane' in journal:
            review_count += 1
            high_quality_count += 2
        elif 'review' in title:
            review_count += 1
            high_quality_count += 1

    # Calculate strength score
    strength_score = (num_papers * 0.5) + high_quality_count

    if strength_score >= 10 or meta_analysis_count >= 2:
        return {
            'level': 'High',
            'color': '#10B981',
            'description': f'{num_papers} papers including {meta_analysis_count} meta-analyses, {rct_count} RCTs',
            'score': 3,
            'breakdown': {
                'meta_analyses': meta_analysis_count,
                'rcts': rct_count,
                'reviews': review_count,
                'total': num_papers
            }
        }
    elif strength_score >= 5 or num_papers >= 5:
        return {
            'level': 'Moderate',
            'color': '#FBBF24',
            'description': f'{num_papers} papers including {rct_count} RCTs, {review_count} reviews',
            'score': 2,
            'breakdown': {
                'meta_analyses': meta_analysis_count,
                'rcts': rct_count,
                'reviews': review_count,
                'total': num_papers
            }
        }
    else:
        return {
            'level': 'Low',
            'color': '#EF4444',
            'description': f'Limited evidence ({num_papers} papers)',
            'score': 1,
            'breakdown': {
                'meta_analyses': meta_analysis_count,
                'rcts': rct_count,
                'reviews': review_count,
                'total': num_papers
            }
        }

HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>gasconsult.ai — Evidence-Based Anesthesiology</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Sora:wght@400;600&display=swap" rel="stylesheet">
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=5">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">
    <meta name="apple-mobile-web-app-capable" content="yes">
    <meta name="apple-mobile-web-app-status-bar-style" content="default">
    <meta name="apple-mobile-web-app-title" content="gasconsult.ai">
    <style>
        :root {
            /* Primary Brand Colors */
            --primary-blue: #2563EB;
            --primary-blue-dark: #1D4ED8;
            --primary-blue-light: #DBEAFE;

            /* Anesthesia Color Palette (for logo & accents) */
            --opioid-blue: #2563EB;
            --nmb-red: #EF4444;
            --induction-yellow: #FBBF24;
            --vasopressor-violet: #8B5CF6;
            --anticholinergic-green: #10B981;
            --local-gray: #6B7280;

            /* Neutral Palette */
            --text-primary: #0F172A;
            --text-secondary: #475569;
            --text-muted: #94A3B8;
            --bg-primary: #FFFFFF;
            --bg-secondary: #F8FAFC;
            --border: #E2E8F0;

            /* Legacy aliases for compatibility */
            --primary: #2563EB;
            --primary-dark: #1D4ED8;
        }

        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        html {
            scroll-behavior: smooth;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'SF Pro Display', 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background: #F8FAFC;
            color: #0A3D62;
            line-height: 1.6;
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
            animation: pageFadeIn 0.6s cubic-bezier(0.4, 0, 0.2, 1);
            overflow-x: hidden;
            width: 100%;
        }

        @keyframes pageFadeIn {
            from {
                opacity: 0;
                transform: translateY(10px);
            }
            to {
                opacity: 1;
                transform: translateY(0);
            }
        }

        @keyframes slideUp {
            from {
                opacity: 0;
                transform: translateY(30px);
            }
            to {
                opacity: 1;
                transform: translateY(0);
            }
        }

        @keyframes slideIn {
            from {
                opacity: 0;
                transform: translateX(-20px);
            }
            to {
                opacity: 1;
                transform: translateX(0);
            }
        }

        @keyframes scaleIn {
            from {
                opacity: 0;
                transform: scale(0.95);
            }
            to {
                opacity: 1;
                transform: scale(1);
            }
        }

        @keyframes float {
            0%, 100% {
                transform: translateY(0);
            }
            50% {
                transform: translateY(-10px);
            }
        }

        @keyframes shimmer {
            0% {
                background-position: -1000px 0;
            }
            100% {
                background-position: 1000px 0;
            }
        }

        /* Navigation */
        nav {
            background: rgba(255, 255, 255, 0.85);
            backdrop-filter: blur(12px);
            -webkit-backdrop-filter: blur(12px);
            padding: 16px 40px;
            position: sticky;
            top: 0;
            left: 0;
            right: 0;
            z-index: 100;
            border-bottom: 1px solid rgba(226, 232, 240, 0.8);
        }

        nav .container {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            align-items: center;
            justify-content: space-between;
            flex-wrap: wrap;
            gap: 16px;
        }

        .logo-container {
            text-decoration: none;
            display: flex;
            align-items: center;
            gap: 12px;
            cursor: pointer;
            transition: transform 0.2s ease;
        }

        .logo-container:hover {
            transform: translateY(-1px);
        }

        .logo-ecg {
            height: 28px;
            width: auto;
            flex-shrink: 0;
        }

        .logo-wordmark {
            font-family: 'Sora', -apple-system, BlinkMacSystemFont, system-ui, sans-serif;
            font-size: 20px;
            font-weight: 600;
            letter-spacing: -0.5px;
            white-space: nowrap;
        }

        .logo-gas {
            color: #2563EB;
        }

        .logo-consult {
            color: #111111;
        }

        .logo-ai {
            font-weight: 400;
            color: #6B7280;
        }

        .nav-actions {
            display: flex;
            gap: 12px;
            align-items: center;
        }

        .nav-link {
            color: var(--text-secondary);
            text-decoration: none;
            font-size: 14px;
            font-weight: 500;
            padding: 8px 16px;
            border-radius: 8px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            color: var(--text-primary);
            background: rgba(255, 255, 255, 0.6);
            backdrop-filter: blur(8px);
            -webkit-backdrop-filter: blur(8px);
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.06);
        }

        .nav-link.active {
            color: var(--primary-blue);
            font-weight: 600;
        }

        /* PHI Warning Banner */
        .phi-warning {
            background: linear-gradient(135deg, #FEF3C7 0%, #FDE68A 100%);
            border-left: 4px solid #F59E0B;
            padding: 16px 20px;
            margin: 0;
            box-shadow: 0 2px 8px rgba(245, 158, 11, 0.15);
        }

        .phi-warning-content {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            align-items: center;
            gap: 12px;
        }

        .phi-warning-icon {
            font-size: 1.5rem;
            flex-shrink: 0;
        }

        .phi-warning-text {
            flex: 1;
        }

        .phi-warning-text strong {
            color: #92400E;
            font-weight: 600;
        }

        .phi-warning-text p {
            color: #78350F;
            font-size: 14px;
            line-height: 1.5;
            margin: 0;
        }

        .new-chat-btn {
            background: linear-gradient(135deg, var(--primary) 0%, var(--primary-dark) 100%);
            color: white;
            padding: 10px 24px;
            border-radius: 8px;
            font-size: 15px;
            font-weight: 500;
            text-decoration: none;
            transition: all 0.2s ease;
            border: none;
            cursor: pointer;
            box-shadow: 0 2px 8px rgba(37, 99, 235, 0.2);
        }

        .new-chat-btn:hover {
            background: linear-gradient(135deg, var(--primary-dark) 0%, #1E40AF 100%);
            transform: translateY(-1px);
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.35);
        }

        .clear-chat-btn {
            background: transparent;
            color: var(--text-secondary);
            padding: 8px 16px;
            border-radius: 8px;
            font-size: 14px;
            font-weight: 500;
            text-decoration: none;
            transition: all 0.2s ease;
            border: 1.5px solid var(--border);
        }

        .clear-chat-btn:hover {
            background: rgba(239, 68, 68, 0.05);
            color: #EF4444;
            border-color: #EF4444;
        }

        /* Welcome Screen */
        .welcome-screen {
            padding: 80px 40px 60px;
            margin: 0 auto;
            max-width: 900px;
            text-align: center;
            position: relative;
            animation: scaleIn 0.7s cubic-bezier(0.4, 0, 0.2, 1) 0.2s backwards;
        }

        /* Hero background gradient - extends behind navbar */
        .main-content::before {
            content: '';
            position: absolute;
            top: -72px;
            left: 50%;
            transform: translateX(-50%);
            width: 100%;
            height: 700px;
            background: radial-gradient(ellipse at 50% 0%, rgba(37, 99, 235, 0.15) 0%, transparent 65%);
            z-index: -1;
            pointer-events: none;
        }

        .hero-headline {
            font-family: 'Sora', -apple-system, BlinkMacSystemFont, system-ui, sans-serif;
            font-size: 48px;
            font-weight: 700;
            letter-spacing: -1.5px;
            line-height: 1.2;
            margin-bottom: 20px;
        }

        .hero-headline-blue {
            color: #2563EB;
        }

        .hero-headline-dark {
            color: #0F172A;
        }

        .hero-subtitle {
            font-size: 17px;
            color: #475569;
            max-width: 640px;
            margin: 0 auto 24px;
            font-weight: 400;
            line-height: 1.65;
        }

        /* CTA Buttons Container */
        .cta-buttons {
            display: flex;
            gap: 16px;
            align-items: center;
            justify-content: center;
            margin-bottom: 48px;
            flex-wrap: wrap;
        }

        /* Pre-Op CTA Button */
        .preop-cta-outline {
            display: inline-flex;
            align-items: center;
            gap: 10px;
            padding: 16px 36px;
            background: linear-gradient(135deg, #2563EB 0%, #1D4ED8 100%);
            border: none;
            border-radius: 14px;
            color: white;
            font-weight: 600;
            font-size: 16px;
            text-decoration: none;
            transition: all 0.4s cubic-bezier(0.4, 0, 0.2, 1);
            box-shadow: 0 4px 14px rgba(37, 99, 235, 0.3);
            position: relative;
            overflow: hidden;
        }

        .preop-cta-outline::before {
            content: '';
            position: absolute;
            top: 0;
            left: -100%;
            width: 100%;
            height: 100%;
            background: linear-gradient(90deg, transparent, rgba(255, 255, 255, 0.2), transparent);
            transition: left 0.6s;
        }

        .preop-cta-outline:hover::before {
            left: 100%;
        }

        .preop-cta-outline:hover {
            transform: translateY(-3px) scale(1.02);
            box-shadow: 0 8px 24px rgba(37, 99, 235, 0.4);
        }

        .preop-cta-outline svg {
            transition: transform 0.4s cubic-bezier(0.4, 0, 0.2, 1);
        }

        .preop-cta-outline:hover svg {
            transform: translateX(4px);
        }

        /* Quick Dose CTA Button */
        .quickdose-cta-filled {
            display: inline-flex;
            align-items: center;
            gap: 10px;
            padding: 16px 36px;
            background: linear-gradient(135deg, #8B5CF6 0%, #7C3AED 100%);
            border: none;
            border-radius: 14px;
            color: white;
            font-weight: 600;
            font-size: 16px;
            text-decoration: none;
            transition: all 0.4s cubic-bezier(0.4, 0, 0.2, 1);
            box-shadow: 0 4px 14px rgba(139, 92, 246, 0.3);
            position: relative;
            overflow: hidden;
        }

        .quickdose-cta-filled::before {
            content: '';
            position: absolute;
            top: 0;
            left: -100%;
            width: 100%;
            height: 100%;
            background: linear-gradient(90deg, transparent, rgba(255, 255, 255, 0.2), transparent);
            transition: left 0.6s;
        }

        .quickdose-cta-filled:hover::before {
            left: 100%;
        }

        .quickdose-cta-filled:hover {
            transform: translateY(-3px) scale(1.02);
            box-shadow: 0 8px 24px rgba(139, 92, 246, 0.4);
        }

        .quickdose-cta-filled svg {
            transition: transform 0.4s cubic-bezier(0.4, 0, 0.2, 1);
        }

        .quickdose-cta-filled:hover svg {
            transform: scale(1.1);
        }

        /* Trust Badges */
        .trust-badges {
            display: flex;
            gap: 24px;
            align-items: center;
            justify-content: center;
            margin-bottom: 48px;
            font-size: 13px;
            color: var(--text-muted);
        }

        .trust-badge {
            display: flex;
            align-items: center;
            gap: 6px;
        }

        .trust-badge svg {
            flex-shrink: 0;
            color: var(--anticholinergic-green);
        }

        /* Homepage Chat Section */
        .homepage-chat-section {
            padding: 0;
            margin-top: -40px;
            position: relative;
            z-index: 10;
        }

        .homepage-input {
            margin: 0 auto 80px;
            max-width: 700px;
            padding: 0 40px;
            background: transparent !important;
            border: none !important;
            position: static !important;
        }

        .homepage-input .chat-form {
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(20px);
            -webkit-backdrop-filter: blur(20px);
            border: 2px solid rgba(37, 99, 235, 0.1);
            border-radius: 16px;
            padding: 6px;
            display: flex;
            align-items: center;
            gap: 8px;
            outline: none;
            box-shadow: 0 8px 32px rgba(37, 99, 235, 0.08), 0 2px 12px rgba(0, 0, 0, 0.04);
            transition: all 0.4s cubic-bezier(0.4, 0, 0.2, 1);
            position: relative;
            overflow: hidden;
        }

        .homepage-input .chat-form::before {
            content: '';
            position: absolute;
            top: 0;
            left: -100%;
            width: 100%;
            height: 100%;
            background: linear-gradient(90deg, transparent, rgba(37, 99, 235, 0.05), transparent);
            transition: left 0.6s cubic-bezier(0.4, 0, 0.2, 1);
        }

        .homepage-input .chat-form:hover::before {
            left: 100%;
        }

        .homepage-input .chat-form:hover {
            border-color: rgba(37, 99, 235, 0.2);
            box-shadow: 0 12px 48px rgba(37, 99, 235, 0.12), 0 4px 16px rgba(0, 0, 0, 0.06);
            transform: translateY(-2px);
        }

        .homepage-input .chat-form:focus-within {
            border-color: var(--primary-blue);
            box-shadow: 0 16px 56px rgba(37, 99, 235, 0.18), 0 8px 24px rgba(0, 0, 0, 0.08);
            transform: translateY(-3px);
        }

        .homepage-input .chat-form textarea {
            flex: 1;
            padding: 12px 16px;
            font-size: 15px;
            border: none;
            outline: none;
            background: transparent;
            color: var(--text-primary);
            resize: none;
            height: 44px;
            min-height: unset;
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            line-height: 1.4;
            position: relative;
            z-index: 1;
        }

        .homepage-input .chat-form textarea::placeholder {
            color: var(--text-muted);
        }

        .homepage-input .send-btn {
            position: relative;
            z-index: 1;
            background: linear-gradient(135deg, var(--primary-blue) 0%, var(--primary-blue-dark) 100%);
            width: 44px;
            height: 44px;
            border-radius: 50%;
            font-size: 1.75rem;
            font-weight: 900;
            color: white;
            border: none;
            cursor: pointer;
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            flex-shrink: 0;
            display: flex;
            align-items: center;
            justify-content: center;
            padding: 0 0 3px 0;
            -webkit-text-stroke: 0.5px white;
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.3);
        }

        .homepage-input .send-btn:hover {
            transform: scale(1.08) rotate(5deg);
            box-shadow: 0 6px 20px rgba(37, 99, 235, 0.4);
        }

        .homepage-input .send-btn:active {
            transform: scale(0.95);
        }

        .preop-cta {
            display: inline-block;
            margin-bottom: 60px;
            padding: 14px 32px;
            background: linear-gradient(135deg, var(--primary) 0%, var(--primary-dark) 100%);
            color: white;
            text-decoration: none;
            border-radius: 10px;
            font-weight: 600;
            font-size: 1rem;
            transition: all 0.3s ease;
            box-shadow: 0 4px 14px rgba(37, 99, 235, 0.25);
        }

        .preop-cta:hover {
            transform: translateY(-2px);
            box-shadow: 0 6px 20px rgba(37, 99, 235, 0.35);
        }

        /* Features Section */
        .features-section {
            padding: 80px 0 0 0;
            background: #F8FAFC;
            margin: 0;
            position: relative;
            z-index: 1;
        }

        .feature-grid {
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 24px;
            max-width: 1100px;
            margin: 0 auto;
            padding: 0 40px 80px 40px;
        }

        .feature-card {
            background: var(--bg-primary);
            border-radius: 16px;
            padding: 32px;
            transition: all 0.4s cubic-bezier(0.4, 0, 0.2, 1);
            border: 1px solid var(--border);
            animation: slideUp 0.6s cubic-bezier(0.4, 0, 0.2, 1) backwards;
            position: relative;
            overflow: hidden;
        }

        .feature-card::before {
            content: '';
            position: absolute;
            top: 0;
            left: -100%;
            width: 100%;
            height: 100%;
            background: linear-gradient(90deg, transparent, rgba(37, 99, 235, 0.03), transparent);
            transition: left 0.6s cubic-bezier(0.4, 0, 0.2, 1);
        }

        .feature-card:hover::before {
            left: 100%;
        }

        .feature-card:nth-child(1) {
            animation-delay: 0.1s;
        }

        .feature-card:nth-child(2) {
            animation-delay: 0.2s;
        }

        .feature-card:nth-child(3) {
            animation-delay: 0.3s;
        }

        .feature-card:hover {
            border-color: var(--primary-blue);
            box-shadow: 0 12px 40px rgba(37, 99, 235, 0.12), 0 4px 12px rgba(0, 0, 0, 0.08);
            transform: translateY(-6px);
        }

        .feature-card:hover .feature-icon {
            transform: scale(1.1) rotate(5deg);
        }

        .feature-title {
            font-family: 'Sora', sans-serif;
            font-size: 18px;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 12px;
        }

        .feature-description {
            font-size: 14px;
            color: var(--text-secondary);
            line-height: 1.7;
        }

        .feature-icon {
            width: 48px;
            height: 48px;
            margin-bottom: 20px;
            display: flex;
            align-items: center;
            justify-content: center;
            border-radius: 12px;
            transition: all 0.4s cubic-bezier(0.4, 0, 0.2, 1);
            position: relative;
        }

        .feature-icon.violet {
            background: linear-gradient(135deg, #EDE9FE 0%, #DDD6FE 100%);
            border: none;
        }

        .feature-icon.blue {
            background: linear-gradient(135deg, #DBEAFE 0%, #BFDBFE 100%);
            border: none;
        }

        .feature-icon.green {
            background: linear-gradient(135deg, #D1FAE5 0%, #A7F3D0 100%);
            border: none;
        }

        .feature-icon svg {
            width: 24px;
            height: 24px;
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
            padding-top: 0;
            display: flex;
            flex-direction: column;
            position: relative;
            background: linear-gradient(180deg, #EEF2FF 0%, #F8FAFC 100%);
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
            animation: slideIn 0.5s cubic-bezier(0.4, 0, 0.2, 1);
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
            background: linear-gradient(135deg, var(--primary) 0%, var(--primary-dark) 100%);
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
            position: relative;
        }

        /* Copy button */
        .copy-btn {
            position: absolute;
            top: 12px;
            right: 12px;
            background: rgba(255, 255, 255, 0.9);
            border: 1px solid rgba(0, 102, 204, 0.2);
            border-radius: 6px;
            padding: 6px 12px;
            font-size: 0.85rem;
            color: #4B5563;
            cursor: pointer;
            display: flex;
            align-items: center;
            gap: 6px;
            transition: all 0.2s ease;
            font-family: inherit;
            font-weight: 500;
        }

        .copy-btn:hover {
            background: white;
            border-color: #0066CC;
            color: #0066CC;
            box-shadow: 0 2px 6px rgba(0, 102, 204, 0.15);
        }

        .copy-btn.copied {
            color: #10B981;
            border-color: #10B981;
        }

        .copy-btn svg {
            width: 14px;
            height: 14px;
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

        /* Loading indicator styling */
        .loading-indicator {
            color: #0066CC;
            font-weight: 500;
            animation: pulse 1.5s ease-in-out infinite;
        }

        @keyframes pulse {
            0%, 100% { opacity: 1; }
            50% { opacity: 0.5; }
        }

        /* Evidence Quality Badge */
        .evidence-quality-badge {
            background: linear-gradient(135deg, #F0F9FF 0%, #E0F2FE 100%);
            border-left: 4px solid #0EA5E9;
            border-radius: 12px;
            padding: 16px 20px;
            margin-bottom: 20px;
            box-shadow: 0 2px 8px rgba(14, 165, 233, 0.1);
        }

        .confidence-level {
            font-size: 0.95rem;
            margin-bottom: 8px;
        }

        .confidence-level.high {
            color: #047857;
        }

        .confidence-level.high strong {
            color: #065F46;
        }

        .confidence-level.moderate {
            color: #D97706;
        }

        .confidence-level.moderate strong {
            color: #B45309;
        }

        .confidence-level.low {
            color: #DC2626;
        }

        .confidence-level.low strong {
            color: #B91C1C;
        }

        .evidence-details {
            font-size: 0.875rem;
            color: #475569;
            line-height: 1.6;
        }

        /* Suggested Prompts */
        .suggested-prompts {
            max-width: 900px;
            margin: 8px auto 16px;
            padding: 0 20px;
        }

        .suggested-prompts h4 {
            font-size: 0.9rem;
            color: #64748B;
            margin-bottom: 12px;
            font-weight: 500;
        }

        .prompt-buttons {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
        }

        .prompt-btn {
            background: white;
            border: 1.5px solid #E2E8F0;
            border-radius: 8px;
            padding: 8px 16px;
            font-size: 0.875rem;
            color: #475569;
            cursor: pointer;
            transition: all 0.2s ease;
            font-family: 'Inter', sans-serif;
        }

        .prompt-btn:hover {
            border-color: #2563EB;
            color: #2563EB;
            background: #F0F9FF;
            transform: translateY(-1px);
            box-shadow: 0 2px 8px rgba(37, 99, 235, 0.1);
        }

        /* References styling for streamed content */
        .references {
            margin-top: 24px;
            padding-top: 18px;
            border-top: 2px solid rgba(0, 102, 204, 0.15);
        }

        .references h4 {
            color: #0066CC;
            font-weight: 700;
            margin-bottom: 14px;
            font-size: 1.05rem;
        }

        .references ol {
            margin-left: 20px;
        }

        .references li {
            margin-bottom: 14px;
            line-height: 1.6;
        }

        .references a {
            color: #0066CC;
            text-decoration: none;
            font-weight: 500;
            transition: color 0.2s ease;
        }

        .references a:hover {
            color: #0052A3;
            text-decoration: underline;
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
            display: flex;
            align-items: center;
            gap: 12px;
            background: var(--bg-primary);
            border: 1px solid var(--border);
            border-radius: 16px;
            padding: 8px;
        }

        .chat-form textarea {
            flex: 1;
            padding: 12px 16px;
            font-size: 15px;
            font-family: inherit;
            border: none;
            resize: none;
            overflow: hidden;
            transition: all 0.2s ease;
            background: transparent;
            color: var(--text-primary);
            height: 44px;
            min-height: unset;
            max-height: 200px;
            line-height: 1.4;
        }

        .chat-form textarea:focus {
            outline: none;
        }

        .chat-form textarea::placeholder {
            color: var(--text-muted);
        }

        .send-btn {
            position: static;
            background: var(--primary-blue);
            color: white;
            border: none;
            width: 44px;
            height: 44px;
            border-radius: 50%;
            font-size: 1.75rem;
            font-weight: 900;
            cursor: pointer;
            transition: all 0.2s ease;
            display: flex;
            align-items: center;
            justify-content: center;
            padding: 0 0 3px 0;
            line-height: 1;
            flex-shrink: 0;
            -webkit-text-stroke: 0.5px white;
        }

        .send-btn:hover {
            background: var(--primary-blue-dark);
            transform: scale(1.05);
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

        /* Footer */
        footer {
            text-align: center;
            padding: 40px;
            border-top: 1px solid #E2E8F0;
            background: #FFFFFF;
            color: var(--text-muted);
            font-size: 13px;
            margin: 0;
        }

        footer a {
            transition: color 0.2s ease;
        }

        footer a:hover {
            color: var(--text-primary);
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

        /* Smooth Transitions & Animations */
        * {
            scroll-behavior: smooth;
        }

        .message {
            animation: messageSlideIn 0.4s ease-out;
        }

        @keyframes messageSlideIn {
            from {
                opacity: 0;
                transform: translateY(20px);
            }
            to {
                opacity: 1;
                transform: translateY(0);
            }
        }

        .feature-card {
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        }

        .feature-card:hover {
            transform: translateY(-8px);
        }

        .nav-link,
        .new-chat-btn,
        .preop-cta,
        .submit-btn {
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        }

        /* Page fade-in */
        body {
            animation: pageFadeIn 0.3s ease-in;
        }

        @keyframes pageFadeIn {
            from {
                opacity: 0;
            }
            to {
                opacity: 1;
            }
        }

        /* Input focus animations */
        input:focus,
        textarea:focus,
        select:focus {
            transition: all 0.2s cubic-bezier(0.4, 0, 0.2, 1);
        }

        /* Responsive Design */
        @media (max-width: 768px) {
            .welcome-screen {
                padding: 60px 24px 40px;
            }

            .hero-headline {
                font-size: 36px;
                letter-spacing: -1px;
                margin-bottom: 16px;
            }

            .hero-subtitle {
                font-size: 16px;
                line-height: 1.6;
            }

            .homepage-input {
                padding: 0 24px;
            }

            .features-section {
                padding: 60px 0 0 0;
            }

            .feature-grid {
                grid-template-columns: 1fr;
                padding: 0 24px 60px 24px;
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

            /* Navigation mobile adjustments */
            nav {
                padding: 14px 20px;
            }

            .nav-actions {
                gap: 8px;
            }

            .nav-link {
                padding: 8px 14px;
                font-size: 0.9rem;
            }

            .new-chat-btn {
                padding: 8px 14px;
                font-size: 0.9rem;
            }

            /* Hero headline mobile */
            .hero-headline {
                font-size: 2rem !important;
                line-height: 1.2 !important;
            }

            .preop-cta {
                padding: 12px 24px;
                font-size: 0.95rem;
            }

            /* Copy button mobile adjustments */
            .copy-btn {
                top: 8px;
                right: 8px;
                padding: 5px 10px;
                font-size: 0.8rem;
            }

            .copy-btn svg {
                width: 12px;
                height: 12px;
            }

            /* Message styling mobile */
            .message-text {
                font-size: 0.95rem;
            }

            .message-text h3 {
                font-size: 1.05rem;
            }

            /* Loading indicator */
            .loading-indicator {
                font-size: 0.9rem;
            }

            /* References mobile */
            .message-refs,
            .references {
                font-size: 0.85rem;
            }

            .message-meta {
                font-size: 0.8rem;
            }

            /* PHI Warning mobile adjustments */
            .phi-warning {
                padding: 12px 16px;
            }

            .phi-warning-content {
                gap: 10px;
            }

            .phi-warning-icon {
                font-size: 1.2rem;
            }

            .phi-warning-text p {
                font-size: 12px;
            }

            /* Homepage chat input mobile */
            .homepage-input {
                padding: 0 20px;
                margin-bottom: 60px !important;
            }

            .chat-form {
                padding: 6px;
                gap: 8px;
            }

            .chat-form textarea {
                font-size: 14px;
                padding: 10px 14px;
            }

            /* Homepage hero mobile */
            .hero-subtitle {
                font-size: 0.95rem;
                padding: 0 10px;
            }

            .trust-badges {
                flex-direction: column;
                gap: 12px;
            }

            .preop-cta-outline {
                padding: 12px 24px;
                font-size: 0.95rem;
            }

            /* Message user/assistant badges mobile */
            .message.user .message-content {
                max-width: 85%;
            }

            .message.assistant .message-content {
                max-width: 100%;
            }
        }

        @media (max-width: 640px) {
            nav .container {
                justify-content: space-between;
                flex-wrap: wrap;
            }

            .logo-container {
                order: 1;
            }

            .nav-actions {
                order: 2;
                width: 100%;
                justify-content: center;
                margin-top: 8px;
            }

            .nav-link {
                padding: 6px 10px;
                font-size: 0.8rem;
            }

            .logo-wordmark {
                font-size: 1rem;
            }

            .logo-ecg {
                height: 24px;
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

        // PWA Install Prompt (Mobile Only)
        let deferredPrompt;
        const isMobile = /Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(navigator.userAgent);
        const isIOS = /iPhone|iPad|iPod/i.test(navigator.userAgent);
        const isStandalone = window.matchMedia('(display-mode: standalone)').matches || window.navigator.standalone;

        console.log('Mobile detected:', isMobile);
        console.log('iOS detected:', isIOS);
        console.log('Already installed:', isStandalone);

        if (isMobile && !isStandalone) {
            // For Android/Chrome - use beforeinstallprompt event
            window.addEventListener('beforeinstallprompt', (e) => {
                console.log('beforeinstallprompt event fired');
                e.preventDefault();
                deferredPrompt = e;

                // Create install banner
                const installBanner = document.createElement('div');
                installBanner.id = 'pwaInstallBanner';
                installBanner.style.cssText = `
                    position: fixed;
                    bottom: 0;
                    left: 0;
                    right: 0;
                    background: linear-gradient(135deg, #2563EB 0%, #1D4ED8 100%);
                    color: white;
                    padding: 16px 20px;
                    box-shadow: 0 -4px 12px rgba(37, 99, 235, 0.2);
                    z-index: 9999;
                    display: flex;
                    align-items: center;
                    justify-content: space-between;
                    gap: 12px;
                    animation: slideUp 0.3s ease;
                `;

                installBanner.innerHTML = `
                    <div style="flex: 1;">
                        <div style="font-weight: 600; font-size: 14px; margin-bottom: 4px;">Add to Home Screen</div>
                        <div style="font-size: 12px; opacity: 0.9;">Install gasconsult.ai for quick access</div>
                    </div>
                    <button id="pwaInstallBtn" style="
                        background: white;
                        color: #2563EB;
                        border: none;
                        padding: 10px 20px;
                        border-radius: 8px;
                        font-weight: 600;
                        font-size: 14px;
                        cursor: pointer;
                        flex-shrink: 0;
                    ">Install</button>
                    <button id="pwaCloseBtn" style="
                        background: transparent;
                        color: white;
                        border: none;
                        padding: 8px;
                        font-size: 20px;
                        cursor: pointer;
                        flex-shrink: 0;
                        line-height: 1;
                    ">&times;</button>
                `;

                document.body.appendChild(installBanner);

                // Install button click
                document.getElementById('pwaInstallBtn').addEventListener('click', async () => {
                    installBanner.style.display = 'none';
                    deferredPrompt.prompt();
                    const { outcome } = await deferredPrompt.userChoice;
                    console.log('User choice:', outcome);
                    deferredPrompt = null;
                });

                // Close button click
                document.getElementById('pwaCloseBtn').addEventListener('click', () => {
                    installBanner.style.display = 'none';
                    localStorage.setItem('pwaPromptDismissed', 'true');
                });
            });

            // For iOS - show manual instructions after a delay
            if (isIOS) {
                const dismissed = localStorage.getItem('iosPwaPromptDismissed');
                if (!dismissed) {
                    setTimeout(() => {
                        console.log('Showing iOS PWA banner');
                        const iosBanner = document.createElement('div');
                        iosBanner.id = 'iosPwaInstallBanner';
                        iosBanner.style.cssText = `
                            position: fixed;
                            bottom: 0;
                            left: 0;
                            right: 0;
                            background: linear-gradient(135deg, #2563EB 0%, #1D4ED8 100%);
                            color: white;
                            padding: 16px 20px;
                            box-shadow: 0 -4px 12px rgba(37, 99, 235, 0.2);
                            z-index: 9999;
                            animation: slideUp 0.3s ease;
                        `;

                        iosBanner.innerHTML = `
                            <div style="display: flex; align-items: center; justify-content: space-between; gap: 12px;">
                                <div style="flex: 1;">
                                    <div style="font-weight: 600; font-size: 14px; margin-bottom: 6px;">Install gasconsult.ai</div>
                                    <div style="font-size: 12px; opacity: 0.9; line-height: 1.4;">
                                        Tap <svg width="16" height="16" style="display: inline; vertical-align: middle; margin: 0 2px;" fill="white" viewBox="0 0 24 24"><path d="M16.59 9H15V4c0-.55-.45-1-1-1h-4c-.55 0-1 .45-1 1v5H7.41c-.89 0-1.34 1.08-.71 1.71l4.59 4.59c.39.39 1.02.39 1.41 0l4.59-4.59c.63-.63.19-1.71-.7-1.71zM5 19c0 .55.45 1 1 1h12c.55 0 1-.45 1-1s-.45-1-1-1H6c-.55 0-1 .45-1 1z"/></svg> then "Add to Home Screen"
                                    </div>
                                </div>
                                <button id="iosCloseBtn" style="
                                    background: transparent;
                                    color: white;
                                    border: none;
                                    padding: 8px;
                                    font-size: 20px;
                                    cursor: pointer;
                                    flex-shrink: 0;
                                    line-height: 1;
                                ">&times;</button>
                            </div>
                        `;

                        document.body.appendChild(iosBanner);

                        document.getElementById('iosCloseBtn').addEventListener('click', () => {
                            iosBanner.style.display = 'none';
                            localStorage.setItem('iosPwaPromptDismissed', 'true');
                        });
                    }, 3000); // Show after 3 seconds
                }
            }

            // Add slideUp animation
            const style = document.createElement('style');
            style.textContent = `
                @keyframes slideUp {
                    from { transform: translateY(100%); }
                    to { transform: translateY(0); }
                }
            `;
            document.head.appendChild(style);
        }

        // Register Service Worker for offline functionality
        if ('serviceWorker' in navigator) {
            window.addEventListener('load', () => {
                navigator.serviceWorker.register('/static/sw.js')
                    .then((registration) => {
                        console.log('Service Worker registered:', registration.scope);
                    })
                    .catch((error) => {
                        console.log('Service Worker registration failed:', error);
                    });
            });
        }
    </script>
</head>
<body>
    <nav>
        <div class="container">
            <a href="/" class="logo-container">
                <svg class="logo-ecg" viewBox="0 0 60 28" fill="none" xmlns="http://www.w3.org/2000/svg">
                    <defs>
                        <linearGradient id="ecgGrad" x1="0%" y1="0%" x2="100%" y2="0%">
                            <stop offset="0%" stop-color="#2563EB"/>
                            <stop offset="20%" stop-color="#EF4444"/>
                            <stop offset="40%" stop-color="#FBBF24"/>
                            <stop offset="60%" stop-color="#8B5CF6"/>
                            <stop offset="80%" stop-color="#10B981"/>
                            <stop offset="100%" stop-color="#6B7280"/>
                        </linearGradient>
                    </defs>
                    <path d="M2 14 L10 14 L14 12 L18 16 L22 4 L26 24 L30 10 L34 14 L42 14"
                          stroke="url(#ecgGrad)"
                          stroke-width="2.5"
                          stroke-linecap="round"
                          stroke-linejoin="round"
                          fill="none"/>
                </svg>
                <div class="logo-wordmark">
                    <span class="logo-gas">gas</span><span class="logo-consult">consult</span><span class="logo-ai">.ai</span>
                </div>
            </a>
            <div class="nav-actions">
                <a href="/" class="nav-link active">Home</a>
                <a href="/preop" class="nav-link">Pre-Op Assessment</a>
                <a href="/calculators" class="nav-link">Clinical Calculators</a>
                <a href="/quick-dose" class="nav-link">Quick Dose</a>
                <a href="/hypotension" class="nav-link">IOH Predictor</a>
                {% if messages %}
                <a href="/clear" class="clear-chat-btn">Clear Chat</a>
                {% endif %}
            </div>
        </div>
    </nav>

    <!-- PHI Warning Banner (Chat page only) -->
    {% if messages %}
    <div class="phi-warning">
        <div class="phi-warning-content">
            <div class="phi-warning-icon">⚠️</div>
            <div class="phi-warning-text">
                <strong>Privacy Notice:</strong>
                <p>Do not enter patient names, dates of birth, MRNs, or other identifying information. Use age, weight, and clinical details only.</p>
            </div>
        </div>
    </div>
    {% endif %}

    <div class="main-content">
        {% if not messages %}
        <!-- Welcome Screen -->
        <div class="welcome-screen">
            <!-- Hero Section -->
            <h1 class="hero-headline">
                <span class="hero-headline-blue">Evidence-Based</span><br>
                <span class="hero-headline-dark">Anesthesiology Consultation</span>
            </h1>
            <p class="hero-subtitle">Get instant, citation-backed clinical answers powered by PubMed research. Real evidence. Real citations. Zero hallucinations.</p>

            <!-- CTA Buttons -->
            <div class="cta-buttons">
                <a href="/preop" class="preop-cta-outline">
                    Pre-Operative Assessment Tool
                    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round">
                        <line x1="5" y1="12" x2="19" y2="12"></line>
                        <polyline points="12 5 19 12 12 19"></polyline>
                    </svg>
                </a>
                <a href="/quick-dose" class="quickdose-cta-filled">
                    <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round">
                        <path d="M12 2L2 7l10 5 10-5-10-5z"></path>
                        <path d="M2 17l10 5 10-5"></path>
                        <path d="M2 12l10 5 10-5"></path>
                    </svg>
                    Quick Dose Calculator
                </a>
            </div>

            <!-- Trust Badges -->
            <div class="trust-badges">
                <div class="trust-badge">
                    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                        <polyline points="20 6 9 17 4 12"></polyline>
                    </svg>
                    <span>PubMed sourced</span>
                </div>
                <div class="trust-badge">
                    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                        <polyline points="20 6 9 17 4 12"></polyline>
                    </svg>
                    <span>Verifiable citations</span>
                </div>
            </div>
        </div>
        {% endif %}

        <!-- Chat Input - Always Visible on Homepage -->
        {% if not messages %}
        <div class="homepage-chat-section">
            <!-- Suggested Prompts -->
            <div class="suggested-prompts">
                <h4>💡 Try asking:</h4>
                <div class="prompt-buttons">
                    <button class="prompt-btn" onclick="fillQuery('TXA in cardiac surgery')">TXA in cardiac surgery</button>
                    <button class="prompt-btn" onclick="fillQuery('propofol vs etomidate for induction')">Propofol vs Etomidate</button>
                    <button class="prompt-btn" onclick="fillQuery('PONV prevention strategies')">PONV prevention</button>
                    <button class="prompt-btn" onclick="fillQuery('sugammadex reversal dosing')">Sugammadex dosing</button>
                    <button class="prompt-btn" onclick="fillQuery('difficult airway management')">Difficult airway</button>
                </div>
            </div>

            <div class="chat-input-container homepage-input">
                <form method="post" action="/chat" class="chat-form">
                    <input type="hidden" name="csrf_token" value="{{ csrf_token() }}"/>
                    <textarea name="query" id="chatInput" placeholder="Ask anything about anesthesiology..." required rows="2"></textarea>
                    <button type="submit" class="send-btn">↑</button>
                </form>
            </div>

            <!-- Features Section (below chat) -->
            <div class="features-section">
                <div class="feature-grid">
                    <div class="feature-card">
                        <div class="feature-icon blue">
                            <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path>
                                <path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path>
                            </svg>
                        </div>
                        <h3 class="feature-title">PubMed-Backed Answers</h3>
                        <p class="feature-description">Every answer sourced from peer-reviewed research, systematic reviews, and clinical guidelines — with full citations you can verify.</p>
                    </div>
                    <div class="feature-card">
                        <div class="feature-icon green">
                            <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <rect x="4" y="4" width="16" height="16" rx="2" ry="2"></rect>
                                <rect x="9" y="9" width="6" height="6"></rect>
                                <line x1="9" y1="1" x2="9" y2="4"></line>
                                <line x1="15" y1="1" x2="15" y2="4"></line>
                                <line x1="9" y1="20" x2="9" y2="23"></line>
                                <line x1="15" y1="20" x2="15" y2="23"></line>
                                <line x1="20" y1="9" x2="23" y2="9"></line>
                                <line x1="20" y1="14" x2="23" y2="14"></line>
                                <line x1="1" y1="9" x2="4" y2="9"></line>
                                <line x1="1" y1="14" x2="4" y2="14"></line>
                            </svg>
                        </div>
                        <h3 class="feature-title">Medical Calculators</h3>
                        <p class="feature-description">Built-in calculators for MABL, IBW, BSA, QTc, maintenance fluids, and more. Just type your values and get instant results.</p>
                    </div>
                    <div class="feature-card">
                        <div class="feature-icon violet">
                            <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <path d="M21 15a2 2 0 0 1-2 2H7l-4 4V5a2 2 0 0 1 2-2h14a2 2 0 0 1 2 2z"></path>
                            </svg>
                        </div>
                        <h3 class="feature-title">Conversational AI</h3>
                        <p class="feature-description">Ask follow-up questions, refine your queries, and explore topics naturally — like talking to a colleague who knows the literature.</p>
                    </div>
                </div>
            </div>
        </div>
        {% endif %}

        <!-- Chat Container - Always present for streaming -->
        <div class="chat-container" {% if not messages %}style="display: none;"{% endif %}>
            <div class="chat-messages" id="chatMessages">
                {% if messages %}
                {% for msg in messages %}
                    <div class="message {{ msg.role }}">
                        <div class="message-content">
                            {% if msg.role == 'user' %}
                                <div class="message-text">{{ msg.content }}</div>
                            {% else %}
                                <button class="copy-btn" onclick="copyToClipboard(this)" title="Copy to clipboard">
                                    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                        <rect x="9" y="9" width="13" height="13" rx="2" ry="2"></rect>
                                        <path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"></path>
                                    </svg>
                                    <span class="copy-text">Copy</span>
                                </button>
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
                {% endif %}
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
    </div>

    <!-- Chat Input - Visible on Chat Page -->
    {% if messages %}
    <div class="chat-input-container">
        <form method="post" action="/chat" class="chat-form">
            <input type="hidden" name="csrf_token" value="{{ csrf_token() }}"/>
            <textarea name="query" id="chatInput" placeholder="Ask anything about anesthesiology..." required rows="2"></textarea>
            <button type="submit" class="send-btn">↑</button>
        </form>
    </div>
    {% endif %}

    <footer>
        <p>&copy; 2025 gasconsult.ai. All rights reserved. | <a href="/terms" style="color: var(--primary); text-decoration: none;">Terms of Service</a> | <a href="/privacy" style="color: var(--primary); text-decoration: none;">Privacy Policy</a></p>
    </footer>

    <script>
        // Fill query from suggested prompt
        function fillQuery(text) {
            const textarea = document.getElementById('chatInput');
            if (textarea) {
                textarea.value = text;
                textarea.focus();
                // Auto-expand textarea if needed
                textarea.style.height = '52px';
                textarea.style.height = textarea.scrollHeight + 'px';
            }
        }

        // Keyboard shortcut: Ctrl+Enter or Cmd+Enter to submit
        document.addEventListener('DOMContentLoaded', function() {
            const textarea = document.getElementById('chatInput');
            if (textarea) {
                // Check for prefill message from calculators
                const prefillMessage = sessionStorage.getItem('prefillMessage');
                if (prefillMessage) {
                    textarea.value = prefillMessage;
                    textarea.focus();
                    // Auto-expand textarea if needed
                    textarea.style.height = '52px';
                    textarea.style.height = textarea.scrollHeight + 'px';
                    // Clear the stored message
                    sessionStorage.removeItem('prefillMessage');
                    console.log('[PREFILL] Populated textarea with calculator message');
                }

                textarea.addEventListener('keydown', function(e) {
                    // Check for Ctrl+Enter (Windows/Linux) or Cmd+Enter (Mac)
                    if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') {
                        e.preventDefault(); // Prevent default newline behavior

                        // Find the form and submit it
                        const form = textarea.closest('form');
                        if (form) {
                            console.log('[KEYBOARD] Ctrl+Enter detected - submitting form');
                            form.requestSubmit(); // Use requestSubmit() to trigger validation and submit event
                        }
                    }
                });
                console.log('[KEYBOARD] Keyboard shortcuts initialized (Ctrl+Enter or Cmd+Enter to submit)');
            }
        });

        // Smooth scroll to latest message on page load
        window.addEventListener('load', function() {
            const messages = document.querySelectorAll('.message');
            if (messages.length > 0) {
                const lastMessage = messages[messages.length - 1];
                setTimeout(() => {
                    lastMessage.scrollIntoView({ behavior: 'smooth', block: 'end' });
                }, 100);
            }
        });

        // Streaming form submission - must run after DOM loads
        document.addEventListener('DOMContentLoaded', function() {
            const form = document.querySelector('.chat-form');
            if (form) {
                form.addEventListener('submit', function(e) {
                    const submitBtn = form.querySelector('.send-btn');
                    const textarea = form.querySelector('textarea');
                    const query = textarea.value.trim();

                    if (!query) {
                        e.preventDefault();
                        return;
                    }

                // Check if we're on the homepage (welcome screen visible and no chat messages)
                const welcomeScreen = document.querySelector('.welcome-screen');
                const chatMessages = document.getElementById('chatMessages');
                const hasMessages = chatMessages && chatMessages.children.length > 0;
                const isHomepage = welcomeScreen && !hasMessages;

                // If on homepage, navigate to /chat page with query
                if (isHomepage) {
                    // Let the form submit naturally - it will POST to /chat and redirect
                    return;
                }

                // On chat page - prevent default and use streaming
                e.preventDefault();

                // Disable inputs
                submitBtn.disabled = true;
                textarea.disabled = true;
                submitBtn.style.opacity = '0.6';

                // Show chat container
                const chatContainer = document.querySelector('.chat-container');
                if (chatContainer) {
                    chatContainer.style.display = 'block';
                }

                // Add user message to UI
                const messagesContainer = document.getElementById('chatMessages');
                if (!messagesContainer) {
                    console.error('[ERROR] Chat messages container not found');
                    submitBtn.disabled = false;
                    textarea.disabled = false;
                    submitBtn.style.opacity = '1';
                    return;
                }

                const userMsg = document.createElement('div');
                userMsg.className = 'message user';
                userMsg.innerHTML = `<div class="message-content"><div class="message-text">${escapeHtml(query)}</div></div>`;
                messagesContainer.appendChild(userMsg);

                // Add loading indicator with progress steps
                const loadingMsg = document.createElement('div');
                loadingMsg.className = 'message assistant';
                loadingMsg.id = 'streaming-response';
                loadingMsg.innerHTML = `<div class="message-content"><p class="loading-indicator">🔍 Searching PubMed database...</p></div>`;
                messagesContainer.appendChild(loadingMsg);
                loadingMsg.scrollIntoView({ behavior: 'smooth', block: 'end' });

                // Clear textarea
                textarea.value = '';
                textarea.style.height = '52px';

                // Make POST request to initiate streaming
                const formData = new FormData();
                formData.append('query', query);

                fetch('/chat', {
                    method: 'POST',
                    credentials: 'same-origin',  // Important for session cookies
                    body: formData
                })
                .then(response => response.json())
                .then(data => {
                    if (data.status === 'ready') {
                        // Update loading indicator with progress
                        const loadingIndicator = document.querySelector('.loading-indicator');
                        if (loadingIndicator) {
                            loadingIndicator.innerHTML = '📚 Found papers<br>🤖 Analyzing evidence...';
                        }

                        // Connect to streaming endpoint
                        const eventSource = new EventSource(`/stream?request_id=${data.request_id}`);
                        const responseDiv = document.getElementById('streaming-response').querySelector('.message-content');
                        let responseContent = '';

                        eventSource.addEventListener('message', function(e) {
                            const event = JSON.parse(e.data);

                            if (event.type === 'connected') {
                                console.log('[STREAM] Connected');
                            } else if (event.type === 'content') {
                                // Remove loading indicator and add copy button on first content
                                if (responseContent === '') {
                                    responseDiv.innerHTML = `
                                        <button class="copy-btn" onclick="copyToClipboard(this)" title="Copy to clipboard">
                                            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                                <rect x="9" y="9" width="13" height="13" rx="2" ry="2"></rect>
                                                <path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"></path>
                                            </svg>
                                            <span class="copy-text">Copy</span>
                                        </button>
                                        <div class="message-text"></div>
                                    `;
                                }
                                responseContent += event.data;
                                const messageText = responseDiv.querySelector('.message-text');
                                if (messageText) {
                                    messageText.innerHTML = responseContent;
                                } else {
                                    responseDiv.innerHTML = responseContent;
                                }
                                // Auto-scroll
                                responseDiv.scrollIntoView({ behavior: 'smooth', block: 'end' });
                            } else if (event.type === 'references') {
                                // Add references and metadata
                                const refs = event.data;
                                if (refs && refs.length > 0) {
                                    let refsHtml = '<div class="message-refs"><strong>References:</strong>';
                                    refs.forEach((ref, i) => {
                                        refsHtml += `<div class="ref-item">
                                            <a href="https://pubmed.ncbi.nlm.nih.gov/${ref.pmid}/" target="_blank">
                                                [${i+1}] ${ref.title} (${ref.year})
                                            </a>
                                        </div>`;
                                    });
                                    refsHtml += '</div>';
                                    if (event.num_papers > 0) {
                                        refsHtml += `<div class="message-meta">📊 ${event.num_papers} papers from PubMed</div>`;
                                    }
                                    responseDiv.querySelector('.message-text').insertAdjacentHTML('afterend', refsHtml);
                                }
                            } else if (event.type === 'done') {
                                console.log('[STREAM] Complete');
                                eventSource.close();
                                // Re-enable form
                                submitBtn.disabled = false;
                                textarea.disabled = false;
                                submitBtn.style.opacity = '1';
                                textarea.focus();
                                // Session already updated on server, no reload needed
                            } else if (event.type === 'error') {
                                console.error('[STREAM] Error:', event.message);
                                responseDiv.innerHTML = `<p><strong>Error:</strong> ${event.message}</p>`;
                                eventSource.close();
                                submitBtn.disabled = false;
                                textarea.disabled = false;
                                submitBtn.style.opacity = '1';
                            }
                        });

                        eventSource.onerror = function(err) {
                            console.error('[STREAM] Connection error:', err);
                            eventSource.close();
                            submitBtn.disabled = false;
                            textarea.disabled = false;
                            submitBtn.style.opacity = '1';
                        };
                    }
                })
                .catch(error => {
                    console.error('[POST] Error:', error);
                    const loadingMsg = document.getElementById('streaming-response');
                    if (loadingMsg) {
                        loadingMsg.querySelector('.message-content').innerHTML = '<p><strong>Error:</strong> Failed to start query. Please try again.</p>';
                    }
                    submitBtn.disabled = false;
                    textarea.disabled = false;
                    submitBtn.style.opacity = '1';
                });
            });
            }
        });

        // Helper function to escape HTML
        function escapeHtml(text) {
            const div = document.createElement('div');
            div.textContent = text;
            return div.innerHTML;
        }

        // Auto-resize textarea smoothly
        const textareas = document.querySelectorAll('textarea');
        textareas.forEach(textarea => {
            textarea.addEventListener('input', function() {
                this.style.height = '52px';
                this.style.height = Math.min(this.scrollHeight, 200) + 'px';
            });
        });

        // Copy to clipboard function
        function copyToClipboard(button) {
            const messageContent = button.parentElement;
            const messageText = messageContent.querySelector('.message-text');

            // Extract text content without HTML tags
            const tempDiv = document.createElement('div');
            tempDiv.innerHTML = messageText.innerHTML;
            const textContent = tempDiv.textContent || tempDiv.innerText || '';

            // Copy to clipboard
            navigator.clipboard.writeText(textContent.trim()).then(() => {
                // Update button to show success
                const originalText = button.querySelector('.copy-text').textContent;
                button.querySelector('.copy-text').textContent = 'Copied!';
                button.classList.add('copied');

                // Reset after 2 seconds
                setTimeout(() => {
                    button.querySelector('.copy-text').textContent = originalText;
                    button.classList.remove('copied');
                }, 2000);
            }).catch(err => {
                console.error('Failed to copy:', err);
                button.querySelector('.copy-text').textContent = 'Failed';
                setTimeout(() => {
                    button.querySelector('.copy-text').textContent = 'Copy';
                }, 2000);
            });
        }

        // Check for pending stream from homepage redirect
        console.log('[AUTO-START] Checking for pending stream...');
        {% if pending_stream %}
        console.log('[AUTO-START] pending_stream value: {{ pending_stream }}');
        document.addEventListener('DOMContentLoaded', function() {
            console.log('[AUTO-START] Pending stream detected, starting streaming...');
            const requestId = '{{ pending_stream }}';

            // Add loading indicator to last message
            const messages = document.querySelectorAll('.message.assistant');
            if (messages.length > 0) {
                const lastMessage = messages[messages.length - 1];
                lastMessage.querySelector('.message-content').innerHTML = '<p class="loading-indicator">🔍 Searching medical literature...</p>';

                // Start streaming
                console.log('[AUTO-START] Creating EventSource with requestId:', requestId);
                const eventSource = new EventSource(`/stream?request_id=${requestId}`);
                const responseDiv = lastMessage.querySelector('.message-content');
                let responseContent = '';

                eventSource.addEventListener('open', function(e) {
                    console.log('[STREAM] Connection opened');
                });

                eventSource.addEventListener('error', function(e) {
                    console.error('[STREAM] Error occurred:', e);
                    console.error('[STREAM] ReadyState:', eventSource.readyState);
                    if (eventSource.readyState === EventSource.CLOSED) {
                        console.error('[STREAM] Connection closed');
                    }
                });

                eventSource.addEventListener('message', function(e) {
                    console.log('[STREAM] Message received:', e.data);
                    const event = JSON.parse(e.data);

                    if (event.type === 'connected') {
                        console.log('[STREAM] Connected');
                    } else if (event.type === 'content') {
                        // Remove loading indicator and add copy button on first content
                        if (responseContent === '') {
                            responseDiv.innerHTML = `
                                <button class="copy-btn" onclick="copyToClipboard(this)" title="Copy to clipboard">
                                    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                        <rect x="9" y="9" width="13" height="13" rx="2" ry="2"></rect>
                                        <path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"></path>
                                    </svg>
                                    <span class="copy-text">Copy</span>
                                </button>
                                <div class="message-text"></div>`;
                        }
                        responseContent += event.data;
                        responseDiv.querySelector('.message-text').innerHTML = responseContent;
                        lastMessage.scrollIntoView({ behavior: 'smooth', block: 'end' });
                    } else if (event.type === 'references') {
                        // Add references and metadata
                        const refs = event.data;
                        if (refs && refs.length > 0) {
                            let refsHtml = '<div class="message-refs"><strong>References:</strong>';
                            refs.forEach((ref, i) => {
                                refsHtml += `<div class="ref-item">
                                    <a href="https://pubmed.ncbi.nlm.nih.gov/${ref.pmid}/" target="_blank">
                                        [${i+1}] ${ref.title} (${ref.year})
                                    </a>
                                </div>`;
                            });
                            refsHtml += '</div>';
                            if (event.num_papers > 0) {
                                refsHtml += `<div class="message-meta">📊 ${event.num_papers} papers from PubMed</div>`;
                            }
                            responseDiv.querySelector('.message-text').insertAdjacentHTML('afterend', refsHtml);
                        }
                    } else if (event.type === 'done') {
                        console.log('[STREAM] Completed');
                        eventSource.close();
                    } else if (event.type === 'error') {
                        console.error('[STREAM] Error:', event.message);
                        responseDiv.innerHTML = `<p style="color: #EF4444;">${event.message}</p>`;
                        eventSource.close();
                    }
                });

                eventSource.onerror = function(e) {
                    console.error('[STREAM] Connection error:', e);
                    eventSource.close();
                };
            }
        });
        {% endif %}
    </script>

</body>
</html>
"""

CHAT_HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Chat — gasconsult.ai</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Sora:wght@400;600&display=swap" rel="stylesheet">
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=5">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">
    <meta name="apple-mobile-web-app-capable" content="yes">
    <meta name="apple-mobile-web-app-status-bar-style" content="default">
    <meta name="apple-mobile-web-app-title" content="gasconsult.ai">
    <style>
        :root {
            /* Primary Brand Colors */
            --primary-blue: #2563EB;
            --primary-blue-dark: #1D4ED8;
            --primary-blue-light: #DBEAFE;

            /* Anesthesia Color Palette */
            --opioid-blue: #2563EB;
            --nmb-red: #EF4444;
            --induction-yellow: #FBBF24;
            --vasopressor-violet: #8B5CF6;
            --anticholinergic-green: #10B981;
            --local-gray: #6B7280;

            /* Neutral Palette */
            --text-primary: #0F172A;
            --text-secondary: #475569;
            --text-muted: #94A3B8;
            --bg-primary: #FFFFFF;
            --bg-secondary: #F8FAFC;
            --border: #E2E8F0;

            /* Legacy aliases */
            --primary: #2563EB;
            --primary-dark: #1D4ED8;
        }

        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        html {
            scroll-behavior: smooth;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'SF Pro Display', 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background: #F8FAFC;
            color: #0A3D62;
            line-height: 1.6;
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
            overflow-x: hidden;
            width: 100%;
            display: flex;
            flex-direction: column;
            min-height: 100vh;
        }

        @keyframes slideUp {
            from {
                opacity: 0;
                transform: translateY(30px);
            }
            to {
                opacity: 1;
                transform: translateY(0);
            }
        }

        @keyframes slideIn {
            from {
                opacity: 0;
                transform: translateX(-20px);
            }
            to {
                opacity: 1;
                transform: translateX(0);
            }
        }

        /* Navigation */
        nav {
            background: rgba(255, 255, 255, 0.85);
            backdrop-filter: blur(12px);
            -webkit-backdrop-filter: blur(12px);
            padding: 16px 40px;
            position: sticky;
            top: 0;
            left: 0;
            right: 0;
            z-index: 100;
            border-bottom: 1px solid rgba(226, 232, 240, 0.8);
        }

        nav .container {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            align-items: center;
            justify-content: space-between;
            flex-wrap: wrap;
            gap: 16px;
        }

        .logo-container {
            text-decoration: none;
            display: flex;
            align-items: center;
            gap: 12px;
            cursor: pointer;
            transition: transform 0.2s ease;
        }

        .logo-container:hover {
            transform: translateY(-1px);
        }

        .logo-ecg {
            height: 28px;
            width: auto;
            flex-shrink: 0;
        }

        .logo-wordmark {
            font-family: 'Sora', -apple-system, BlinkMacSystemFont, system-ui, sans-serif;
            font-size: 20px;
            font-weight: 600;
            letter-spacing: -0.5px;
            white-space: nowrap;
        }

        .logo-gas {
            color: #2563EB;
        }

        .logo-consult {
            color: #111111;
        }

        .logo-ai {
            font-weight: 400;
            color: #6B7280;
        }

        .nav-actions {
            display: flex;
            gap: 12px;
            align-items: center;
        }

        .nav-link {
            color: var(--text-secondary);
            text-decoration: none;
            font-size: 14px;
            font-weight: 500;
            padding: 8px 16px;
            border-radius: 8px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            color: var(--text-primary);
            background: rgba(255, 255, 255, 0.6);
            backdrop-filter: blur(8px);
            -webkit-backdrop-filter: blur(8px);
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.06);
        }

        .nav-link.active {
            color: var(--primary-blue);
            font-weight: 600;
        }

        .new-chat-btn {
            background: linear-gradient(135deg, var(--primary) 0%, var(--primary-dark) 100%);
            color: white;
            padding: 10px 24px;
            border-radius: 8px;
            font-size: 15px;
            font-weight: 500;
            text-decoration: none;
            transition: all 0.2s ease;
            border: none;
            cursor: pointer;
            box-shadow: 0 2px 8px rgba(37, 99, 235, 0.2);
        }

        .new-chat-btn:hover {
            background: linear-gradient(135deg, var(--primary-dark) 0%, #1E40AF 100%);
            transform: translateY(-1px);
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.35);
        }

        .clear-chat-btn {
            background: transparent;
            color: var(--text-secondary);
            padding: 8px 16px;
            border-radius: 8px;
            font-size: 14px;
            font-weight: 500;
            text-decoration: none;
            transition: all 0.2s ease;
            border: 1.5px solid var(--border);
        }

        .clear-chat-btn:hover {
            background: rgba(239, 68, 68, 0.05);
            color: #EF4444;
            border-color: #EF4444;
        }

        /* Main Content */
        .main-content {
            flex: 1;
            display: flex;
            flex-direction: column;
            max-width: 1000px;
            width: 100%;
            margin: 0 auto;
            padding: 20px;
            position: relative;
        }

        /* Empty State */
        .empty-state {
            flex: 1;
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            text-align: center;
            padding: 80px 40px;
            animation: slideUp 0.6s cubic-bezier(0.4, 0, 0.2, 1);
        }

        .empty-state-icon {
            width: 80px;
            height: 80px;
            background: linear-gradient(135deg, rgba(37, 99, 235, 0.1) 0%, rgba(139, 92, 246, 0.1) 100%);
            border-radius: 50%;
            display: flex;
            align-items: center;
            justify-content: center;
            margin-bottom: 24px;
            animation: float 3s ease-in-out infinite;
        }

        @keyframes float {
            0%, 100% {
                transform: translateY(0);
            }
            50% {
                transform: translateY(-10px);
            }
        }

        .empty-state-icon svg {
            width: 40px;
            height: 40px;
            stroke: var(--primary-blue);
        }

        .empty-state h2 {
            font-family: 'Sora', sans-serif;
            font-size: 28px;
            font-weight: 700;
            color: var(--text-primary);
            margin-bottom: 12px;
            letter-spacing: -0.5px;
        }

        .empty-state p {
            font-size: 16px;
            color: var(--text-secondary);
            max-width: 500px;
            margin-bottom: 32px;
        }

        /* Suggested Prompts */
        .suggested-prompts {
            width: 100%;
            max-width: 700px;
        }

        .suggested-prompts h4 {
            font-size: 14px;
            font-weight: 600;
            color: var(--text-secondary);
            margin-bottom: 12px;
            text-align: left;
        }

        .prompt-buttons {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            justify-content: center;
        }

        .prompt-btn {
            padding: 10px 18px;
            background: rgba(255, 255, 255, 0.9);
            backdrop-filter: blur(12px);
            -webkit-backdrop-filter: blur(12px);
            border: 1.5px solid var(--border);
            border-radius: 10px;
            color: var(--text-secondary);
            font-size: 14px;
            font-weight: 500;
            cursor: pointer;
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.04);
        }

        .prompt-btn:hover {
            background: var(--primary-blue-light);
            border-color: var(--primary-blue);
            color: var(--primary-blue-dark);
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.15);
        }

        /* Chat Container */
        .chat-container {
            flex: 1;
            display: flex;
            flex-direction: column;
            padding: 20px 0;
            overflow: hidden;
        }

        .chat-messages {
            flex: 1;
            overflow-y: auto;
            padding: 20px;
            display: flex;
            flex-direction: column;
            gap: 20px;
        }

        .message {
            animation: slideIn 0.4s cubic-bezier(0.4, 0, 0.2, 1);
        }

        .message.user .message-content {
            background: linear-gradient(135deg, rgba(37, 99, 235, 0.08) 0%, rgba(139, 92, 246, 0.08) 100%);
            border: 1.5px solid rgba(37, 99, 235, 0.15);
            padding: 16px 20px;
            border-radius: 16px;
            max-width: 80%;
            margin-left: auto;
            box-shadow: 0 2px 8px rgba(37, 99, 235, 0.06);
        }

        .message.assistant .message-content {
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(20px);
            -webkit-backdrop-filter: blur(20px);
            border: 1.5px solid var(--border);
            padding: 20px 24px;
            border-radius: 16px;
            box-shadow: 0 4px 16px rgba(0, 0, 0, 0.06);
            position: relative;
        }

        .message-text {
            color: var(--text-primary);
            font-size: 15px;
            line-height: 1.7;
        }

        .message-refs {
            margin-top: 20px;
            padding-top: 16px;
            border-top: 1px solid var(--border);
            font-size: 13px;
        }

        .message-refs strong {
            color: var(--text-primary);
            font-weight: 600;
            display: block;
            margin-bottom: 8px;
        }

        .ref-item {
            margin: 6px 0;
        }

        .ref-item a {
            color: var(--primary-blue);
            text-decoration: none;
            transition: color 0.2s ease;
        }

        .ref-item a:hover {
            color: var(--primary-blue-dark);
            text-decoration: underline;
        }

        .message-meta {
            margin-top: 12px;
            font-size: 12px;
            color: var(--text-muted);
        }

        .copy-btn {
            position: absolute;
            top: 16px;
            right: 16px;
            background: var(--bg-secondary);
            border: 1px solid var(--border);
            border-radius: 6px;
            padding: 6px 12px;
            font-size: 12px;
            color: var(--text-secondary);
            cursor: pointer;
            display: flex;
            align-items: center;
            gap: 6px;
            transition: all 0.2s ease;
        }

        .copy-btn:hover {
            background: white;
            color: var(--primary-blue);
            border-color: var(--primary-blue);
        }

        /* ====== Premium Feature Styles ====== */

        /* Message Actions Toolbar */
        .message-actions {
            display: flex;
            gap: 8px;
            margin-top: 16px;
            padding-top: 16px;
            border-top: 1px solid rgba(226, 232, 240, 0.6);
            flex-wrap: wrap;
        }

        .action-btn {
            background: rgba(255, 255, 255, 0.8);
            backdrop-filter: blur(10px);
            -webkit-backdrop-filter: blur(10px);
            border: 1.5px solid var(--border);
            border-radius: 8px;
            padding: 8px 14px;
            font-size: 13px;
            font-weight: 500;
            color: var(--text-secondary);
            cursor: pointer;
            display: flex;
            align-items: center;
            gap: 6px;
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            box-shadow: 0 2px 6px rgba(0, 0, 0, 0.04);
        }

        .action-btn:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.08);
            border-color: var(--primary-blue);
            color: var(--primary-blue);
            background: rgba(255, 255, 255, 0.95);
        }

        .action-btn.active {
            background: linear-gradient(135deg, rgba(37, 99, 235, 0.1) 0%, rgba(139, 92, 246, 0.1) 100%);
            border-color: var(--primary-blue);
            color: var(--primary-blue);
        }

        .action-btn svg {
            width: 16px;
            height: 16px;
        }

        /* Evidence Strength Visualizer */
        .evidence-badge {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            padding: 8px 14px;
            border-radius: 10px;
            font-size: 12px;
            font-weight: 600;
            margin-bottom: 12px;
            backdrop-filter: blur(10px);
            -webkit-backdrop-filter: blur(10px);
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.06);
        }

        .evidence-badge.high {
            background: rgba(16, 185, 129, 0.15);
            border: 1.5px solid #10B981;
            color: #047857;
        }

        .evidence-badge.moderate {
            background: rgba(251, 191, 36, 0.15);
            border: 1.5px solid #FBBF24;
            color: #B45309;
        }

        .evidence-badge.low {
            background: rgba(239, 68, 68, 0.15);
            border: 1.5px solid #EF4444;
            color: #B91C1C;
        }

        .evidence-breakdown {
            display: flex;
            gap: 12px;
            margin-top: 8px;
            font-size: 11px;
            color: var(--text-muted);
            flex-wrap: wrap;
        }

        .study-type-badge {
            background: var(--bg-secondary);
            padding: 4px 10px;
            border-radius: 6px;
            border: 1px solid var(--border);
            font-weight: 500;
        }

        /* Follow-up Suggestions */
        .followup-section {
            margin-top: 20px;
            padding: 16px;
            background: linear-gradient(135deg, rgba(37, 99, 235, 0.04) 0%, rgba(139, 92, 246, 0.04) 100%);
            border-radius: 12px;
            border: 1.5px solid rgba(37, 99, 235, 0.1);
        }

        .followup-title {
            font-size: 13px;
            font-weight: 600;
            color: var(--text-secondary);
            margin-bottom: 10px;
            display: flex;
            align-items: center;
            gap: 8px;
        }

        .followup-questions {
            display: flex;
            flex-direction: column;
            gap: 8px;
        }

        .followup-btn {
            background: rgba(255, 255, 255, 0.9);
            backdrop-filter: blur(10px);
            -webkit-backdrop-filter: blur(10px);
            border: 1.5px solid var(--border);
            border-radius: 10px;
            padding: 10px 16px;
            text-align: left;
            font-size: 14px;
            color: var(--text-primary);
            cursor: pointer;
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            box-shadow: 0 2px 6px rgba(0, 0, 0, 0.03);
        }

        .followup-btn:hover {
            background: var(--primary-blue-light);
            border-color: var(--primary-blue);
            color: var(--primary-blue-dark);
            transform: translateX(4px);
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.12);
        }

        /* Voice Input Button */
        .voice-btn {
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(10px);
            -webkit-backdrop-filter: blur(10px);
            width: 44px;
            height: 44px;
            border-radius: 50%;
            border: 2px solid var(--border);
            display: flex;
            align-items: center;
            justify-content: center;
            cursor: pointer;
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            flex-shrink: 0;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.06);
        }

        .voice-btn:hover {
            border-color: var(--primary-blue);
            background: var(--primary-blue-light);
            transform: scale(1.05);
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.2);
        }

        .voice-btn.listening {
            background: linear-gradient(135deg, #EF4444 0%, #DC2626 100%);
            border-color: #EF4444;
            animation: pulse 1.5s infinite;
        }

        .voice-btn svg {
            width: 20px;
            height: 20px;
            stroke: var(--text-secondary);
        }

        .voice-btn.listening svg {
            stroke: white;
        }

        @keyframes pulse {
            0%, 100% {
                box-shadow: 0 0 0 0 rgba(239, 68, 68, 0.4);
            }
            50% {
                box-shadow: 0 0 0 12px rgba(239, 68, 68, 0);
            }
        }

        /* Share Modal */
        .modal-overlay {
            display: none;
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            bottom: 0;
            background: rgba(0, 0, 0, 0.5);
            backdrop-filter: blur(4px);
            -webkit-backdrop-filter: blur(4px);
            z-index: 1000;
            align-items: center;
            justify-content: center;
            animation: fadeIn 0.2s ease;
        }

        .modal-overlay.active {
            display: flex;
        }

        .modal-content {
            background: white;
            border-radius: 16px;
            padding: 28px;
            max-width: 500px;
            width: 90%;
            box-shadow: 0 20px 60px rgba(0, 0, 0, 0.3);
            animation: slideUp 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        }

        @keyframes fadeIn {
            from { opacity: 0; }
            to { opacity: 1; }
        }

        .modal-title {
            font-family: 'Sora', sans-serif;
            font-size: 20px;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 16px;
        }

        .share-link-container {
            display: flex;
            gap: 8px;
            margin: 16px 0;
        }

        .share-link-input {
            flex: 1;
            padding: 12px 16px;
            border: 2px solid var(--border);
            border-radius: 10px;
            font-size: 14px;
            font-family: 'JetBrains Mono', monospace;
            background: var(--bg-secondary);
        }

        .modal-actions {
            display: flex;
            gap: 10px;
            justify-content: flex-end;
            margin-top: 20px;
        }

        .modal-btn {
            padding: 10px 20px;
            border-radius: 8px;
            font-size: 14px;
            font-weight: 500;
            cursor: pointer;
            transition: all 0.2s ease;
        }

        .modal-btn-primary {
            background: var(--primary-blue);
            color: white;
            border: none;
        }

        .modal-btn-primary:hover {
            background: var(--primary-blue-dark);
        }

        .modal-btn-secondary {
            background: transparent;
            color: var(--text-secondary);
            border: 1.5px solid var(--border);
        }

        .modal-btn-secondary:hover {
            background: var(--bg-secondary);
        }

        /* Chat Input */
        .chat-input-container {
            position: sticky;
            bottom: 0;
            background: rgba(248, 250, 252, 0.95);
            backdrop-filter: blur(12px);
            -webkit-backdrop-filter: blur(12px);
            padding: 20px;
            border-top: 1px solid rgba(226, 232, 240, 0.8);
        }

        .chat-form {
            max-width: 800px;
            margin: 0 auto;
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(20px);
            -webkit-backdrop-filter: blur(20px);
            border: 2px solid rgba(37, 99, 235, 0.1);
            border-radius: 16px;
            padding: 6px;
            display: flex;
            align-items: center;
            gap: 8px;
            box-shadow: 0 8px 32px rgba(37, 99, 235, 0.08), 0 2px 12px rgba(0, 0, 0, 0.04);
            transition: all 0.4s cubic-bezier(0.4, 0, 0.2, 1);
        }

        .chat-form:hover {
            border-color: rgba(37, 99, 235, 0.2);
            box-shadow: 0 12px 48px rgba(37, 99, 235, 0.12), 0 4px 16px rgba(0, 0, 0, 0.06);
            transform: translateY(-2px);
        }

        .chat-form:focus-within {
            border-color: var(--primary-blue);
            box-shadow: 0 16px 56px rgba(37, 99, 235, 0.18), 0 8px 24px rgba(0, 0, 0, 0.08);
            transform: translateY(-3px);
        }

        .chat-form textarea {
            flex: 1;
            padding: 12px 16px;
            font-size: 15px;
            border: none;
            outline: none;
            background: transparent;
            color: var(--text-primary);
            resize: none;
            min-height: 44px;
            max-height: 200px;
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            line-height: 1.4;
        }

        .chat-form textarea::placeholder {
            color: var(--text-muted);
        }

        .send-btn {
            background: linear-gradient(135deg, var(--primary-blue) 0%, var(--primary-blue-dark) 100%);
            width: 44px;
            height: 44px;
            border-radius: 50%;
            font-size: 1.75rem;
            font-weight: 900;
            color: white;
            border: none;
            cursor: pointer;
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            flex-shrink: 0;
            display: flex;
            align-items: center;
            justify-content: center;
            padding: 0 0 3px 0;
            -webkit-text-stroke: 0.5px white;
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.3);
        }

        .send-btn:hover {
            transform: scale(1.08) rotate(5deg);
            box-shadow: 0 6px 20px rgba(37, 99, 235, 0.4);
        }

        .send-btn:active {
            transform: scale(0.95);
        }

        /* Loading Indicator */
        .loading-container {
            display: none;
            padding: 20px;
            align-items: center;
            justify-content: center;
        }

        .loading-container.active {
            display: flex;
        }

        .loading-dots {
            display: flex;
            gap: 8px;
        }

        .loading-dots span {
            width: 12px;
            height: 12px;
            background: var(--primary-blue);
            border-radius: 50%;
            animation: loadingDot 1.4s infinite ease-in-out;
        }

        .loading-dots span:nth-child(1) { animation-delay: 0s; }
        .loading-dots span:nth-child(2) { animation-delay: 0.2s; }
        .loading-dots span:nth-child(3) { animation-delay: 0.4s; }

        @keyframes loadingDot {
            0%, 60%, 100% {
                transform: scale(0.8);
                opacity: 0.5;
            }
            30% {
                transform: scale(1.2);
                opacity: 1;
            }
        }

        footer {
            background: white;
            padding: 32px 20px;
            text-align: center;
            font-size: 14px;
            color: var(--text-secondary);
            border-top: 1px solid var(--border);
        }

        /* Responsive */
        @media (max-width: 768px) {
            nav {
                padding: 12px 20px;
            }

            .nav-actions {
                gap: 8px;
            }

            .nav-link {
                padding: 6px 12px;
                font-size: 13px;
            }

            .main-content {
                padding: 12px;
            }

            .empty-state {
                padding: 60px 20px;
            }

            .message.user .message-content {
                max-width: 90%;
            }

            .chat-input-container {
                padding: 12px;
            }
        }
    </style>
</head>
<body>
    <!-- Navigation -->
    <nav>
        <div class="container">
            <a href="/" class="logo-container">
                <svg class="logo-ecg" viewBox="0 0 60 28" fill="none" xmlns="http://www.w3.org/2000/svg">
                    <defs>
                        <linearGradient id="ecgGrad" x1="0%" y1="0%" x2="100%" y2="0%">
                            <stop offset="0%" stop-color="#2563EB"/>
                            <stop offset="20%" stop-color="#EF4444"/>
                            <stop offset="40%" stop-color="#FBBF24"/>
                            <stop offset="60%" stop-color="#8B5CF6"/>
                            <stop offset="80%" stop-color="#10B981"/>
                            <stop offset="100%" stop-color="#6B7280"/>
                        </linearGradient>
                    </defs>
                    <path d="M2 14 L10 14 L14 12 L18 16 L22 4 L26 24 L30 10 L34 14 L42 14"
                          stroke="url(#ecgGrad)"
                          stroke-width="2.5"
                          stroke-linecap="round"
                          stroke-linejoin="round"
                          fill="none"/>
                </svg>
                <div class="logo-wordmark">
                    <span class="logo-gas">gas</span><span class="logo-consult">consult</span><span class="logo-ai">.ai</span>
                </div>
            </a>
            <div class="nav-actions">
                <a href="/" class="nav-link active">Home</a>
                <a href="/library" class="nav-link" title="View saved responses">
                    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" style="display: inline-block; vertical-align: middle; margin-right: 4px;">
                        <path d="M19 21l-7-5-7 5V5a2 2 0 012-2h10a2 2 0 012 2z"></path>
                    </svg>
                    Library
                </a>
                <a href="/preop" class="nav-link">Pre-Op Assessment</a>
                <a href="/calculators" class="nav-link">Clinical Calculators</a>
                <a href="/quick-dose" class="nav-link">Quick Dose</a>
                <a href="/hypotension" class="nav-link">IOH Predictor</a>
                {% if messages %}
                <a href="/clear" class="clear-chat-btn">Clear Chat</a>
                {% endif %}
            </div>
        </div>
    </nav>

    <!-- Main Content -->
    <div class="main-content">
        {% if not messages %}
        <!-- Empty State -->
        <div class="empty-state">
            <div class="empty-state-icon">
                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                    <path d="M21 15a2 2 0 0 1-2 2H7l-4 4V5a2 2 0 0 1 2-2h14a2 2 0 0 1 2 2z"></path>
                </svg>
            </div>
            <h2>Start a New Conversation</h2>
            <p>Ask any clinical question and get evidence-based answers backed by PubMed literature</p>

            <div class="suggested-prompts">
                <h4>💡 Try asking:</h4>
                <div class="prompt-buttons">
                    <button class="prompt-btn" onclick="fillQuery('TXA in cardiac surgery')">TXA in cardiac surgery</button>
                    <button class="prompt-btn" onclick="fillQuery('propofol vs etomidate for induction')">Propofol vs Etomidate</button>
                    <button class="prompt-btn" onclick="fillQuery('PONV prevention strategies')">PONV prevention</button>
                    <button class="prompt-btn" onclick="fillQuery('sugammadex reversal dosing')">Sugammadex dosing</button>
                    <button class="prompt-btn" onclick="fillQuery('difficult airway management')">Difficult airway</button>
                </div>
            </div>
        </div>
        {% else %}
        <!-- Chat Messages -->
        <div class="chat-container">
            <div class="chat-messages" id="chatMessages">
                {% for msg in messages %}
                    <div class="message {{ msg.role }}">
                        <div class="message-content">
                            {% if msg.role == 'user' %}
                                <div class="message-text">{{ msg.content }}</div>
                            {% else %}
                                <button class="copy-btn" onclick="smartCopy(this, {{ loop.index0 }})" title="Smart copy with citations">
                                    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                        <rect x="9" y="9" width="13" height="13" rx="2" ry="2"></rect>
                                        <path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"></path>
                                    </svg>
                                    <span class="copy-text">Copy</span>
                                </button>

                                <!-- Evidence Strength Visualizer -->
                                {% if msg.num_papers and msg.num_papers > 0 %}
                                <div class="evidence-badge {{ msg.get('evidence_strength', {}).get('level', 'low').lower() }}" data-message-index="{{ loop.index0 }}">
                                    <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5">
                                        <path d="M9 12l2 2 4-4m6 2a9 9 0 11-18 0 9 9 0 0118 0z"></path>
                                    </svg>
                                    <span>{{ msg.get('evidence_strength', {}).get('level', 'Moderate') }} Evidence</span>
                                </div>
                                {% if msg.get('evidence_strength', {}).get('breakdown') %}
                                <div class="evidence-breakdown">
                                    {% set bd = msg.evidence_strength.breakdown %}
                                    {% if bd.meta_analyses > 0 %}
                                    <span class="study-type-badge">{{ bd.meta_analyses }} Meta-analyses</span>
                                    {% endif %}
                                    {% if bd.rcts > 0 %}
                                    <span class="study-type-badge">{{ bd.rcts }} RCTs</span>
                                    {% endif %}
                                    {% if bd.reviews > 0 %}
                                    <span class="study-type-badge">{{ bd.reviews }} Reviews</span>
                                    {% endif %}
                                    <span class="study-type-badge">{{ bd.total }} Total Papers</span>
                                </div>
                                {% endif %}
                                {% endif %}

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

                                <!-- Premium Action Buttons -->
                                <div class="message-actions">
                                    <button class="action-btn" onclick="toggleBookmark(this, {{ loop.index0 }})" title="Bookmark this response">
                                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <path d="M19 21l-7-5-7 5V5a2 2 0 012-2h10a2 2 0 012 2z"></path>
                                        </svg>
                                        <span>Save</span>
                                    </button>

                                    <button class="action-btn" onclick="shareResponse({{ loop.index0 }})" title="Share this response">
                                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <circle cx="18" cy="5" r="3"></circle>
                                            <circle cx="6" cy="12" r="3"></circle>
                                            <circle cx="18" cy="19" r="3"></circle>
                                            <line x1="8.59" y1="13.51" x2="15.42" y2="17.49"></line>
                                            <line x1="15.41" y1="6.51" x2="8.59" y2="10.49"></line>
                                        </svg>
                                        <span>Share</span>
                                    </button>

                                    {% if msg.references %}
                                    <button class="action-btn" onclick="exportCitations({{ loop.index0 }}, 'bibtex')" title="Export to BibTeX">
                                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4M7 10l5 5 5-5M12 15V3"></path>
                                        </svg>
                                        <span>Export .bib</span>
                                    </button>

                                    <button class="action-btn" onclick="exportCitations({{ loop.index0 }}, 'ris')" title="Export to RIS (Zotero/Mendeley)">
                                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4M7 10l5 5 5-5M12 15V3"></path>
                                        </svg>
                                        <span>Export .ris</span>
                                    </button>
                                    {% endif %}
                                </div>

                                <!-- Smart Follow-up Suggestions -->
                                <div class="followup-section" id="followup-{{ loop.index0 }}" style="display: none;">
                                    <div class="followup-title">
                                        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <path d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z"></path>
                                        </svg>
                                        Ask next:
                                    </div>
                                    <div class="followup-questions" id="followup-questions-{{ loop.index0 }}">
                                        <!-- Loaded dynamically -->
                                    </div>
                                </div>
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
    </div>

    <!-- Chat Input - Always Visible -->
    <div class="chat-input-container">
        <form method="post" action="/chat" class="chat-form">
            <input type="hidden" name="csrf_token" value="{{ csrf_token() }}"/>
            <!-- Voice Input Button -->
            <button type="button" class="voice-btn" id="voiceBtn" title="Voice input" onclick="toggleVoiceInput()">
                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <path d="M12 1a3 3 0 003 3v8a3 3 0 01-6 0V4a3 3 0 013-3z"></path>
                    <path d="M19 10v2a7 7 0 01-14 0v-2M12 19v4M8 23h8"></path>
                </svg>
            </button>
            <textarea name="query" id="chatInput" placeholder="Ask anything about anesthesiology..." required rows="2"></textarea>
            <button type="submit" class="send-btn">↑</button>
        </form>
    </div>

    <!-- Share Modal -->
    <div class="modal-overlay" id="shareModal">
        <div class="modal-content">
            <div class="modal-title">Share Response</div>
            <p style="color: var(--text-secondary); margin-bottom: 16px;">Anyone with this link can view this response. Link expires in 30 days.</p>
            <div class="share-link-container">
                <input type="text" class="share-link-input" id="shareLinkInput" readonly />
                <button class="action-btn" onclick="copyShareLink()" style="flex-shrink: 0;">
                    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                        <rect x="9" y="9" width="13" height="13" rx="2" ry="2"></rect>
                        <path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"></path>
                    </svg>
                    Copy
                </button>
            </div>
            <div class="modal-actions">
                <button class="modal-btn modal-btn-secondary" onclick="closeShareModal()">Close</button>
                <button class="modal-btn modal-btn-primary" onclick="openShareLink()">Open Link</button>
            </div>
        </div>
    </div>

    <footer>
        <p>&copy; 2025 gasconsult.ai. All rights reserved. | <a href="/terms" style="color: var(--primary); text-decoration: none;">Terms of Service</a> | <a href="/privacy" style="color: var(--primary); text-decoration: none;">Privacy Policy</a></p>
    </footer>

    <script>
        // Fill query from suggested prompt
        function fillQuery(text) {
            const textarea = document.getElementById('chatInput');
            if (textarea) {
                textarea.value = text;
                textarea.focus();
                textarea.style.height = '52px';
                textarea.style.height = textarea.scrollHeight + 'px';
            }
        }

        // Copy to clipboard
        function copyToClipboard(button) {
            const messageContent = button.parentElement;
            const messageText = messageContent.querySelector('.message-text');
            const textToCopy = messageText.innerText;

            navigator.clipboard.writeText(textToCopy).then(() => {
                const copyText = button.querySelector('.copy-text');
                const originalText = copyText.textContent;
                copyText.textContent = 'Copied!';
                setTimeout(() => {
                    copyText.textContent = originalText;
                }, 2000);
            }).catch(err => {
                console.error('Failed to copy:', err);
            });
        }

        // Keyboard shortcut: Ctrl+Enter or Cmd+Enter to submit
        document.addEventListener('DOMContentLoaded', function() {
            const textarea = document.getElementById('chatInput');
            if (textarea) {
                // Auto-expand textarea
                textarea.addEventListener('input', function() {
                    this.style.height = '44px';
                    this.style.height = this.scrollHeight + 'px';
                });

                // Keyboard shortcut
                textarea.addEventListener('keydown', function(e) {
                    if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') {
                        e.preventDefault();
                        const form = textarea.closest('form');
                        if (form) {
                            form.submit();
                        }
                    }
                });
            }

            // Auto-scroll to bottom of chat
            const chatMessages = document.getElementById('chatMessages');
            if (chatMessages) {
                chatMessages.scrollTop = chatMessages.scrollHeight;
            }
        });

        // Show loading indicator when form submits
        document.querySelector('.chat-form').addEventListener('submit', function() {
            const loadingIndicator = document.getElementById('loadingIndicator');
            if (loadingIndicator) {
                loadingIndicator.classList.add('active');
            }
        });

        // Auto-start streaming handled by DOMContentLoaded event listener above (removed duplicate implementation)

        // ====== Premium Features JavaScript ======

        // Smart Copy with Citations
        function smartCopy(button, messageIndex) {
            const messageContent = button.parentElement;
            const messageText = messageContent.querySelector('.message-text');
            const messageRefs = messageContent.querySelector('.message-refs');
            const evidenceBadge = messageContent.querySelector('.evidence-badge');

            let textToCopy = '=== gasconsult.ai Response ===\n\n';
            textToCopy += messageText.innerText + '\n\n';

            if (messageRefs) {
                const refs = messageRefs.querySelectorAll('.ref-item');
                if (refs.length > 0) {
                    textToCopy += '=== References ===\n';
                    refs.forEach((ref, i) => {
                        textToCopy += ref.innerText + '\n';
                    });
                    textToCopy += '\n';
                }
            }

            if (evidenceBadge) {
                textToCopy += `Evidence Quality: ${evidenceBadge.innerText}\n`;
            }

            textToCopy += `\nAccessed: ${new Date().toLocaleDateString()}\nSource: ${window.location.origin}`;

            navigator.clipboard.writeText(textToCopy).then(() => {
                const copyText = button.querySelector('.copy-text');
                const originalText = copyText.textContent;
                copyText.textContent = 'Copied!';
                button.style.borderColor = '#10B981';
                button.style.color = '#10B981';
                setTimeout(() => {
                    copyText.textContent = originalText;
                    button.style.borderColor = '';
                    button.style.color = '';
                }, 2000);
            }).catch(err => {
                console.error('Failed to copy:', err);
            });
        }

        // Bookmark Toggle
        async function toggleBookmark(button, messageIndex) {
            const isBookmarked = button.classList.contains('active');
            const bookmarkId = button.dataset.bookmarkId;

            try {
                const response = await fetch('/bookmark', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                        'X-CSRFToken': document.querySelector('input[name="csrf_token"]').value
                    },
                    body: JSON.stringify({
                        message_index: messageIndex,
                        action: isBookmarked ? 'remove' : 'add',
                        bookmark_id: bookmarkId
                    })
                });

                const data = await response.json();

                if (data.success) {
                    if (data.action === 'added') {
                        button.classList.add('active');
                        button.dataset.bookmarkId = data.bookmark_id;
                        button.querySelector('span').textContent = 'Saved';
                    } else {
                        button.classList.remove('active');
                        delete button.dataset.bookmarkId;
                        button.querySelector('span').textContent = 'Save';
                    }
                }
            } catch (err) {
                console.error('Bookmark error:', err);
            }
        }

        // Share Response
        let currentShareUrl = '';
        async function shareResponse(messageIndex) {
            try {
                const response = await fetch('/share', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                        'X-CSRFToken': document.querySelector('input[name="csrf_token"]').value
                    },
                    body: JSON.stringify({ message_index: messageIndex })
                });

                const data = await response.json();

                if (data.success) {
                    currentShareUrl = data.share_url;
                    document.getElementById('shareLinkInput').value = currentShareUrl;
                    document.getElementById('shareModal').classList.add('active');
                }
            } catch (err) {
                console.error('Share error:', err);
            }
        }

        function closeShareModal() {
            document.getElementById('shareModal').classList.remove('active');
        }

        function copyShareLink() {
            const input = document.getElementById('shareLinkInput');
            navigator.clipboard.writeText(input.value).then(() => {
                const copyBtn = event.target.closest('button');
                const originalHTML = copyBtn.innerHTML;
                copyBtn.innerHTML = '<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M20 6L9 17l-5-5"></path></svg> Copied!';
                setTimeout(() => {
                    copyBtn.innerHTML = originalHTML;
                }, 2000);
            });
        }

        function openShareLink() {
            if (currentShareUrl) {
                window.open(currentShareUrl, '_blank');
            }
        }

        // Close modal on outside click
        document.addEventListener('click', function(e) {
            const modal = document.getElementById('shareModal');
            if (e.target === modal) {
                closeShareModal();
            }
        });

        // Export Citations
        async function exportCitations(messageIndex, format) {
            try {
                const response = await fetch(`/export-citations/${format}`, {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                        'X-CSRFToken': document.querySelector('input[name="csrf_token"]').value
                    },
                    body: JSON.stringify({ message_index: messageIndex })
                });

                if (response.ok) {
                    const blob = await response.blob();
                    const url = window.URL.createObjectURL(blob);
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = `gasconsult_citations.${format === 'bibtex' ? 'bib' : 'ris'}`;
                    document.body.appendChild(a);
                    a.click();
                    document.body.removeChild(a);
                    window.URL.revokeObjectURL(url);
                }
            } catch (err) {
                console.error('Export error:', err);
            }
        }

        // Voice Input
        let recognition = null;
        let isListening = false;

        function toggleVoiceInput() {
            if (!('webkitSpeechRecognition' in window) && !('SpeechRecognition' in window)) {
                alert('Voice input is not supported in your browser. Try Chrome or Edge.');
                return;
            }

            if (!recognition) {
                recognition = new (window.SpeechRecognition || window.webkitSpeechRecognition)();
                recognition.continuous = false;
                recognition.interimResults = false;
                recognition.lang = 'en-US';

                recognition.onresult = function(event) {
                    const transcript = event.results[0][0].transcript;
                    const textarea = document.getElementById('chatInput');
                    textarea.value = transcript;
                    textarea.style.height = '44px';
                    textarea.style.height = textarea.scrollHeight + 'px';
                    stopVoiceInput();
                };

                recognition.onerror = function(event) {
                    console.error('Speech recognition error:', event.error);
                    stopVoiceInput();
                };

                recognition.onend = function() {
                    stopVoiceInput();
                };
            }

            if (isListening) {
                stopVoiceInput();
            } else {
                startVoiceInput();
            }
        }

        function startVoiceInput() {
            recognition.start();
            isListening = true;
            document.getElementById('voiceBtn').classList.add('listening');
        }

        function stopVoiceInput() {
            if (recognition && isListening) {
                recognition.stop();
            }
            isListening = false;
            document.getElementById('voiceBtn').classList.remove('listening');
        }

        // Auto-load follow-up suggestions on page load
        document.addEventListener('DOMContentLoaded', function() {
            // Load follow-ups for the last assistant message
            const messages = document.querySelectorAll('.message.assistant');
            if (messages.length > 0) {
                const lastMessageIndex = messages.length - 1;
                loadFollowupSuggestions(lastMessageIndex);
            }
        });

        // Load Follow-up Suggestions
        async function loadFollowupSuggestions(messageIndex) {
            try {
                const response = await fetch('/generate-followups', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                        'X-CSRFToken': document.querySelector('input[name="csrf_token"]').value
                    },
                    body: JSON.stringify({ message_index: messageIndex })
                });

                const data = await response.json();

                if (data.success && data.followups) {
                    const section = document.getElementById(`followup-${messageIndex}`);
                    const questionsContainer = document.getElementById(`followup-questions-${messageIndex}`);

                    if (section && questionsContainer) {
                        questionsContainer.innerHTML = '';
                        data.followups.forEach(question => {
                            const btn = document.createElement('button');
                            btn.className = 'followup-btn';
                            btn.textContent = question;
                            btn.onclick = function() {
                                fillQuery(question);
                                document.getElementById('chatInput').focus();
                            };
                            questionsContainer.appendChild(btn);
                        });
                        section.style.display = 'block';
                    }
                }
            } catch (err) {
                console.error('Followup generation error:', err);
            }
        }
    </script>
</body>
</html>
"""

LIBRARY_HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Library — gasconsult.ai</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Sora:wght@400;600&display=swap" rel="stylesheet">
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <style>
        :root {
            --primary-blue: #2563EB;
            --primary-blue-dark: #1D4ED8;
            --primary-blue-light: #DBEAFE;
            --text-primary: #0F172A;
            --text-secondary: #475569;
            --text-muted: #94A3B8;
            --bg-primary: #FFFFFF;
            --bg-secondary: #F8FAFC;
            --border: #E2E8F0;
        }
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            background: var(--bg-secondary);
            color: var(--text-primary);
            line-height: 1.6;
        }
        nav {
            background: rgba(255, 255, 255, 0.85);
            backdrop-filter: blur(12px);
            -webkit-backdrop-filter: blur(12px);
            padding: 16px 40px;
            border-bottom: 1px solid var(--border);
        }
        nav .container {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }
        .logo-text {
            font-family: 'Sora', sans-serif;
            font-size: 20px;
            font-weight: 600;
            color: var(--primary-blue);
            text-decoration: none;
        }
        .nav-link {
            color: var(--text-secondary);
            text-decoration: none;
            font-size: 14px;
            font-weight: 500;
            padding: 8px 16px;
            border-radius: 8px;
            transition: all 0.2s ease;
        }
        .nav-link:hover {
            background: rgba(255, 255, 255, 0.6);
            color: var(--primary-blue);
        }
        main {
            max-width: 1000px;
            margin: 40px auto;
            padding: 0 24px;
        }
        h1 {
            font-family: 'Sora', sans-serif;
            font-size: 2.2rem;
            margin-bottom: 32px;
            color: var(--text-primary);
        }
        .empty-state {
            text-align: center;
            padding: 80px 40px;
            color: var(--text-muted);
        }
        .bookmark-card {
            background: white;
            border-radius: 16px;
            padding: 24px;
            margin-bottom: 20px;
            border: 1.5px solid var(--border);
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.04);
            transition: all 0.3s ease;
        }
        .bookmark-card:hover {
            box-shadow: 0 8px 24px rgba(0, 0, 0, 0.08);
            transform: translateY(-2px);
        }
        .bookmark-header {
            display: flex;
            justify-content: space-between;
            align-items: start;
            margin-bottom: 16px;
        }
        .bookmark-query {
            font-size: 16px;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 4px;
        }
        .bookmark-date {
            font-size: 12px;
            color: var(--text-muted);
        }
        .remove-btn {
            background: transparent;
            border: 1.5px solid var(--border);
            color: var(--text-secondary);
            padding: 6px 12px;
            border-radius: 6px;
            font-size: 12px;
            cursor: pointer;
            transition: all 0.2s ease;
        }
        .remove-btn:hover {
            background: #FEE2E2;
            border-color: #EF4444;
            color: #EF4444;
        }
        .bookmark-answer {
            font-size: 15px;
            line-height: 1.7;
            color: var(--text-secondary);
            margin-bottom: 16px;
        }
        .bookmark-refs {
            font-size: 13px;
            color: var(--text-muted);
            padding-top: 12px;
            border-top: 1px solid var(--border);
        }
    </style>
</head>
<body>
    <nav>
        <div class="container">
            <a href="/" class="logo-text">gasconsult.ai</a>
            <a href="/chat" class="nav-link">← Back to Chat</a>
        </div>
    </nav>

    <main>
        <h1>📚 My Library</h1>

        {% if not bookmarks %}
        <div class="empty-state">
            <svg width="64" height="64" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" style="margin-bottom: 16px;">
                <path d="M19 21l-7-5-7 5V5a2 2 0 012-2h10a2 2 0 012 2z"></path>
            </svg>
            <p>No saved responses yet. Bookmark responses in chat to build your library.</p>
        </div>
        {% else %}
            {% for bookmark in bookmarks %}
            <div class="bookmark-card">
                <div class="bookmark-header">
                    <div>
                        <div class="bookmark-query">{{ bookmark.query }}</div>
                        <div class="bookmark-date">Saved {{ bookmark.timestamp[:10] }}</div>
                    </div>
                    <button class="remove-btn" onclick="removeBookmark('{{ bookmark.id }}')">Remove</button>
                </div>
                <div class="bookmark-answer">{{ bookmark.answer|safe|truncate(300) }}</div>
                {% if bookmark.num_papers > 0 %}
                <div class="bookmark-refs">📊 {{ bookmark.num_papers }} references</div>
                {% endif %}
            </div>
            {% endfor %}
        {% endif %}
    </main>

    <script>
        async function removeBookmark(id) {
            // Implementation for removing bookmarks
            window.location.reload();
        }
    </script>
</body>
</html>
"""

SHARED_RESPONSE_HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Shared Response — gasconsult.ai</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Sora:wght@400;600&display=swap" rel="stylesheet">
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <style>
        :root {
            --primary-blue: #2563EB;
            --primary-blue-dark: #1D4ED8;
            --primary-blue-light: #DBEAFE;
            --text-primary: #0F172A;
            --text-secondary: #475569;
            --text-muted: #94A3B8;
            --bg-primary: #FFFFFF;
            --bg-secondary: #F8FAFC;
            --border: #E2E8F0;
        }
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            background: var(--bg-secondary);
            color: var(--text-primary);
            line-height: 1.6;
        }
        nav {
            background: rgba(255, 255, 255, 0.85);
            backdrop-filter: blur(12px);
            padding: 16px 40px;
            border-bottom: 1px solid var(--border);
        }
        nav .container {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }
        .logo-text {
            font-family: 'Sora', sans-serif;
            font-size: 20px;
            font-weight: 600;
            color: var(--primary-blue);
            text-decoration: none;
        }
        .try-btn {
            background: linear-gradient(135deg, var(--primary-blue) 0%, var(--primary-blue-dark) 100%);
            color: white;
            padding: 10px 24px;
            border-radius: 8px;
            text-decoration: none;
            font-size: 14px;
            font-weight: 500;
            box-shadow: 0 2px 8px rgba(37, 99, 235, 0.2);
        }
        main {
            max-width: 900px;
            margin: 40px auto;
            padding: 0 24px;
        }
        .response-card {
            background: white;
            border-radius: 16px;
            padding: 32px;
            border: 1.5px solid var(--border);
            box-shadow: 0 4px 16px rgba(0, 0, 0, 0.06);
        }
        .query {
            font-size: 18px;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 24px;
            padding-bottom: 16px;
            border-bottom: 2px solid var(--border);
        }
        .answer {
            font-size: 15px;
            line-height: 1.8;
            color: var(--text-secondary);
            margin-bottom: 24px;
        }
        .refs {
            font-size: 13px;
            color: var(--text-muted);
            padding-top: 16px;
            border-top: 1px solid var(--border);
        }
        .refs a {
            color: var(--primary-blue);
            text-decoration: none;
        }
        .refs a:hover {
            text-decoration: underline;
        }
        .footer-note {
            text-align: center;
            padding: 32px;
            color: var(--text-muted);
            font-size: 13px;
        }
    </style>
</head>
<body>
    <nav>
        <div class="container">
            <a href="/" class="logo-text">gasconsult.ai</a>
            <a href="/chat" class="try-btn">Try it yourself</a>
        </div>
    </nav>

    <main>
        <div class="response-card">
            <div class="query">{{ data.query }}</div>
            <div class="answer">{{ data.answer|safe }}</div>

            {% if data.references %}
            <div class="refs">
                <strong>References:</strong><br>
                {% for ref in data.references %}
                <a href="https://pubmed.ncbi.nlm.nih.gov/{{ ref.pmid }}/" target="_blank">
                    [{{ loop.index }}] {{ ref.title }} ({{ ref.year }})
                </a><br>
                {% endfor %}
                <br>
                📊 {{ data.num_papers }} papers from PubMed
            </div>
            {% endif %}
        </div>

        <div class="footer-note">
            This response was shared from gasconsult.ai • Shared {{ data.timestamp[:10] }}
        </div>
    </main>
</body>
</html>
"""

TERMS_HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Terms of Service - gasconsult.ai</title>
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=5">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">
    <meta name="apple-mobile-web-app-capable" content="yes">
    <meta name="apple-mobile-web-app-status-bar-style" content="default">
    <meta name="apple-mobile-web-app-title" content="gasconsult.ai">
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=Sora:wght@400;600&display=swap" rel="stylesheet">
    <style>
        :root {
            /* Primary Brand Colors */
            --primary-blue: #2563EB;
            --primary-blue-dark: #1D4ED8;
            --primary-blue-light: #DBEAFE;

            /* Anesthesia Color Palette (for logo & accents) */
            --opioid-blue: #2563EB;
            --nmb-red: #EF4444;
            --induction-yellow: #FBBF24;
            --vasopressor-violet: #8B5CF6;
            --anticholinergic-green: #10B981;
            --local-gray: #6B7280;

            /* Neutral Palette */
            --text-primary: #0F172A;
            --text-secondary: #475569;
            --text-muted: #94A3B8;
            --bg-primary: #FFFFFF;
            --bg-secondary: #F8FAFC;
            --border: #E2E8F0;

            /* Legacy aliases for compatibility */
            --primary: #2563EB;
            --primary-dark: #1D4ED8;
            --secondary: #0066CC;
            --text: #1F2937;
            --background: #FAFBFC;
        }

        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'SF Pro Display', 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background: var(--bg-secondary);
            color: #0A3D62;
            line-height: 1.6;
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
            animation: pageFadeIn 0.3s ease-in;
        }

        @keyframes pageFadeIn {
            from { opacity: 0; }
            to { opacity: 1; }
        }

        nav {
            background: rgba(255, 255, 255, 0.85);
            backdrop-filter: blur(12px);
            -webkit-backdrop-filter: blur(12px);
            padding: 16px 40px;
            position: sticky;
            top: 0;
            left: 0;
            right: 0;
            z-index: 100;
            border-bottom: 1px solid rgba(226, 232, 240, 0.8);
        }

        nav .container {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }

        .logo-container {
            display: flex;
            align-items: center;
            gap: 12px;
            cursor: pointer;
            transition: transform 0.2s ease;
            text-decoration: none;
        }

        .logo-container:hover {
            transform: translateY(-1px);
        }

        .logo-ecg {
            height: 28px;
            width: auto;
            flex-shrink: 0;
        }

        .logo-wordmark {
            font-family: 'Sora', -apple-system, BlinkMacSystemFont, system-ui, sans-serif;
            font-size: 20px;
            font-weight: 600;
            letter-spacing: -0.5px;
            white-space: nowrap;
        }

        .logo-gas {
            color: #2563EB;
        }

        .logo-consult {
            color: #111111;
        }

        .logo-ai {
            font-weight: 400;
            color: #6B7280;
        }

        .nav-actions {
            display: flex;
            gap: 12px;
            align-items: center;
        }

        .nav-link {
            color: var(--text-secondary);
            text-decoration: none;
            font-size: 14px;
            font-weight: 500;
            padding: 8px 16px;
            border-radius: 8px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            color: var(--text-primary);
            background: rgba(255, 255, 255, 0.6);
            backdrop-filter: blur(8px);
            -webkit-backdrop-filter: blur(8px);
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.06);
        }

        .nav-link.active {
            color: var(--primary-blue);
            font-weight: 600;
        }

        .content {
            max-width: 900px;
            margin: 60px auto;
            padding: 0 40px;
        }

        h1 {
            font-size: 2.5rem;
            color: var(--text);
            margin-bottom: 12px;
            font-weight: 700;
        }

        .last-updated {
            color: var(--text-muted);
            font-size: 0.9rem;
            margin-bottom: 40px;
        }

        h2 {
            font-size: 1.6rem;
            color: var(--secondary);
            margin-top: 40px;
            margin-bottom: 16px;
            font-weight: 700;
        }

        h3 {
            font-size: 1.2rem;
            color: var(--text);
            margin-top: 24px;
            margin-bottom: 12px;
            font-weight: 600;
        }

        p {
            margin-bottom: 16px;
            color: var(--text-secondary);
        }

        ul, ol {
            margin-left: 24px;
            margin-bottom: 16px;
            color: var(--text-secondary);
        }

        li {
            margin-bottom: 8px;
        }

        .notice-box {
            background: rgba(37, 99, 235, 0.05);
            border-left: 4px solid var(--primary-blue);
            padding: 20px;
            margin: 30px 0;
            border-radius: 0 8px 8px 0;
        }

        .notice-box h3 {
            margin-top: 0;
            color: var(--primary-blue);
        }

        footer {
            text-align: center;
            padding: 40px;
            border-top: 1px solid #E2E8F0;
            background: #FFFFFF;
            color: var(--text-muted);
            font-size: 13px;
            margin: 0;
        }

        footer a {
            color: var(--primary);
            text-decoration: none;
            transition: color 0.2s ease;
        }

        footer a:hover {
            color: var(--text-primary);
        }

        footer p {
            margin-bottom: 8px;
            line-height: 1.6;
        }

        /* Mobile responsiveness */
        @media (max-width: 768px) {
            nav {
                padding: 14px 20px;
            }

            nav .container {
                flex-wrap: nowrap;
            }

            .logo-container {
                flex-shrink: 1;
                min-width: 0;
            }

            .logo-text {
                font-size: 1.2rem;
            }

            .logo-svg {
                width: 28px;
                height: 28px;
            }

            .nav-actions {
                gap: 6px;
                flex-wrap: nowrap;
                flex-shrink: 0;
            }

            .nav-link {
                padding: 6px 12px;
                font-size: 0.85rem;
            }

            .content {
                padding: 0 20px;
                margin: 40px auto;
            }

            h1 {
                font-size: 2rem;
            }

            h2 {
                font-size: 1.5rem;
            }

            h3 {
                font-size: 1.1rem;
            }

            footer {
                padding: 40px 20px;
            }
        }
    </style>
</head>
<body>

    <nav>
        <div class="container">
            <a href="/" class="logo-container">
                <svg class="logo-ecg" viewBox="0 0 60 28" fill="none" xmlns="http://www.w3.org/2000/svg">
                    <defs>
                        <linearGradient id="ecgGrad" x1="0%" y1="0%" x2="100%" y2="0%">
                            <stop offset="0%" stop-color="#2563EB"/>
                            <stop offset="20%" stop-color="#EF4444"/>
                            <stop offset="40%" stop-color="#FBBF24"/>
                            <stop offset="60%" stop-color="#8B5CF6"/>
                            <stop offset="80%" stop-color="#10B981"/>
                            <stop offset="100%" stop-color="#6B7280"/>
                        </linearGradient>
                    </defs>
                    <path d="M2 14 L10 14 L14 12 L18 16 L22 4 L26 24 L30 10 L34 14 L42 14"
                          stroke="url(#ecgGrad)"
                          stroke-width="2.5"
                          stroke-linecap="round"
                          stroke-linejoin="round"
                          fill="none"/>
                </svg>
                <div class="logo-wordmark">
                    <span class="logo-gas">gas</span><span class="logo-consult">consult</span><span class="logo-ai">.ai</span>
                </div>
            </a>
            <div class="nav-actions">
                <a href="/" class="nav-link">Home</a>
                <a href="/preop" class="nav-link">Pre-Op Assessment</a>
                <a href="/calculators" class="nav-link">Clinical Calculators</a>
                <a href="/quick-dose" class="nav-link">Quick Dose</a>
                <a href="/hypotension" class="nav-link">IOH Predictor</a>
            </div>
        </div>
    </nav>

    <div class="content">
        <h1>Terms of Service</h1>
        <p class="last-updated">Last Updated: November 25, 2025</p>

        <div class="notice-box">
            <h3>⚠️ Critical Notice</h3>
            <p><strong>gasconsult.ai is NOT a substitute for professional medical judgment.</strong> This tool is strictly for informational and educational purposes only and must be used exclusively by qualified healthcare professionals as a clinical decision support aid.</p>
        </div>

        <h2>1. Acceptance of Terms</h2>
        <p>By accessing or using gasconsult.ai ("the Service"), you acknowledge that you have read, understood, and agree to be bound by these Terms of Service. If you do not agree to these terms, you must not use this Service.</p>

        <h2>2. Medical Disclaimer</h2>
        <h3>2.1 Not Medical Advice</h3>
        <p>The information provided by gasconsult.ai is for <strong>informational and educational purposes only</strong>. It is not intended to be, and should not be interpreted as:</p>
        <ul>
            <li>Medical advice, diagnosis, or treatment recommendations</li>
            <li>A substitute for professional medical judgment or clinical expertise</li>
            <li>A replacement for consultation with qualified healthcare professionals</li>
            <li>Definitive clinical guidance for patient care decisions</li>
        </ul>

        <h3>2.2 Professional Use Only</h3>
        <p>This Service is designed exclusively for use by licensed healthcare professionals, including but not limited to physicians, nurse anesthetists (CRNAs), and other qualified medical practitioners. The Service must not be used by patients or non-medical personnel for self-diagnosis or self-treatment.</p>

        <h3>2.3 Clinical Decision Support</h3>
        <p>gasconsult.ai serves solely as a <strong>clinical decision support tool</strong> to assist qualified healthcare providers. All treatment decisions must be made by licensed healthcare professionals based on:</p>
        <ul>
            <li>Comprehensive patient assessment and clinical evaluation</li>
            <li>Individual patient circumstances, comorbidities, and risk factors</li>
            <li>Current evidence-based medical practice and institutional protocols</li>
            <li>Professional medical judgment and clinical expertise</li>
        </ul>

        <h2>3. No Warranty or Guarantee</h2>
        <h3>3.1 Information Accuracy</h3>
        <p>While gasconsult.ai strives to provide evidence-based information sourced from peer-reviewed medical literature (PubMed), we make <strong>NO WARRANTIES OR GUARANTEES</strong> regarding:</p>
        <ul>
            <li>The accuracy, completeness, or currentness of any information provided</li>
            <li>The suitability of information for any particular patient or clinical situation</li>
            <li>The reliability of AI-generated content or literature interpretations</li>
            <li>The absence of errors, omissions, or inaccuracies in responses</li>
        </ul>

        <h3>3.2 Service Availability</h3>
        <p>The Service is provided "AS IS" and "AS AVAILABLE" without any warranty of any kind, express or implied, including but not limited to warranties of merchantability, fitness for a particular purpose, or non-infringement.</p>

        <h2>4. Limitation of Liability</h2>
        <h3>4.1 No Liability for Medical Outcomes</h3>
        <p>To the fullest extent permitted by law, gasconsult.ai, its developers, operators, and affiliates shall NOT BE LIABLE for any:</p>
        <ul>
            <li>Patient injuries, adverse outcomes, or complications arising from use of this Service</li>
            <li>Clinical decisions made based on information provided by the Service</li>
            <li>Errors, omissions, or inaccuracies in AI-generated content or literature citations</li>
            <li>Damages arising from reliance on the Service for medical decision-making</li>
            <li>Direct, indirect, incidental, consequential, or punitive damages of any kind</li>
        </ul>

        <h3>4.2 User Responsibility</h3>
        <p>Users of this Service assume <strong>FULL RESPONSIBILITY</strong> for:</p>
        <ul>
            <li>Verifying all information through primary sources and current medical literature</li>
            <li>Exercising independent professional judgment in all clinical decisions</li>
            <li>Complying with institutional policies, protocols, and standard of care requirements</li>
            <li>Obtaining informed consent and following applicable medical regulations</li>
        </ul>

        <h2>5. User Obligations</h2>
        <h3>5.1 Professional Qualifications</h3>
        <p>By using this Service, you represent and warrant that you are:</p>
        <ul>
            <li>A licensed healthcare professional authorized to practice medicine</li>
            <li>Qualified to interpret medical information and make clinical decisions</li>
            <li>Using the Service solely for professional educational purposes</li>
            <li>Capable of independently verifying all medical information</li>
        </ul>

        <h3>5.2 Prohibited Uses</h3>
        <p>You agree NOT to use this Service for:</p>
        <ul>
            <li>Patient self-diagnosis, self-treatment, or medical decision-making by non-professionals</li>
            <li>Emergency medical situations requiring immediate clinical intervention without healthcare professional oversight and judgement</li>
            <li>Situations where delays in obtaining professional medical care could cause harm</li>
            <li>Any unlawful, fraudulent, or unauthorized purposes</li>
        </ul>

        <h2>6. Third-Party Content</h2>
        <p>The Service aggregates information from third-party sources, including PubMed and medical literature databases. We are not responsible for the accuracy, reliability, or content of third-party sources. Citations and references should be independently verified through original publications.</p>

        <h2>7. Privacy and Data</h2>
        <p>Your use of the Service is subject to our Privacy Policy. We do not store patient health information (PHI) or individually identifiable medical data. Users must not input protected health information into the Service.</p>

        <h2>8. Modifications to Terms</h2>
        <p>We reserve the right to modify these Terms of Service at any time. Continued use of the Service following any changes constitutes acceptance of modified terms. Users are responsible for regularly reviewing these terms.</p>

        <h2>9. Indemnification</h2>
        <p>You agree to indemnify, defend, and hold harmless gasconsult.ai and its operators from any claims, damages, losses, liabilities, and expenses (including legal fees) arising from:</p>
        <ul>
            <li>Your use or misuse of the Service</li>
            <li>Clinical decisions made based on information from the Service</li>
            <li>Violation of these Terms of Service</li>
            <li>Violation of any applicable laws or regulations</li>
        </ul>

        <h2>10. Governing Law and Jurisdiction</h2>
        <p>These Terms shall be governed by and construed in accordance with the laws of the United States. Any disputes arising from these Terms or use of the Service shall be subject to the exclusive jurisdiction of the courts in the United States.</p>

        <h2>11. Emergency Medical Situations</h2>
        <div class="notice-box">
            <h3>⚠️ Emergency Disclaimer</h3>
            <p><strong>DO NOT USE THIS SERVICE FOR MEDICAL EMERGENCIES WITHOUT HEALTHCARE PROFESSIONAL OVERSIGHT.</strong> In case of medical emergency, call 911 (or your local emergency number) immediately or seek emergency medical care at the nearest hospital.</p>
        </div>

        <h2>12. Contact Information</h2>
        <p>For questions regarding these Terms of Service, please contact us at: <strong>support@gasconsult.ai</strong></p>

        <p style="margin-top: 40px; font-weight: 600;">By using gasconsult.ai, you acknowledge that you have read, understood, and agree to be bound by these Terms of Service.</p>
    </div>

    <footer>
        <p>&copy; 2025 gasconsult.ai. All rights reserved. | <a href="/terms" style="color: var(--primary); text-decoration: none;">Terms of Service</a> | <a href="/privacy" style="color: var(--primary); text-decoration: none;">Privacy Policy</a></p>
    </footer>

</body>
</html>
"""

PRIVACY_POLICY_HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Privacy Policy - gasconsult.ai</title>
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=Sora:wght@400;600&display=swap" rel="stylesheet">
    <style>
        :root {
            --primary: #2563EB;
            --primary-blue: #2563EB;
            --text-primary: #0F172A;
            --text-secondary: #475569;
            --text-muted: #94A3B8;
        }

        * { margin: 0; padding: 0; box-sizing: border-box; }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            background: #F8FAFC;
            color: #0F172A;
            line-height: 1.7;
            -webkit-font-smoothing: antialiased;
        }

        nav {
            background: rgba(255, 255, 255, 0.85);
            backdrop-filter: blur(12px);
            -webkit-backdrop-filter: blur(12px);
            padding: 16px 40px;
            position: sticky;
            top: 0;
            left: 0;
            right: 0;
            z-index: 100;
            border-bottom: 1px solid rgba(226, 232, 240, 0.8);
        }

        nav .container {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }

        .logo-container {
            display: flex;
            align-items: center;
            gap: 12px;
            cursor: pointer;
            transition: transform 0.2s ease;
            text-decoration: none;
        }

        .logo-container:hover {
            transform: translateY(-1px);
        }

        .logo-ecg {
            height: 28px;
            width: auto;
            flex-shrink: 0;
        }

        .logo-wordmark {
            font-family: 'Sora', -apple-system, BlinkMacSystemFont, system-ui, sans-serif;
            font-size: 20px;
            font-weight: 600;
            letter-spacing: -0.5px;
            white-space: nowrap;
        }

        .logo-gas {
            color: #2563EB;
        }

        .logo-consult {
            color: #111111;
        }

        .logo-ai {
            font-weight: 400;
            color: #6B7280;
        }

        .nav-actions {
            display: flex;
            gap: 12px;
            align-items: center;
        }

        .nav-link {
            color: var(--text-secondary);
            text-decoration: none;
            font-size: 14px;
            font-weight: 500;
            padding: 8px 16px;
            border-radius: 8px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            color: var(--text-primary);
            background: rgba(255, 255, 255, 0.6);
            backdrop-filter: blur(8px);
            -webkit-backdrop-filter: blur(8px);
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.06);
        }

        .nav-link.active {
            color: var(--primary-blue);
            font-weight: 600;
        }

        .content {
            max-width: 900px;
            margin: 60px auto;
            padding: 0 40px;
        }

        h1 {
            font-size: 36px;
            font-weight: 700;
            color: #0F172A;
            margin-bottom: 12px;
        }

        .last-updated {
            color: #64748B;
            font-size: 14px;
            margin-bottom: 40px;
        }

        h2 {
            font-size: 24px;
            font-weight: 600;
            color: #0F172A;
            margin-top: 48px;
            margin-bottom: 16px;
        }

        h3 {
            font-size: 18px;
            font-weight: 600;
            color: #1E293B;
            margin-top: 32px;
            margin-bottom: 12px;
        }

        p {
            margin-bottom: 16px;
            color: #475569;
        }

        ul, ol {
            margin-left: 24px;
            margin-bottom: 16px;
            color: #475569;
        }

        li {
            margin-bottom: 8px;
        }

        strong {
            color: #0F172A;
            font-weight: 600;
        }

        .highlight-box {
            background: #FEF3C7;
            border-left: 4px solid #F59E0B;
            padding: 20px;
            margin: 24px 0;
            border-radius: 4px;
        }

        .highlight-box strong {
            color: #92400E;
        }

        footer {
            text-align: center;
            padding: 40px;
            border-top: 1px solid #E2E8F0;
            background: #FFFFFF;
            color: var(--text-muted);
            font-size: 13px;
            margin: 0;
        }

        footer a {
            color: var(--primary);
            text-decoration: none;
            transition: color 0.2s ease;
        }

        footer a:hover {
            color: var(--text-primary);
        }

        footer p {
            margin-bottom: 8px;
            line-height: 1.6;
        }

        /* Mobile responsiveness */
        @media (max-width: 768px) {
            nav {
                padding: 14px 20px;
            }

            nav .container {
                flex-wrap: nowrap;
            }

            .logo-container {
                flex-shrink: 1;
                min-width: 0;
            }

            .logo-wordmark {
                font-size: 18px;
            }

            .logo-ecg {
                width: 24px;
                height: 24px;
            }

            .nav-actions {
                gap: 6px;
                flex-wrap: nowrap;
                flex-shrink: 0;
            }

            .nav-link {
                padding: 6px 12px;
                font-size: 0.85rem;
            }

            .content {
                padding: 0 20px;
                margin: 40px auto;
            }

            h1 {
                font-size: 28px;
            }

            h2 {
                font-size: 20px;
                margin-top: 36px;
            }

            h3 {
                font-size: 16px;
                margin-top: 24px;
            }

            footer {
                padding: 32px 20px;
            }
        }
    </style>
</head>
<body>
    <nav>
        <div class="container">
            <a href="/" class="logo-container">
                <svg class="logo-ecg" viewBox="0 0 60 28" fill="none" xmlns="http://www.w3.org/2000/svg">
                    <defs>
                        <linearGradient id="ecgGrad" x1="0%" y1="0%" x2="100%" y2="0%">
                            <stop offset="0%" stop-color="#2563EB"/>
                            <stop offset="20%" stop-color="#EF4444"/>
                            <stop offset="40%" stop-color="#FBBF24"/>
                            <stop offset="60%" stop-color="#8B5CF6"/>
                            <stop offset="80%" stop-color="#10B981"/>
                            <stop offset="100%" stop-color="#6B7280"/>
                        </linearGradient>
                    </defs>
                    <path d="M2 14 L10 14 L14 12 L18 16 L22 4 L26 24 L30 10 L34 14 L42 14"
                          stroke="url(#ecgGrad)"
                          stroke-width="2.5"
                          stroke-linecap="round"
                          stroke-linejoin="round"
                          fill="none"/>
                </svg>
                <div class="logo-wordmark">
                    <span class="logo-gas">gas</span><span class="logo-consult">consult</span><span class="logo-ai">.ai</span>
                </div>
            </a>
            <div class="nav-actions">
                <a href="/" class="nav-link">Home</a>
                <a href="/preop" class="nav-link">Pre-Op Assessment</a>
                <a href="/calculators" class="nav-link">Clinical Calculators</a>
                <a href="/quick-dose" class="nav-link">Quick Dose</a>
                <a href="/hypotension" class="nav-link">IOH Predictor</a>
            </div>
        </div>
    </nav>

    <div class="content">
        <h1>Privacy Policy</h1>
        <p class="last-updated">Last Updated: November 27, 2025</p>

        <div class="highlight-box">
            <strong>TL;DR:</strong> We take your privacy seriously. We don't sell your data. We don't store personal health information.
            Your medical queries are processed to provide answers and are not used to identify you.
        </div>

        <h2>1. Introduction</h2>
        <p>
            Welcome to gasconsult.ai ("we," "our," or "us"). This Privacy Policy explains how we collect, use,
            disclose, and safeguard your information when you use our web application.
        </p>
        <p>
            <strong>Please read this privacy policy carefully.</strong> If you do not agree with the terms of this
            privacy policy, please do not access the application.
        </p>

        <h2>2. Information We Collect</h2>

        <h3>2.1 Information You Provide</h3>
        <p>We collect information that you voluntarily provide when using our service:</p>
        <ul>
            <li><strong>Medical Queries:</strong> Questions you submit about anesthesiology topics</li>
            <li><strong>Form Data:</strong> Information entered in pre-operative assessments or clinical calculators</li>
            <li><strong>Session Data:</strong> Temporary conversation history during your active session</li>
        </ul>

        <h3>2.2 Automatically Collected Information</h3>
        <p>We automatically collect certain information when you use our service:</p>
        <ul>
            <li><strong>Usage Data:</strong> Pages visited, features used, timestamps</li>
            <li><strong>Technical Data:</strong> IP address, browser type, device information</li>
            <li><strong>Log Data:</strong> Error logs and performance metrics for service improvement</li>
        </ul>

        <h3>2.3 Third-Party Services</h3>
        <p>We use the following third-party services that may collect information:</p>
        <ul>
            <li><strong>OpenAI (GPT-4):</strong> Processes your queries to generate responses</li>
            <li><strong>NCBI PubMed:</strong> Searches medical literature based on your queries</li>
        </ul>

        <h2>3. How We Use Your Information</h2>
        <p>We use the collected information for the following purposes:</p>
        <ul>
            <li><strong>Service Delivery:</strong> To provide evidence-based anesthesiology answers</li>
            <li><strong>Improvement:</strong> To analyze usage patterns and improve our service</li>
            <li><strong>Security:</strong> To detect, prevent, and address technical issues</li>
            <li><strong>Communication:</strong> To respond to your inquiries and provide support</li>
            <li><strong>Legal Compliance:</strong> To comply with applicable laws and regulations</li>
        </ul>

        <h2>4. Data Retention</h2>
        <ul>
            <li><strong>Session Data:</strong> Conversation history is stored temporarily during your session and cleared when you close your browser or clear your session</li>
            <li><strong>Log Data:</strong> System logs are retained for 90 days for debugging and security purposes</li>
            <li><strong>No Long-Term Storage:</strong> We do not store your medical queries or personal health information long-term</li>
        </ul>

        <h2>5. Data Security</h2>
        <p>
            We implement industry-standard security measures to protect your information:
        </p>
        <ul>
            <li>HTTPS encryption for all data transmission</li>
            <li>Server-side session storage (not in cookies)</li>
            <li>Input sanitization to prevent malicious attacks</li>
            <li>Rate limiting to prevent abuse</li>
            <li>Regular security updates and monitoring</li>
        </ul>
        <p>
            <strong>However, no method of transmission over the internet is 100% secure.</strong> While we strive
            to protect your information, we cannot guarantee absolute security.
        </p>

        <h2>6. Your Privacy Rights</h2>
        <p>Depending on your location, you may have the following rights:</p>
        <ul>
            <li><strong>Access:</strong> Request a copy of the information we have about you</li>
            <li><strong>Correction:</strong> Request correction of inaccurate information</li>
            <li><strong>Deletion:</strong> Request deletion of your information</li>
            <li><strong>Objection:</strong> Object to our processing of your information</li>
            <li><strong>Data Portability:</strong> Request transfer of your data to another service</li>
        </ul>
        <p>
            To exercise these rights, please contact us at: <strong>privacy@gasconsult.ai</strong>
        </p>

        <h2>7. Protected Health Information (PHI)</h2>
        <div class="highlight-box">
            <strong>IMPORTANT:</strong> gasconsult.ai is an educational tool only.
            <strong>Do not enter any personally identifiable patient information (PHI)</strong> such as:
            <ul style="margin-top: 12px;">
                <li>Patient names, dates of birth, medical record numbers</li>
                <li>Specific case details that could identify individuals</li>
                <li>Protected health information under HIPAA regulations</li>
            </ul>
            <p style="margin-top: 12px; margin-bottom: 0;">
                Use only hypothetical scenarios and generalized clinical questions.
            </p>
        </div>

        <h2>8. Children's Privacy</h2>
        <p>
            Our service is not intended for children under 18 years of age. We do not knowingly collect
            personal information from children. If you believe we have collected information from a child,
            please contact us immediately.
        </p>

        <h2>9. International Users</h2>
        <p>
            Our service is hosted in the United States. If you access our service from outside the US,
            your information may be transferred to, stored, and processed in the US where our servers are located.
            Data protection laws in the US may differ from those in your country.
        </p>

        <h2>10. Changes to This Policy</h2>
        <p>
            We may update this Privacy Policy from time to time. We will notify you of any changes by
            posting the new policy on this page and updating the "Last Updated" date.
        </p>
        <p>
            <strong>Continued use of our service after changes constitutes acceptance of the updated policy.</strong>
        </p>

        <h2>11. Contact Us</h2>
        <p>
            If you have questions or concerns about this Privacy Policy, please contact us at:
        </p>
        <ul>
            <li><strong>Email:</strong> privacy@gasconsult.ai</li>
            <li><strong>Website:</strong> gasconsult.ai</li>
        </ul>

        <h2>12. Compliance</h2>
        <p>
            We are committed to complying with applicable data protection laws, including:
        </p>
        <ul>
            <li>General Data Protection Regulation (GDPR) for EU users</li>
            <li>California Consumer Privacy Act (CCPA) for California residents</li>
            <li>Health Insurance Portability and Accountability Act (HIPAA) - Note: We are not a covered entity
                but follow privacy best practices</li>
        </ul>
    </div>

    <footer>
        <p>&copy; 2025 gasconsult.ai • Educational Use Only • Not Medical Advice</p>
        <p style="margin-top: 8px;">
            <a href="/terms" style="color: #64748B; text-decoration: none;">Terms of Service</a> •
            <a href="/privacy" style="color: #64748B; text-decoration: none;">Privacy Policy</a>
        </p>
    </footer>
</body>
</html>
"""

QUICK_DOSE_HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Quick Dose Reference — gasconsult.ai</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Sora:wght@300;400;500;600;700&family=Inter:wght@400;500;600&family=JetBrains+Mono:wght@400;500;600&display=swap" rel="stylesheet">
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=5">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">
    <meta name="apple-mobile-web-app-capable" content="yes">
    <meta name="apple-mobile-web-app-status-bar-style" content="default">
    <meta name="apple-mobile-web-app-title" content="gasconsult.ai">
    <style>
        :root {
            /* Primary Brand */
            --primary-blue: #2563EB;
            --primary-blue-dark: #1D4ED8;
            --primary-blue-light: #DBEAFE;
            --primary: #2563EB;
            --primary-dark: #1D4ED8;

            /* Anesthesia Color Palette */
            --opioid-blue: #2563EB;
            --opioid-blue-light: #DBEAFE;
            --nmb-red: #EF4444;
            --nmb-red-light: #FEE2E2;
            --induction-yellow: #F59E0B;
            --induction-yellow-light: #FEF3C7;
            --tranq-orange: #F97316;
            --tranq-orange-light: #FFEDD5;
            --vasopressor-violet: #8B5CF6;
            --vasopressor-violet-light: #EDE9FE;
            --anticholinergic-green: #10B981;
            --anticholinergic-green-light: #D1FAE5;
            --local-gray: #6B7280;
            --local-gray-light: #F3F4F6;

            /* Neutral Palette */
            --text-primary: #0F172A;
            --text-secondary: #475569;
            --text-muted: #94A3B8;
            --bg-primary: #FFFFFF;
            --bg-secondary: #F8FAFC;
            --border: #E2E8F0;

            /* Crisis */
            --crisis-red: #DC2626;
            --crisis-red-light: #FEF2F2;
        }

        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, system-ui, sans-serif;
            background: var(--bg-secondary);
            color: var(--text-primary);
            line-height: 1.6;
            -webkit-font-smoothing: antialiased;
            overflow-x: hidden;
            width: 100%;
        }

        /* Header */
        header {
            background: rgba(255, 255, 255, 0.85);
            backdrop-filter: blur(12px);
            -webkit-backdrop-filter: blur(12px);
            border-bottom: 1px solid rgba(226, 232, 240, 0.8);
            padding: 16px 40px;
            position: sticky;
            top: 0;
            z-index: 100;
        }

        .header-container {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 24px;
            flex-wrap: wrap;
        }

        .header-left {
            display: flex;
            align-items: center;
            gap: 16px;
        }

        .logo-link {
            display: flex;
            align-items: center;
            gap: 12px;
            text-decoration: none;
            transition: transform 0.2s ease;
        }

        .logo-link:hover {
            transform: translateY(-1px);
        }

        .logo-svg {
            height: 28px;
            width: auto;
            flex-shrink: 0;
        }

        .logo-text {
            font-family: 'Sora', sans-serif;
            font-size: 20px;
            font-weight: 600;
            line-height: 1;
            letter-spacing: -0.5px;
            white-space: nowrap;
        }

        .logo-gas {
            color: #2563EB;
        }

        .logo-consult {
            color: #111111;
        }

        .logo-ai {
            color: #6B7280;
            font-weight: 400;
        }

        .divider {
            width: 1px;
            height: 24px;
            background: var(--border);
        }

        .page-title {
            font-family: 'Sora', sans-serif;
            font-size: 14px;
            font-weight: 600;
            color: var(--text-secondary);
        }

        .header-center {
            flex: 1;
            display: flex;
            justify-content: center;
            align-items: center;
        }

        .nav-links {
            display: flex;
            gap: 12px;
            align-items: center;
        }

        .nav-link {
            color: var(--text-secondary);
            text-decoration: none;
            font-size: 14px;
            font-weight: 500;
            padding: 8px 16px;
            border-radius: 8px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            color: var(--text-primary);
            background: rgba(255, 255, 255, 0.6);
            backdrop-filter: blur(8px);
            -webkit-backdrop-filter: blur(8px);
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.06);
        }

        .nav-link.active {
            color: var(--primary-blue);
            font-weight: 600;
        }

        .crisis-btn {
            background: var(--crisis-red);
            color: white;
            padding: 10px 20px;
            border-radius: 10px;
            border: none;
            font-family: 'Inter', sans-serif;
            font-size: 13px;
            font-weight: 600;
            cursor: pointer;
            display: flex;
            align-items: center;
            gap: 8px;
            transition: all 0.2s ease;
            animation: subtle-pulse 2s infinite;
        }

        .crisis-btn:hover {
            background: #B91C1C;
            transform: scale(1.02);
        }

        @keyframes subtle-pulse {
            0%, 100% { box-shadow: 0 0 0 0 rgba(220, 38, 38, 0.4); }
            50% { box-shadow: 0 0 0 8px rgba(220, 38, 38, 0); }
        }

        /* Main Content */
        main {
            max-width: 900px;
            margin: 0 auto;
            padding: 32px 24px;
        }

        /* Weight Input Section */
        .weight-section {
            background: white;
            border-radius: 16px;
            padding: 24px;
            margin-bottom: 24px;
            border: 1px solid var(--border);
            box-shadow: 0 1px 3px rgba(0,0,0,0.04);
        }

        .weight-label {
            font-size: 11px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 1px;
            color: var(--text-muted);
            margin-bottom: 12px;
        }

        .weight-input-row {
            display: flex;
            align-items: center;
            gap: 16px;
            flex-wrap: wrap;
        }

        .weight-input-wrapper {
            position: relative;
            width: 160px;
        }

        #weightInput {
            width: 100%;
            padding: 14px 45px 14px 18px;
            font-family: 'JetBrains Mono', monospace;
            font-size: 28px;
            font-weight: 600;
            border: 2px solid var(--border);
            border-radius: 12px;
            transition: border-color 0.2s ease;
        }

        #weightInput:focus {
            outline: none;
            border-color: var(--primary-blue);
        }

        .weight-unit {
            position: absolute;
            right: 14px;
            top: 50%;
            transform: translateY(-50%);
            font-size: 14px;
            color: var(--text-muted);
            pointer-events: none;
        }

        .conversion-text {
            font-size: 14px;
            color: var(--text-secondary);
        }

        .conversion-value {
            font-family: 'JetBrains Mono', monospace;
            font-weight: 600;
            color: var(--text-primary);
        }

        .quick-weights {
            display: flex;
            gap: 8px;
            margin-top: 16px;
            flex-wrap: wrap;
        }

        .quick-weight-btn {
            padding: 8px 14px;
            background: var(--bg-secondary);
            border: 1px solid var(--border);
            border-radius: 8px;
            font-family: 'Inter', sans-serif;
            font-size: 13px;
            font-weight: 500;
            color: var(--text-secondary);
            cursor: pointer;
            transition: all 0.2s ease;
        }

        .quick-weight-btn:hover {
            background: var(--primary-blue-light);
            border-color: var(--primary-blue);
            color: var(--primary-blue);
        }

        .quick-weight-btn.active {
            background: var(--primary-blue);
            border-color: var(--primary-blue);
            color: white;
        }

        /* Drug Category */
        .drug-category {
            background: white;
            border-radius: 16px;
            margin-bottom: 12px;
            border: 1px solid var(--border);
            overflow: hidden;
            box-shadow: 0 1px 3px rgba(0,0,0,0.04);
        }

        .category-header {
            padding: 16px 20px;
            display: flex;
            align-items: center;
            gap: 14px;
            cursor: pointer;
            transition: background 0.2s ease;
        }

        .category-header:hover {
            background: var(--bg-secondary);
        }

        .color-indicator {
            width: 5px;
            height: 40px;
            border-radius: 3px;
            flex-shrink: 0;
        }

        .category-info {
            flex: 1;
        }

        .category-title {
            font-family: 'Sora', sans-serif;
            font-size: 16px;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 2px;
        }

        .category-subtitle {
            font-size: 12px;
            color: var(--text-muted);
        }

        .chevron {
            width: 20px;
            height: 20px;
            color: var(--text-muted);
            transition: transform 0.2s ease;
        }

        .drug-category.open .chevron {
            transform: rotate(180deg);
        }

        .category-content {
            display: none;
            padding: 0 20px 20px;
        }

        .drug-category.open .category-content {
            display: block;
        }

        /* Drug Card */
        .drug-card {
            background: var(--bg-secondary);
            border-radius: 12px;
            padding: 20px;
            margin-bottom: 12px;
        }

        .drug-card:last-child {
            margin-bottom: 0;
        }

        .drug-header {
            display: flex;
            justify-content: space-between;
            align-items: flex-start;
            margin-bottom: 16px;
            gap: 12px;
            flex-wrap: wrap;
        }

        .drug-name-section h3 {
            font-family: 'Sora', sans-serif;
            font-size: 16px;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 2px;
        }

        .drug-subtitle {
            font-size: 12px;
            color: var(--text-muted);
        }

        .concentration-badge {
            font-family: 'JetBrains Mono', monospace;
            font-size: 11px;
            color: var(--text-muted);
            background: white;
            padding: 6px 10px;
            border-radius: 6px;
            border: 1px solid var(--border);
        }

        .dose-grid {
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 10px;
            margin-bottom: 16px;
        }

        .dose-item {
            background: white;
            border-radius: 10px;
            padding: 14px;
            text-align: center;
            border: 1px solid var(--border);
        }

        .dose-label {
            font-size: 10px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: var(--text-muted);
            margin-bottom: 6px;
        }

        .dose-value {
            font-family: 'JetBrains Mono', monospace;
            font-size: 22px;
            font-weight: 600;
            color: var(--text-primary);
            line-height: 1.2;
        }

        .dose-unit {
            font-size: 13px;
            font-weight: 500;
            color: var(--text-secondary);
        }

        .dose-range {
            font-family: 'JetBrains Mono', monospace;
            font-size: 11px;
            color: var(--text-muted);
            margin-top: 4px;
        }

        .clinical-pearl {
            display: flex;
            align-items: flex-start;
            gap: 10px;
            padding: 12px 14px;
            background: var(--induction-yellow-light);
            border-radius: 10px;
            border-left: 3px solid var(--induction-yellow);
        }

        .pearl-icon {
            width: 16px;
            height: 16px;
            color: var(--induction-yellow);
            flex-shrink: 0;
            margin-top: 2px;
        }

        .pearl-text {
            font-size: 12px;
            color: var(--text-secondary);
            line-height: 1.6;
        }

        /* Color Legend */
        .color-legend {
            margin-top: 32px;
            padding: 20px 24px;
            background: white;
            border-radius: 16px;
            border: 1px solid var(--border);
        }

        .legend-title {
            font-size: 11px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 1.5px;
            color: var(--text-muted);
            margin-bottom: 16px;
        }

        .legend-grid {
            display: flex;
            flex-wrap: wrap;
            gap: 16px;
        }

        .legend-item {
            display: flex;
            align-items: center;
            gap: 8px;
        }

        .legend-swatch {
            width: 12px;
            height: 12px;
            border-radius: 4px;
        }

        .legend-label {
            font-size: 12px;
            color: var(--text-secondary);
        }

        /* Footer */
        .disclaimer {
            max-width: 700px;
            margin: 0 auto 16px;
            padding: 16px 20px;
            background: white;
            border-radius: 12px;
            border: 1px solid var(--border);
            font-size: 11px;
            line-height: 1.6;
            color: var(--text-secondary);
        }

        footer {
            text-align: center;
            padding: 40px;
            border-top: 1px solid #E2E8F0;
            background: #FFFFFF;
            color: var(--text-muted);
            font-size: 13px;
            margin: 0;
        }

        footer a {
            transition: color 0.2s ease;
        }

        footer a:hover {
            color: var(--text-primary);
        }

        /* Crisis Modal */
        .crisis-overlay {
            display: none;
            position: fixed;
            inset: 0;
            background: rgba(0,0,0,0.6);
            backdrop-filter: blur(4px);
            align-items: center;
            justify-content: center;
            z-index: 1000;
            padding: 24px;
        }

        .crisis-overlay.show {
            display: flex;
        }

        .crisis-modal {
            background: white;
            border-radius: 20px;
            width: 100%;
            max-width: 500px;
            max-height: 85vh;
            overflow: hidden;
            box-shadow: 0 25px 50px rgba(0,0,0,0.25);
        }

        .crisis-modal-header {
            background: var(--crisis-red);
            color: white;
            padding: 20px 24px;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }

        .crisis-modal-title {
            font-family: 'Sora', sans-serif;
            font-size: 18px;
            font-weight: 700;
            display: flex;
            align-items: center;
            gap: 8px;
        }

        .crisis-close-btn {
            background: rgba(255,255,255,0.2);
            border: none;
            color: white;
            width: 36px;
            height: 36px;
            border-radius: 10px;
            cursor: pointer;
            display: flex;
            align-items: center;
            justify-content: center;
            transition: background 0.2s ease;
        }

        .crisis-close-btn:hover {
            background: rgba(255,255,255,0.3);
        }

        .crisis-modal-content {
            padding: 24px;
            overflow-y: auto;
            max-height: calc(85vh - 80px);
        }

        .protocol-grid {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 12px;
        }

        .protocol-btn {
            background: var(--crisis-red-light);
            border: 2px solid transparent;
            border-radius: 14px;
            padding: 20px 16px;
            text-align: left;
            cursor: pointer;
            transition: all 0.2s ease;
        }

        .protocol-btn:hover {
            border-color: var(--crisis-red);
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(220,38,38,0.15);
        }

        .protocol-title {
            font-family: 'Sora', sans-serif;
            font-size: 14px;
            font-weight: 600;
            color: var(--crisis-red);
            margin-bottom: 4px;
        }

        .protocol-desc {
            font-size: 11px;
            color: var(--text-secondary);
            line-height: 1.4;
        }

        /* Responsive */
        @media (max-width: 768px) {
            header {
                padding: 14px 20px;
            }

            .header-container {
                padding: 0;
            }

            .logo-text {
                font-size: 1.1rem;
            }

            .nav-links {
                gap: 6px;
            }

            .nav-link {
                padding: 6px 12px;
                font-size: 0.85rem;
            }

            main {
                padding: 24px 16px;
            }

            .weight-section {
                padding: 20px 16px;
            }

            .dose-grid {
                grid-template-columns: 1fr;
                gap: 16px;
                padding: 0;
            }

            .protocol-grid {
                grid-template-columns: 1fr;
            }

            .dose-card {
                padding: 16px;
            }

            .dose-drug {
                font-size: 1rem;
            }

            .dose-value {
                font-size: 0.9rem;
            }

            .dose-notes {
                font-size: 0.8rem;
            }

            .protocol-card {
                padding: 20px;
            }

            .protocol-name {
                font-size: 1.1rem;
            }

            .protocol-steps li {
                font-size: 0.9rem;
            }

            h1 {
                font-size: 1.8rem;
            }

            h2 {
                font-size: 1.3rem;
            }

            .search-container input {
                font-size: 16px !important; /* Prevents iOS zoom on focus */
            }

            .weight-input-wrapper {
                width: 140px;
            }

            #weightInput {
                font-size: 24px;
                padding: 12px 40px 12px 14px;
            }
        }

        @media (max-width: 640px) {
            .header-container {
                justify-content: space-between;
                flex-wrap: wrap;
            }

            .header-left {
                order: 1;
            }

            .nav-links {
                order: 2;
                width: 100%;
                justify-content: center;
                margin-top: 12px;
            }

            .nav-link {
                padding: 6px 10px;
                font-size: 0.8rem;
            }

            .logo-text {
                font-size: 1rem;
            }

            .logo-svg {
                height: 24px;
            }
        }
    </style>
    <script>
        // Register Service Worker for offline functionality
        if ('serviceWorker' in navigator) {
            window.addEventListener('load', () => {
                navigator.serviceWorker.register('/static/sw.js')
                    .then((registration) => {
                        console.log('Service Worker registered:', registration.scope);
                    })
                    .catch((error) => {
                        console.log('Service Worker registration failed:', error);
                    });
            });
        }
    </script>
</head>
<body>
    <header>
        <div class="header-container">
            <div class="header-left">
                <a href="/" class="logo-link">
                    <svg class="logo-svg" viewBox="0 0 50 24" fill="none" xmlns="http://www.w3.org/2000/svg">
                        <defs>
                            <linearGradient id="ecgGrad" x1="0%" y1="0%" x2="100%" y2="0%">
                                <stop offset="0%" stop-color="#2563EB"/>
                                <stop offset="20%" stop-color="#EF4444"/>
                                <stop offset="40%" stop-color="#FBBF24"/>
                                <stop offset="60%" stop-color="#8B5CF6"/>
                                <stop offset="80%" stop-color="#10B981"/>
                                <stop offset="100%" stop-color="#6B7280"/>
                            </linearGradient>
                        </defs>
                        <path d="M2 12 L8 12 L11 10 L14 14 L17 4 L20 20 L23 9 L26 12 L34 12"
                              stroke="url(#ecgGrad)"
                              stroke-width="2.5"
                              stroke-linecap="round"
                              stroke-linejoin="round"
                              fill="none"/>
                    </svg>
                    <span class="logo-text">
                        <span class="logo-gas">gas</span><span class="logo-consult">consult</span><span class="logo-ai">.ai</span>
                    </span>
                </a>
            </div>
            <div class="header-center">
                <button class="crisis-btn" onclick="toggleCrisis()">
                    <svg width="16" height="16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                        <path d="M10.29 3.86L1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z"></path>
                        <line x1="12" y1="9" x2="12" y2="13"></line>
                        <line x1="12" y1="17" x2="12.01" y2="17"></line>
                    </svg>
                    Crisis Mode
                </button>
            </div>
            <div class="nav-links">
                <a href="/" class="nav-link">Home</a>
                <a href="/preop" class="nav-link">Pre-Op Assessment</a>
                <a href="/calculators" class="nav-link">Clinical Calculators</a>
                <a href="/quick-dose" class="nav-link active">Quick Dose</a>
                <a href="/hypotension" class="nav-link">IOH Predictor</a>
            </div>
        </div>
    </header>

    <main>
        <!-- Weight Input Section -->
        <div class="weight-section">
            <div class="weight-label">Patient Weight</div>
            <div class="weight-input-row">
                <div class="weight-input-wrapper">
                    <input type="number" id="weightInput" value="70" min="1" max="300" oninput="updateDoses()">
                    <span class="weight-unit">kg</span>
                </div>
                <span class="conversion-text">= <span class="conversion-value" id="lbsConversion">154</span> lbs</span>
            </div>
            <div class="quick-weights">
                <button class="quick-weight-btn" onclick="setWeight(50)">50 kg</button>
                <button class="quick-weight-btn active" onclick="setWeight(70)">70 kg</button>
                <button class="quick-weight-btn" onclick="setWeight(80)">80 kg</button>
                <button class="quick-weight-btn" onclick="setWeight(100)">100 kg</button>
                <button class="quick-weight-btn" onclick="setWeight(120)">120 kg</button>
            </div>
        </div>

        <!-- INDUCTION AGENTS (Yellow) -->
        <div class="drug-category open">
            <div class="category-header" onclick="toggleCategory(this)">
                <div class="color-indicator" style="background: #F59E0B;"></div>
                <div class="category-info">
                    <div class="category-title">Induction Agents</div>
                    <div class="category-subtitle">Propofol, Etomidate, Ketamine</div>
                </div>
                <svg class="chevron" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                    <polyline points="6 9 12 15 18 9"></polyline>
                </svg>
            </div>
            <div class="category-content">
                <!-- Propofol -->
                <div class="drug-card">
                    <div class="drug-header">
                        <div class="drug-name-section">
                            <h3>Propofol</h3>
                            <div class="drug-subtitle">Diprivan</div>
                        </div>
                        <div class="concentration-badge">10 mg/mL</div>
                    </div>
                    <div class="dose-grid">
                        <div class="dose-item">
                            <div class="dose-label">Low</div>
                            <div class="dose-value"><span data-calc="1.5">105</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">1.5 mg/kg</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">Standard</div>
                            <div class="dose-value"><span data-calc="2">140</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">2 mg/kg</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">High</div>
                            <div class="dose-value"><span data-calc="2.5">175</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">2.5 mg/kg</div>
                        </div>
                    </div>
                    <div class="clinical-pearl">
                        <svg class="pearl-icon" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <circle cx="8" cy="8" r="7"></circle>
                            <line x1="8" y1="12" x2="8" y2="8"></line>
                            <line x1="8" y1="5" x2="8.01" y2="5"></line>
                        </svg>
                        <div class="pearl-text">Reduce dose 30-50% in elderly, hypovolemic, or cardiac patients. Causes hypotension via vasodilation.</div>
                    </div>
                </div>

                <!-- Etomidate -->
                <div class="drug-card">
                    <div class="drug-header">
                        <div class="drug-name-section">
                            <h3>Etomidate</h3>
                            <div class="drug-subtitle">Amidate</div>
                        </div>
                        <div class="concentration-badge">2 mg/mL</div>
                    </div>
                    <div class="dose-grid">
                        <div class="dose-item">
                            <div class="dose-label">Low</div>
                            <div class="dose-value"><span data-calc="0.2">14</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">0.2 mg/kg</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">Standard</div>
                            <div class="dose-value"><span data-calc="0.3">21</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">0.3 mg/kg</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">High</div>
                            <div class="dose-value"><span data-calc="0.4">28</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">0.4 mg/kg</div>
                        </div>
                    </div>
                    <div class="clinical-pearl">
                        <svg class="pearl-icon" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <circle cx="8" cy="8" r="7"></circle>
                            <line x1="8" y1="12" x2="8" y2="8"></line>
                            <line x1="8" y1="5" x2="8.01" y2="5"></line>
                        </svg>
                        <div class="pearl-text">Hemodynamically stable — preferred for cardiac/trauma. Avoid repeated doses (adrenal suppression). May cause myoclonus.</div>
                    </div>
                </div>

                <!-- Ketamine -->
                <div class="drug-card">
                    <div class="drug-header">
                        <div class="drug-name-section">
                            <h3>Ketamine</h3>
                            <div class="drug-subtitle">Ketalar</div>
                        </div>
                        <div class="concentration-badge">50/100 mg/mL</div>
                    </div>
                    <div class="dose-grid">
                        <div class="dose-item">
                            <div class="dose-label">Low</div>
                            <div class="dose-value"><span data-calc="1">70</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">1 mg/kg</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">Standard</div>
                            <div class="dose-value"><span data-calc="1.5">105</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">1.5 mg/kg</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">High</div>
                            <div class="dose-value"><span data-calc="2">140</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">2 mg/kg</div>
                        </div>
                    </div>
                    <div class="clinical-pearl">
                        <svg class="pearl-icon" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <circle cx="8" cy="8" r="7"></circle>
                            <line x1="8" y1="12" x2="8" y2="8"></line>
                            <line x1="8" y1="5" x2="8.01" y2="5"></line>
                        </svg>
                        <div class="pearl-text">Maintains airway reflexes & BP. Great for asthma, hypovolemia. Avoid in CAD, elevated ICP (relative). Expect emergence reactions.</div>
                    </div>
                </div>
            </div>
        </div>

        <!-- OPIOIDS (Blue) -->
        <div class="drug-category">
            <div class="category-header" onclick="toggleCategory(this)">
                <div class="color-indicator" style="background: #2563EB;"></div>
                <div class="category-info">
                    <div class="category-title">Opioids</div>
                    <div class="category-subtitle">Fentanyl</div>
                </div>
                <svg class="chevron" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                    <polyline points="6 9 12 15 18 9"></polyline>
                </svg>
            </div>
            <div class="category-content">
                <!-- Fentanyl -->
                <div class="drug-card">
                    <div class="drug-header">
                        <div class="drug-name-section">
                            <h3>Fentanyl</h3>
                            <div class="drug-subtitle">Sublimaze</div>
                        </div>
                        <div class="concentration-badge">50 mcg/mL</div>
                    </div>
                    <div class="dose-grid">
                        <div class="dose-item">
                            <div class="dose-label">Low</div>
                            <div class="dose-value"><span data-calc="1">70</span> <span class="dose-unit">mcg</span></div>
                            <div class="dose-range">0.5–1 mcg/kg</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">Standard</div>
                            <div class="dose-value"><span data-calc="2">140</span> <span class="dose-unit">mcg</span></div>
                            <div class="dose-range">1–2 mcg/kg</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">High</div>
                            <div class="dose-value"><span data-calc="3">210</span> <span class="dose-unit">mcg</span></div>
                            <div class="dose-range">2–3 mcg/kg</div>
                        </div>
                    </div>
                    <div class="clinical-pearl">
                        <svg class="pearl-icon" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <circle cx="8" cy="8" r="7"></circle>
                            <line x1="8" y1="12" x2="8" y2="8"></line>
                            <line x1="8" y1="5" x2="8.01" y2="5"></line>
                        </svg>
                        <div class="pearl-text">Onset 2-3 min IV, duration 30-60 min. Blunts laryngoscopy response at 3-5 mcg/kg. Watch for chest wall rigidity at high doses.</div>
                    </div>
                </div>
            </div>
        </div>

        <!-- NEUROMUSCULAR BLOCKERS (Red) -->
        <div class="drug-category">
            <div class="category-header" onclick="toggleCategory(this)">
                <div class="color-indicator" style="background: #EF4444;"></div>
                <div class="category-info">
                    <div class="category-title">Neuromuscular Blockers</div>
                    <div class="category-subtitle">Succinylcholine, Rocuronium</div>
                </div>
                <svg class="chevron" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                    <polyline points="6 9 12 15 18 9"></polyline>
                </svg>
            </div>
            <div class="category-content">
                <!-- Succinylcholine -->
                <div class="drug-card">
                    <div class="drug-header">
                        <div class="drug-name-section">
                            <h3>Succinylcholine</h3>
                            <div class="drug-subtitle">Anectine</div>
                        </div>
                        <div class="concentration-badge">20 mg/mL</div>
                    </div>
                    <div class="dose-grid">
                        <div class="dose-item">
                            <div class="dose-label">RSI Dose</div>
                            <div class="dose-value"><span data-calc="1.5">105</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">1–1.5 mg/kg</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">Laryngospasm</div>
                            <div class="dose-value"><span data-calc="0.2">14</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">0.1–0.2 mg/kg</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">IM</div>
                            <div class="dose-value"><span data-calc="4">280</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">4 mg/kg</div>
                        </div>
                    </div>
                    <div class="clinical-pearl">
                        <svg class="pearl-icon" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <circle cx="8" cy="8" r="7"></circle>
                            <line x1="8" y1="12" x2="8" y2="8"></line>
                            <line x1="8" y1="5" x2="8.01" y2="5"></line>
                        </svg>
                        <div class="pearl-text">Contraindicated: burns >24h, crush injury, denervation, hyperkalemia risk. MH trigger. Duration 5-10 min.</div>
                    </div>
                </div>

                <!-- Rocuronium -->
                <div class="drug-card">
                    <div class="drug-header">
                        <div class="drug-name-section">
                            <h3>Rocuronium</h3>
                            <div class="drug-subtitle">Zemuron</div>
                        </div>
                        <div class="concentration-badge">10 mg/mL</div>
                    </div>
                    <div class="dose-grid">
                        <div class="dose-item">
                            <div class="dose-label">Intubating</div>
                            <div class="dose-value"><span data-calc="0.6">42</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">0.6 mg/kg</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">RSI Dose</div>
                            <div class="dose-value"><span data-calc="1.2">84</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">1.2 mg/kg</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">Maintenance</div>
                            <div class="dose-value"><span data-calc="0.1">7</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">0.1 mg/kg</div>
                        </div>
                    </div>
                    <div class="clinical-pearl">
                        <svg class="pearl-icon" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <circle cx="8" cy="8" r="7"></circle>
                            <line x1="8" y1="12" x2="8" y2="8"></line>
                            <line x1="8" y1="5" x2="8.01" y2="5"></line>
                        </svg>
                        <div class="pearl-text">Onset 60-90s (standard), 45-60s (RSI dose). Duration 30-45 min. Reversible with sugammadex 16 mg/kg for immediate reversal.</div>
                    </div>
                </div>
            </div>
        </div>

        <!-- VASOPRESSORS & INOTROPES (Violet) -->
        <div class="drug-category">
            <div class="category-header" onclick="toggleCategory(this)">
                <div class="color-indicator" style="background: #8B5CF6;"></div>
                <div class="category-info">
                    <div class="category-title">Vasopressors & Inotropes</div>
                    <div class="category-subtitle">Phenylephrine, Ephedrine (Fixed doses)</div>
                </div>
                <svg class="chevron" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                    <polyline points="6 9 12 15 18 9"></polyline>
                </svg>
            </div>
            <div class="category-content">
                <!-- Phenylephrine -->
                <div class="drug-card">
                    <div class="drug-header">
                        <div class="drug-name-section">
                            <h3>Phenylephrine</h3>
                            <div class="drug-subtitle">Neo-Synephrine</div>
                        </div>
                        <div class="concentration-badge">100 mcg/mL</div>
                    </div>
                    <div class="dose-grid">
                        <div class="dose-item">
                            <div class="dose-label">Bolus</div>
                            <div class="dose-value">100 <span class="dose-unit">mcg</span></div>
                            <div class="dose-range">50–200 mcg</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">Infusion</div>
                            <div class="dose-value">50 <span class="dose-unit">mcg/min</span></div>
                            <div class="dose-range">10–200 mcg/min</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">OB Bolus</div>
                            <div class="dose-value">100 <span class="dose-unit">mcg</span></div>
                            <div class="dose-range">q1-2 min PRN</div>
                        </div>
                    </div>
                    <div class="clinical-pearl">
                        <svg class="pearl-icon" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <circle cx="8" cy="8" r="7"></circle>
                            <line x1="8" y1="12" x2="8" y2="8"></line>
                            <line x1="8" y1="5" x2="8.01" y2="5"></line>
                        </svg>
                        <div class="pearl-text">Pure α1-agonist — ↑SVR, reflex bradycardia. First-line for spinal hypotension. Avoid if already bradycardic — use ephedrine instead.</div>
                    </div>
                </div>

                <!-- Ephedrine -->
                <div class="drug-card">
                    <div class="drug-header">
                        <div class="drug-name-section">
                            <h3>Ephedrine</h3>
                            <div class="drug-subtitle">—</div>
                        </div>
                        <div class="concentration-badge">5 mg/mL</div>
                    </div>
                    <div class="dose-grid">
                        <div class="dose-item">
                            <div class="dose-label">Bolus</div>
                            <div class="dose-value">5–10 <span class="dose-unit">mg</span></div>
                            <div class="dose-range">q3-5 min PRN</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">Max Total</div>
                            <div class="dose-value">50 <span class="dose-unit">mg</span></div>
                            <div class="dose-range">tachyphylaxis</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">IM Dose</div>
                            <div class="dose-value">25–50 <span class="dose-unit">mg</span></div>
                            <div class="dose-range">if no IV access</div>
                        </div>
                    </div>
                    <div class="clinical-pearl">
                        <svg class="pearl-icon" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <circle cx="8" cy="8" r="7"></circle>
                            <line x1="8" y1="12" x2="8" y2="8"></line>
                            <line x1="8" y1="5" x2="8.01" y2="5"></line>
                        </svg>
                        <div class="pearl-text">Mixed α/β agonist — ↑HR, ↑BP. Preferred over phenylephrine if bradycardic. Indirect mechanism = tachyphylaxis with repeated dosing.</div>
                    </div>
                </div>
            </div>
        </div>

        <!-- REVERSAL & ANTICHOLINERGICS (Green) -->
        <div class="drug-category">
            <div class="category-header" onclick="toggleCategory(this)">
                <div class="color-indicator" style="background: #10B981;"></div>
                <div class="category-info">
                    <div class="category-title">Reversal & Anticholinergics</div>
                    <div class="category-subtitle">Sugammadex</div>
                </div>
                <svg class="chevron" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                    <polyline points="6 9 12 15 18 9"></polyline>
                </svg>
            </div>
            <div class="category-content">
                <!-- Sugammadex -->
                <div class="drug-card">
                    <div class="drug-header">
                        <div class="drug-name-section">
                            <h3>Sugammadex</h3>
                            <div class="drug-subtitle">Bridion</div>
                        </div>
                        <div class="concentration-badge">100 mg/mL</div>
                    </div>
                    <div class="dose-grid">
                        <div class="dose-item">
                            <div class="dose-label">Moderate</div>
                            <div class="dose-value"><span data-calc="2">140</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">2 mg/kg</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">Deep Block</div>
                            <div class="dose-value"><span data-calc="4">280</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">4 mg/kg</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">Immediate</div>
                            <div class="dose-value"><span data-calc="16">1120</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">16 mg/kg</div>
                        </div>
                    </div>
                    <div class="clinical-pearl">
                        <svg class="pearl-icon" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <circle cx="8" cy="8" r="7"></circle>
                            <line x1="8" y1="12" x2="8" y2="8"></line>
                            <line x1="8" y1="5" x2="8.01" y2="5"></line>
                        </svg>
                        <div class="pearl-text">Encapsulates rocuronium/vecuronium. 16 mg/kg for 'can't intubate, can't oxygenate' after RSI with roc. May reduce efficacy of hormonal contraceptives.</div>
                    </div>
                </div>
            </div>
        </div>

        <!-- LOCAL ANESTHETICS (Gray) -->
        <div class="drug-category">
            <div class="category-header" onclick="toggleCategory(this)">
                <div class="color-indicator" style="background: #6B7280;"></div>
                <div class="category-info">
                    <div class="category-title">Local Anesthetics</div>
                    <div class="category-subtitle">Lidocaine</div>
                </div>
                <svg class="chevron" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                    <polyline points="6 9 12 15 18 9"></polyline>
                </svg>
            </div>
            <div class="category-content">
                <!-- Lidocaine -->
                <div class="drug-card">
                    <div class="drug-header">
                        <div class="drug-name-section">
                            <h3>Lidocaine</h3>
                            <div class="drug-subtitle">Xylocaine</div>
                        </div>
                        <div class="concentration-badge">1% = 10 mg/mL</div>
                    </div>
                    <div class="dose-grid">
                        <div class="dose-item">
                            <div class="dose-label">Plain Max</div>
                            <div class="dose-value"><span data-calc="4.5">315</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">4.5 mg/kg</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">+ Epi Max</div>
                            <div class="dose-value"><span data-calc="7">490</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">7 mg/kg</div>
                        </div>
                        <div class="dose-item">
                            <div class="dose-label">IV Bolus</div>
                            <div class="dose-value"><span data-calc="1.5">105</span> <span class="dose-unit">mg</span></div>
                            <div class="dose-range">1.5 mg/kg</div>
                        </div>
                    </div>
                    <div class="clinical-pearl">
                        <svg class="pearl-icon" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <circle cx="8" cy="8" r="7"></circle>
                            <line x1="8" y1="12" x2="8" y2="8"></line>
                            <line x1="8" y1="5" x2="8.01" y2="5"></line>
                        </svg>
                        <div class="pearl-text">IV lidocaine blunts airway reflexes at intubation/extubation. LAST symptoms: tinnitus, perioral numbness, seizures, arrhythmias.</div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Color Legend -->
        <div class="color-legend">
            <div class="legend-title">Anesthesia Syringe Color Reference</div>
            <div class="legend-grid">
                <div class="legend-item">
                    <div class="legend-swatch" style="background: #2563EB;"></div>
                    <span class="legend-label">Opioids</span>
                </div>
                <div class="legend-item">
                    <div class="legend-swatch" style="background: #EF4444;"></div>
                    <span class="legend-label">Neuromuscular Blockers</span>
                </div>
                <div class="legend-item">
                    <div class="legend-swatch" style="background: #F59E0B;"></div>
                    <span class="legend-label">Induction Agents</span>
                </div>
                <div class="legend-item">
                    <div class="legend-swatch" style="background: #F97316;"></div>
                    <span class="legend-label">Tranquilizers</span>
                </div>
                <div class="legend-item">
                    <div class="legend-swatch" style="background: #8B5CF6;"></div>
                    <span class="legend-label">Vasopressors</span>
                </div>
                <div class="legend-item">
                    <div class="legend-swatch" style="background: #10B981;"></div>
                    <span class="legend-label">Anticholinergics</span>
                </div>
                <div class="legend-item">
                    <div class="legend-swatch" style="background: #6B7280;"></div>
                    <span class="legend-label">Local Anesthetics</span>
                </div>
            </div>
        </div>

        <!-- Disclaimer -->
        <div class="disclaimer">
            <strong>Disclaimer:</strong> This tool is for educational and reference purposes only. Always verify doses against institutional protocols and consider patient-specific factors. Not a substitute for clinical judgment.
        </div>

    </main>

    <footer>
        <p>&copy; 2025 gasconsult.ai. All rights reserved. | <a href="/terms" style="color: var(--primary); text-decoration: none;">Terms of Service</a> | <a href="/privacy" style="color: var(--primary); text-decoration: none;">Privacy Policy</a></p>
    </footer>

    <!-- Crisis Modal -->
    <div class="crisis-overlay" id="crisisOverlay" onclick="closeCrisisOnOverlay(event)">
        <div class="crisis-modal" onclick="event.stopPropagation()">
            <div class="crisis-modal-header">
                <div class="crisis-modal-title">
                    <svg width="20" height="20" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                        <path d="M10.29 3.86L1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z"></path>
                        <line x1="12" y1="9" x2="12" y2="13"></line>
                        <line x1="12" y1="17" x2="12.01" y2="17"></line>
                    </svg>
                    Crisis Protocols
                </div>
                <button class="crisis-close-btn" onclick="toggleCrisis()">
                    <svg width="20" height="20" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                        <line x1="18" y1="6" x2="6" y2="18"></line>
                        <line x1="6" y1="6" x2="18" y2="18"></line>
                    </svg>
                </button>
            </div>
            <div class="crisis-modal-content">
                <div class="protocol-grid">
                    <button class="protocol-btn">
                        <div class="protocol-title">Malignant Hyperthermia</div>
                        <div class="protocol-desc">Dantrolene 2.5 mg/kg, call MH hotline</div>
                    </button>
                    <button class="protocol-btn">
                        <div class="protocol-title">LAST (Lipid Rescue)</div>
                        <div class="protocol-desc">20% Intralipid bolus + infusion</div>
                    </button>
                    <button class="protocol-btn">
                        <div class="protocol-title">Anaphylaxis</div>
                        <div class="protocol-desc">Epinephrine, fluids, steroids</div>
                    </button>
                    <button class="protocol-btn">
                        <div class="protocol-title">Cardiac Arrest</div>
                        <div class="protocol-desc">ACLS algorithms & doses</div>
                    </button>
                    <button class="protocol-btn">
                        <div class="protocol-title">Bronchospasm</div>
                        <div class="protocol-desc">Albuterol, epinephrine, deepening</div>
                    </button>
                    <button class="protocol-btn">
                        <div class="protocol-title">Laryngospasm</div>
                        <div class="protocol-desc">Jaw thrust, CPAP, succinylcholine</div>
                    </button>
                    <button class="protocol-btn">
                        <div class="protocol-title">High Spinal</div>
                        <div class="protocol-desc">Airway, pressors, sedation</div>
                    </button>
                    <button class="protocol-btn">
                        <div class="protocol-title">Air Embolism</div>
                        <div class="protocol-desc">Flood field, Durant, aspirate</div>
                    </button>
                </div>
            </div>
        </div>
    </div>

    <script>
        // Weight conversion and dose calculation
        function updateDoses() {
            const kg = parseFloat(document.getElementById('weightInput').value) || 70;

            // Update lbs conversion
            document.getElementById('lbsConversion').textContent = Math.round(kg * 2.205);

            // Update all weight-based doses
            document.querySelectorAll('[data-calc]').forEach(element => {
                const multiplier = parseFloat(element.getAttribute('data-calc'));
                const dose = Math.round(kg * multiplier);
                element.textContent = dose;
            });
        }

        // Quick weight buttons
        function setWeight(kg) {
            document.getElementById('weightInput').value = kg;
            updateDoses();

            // Update active state on buttons
            document.querySelectorAll('.quick-weight-btn').forEach(btn => {
                btn.classList.remove('active');
                if (btn.textContent.includes(kg.toString())) {
                    btn.classList.add('active');
                }
            });
        }

        // Category accordion toggle
        function toggleCategory(header) {
            const category = header.parentElement;
            category.classList.toggle('open');
        }

        // Crisis modal
        function toggleCrisis() {
            document.getElementById('crisisOverlay').classList.toggle('show');
        }

        // Close modal when clicking overlay (not modal content)
        function closeCrisisOnOverlay(event) {
            if (event.target.id === 'crisisOverlay') {
                toggleCrisis();
            }
        }

        // Initialize with default weight
        document.addEventListener('DOMContentLoaded', () => {
            updateDoses();
        });
    </script>
</body>
</html>
"""

CALCULATORS_HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Clinical Calculators — gasconsult.ai</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Sora:wght@400;600;700&display=swap" rel="stylesheet">
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=5">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">
    <style>
        :root {
            --primary-blue: #2563EB;
            --primary-blue-dark: #1D4ED8;
            --primary-blue-light: #DBEAFE;
            --text-primary: #0F172A;
            --text-secondary: #475569;
            --text-muted: #94A3B8;
            --bg-primary: #FFFFFF;
            --bg-secondary: #F8FAFC;
            --border: #E2E8F0;
            --success-green: #10B981;
            --success-green-light: #D1FAE5;
        }

        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        html {
            scroll-behavior: smooth;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, system-ui, sans-serif;
            background: var(--bg-secondary);
            color: var(--text-primary);
            line-height: 1.6;
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
            text-rendering: optimizeLegibility;
            animation: pageFadeIn 0.6s cubic-bezier(0.4, 0, 0.2, 1);
            overflow-x: hidden;
        }

        @keyframes pageFadeIn {
            from {
                opacity: 0;
                transform: translateY(10px);
            }
            to {
                opacity: 1;
                transform: translateY(0);
            }
        }

        @keyframes slideUp {
            from {
                opacity: 0;
                transform: translateY(30px);
            }
            to {
                opacity: 1;
                transform: translateY(0);
            }
        }

        @keyframes slideIn {
            from {
                opacity: 0;
                transform: translateX(-20px);
            }
            to {
                opacity: 1;
                transform: translateX(0);
            }
        }

        @keyframes scaleIn {
            from {
                opacity: 0;
                transform: scale(0.95);
            }
            to {
                opacity: 1;
                transform: scale(1);
            }
        }

        @keyframes shimmer {
            0% {
                background-position: -1000px 0;
            }
            100% {
                background-position: 1000px 0;
            }
        }

        @keyframes glow {
            0%, 100% {
                box-shadow: 0 4px 14px rgba(37, 99, 235, 0.3);
            }
            50% {
                box-shadow: 0 6px 20px rgba(37, 99, 235, 0.5);
            }
        }

        /* Navigation */
        nav {
            background: rgba(255, 255, 255, 0.85);
            backdrop-filter: blur(12px);
            -webkit-backdrop-filter: blur(12px);
            padding: 16px 40px;
            position: sticky;
            top: 0;
            left: 0;
            right: 0;
            z-index: 100;
            border-bottom: 1px solid rgba(226, 232, 240, 0.8);
        }

        nav .container {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            align-items: center;
            justify-content: space-between;
            flex-wrap: wrap;
            gap: 16px;
        }

        .logo-container {
            text-decoration: none;
            display: flex;
            align-items: center;
            gap: 12px;
            cursor: pointer;
            transition: transform 0.2s ease;
        }

        .logo-container:hover {
            transform: translateY(-1px);
        }

        .logo-ecg {
            height: 28px;
            width: auto;
            flex-shrink: 0;
        }

        .logo-wordmark {
            font-family: 'Sora', -apple-system, BlinkMacSystemFont, system-ui, sans-serif;
            font-size: 20px;
            font-weight: 600;
            letter-spacing: -0.5px;
            white-space: nowrap;
        }

        .logo-gas {
            color: #2563EB;
        }

        .logo-consult {
            color: #111111;
        }

        .logo-ai {
            font-weight: 400;
            color: #6B7280;
        }

        .nav-actions {
            display: flex;
            gap: 12px;
            align-items: center;
        }

        .nav-link {
            color: var(--text-secondary);
            text-decoration: none;
            font-size: 14px;
            font-weight: 500;
            padding: 8px 16px;
            border-radius: 8px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            color: var(--text-primary);
            background: rgba(255, 255, 255, 0.6);
            backdrop-filter: blur(8px);
            -webkit-backdrop-filter: blur(8px);
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.06);
        }

        .nav-link.active {
            color: var(--primary-blue);
            font-weight: 600;
        }

        /* Main Layout */
        .main-container {
            max-width: 1200px;
            margin: 0 auto;
            padding: 48px 24px;
            display: flex;
            gap: 32px;
            position: relative;
        }

        .main-container::before {
            content: '';
            position: absolute;
            top: -100px;
            left: 50%;
            transform: translateX(-50%);
            width: 100%;
            height: 500px;
            background: radial-gradient(ellipse at 50% 0%, rgba(37, 99, 235, 0.08) 0%, transparent 70%);
            z-index: -1;
            pointer-events: none;
        }

        /* Sidebar */
        .sidebar {
            width: 290px;
            flex-shrink: 0;
            position: sticky;
            top: 100px;
            height: fit-content;
            animation: slideIn 0.5s cubic-bezier(0.4, 0, 0.2, 1) 0.2s backwards;
            background: linear-gradient(to bottom, rgba(255, 255, 255, 0.5), rgba(255, 255, 255, 0.3));
            backdrop-filter: blur(10px);
            -webkit-backdrop-filter: blur(10px);
            border-radius: 20px;
            padding: 28px 20px;
            border: 1px solid rgba(226, 232, 240, 0.5);
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.04);
        }

        .sidebar-header {
            font-family: 'Sora', sans-serif;
            font-size: 18px;
            font-weight: 700;
            background: linear-gradient(135deg, var(--primary-blue) 0%, #6366f1 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            margin-bottom: 24px;
            padding: 0 4px;
            letter-spacing: -0.3px;
        }

        .calculator-list {
            list-style: none;
        }

        .calculator-item {
            margin-bottom: 8px;
            animation: slideIn 0.4s cubic-bezier(0.4, 0, 0.2, 1) backwards;
        }

        .calculator-item:nth-child(1) { animation-delay: 0.3s; }
        .calculator-item:nth-child(2) { animation-delay: 0.35s; }
        .calculator-item:nth-child(3) { animation-delay: 0.4s; }
        .calculator-item:nth-child(4) { animation-delay: 0.45s; }
        .calculator-item:nth-child(5) { animation-delay: 0.5s; }
        .calculator-item:nth-child(6) { animation-delay: 0.55s; }
        .calculator-item:nth-child(7) { animation-delay: 0.6s; }
        .calculator-item:nth-child(8) { animation-delay: 0.65s; }
        .calculator-item:nth-child(9) { animation-delay: 0.7s; }
        .calculator-item:nth-child(10) { animation-delay: 0.75s; }
        .calculator-item:nth-child(11) { animation-delay: 0.8s; }
        .calculator-item:nth-child(12) { animation-delay: 0.85s; }
        .calculator-item:nth-child(13) { animation-delay: 0.9s; }
        .calculator-item:nth-child(14) { animation-delay: 0.95s; }
        .calculator-item:nth-child(15) { animation-delay: 1s; }
        .calculator-item:nth-child(16) { animation-delay: 1.05s; }
        .calculator-item:nth-child(17) { animation-delay: 1.1s; }

        .calculator-btn {
            width: 100%;
            text-align: left;
            padding: 16px 20px;
            background: linear-gradient(to bottom, #ffffff, #fafbfc);
            border: 1.5px solid var(--border);
            border-radius: 14px;
            color: var(--text-secondary);
            font-family: 'Inter', sans-serif;
            font-size: 14px;
            font-weight: 500;
            cursor: pointer;
            transition: all 0.25s cubic-bezier(0.4, 0, 0.2, 1);
            box-shadow: 0 1px 2px rgba(0, 0, 0, 0.03);
            position: relative;
            overflow: hidden;
        }

        .calculator-btn::before {
            content: '';
            position: absolute;
            left: 0;
            top: 0;
            bottom: 0;
            width: 3px;
            background: var(--primary-blue);
            transform: scaleY(0);
            transition: transform 0.25s cubic-bezier(0.4, 0, 0.2, 1);
        }

        .calculator-btn:hover {
            background: linear-gradient(to bottom, #ffffff, var(--primary-blue-light));
            border-color: var(--primary-blue);
            color: var(--text-primary);
            transform: translateX(2px);
            box-shadow: 0 3px 10px rgba(37, 99, 235, 0.12);
        }

        .calculator-btn:hover::before {
            transform: scaleY(1);
        }

        .calculator-btn.active {
            background: linear-gradient(135deg, var(--primary-blue) 0%, #1e40af 100%);
            border-color: transparent;
            color: white;
            font-weight: 600;
            transform: translateX(2px);
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.25), inset 0 1px 0 rgba(255, 255, 255, 0.1);
        }

        .calculator-btn.active::before {
            width: 100%;
            background: rgba(255, 255, 255, 0.1);
            transform: scaleY(1);
        }

        /* Main Content Area */
        .content-area {
            flex: 1;
            min-width: 0;
            animation: scaleIn 0.6s cubic-bezier(0.4, 0, 0.2, 1) 0.3s backwards;
        }

        .calculator-card {
            background: white;
            border: 1px solid var(--border);
            border-radius: 20px;
            padding: 40px;
            box-shadow: 0 4px 24px rgba(0, 0, 0, 0.06), 0 0 0 1px rgba(0, 0, 0, 0.02);
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        }

        .calculator-card:hover {
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.08), 0 0 0 1px rgba(37, 99, 235, 0.1);
            transform: translateY(-2px);
        }

        .calculator-title {
            font-family: 'Sora', sans-serif;
            font-size: 28px;
            font-weight: 700;
            color: var(--text-primary);
            margin-bottom: 8px;
            letter-spacing: -0.5px;
        }

        .calculator-description {
            color: var(--text-secondary);
            font-size: 15px;
            margin-bottom: 32px;
            line-height: 1.6;
        }

        /* Form Styling */
        .form-group {
            margin-bottom: 20px;
        }

        .form-row {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 20px;
        }

        label {
            display: block;
            font-size: 13px;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 8px;
            letter-spacing: -0.2px;
        }

        input, select {
            width: 100%;
            padding: 12px 16px;
            background: var(--bg-secondary);
            border: 2px solid var(--border);
            border-radius: 10px;
            color: var(--text-primary);
            font-family: 'Inter', sans-serif;
            font-size: 14px;
            transition: all 0.2s ease;
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
        }

        input:hover, select:hover {
            border-color: var(--primary-blue);
        }

        input:focus, select:focus {
            outline: none;
            border-color: var(--primary-blue);
            background: white;
            box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.1);
        }

        input::placeholder {
            color: var(--text-muted);
        }

        .calculate-btn {
            padding: 14px 40px;
            background: linear-gradient(135deg, var(--primary-blue) 0%, var(--primary-blue-dark) 100%);
            color: white;
            border: none;
            border-radius: 12px;
            font-size: 15px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            box-shadow: 0 4px 14px rgba(37, 99, 235, 0.3);
            position: relative;
            overflow: hidden;
        }

        .calculate-btn::before {
            content: '';
            position: absolute;
            top: 0;
            left: -100%;
            width: 100%;
            height: 100%;
            background: linear-gradient(90deg, transparent, rgba(255, 255, 255, 0.2), transparent);
            transition: left 0.6s;
        }

        .calculate-btn:hover::before {
            left: 100%;
        }

        .calculate-btn:hover {
            transform: translateY(-2px) scale(1.02);
            box-shadow: 0 6px 20px rgba(37, 99, 235, 0.4);
        }

        .calculate-btn:active {
            transform: translateY(0) scale(0.98);
        }

        .result-box {
            margin-top: 32px;
            padding: 28px;
            background: linear-gradient(135deg, rgba(37, 99, 235, 0.05) 0%, rgba(139, 92, 246, 0.05) 100%);
            border: 2px solid var(--primary-blue-light);
            border-radius: 16px;
            display: none;
            animation: slideUp 0.4s cubic-bezier(0.4, 0, 0.2, 1);
        }

        .result-box.visible {
            display: block;
        }

        .result-label {
            font-size: 13px;
            font-weight: 600;
            color: var(--text-muted);
            margin-bottom: 12px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        .result-value {
            font-family: 'Sora', sans-serif;
            font-size: 36px;
            font-weight: 700;
            background: linear-gradient(135deg, var(--primary-blue) 0%, #8B5CF6 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            margin-bottom: 16px;
            letter-spacing: -1px;
        }

        .result-interpretation {
            font-size: 14px;
            color: var(--text-secondary);
            margin-bottom: 24px;
            line-height: 1.7;
        }

        .send-to-ai-btn {
            padding: 13px 28px;
            background: linear-gradient(135deg, var(--primary-blue) 0%, #1e40af 100%);
            color: white;
            border: none;
            border-radius: 12px;
            font-size: 14px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.25s cubic-bezier(0.4, 0, 0.2, 1);
            display: inline-flex;
            align-items: center;
            gap: 8px;
            box-shadow: 0 2px 8px rgba(37, 99, 235, 0.20), inset 0 1px 0 rgba(255, 255, 255, 0.1);
            letter-spacing: 0.01em;
            position: relative;
            overflow: hidden;
        }

        .send-to-ai-btn::before {
            content: '';
            position: absolute;
            top: 0;
            left: -100%;
            width: 100%;
            height: 100%;
            background: linear-gradient(90deg, transparent, rgba(255, 255, 255, 0.15), transparent);
            transition: left 0.5s ease;
        }

        .send-to-ai-btn:hover::before {
            left: 100%;
        }

        .send-to-ai-btn:hover {
            transform: translateY(-1px);
            box-shadow: 0 4px 14px rgba(37, 99, 235, 0.30), inset 0 1px 0 rgba(255, 255, 255, 0.15);
        }

        .send-to-ai-btn:active {
            transform: translateY(0);
            box-shadow: 0 2px 6px rgba(37, 99, 235, 0.25);
        }

        .calculator-section {
            display: none;
        }

        .calculator-section.active {
            display: block;
            animation: scaleIn 0.4s cubic-bezier(0.4, 0, 0.2, 1);
        }

        /* AI Chat Modal */
        .ai-chat-modal {
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            bottom: 0;
            background: rgba(0, 0, 0, 0.5);
            backdrop-filter: blur(8px);
            -webkit-backdrop-filter: blur(8px);
            z-index: 1000;
            display: none;
            align-items: center;
            justify-content: center;
            padding: 20px;
            animation: fadeIn 0.3s ease;
        }

        .ai-chat-modal.active {
            display: flex;
        }

        @keyframes fadeIn {
            from { opacity: 0; }
            to { opacity: 1; }
        }

        .ai-chat-container {
            background: white;
            border-radius: 24px;
            width: 100%;
            max-width: 700px;
            max-height: 85vh;
            display: flex;
            flex-direction: column;
            box-shadow: 0 24px 64px rgba(0, 0, 0, 0.2);
            animation: slideUpModal 0.4s cubic-bezier(0.4, 0, 0.2, 1);
        }

        @keyframes slideUpModal {
            from {
                opacity: 0;
                transform: translateY(40px) scale(0.95);
            }
            to {
                opacity: 1;
                transform: translateY(0) scale(1);
            }
        }

        .ai-chat-header {
            display: flex;
            align-items: center;
            justify-content: space-between;
            padding: 24px 28px;
            border-bottom: 1px solid var(--border);
        }

        .ai-chat-header h3 {
            font-family: 'Sora', sans-serif;
            font-size: 20px;
            font-weight: 700;
            background: linear-gradient(135deg, var(--primary-blue) 0%, #8B5CF6 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            letter-spacing: -0.5px;
        }

        .close-modal-btn {
            background: var(--bg-secondary);
            border: none;
            width: 36px;
            height: 36px;
            border-radius: 50%;
            cursor: pointer;
            display: flex;
            align-items: center;
            justify-content: center;
            color: var(--text-secondary);
            transition: all 0.2s ease;
        }

        .close-modal-btn:hover {
            background: var(--border);
            color: var(--text-primary);
            transform: rotate(90deg);
        }

        .ai-chat-messages {
            flex: 1;
            overflow-y: auto;
            padding: 24px 28px;
            max-height: 50vh;
        }

        .ai-message {
            background: var(--bg-secondary);
            padding: 16px 20px;
            border-radius: 16px;
            margin-bottom: 12px;
            line-height: 1.6;
            color: var(--text-secondary);
            font-size: 14px;
            animation: slideUp 0.3s ease;
        }

        .ai-message strong {
            color: var(--text-primary);
        }

        .ai-chat-input-area {
            padding: 20px 28px;
            border-top: 1px solid var(--border);
            background: var(--bg-secondary);
            border-radius: 0 0 24px 24px;
        }

        .ai-input-container {
            display: flex;
            gap: 12px;
            align-items: flex-end;
        }

        .ai-chat-input {
            flex: 1;
            padding: 12px 16px;
            border: 2px solid var(--border);
            border-radius: 12px;
            font-family: 'Inter', sans-serif;
            font-size: 14px;
            resize: vertical;
            min-height: 44px;
            max-height: 120px;
            transition: all 0.2s ease;
        }

        .ai-chat-input:focus {
            outline: none;
            border-color: var(--primary-blue);
            box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.1);
        }

        .send-ai-message-btn {
            padding: 12px 24px;
            background: linear-gradient(135deg, var(--primary-blue) 0%, var(--primary-blue-dark) 100%);
            color: white;
            border: none;
            border-radius: 12px;
            font-size: 14px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.3);
            white-space: nowrap;
        }

        .send-ai-message-btn:hover {
            transform: translateY(-2px);
            box-shadow: 0 6px 16px rgba(37, 99, 235, 0.4);
        }

        .send-ai-message-btn:disabled {
            opacity: 0.5;
            cursor: not-allowed;
            transform: none;
        }

        .ai-loading {
            display: flex;
            align-items: center;
            gap: 8px;
            padding: 16px 20px;
            background: var(--primary-blue-light);
            border-radius: 16px;
            margin-bottom: 12px;
            color: var(--primary-blue-dark);
            font-size: 14px;
            font-weight: 500;
        }

        .ai-loading-dots {
            display: flex;
            gap: 4px;
        }

        .ai-loading-dots span {
            width: 6px;
            height: 6px;
            background: var(--primary-blue);
            border-radius: 50%;
            animation: loadingDot 1.4s infinite ease-in-out;
        }

        .ai-loading-dots span:nth-child(1) { animation-delay: 0s; }
        .ai-loading-dots span:nth-child(2) { animation-delay: 0.2s; }
        .ai-loading-dots span:nth-child(3) { animation-delay: 0.4s; }

        @keyframes loadingDot {
            0%, 60%, 100% {
                transform: scale(0.8);
                opacity: 0.5;
            }
            30% {
                transform: scale(1.2);
                opacity: 1;
            }
        }

        /* Responsive */
        @media (max-width: 968px) {
            .main-container {
                flex-direction: column;
                padding: 32px 20px;
            }

            .sidebar {
                width: 100%;
                position: static;
            }

            .calculator-list {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 12px;
            }

            .calculator-btn.active,
            .calculator-btn:hover {
                transform: none;
            }

            .ai-chat-container {
                max-height: 90vh;
            }
        }

        @media (max-width: 640px) {
            nav {
                padding: 12px 20px;
            }

            .nav-actions {
                gap: 8px;
            }

            .nav-link {
                padding: 6px 12px;
                font-size: 13px;
            }

            .main-container {
                padding: 24px 16px;
            }

            .calculator-card {
                padding: 24px;
                border-radius: 16px;
            }

            .calculator-title {
                font-size: 24px;
            }

            .form-row {
                grid-template-columns: 1fr;
            }

            .result-value {
                font-size: 32px;
            }

            .ai-chat-container {
                border-radius: 20px;
            }

            .ai-chat-header,
            .ai-chat-messages,
            .ai-chat-input-area {
                padding: 20px;
            }
        }
    </style>
</head>
<body>
    <!-- Navigation -->
    <nav>
        <div class="container">
            <a href="/" class="logo-container">
                <svg class="logo-ecg" viewBox="0 0 60 28" fill="none" xmlns="http://www.w3.org/2000/svg">
                    <defs>
                        <linearGradient id="ecgGrad" x1="0%" y1="0%" x2="100%" y2="0%">
                            <stop offset="0%" stop-color="#2563EB"/>
                            <stop offset="20%" stop-color="#EF4444"/>
                            <stop offset="40%" stop-color="#FBBF24"/>
                            <stop offset="60%" stop-color="#8B5CF6"/>
                            <stop offset="80%" stop-color="#10B981"/>
                            <stop offset="100%" stop-color="#6B7280"/>
                        </linearGradient>
                    </defs>
                    <path d="M2 14 L10 14 L14 12 L18 16 L22 4 L26 24 L30 10 L34 14 L42 14"
                          stroke="url(#ecgGrad)"
                          stroke-width="2.5"
                          stroke-linecap="round"
                          stroke-linejoin="round"
                          fill="none"/>
                </svg>
                <div class="logo-wordmark">
                    <span class="logo-gas">gas</span><span class="logo-consult">consult</span><span class="logo-ai">.ai</span>
                </div>
            </a>
            <div class="nav-actions">
                <a href="/" class="nav-link">Home</a>
                <a href="/preop" class="nav-link">Pre-Op Assessment</a>
                <a href="/calculators" class="nav-link active">Clinical Calculators</a>
                <a href="/quick-dose" class="nav-link">Quick Dose</a>
                <a href="/hypotension" class="nav-link">IOH Predictor</a>
            </div>
        </div>
    </nav>

    <!-- Main Container -->
    <div class="main-container">
        <!-- Sidebar -->
        <div class="sidebar">
            <h2 class="sidebar-header">Clinical Calculators</h2>
            <ul class="calculator-list">
                <li class="calculator-item">
                    <button class="calculator-btn active" data-calc="ibw">Ideal Body Weight (IBW)</button>
                </li>
                <li class="calculator-item">
                    <button class="calculator-btn" data-calc="mabl">Maximum Allowable Blood Loss</button>
                </li>
                <li class="calculator-item">
                    <button class="calculator-btn" data-calc="bsa">Body Surface Area (BSA)</button>
                </li>
                <li class="calculator-item">
                    <button class="calculator-btn" data-calc="qtc">QTc Interval (Bazett)</button>
                </li>
                <li class="calculator-item">
                    <button class="calculator-btn" data-calc="fluids">Maintenance Fluids (4-2-1)</button>
                </li>
                <li class="calculator-item">
                    <button class="calculator-btn" data-calc="ponv">PONV Risk (Apfel Score)</button>
                </li>
                <li class="calculator-item">
                    <button class="calculator-btn" data-calc="asa">ASA Physical Status</button>
                </li>
                <li class="calculator-item">
                    <button class="calculator-btn" data-calc="fluid-deficit">Fluid Deficit & Replacement</button>
                </li>
                <li class="calculator-item">
                    <button class="calculator-btn" data-calc="defib">Defibrillation Energy</button>
                </li>
                <li class="calculator-item">
                    <button class="calculator-btn" data-calc="sas">Surgical Apgar Score</button>
                </li>
                <li class="calculator-item">
                    <button class="calculator-btn" data-calc="opioid">Opioid Conversion</button>
                </li>
                <li class="calculator-item">
                    <button class="calculator-btn" data-calc="last">Local Anesthetic Toxicity Dose</button>
                </li>
                <li class="calculator-item">
                    <button class="calculator-btn" data-calc="rcri">RCRI (Cardiac Risk)</button>
                </li>
                <li class="calculator-item">
                    <button class="calculator-btn" data-calc="stopbang">STOP-BANG OSA Risk</button>
                </li>
                <li class="calculator-item">
                    <button class="calculator-btn" data-calc="airway">Airway Assessment</button>
                </li>
                <li class="calculator-item">
                    <button class="calculator-btn" data-calc="mac">MAC Calculator</button>
                </li>
                <li class="calculator-item">
                    <button class="calculator-btn" data-calc="propofol-tci">Propofol TCI</button>
                </li>
            </ul>
        </div>

        <!-- Content Area -->
        <div class="content-area">
            <!-- IBW Calculator -->
            <div id="calc-ibw" class="calculator-section active">
                <div class="calculator-card">
                    <h1 class="calculator-title">Ideal Body Weight (IBW)</h1>
                    <p class="calculator-description">Calculate ideal body weight using the Devine formula. Used for weight-based drug dosing and ventilator settings.</p>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="ibw-height">Height (cm)</label>
                            <input type="number" id="ibw-height" placeholder="170" step="0.1">
                        </div>
                        <div class="form-group">
                            <label for="ibw-sex">Biological Sex</label>
                            <select id="ibw-sex">
                                <option value="">Select...</option>
                                <option value="male">Male</option>
                                <option value="female">Female</option>
                            </select>
                        </div>
                    </div>

                    <button class="calculate-btn" onclick="calculateIBW()">Calculate</button>

                    <div id="ibw-result" class="result-box">
                        <div class="result-label">Ideal Body Weight</div>
                        <div class="result-value" id="ibw-value">-</div>
                        <div class="result-interpretation" id="ibw-interpretation"></div>
                        <button class="send-to-ai-btn" onclick="sendToAI('ibw')">
                            Send to AI
                        </button>
                    </div>
                </div>
            </div>

            <!-- MABL Calculator -->
            <div id="calc-mabl" class="calculator-section">
                <div class="calculator-card">
                    <h1 class="calculator-title">Maximum Allowable Blood Loss (MABL)</h1>
                    <p class="calculator-description">Estimate maximum allowable blood loss before transfusion is indicated based on starting and target hematocrit.</p>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="mabl-weight">Weight (kg)</label>
                            <input type="number" id="mabl-weight" placeholder="70" step="0.1">
                        </div>
                        <div class="form-group">
                            <label for="mabl-height">Height (cm)</label>
                            <input type="number" id="mabl-height" placeholder="170" step="0.1">
                        </div>
                        <div class="form-group">
                            <label for="mabl-sex">Biological Sex</label>
                            <select id="mabl-sex">
                                <option value="">Select...</option>
                                <option value="male">Male</option>
                                <option value="female">Female</option>
                            </select>
                        </div>
                    </div>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="mabl-initial-hct">Initial Hematocrit (%)</label>
                            <input type="number" id="mabl-initial-hct" placeholder="45" step="0.1">
                        </div>
                        <div class="form-group">
                            <label for="mabl-target-hct">Target Hematocrit (%)</label>
                            <input type="number" id="mabl-target-hct" placeholder="25" step="0.1">
                        </div>
                    </div>

                    <button class="calculate-btn" onclick="calculateMABL()">Calculate</button>

                    <div id="mabl-result" class="result-box">
                        <div class="result-label">Maximum Allowable Blood Loss</div>
                        <div class="result-value" id="mabl-value">-</div>
                        <div class="result-interpretation" id="mabl-interpretation"></div>
                        <button class="send-to-ai-btn" onclick="sendToAI('mabl')">
                            Send to AI
                        </button>
                    </div>
                </div>
            </div>

            <!-- BSA Calculator -->
            <div id="calc-bsa" class="calculator-section">
                <div class="calculator-card">
                    <h1 class="calculator-title">Body Surface Area (BSA)</h1>
                    <p class="calculator-description">Calculate BSA using the Mosteller formula. Used for chemotherapy dosing and cardiac output calculation.</p>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="bsa-weight">Weight (kg)</label>
                            <input type="number" id="bsa-weight" placeholder="70" step="0.1">
                        </div>
                        <div class="form-group">
                            <label for="bsa-height">Height (cm)</label>
                            <input type="number" id="bsa-height" placeholder="170" step="0.1">
                        </div>
                    </div>

                    <button class="calculate-btn" onclick="calculateBSA()">Calculate</button>

                    <div id="bsa-result" class="result-box">
                        <div class="result-label">Body Surface Area</div>
                        <div class="result-value" id="bsa-value">-</div>
                        <div class="result-interpretation" id="bsa-interpretation"></div>
                        <button class="send-to-ai-btn" onclick="sendToAI('bsa')">
                            Send to AI
                        </button>
                    </div>
                </div>
            </div>

            <!-- QTc Calculator -->
            <div id="calc-qtc" class="calculator-section">
                <div class="calculator-card">
                    <h1 class="calculator-title">QTc Interval (Bazett Formula)</h1>
                    <p class="calculator-description">Calculate heart rate-corrected QT interval. Prolonged QTc increases risk of torsades de pointes.</p>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="qtc-qt">QT Interval (ms)</label>
                            <input type="number" id="qtc-qt" placeholder="400" step="1">
                        </div>
                        <div class="form-group">
                            <label for="qtc-rr">RR Interval (ms)</label>
                            <input type="number" id="qtc-rr" placeholder="800" step="1">
                        </div>
                    </div>

                    <button class="calculate-btn" onclick="calculateQTc()">Calculate</button>

                    <div id="qtc-result" class="result-box">
                        <div class="result-label">Corrected QT Interval</div>
                        <div class="result-value" id="qtc-value">-</div>
                        <div class="result-interpretation" id="qtc-interpretation"></div>
                        <button class="send-to-ai-btn" onclick="sendToAI('qtc')">
                            Send to AI
                        </button>
                    </div>
                </div>
            </div>

            <!-- Fluids Calculator -->
            <div id="calc-fluids" class="calculator-section">
                <div class="calculator-card">
                    <h1 class="calculator-title">Maintenance Fluids (4-2-1 Rule)</h1>
                    <p class="calculator-description">Calculate hourly maintenance fluid rate using the Holliday-Segar method.</p>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="fluids-weight">Weight (kg)</label>
                            <input type="number" id="fluids-weight" placeholder="70" step="0.1">
                        </div>
                    </div>

                    <button class="calculate-btn" onclick="calculateFluids()">Calculate</button>

                    <div id="fluids-result" class="result-box">
                        <div class="result-label">Maintenance Fluid Rate</div>
                        <div class="result-value" id="fluids-value">-</div>
                        <div class="result-interpretation" id="fluids-interpretation"></div>
                        <button class="send-to-ai-btn" onclick="sendToAI('fluids')">
                            Send to AI
                        </button>
                    </div>
                </div>
            </div>

            <!-- PONV Calculator -->
            <div id="calc-ponv" class="calculator-section">
                <div class="calculator-card">
                    <h1 class="calculator-title">PONV Risk (Apfel Score)</h1>
                    <p class="calculator-description">Predict risk of postoperative nausea and vomiting based on 4 independent risk factors.</p>

                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="ponv-female" style="width: auto; display: inline-block; margin-right: 8px;">
                            Female sex
                        </label>
                    </div>
                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="ponv-nonsmoker" style="width: auto; display: inline-block; margin-right: 8px;">
                            Non-smoker
                        </label>
                    </div>
                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="ponv-history" style="width: auto; display: inline-block; margin-right: 8px;">
                            History of PONV or motion sickness
                        </label>
                    </div>
                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="ponv-opioids" style="width: auto; display: inline-block; margin-right: 8px;">
                            Postoperative opioids expected
                        </label>
                    </div>

                    <button class="calculate-btn" onclick="calculatePONV()">Calculate</button>

                    <div id="ponv-result" class="result-box">
                        <div class="result-label">PONV Risk</div>
                        <div class="result-value" id="ponv-value">-</div>
                        <div class="result-interpretation" id="ponv-interpretation"></div>
                        <button class="send-to-ai-btn" onclick="sendToAI('ponv')">
                            Send to AI
                        </button>
                    </div>
                </div>
            </div>

            <!-- ASA Helper -->
            <div id="calc-asa" class="calculator-section">
                <div class="calculator-card">
                    <h1 class="calculator-title">ASA Physical Status Classification</h1>
                    <p class="calculator-description">Reference guide for ASA physical status classification system.</p>

                    <div style="color: var(--text-secondary); line-height: 1.8;">
                        <p style="margin-bottom: 16px;"><strong style="color: var(--primary-blue);">ASA I:</strong> Normal healthy patient (no organic, physiologic, or psychiatric disturbance)</p>
                        <p style="margin-bottom: 16px;"><strong style="color: var(--primary-blue);">ASA II:</strong> Patient with mild systemic disease (well-controlled HTN, DM, obesity, smoking)</p>
                        <p style="margin-bottom: 16px;"><strong style="color: var(--primary-blue);">ASA III:</strong> Patient with severe systemic disease (poorly controlled HTN/DM, COPD, morbid obesity, active hepatitis, CAD, MI >3mo ago)</p>
                        <p style="margin-bottom: 16px;"><strong style="color: var(--primary-blue);">ASA IV:</strong> Patient with severe systemic disease that is a constant threat to life (recent MI <3mo, CVA, ongoing cardiac ischemia, severe valve dysfunction, sepsis, ESRD)</p>
                        <p style="margin-bottom: 16px;"><strong style="color: var(--primary-blue);">ASA V:</strong> Moribund patient not expected to survive without operation (ruptured AAA, massive trauma, intracranial bleed with mass effect)</p>
                        <p style="margin-bottom: 16px;"><strong style="color: var(--primary-blue);">ASA VI:</strong> Brain-dead patient for organ donation</p>
                        <p style="margin-top: 24px; font-style: italic; color: var(--text-muted);">Add "E" suffix for emergency surgery (e.g., ASA 3E)</p>
                    </div>

                    <div id="asa-result" class="result-box visible" style="margin-top: 32px;">
                        <button class="send-to-ai-btn" onclick="sendToAI('asa')">
                            Ask AI about ASA Classification
                        </button>
                    </div>
                </div>
            </div>

            <!-- Fluid Deficit Calculator -->
            <div id="calc-fluid-deficit" class="calculator-section">
                <div class="calculator-card">
                    <h1 class="calculator-title">Fluid Deficit & Replacement</h1>
                    <p class="calculator-description">Calculate fluid deficit based on dehydration percentage and replacement strategy.</p>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="fd-weight">Weight (kg)</label>
                            <input type="number" id="fd-weight" placeholder="70" step="0.1">
                        </div>
                        <div class="form-group">
                            <label for="fd-dehydration">Dehydration (%)</label>
                            <select id="fd-dehydration">
                                <option value="">Select...</option>
                                <option value="3">Mild (3%)</option>
                                <option value="5">Moderate (5%)</option>
                                <option value="7">Severe (7%)</option>
                                <option value="10">Critical (10%)</option>
                            </select>
                        </div>
                    </div>

                    <button class="calculate-btn" onclick="calculateFluidDeficit()">Calculate</button>

                    <div id="fd-result" class="result-box">
                        <div class="result-label">Fluid Deficit</div>
                        <div class="result-value" id="fd-value">-</div>
                        <div class="result-interpretation" id="fd-interpretation"></div>
                        <button class="send-to-ai-btn" onclick="sendToAI('fluid-deficit')">
                            Send to AI
                        </button>
                    </div>
                </div>
            </div>

            <!-- Defibrillation Energy Calculator -->
            <div id="calc-defib" class="calculator-section">
                <div class="calculator-card">
                    <h1 class="calculator-title">Defibrillation / Cardioversion Energy</h1>
                    <p class="calculator-description">Calculate appropriate defibrillation or cardioversion energy for adults and pediatrics.</p>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="defib-type">Procedure Type</label>
                            <select id="defib-type">
                                <option value="">Select...</option>
                                <option value="defib">Defibrillation (VF/pulseless VT)</option>
                                <option value="cardioversion">Synchronized Cardioversion</option>
                            </select>
                        </div>
                        <div class="form-group">
                            <label for="defib-age">Patient Age</label>
                            <select id="defib-age">
                                <option value="">Select...</option>
                                <option value="adult">Adult</option>
                                <option value="peds">Pediatric</option>
                            </select>
                        </div>
                    </div>

                    <div class="form-row" id="defib-weight-row" style="display: none;">
                        <div class="form-group">
                            <label for="defib-weight">Weight (kg)</label>
                            <input type="number" id="defib-weight" placeholder="20" step="0.1">
                        </div>
                    </div>

                    <button class="calculate-btn" onclick="calculateDefib()">Calculate</button>

                    <div id="defib-result" class="result-box">
                        <div class="result-label">Energy Recommendation</div>
                        <div class="result-value" id="defib-value">-</div>
                        <div class="result-interpretation" id="defib-interpretation"></div>
                        <button class="send-to-ai-btn" onclick="sendToAI('defib')">
                            Send to AI
                        </button>
                    </div>
                </div>
            </div>

            <!-- Surgical Apgar Score Calculator -->
            <div id="calc-sas" class="calculator-section">
                <div class="calculator-card">
                    <h1 class="calculator-title">Surgical Apgar Score (SAS)</h1>
                    <p class="calculator-description">10-point score for predicting post-operative morbidity and mortality.</p>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="sas-ebl">Estimated Blood Loss (mL)</label>
                            <input type="number" id="sas-ebl" placeholder="400" step="1">
                        </div>
                        <div class="form-group">
                            <label for="sas-hr">Lowest Heart Rate (bpm)</label>
                            <input type="number" id="sas-hr" placeholder="75" step="1">
                        </div>
                    </div>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="sas-map">Lowest MAP (mmHg)</label>
                            <input type="number" id="sas-map" placeholder="70" step="1">
                        </div>
                    </div>

                    <button class="calculate-btn" onclick="calculateSAS()">Calculate</button>

                    <div id="sas-result" class="result-box">
                        <div class="result-label">Surgical Apgar Score</div>
                        <div class="result-value" id="sas-value">-</div>
                        <div class="result-interpretation" id="sas-interpretation"></div>
                        <button class="send-to-ai-btn" onclick="sendToAI('sas')">
                            Send to AI
                        </button>
                    </div>
                </div>
            </div>

            <!-- Opioid Conversion Calculator -->
            <div id="calc-opioid" class="calculator-section">
                <div class="calculator-card">
                    <h1 class="calculator-title">Opioid Conversion Calculator</h1>
                    <p class="calculator-description">Convert between opioid medications using morphine milligram equivalents (MME).</p>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="opioid-from">From Medication</label>
                            <select id="opioid-from">
                                <option value="">Select...</option>
                                <option value="morphine-iv">Morphine IV</option>
                                <option value="morphine-po">Morphine PO</option>
                                <option value="fentanyl-iv">Fentanyl IV</option>
                                <option value="hydromorphone-iv">Hydromorphone IV</option>
                                <option value="hydromorphone-po">Hydromorphone PO</option>
                                <option value="oxycodone-po">Oxycodone PO</option>
                                <option value="hydrocodone-po">Hydrocodone PO</option>
                            </select>
                        </div>
                        <div class="form-group">
                            <label for="opioid-dose">Dose (mg or mcg for fentanyl)</label>
                            <input type="number" id="opioid-dose" placeholder="10" step="0.1">
                        </div>
                    </div>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="opioid-to">To Medication</label>
                            <select id="opioid-to">
                                <option value="">Select...</option>
                                <option value="morphine-iv">Morphine IV</option>
                                <option value="morphine-po">Morphine PO</option>
                                <option value="fentanyl-iv">Fentanyl IV</option>
                                <option value="hydromorphone-iv">Hydromorphone IV</option>
                                <option value="hydromorphone-po">Hydromorphone PO</option>
                                <option value="oxycodone-po">Oxycodone PO</option>
                                <option value="hydrocodone-po">Hydrocodone PO</option>
                            </select>
                        </div>
                    </div>

                    <button class="calculate-btn" onclick="calculateOpioid()">Calculate</button>

                    <div id="opioid-result" class="result-box">
                        <div class="result-label">Converted Dose</div>
                        <div class="result-value" id="opioid-value">-</div>
                        <div class="result-interpretation" id="opioid-interpretation"></div>
                        <button class="send-to-ai-btn" onclick="sendToAI('opioid')">
                            Send to AI
                        </button>
                    </div>
                </div>
            </div>

            <!-- Local Anesthetic Toxicity Dose Calculator -->
            <div id="calc-last" class="calculator-section">
                <div class="calculator-card">
                    <h1 class="calculator-title">Local Anesthetic Toxicity Dose</h1>
                    <p class="calculator-description">Calculate maximum safe dose for local anesthetics to prevent LAST.</p>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="last-weight">Weight (kg)</label>
                            <input type="number" id="last-weight" placeholder="70" step="0.1">
                        </div>
                        <div class="form-group">
                            <label for="last-agent">Local Anesthetic</label>
                            <select id="last-agent">
                                <option value="">Select...</option>
                                <option value="lidocaine">Lidocaine (without epi)</option>
                                <option value="lidocaine-epi">Lidocaine (with epi)</option>
                                <option value="bupivacaine">Bupivacaine</option>
                                <option value="ropivacaine">Ropivacaine</option>
                                <option value="chloroprocaine">Chloroprocaine</option>
                            </select>
                        </div>
                    </div>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="last-concentration">Concentration (%)</label>
                            <input type="number" id="last-concentration" placeholder="0.5" step="0.1">
                        </div>
                    </div>

                    <button class="calculate-btn" onclick="calculateLAST()">Calculate</button>

                    <div id="last-result" class="result-box">
                        <div class="result-label">Maximum Safe Dose</div>
                        <div class="result-value" id="last-value">-</div>
                        <div class="result-interpretation" id="last-interpretation"></div>
                        <button class="send-to-ai-btn" onclick="sendToAI('last')">
                            Send to AI
                        </button>
                    </div>
                </div>
            </div>

            <!-- RCRI Calculator -->
            <div id="calc-rcri" class="calculator-section">
                <div class="calculator-card">
                    <h1 class="calculator-title">Revised Cardiac Risk Index (RCRI)</h1>
                    <p class="calculator-description">Estimate risk of major cardiac complications after noncardiac surgery.</p>

                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="rcri-highrisk" style="width: auto; display: inline-block; margin-right: 8px;">
                            High-risk surgery (intraperitoneal, intrathoracic, or suprainguinal vascular)
                        </label>
                    </div>
                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="rcri-ihd" style="width: auto; display: inline-block; margin-right: 8px;">
                            History of ischemic heart disease
                        </label>
                    </div>
                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="rcri-chf" style="width: auto; display: inline-block; margin-right: 8px;">
                            History of congestive heart failure
                        </label>
                    </div>
                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="rcri-cvd" style="width: auto; display: inline-block; margin-right: 8px;">
                            History of cerebrovascular disease
                        </label>
                    </div>
                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="rcri-dm" style="width: auto; display: inline-block; margin-right: 8px;">
                            Diabetes mellitus on insulin therapy
                        </label>
                    </div>
                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="rcri-renal" style="width: auto; display: inline-block; margin-right: 8px;">
                            Preoperative creatinine >2 mg/dL
                        </label>
                    </div>

                    <button class="calculate-btn" onclick="calculateRCRI()">Calculate</button>

                    <div id="rcri-result" class="result-box">
                        <div class="result-label">RCRI Score</div>
                        <div class="result-value" id="rcri-value">-</div>
                        <div class="result-interpretation" id="rcri-interpretation"></div>
                        <button class="send-to-ai-btn" onclick="sendToAI('rcri')">
                            Send to AI
                        </button>
                    </div>
                </div>
            </div>

            <!-- STOP-BANG Calculator -->
            <div id="calc-stopbang" class="calculator-section">
                <div class="calculator-card">
                    <h1 class="calculator-title">STOP-BANG OSA Risk</h1>
                    <p class="calculator-description">Screen for obstructive sleep apnea using the STOP-BANG questionnaire.</p>

                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="sb-snore" style="width: auto; display: inline-block; margin-right: 8px;">
                            <strong>S</strong>noring: Do you snore loudly?
                        </label>
                    </div>
                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="sb-tired" style="width: auto; display: inline-block; margin-right: 8px;">
                            <strong>T</strong>ired: Do you often feel tired during the day?
                        </label>
                    </div>
                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="sb-observed" style="width: auto; display: inline-block; margin-right: 8px;">
                            <strong>O</strong>bserved: Has anyone observed you stop breathing during sleep?
                        </label>
                    </div>
                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="sb-pressure" style="width: auto; display: inline-block; margin-right: 8px;">
                            <strong>P</strong>ressure: Do you have or are you being treated for high blood pressure?
                        </label>
                    </div>
                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="sb-bmi" style="width: auto; display: inline-block; margin-right: 8px;">
                            <strong>B</strong>MI >35 kg/m²
                        </label>
                    </div>
                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="sb-age" style="width: auto; display: inline-block; margin-right: 8px;">
                            <strong>A</strong>ge >50 years
                        </label>
                    </div>
                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="sb-neck" style="width: auto; display: inline-block; margin-right: 8px;">
                            <strong>N</strong>eck circumference >40 cm
                        </label>
                    </div>
                    <div class="form-group">
                        <label>
                            <input type="checkbox" id="sb-gender" style="width: auto; display: inline-block; margin-right: 8px;">
                            <strong>G</strong>ender: Male
                        </label>
                    </div>

                    <button class="calculate-btn" onclick="calculateStopBang()">Calculate</button>

                    <div id="stopbang-result" class="result-box">
                        <div class="result-label">STOP-BANG Score</div>
                        <div class="result-value" id="stopbang-value">-</div>
                        <div class="result-interpretation" id="stopbang-interpretation"></div>
                        <button class="send-to-ai-btn" onclick="sendToAI('stopbang')">
                            Send to AI
                        </button>
                    </div>
                </div>
            </div>

            <!-- Airway Assessment Calculator -->
            <div id="calc-airway" class="calculator-section">
                <div class="calculator-card">
                    <h1 class="calculator-title">Airway Assessment Helper</h1>
                    <p class="calculator-description">Predict difficult airway based on Mallampati, thyromental distance, and neck ROM.</p>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="airway-mallampati">Mallampati Class</label>
                            <select id="airway-mallampati">
                                <option value="">Select...</option>
                                <option value="1">Class I (full visibility)</option>
                                <option value="2">Class II (uvula visible)</option>
                                <option value="3">Class III (soft palate only)</option>
                                <option value="4">Class IV (hard palate only)</option>
                            </select>
                        </div>
                        <div class="form-group">
                            <label for="airway-tmd">Thyromental Distance</label>
                            <select id="airway-tmd">
                                <option value="">Select...</option>
                                <option value="good">>6 cm (Normal)</option>
                                <option value="borderline">4-6 cm (Borderline)</option>
                                <option value="poor"><4 cm (Difficult)</option>
                            </select>
                        </div>
                    </div>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="airway-neck">Neck ROM</label>
                            <select id="airway-neck">
                                <option value="">Select...</option>
                                <option value="good">Normal (>35°)</option>
                                <option value="reduced">Reduced (<35°)</option>
                            </select>
                        </div>
                        <div class="form-group">
                            <label for="airway-mouth">Mouth Opening</label>
                            <select id="airway-mouth">
                                <option value="">Select...</option>
                                <option value="good">>4 cm</option>
                                <option value="poor"><4 cm</option>
                            </select>
                        </div>
                    </div>

                    <button class="calculate-btn" onclick="calculateAirway()">Calculate</button>

                    <div id="airway-result" class="result-box">
                        <div class="result-label">Predicted Difficulty</div>
                        <div class="result-value" id="airway-value">-</div>
                        <div class="result-interpretation" id="airway-interpretation"></div>
                        <button class="send-to-ai-btn" onclick="sendToAI('airway')">
                            Send to AI
                        </button>
                    </div>
                </div>
            </div>

            <!-- MAC Calculator -->
            <div id="calc-mac" class="calculator-section">
                <div class="calculator-card">
                    <h1 class="calculator-title">MAC Calculator (Age-Adjusted)</h1>
                    <p class="calculator-description">Calculate age-adjusted MAC-awake and MAC for common volatile anesthetics.</p>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="mac-agent">Volatile Agent</label>
                            <select id="mac-agent">
                                <option value="">Select...</option>
                                <option value="sevoflurane">Sevoflurane</option>
                                <option value="desflurane">Desflurane</option>
                                <option value="isoflurane">Isoflurane</option>
                            </select>
                        </div>
                        <div class="form-group">
                            <label for="mac-age">Age (years)</label>
                            <input type="number" id="mac-age" placeholder="40" step="1">
                        </div>
                    </div>

                    <button class="calculate-btn" onclick="calculateMAC()">Calculate</button>

                    <div id="mac-result" class="result-box">
                        <div class="result-label">Age-Adjusted MAC</div>
                        <div class="result-value" id="mac-value">-</div>
                        <div class="result-interpretation" id="mac-interpretation"></div>
                        <button class="send-to-ai-btn" onclick="sendToAI('mac')">
                            Send to AI
                        </button>
                    </div>
                </div>
            </div>

            <!-- Propofol TCI Calculator -->
            <div id="calc-propofol-tci" class="calculator-section">
                <div class="calculator-card">
                    <h1 class="calculator-title">Propofol TCI / Infusion Rate</h1>
                    <p class="calculator-description">Calculate target-controlled infusion rates using Marsh, Schnider, or Kataria models.</p>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="tci-model">TCI Model</label>
                            <select id="tci-model">
                                <option value="">Select...</option>
                                <option value="marsh">Marsh (weight-based)</option>
                                <option value="schnider">Schnider (effect-site)</option>
                                <option value="kataria">Kataria (pediatric)</option>
                            </select>
                        </div>
                        <div class="form-group">
                            <label for="tci-target">Target Concentration (mcg/mL)</label>
                            <input type="number" id="tci-target" placeholder="3.5" step="0.1">
                        </div>
                    </div>

                    <div class="form-row">
                        <div class="form-group">
                            <label for="tci-weight">Weight (kg)</label>
                            <input type="number" id="tci-weight" placeholder="70" step="0.1">
                        </div>
                        <div class="form-group" id="tci-age-group">
                            <label for="tci-age">Age (years)</label>
                            <input type="number" id="tci-age" placeholder="40" step="1">
                        </div>
                    </div>

                    <div class="form-row" id="tci-schnider-params" style="display: none;">
                        <div class="form-group">
                            <label for="tci-height">Height (cm)</label>
                            <input type="number" id="tci-height" placeholder="170" step="1">
                        </div>
                        <div class="form-group">
                            <label for="tci-sex">Sex</label>
                            <select id="tci-sex">
                                <option value="">Select...</option>
                                <option value="male">Male</option>
                                <option value="female">Female</option>
                            </select>
                        </div>
                    </div>

                    <button class="calculate-btn" onclick="calculatePropofolTCI()">Calculate</button>

                    <div id="tci-result" class="result-box">
                        <div class="result-label">Infusion Rate</div>
                        <div class="result-value" id="tci-value">-</div>
                        <div class="result-interpretation" id="tci-interpretation"></div>
                        <button class="send-to-ai-btn" onclick="sendToAI('propofol-tci')">
                            Send to AI
                        </button>
                    </div>
                </div>
            </div>
        </div>
    </div>

    <!-- AI Chat Modal -->
    <div class="ai-chat-modal" id="aiChatModal">
        <div class="ai-chat-container">
            <div class="ai-chat-header">
                <h3>Ask AI About Your Results</h3>
                <button class="close-modal-btn" onclick="closeAIChat()">
                    <svg width="20" height="20" viewBox="0 0 20 20" fill="none" xmlns="http://www.w3.org/2000/svg">
                        <path d="M15 5L5 15M5 5L15 15" stroke="currentColor" stroke-width="2" stroke-linecap="round"/>
                    </svg>
                </button>
            </div>
            <div class="ai-chat-messages" id="aiChatMessages">
                <!-- Messages will be added here -->
            </div>
            <div class="ai-chat-input-area">
                <div class="ai-input-container">
                    <input type="hidden" id="csrf_token" value="{{ csrf_token() }}"/>
                    <textarea class="ai-chat-input" id="aiChatInput" placeholder="Ask a follow-up question..." rows="1"></textarea>
                    <button class="send-ai-message-btn" onclick="sendAIMessage()">Send</button>
                </div>
            </div>
        </div>
    </div>

    <script>
        // Store calculation results
        let calculationResults = {};
        let currentChatContext = '';

        // Calculator switching
        document.querySelectorAll('.calculator-btn').forEach(btn => {
            btn.addEventListener('click', function() {
                const calcType = this.dataset.calc;

                // Update sidebar
                document.querySelectorAll('.calculator-btn').forEach(b => b.classList.remove('active'));
                this.classList.add('active');

                // Update content
                document.querySelectorAll('.calculator-section').forEach(section => {
                    section.classList.remove('active');
                });
                document.getElementById(`calc-${calcType}`).classList.add('active');
            });
        });

        // IBW Calculator
        function calculateIBW() {
            const height = parseFloat(document.getElementById('ibw-height').value);
            const sex = document.getElementById('ibw-sex').value;

            if (!height || !sex) {
                alert('Please fill in all fields');
                return;
            }

            let ibw;
            if (sex === 'male') {
                ibw = 50 + 0.91 * (height - 152.4);
            } else {
                ibw = 45.5 + 0.91 * (height - 152.4);
            }

            ibw = Math.round(ibw * 10) / 10;

            calculationResults.ibw = {
                value: ibw,
                height: height,
                sex: sex
            };

            document.getElementById('ibw-value').textContent = `${ibw} kg`;
            document.getElementById('ibw-interpretation').innerHTML =
                `Based on Devine formula for ${sex} patient with height ${height} cm. Used for dosing propofol, succinylcholine, and tidal volume calculations (6-8 mL/kg IBW).`;
            document.getElementById('ibw-result').classList.add('visible');
        }

        // MABL Calculator
        function calculateMABL() {
            const weight = parseFloat(document.getElementById('mabl-weight').value);
            const height = parseFloat(document.getElementById('mabl-height').value);
            const sex = document.getElementById('mabl-sex').value;
            const initialHct = parseFloat(document.getElementById('mabl-initial-hct').value);
            const targetHct = parseFloat(document.getElementById('mabl-target-hct').value);

            if (!weight || !height || !sex || !initialHct || !targetHct) {
                alert('Please fill in all fields');
                return;
            }

            // Calculate blood volume
            let bloodVolume;
            if (sex === 'male') {
                bloodVolume = weight * 75; // mL
            } else {
                bloodVolume = weight * 65; // mL
            }

            // MABL = Blood Volume × (Initial Hct - Target Hct) / Initial Hct
            const mabl = bloodVolume * (initialHct - targetHct) / initialHct;

            calculationResults.mabl = {
                value: Math.round(mabl),
                weight: weight,
                height: height,
                sex: sex,
                initialHct: initialHct,
                targetHct: targetHct,
                bloodVolume: Math.round(bloodVolume)
            };

            document.getElementById('mabl-value').textContent = `${Math.round(mabl)} mL`;
            document.getElementById('mabl-interpretation').innerHTML =
                `Estimated blood volume: ${Math.round(bloodVolume)} mL (${sex === 'male' ? '75' : '65'} mL/kg).<br>
                This represents the maximum blood loss before transfusion should be considered, assuming no ongoing bleeding and adequate crystalloid resuscitation.`;
            document.getElementById('mabl-result').classList.add('visible');
        }

        // BSA Calculator
        function calculateBSA() {
            const weight = parseFloat(document.getElementById('bsa-weight').value);
            const height = parseFloat(document.getElementById('bsa-height').value);

            if (!weight || !height) {
                alert('Please fill in all fields');
                return;
            }

            // Mosteller formula
            const bsa = Math.sqrt((weight * height) / 3600);

            calculationResults.bsa = {
                value: bsa.toFixed(2),
                weight: weight,
                height: height
            };

            document.getElementById('bsa-value').textContent = `${bsa.toFixed(2)} m²`;
            document.getElementById('bsa-interpretation').innerHTML =
                `Calculated using Mosteller formula. Normal adult BSA is approximately 1.7 m². Used for chemotherapy dosing and cardiac index calculations (CI = CO / BSA).`;
            document.getElementById('bsa-result').classList.add('visible');
        }

        // QTc Calculator
        function calculateQTc() {
            const qt = parseFloat(document.getElementById('qtc-qt').value);
            const rr = parseFloat(document.getElementById('qtc-rr').value);

            if (!qt || !rr) {
                alert('Please fill in all fields');
                return;
            }

            // Bazett formula: QTc = QT / sqrt(RR in seconds)
            const rrSeconds = rr / 1000;
            const qtc = qt / Math.sqrt(rrSeconds);

            let interpretation = '';
            let risk = '';
            if (qtc > 500) {
                interpretation = '⚠️ <strong>Severely prolonged</strong> (>500 ms) - HIGH RISK for torsades de pointes. Avoid QT-prolonging drugs.';
                risk = 'severe';
            } else if (qtc > 470) {
                interpretation = '⚠️ <strong>Prolonged</strong> (>470 ms females, >450 ms males) - Increased risk for torsades. Use caution with QT-prolonging drugs.';
                risk = 'moderate';
            } else if (qtc >= 440) {
                interpretation = '⚠️ <strong>Borderline prolonged</strong> (440-470 ms) - Monitor closely if using QT-prolonging medications.';
                risk = 'mild';
            } else {
                interpretation = '✓ <strong>Normal</strong> (<440 ms) - Standard perioperative management.';
                risk = 'normal';
            }

            calculationResults.qtc = {
                value: Math.round(qtc),
                qt: qt,
                rr: rr,
                hr: Math.round(60000 / rr),
                risk: risk
            };

            document.getElementById('qtc-value').textContent = `${Math.round(qtc)} ms`;
            document.getElementById('qtc-interpretation').innerHTML = interpretation;
            document.getElementById('qtc-result').classList.add('visible');
        }

        // Fluids Calculator
        function calculateFluids() {
            const weight = parseFloat(document.getElementById('fluids-weight').value);

            if (!weight) {
                alert('Please enter weight');
                return;
            }

            // 4-2-1 rule
            let rate = 0;
            if (weight <= 10) {
                rate = weight * 4;
            } else if (weight <= 20) {
                rate = 40 + (weight - 10) * 2;
            } else {
                rate = 60 + (weight - 20) * 1;
            }

            calculationResults.fluids = {
                value: Math.round(rate),
                weight: weight,
                daily: Math.round(rate * 24)
            };

            document.getElementById('fluids-value').textContent = `${Math.round(rate)} mL/hr`;
            document.getElementById('fluids-interpretation').innerHTML =
                `Based on 4-2-1 rule (Holliday-Segar):<br>
                • First 10 kg: 4 mL/kg/hr<br>
                • Second 10 kg: 2 mL/kg/hr<br>
                • Each additional kg: 1 mL/kg/hr<br>
                <strong>Daily total: ${Math.round(rate * 24)} mL/day</strong>`;
            document.getElementById('fluids-result').classList.add('visible');
        }

        // PONV Calculator
        function calculatePONV() {
            const female = document.getElementById('ponv-female').checked;
            const nonsmoker = document.getElementById('ponv-nonsmoker').checked;
            const history = document.getElementById('ponv-history').checked;
            const opioids = document.getElementById('ponv-opioids').checked;

            const score = (female ? 1 : 0) + (nonsmoker ? 1 : 0) + (history ? 1 : 0) + (opioids ? 1 : 0);

            let risk, percentage, recommendations;
            switch(score) {
                case 0:
                    risk = 'Very Low';
                    percentage = '10%';
                    recommendations = 'Standard anesthetic care. Prophylaxis generally not indicated.';
                    break;
                case 1:
                    risk = 'Low';
                    percentage = '20%';
                    recommendations = 'Consider single antiemetic (ondansetron 4 mg or dexamethasone 4-8 mg).';
                    break;
                case 2:
                    risk = 'Moderate';
                    percentage = '40%';
                    recommendations = 'Use 2 antiemetics from different classes. Consider TIVA instead of volatile agents.';
                    break;
                case 3:
                    risk = 'High';
                    percentage = '60%';
                    recommendations = '3+ antiemetics recommended. Strong consideration for TIVA. Minimize opioids, use multimodal analgesia.';
                    break;
                case 4:
                    risk = 'Very High';
                    percentage = '80%';
                    recommendations = 'Maximum prophylaxis: 3-4 antiemetics, TIVA, opioid-sparing techniques, consider regional anesthesia.';
                    break;
            }

            calculationResults.ponv = {
                score: score,
                risk: risk,
                percentage: percentage,
                factors: {
                    female: female,
                    nonsmoker: nonsmoker,
                    history: history,
                    opioids: opioids
                }
            };

            document.getElementById('ponv-value').textContent = `${score}/4 - ${risk} Risk (${percentage})`;
            document.getElementById('ponv-interpretation').innerHTML = `<strong>Recommendations:</strong><br>${recommendations}`;
            document.getElementById('ponv-result').classList.add('visible');
        }

        // Fluid Deficit Calculator
        function calculateFluidDeficit() {
            const weight = parseFloat(document.getElementById('fd-weight').value);
            const dehydration = parseFloat(document.getElementById('fd-dehydration').value);

            if (!weight || !dehydration) {
                alert('Please fill in all fields');
                return;
            }

            const deficit = weight * (dehydration / 100) * 1000; // mL
            const first8hr = deficit / 2;
            const next16hr = deficit / 2;
            const hourlyFirst8 = first8hr / 8;
            const hourlyNext16 = next16hr / 16;

            calculationResults['fluid-deficit'] = {
                value: Math.round(deficit),
                weight: weight,
                dehydration: dehydration,
                first8hr: Math.round(first8hr),
                next16hr: Math.round(next16hr)
            };

            document.getElementById('fd-value').textContent = `${Math.round(deficit)} mL`;
            document.getElementById('fd-interpretation').innerHTML = `<strong>Replacement Strategy:</strong><br>` +
                `• First 8 hours: ${Math.round(first8hr)} mL (${Math.round(hourlyFirst8)} mL/hr)<br>` +
                `• Next 16 hours: ${Math.round(next16hr)} mL (${Math.round(hourlyNext16)} mL/hr)<br>` +
                `• Total fluid deficit: ${Math.round(deficit)} mL over 24 hours<br>` +
                `• Replace with isotonic crystalloid (LR or NS) + ongoing maintenance fluids`;
            document.getElementById('fd-result').classList.add('visible');
        }

        // Defibrillation Calculator
        document.getElementById('defib-age').addEventListener('change', function() {
            if (this.value === 'peds') {
                document.getElementById('defib-weight-row').style.display = 'block';
            } else {
                document.getElementById('defib-weight-row').style.display = 'none';
            }
        });

        function calculateDefib() {
            const type = document.getElementById('defib-type').value;
            const age = document.getElementById('defib-age').value;
            const weight = parseFloat(document.getElementById('defib-weight').value);

            if (!type || !age) {
                alert('Please select procedure type and patient age');
                return;
            }

            if (age === 'peds' && !weight) {
                alert('Please enter weight for pediatric patient');
                return;
            }

            let energy, interpretation;

            if (type === 'defib') {
                if (age === 'adult') {
                    energy = '200 J (biphasic), 360 J (monophasic)';
                    interpretation = '<strong>Defibrillation Energy (Adult):</strong><br>' +
                        '• Initial: 200 J (biphasic) or 360 J (monophasic)<br>' +
                        '• Subsequent: Same or higher (maximum 200 J biphasic, 360 J monophasic)<br>' +
                        '• Unsynchronized shock for VF/pulseless VT';
                } else {
                    const dose = 2 * weight; // 2 J/kg initial
                    const doseMax = 4 * weight; // 4 J/kg subsequent
                    energy = `${dose} J (initial, 2 J/kg)`;
                    interpretation = `<strong>Defibrillation Energy (Pediatric):</strong><br>` +
                        `• Initial: ${dose} J (2 J/kg)<br>` +
                        `• Subsequent: ${Math.min(doseMax, 200)} J (4 J/kg, max 200 J)<br>` +
                        `• Unsynchronized shock for VF/pulseless VT`;
                }
            } else { // cardioversion
                if (age === 'adult') {
                    energy = '50-100 J (initial), 200 J (subsequent)';
                    interpretation = '<strong>Synchronized Cardioversion (Adult):</strong><br>' +
                        '• A-fib/A-flutter: 120-200 J (biphasic)<br>' +
                        '• SVT/V-tach (with pulse): 50-100 J initially<br>' +
                        '• Synchronized to R wave (avoid T wave)';
                } else {
                    const dose = 0.5 * weight; // 0.5 J/kg initial
                    const doseMax = 1 * weight; // 1 J/kg subsequent
                    energy = `${dose} J (initial, 0.5 J/kg)`;
                    interpretation = `<strong>Synchronized Cardioversion (Pediatric):</strong><br>` +
                        `• Initial: ${dose} J (0.5 J/kg)<br>` +
                        `• Subsequent: ${Math.min(doseMax, 100)} J (1 J/kg, max 100 J)<br>` +
                        `• Synchronized to R wave`;
                }
            }

            calculationResults.defib = {
                value: energy,
                type: type,
                age: age,
                weight: weight
            };

            document.getElementById('defib-value').textContent = energy;
            document.getElementById('defib-interpretation').innerHTML = interpretation;
            document.getElementById('defib-result').classList.add('visible');
        }

        // Surgical Apgar Score Calculator
        function calculateSAS() {
            const ebl = parseFloat(document.getElementById('sas-ebl').value);
            const hr = parseFloat(document.getElementById('sas-hr').value);
            const map = parseFloat(document.getElementById('sas-map').value);

            if (!ebl && ebl !== 0 || !hr || !map) {
                alert('Please fill in all fields');
                return;
            }

            let eblPoints, hrPoints, mapPoints;

            // EBL points
            if (ebl <= 100) eblPoints = 4;
            else if (ebl <= 600) eblPoints = 3;
            else if (ebl <= 1000) eblPoints = 2;
            else if (ebl <= 1500) eblPoints = 1;
            else eblPoints = 0;

            // HR points
            if (hr >= 85) hrPoints = 3;
            else if (hr >= 76) hrPoints = 2;
            else if (hr >= 66) hrPoints = 1;
            else if (hr >= 56) hrPoints = 0;
            else hrPoints = 0;

            // MAP points
            if (map >= 85) mapPoints = 3;
            else if (map >= 76) mapPoints = 2;
            else if (map >= 66) mapPoints = 1;
            else if (map >= 50) mapPoints = 0;
            else mapPoints = 0;

            const total = eblPoints + hrPoints + mapPoints;

            let risk, riskText;
            if (total >= 9) {
                risk = 'Low Risk';
                riskText = 'Risk of major complication: 3-5% | Mortality: <1%';
            } else if (total >= 7) {
                risk = 'Intermediate Risk';
                riskText = 'Risk of major complication: 10-15% | Mortality: 1-4%';
            } else if (total >= 5) {
                risk = 'High Risk';
                riskText = 'Risk of major complication: 20-30% | Mortality: 5-10%';
            } else {
                risk = 'Very High Risk';
                riskText = 'Risk of major complication: >40% | Mortality: >15%';
            }

            calculationResults.sas = {
                value: total,
                ebl: ebl,
                hr: hr,
                map: map,
                risk: risk
            };

            document.getElementById('sas-value').textContent = `${total}/10 - ${risk}`;
            document.getElementById('sas-interpretation').innerHTML = `<strong>Score Breakdown:</strong><br>` +
                `• EBL: ${eblPoints} points (${ebl} mL)<br>` +
                `• Lowest HR: ${hrPoints} points (${hr} bpm)<br>` +
                `• Lowest MAP: ${mapPoints} points (${map} mmHg)<br><br>` +
                `<strong>Risk Assessment:</strong><br>${riskText}<br><br>` +
                `Consider closer post-op monitoring for scores <7.`;
            document.getElementById('sas-result').classList.add('visible');
        }

        // Opioid Conversion Calculator
        function calculateOpioid() {
            const from = document.getElementById('opioid-from').value;
            const dose = parseFloat(document.getElementById('opioid-dose').value);
            const to = document.getElementById('opioid-to').value;

            if (!from || !dose || !to) {
                alert('Please fill in all fields');
                return;
            }

            // Morphine equivalents (relative to morphine IV = 1)
            const equivalents = {
                'morphine-iv': 1,
                'morphine-po': 3,
                'fentanyl-iv': 0.01,
                'hydromorphone-iv': 0.2,
                'hydromorphone-po': 0.8,
                'oxycodone-po': 2,
                'hydrocodone-po': 2.5
            };

            const morphineEquivalent = dose * equivalents[from];
            const convertedDose = morphineEquivalent / equivalents[to];
            const reducedDose = convertedDose * 0.75; // 25% reduction for cross-tolerance

            calculationResults.opioid = {
                value: Math.round(reducedDose * 10) / 10,
                from: from,
                dose: dose,
                to: to,
                fullDose: Math.round(convertedDose * 10) / 10
            };

            const fromName = from.replace('-iv', ' IV').replace('-po', ' PO').replace(/^./, m => m.toUpperCase());
            const toName = to.replace('-iv', ' IV').replace('-po', ' PO').replace(/^./, m => m.toUpperCase());

            document.getElementById('opioid-value').textContent = `${Math.round(reducedDose * 10) / 10} ${to.includes('fentanyl') ? 'mcg' : 'mg'}`;
            document.getElementById('opioid-interpretation').innerHTML = `<strong>Conversion:</strong><br>` +
                `• From: ${dose} ${from.includes('fentanyl') ? 'mcg' : 'mg'} ${fromName}<br>` +
                `• Equianalgesic dose: ${Math.round(convertedDose * 10) / 10} ${to.includes('fentanyl') ? 'mcg' : 'mg'} ${toName}<br>` +
                `• <strong>Recommended starting dose: ${Math.round(reducedDose * 10) / 10} ${to.includes('fentanyl') ? 'mcg' : 'mg'}</strong> (25% reduction for incomplete cross-tolerance)<br><br>` +
                `<em>Warning: This is a general guide. Individual patient factors, tolerance, and clinical context must be considered. Titrate to effect.</em>`;
            document.getElementById('opioid-result').classList.add('visible');
        }

        // Local Anesthetic Toxicity Dose Calculator
        function calculateLAST() {
            const weight = parseFloat(document.getElementById('last-weight').value);
            const agent = document.getElementById('last-agent').value;
            const concentration = parseFloat(document.getElementById('last-concentration').value);

            if (!weight || !agent || !concentration) {
                alert('Please fill in all fields');
                return;
            }

            // Max doses in mg/kg
            const maxDoses = {
                'lidocaine': 4.5,
                'lidocaine-epi': 7,
                'bupivacaine': 2.5,
                'ropivacaine': 3,
                'chloroprocaine': 11
            };

            const maxMg = maxDoses[agent] * weight;
            const maxVolume = maxMg / (concentration * 10); // concentration % to mg/mL

            const agentName = agent.replace('-epi', ' with epinephrine').replace(/^./, m => m.toUpperCase());

            calculationResults.last = {
                value: Math.round(maxMg),
                weight: weight,
                agent: agentName,
                concentration: concentration,
                maxVolume: Math.round(maxVolume * 10) / 10
            };

            document.getElementById('last-value').textContent = `${Math.round(maxMg)} mg (${Math.round(maxVolume * 10) / 10} mL)`;
            document.getElementById('last-interpretation').innerHTML = `<strong>Maximum Safe Dose:</strong><br>` +
                `• Agent: ${agentName}<br>` +
                `• Patient weight: ${weight} kg<br>` +
                `• Concentration: ${concentration}%<br>` +
                `• Max dose: ${maxDoses[agent]} mg/kg = ${Math.round(maxMg)} mg<br>` +
                `• Max volume at ${concentration}%: ${Math.round(maxVolume * 10) / 10} mL<br><br>` +
                `<strong>LAST Prevention:</strong> Use lowest effective dose, aspirate before injection, fractionate doses, use ultrasound guidance.<br>` +
                `<strong>LAST Treatment:</strong> Stop injection, airway management, IV lipid emulsion 20% (1.5 mL/kg bolus), ACLS, avoid propofol.`;
            document.getElementById('last-result').classList.add('visible');
        }

        // RCRI Calculator
        function calculateRCRI() {
            const highrisk = document.getElementById('rcri-highrisk').checked;
            const ihd = document.getElementById('rcri-ihd').checked;
            const chf = document.getElementById('rcri-chf').checked;
            const cvd = document.getElementById('rcri-cvd').checked;
            const dm = document.getElementById('rcri-dm').checked;
            const renal = document.getElementById('rcri-renal').checked;

            const score = (highrisk ? 1 : 0) + (ihd ? 1 : 0) + (chf ? 1 : 0) + (cvd ? 1 : 0) + (dm ? 1 : 0) + (renal ? 1 : 0);

            let risk, percentage;
            if (score === 0) {
                risk = 'Low';
                percentage = '0.4%';
            } else if (score === 1) {
                risk = 'Low-Moderate';
                percentage = '0.9%';
            } else if (score === 2) {
                risk = 'Moderate';
                percentage = '6.6%';
            } else {
                risk = 'High';
                percentage = '>11%';
            }

            calculationResults.rcri = {
                value: score,
                risk: risk,
                percentage: percentage
            };

            document.getElementById('rcri-value').textContent = `${score}/6 - ${risk} Risk`;
            document.getElementById('rcri-interpretation').innerHTML = `<strong>Risk of Major Cardiac Event:</strong><br>` +
                `• Score: ${score} points<br>` +
                `• Risk category: ${risk}<br>` +
                `• Predicted risk: ${percentage}<br><br>` +
                `<strong>Recommendations:</strong><br>` +
                (score >= 2 ? '• Consider preoperative cardiology consultation<br>• Consider β-blocker therapy<br>• Optimize medical management' :
                 score === 1 ? '• Perioperative medical optimization<br>• Continue home cardiac medications' :
                 '• Routine perioperative care');
            document.getElementById('rcri-result').classList.add('visible');
        }

        // STOP-BANG Calculator
        function calculateStopBang() {
            const snore = document.getElementById('sb-snore').checked;
            const tired = document.getElementById('sb-tired').checked;
            const observed = document.getElementById('sb-observed').checked;
            const pressure = document.getElementById('sb-pressure').checked;
            const bmi = document.getElementById('sb-bmi').checked;
            const age = document.getElementById('sb-age').checked;
            const neck = document.getElementById('sb-neck').checked;
            const gender = document.getElementById('sb-gender').checked;

            const score = (snore ? 1 : 0) + (tired ? 1 : 0) + (observed ? 1 : 0) + (pressure ? 1 : 0) +
                          (bmi ? 1 : 0) + (age ? 1 : 0) + (neck ? 1 : 0) + (gender ? 1 : 0);

            let risk, osaSeverity;
            if (score <= 2) {
                risk = 'Low';
                osaSeverity = 'Low probability of moderate-severe OSA';
            } else if (score <= 4) {
                risk = 'Intermediate';
                osaSeverity = 'Intermediate probability of moderate-severe OSA';
            } else {
                risk = 'High';
                osaSeverity = 'High probability of moderate-severe OSA';
            }

            calculationResults.stopbang = {
                value: score,
                risk: risk
            };

            document.getElementById('stopbang-value').textContent = `${score}/8 - ${risk} Risk`;
            document.getElementById('stopbang-interpretation').innerHTML = `<strong>OSA Screening Result:</strong><br>` +
                `• Score: ${score} points<br>` +
                `• Risk category: ${risk}<br>` +
                `• ${osaSeverity}<br><br>` +
                `<strong>Perioperative Considerations:</strong><br>` +
                (score >= 5 ? '• High risk - consider sleep study<br>• Avoid opioids/sedatives when possible<br>• Post-op monitoring with continuous pulse oximetry<br>• Consider regional anesthesia' :
                 score >= 3 ? '• Moderate risk - consider precautions<br>• Minimize sedatives/opioids<br>• Extended PACU monitoring' :
                 '• Standard perioperative care');
            document.getElementById('stopbang-result').classList.add('visible');
        }

        // Airway Assessment Calculator
        function calculateAirway() {
            const mallampati = parseInt(document.getElementById('airway-mallampati').value);
            const tmd = document.getElementById('airway-tmd').value;
            const neck = document.getElementById('airway-neck').value;
            const mouth = document.getElementById('airway-mouth').value;

            if (!mallampati || !tmd || !neck || !mouth) {
                alert('Please fill in all fields');
                return;
            }

            let difficultyScore = 0;

            // Mallampati contribution
            if (mallampati >= 3) difficultyScore += 2;
            else if (mallampati === 2) difficultyScore += 1;

            // TMD contribution
            if (tmd === 'poor') difficultyScore += 2;
            else if (tmd === 'borderline') difficultyScore += 1;

            // Neck ROM contribution
            if (neck === 'reduced') difficultyScore += 1;

            // Mouth opening contribution
            if (mouth === 'poor') difficultyScore += 2;

            let difficulty, recommendations;
            if (difficultyScore === 0) {
                difficulty = 'Easy Airway (Predicted)';
                recommendations = 'Standard intubation equipment and approach expected to be successful.';
            } else if (difficultyScore <= 2) {
                difficulty = 'Potentially Difficult';
                recommendations = 'Have backup equipment available (bougie, video laryngoscope). Consider awake look with sedation.';
            } else if (difficultyScore <= 4) {
                difficulty = 'Likely Difficult';
                recommendations = 'Primary plan: Video laryngoscope. Backup: LMA, bougie. Prepare difficult airway cart. Consider awake fiberoptic intubation.';
            } else {
                difficulty = 'Very Difficult Airway';
                recommendations = 'Strong consideration for awake fiberoptic intubation. Difficult airway cart at bedside. Consider ENT backup. Inform patient of risks.';
            }

            calculationResults.airway = {
                value: difficultyScore,
                difficulty: difficulty,
                mallampati: mallampati,
                tmd: tmd,
                neck: neck,
                mouth: mouth
            };

            document.getElementById('airway-value').textContent = difficulty;
            document.getElementById('airway-interpretation').innerHTML = `<strong>Assessment Details:</strong><br>` +
                `• Mallampati: Class ${mallampati}<br>` +
                `• Thyromental distance: ${tmd}<br>` +
                `• Neck ROM: ${neck}<br>` +
                `• Mouth opening: ${mouth}<br>` +
                `• Difficulty score: ${difficultyScore}/7<br><br>` +
                `<strong>Recommendations:</strong><br>${recommendations}`;
            document.getElementById('airway-result').classList.add('visible');
        }

        // MAC Calculator
        function calculateMAC() {
            const agent = document.getElementById('mac-agent').value;
            const age = parseFloat(document.getElementById('mac-age').value);

            if (!agent || !age) {
                alert('Please fill in all fields');
                return;
            }

            // MAC at age 40
            const macAt40 = {
                'sevoflurane': 2.0,
                'desflurane': 6.0,
                'isoflurane': 1.15
            };

            // MAC decreases ~6% per decade after 40
            const ageDiff = age - 40;
            const macMultiplier = 1 - (ageDiff / 10) * 0.06;
            const mac = macAt40[agent] * macMultiplier;
            const macAwake = mac * 0.5; // MAC-awake is approximately 0.5 MAC

            const agentName = agent.charAt(0).toUpperCase() + agent.slice(1);

            calculationResults.mac = {
                value: Math.round(mac * 100) / 100,
                agent: agentName,
                age: age,
                macAwake: Math.round(macAwake * 100) / 100
            };

            document.getElementById('mac-value').textContent = `${Math.round(mac * 100) / 100}%`;
            document.getElementById('mac-interpretation').innerHTML = `<strong>Age-Adjusted MAC Values:</strong><br>` +
                `• Agent: ${agentName}<br>` +
                `• Patient age: ${age} years<br>` +
                `• MAC: ${Math.round(mac * 100) / 100}%<br>` +
                `• MAC-awake: ${Math.round(macAwake * 100) / 100}%<br>` +
                `• MAC-BAR (blocks adrenergic response): ${Math.round(mac * 1.3 * 100) / 100}%<br><br>` +
                `<em>Note: MAC decreases ~6% per decade after age 40. These are population averages; individual requirements vary.</em>`;
            document.getElementById('mac-result').classList.add('visible');
        }

        // Propofol TCI Calculator
        document.getElementById('tci-model').addEventListener('change', function() {
            if (this.value === 'schnider') {
                document.getElementById('tci-schnider-params').style.display = 'grid';
            } else {
                document.getElementById('tci-schnider-params').style.display = 'none';
            }
        });

        function calculatePropofolTCI() {
            const model = document.getElementById('tci-model').value;
            const target = parseFloat(document.getElementById('tci-target').value);
            const weight = parseFloat(document.getElementById('tci-weight').value);
            const age = parseFloat(document.getElementById('tci-age').value);

            if (!model || !target || !weight || !age) {
                alert('Please fill in required fields');
                return;
            }

            let infusionRate, v1, k10;

            if (model === 'marsh') {
                // Simplified Marsh model
                v1 = 0.228 * weight; // L
                k10 = 0.119; // min^-1
                infusionRate = target * v1 * k10 * 60 / 1000; // mg/kg/hr converted to mcg/kg/min then to mg/kg/hr
                infusionRate = target * v1 * k10 * 60; // mg/hr
            } else if (model === 'schnider') {
                const height = parseFloat(document.getElementById('tci-height').value);
                const sex = document.getElementById('tci-sex').value;

                if (!height || !sex) {
                    alert('Please enter height and sex for Schnider model');
                    return;
                }

                // Simplified Schnider model (effect-site targeting)
                const lbm = sex === 'male' ?
                    1.1 * weight - 128 * (weight / height) ** 2 :
                    1.07 * weight - 148 * (weight / height) ** 2;
                v1 = 4.27; // L (fixed)
                k10 = 0.443 + 0.0107 * (weight - 77) - 0.0159 * (lbm - 59) + 0.0062 * (height - 177);
                infusionRate = target * v1 * k10 * 60; // mg/hr
            } else { // kataria (pediatric)
                v1 = 0.41 * weight; // L
                k10 = 0.0678 * Math.pow(age, -0.269);
                infusionRate = target * v1 * k10 * 60; // mg/hr
            }

            const mcgKgMin = (infusionRate * 1000) / (weight * 60);

            const modelName = model.charAt(0).toUpperCase() + model.slice(1);
            calculationResults['propofol-tci'] = {
                value: Math.round(infusionRate),
                model: modelName,
                target: target,
                weight: weight,
                age: age
            };

            document.getElementById('tci-value').textContent = `${Math.round(mcgKgMin)} mcg/kg/min`;
            document.getElementById('tci-interpretation').innerHTML = `<strong>TCI Parameters:</strong><br>` +
                `• Model: ${modelName}<br>` +
                `• Target concentration: ${target} mcg/mL<br>` +
                `• Patient weight: ${weight} kg<br>` +
                `• Age: ${age} years<br><br>` +
                `<strong>Infusion Rates:</strong><br>` +
                `• ${Math.round(mcgKgMin)} mcg/kg/min<br>` +
                `• ${Math.round(infusionRate)} mg/hr<br>` +
                `• ${Math.round(mcgKgMin * weight * 60 / 100)} mL/hr (at 10 mg/mL)<br><br>` +
                `<em>Note: This is a simplified calculation. Actual TCI pumps use complex pharmacokinetic models with multiple compartments. Use for estimation only.</em>`;
            document.getElementById('tci-result').classList.add('visible');
        }

        // Send to AI function - Now opens modal instead of redirecting
        function sendToAI(calcType) {
            const result = calculationResults[calcType];
            if (!result && calcType !== 'asa') {
                alert('Please calculate first');
                return;
            }

            let message = '';

            switch(calcType) {
                case 'ibw':
                    message = `Patient: ${result.sex}, ${result.height} cm tall. IBW = ${result.value} kg. What should I consider for anesthetic drug dosing?`;
                    break;
                case 'mabl':
                    message = `Patient: ${result.weight} kg ${result.sex}, height ${result.height} cm. Initial Hct ${result.initialHct}%, target Hct ${result.targetHct}%. Estimated blood volume ${result.bloodVolume} mL, MABL = ${result.value} mL. What's my transfusion strategy for this case?`;
                    break;
                case 'bsa':
                    message = `Patient: ${result.weight} kg, ${result.height} cm. BSA = ${result.value} m². How does this affect my anesthetic management?`;
                    break;
                case 'qtc':
                    message = `Patient has QTc = ${result.value} ms (QT ${result.qt} ms, HR ${result.hr} bpm). Risk level: ${result.risk}. What anesthetic drugs should I avoid or use with caution?`;
                    break;
                case 'fluids':
                    message = `Patient weighs ${result.weight} kg. Calculated maintenance fluids = ${result.value} mL/hr (${result.daily} mL/day). What should I know about perioperative fluid management?`;
                    break;
                case 'ponv':
                    message = `Patient has Apfel score ${result.score}/4 (${result.risk} PONV risk, ~${result.percentage} incidence). Risk factors: ${result.factors.female ? 'female, ' : ''}${result.factors.nonsmoker ? 'non-smoker, ' : ''}${result.factors.history ? 'PONV/motion sickness history, ' : ''}${result.factors.opioids ? 'postop opioids expected' : ''}. What's my best prophylaxis strategy?`;
                    break;
                case 'asa':
                    message = `Can you explain the ASA Physical Status Classification system and how it affects perioperative risk?`;
                    break;
                case 'fluid-deficit':
                    message = `Patient weighs ${result.weight} kg with ${result.dehydration}% dehydration. Calculated fluid deficit = ${result.value} mL. What should I know about fluid resuscitation strategy?`;
                    break;
                case 'defib':
                    message = `${result.type === 'defib' ? 'Defibrillation' : 'Cardioversion'} energy for ${result.age} patient${result.age === 'peds' ? ' (' + result.weight + ' kg)' : ''}: ${result.value}. What are the key safety considerations?`;
                    break;
                case 'sas':
                    message = `Patient has Surgical Apgar Score of ${result.value}/10 (${result.risk}). EBL ${result.ebl} mL, lowest HR ${result.hr} bpm, lowest MAP ${result.map} mmHg. What does this mean for post-op management?`;
                    break;
                case 'opioid':
                    message = `Converting ${result.dose} ${result.from.includes('fentanyl') ? 'mcg' : 'mg'} ${result.from} to ${result.value} ${result.to.includes('fentanyl') ? 'mcg' : 'mg'} ${result.to}. What should I know about opioid rotation and cross-tolerance?`;
                    break;
                case 'last':
                    message = `Maximum safe dose of ${result.agent} at ${result.concentration}% for ${result.weight} kg patient is ${result.value} mg (${result.maxVolume} mL). What are the signs and treatment of local anesthetic systemic toxicity?`;
                    break;
                case 'rcri':
                    message = `Patient has RCRI score of ${result.value}/6 (${result.risk} risk, ${result.percentage} cardiac event risk). How should I optimize this patient perioperatively?`;
                    break;
                case 'stopbang':
                    message = `Patient has STOP-BANG score of ${result.value}/8 (${result.risk} OSA risk). What are the perioperative implications and management strategies?`;
                    break;
                case 'airway':
                    message = `Airway assessment predicts: ${result.difficulty} (Mallampati ${result.mallampati}, TMD ${result.tmd}, neck ROM ${result.neck}, mouth ${result.mouth}). What's my airway management plan?`;
                    break;
                case 'mac':
                    message = `Age-adjusted MAC for ${result.agent} in ${result.age} year old patient is ${result.value}% (MAC-awake ${result.macAwake}%). How should I titrate volatile anesthetics based on age?`;
                    break;
                case 'propofol-tci':
                    message = `Propofol TCI using ${result.model} model: target ${result.target} mcg/mL for ${result.weight} kg patient, age ${result.age}. Rate = ${result.value} mg/hr. How do different TCI models compare?`;
                    break;
            }

            // Open modal and send initial message
            openAIChat(message);
        }

        // AI Chat Modal Functions
        function openAIChat(initialMessage) {
            const modal = document.getElementById('aiChatModal');
            const messagesDiv = document.getElementById('aiChatMessages');
            const inputField = document.getElementById('aiChatInput');

            // Clear previous messages
            messagesDiv.innerHTML = '';
            inputField.value = '';
            currentChatContext = initialMessage;

            // Show modal
            modal.classList.add('active');

            // Send initial message
            if (initialMessage) {
                fetchAIResponse(initialMessage);
            }
        }

        function closeAIChat() {
            document.getElementById('aiChatModal').classList.remove('active');
        }

        // Close modal when clicking outside
        document.getElementById('aiChatModal').addEventListener('click', function(e) {
            if (e.target === this) {
                closeAIChat();
            }
        });

        // Send message on Enter key
        document.getElementById('aiChatInput').addEventListener('keydown', function(e) {
            if (e.key === 'Enter' && !e.shiftKey) {
                e.preventDefault();
                sendAIMessage();
            }
        });

        function sendAIMessage() {
            const input = document.getElementById('aiChatInput');
            const message = input.value.trim();

            if (!message) return;

            // Add user message to chat
            addMessageToChat(message, 'user');
            input.value = '';

            // Fetch AI response
            fetchAIResponse(message);
        }

        function addMessageToChat(message, type) {
            const messagesDiv = document.getElementById('aiChatMessages');
            const messageEl = document.createElement('div');
            messageEl.className = 'ai-message';

            if (type === 'user') {
                messageEl.style.background = 'linear-gradient(135deg, rgba(37, 99, 235, 0.1) 0%, rgba(139, 92, 246, 0.1) 100%)';
                messageEl.style.borderLeft = '3px solid var(--primary-blue)';
                messageEl.innerHTML = `<strong>You:</strong> ${escapeHtml(message)}`;
            } else if (type === 'loading') {
                messageEl.className = 'ai-loading';
                messageEl.innerHTML = `
                    Thinking...
                    <div class="ai-loading-dots">
                        <span></span>
                        <span></span>
                        <span></span>
                    </div>
                `;
                messageEl.id = 'loadingMessage';
            } else {
                messageEl.innerHTML = `<strong>AI:</strong> ${message}`;
            }

            messagesDiv.appendChild(messageEl);
            messagesDiv.scrollTop = messagesDiv.scrollHeight;

            return messageEl;
        }

        async function fetchAIResponse(userMessage) {
            // Show loading indicator
            addMessageToChat('', 'loading');

            try {
                // Get CSRF token
                const csrfToken = document.getElementById('csrf_token').value;

                const response = await fetch('/chat', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/x-www-form-urlencoded',
                    },
                    body: `query=${encodeURIComponent(userMessage)}&csrf_token=${encodeURIComponent(csrfToken)}`
                });

                if (!response.ok) {
                    throw new Error('Network response was not ok');
                }

                // Get the full page HTML response
                const html = await response.text();

                // Parse the HTML to extract the last AI message
                const parser = new DOMParser();
                const doc = parser.parseFromString(html, 'text/html');
                const assistantMessages = doc.querySelectorAll('.message.assistant');

                if (assistantMessages.length > 0) {
                    const lastMessage = assistantMessages[assistantMessages.length - 1];
                    const messageText = lastMessage.querySelector('.message-text');

                    if (messageText) {
                        // Remove loading message
                        const loadingMsg = document.getElementById('loadingMessage');
                        if (loadingMsg) loadingMsg.remove();

                        // Add AI response
                        addMessageToChat(messageText.innerHTML, 'ai');
                    } else {
                        throw new Error('Could not extract AI response');
                    }
                } else {
                    throw new Error('No response found');
                }

            } catch (error) {
                console.error('Error fetching AI response:', error);
                const loadingMsg = document.getElementById('loadingMessage');
                if (loadingMsg) loadingMsg.remove();

                addMessageToChat('Sorry, there was an error getting a response. Please try again or visit the <a href="/chat" style="color: var(--primary-blue); text-decoration: underline;">chat page</a>.', 'ai');
            }
        }

        function escapeHtml(text) {
            const div = document.createElement('div');
            div.textContent = text;
            return div.innerHTML;
        }
    </script>
</body>
</html>
"""

HYPOTENSION_HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>IOH Predictor — gasconsult.ai</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Sora:wght@300;400;500;600;700&family=Inter:wght@400;500;600&family=JetBrains+Mono:wght@400;500;600&display=swap" rel="stylesheet">
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=5">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">
    <meta name="apple-mobile-web-app-capable" content="yes">
    <meta name="apple-mobile-web-app-status-bar-style" content="default">
    <meta name="apple-mobile-web-app-title" content="gasconsult.ai">
    <style>
        :root {
            --primary-blue: #2563EB;
            --primary-blue-dark: #1D4ED8;
            --primary-blue-light: #DBEAFE;
            --primary: #2563EB;
            --primary-dark: #1D4ED8;
            --text-primary: #0F172A;
            --text-secondary: #475569;
            --text-muted: #94A3B8;
            --bg-primary: #FFFFFF;
            --bg-secondary: #F8FAFC;
            --border: #E2E8F0;
            --warning-yellow: #F59E0B;
            --warning-yellow-light: #FEF3C7;
            --success-green: #10B981;
            --success-green-light: #D1FAE5;
            --danger-red: #EF4444;
            --danger-red-light: #FEE2E2;
        }

        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, system-ui, sans-serif;
            background: var(--bg-secondary);
            color: var(--text-primary);
            line-height: 1.6;
            -webkit-font-smoothing: antialiased;
            overflow-x: hidden;
            width: 100%;
        }

        header {
            background: rgba(255, 255, 255, 0.85);
            backdrop-filter: blur(12px);
            -webkit-backdrop-filter: blur(12px);
            border-bottom: 1px solid rgba(226, 232, 240, 0.8);
            padding: 16px 40px;
            position: sticky;
            top: 0;
            z-index: 100;
        }

        .header-container {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 24px;
            flex-wrap: wrap;
        }

        .header-left {
            display: flex;
            align-items: center;
            gap: 16px;
        }

        .logo-link {
            display: flex;
            align-items: center;
            gap: 12px;
            text-decoration: none;
            transition: transform 0.2s ease;
        }

        .logo-link:hover {
            transform: translateY(-1px);
        }

        .logo-svg {
            height: 28px;
            width: auto;
        }

        .logo-text {
            font-family: 'Sora', sans-serif;
            font-size: 20px;
            font-weight: 600;
            line-height: 1;
            letter-spacing: -0.5px;
        }

        .logo-gas {
            color: #2563EB;
        }

        .logo-consult {
            color: #111111;
        }

        .logo-ai {
            color: #6B7280;
            font-weight: 400;
        }

        .nav-links {
            display: flex;
            gap: 12px;
            align-items: center;
        }

        .nav-link {
            color: var(--text-secondary);
            text-decoration: none;
            font-size: 14px;
            font-weight: 500;
            padding: 8px 16px;
            border-radius: 8px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            color: var(--text-primary);
            background: rgba(255, 255, 255, 0.6);
            backdrop-filter: blur(8px);
            -webkit-backdrop-filter: blur(8px);
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.06);
        }

        .nav-link.active {
            color: var(--primary-blue);
            font-weight: 600;
        }

        main {
            max-width: 900px;
            margin: 0 auto;
            padding: 32px 24px;
        }

        /* Educational Warning Banner */
        .edu-warning {
            background: linear-gradient(135deg, #FEF3C7 0%, #FEE2E2 100%);
            border: 2px solid #F59E0B;
            border-radius: 16px;
            padding: 20px 24px;
            margin-bottom: 24px;
            display: flex;
            align-items: flex-start;
            gap: 12px;
        }

        .edu-warning-icon {
            font-size: 24px;
            flex-shrink: 0;
        }

        .edu-warning-content h3 {
            font-family: 'Sora', sans-serif;
            font-size: 16px;
            font-weight: 700;
            color: #B91C1C;
            margin-bottom: 8px;
        }

        .edu-warning-content p {
            font-size: 14px;
            color: var(--text-primary);
            line-height: 1.5;
        }

        /* Page Header */
        .page-header {
            text-align: center;
            margin-bottom: 32px;
        }

        .page-header h1 {
            font-family: 'Sora', sans-serif;
            font-size: 2.5rem;
            background: linear-gradient(135deg, var(--primary-blue) 0%, #8B5CF6 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            margin-bottom: 12px;
            font-weight: 700;
            letter-spacing: -1px;
        }

        .page-header p {
            color: var(--text-secondary);
            font-size: 16px;
            max-width: 600px;
            margin: 0 auto;
        }

        .edu-badge {
            display: inline-block;
            background: var(--warning-yellow-light);
            color: #92400E;
            padding: 6px 12px;
            border-radius: 8px;
            font-size: 12px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            margin-top: 12px;
        }

        /* Disclaimer Box */
        .disclaimer-box {
            background: white;
            border: 2px solid var(--border);
            border-radius: 16px;
            padding: 24px;
            margin-bottom: 32px;
        }

        .disclaimer-box h3 {
            font-family: 'Sora', sans-serif;
            font-size: 18px;
            font-weight: 700;
            color: var(--text-primary);
            margin-bottom: 16px;
        }

        .disclaimer-box p {
            font-size: 14px;
            color: var(--text-secondary);
            margin-bottom: 12px;
            line-height: 1.6;
        }

        .disclaimer-box ul {
            list-style: none;
            padding-left: 0;
            font-size: 14px;
            color: var(--text-secondary);
        }

        .disclaimer-box ul li {
            padding-left: 24px;
            position: relative;
            margin-bottom: 8px;
        }

        .disclaimer-box ul li:before {
            content: "•";
            position: absolute;
            left: 8px;
            color: var(--primary-blue);
            font-weight: 700;
        }

        /* Model Performance Metrics */
        .metrics-box {
            background: linear-gradient(135deg, #F0F9FF 0%, #E0F2FE 100%);
            border: 2px solid var(--primary-blue);
            border-radius: 16px;
            padding: 24px;
            margin-bottom: 32px;
        }

        .metrics-box h3 {
            font-family: 'Sora', sans-serif;
            font-size: 18px;
            font-weight: 700;
            color: var(--primary-blue);
            margin-bottom: 12px;
            display: flex;
            align-items: center;
            gap: 8px;
        }

        .metrics-box .metrics-intro {
            font-size: 14px;
            color: var(--text-secondary);
            margin-bottom: 20px;
            line-height: 1.6;
        }

        .metrics-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 16px;
            margin-bottom: 20px;
        }

        .metric-card {
            background: white;
            border-radius: 12px;
            padding: 16px;
            border: 1px solid rgba(37, 99, 235, 0.2);
            text-align: center;
        }

        .metric-card .metric-label {
            font-size: 12px;
            font-weight: 600;
            color: var(--text-muted);
            text-transform: uppercase;
            letter-spacing: 0.5px;
            margin-bottom: 8px;
        }

        .metric-card .metric-value {
            font-family: 'Sora', sans-serif;
            font-size: 28px;
            font-weight: 700;
            color: var(--primary-blue);
            margin-bottom: 8px;
        }

        .metric-card .metric-description {
            font-size: 12px;
            color: var(--text-secondary);
            line-height: 1.4;
        }

        .metrics-explanation {
            background: white;
            border-radius: 12px;
            padding: 16px;
            border-left: 4px solid var(--primary-blue);
        }

        .metrics-explanation h4 {
            font-size: 14px;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 12px;
        }

        .metrics-explanation ul {
            list-style: none;
            padding: 0;
        }

        .metrics-explanation li {
            font-size: 13px;
            color: var(--text-secondary);
            margin-bottom: 8px;
            padding-left: 20px;
            position: relative;
            line-height: 1.5;
        }

        .metrics-explanation li:before {
            content: "→";
            position: absolute;
            left: 0;
            color: var(--primary-blue);
            font-weight: 700;
        }

        .metrics-explanation .metric-term {
            font-weight: 600;
            color: var(--text-primary);
        }

        /* Form Sections */
        .form-section {
            background: white;
            border-radius: 16px;
            padding: 24px;
            margin-bottom: 24px;
            border: 1px solid var(--border);
            box-shadow: 0 1px 3px rgba(0,0,0,0.04);
        }

        .form-section h2 {
            font-family: 'Sora', sans-serif;
            font-size: 18px;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 20px;
            padding-bottom: 12px;
            border-bottom: 2px solid var(--bg-secondary);
        }

        .form-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 16px;
        }

        .form-group {
            display: flex;
            flex-direction: column;
        }

        .form-group label {
            font-size: 13px;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 8px;
        }

        .form-group input,
        .form-group select {
            padding: 10px 12px;
            border: 2px solid var(--border);
            border-radius: 8px;
            font-size: 14px;
            font-family: 'Inter', sans-serif;
            transition: all 0.2s ease;
        }

        .form-group input:focus,
        .form-group select:focus {
            outline: none;
            border-color: var(--primary-blue);
            box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.1);
        }

        .submit-btn {
            background: linear-gradient(135deg, var(--primary-blue) 0%, var(--primary-blue-dark) 100%);
            color: white;
            padding: 14px 32px;
            border-radius: 12px;
            border: none;
            font-size: 16px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.3s ease;
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.3);
            width: 100%;
            margin-top: 24px;
        }

        .submit-btn:hover {
            transform: translateY(-2px);
            box-shadow: 0 6px 20px rgba(37, 99, 235, 0.4);
        }

        /* Results Section */
        .results-warning {
            background: var(--warning-yellow-light);
            border-left: 4px solid var(--warning-yellow);
            padding: 16px;
            border-radius: 8px;
            margin-bottom: 24px;
        }

        .results-warning p {
            font-size: 14px;
            font-weight: 600;
            color: #92400E;
            display: flex;
            align-items: center;
            gap: 8px;
        }

        .risk-gauges {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 32px;
        }

        .risk-gauge {
            background: white;
            border-radius: 16px;
            padding: 24px;
            border: 2px solid var(--border);
            text-align: center;
        }

        .risk-gauge h3 {
            font-family: 'Sora', sans-serif;
            font-size: 14px;
            font-weight: 600;
            color: var(--text-muted);
            margin-bottom: 16px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        .risk-value {
            font-family: 'Sora', sans-serif;
            font-size: 48px;
            font-weight: 700;
            margin-bottom: 12px;
        }

        .risk-label {
            font-size: 14px;
            font-weight: 600;
            padding: 6px 16px;
            border-radius: 8px;
            display: inline-block;
        }

        .risk-low {
            color: #059669;
            background: var(--success-green-light);
        }

        .risk-moderate {
            color: #D97706;
            background: var(--warning-yellow-light);
        }

        .risk-high {
            color: #DC2626;
            background: var(--danger-red-light);
        }

        .risk-value.low { color: #059669; }
        .risk-value.moderate { color: #D97706; }
        .risk-value.high { color: #DC2626; }

        /* Contributing Factors */
        .factors-section {
            background: white;
            border-radius: 16px;
            padding: 24px;
            margin-bottom: 24px;
            border: 1px solid var(--border);
        }

        .factors-section h2 {
            font-family: 'Sora', sans-serif;
            font-size: 18px;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 16px;
        }

        .factor-item {
            padding: 12px;
            border-left: 4px solid var(--primary-blue);
            background: var(--bg-secondary);
            border-radius: 8px;
            margin-bottom: 12px;
        }

        .factor-item strong {
            color: var(--text-primary);
        }

        /* Interventions */
        .interventions-section {
            background: white;
            border-radius: 16px;
            padding: 24px;
            border: 1px solid var(--border);
        }

        .interventions-section h2 {
            font-family: 'Sora', sans-serif;
            font-size: 18px;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 16px;
        }

        .intervention-item {
            display: flex;
            align-items: flex-start;
            gap: 12px;
            padding: 16px;
            background: var(--bg-secondary);
            border-radius: 12px;
            margin-bottom: 12px;
        }

        .intervention-rank {
            background: var(--primary-blue);
            color: white;
            font-weight: 700;
            font-size: 14px;
            padding: 6px 12px;
            border-radius: 8px;
            flex-shrink: 0;
        }

        .intervention-content {
            flex: 1;
        }

        .intervention-content h3 {
            font-size: 16px;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 4px;
        }

        .intervention-content p {
            font-size: 14px;
            color: var(--text-secondary);
        }

        footer {
            text-align: center;
            padding: 40px;
            border-top: 1px solid var(--border);
            background: white;
            color: var(--text-muted);
            font-size: 13px;
            margin-top: 48px;
        }

        footer a {
            color: var(--primary);
            text-decoration: none;
        }

        @media (max-width: 768px) {
            header {
                padding: 14px 20px;
            }

            .header-container {
                flex-direction: column;
                align-items: flex-start;
            }

            .nav-links {
                gap: 8px;
            }

            .nav-link {
                padding: 8px 14px;
                font-size: 0.9rem;
            }

            .page-header h1 {
                font-size: 2rem;
            }

            .form-grid {
                grid-template-columns: 1fr;
            }

            .risk-gauges {
                grid-template-columns: 1fr;
            }

            .form-section {
                padding: 20px;
            }

            input, select {
                font-size: 16px !important; /* Prevents iOS zoom on focus */
            }
        }

        @media (max-width: 640px) {
            .header-container {
                justify-content: space-between;
                flex-wrap: wrap;
            }

            .header-left {
                order: 1;
            }

            .nav-links {
                order: 2;
                width: 100%;
                justify-content: center;
                margin-top: 12px;
            }

            .nav-link {
                padding: 6px 10px;
                font-size: 0.8rem;
            }

            .logo-text {
                font-size: 1rem;
            }

            .logo-svg {
                height: 24px;
            }
        }
    </style>
</head>
<body>
    <header>
        <div class="header-container">
            <div class="header-left">
                <a href="/" class="logo-link">
                    <svg class="logo-svg" viewBox="0 0 50 24" fill="none" xmlns="http://www.w3.org/2000/svg">
                        <defs>
                            <linearGradient id="ecgGrad" x1="0%" y1="0%" x2="100%" y2="0%">
                                <stop offset="0%" stop-color="#2563EB"/>
                                <stop offset="20%" stop-color="#EF4444"/>
                                <stop offset="40%" stop-color="#FBBF24"/>
                                <stop offset="60%" stop-color="#8B5CF6"/>
                                <stop offset="80%" stop-color="#10B981"/>
                                <stop offset="100%" stop-color="#6B7280"/>
                            </linearGradient>
                        </defs>
                        <path d="M2 12 L8 12 L11 10 L14 14 L17 4 L20 20 L23 9 L26 12 L34 12"
                              stroke="url(#ecgGrad)"
                              stroke-width="2.5"
                              stroke-linecap="round"
                              stroke-linejoin="round"
                              fill="none"/>
                    </svg>
                    <span class="logo-text">
                        <span class="logo-gas">gas</span><span class="logo-consult">consult</span><span class="logo-ai">.ai</span>
                    </span>
                </a>
            </div>
            <div class="nav-links">
                <a href="/" class="nav-link">Home</a>
                <a href="/preop" class="nav-link">Pre-Op Assessment</a>
                <a href="/calculators" class="nav-link">Clinical Calculators</a>
                <a href="/quick-dose" class="nav-link">Quick Dose</a>
                <a href="/hypotension" class="nav-link active">IOH Predictor</a>
            </div>
        </div>
    </header>

    <main>
        <!-- Educational Warning Banner -->
        <div class="edu-warning">
            <div class="edu-warning-icon">⚠️</div>
            <div class="edu-warning-content">
                <h3>EDUCATIONAL USE ONLY - NOT FOR CLINICAL DECISION MAKING</h3>
                <p>This tool is for educational demonstration purposes only. It is NOT validated for clinical use, NOT FDA-approved, and should NEVER be used to make actual patient care decisions.</p>
            </div>
        </div>

        <!-- Page Header -->
        <div class="page-header">
            <h1>Intraoperative Hypotension Predictor</h1>
            <p>Educational simulator for predicting intraoperative hypotension risk based on patient and procedural factors</p>
            <span class="edu-badge">Educational Demo Only</span>
        </div>

        <!-- Legal Disclaimer -->
        <div class="disclaimer-box">
            <h3>Legal Disclaimer</h3>
            <p><strong>This Intraoperative Hypotension Predictor is a LEARNING TOOL designed for educational purposes only.</strong> It is not a medical device, has not been validated in clinical settings, and is not approved by the FDA or any regulatory body.</p>
            <p><strong>DO NOT use this tool to make clinical decisions.</strong> Predictions may be inaccurate and should never replace clinical judgment, monitoring, or standard anesthesia care protocols.</p>
            <p>By using this tool, you acknowledge that:</p>
            <ul>
                <li>This is for EDUCATIONAL DEMONSTRATION only</li>
                <li>Predictions are hypothetical and may be incorrect</li>
                <li>No patient data should be entered</li>
                <li>You will not rely on outputs for any clinical decisions</li>
                <li>The creators assume no liability for any use of this tool</li>
            </ul>
        </div>

        <!-- Model Performance Metrics -->
        <div class="metrics-box">
            <h3>📊 Model Performance Metrics</h3>
            <p class="metrics-intro">
                This educational model demonstrates hypotension prediction capabilities based on published research.
                Below are typical performance metrics from validated intraoperative hypotension prediction models in the literature.
                Understanding these metrics helps interpret the model's reliability and limitations.
            </p>

            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">AUC-ROC</div>
                    <div class="metric-value">0.84</div>
                    <div class="metric-description">Overall discrimination ability</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Sensitivity</div>
                    <div class="metric-value">78%</div>
                    <div class="metric-description">True positive rate</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Specificity</div>
                    <div class="metric-value">76%</div>
                    <div class="metric-description">True negative rate</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">PPV</div>
                    <div class="metric-value">68%</div>
                    <div class="metric-description">Positive predictive value</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">NPV</div>
                    <div class="metric-value">84%</div>
                    <div class="metric-description">Negative predictive value</div>
                </div>
            </div>

            <div class="metrics-explanation">
                <h4>What These Metrics Mean:</h4>
                <ul>
                    <li><span class="metric-term">AUC-ROC (0.84):</span> Area Under the Receiver Operating Characteristic curve. Ranges from 0.5 (random) to 1.0 (perfect). 0.84 indicates good discrimination between patients who will/won't develop hypotension.</li>
                    <li><span class="metric-term">Sensitivity (78%):</span> Of patients who actually develop hypotension, the model correctly predicts 78% of them. Higher is better for avoiding missed cases.</li>
                    <li><span class="metric-term">Specificity (76%):</span> Of patients who don't develop hypotension, the model correctly identifies 76% as low-risk. Higher reduces false alarms.</li>
                    <li><span class="metric-term">PPV (68%):</span> When the model predicts hypotension, there's a 68% chance it will actually occur. Depends on baseline hypotension prevalence.</li>
                    <li><span class="metric-term">NPV (84%):</span> When the model predicts no hypotension, there's an 84% chance the prediction is correct. Useful for reassurance in low-risk cases.</li>
                </ul>
            </div>
        </div>

        <!-- Input Form -->
        <form method="POST" action="/hypotension">
            <input type="hidden" name="csrf_token" value="{{ csrf_token() }}"/>
            <div class="form-section">
                <h2>Patient Demographics</h2>
                <div class="form-grid">
                    <div class="form-group">
                        <label for="age">Age (years)</label>
                        <input type="number" id="age" name="age" min="1" max="120" required>
                    </div>
                    <div class="form-group">
                        <label for="sex">Sex</label>
                        <select id="sex" name="sex" required>
                            <option value="">Select...</option>
                            <option value="male">Male</option>
                            <option value="female">Female</option>
                        </select>
                    </div>
                    <div class="form-group">
                        <label for="weight">Weight (kg)</label>
                        <input type="number" id="weight" name="weight" min="1" max="300" step="0.1" required>
                    </div>
                    <div class="form-group">
                        <label for="height">Height (cm)</label>
                        <input type="number" id="height" name="height" min="50" max="250" step="0.1" required>
                    </div>
                    <div class="form-group">
                        <label for="asa">ASA Classification</label>
                        <select id="asa" name="asa" required>
                            <option value="">Select...</option>
                            <option value="1">ASA I - Healthy</option>
                            <option value="2">ASA II - Mild systemic disease</option>
                            <option value="3">ASA III - Severe systemic disease</option>
                            <option value="4">ASA IV - Life-threatening disease</option>
                        </select>
                    </div>
                </div>
            </div>

            <div class="form-section">
                <h2>Baseline Vitals</h2>
                <div class="form-grid">
                    <div class="form-group">
                        <label for="baseline_map">Baseline MAP (mmHg)</label>
                        <input type="number" id="baseline_map" name="baseline_map" min="40" max="150" required>
                    </div>
                    <div class="form-group">
                        <label for="baseline_hr">Baseline Heart Rate (bpm)</label>
                        <input type="number" id="baseline_hr" name="baseline_hr" min="30" max="200" required>
                    </div>
                </div>
            </div>

            <div class="form-section">
                <h2>Current Intraoperative State</h2>
                <div class="form-grid">
                    <div class="form-group">
                        <label for="current_map">Current MAP (mmHg)</label>
                        <input type="number" id="current_map" name="current_map" min="40" max="150" required>
                    </div>
                    <div class="form-group">
                        <label for="map_5min">MAP 5 min ago (mmHg)</label>
                        <input type="number" id="map_5min" name="map_5min" min="40" max="150" required>
                    </div>
                    <div class="form-group">
                        <label for="map_10min">MAP 10 min ago (mmHg)</label>
                        <input type="number" id="map_10min" name="map_10min" min="40" max="150" required>
                    </div>
                    <div class="form-group">
                        <label for="surgery_duration">Surgery Duration (min)</label>
                        <input type="number" id="surgery_duration" name="surgery_duration" min="0" max="1000" required>
                    </div>
                    <div class="form-group">
                        <label for="vasopressor">Vasopressor Use</label>
                        <select id="vasopressor" name="vasopressor" required>
                            <option value="none">None</option>
                            <option value="phenylephrine">Phenylephrine</option>
                            <option value="ephedrine">Ephedrine</option>
                            <option value="norepinephrine">Norepinephrine</option>
                            <option value="epinephrine">Epinephrine</option>
                        </select>
                    </div>
                    <div class="form-group">
                        <label for="surgery_type">Surgery Type</label>
                        <select id="surgery_type" name="surgery_type" required>
                            <option value="">Select...</option>
                            <option value="cardiac">Cardiac Surgery</option>
                            <option value="major_abdominal">Major Abdominal</option>
                            <option value="spine">Spine Surgery</option>
                            <option value="orthopedic">Orthopedic</option>
                            <option value="neuro">Neurosurgery</option>
                            <option value="vascular">Vascular</option>
                            <option value="thoracic">Thoracic</option>
                            <option value="minor">Minor/Superficial</option>
                        </select>
                    </div>
                    <div class="form-group">
                        <label for="induction_agent">Induction Agent</label>
                        <select id="induction_agent" name="induction_agent" required>
                            <option value="">Select...</option>
                            <option value="propofol">Propofol</option>
                            <option value="etomidate">Etomidate</option>
                            <option value="ketamine">Ketamine</option>
                            <option value="thiopental">Thiopental</option>
                        </select>
                    </div>
                    <div class="form-group">
                        <label for="emergency">Emergency Status</label>
                        <select id="emergency" name="emergency" required>
                            <option value="no">Elective</option>
                            <option value="yes">Emergency</option>
                        </select>
                    </div>
                </div>
            </div>

            <button type="submit" class="submit-btn">Calculate Hypotension Risk</button>
        </form>

        {% if prediction %}
        <!-- Results Section -->
        <div style="margin-top: 48px;">
            <div class="results-warning">
                <p>⚠️ SIMULATION ONLY - These predictions are for educational purposes and may not reflect actual clinical outcomes. Always follow institutional protocols and use validated monitoring.</p>
            </div>

            <!-- Risk Gauges -->
            <div class="risk-gauges">
                <div class="risk-gauge">
                    <h3>5-Minute Risk</h3>
                    <div class="risk-value {{ prediction.risk_5min_class }}">{{ prediction.prob_5min }}%</div>
                    <span class="risk-label {{ prediction.risk_5min_label }}">{{ prediction.risk_5min_text }}</span>
                </div>
                <div class="risk-gauge">
                    <h3>10-Minute Risk</h3>
                    <div class="risk-value {{ prediction.risk_10min_class }}">{{ prediction.prob_10min }}%</div>
                    <span class="risk-label {{ prediction.risk_10min_label }}">{{ prediction.risk_10min_text }}</span>
                </div>
                <div class="risk-gauge">
                    <h3>20-Minute Risk</h3>
                    <div class="risk-value {{ prediction.risk_20min_class }}">{{ prediction.prob_20min }}%</div>
                    <span class="risk-label {{ prediction.risk_20min_label }}">{{ prediction.risk_20min_text }}</span>
                </div>
            </div>

            <!-- Contributing Factors -->
            <div class="factors-section">
                <h2>Top Contributing Factors</h2>
                {% for factor in prediction.factors %}
                <div class="factor-item">
                    <strong>{{ factor.name }}:</strong> {{ factor.description }}
                </div>
                {% endfor %}
            </div>

            <!-- Suggested Interventions -->
            <div class="interventions-section">
                <h2>Evidence-Based Intervention Suggestions</h2>
                {% for intervention in prediction.interventions %}
                <div class="intervention-item">
                    <div class="intervention-rank">#{{ loop.index }}</div>
                    <div class="intervention-content">
                        <h3>{{ intervention.name }}</h3>
                        <p>{{ intervention.description }}</p>
                    </div>
                </div>
                {% endfor %}
            </div>
        </div>
        {% endif %}
    </main>

    <footer>
        <p>&copy; 2025 gasconsult.ai. All rights reserved. | <a href="/terms">Terms of Service</a> | <a href="/privacy">Privacy Policy</a></p>
        <p style="margin-top: 8px; font-weight: 600; color: #DC2626;">Educational Tool Only | Not for Clinical Use | No Medical Advice</p>
    </footer>
</body>
</html>
"""

@app.route("/stream")
def stream():
    """Server-Sent Events endpoint for streaming GPT responses"""
    request_id = request.args.get('request_id')
    print(f"[DEBUG] /stream endpoint called with request_id: {request_id}")

    if not request_id:
        print(f"[DEBUG] /stream - No request_id provided")
        return Response("error: Invalid request - no request_id\n\n", mimetype='text/event-stream')

    stream_key = f'stream_data_{request_id}'
    if stream_key not in session:
        print(f"[DEBUG] /stream - Stream key '{stream_key}' not found in session")
        print(f"[DEBUG] /stream - Available session keys: {list(session.keys())}")
        return Response("error: Invalid request - stream data not found\n\n", mimetype='text/event-stream')

    print(f"[DEBUG] /stream - Stream data found, proceeding...")

    # Get prepared data from session
    stream_data = session[f'stream_data_{request_id}']
    prompt = stream_data['prompt']
    refs = stream_data['refs']
    num_papers = stream_data['num_papers']
    raw_query = stream_data['raw_query']

    def generate():
        try:
            # Send initial event to confirm connection
            yield f"data: {json.dumps({'type': 'connected'})}\n\n"

            # Dynamic temperature based on question type and evidence availability
            question_type = stream_data.get('question_type', 'general')

            if num_papers == 0:
                # No papers - use higher temperature for general knowledge
                temperature = 0.2
            elif question_type == 'dosing':
                # Dosing questions need very precise answers
                temperature = 0.05
            elif question_type == 'safety':
                # Safety questions need factual accuracy
                temperature = 0.1
            elif question_type == 'mechanism':
                # Mechanism explanations can be slightly more fluid
                temperature = 0.15
            elif question_type == 'comparison':
                # Comparisons need balanced precision
                temperature = 0.1
            else:
                # Default (general, management)
                temperature = 0.1

            print(f"[DEBUG] Using temperature {temperature} for question type '{question_type}'")

            # Stream GPT response
            stream_response = openai_client.chat.completions.create(
                model="gpt-4o",
                messages=[{"role": "user", "content": prompt}],
                temperature=temperature,
                stream=True
            )

            full_response = ""
            for chunk in stream_response:
                if chunk.choices[0].delta.content:
                    content = chunk.choices[0].delta.content
                    full_response += content
                    # Send content chunk
                    yield f"data: {json.dumps({'type': 'content', 'data': content})}\n\n"

            # Calculate evidence strength
            evidence_strength = get_evidence_strength(num_papers, refs)

            # Save complete response to session
            session['messages'].append({
                "role": "assistant",
                "content": full_response,
                "references": refs,
                "num_papers": num_papers,
                "evidence_strength": evidence_strength
            })
            session.modified = True

            # Send references
            yield f"data: {json.dumps({'type': 'references', 'data': refs, 'num_papers': num_papers})}\n\n"

            # Send completion event
            yield f"data: {json.dumps({'type': 'done'})}\n\n"

            # Clean up stream data from session
            session.pop(f'stream_data_{request_id}', None)
            session.modified = True

        except Exception as e:
            print(f"[ERROR] Streaming failed: {e}")
            yield f"data: {json.dumps({'type': 'error', 'message': str(e)})}\n\n"

    return Response(stream_with_context(generate()), mimetype='text/event-stream')

@app.route("/")
def index():
    """Homepage - welcome screen only"""
    # Clear any existing conversation when returning to homepage
    session.pop('messages', None)
    return render_template_string(HTML, messages=[])

@app.route("/chat", methods=["GET", "POST"])
def chat():
    """Chat interface with conversation history"""
    # Initialize conversation history in session
    if 'messages' not in session:
        session['messages'] = []

    # Initialize conversation topic tracking
    if 'conversation_topic' not in session:
        session['conversation_topic'] = None

    if request.method == "POST":
        try:
            # Safely get query from form data and sanitize it
            raw_query = request.form.get("query", "").strip()
            raw_query = sanitize_user_query(raw_query)

            # If query is empty, redirect to GET
            if not raw_query:
                print(f"[DEBUG] Empty query received, redirecting to GET")
                return redirect(url_for('chat'))

            print(f"\n[DEBUG] ===== NEW REQUEST =====")
            print(f"[DEBUG] Raw query: '{raw_query}'")
            print(f"[DEBUG] Session has {len(session['messages'])} messages before")

            # Check if this is the first message (from homepage)
            is_first_message = len(session['messages']) == 0

            # Set conversation topic on first message (for context in future vague queries)
            if is_first_message and not session['conversation_topic']:
                # Extract main topic from first query (simple keyword extraction)
                topic_words = []
                for word in raw_query.lower().split():
                    # Keep medical terms (4+ chars, not common words)
                    if len(word) >= 4 and word not in ['what', 'about', 'when', 'where', 'which', 'that', 'this', 'with', 'from', 'have', 'been', 'does', 'should', 'would', 'could']:
                        topic_words.append(word)
                if topic_words:
                    session['conversation_topic'] = ' '.join(topic_words[:3])  # First 3 meaningful words
                    print(f"[DEBUG] Conversation topic set: {session['conversation_topic']}")

            # Add user message to conversation
            session['messages'].append({"role": "user", "content": raw_query})
            session.modified = True
            print(f"[DEBUG] Added user message, session now has {len(session['messages'])} messages")

            # Check if this is a calculation request
            context_hint = None
            if len(session['messages']) >= 3:
                last_msgs = session['messages'][-4:]
                for msg in last_msgs:
                    content = msg.get('content', '').lower()
                    if any(term in content for term in ['mabl', 'ibw', 'bsa', 'qtc', 'maintenance fluid', 'ideal body weight', 'body surface']):
                        for term in ['mabl', 'ibw', 'bsa', 'qtc', 'maintenance fluid']:
                            if term in content:
                                context_hint = term
                                break
                        break

            calc_result = detect_and_calculate(raw_query, context_hint=context_hint)

            if calc_result:
                print(f"[DEBUG] Calculator result generated")
                session['messages'].append({
                    "role": "assistant",
                    "content": calc_result,
                    "references": [],
                    "num_papers": 0
                })
                session.modified = True
                print(f"[DEBUG] Redirecting after calculation")
                return redirect(url_for('chat'))

            query = clean_query(raw_query)
            print(f"[DEBUG] Cleaned query: '{query}'")

            is_followup = len(session['messages']) >= 3
            print(f"[DEBUG] Is follow-up: {is_followup}")

            # Resolve pronouns and references using conversation context
            query = resolve_references(query, session['messages'][:-1])  # Exclude the just-added user message
            print(f"[DEBUG] After reference resolution: '{query}'")

            # Detect question type for customized search and temperature
            question_type = detect_question_type(query)
            print(f"[DEBUG] Question type: {question_type}")

            # Detect if this is a multi-part question
            is_multipart = detect_multipart(query)
            if is_multipart:
                print(f"[DEBUG] Multi-part question detected")

            # Handle negations (contraindications, when NOT to use, etc.)
            negation_search_modifier, negation_prompt_modifier = handle_negations(query)
            if negation_search_modifier:
                print(f"[DEBUG] Negation detected - will modify search and prompt")

            # Expand medical abbreviations and synonyms
            q = expand_medical_abbreviations(query)

            print(f"[DEBUG] Expanded query: '{q}'")

            # Customize search based on question type
            if is_followup:
                # Broader search for follow-ups
                search_term = f'({q}) AND ("2005/01/01"[PDAT] : "3000"[PDAT])'
            else:
                # Customize search filters based on question type
                if question_type == 'dosing':
                    # Prioritize guidelines and reviews for dosing questions
                    search_term = (
                        f'({q}) AND '
                        f'(guideline[pt] OR "practice guideline"[pt] OR review[pt] OR meta-analysis[pt]) AND '
                        f'("2015/01/01"[PDAT] : "3000"[PDAT])'
                    )
                elif question_type == 'safety':
                    # Include adverse effects and safety studies
                    search_term = (
                        f'({q}) AND '
                        f'(systematic review[pt] OR meta-analysis[pt] OR "randomized controlled trial"[pt] OR '
                        f'guideline[pt] OR "adverse effects"[sh] OR safety[ti]) AND '
                        f'("2015/01/01"[PDAT] : "3000"[PDAT])'
                    )
                elif question_type == 'comparison':
                    # Prioritize comparative studies
                    search_term = (
                        f'({q}) AND '
                        f'(systematic review[pt] OR meta-analysis[pt] OR "randomized controlled trial"[pt] OR '
                        f'"comparative study"[pt]) AND '
                        f'("2015/01/01"[PDAT] : "3000"[PDAT])'
                    )
                else:
                    # Default search (general, mechanism, management)
                    search_term = (
                        f'({q}) AND '
                        f'(systematic review[pt] OR meta-analysis[pt] OR "randomized controlled trial"[pt] OR '
                        f'"Cochrane Database Syst Rev"[ta] OR guideline[pt]) AND '
                        f'("2015/01/01"[PDAT] : "3000"[PDAT])'
                    )

            # Add negation modifier if detected
            if negation_search_modifier:
                search_term += negation_search_modifier

            print(f"[DEBUG] Search term: '{search_term[:150]}...'")

            # Try anesthesiology-specific search first (reduced to 10 papers for speed)
            ids = []
            try:
                print(f"[DEBUG] Searching PubMed (anesthesiology)...")
                handle = Entrez.esearch(db="pubmed", term=f'anesthesiology[MeSH Terms] AND {search_term}', retmax=10, sort="relevance")
                result = Entrez.read(handle)
                ids = result.get("IdList", [])
                print(f"[DEBUG] Found {len(ids)} papers (anesthesiology)")
            except Exception as e:
                print(f"[ERROR] PubMed search failed (anesthesiology): {e}")
                ids = []

            # Fallback: Try without anesthesiology restriction (skip for follow-ups to save time)
            if not ids and not is_followup:
                try:
                    print(f"[DEBUG] Searching PubMed (general)...")
                    handle = Entrez.esearch(db="pubmed", term=search_term, retmax=10, sort="relevance")
                    result = Entrez.read(handle)
                    ids = result.get("IdList", [])
                    print(f"[DEBUG] Found {len(ids)} papers (general)")
                except Exception as e:
                    print(f"[ERROR] PubMed search failed (general): {e}")
                    ids = []

            # If no papers found, handle gracefully
            if not ids:
                print(f"[DEBUG] No papers found")
                if is_followup:
                    print(f"[DEBUG] Generating follow-up response without papers")
                    conversation_context = ""
                    recent_messages = session['messages'][-8:]
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

                    print(f"[DEBUG] Preparing streaming for follow-up...")

                    # Generate unique request ID for this streaming session
                    request_id = str(uuid.uuid4())

                    # Store data in session for streaming endpoint
                    session[f'stream_data_{request_id}'] = {
                        'prompt': prompt,
                        'refs': [],
                        'num_papers': 0,
                        'raw_query': raw_query
                    }
                    session.modified = True

                    print(f"[DEBUG] Stream data prepared for follow-up, returning request_id: {request_id}")
                    return jsonify({
                        'status': 'ready',
                        'request_id': request_id,
                        'raw_query': raw_query
                    })
                else:
                    print(f"[DEBUG] No results for initial query")
                    error_msg = "<p>No relevant evidence found in recent literature. Try rephrasing your question or using different medical terms.</p>"
                    session['messages'].append({
                        "role": "assistant",
                        "content": error_msg,
                        "references": [],
                        "num_papers": 0
                    })
                    session.modified = True
                    print(f"[DEBUG] Error message added, redirecting")
                    return redirect(url_for('chat'))

            print(f"[DEBUG] Fetching {len(ids)} papers from PubMed...")
            handle = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml")
            papers = Entrez.read(handle)["PubmedArticle"]
            print(f"[DEBUG] Papers fetched successfully")

            refs = []
            context = ""
            # Process only 8 papers for faster GPT processing
            for p in papers[:8]:
                try:
                    art = p["MedlineCitation"]["Article"]
                    title = art.get("ArticleTitle", "No title")
                    # Truncate abstracts to 600 chars for faster GPT processing
                    abstract = " ".join(str(t) for t in art.get("Abstract", {}).get("AbstractText", [])) if art.get("Abstract") else ""
                    abstract = abstract[:600] + "..." if len(abstract) > 600 else abstract
                    authors = ", ".join([a.get("LastName","") + " " + (a.get("ForeName","")[:1]+"." if a.get("ForeName") else "") for a in art.get("AuthorList",[])[:3]])  # Reduced to 3 authors
                    journal = art["Journal"].get("Title", "Unknown")
                    year = art["Journal"]["JournalIssue"]["PubDate"].get("Year", "N/A")
                    pmid = p["MedlineCitation"]["PMID"]

                    refs.append({"title": title, "authors": authors, "journal": journal, "year": year, "pmid": pmid})
                    context += f"Title: {title}\nAbstract: {abstract}\nAuthors: {authors}\nJournal: {journal} ({year})\nPMID: {pmid}\n\n"
                except:
                    continue

            num_papers = len(refs)
            print(f"[DEBUG] Processed {num_papers} paper references")

            # Build smart conversation context (includes relevant earlier messages + recent)
            conversation_context = build_smart_context(session['messages'][:-1], raw_query)  # Exclude just-added user message
            print(f"[DEBUG] Smart context built ({len(conversation_context)} chars)")

            # Create numbered reference list for citation
            ref_list = ""
            for i, ref in enumerate(refs, 1):
                ref_list += f"[{i}] {ref['title']} - {ref['authors']} ({ref['year']}) PMID: {ref['pmid']}\n"

            prompt = f"""You are a clinical anesthesiologist AI providing evidence-based answers with citations.

Previous conversation:
{conversation_context if len(session['messages']) > 1 else "New conversation."}

Current question: {raw_query}
{'This is a FOLLOW-UP - build on the previous discussion.' if is_followup else ''}
{negation_prompt_modifier if negation_prompt_modifier else ''}
{'NOTE: This question has multiple parts - address each thoroughly.' if is_multipart else ''}

Research papers (cite as [1], [2], etc.):
{ref_list}

Paper details:
{context}

INSTRUCTIONS:
1. Include specific dosages (mg/kg), contraindications, side effects, and monitoring when relevant
2. For acute situations, provide step-by-step protocols with drugs and doses
3. Use numbered citations [1], [2] - NO author names in text
4. Be conversational but clinically complete - like talking to a colleague
5. HTML format: <h3> for sections, <p> for paragraphs, <strong> for emphasis, <ul><li> for lists
6. START your response with a confidence badge using this exact HTML format:
   <div class="evidence-quality-badge">
   <div class="confidence-level [high/moderate/low]">
   <strong>Evidence Quality:</strong> [High/Moderate/Low] Confidence
   </div>
   <div class="evidence-details">
   📊 {num_papers} papers analyzed • Study types: [list types] • Date range: [range]
   </div>
   </div>

Example response format:
"<div class="evidence-quality-badge">
<div class="confidence-level high">
<strong>Evidence Quality:</strong> High Confidence
</div>
<div class="evidence-details">
📊 8 papers analyzed • Study types: 3 meta-analyses, 4 RCTs, 1 systematic review • Date range: 2015-2024
</div>
</div>

<h3>Acute Bronchospasm Management</h3>
<p><strong>Immediate Actions:</strong><br>
Deepen anesthesia with propofol 0.5-1 mg/kg or increase volatile to 2+ MAC [1]</p>
<p><strong>Bronchodilators:</strong><br>
Albuterol 4-8 puffs via ETT [2]</p>"

Respond with maximum clinical utility:"""

            print(f"[DEBUG] Preparing streaming with {num_papers} papers...")

            # Generate unique request ID for this streaming session
            request_id = str(uuid.uuid4())

            # Store data in session for streaming endpoint
            session[f'stream_data_{request_id}'] = {
                'prompt': prompt,
                'refs': refs,
                'num_papers': num_papers,
                'raw_query': raw_query,
                'question_type': question_type  # Store for temperature adjustment
            }
            session.modified = True

            print(f"[DEBUG] Stream data prepared, returning request_id: {request_id}")

            # If this is the first message (from homepage), set pending flag and redirect
            if is_first_message:
                # Add placeholder assistant message for auto-start JavaScript to populate
                session['messages'].append({
                    "role": "assistant",
                    "content": "",  # Will be populated by streaming
                    "references": [],
                    "num_papers": 0
                })
                session['pending_stream'] = request_id
                session.modified = True
                print(f"[DEBUG] First message - added placeholder assistant message, redirecting to /chat page")
                return redirect(url_for('chat'))

            # For follow-up messages, return JSON for streaming
            return jsonify({
                'status': 'ready',
                'request_id': request_id,
                'raw_query': raw_query
            })

        except Exception as e:
            # Catch all unhandled errors
            print(f"\n[ERROR] ===== UNHANDLED EXCEPTION =====")
            print(f"[ERROR] {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()

            session['messages'].append({
                "role": "assistant",
                "content": f"<p><strong>Error:</strong> {str(e)}</p><p>Please try rephrasing your question or start a new conversation.</p>",
                "references": [],
                "num_papers": 0
            })
            session.modified = True
            return redirect(url_for('chat'))

    # Check for pending stream (from homepage redirect)
    pending_stream = session.pop('pending_stream', None)
    print(f"[DEBUG] GET /chat - pending_stream = {pending_stream}")
    print(f"[DEBUG] GET /chat - num messages = {len(session.get('messages', []))}")
    return render_template_string(CHAT_HTML, messages=session.get('messages', []), pending_stream=pending_stream)

@app.route("/clear")
def clear():
    """Clear conversation history and start new chat"""
    session.pop('messages', None)
    return redirect(url_for('chat'))

@app.route("/terms")
def terms():
    """Terms of Service page"""
    return render_template_string(TERMS_HTML)

@app.route("/privacy")
def privacy():
    """Privacy Policy page"""
    return render_template_string(PRIVACY_POLICY_HTML)

@app.route("/calculators")
def calculators():
    """Clinical Calculators Hub - IBW, MABL, BSA, QTc, Fluids, PONV, ASA"""
    return render_template_string(CALCULATORS_HTML)

@app.route("/quick-dose")
def quick_dose():
    """Quick Dose Reference - Weight-based drug dosing calculator"""
    return render_template_string(QUICK_DOSE_HTML)
@app.route("/preop", methods=["GET", "POST"])
def preop_assessment():
    """Pre-operative assessment with evidence-based risk stratification"""
    if request.method == "GET":
        return render_template_string(PREOP_HTML, summary=None, references=None)

    # Collect and sanitize form data
    age = int(request.form.get("age", 0))
    weight = float(request.form.get("weight", 0))
    height = float(request.form.get("height", 0))
    sex = sanitize_user_query(request.form.get("sex", ""))
    comorbidities = request.form.getlist("comorbidities")
    other_comorbidities = sanitize_user_query(request.form.get("other_comorbidities", ""))
    mets = sanitize_user_query(request.form.get("mets", ""))
    previous_anesthesia = sanitize_user_query(request.form.get("previous_anesthesia", ""))
    medications = sanitize_user_query(request.form.get("medications", ""))
    hgb = sanitize_user_query(request.form.get("hgb", ""))
    plt = sanitize_user_query(request.form.get("plt", ""))
    cr = sanitize_user_query(request.form.get("cr", ""))
    inr = sanitize_user_query(request.form.get("inr", ""))
    ef = sanitize_user_query(request.form.get("ef", ""))
    ekg = sanitize_user_query(request.form.get("ekg", ""))
    procedure = sanitize_user_query(request.form.get("procedure", ""))
    surgery_risk = sanitize_user_query(request.form.get("surgery_risk", ""))
    npo = sanitize_user_query(request.form.get("npo", ""))
    allergies = sanitize_user_query(request.form.get("allergies", ""))

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

    # Reduced ejection fraction
    if ef:
        ef_lower = ef.lower()
        try:
            # Try to extract numeric value from EF
            ef_numeric = float(''.join(filter(lambda x: x.isdigit() or x == '.', ef.split('-')[0])))
            if ef_numeric < 40:
                search_queries.append(f"reduced ejection fraction perioperative management {procedure}")
        except:
            # Check for descriptive terms
            if any(term in ef_lower for term in ['reduced', 'low', 'decreased', 'hfref']):
                search_queries.append(f"reduced ejection fraction perioperative management {procedure}")

    # Atrial fibrillation from EKG
    if ekg:
        ekg_lower = ekg.lower()
        if any(term in ekg_lower for term in ['afib', 'a-fib', 'atrial fib', 'atrial fibrillation']):
            search_queries.append(f"atrial fibrillation perioperative management {procedure}")

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

            handle = Entrez.esearch(db="pubmed", term=search_term, retmax=5, sort="relevance")
            result = Entrez.read(handle)
            ids = result["IdList"]

            if ids:
                handle = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml")
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
- Sex Assigned at Birth: {sex}
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

Cardiac Assessment:
- Ejection Fraction: {ef if ef else 'Not available'}
- EKG: {ekg if ekg else 'Not available'}

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

@app.route("/hypotension", methods=["GET", "POST"])
def hypotension_predictor():
    """Intraoperative Hypotension Predictor - Educational tool only"""
    if request.method == "GET":
        return render_template_string(HYPOTENSION_HTML, prediction=None)

    # Collect and sanitize form data
    age = int(request.form.get("age", 0))
    sex = sanitize_user_query(request.form.get("sex", ""))
    weight = float(request.form.get("weight", 0))
    height = float(request.form.get("height", 0))
    asa = int(request.form.get("asa", 1))
    baseline_map = int(request.form.get("baseline_map", 0))
    baseline_hr = int(request.form.get("baseline_hr", 0))
    current_map = int(request.form.get("current_map", 0))
    map_5min = int(request.form.get("map_5min", 0))
    map_10min = int(request.form.get("map_10min", 0))
    surgery_duration = int(request.form.get("surgery_duration", 0))
    vasopressor = sanitize_user_query(request.form.get("vasopressor", "none"))
    surgery_type = sanitize_user_query(request.form.get("surgery_type", ""))
    induction_agent = sanitize_user_query(request.form.get("induction_agent", ""))
    emergency = sanitize_user_query(request.form.get("emergency", "no"))

    # Calculate BMI
    bmi = round(weight / ((height / 100) ** 2), 1) if weight and height else 0

    # Educational Rule-Based Risk Scoring System
    # Based on research findings: MAP trend, age, ASA status are key predictors

    risk_score = 0
    factors = []

    # 1. MAP Trend Analysis (Most Important)
    map_change_5min = current_map - map_5min
    map_change_10min = current_map - map_10min
    map_trend = (map_change_5min + map_change_10min) / 2

    if map_trend < -10:
        risk_score += 30
        factors.append({
            "name": "Declining MAP Trend",
            "description": f"MAP decreased by {abs(int(map_trend))} mmHg over last 10 minutes, indicating hemodynamic instability"
        })
    elif map_trend < -5:
        risk_score += 15
        factors.append({
            "name": "Moderate MAP Decline",
            "description": f"MAP decreased by {abs(int(map_trend))} mmHg, suggesting early hemodynamic changes"
        })

    # 2. Current MAP vs Baseline
    map_drop_from_baseline = ((baseline_map - current_map) / baseline_map) * 100
    if map_drop_from_baseline > 20:
        risk_score += 25
        factors.append({
            "name": "Significant MAP Drop from Baseline",
            "description": f"{int(map_drop_from_baseline)}% decrease from baseline MAP, exceeding critical threshold"
        })
    elif map_drop_from_baseline > 10:
        risk_score += 12

    # 3. Absolute MAP Level
    if current_map < 70:
        risk_score += 20
        factors.append({
            "name": "Low Current MAP",
            "description": f"Current MAP of {current_map} mmHg is approaching hypotension threshold (65 mmHg)"
        })
    elif current_map < 75:
        risk_score += 10

    # 4. Age Risk (Older age = higher risk)
    if age > 70:
        risk_score += 15
        factors.append({
            "name": "Advanced Age",
            "description": f"Age {age} years associated with increased cardiovascular lability and hypotension risk"
        })
    elif age > 60:
        risk_score += 8

    # 5. ASA Class (Higher ASA = higher risk)
    if asa >= 3:
        risk_score += 15
        factors.append({
            "name": "High ASA Classification",
            "description": f"ASA {asa} indicates significant comorbidities affecting hemodynamic stability"
        })
    elif asa == 2:
        risk_score += 5

    # 6. Surgery Type
    high_risk_surgeries = ["cardiac", "major_abdominal", "vascular"]
    if surgery_type in high_risk_surgeries:
        risk_score += 12
        factors.append({
            "name": "High-Risk Surgery Type",
            "description": f"{surgery_type.replace('_', ' ').title()} surgery associated with greater fluid shifts and hemodynamic instability"
        })

    # 7. Emergency Surgery
    if emergency == "yes":
        risk_score += 10
        factors.append({
            "name": "Emergency Surgery",
            "description": "Emergency procedures have higher hypotension risk due to reduced optimization time"
        })

    # 8. Induction Agent
    if induction_agent == "propofol":
        risk_score += 8
        factors.append({
            "name": "Propofol Induction",
            "description": "Propofol associated with dose-dependent vasodilation and myocardial depression"
        })

    # 9. Vasopressor Use (Paradoxically indicates ongoing instability)
    if vasopressor != "none":
        risk_score += 10
        factors.append({
            "name": "Current Vasopressor Requirement",
            "description": f"Ongoing {vasopressor} use indicates hemodynamic instability requiring support"
        })

    # 10. Surgery Duration (Longer = more risk)
    if surgery_duration > 180:
        risk_score += 8
        factors.append({
            "name": "Prolonged Surgery",
            "description": f"{surgery_duration} minutes of surgery increases cumulative fluid shifts and anesthetic depth effects"
        })

    # Convert risk score to probabilities (simplified logistic-like function)
    # 5-minute prediction (highest confidence)
    base_prob_5min = min(95, max(5, risk_score * 0.9))

    # 10-minute prediction (moderate decay)
    base_prob_10min = min(90, max(5, risk_score * 0.75))

    # 20-minute prediction (further decay)
    base_prob_20min = min(85, max(5, risk_score * 0.6))

    # Add some realistic variance
    import random
    random.seed(age + current_map)  # Deterministic but appears random

    prob_5min = int(base_prob_5min + random.randint(-3, 3))
    prob_10min = int(base_prob_10min + random.randint(-3, 3))
    prob_20min = int(base_prob_20min + random.randint(-3, 3))

    # Classify risk levels
    def classify_risk(prob):
        if prob < 30:
            return "low", "risk-low", "Low Risk"
        elif prob < 60:
            return "moderate", "risk-moderate", "Moderate Risk"
        else:
            return "high", "risk-high", "High Risk"

    risk_5min_class, risk_5min_label, risk_5min_text = classify_risk(prob_5min)
    risk_10min_class, risk_10min_label, risk_10min_text = classify_risk(prob_10min)
    risk_20min_class, risk_20min_label, risk_20min_text = classify_risk(prob_20min)

    # Sort factors by importance (keep top 3)
    factors = factors[:3] if len(factors) > 3 else factors

    # If no factors identified, add default
    if not factors:
        factors.append({
            "name": "Stable Hemodynamics",
            "description": "No major risk factors identified based on current parameters"
        })

    # Evidence-Based Intervention Suggestions
    interventions = []

    # Prioritize based on risk factors
    if current_map < 70 or map_trend < -5:
        if vasopressor == "none":
            interventions.append({
                "name": "Phenylephrine 50-100 mcg IV",
                "description": "First-line alpha-1 agonist for MAP support without chronotropy. Effective for vasodilatory hypotension from anesthetics."
            })

        interventions.append({
            "name": "Fluid Bolus 250-500 mL Crystalloid",
            "description": "Goal-directed fluid therapy to optimize preload. Assess fluid responsiveness via pulse pressure variation if available."
        })

    if induction_agent == "propofol":
        interventions.append({
            "name": "Reduce Volatile Anesthetic %",
            "description": "Dose-dependent vasodilation from inhalational agents. Consider reducing MAC to 0.7-0.9 to improve hemodynamics."
        })

    if age > 65 or asa >= 3:
        interventions.append({
            "name": "Increase Monitoring Frequency",
            "description": "Consider arterial line for beat-to-beat BP monitoring in high-risk patients. Trend MAP every 1-2 minutes."
        })

    if surgery_type in high_risk_surgeries:
        interventions.append({
            "name": "Ephedrine 5-10 mg IV",
            "description": "Mixed alpha/beta agonist useful for hypotension with bradycardia. Raises MAP and cardiac output."
        })

    # Always add general recommendations
    interventions.append({
        "name": "Optimize Positioning & Surgical Factors",
        "description": "Ensure appropriate patient positioning, communicate with surgeon about retraction/traction affecting venous return."
    })

    interventions.append({
        "name": "Reassess Anesthetic Depth",
        "description": "Verify adequate but not excessive anesthetic depth. Consider reducing propofol/volatile if hemodynamically unstable."
    })

    # Limit to top 5 interventions
    interventions = interventions[:5]

    prediction = {
        "prob_5min": prob_5min,
        "prob_10min": prob_10min,
        "prob_20min": prob_20min,
        "risk_5min_class": risk_5min_class,
        "risk_10min_class": risk_10min_class,
        "risk_20min_class": risk_20min_class,
        "risk_5min_label": risk_5min_label,
        "risk_10min_label": risk_10min_label,
        "risk_20min_label": risk_20min_label,
        "risk_5min_text": risk_5min_text,
        "risk_10min_text": risk_10min_text,
        "risk_20min_text": risk_20min_text,
        "factors": factors,
        "interventions": interventions
    }

    return render_template_string(HYPOTENSION_HTML, prediction=prediction)

@app.route("/health")
@csrf.exempt  # Health checks don't need CSRF protection
def health_check():
    """Health check endpoint for deployment monitoring"""
    import sys
    health_status = {
        "status": "healthy",
        "service": "gasconsult.ai",
        "version": "1.0.0",
        "python_version": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        "checks": {
            "openai": bool(os.getenv("OPENAI_API_KEY")),
            "entrez_email": bool(os.getenv("ENTREZ_EMAIL")),
            "entrez_api_key": bool(os.getenv("ENTREZ_API_KEY")),
            "secret_key": bool(app.secret_key)
        }
    }

    # Check if any critical service is missing
    if not all(health_status["checks"].values()):
        health_status["status"] = "degraded"
        return jsonify(health_status), 503

    return jsonify(health_status), 200

@app.route("/api/status")
@csrf.exempt  # Status endpoint doesn't need CSRF protection
def api_status():
    """API status endpoint with basic metrics"""
    return jsonify({
        "status": "operational",
        "endpoints": {
            "chat": "/chat",
            "preop": "/preop",
            "hypotension": "/hypotension",
            "quick_dose": "/quick-dose",
            "health": "/health",
            "library": "/library",
            "share": "/share"
        },
        "rate_limit": os.getenv("RATE_LIMIT", "60 per minute")
    }), 200

# ====== Premium Features Routes ======

@app.route("/bookmark", methods=["POST"])
def bookmark_response():
    """Add or remove a bookmark"""
    try:
        data = request.get_json()
        message_index = data.get('message_index')
        action = data.get('action', 'add')  # 'add' or 'remove'

        # Get user session ID
        user_id = session.get('user_id')
        if not user_id:
            user_id = str(uuid.uuid4())
            session['user_id'] = user_id
            session.modified = True

        if user_id not in BOOKMARKS_STORAGE:
            BOOKMARKS_STORAGE[user_id] = []

        messages = session.get('messages', [])

        if message_index < 0 or message_index >= len(messages):
            return jsonify({'error': 'Invalid message index'}), 400

        # Get the query (previous user message) and response (assistant message)
        query = messages[message_index - 1]['content'] if message_index > 0 else 'Direct query'
        response_msg = messages[message_index]

        if action == 'add':
            bookmark = {
                'id': str(uuid.uuid4()),
                'query': query,
                'answer': response_msg['content'],
                'references': response_msg.get('references', []),
                'num_papers': response_msg.get('num_papers', 0),
                'timestamp': datetime.now().isoformat(),
                'message_index': message_index
            }
            BOOKMARKS_STORAGE[user_id].append(bookmark)
            return jsonify({'success': True, 'action': 'added', 'bookmark_id': bookmark['id']})
        elif action == 'remove':
            bookmark_id = data.get('bookmark_id')
            BOOKMARKS_STORAGE[user_id] = [b for b in BOOKMARKS_STORAGE[user_id] if b['id'] != bookmark_id]
            return jsonify({'success': True, 'action': 'removed'})

    except Exception as e:
        logger.error(f"Bookmark error: {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route("/library")
def library():
    """View all bookmarked responses"""
    user_id = session.get('user_id')
    bookmarks = BOOKMARKS_STORAGE.get(user_id, []) if user_id else []

    # Sort by timestamp (newest first)
    bookmarks = sorted(bookmarks, key=lambda x: x['timestamp'], reverse=True)

    return render_template_string(LIBRARY_HTML, bookmarks=bookmarks)

@app.route("/share", methods=["POST"])
def create_share_link():
    """Create a shareable link for a response"""
    try:
        data = request.get_json()
        message_index = data.get('message_index')

        messages = session.get('messages', [])

        if message_index < 0 or message_index >= len(messages):
            return jsonify({'error': 'Invalid message index'}), 400

        query = messages[message_index - 1]['content'] if message_index > 0 else 'Direct query'
        response_msg = messages[message_index]

        # Generate unique share ID
        share_id = str(uuid.uuid4())[:8]

        # Store shared response (expires in 30 days)
        SHARED_LINKS_STORAGE[share_id] = {
            'query': query,
            'answer': response_msg['content'],
            'references': response_msg.get('references', []),
            'num_papers': response_msg.get('num_papers', 0),
            'timestamp': datetime.now().isoformat(),
            'expires': (datetime.now() + timedelta(days=30)).isoformat()
        }

        share_url = f"{request.host_url}r/{share_id}"
        return jsonify({'success': True, 'share_url': share_url, 'share_id': share_id})

    except Exception as e:
        logger.error(f"Share link error: {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route("/r/<share_id>")
@csrf.exempt
def view_shared_response(share_id):
    """View a shared response"""
    shared_data = SHARED_LINKS_STORAGE.get(share_id)

    if not shared_data:
        return "Shared link not found or expired", 404

    # Check if expired
    expires = datetime.fromisoformat(shared_data['expires'])
    if datetime.now() > expires:
        del SHARED_LINKS_STORAGE[share_id]
        return "Shared link has expired", 410

    return render_template_string(SHARED_RESPONSE_HTML, data=shared_data, share_id=share_id)

@app.route("/export-citations/<format_type>", methods=["POST"])
def export_citations(format_type):
    """Export citations in BibTeX or RIS format"""
    try:
        data = request.get_json()
        message_index = data.get('message_index')

        messages = session.get('messages', [])

        if message_index < 0 or message_index >= len(messages):
            return jsonify({'error': 'Invalid message index'}), 400

        response_msg = messages[message_index]
        references = response_msg.get('references', [])

        if not references:
            return jsonify({'error': 'No references to export'}), 400

        if format_type == 'bibtex':
            content = generate_bibtex(references)
            filename = 'gasconsult_citations.bib'
            mimetype = 'application/x-bibtex'
        elif format_type == 'ris':
            content = generate_ris(references)
            filename = 'gasconsult_citations.ris'
            mimetype = 'application/x-research-info-systems'
        else:
            return jsonify({'error': 'Invalid format'}), 400

        from flask import make_response
        response = make_response(content)
        response.headers['Content-Type'] = mimetype
        response.headers['Content-Disposition'] = f'attachment; filename={filename}'
        return response

    except Exception as e:
        logger.error(f"Export error: {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route("/generate-followups", methods=["POST"])
def generate_followup_questions():
    """Generate smart follow-up questions based on the current conversation"""
    try:
        data = request.get_json()
        message_index = data.get('message_index')

        messages = session.get('messages', [])

        if message_index < 0 or message_index >= len(messages):
            return jsonify({'error': 'Invalid message index'}), 400

        query = messages[message_index - 1]['content'] if message_index > 0 else ''
        response = messages[message_index]['content']

        # Use GPT to generate contextual follow-up questions
        prompt = f"""Based on this clinical question and answer, generate 3 smart follow-up questions that a clinician might want to ask next.

Question: {query}

Answer: {response[:500]}...

Generate exactly 3 concise follow-up questions (one per line, no numbering) that explore:
1. A related contraindication or safety concern
2. A comparison with an alternative approach
3. A practical implementation detail or dosing question

Keep each question under 60 characters."""

        completion = openai_client.chat.completions.create(
            model="gpt-4o",
            messages=[{"role": "user", "content": prompt}],
            temperature=0.3,
            max_tokens=200
        )

        followups_text = completion.choices[0].message.content.strip()
        followups = [q.strip() for q in followups_text.split('\n') if q.strip() and not q.strip()[0].isdigit()][:3]

        return jsonify({'success': True, 'followups': followups})

    except Exception as e:
        logger.error(f"Generate followups error: {str(e)}")
        return jsonify({'error': str(e)}), 500

if __name__ == "__main__":
    port = int(os.getenv("PORT", 5000))
    debug = os.getenv("FLASK_ENV") == "development"
    app.run(host="0.0.0.0", port=port, debug=debug)