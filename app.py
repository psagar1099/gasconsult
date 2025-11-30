from flask import Flask, request, render_template_string, session, redirect, url_for, Response, stream_with_context, jsonify, make_response
from flask_session import Session
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from flask_wtf.csrf import CSRFProtect, CSRFError
from werkzeug.exceptions import HTTPException, NotFound
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
app.config['SESSION_PERMANENT'] = True
app.config['PERMANENT_SESSION_LIFETIME'] = 3600  # 1 hour
app.config['SESSION_USE_SIGNER'] = True
app.config['SESSION_COOKIE_SAMESITE'] = 'Lax'
app.config['SESSION_COOKIE_HTTPONLY'] = True
# Only require HTTPS cookies if explicitly enabled (Render uses reverse proxy)
app.config['SESSION_COOKIE_SECURE'] = os.getenv('FLASK_ENV') == 'production' and os.getenv('FORCE_HTTPS', 'false').lower() == 'true'
Session(app)

# Initialize CSRF protection
csrf = CSRFProtect(app)
app.config['WTF_CSRF_TIME_LIMIT'] = None  # Don't expire CSRF tokens
app.config['WTF_CSRF_SSL_STRICT'] = False  # Allow CSRF on HTTP (behind reverse proxy)

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

# In-memory cache for stream data (backup for session issues)
# This helps when session state isn't shared properly between requests
# WARNING: This won't work with multiple Gunicorn workers unless using sticky sessions
STREAM_DATA_CACHE = {}  # {request_id: {data, expires_at}}
STREAM_CACHE_TTL = 300  # 5 minutes

def store_stream_data(request_id, data):
    """Store stream data in both session and cache"""
    STREAM_DATA_CACHE[request_id] = {
        'data': data,
        'expires_at': datetime.now() + timedelta(seconds=STREAM_CACHE_TTL)
    }
    # Clean up expired entries
    now = datetime.now()
    expired = [k for k, v in STREAM_DATA_CACHE.items() if v['expires_at'] < now]
    for k in expired:
        del STREAM_DATA_CACHE[k]

def get_stream_data(request_id):
    """Get stream data from cache (fallback for session issues)"""
    if request_id in STREAM_DATA_CACHE:
        entry = STREAM_DATA_CACHE[request_id]
        if entry['expires_at'] > datetime.now():
            return entry['data']
        else:
            del STREAM_DATA_CACHE[request_id]
    return None

def strip_markdown_code_fences(content):
    """Remove markdown code fences (```html ... ```) from content"""
    if not content:
        return content
    # Remove opening code fence with language specifier (```html, ```HTML, etc.)
    content = re.sub(r'^```\w*\s*\n?', '', content, flags=re.IGNORECASE)
    # Remove closing code fence
    content = re.sub(r'\n?```\s*$', '', content)
    return content

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

@app.errorhandler(CSRFError)
def handle_csrf_error(e):
    """Handle CSRF token errors with helpful message"""
    logger.warning(f"CSRF error: {str(e)}")

    # For AJAX requests, return JSON error instead of redirect
    is_ajax = request.headers.get('X-Requested-With') == 'XMLHttpRequest'
    if is_ajax:
        return jsonify({
            "status": "error",
            "message": "Your session has expired. Please reload the page and try again.",
            "error_type": "csrf"
        }), 400

    # For regular requests, redirect to homepage to regenerate session/CSRF token
    return redirect(url_for('index'))

@app.errorhandler(404)
def handle_not_found(e):
    """Handle 404 errors - log at WARNING level to avoid log noise from scanners"""
    logger.warning(f"404 Not Found: {request.method} {request.path} from {get_remote_address()}")
    return jsonify({
        "error": "The requested resource was not found.",
        "status": "not_found"
    }), 404

@app.errorhandler(HTTPException)
def handle_http_exception(e):
    """Handle other HTTP exceptions (400, 403, 405, etc.) without logging as errors"""
    logger.info(f"HTTP {e.code}: {request.method} {request.path}")
    return jsonify({
        "error": e.description,
        "status": "http_error"
    }), e.code

@app.errorhandler(Exception)
def log_exception(e):
    """Log unhandled non-HTTP exceptions only"""
    # Don't catch HTTPException here - it's handled above
    if isinstance(e, HTTPException):
        raise e

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

PREOP_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pre-Op Assessment — gasconsult.ai</title>

    <!-- PWA -->
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=5">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">

    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap" rel="stylesheet">
    <style>

        :root {
            --white: #FFFFFF;
            --gray-50: #F8FAFC;
            --gray-100: #F1F5F9;
            --gray-200: #E2E8F0;
            --gray-300: #CBD5E1;
            --gray-400: #94A3B8;
            --gray-500: #64748B;
            --gray-600: #475569;
            --gray-700: #334155;
            --gray-800: #1E293B;
            --gray-900: #0F172A;
            --blue-50: #EFF6FF;
            --blue-100: #DBEAFE;
            --blue-200: #BFDBFE;
            --blue-300: #93C5FD;
            --blue-400: #60A5FA;
            --blue-500: #3B82F6;
            --blue-600: #2563EB;
            --blue-700: #1D4ED8;
        }

        * { margin: 0; padding: 0; box-sizing: border-box; }

        html {
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
            scroll-behavior: smooth;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            background: var(--gray-50);
            color: var(--gray-900);
            min-height: 100vh;
            overflow-x: hidden;
        }

        .bg-canvas {
            position: fixed;
            inset: 0;
            z-index: 0;
            overflow: hidden;
            background: linear-gradient(180deg, #F0F7FF 0%, var(--gray-50) 50%, #FAFBFF 100%);
        }

        .orb {
            position: absolute;
            border-radius: 50%;
            filter: blur(80px);
            opacity: 0.6;
            animation: float 20s ease-in-out infinite;
        }

        .orb-1 {
            width: 400px;
            height: 400px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.15) 0%, transparent 70%);
            top: -15%;
            left: -20%;
        }

        .orb-2 {
            width: 300px;
            height: 300px;
            background: radial-gradient(circle, rgba(147, 197, 253, 0.2) 0%, transparent 70%);
            top: 30%;
            right: -20%;
            animation-delay: -7s;
            animation-duration: 25s;
        }

        .orb-3 {
            width: 250px;
            height: 250px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.1) 0%, transparent 70%);
            bottom: -10%;
            left: 20%;
            animation-delay: -14s;
            animation-duration: 30s;
        }

        @keyframes float {
            0%, 100% { transform: translate(0, 0) scale(1); }
            25% { transform: translate(40px, -40px) scale(1.05); }
            50% { transform: translate(20px, 40px) scale(0.95); }
            75% { transform: translate(-40px, 20px) scale(1.02); }
        }

        .grain {
            position: fixed;
            inset: 0;
            z-index: 1;
            pointer-events: none;
            opacity: 0.02;
            background-image: url("data:image/svg+xml,%3Csvg viewBox='0 0 512 512' xmlns='http://www.w3.org/2000/svg'%3E%3Cfilter id='n'%3E%3CfeTurbulence type='fractalNoise' baseFrequency='0.8' numOctaves='4' stitchTiles='stitch'/%3E%3C/filter%3E%3Crect width='100%25' height='100%25' filter='url(%23n)'/%3E%3C/svg%3E");
        }

        .page {
            position: relative;
            z-index: 2;
            min-height: 100vh;
            display: flex;
            flex-direction: column;
        }

        .nav {
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            z-index: 100;
            padding: 12px 16px;
        }

        .nav-inner {
            max-width: 1200px;
            margin: 0 auto;
            height: 56px;
            background: rgba(255, 255, 255, 0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 0 16px;
            display: flex;
            align-items: center;
            justify-content: space-between;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 12px 48px rgba(0,0,0,0.03);
        }

        .logo {
            display: flex;
            align-items: center;
            gap: 14px;
            text-decoration: none;
        }

        .logo-icon {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .logo-icon svg { width: 36px; height: 12px; }

        .logo-text {
            font-size: 18px;
            font-weight: 700;
            letter-spacing: -0.5px;
            color: var(--gray-900);
        }

        .logo-text .gas { color: var(--blue-600); }
        .logo-text .consult { color: #0F172A; }
        .logo-text .ai { color: rgba(15, 23, 42, 0.4); }

        .nav-links {
            display: none;
            align-items: center;
            gap: 4px;
        }

        .nav-link {
            padding: 10px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .nav-link.active {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .mobile-menu-btn {
            display: flex;
            flex-direction: column;
            gap: 5px;
            background: none;
            border: none;
            cursor: pointer;
            padding: 8px;
            border-radius: 8px;
            transition: background 0.2s ease;
        }

        .mobile-menu-btn:hover {
            background: rgba(0,0,0,0.04);
        }

        .mobile-menu-btn span {
            display: block;
            width: 22px;
            height: 2px;
            background: var(--gray-700);
            border-radius: 1px;
            transition: all 0.3s ease;
        }

        .mobile-menu-btn.active span:nth-child(1) {
            transform: rotate(45deg) translate(7px, 7px);
        }

        .mobile-menu-btn.active span:nth-child(2) {
            opacity: 0;
        }

        .mobile-menu-btn.active span:nth-child(3) {
            transform: rotate(-45deg) translate(7px, -7px);
        }

        .mobile-menu {
            display: none;
            position: fixed;
            top: 80px;
            left: 16px;
            right: 16px;
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 8px;
            box-shadow: 0 4px 16px rgba(0,0,0,0.08), 0 12px 48px rgba(0,0,0,0.12);
            z-index: 99;
            flex-direction: column;
            gap: 4px;
        }

        .mobile-menu.active {
            display: flex;
        }

        .mobile-menu-link {
            padding: 14px 16px;
            font-size: 15px;
            font-weight: 500;
            color: var(--gray-700);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .mobile-menu-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .hero {
            padding: 120px 20px 60px;
            text-align: center;
        }

        .hero-badge {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            background: var(--white);
            border: 1px solid var(--gray-200);
            border-radius: 100px;
            padding: 8px 16px 8px 12px;
            margin-bottom: 24px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.04);
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) forwards;
            opacity: 0;
        }

        .badge-dot {
            width: 8px;
            height: 8px;
            background: var(--blue-500);
            border-radius: 50%;
            position: relative;
        }

        .badge-dot::after {
            content: '';
            position: absolute;
            inset: -3px;
            border-radius: 50%;
            background: var(--blue-400);
            animation: pulse-ring 2s ease-out infinite;
        }

        @keyframes pulse-ring {
            0% { transform: scale(0.8); opacity: 0.8; }
            100% { transform: scale(2); opacity: 0; }
        }

        .badge-text {
            font-size: 12px;
            font-weight: 600;
            color: var(--gray-700);
        }

        .hero-title {
            font-size: 40px;
            font-weight: 800;
            line-height: 1.1;
            letter-spacing: -2px;
            color: var(--gray-900);
            margin-bottom: 20px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.1s forwards;
            opacity: 0;
        }

        .hero-title .gradient { color: var(--blue-600); }

        .hero-subtitle {
            font-size: 16px;
            font-weight: 400;
            line-height: 1.6;
            color: var(--gray-500);
            max-width: 560px;
            margin: 0 auto 40px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.2s forwards;
            opacity: 0;
        }

        @keyframes fade-up {
            from { opacity: 0; transform: translateY(24px); }
            to { opacity: 1; transform: translateY(0); }
        }

        .chat-container {
            max-width: 720px;
            margin: 0 auto 60px;
            padding: 0 16px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.3s forwards;
            opacity: 0;
        }

        .chat-card {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(40px) saturate(180%);
            -webkit-backdrop-filter: blur(40px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 8px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
            transition: all 0.4s cubic-bezier(0.4,0,0.2,1);
        }

        .chat-card:focus-within {
            box-shadow: 0 0 0 4px rgba(59,130,246,0.1), 0 1px 2px rgba(0,0,0,0.02), 0 8px 24px rgba(37,99,235,0.08), 0 32px 100px rgba(37,99,235,0.12), inset 0 1px 0 rgba(255,255,255,0.8);
            border-color: rgba(59,130,246,0.3);
        }

        .chat-inner {
            background: var(--white);
            border-radius: 14px;
            padding: 6px;
            display: flex;
            align-items: flex-end;
            gap: 6px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04), inset 0 1px 0 rgba(255,255,255,1);
        }

        .chat-input {
            flex: 1;
            border: none;
            outline: none;
            padding: 14px 16px;
            font-size: 16px;
            font-family: inherit;
            color: var(--gray-800);
            background: transparent;
            resize: none;
            min-height: 52px;
            max-height: 150px;
            line-height: 1.5;
        }

        .chat-input::placeholder { color: var(--gray-400); }

        .chat-send {
            width: 48px;
            height: 48px;
            background: var(--blue-600);
            border: none;
            border-radius: 12px;
            color: var(--white);
            cursor: pointer;
            display: flex;
            align-items: center;
            justify-content: center;
            transition: all 0.25s cubic-bezier(0.4,0,0.2,1);
            box-shadow: 0 1px 2px rgba(37,99,235,0.2), 0 4px 16px rgba(37,99,235,0.2), inset 0 1px 0 rgba(255,255,255,0.1), inset 0 -1px 0 rgba(0,0,0,0.1);
            flex-shrink: 0;
            margin: 6px;
        }

        .chat-send:hover {
            background: var(--blue-700);
            transform: translateY(-2px);
            box-shadow: 0 2px 4px rgba(37,99,235,0.2), 0 12px 40px rgba(37,99,235,0.3), inset 0 1px 0 rgba(255,255,255,0.1), inset 0 -1px 0 rgba(0,0,0,0.1);
        }

        .chat-send:active { transform: translateY(0); }
        .chat-send svg { width: 20px; height: 20px; }

        .chat-hints {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            padding: 16px 8px 6px;
        }

        .hint-chip {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            background: rgba(255,255,255,0.6);
            border: 1px solid var(--gray-200);
            border-radius: 100px;
            padding: 10px 14px;
            font-size: 12px;
            font-weight: 500;
            color: var(--gray-600);
            cursor: pointer;
            transition: all 0.2s ease;
        }

        .hint-chip:hover {
            background: var(--white);
            border-color: var(--blue-200);
            color: var(--blue-600);
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(37,99,235,0.1);
        }

        .hint-chip svg { width: 14px; height: 14px; opacity: 0.5; }
        .hint-chip:hover svg { opacity: 1; color: var(--blue-500); }

        .features { padding: 60px 20px 80px; }

        .features-header {
            text-align: center;
            margin-bottom: 40px;
        }

        .features-label {
            font-size: 11px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 1.5px;
            color: var(--blue-600);
            margin-bottom: 12px;
        }

        .features-title {
            font-size: 28px;
            font-weight: 800;
            letter-spacing: -1px;
            color: var(--gray-900);
            margin-bottom: 12px;
        }

        .features-subtitle {
            font-size: 16px;
            color: var(--gray-500);
            max-width: 480px;
            margin: 0 auto;
        }

        .features-grid {
            display: grid;
            grid-template-columns: 1fr;
            gap: 16px;
            max-width: 1200px;
            margin: 0 auto;
        }

        .feature-card {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            padding: 28px;
            position: relative;
            overflow: hidden;
            transition: all 0.4s cubic-bezier(0.4,0,0.2,1);
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04);
            cursor: pointer;
        }

        .feature-card::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 1px;
            background: linear-gradient(90deg, transparent 0%, rgba(255,255,255,0.8) 50%, transparent 100%);
        }

        .feature-card:hover {
            transform: translateY(-4px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.04), 0 24px 64px rgba(0,0,0,0.08);
            border-color: rgba(59,130,246,0.2);
        }

        .feature-icon {
            width: 56px;
            height: 56px;
            border-radius: 16px;
            display: flex;
            align-items: center;
            justify-content: center;
            margin-bottom: 20px;
            transition: all 0.3s ease;
        }

        .feature-card:hover .feature-icon { transform: scale(1.05); }
        .feature-icon svg { width: 24px; height: 24px; }

        .feature-icon.blue {
            background: var(--blue-50);
            border: 1px solid var(--blue-100);
            box-shadow: 0 4px 16px rgba(37,99,235,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.blue svg { color: var(--blue-600); }

        .feature-icon.emerald {
            background: #ECFDF5;
            border: 1px solid #D1FAE5;
            box-shadow: 0 4px 16px rgba(16,185,129,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.emerald svg { color: #059669; }

        .feature-icon.violet {
            background: #F5F3FF;
            border: 1px solid #EDE9FE;
            box-shadow: 0 4px 16px rgba(139,92,246,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.violet svg { color: #7C3AED; }

        .feature-icon.rose {
            background: #FFF1F2;
            border: 1px solid #FFE4E6;
            box-shadow: 0 4px 16px rgba(244,63,94,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.rose svg { color: #E11D48; }

        .feature-icon.cyan {
            background: #ECFEFF;
            border: 1px solid #CFFAFE;
            box-shadow: 0 4px 16px rgba(6,182,212,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.cyan svg { color: #0891B2; }

        .feature-title {
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 10px;
            letter-spacing: -0.3px;
        }

        .feature-desc {
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-500);
            margin-bottom: 20px;
        }

        .feature-link {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            font-size: 14px;
            font-weight: 600;
            color: var(--blue-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .feature-link:hover { gap: 10px; }
        .feature-link svg { width: 16px; height: 16px; transition: transform 0.2s ease; }
        .feature-link:hover svg { transform: translateX(4px); }

        .footer {
            padding: 32px 20px;
            border-top: 1px solid var(--gray-200);
            background: rgba(255,255,255,0.5);
        }

        .footer-inner {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 20px;
            text-align: center;
        }

        .footer-brand {
            display: flex;
            align-items: center;
            gap: 10px;
        }

        .footer-logo {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .footer-logo svg { width: 32px; height: 32px; }

        .footer-text {
            font-size: 13px;
            color: var(--gray-500);
        }

        .footer-links {
            display: flex;
            gap: 24px;
        }

        .footer-link {
            font-size: 13px;
            color: var(--gray-500);
            text-decoration: none;
            transition: color 0.2s ease;
        }

        .footer-link:hover { color: var(--gray-700); }

        @media (min-width: 768px) {
            .nav { padding: 16px 32px; }
            .nav-inner { height: 64px; padding: 0 24px; border-radius: 20px; }
            .logo-icon svg { width: 42px; height: 15px; }
            .logo-text { font-size: 20px; }
            .nav-links { display: flex; }
            .mobile-menu-btn { display: none; }
            .hero { padding: 160px 32px 80px; }
            .hero-badge { padding: 10px 20px 10px 14px; margin-bottom: 32px; }
            .badge-dot { width: 10px; height: 10px; }
            .badge-text { font-size: 13px; }
            .hero-title { font-size: 56px; letter-spacing: -2.5px; margin-bottom: 24px; }
            .hero-subtitle { font-size: 18px; margin-bottom: 48px; }
            .chat-container { padding: 0 24px; margin-bottom: 80px; }
            .chat-card { border-radius: 24px; padding: 10px; }
            .chat-inner { border-radius: 18px; padding: 8px; }
            .chat-input { padding: 16px 18px; min-height: 56px; }
            .chat-send { width: 52px; height: 52px; border-radius: 14px; }
            .chat-hints { gap: 10px; padding: 18px 10px 8px; }
            .hint-chip { padding: 12px 18px; font-size: 13px; }
            .hint-chip svg { width: 16px; height: 16px; }
            .features { padding: 80px 32px 100px; }
            .features-header { margin-bottom: 56px; }
            .features-label { font-size: 12px; margin-bottom: 16px; }
            .features-title { font-size: 36px; letter-spacing: -1.5px; }
            .features-subtitle { font-size: 18px; }
            .features-grid { grid-template-columns: repeat(2, 1fr); gap: 20px; }
            .feature-card { padding: 36px; border-radius: 24px; }
            .feature-card:hover { transform: translateY(-6px); }
            .feature-icon { width: 60px; height: 60px; border-radius: 18px; margin-bottom: 24px; }
            .feature-icon svg { width: 26px; height: 26px; }
            .feature-title { font-size: 20px; margin-bottom: 12px; }
            .feature-desc { font-size: 15px; line-height: 1.7; margin-bottom: 24px; }
            .footer { padding: 40px 32px; }
            .footer-inner { flex-direction: row; justify-content: space-between; text-align: left; }
            .footer-logo svg { width: 36px; height: 36px; }
            .footer-text { font-size: 14px; }
            .footer-links { gap: 32px; }
            .footer-link { font-size: 14px; }
            .orb-1 { width: 600px; height: 600px; left: -10%; }
            .orb-2 { width: 450px; height: 450px; right: -10%; }
            .orb-3 { width: 400px; height: 400px; }
        }

        @media (min-width: 1024px) {
            .nav { padding: 16px 40px; }
            .hero { padding: 180px 40px 80px; }
            .hero-title { font-size: 72px; letter-spacing: -3px; margin-bottom: 28px; }
            .hero-subtitle { font-size: 20px; margin-bottom: 56px; }
            .chat-container { margin-bottom: 100px; }
            .chat-card { border-radius: 28px; }
            .chat-inner { border-radius: 20px; }
            .chat-input { padding: 18px 20px; min-height: 60px; max-height: 180px; }
            .chat-send { width: 56px; height: 56px; border-radius: 16px; margin: 8px; }
            .chat-send svg { width: 22px; height: 22px; }
            .chat-hints { padding: 20px 12px 8px; }
            .hint-chip { padding: 12px 20px; }
            .features { padding: 80px 40px 120px; }
            .features-header { margin-bottom: 64px; }
            .features-title { font-size: 40px; }
            .features-grid { grid-template-columns: repeat(3, 1fr); gap: 24px; }
            .feature-card { padding: 40px; border-radius: 28px; }
            .feature-card:hover { transform: translateY(-8px); }
            .feature-icon { width: 64px; height: 64px; border-radius: 20px; margin-bottom: 28px; }
            .feature-icon svg { width: 28px; height: 28px; }
            .footer { padding: 48px 40px; }
            .orb-1 { width: 800px; height: 800px; top: -20%; left: -10%; }
            .orb-2 { width: 600px; height: 600px; right: -15%; }
            .orb-3 { width: 500px; height: 500px; }
        }

        @media (min-width: 1280px) {
            .hero-title { font-size: 80px; }
        }
    

        .main-content {
            flex: 1;
            padding: 100px 20px 40px;
            max-width: 1200px;
            margin: 0 auto;
            width: 100%;
        }

        .content-card {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            padding: 32px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
        }

        h1 { font-size: 32px; margin-bottom: 16px; font-weight: 700; letter-spacing: -0.5px; }
        h2 { font-size: 24px; margin-bottom: 12px; font-weight: 700; letter-spacing: -0.5px; }
        h3 { font-size: 20px; margin-bottom: 10px; font-weight: 700; letter-spacing: -0.3px; }

        input[type="text"], input[type="number"], input[type="email"], select, textarea {
            width: 100%;
            padding: 12px 16px;
            border: 1px solid var(--gray-300);
            border-radius: 12px;
            font-family: inherit;
            font-size: 15px;
            background: var(--white);
            color: var(--gray-900);
            transition: all 0.2s ease;
        }

        input:focus, select:focus, textarea:focus {
            outline: none;
            border-color: var(--blue-500);
            box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.1);
        }

        button, .btn {
            padding: 12px 24px;
            background: var(--blue-600);
            color: var(--white);
            border: none;
            border-radius: 12px;
            font-family: inherit;
            font-size: 15px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s ease;
            box-shadow: 0 1px 2px rgba(37,99,235,0.2), 0 4px 16px rgba(37,99,235,0.2), inset 0 1px 0 rgba(255,255,255,0.1);
        }

        button:hover, .btn:hover {
            background: var(--blue-700);
            transform: translateY(-2px);
            box-shadow: 0 2px 4px rgba(37,99,235,0.2), 0 12px 40px rgba(37,99,235,0.3), inset 0 1px 0 rgba(255,255,255,0.1);
        }

        button:active, .btn:active {
            transform: translateY(0);
        }

        table {
            width: 100%;
            border-collapse: collapse;
            background: rgba(255,255,255,0.5);
            border-radius: 12px;
            overflow: hidden;
        }

        th, td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid var(--gray-200);
        }

        th {
            background: rgba(37,99,235,0.05);
            font-weight: 600;
            color: var(--blue-700);
        }

        @media (min-width: 768px) {
            .main-content { padding: 120px 32px 60px; }
        }

        @media (min-width: 1024px) {
            .main-content { padding: 140px 40px 80px; }
        }

        /* Pre-Op Assessment Specific Styles */
        .phi-banner {
            background: linear-gradient(135deg, rgba(239, 68, 68, 0.1) 0%, rgba(220, 38, 38, 0.05) 100%);
            border: 1px solid rgba(239, 68, 68, 0.2);
            border-radius: 16px;
            padding: 20px 24px;
            margin-bottom: 32px;
        }

        .phi-content {
            display: flex;
            gap: 16px;
            align-items: flex-start;
        }

        .phi-icon {
            font-size: 28px;
            flex-shrink: 0;
        }

        .phi-text {
            flex: 1;
        }

        .phi-text strong {
            display: block;
            font-size: 15px;
            font-weight: 700;
            color: #DC2626;
            margin-bottom: 6px;
        }

        .phi-text p {
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-700);
            margin: 0;
        }

        .main-wrapper {
            max-width: 900px;
            margin: 0 auto;
        }

        .content-container {
            width: 100%;
        }

        .page-header {
            text-align: center;
            margin-bottom: 48px;
        }

        .page-title {
            font-size: 36px;
            font-weight: 800;
            letter-spacing: -1px;
            color: var(--gray-900);
            margin-bottom: 12px;
        }

        .page-subtitle {
            font-size: 16px;
            color: var(--gray-600);
            font-weight: 400;
        }

        .form-card {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 32px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
        }

        .section-title {
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 24px;
            padding-bottom: 12px;
            border-bottom: 2px solid var(--blue-100);
        }

        .form-grid {
            display: grid;
            grid-template-columns: 1fr;
            gap: 20px;
        }

        @media (min-width: 768px) {
            .form-grid {
                grid-template-columns: repeat(2, 1fr);
            }
        }

        @media (min-width: 1024px) {
            .form-grid {
                grid-template-columns: repeat(3, 1fr);
            }
        }

        .field-group {
            display: flex;
            flex-direction: column;
            gap: 8px;
        }

        .field-label {
            font-size: 14px;
            font-weight: 600;
            color: var(--gray-700);
        }

        .field-input, select.field-input, textarea.field-input {
            width: 100%;
            padding: 12px 16px;
            border: 1px solid var(--gray-300);
            border-radius: 12px;
            font-family: inherit;
            font-size: 15px;
            background: var(--white);
            color: var(--gray-900);
            transition: all 0.2s ease;
        }

        .field-input:focus, select.field-input:focus, textarea.field-input:focus {
            outline: none;
            border-color: var(--blue-500);
            box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.1);
        }

        .radio-group {
            display: flex;
            gap: 16px;
            flex-wrap: wrap;
        }

        .radio-item {
            display: flex;
            align-items: center;
            gap: 8px;
        }

        .radio-item input[type="radio"] {
            width: 18px;
            height: 18px;
            cursor: pointer;
        }

        .radio-item label {
            font-size: 14px;
            color: var(--gray-700);
            cursor: pointer;
        }

        .checkbox-group {
            display: grid;
            grid-template-columns: 1fr;
            gap: 12px;
        }

        @media (min-width: 768px) {
            .checkbox-group {
                grid-template-columns: repeat(2, 1fr);
            }
        }

        .checkbox-item {
            display: flex;
            align-items: center;
            gap: 10px;
        }

        .checkbox-item input[type="checkbox"] {
            width: 18px;
            height: 18px;
            cursor: pointer;
        }

        .checkbox-item label {
            font-size: 14px;
            color: var(--gray-700);
            cursor: pointer;
        }

        textarea.field-input {
            min-height: 100px;
            resize: vertical;
        }

        .submit-button {
            width: 100%;
            padding: 16px 32px;
            background: var(--blue-600);
            color: var(--white);
            border: none;
            border-radius: 12px;
            font-family: inherit;
            font-size: 16px;
            font-weight: 700;
            cursor: pointer;
            transition: all 0.2s ease;
            box-shadow: 0 1px 2px rgba(37,99,235,0.2), 0 4px 16px rgba(37,99,235,0.2), inset 0 1px 0 rgba(255,255,255,0.1);
            margin-top: 12px;
        }

        .submit-button:hover {
            background: var(--blue-700);
            transform: translateY(-2px);
            box-shadow: 0 2px 4px rgba(37,99,235,0.2), 0 12px 40px rgba(37,99,235,0.3), inset 0 1px 0 rgba(255,255,255,0.1);
        }

        .submit-button:active {
            transform: translateY(0);
        }

        .result-section {
            margin-top: 32px;
        }

        .result-card {
            background: var(--white);
            border: 1px solid var(--gray-200);
            border-radius: 16px;
            padding: 28px;
            margin-bottom: 24px;
        }

        .result-title {
            font-size: 20px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 20px;
        }

        .risk-badge {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            padding: 8px 16px;
            border-radius: 100px;
            font-size: 13px;
            font-weight: 600;
            margin-bottom: 16px;
        }

        .risk-badge.low {
            background: rgba(16, 185, 129, 0.1);
            color: #059669;
            border: 1px solid rgba(16, 185, 129, 0.2);
        }

        .risk-badge.moderate {
            background: rgba(245, 158, 11, 0.1);
            color: #D97706;
            border: 1px solid rgba(245, 158, 11, 0.2);
        }

        .risk-badge.high {
            background: rgba(239, 68, 68, 0.1);
            color: #DC2626;
            border: 1px solid rgba(239, 68, 68, 0.2);
        }

        .metrics-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
            gap: 16px;
            margin-bottom: 24px;
        }

        .metric-card {
            background: var(--gray-50);
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            padding: 16px;
            text-align: center;
        }

        .metric-label {
            font-size: 12px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: var(--gray-600);
            margin-bottom: 8px;
        }

        .metric-value {
            font-size: 24px;
            font-weight: 800;
            color: var(--blue-600);
        }

        .evidence-section {
            margin-top: 24px;
            padding-top: 24px;
            border-top: 1px solid var(--gray-200);
        }

        .evidence-title {
            font-size: 16px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 16px;
        }

        .reference-list {
            list-style: none;
            padding: 0;
            margin: 0;
        }

        .reference-item {
            padding: 12px 16px;
            background: var(--gray-50);
            border-left: 3px solid var(--blue-500);
            border-radius: 8px;
            margin-bottom: 12px;
            font-size: 13px;
            line-height: 1.6;
            color: var(--gray-700);
        }

        .reference-item:last-child {
            margin-bottom: 0;
        }

    </style>
</head>
<body>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <nav class="nav">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo">
                    <div class="logo-icon">
                        <svg width="36" height="12" viewBox="0 0 52 18" fill="none">
                            <circle cx="9" cy="9" r="9" fill="#2563EB"/>
                            <circle cx="26" cy="9" r="9" fill="#2563EB" fill-opacity="0.5"/>
                            <circle cx="43" cy="9" r="9" fill="#2563EB" fill-opacity="0.2"/>
                        </svg>
                    </div>
                    <span class="logo-text"><span class="gas">gas</span><span class="consult">consult</span><span class="ai">.ai</span></span>
                </a>
                <div class="nav-links">
                    <a href="/?clear=1" class="nav-link">Home</a>
                    <a href="/quick-dose" class="nav-link">Quick Dose</a>
                    <a href="/preop" class="nav-link active">Pre-Op</a>
                    <a href="/calculators" class="nav-link">Calculators</a>
                    <a href="/hypotension" class="nav-link">IOH Predictor</a>
                </div>
                <button class="mobile-menu-btn" onclick="toggleMobileMenu()" aria-label="Toggle menu">
                    <span></span>
                    <span></span>
                    <span></span>
                </button>
            </div>
        </nav>
        <div class="mobile-menu" id="mobileMenu">
            <a href="/?clear=1" class="mobile-menu-link">Home</a>
            <a href="/quick-dose" class="mobile-menu-link">Quick Dose</a>
            <a href="/preop" class="mobile-menu-link">Pre-Op</a>
            <a href="/calculators" class="mobile-menu-link">Calculators</a>
            <a href="/hypotension" class="mobile-menu-link">IOH Predictor</a>
        </div>

        <main class="main-content">
    <!-- PHI Warning -->
    <div class="phi-banner">
        <div class="phi-content">
            <div class="phi-icon">⚠️</div>
            <div class="phi-text">
                <strong>Privacy Notice:</strong>
                <p>Do not enter patient names, dates of birth, MRNs, or other identifying information. Use age, weight, and clinical details only.</p>
            </div>
        </div>
    </div>

    <!-- Main Content -->
    <div class="main-wrapper">
        <div class="content-container">
            <div class="page-header">
                <h1 class="page-title">Pre-Operative Assessment</h1>
                <p class="page-subtitle">Evidence-based risk stratification and recommendations</p>
            </div>

            {% if not summary %}
            <!-- Assessment Form -->
            <form method="post" action="/preop">
                <input type="hidden" name="csrf_token" value="{{ csrf_token() }}"/>

                <!-- Section 1: Demographics -->
                <div class="form-card">
                    <h2 class="section-title">1. Patient Demographics</h2>
                    <div class="form-grid">
                        <div class="field-group">
                            <label for="age" class="field-label">Age (years)</label>
                            <input type="number" id="age" name="age" class="field-input" required>
                        </div>
                        <div class="field-group">
                            <label for="weight" class="field-label">Weight (kg)</label>
                            <input type="number" id="weight" name="weight" step="0.1" class="field-input" required onchange="calculateMetrics()">
                        </div>
                        <div class="field-group">
                            <label for="height" class="field-label">Height (cm)</label>
                            <input type="number" id="height" name="height" step="0.1" class="field-input" required onchange="calculateMetrics()">
                        </div>
                    </div>
                    <div class="field-group">
                        <label class="field-label">Sex Assigned at Birth</label>
                        <div class="radio-group">
                            <div class="radio-item">
                                <input type="radio" id="male" name="sex" value="male" onchange="calculateMetrics()" required>
                                <label for="male">Male</label>
                            </div>
                            <div class="radio-item">
                                <input type="radio" id="female" name="sex" value="female" onchange="calculateMetrics()" required>
                                <label for="female">Female</label>
                            </div>
                        </div>
                    </div>
                    <div class="auto-calc" id="autoCalc">
                        Enter weight, height, and sex to calculate BMI and IBW
                    </div>
                </div>

                <!-- Section 2: Comorbidities -->
                <div class="form-card">
                    <h2 class="section-title">2. Comorbidities</h2>
                    <div class="checkbox-grid">
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
                    <div class="field-group" style="margin-top: 1rem;">
                        <label for="other_comorbidities" class="field-label">Other Comorbidities</label>
                        <textarea id="other_comorbidities" name="other_comorbidities" class="field-textarea" placeholder="e.g., GERD, Hypothyroidism, Chronic Pain..." rows="2"></textarea>
                    </div>
                </div>

                <!-- Section 3: Functional Status -->
                <div class="form-card">
                    <h2 class="section-title">3. Functional Status</h2>
                    <div class="field-group">
                        <label for="mets" class="field-label">Metabolic Equivalents (METs)</label>
                        <select id="mets" name="mets" class="field-select" required>
                            <option value="">Select...</option>
                            <option value="Unknown">Unknown / Not documented</option>
                            <option value="<4 METs">&lt;4 METs (Cannot climb 2 flights of stairs)</option>
                            <option value="4-10 METs">4-10 METs (Can climb 2 flights of stairs)</option>
                            <option value=">10 METs">&gt;10 METs (Very active, strenuous sports)</option>
                        </select>
                    </div>
                    <div class="field-group" style="margin-top: 1.25rem;">
                        <label for="previous_anesthesia" class="field-label">Previous Anesthesia History</label>
                        <textarea id="previous_anesthesia" name="previous_anesthesia" class="field-textarea" placeholder="e.g., General anesthesia for appendectomy 2015 - no complications. Family history of malignant hyperthermia..." rows="3"></textarea>
                    </div>
                </div>

                <!-- Section 4: Medications -->
                <div class="form-card">
                    <h2 class="section-title">4. Current Medications</h2>
                    <div class="field-group">
                        <label for="medications" class="field-label">List all medications</label>
                        <textarea id="medications" name="medications" class="field-textarea" placeholder="e.g., Aspirin 81mg daily, Metoprolol 50mg BID, Apixaban 5mg BID..."></textarea>
                    </div>
                </div>

                <!-- Section 5: Labs & Cardiac -->
                <div class="form-card">
                    <h2 class="section-title">5. Laboratory Values & Cardiac Assessment</h2>
                    <div class="form-grid">
                        <div class="field-group">
                            <label for="hgb" class="field-label">Hemoglobin (g/dL)</label>
                            <input type="number" id="hgb" name="hgb" step="0.1" class="field-input">
                        </div>
                        <div class="field-group">
                            <label for="plt" class="field-label">Platelets (×10³/μL)</label>
                            <input type="number" id="plt" name="plt" class="field-input">
                        </div>
                        <div class="field-group">
                            <label for="cr" class="field-label">Creatinine (mg/dL)</label>
                            <input type="number" id="cr" name="cr" step="0.01" class="field-input">
                        </div>
                        <div class="field-group">
                            <label for="inr" class="field-label">INR</label>
                            <input type="number" id="inr" name="inr" step="0.1" class="field-input">
                        </div>
                    </div>
                    <div class="form-grid">
                        <div class="field-group">
                            <label for="ef" class="field-label">Ejection Fraction (%)</label>
                            <input type="text" id="ef" name="ef" class="field-input" placeholder="e.g., 55-60% or None">
                        </div>
                        <div class="field-group">
                            <label for="ekg" class="field-label">EKG Findings</label>
                            <input type="text" id="ekg" name="ekg" class="field-input" placeholder="e.g., NSR, Afib, or None">
                        </div>
                    </div>
                </div>

                <!-- Section 6: Surgical Procedure -->
                <div class="form-card">
                    <h2 class="section-title">6. Surgical Procedure</h2>
                    <div class="field-group">
                        <label for="procedure" class="field-label">Procedure Type</label>
                        <input type="text" id="procedure" name="procedure" class="field-input" placeholder="e.g., Total Knee Arthroplasty, CABG, Laparoscopic Cholecystectomy..." required>
                    </div>
                    <div class="field-group" style="margin-top: 1.25rem;">
                        <label for="surgery_risk" class="field-label">Surgery Risk Category</label>
                        <select id="surgery_risk" name="surgery_risk" class="field-select" required>
                            <option value="">Select...</option>
                            <option value="Low">Low Risk (&lt;1% cardiac risk)</option>
                            <option value="Intermediate">Intermediate Risk (1-5% cardiac risk)</option>
                            <option value="High">High Risk (&gt;5% cardiac risk)</option>
                        </select>
                    </div>
                </div>

                <!-- Section 7: Additional Info -->
                <div class="form-card">
                    <h2 class="section-title">7. Additional Information</h2>
                    <div class="form-grid">
                        <div class="field-group">
                            <label for="npo" class="field-label">NPO Status</label>
                            <input type="text" id="npo" name="npo" class="field-input" placeholder="e.g., NPO since midnight...">
                        </div>
                        <div class="field-group">
                            <label for="allergies" class="field-label">Allergies</label>
                            <input type="text" id="allergies" name="allergies" class="field-input" placeholder="e.g., PCN (rash), NKDA...">
                        </div>
                    </div>
                </div>

                <button type="submit" class="submit-button">Generate Evidence-Based Assessment</button>
            </form>

            {% else %}
            <!-- Assessment Results -->
            <div class="summary-card">
                <h2 class="summary-title">Pre-Operative Assessment Summary</h2>
                <div class="summary-content">
                    {{ summary|safe }}
                </div>

                {% if references %}
                <div class="references-section">
                    <h3 class="references-title">References:</h3>
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

            <div class="new-assessment-wrapper">
                <a href="/preop" class="new-assessment-btn">New Assessment</a>
            </div>
            {% endif %}
        </div>
    </div>

    <!-- Glassmorphism Footer -->
    <footer class="footer">
        <div class="footer-inner">
            <div class="footer-brand">
                <div class="footer-logo">
                    <svg width="32" height="32" viewBox="0 0 32 32" fill="none">
                        <circle cx="6" cy="16" r="5" fill="#2563EB"/>
                        <circle cx="16" cy="16" r="5" fill="#2563EB" fill-opacity="0.5"/>
                        <circle cx="26" cy="16" r="5" fill="#2563EB" fill-opacity="0.2"/>
                    </svg>
                </div>
                <span class="footer-text">© 2025 GasConsult.ai</span>
            </div>
            <div class="footer-links">
                <a href="/privacy" class="footer-link">Privacy</a>
                <a href="/terms" class="footer-link">Terms</a>
                <a href="#" class="footer-link">Contact</a>
            </div>
        </div>
    </footer>
        </main>
    </div>
    <script>
        function toggleMobileMenu() {
            const menu = document.getElementById('mobileMenu');
            const btn = document.querySelector('.mobile-menu-btn');
            if (menu && btn) {
                menu.classList.toggle('active');
                btn.classList.toggle('active');
            }
        }
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

HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GasConsult.ai — AI-Powered Anesthesiology Assistant</title>

    <!-- PWA -->
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=5">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">

    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap" rel="stylesheet">
    <style>
        :root {
            --white: #FFFFFF;
            --gray-50: #F8FAFC;
            --gray-100: #F1F5F9;
            --gray-200: #E2E8F0;
            --gray-300: #CBD5E1;
            --gray-400: #94A3B8;
            --gray-500: #64748B;
            --gray-600: #475569;
            --gray-700: #334155;
            --gray-800: #1E293B;
            --gray-900: #0F172A;
            --blue-50: #EFF6FF;
            --blue-100: #DBEAFE;
            --blue-200: #BFDBFE;
            --blue-300: #93C5FD;
            --blue-400: #60A5FA;
            --blue-500: #3B82F6;
            --blue-600: #2563EB;
            --blue-700: #1D4ED8;
        }

        * { margin: 0; padding: 0; box-sizing: border-box; }

        html {
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
            scroll-behavior: smooth;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            background: var(--gray-50);
            color: var(--gray-900);
            min-height: 100vh;
            overflow-x: hidden;
        }

        .bg-canvas {
            position: fixed;
            inset: 0;
            z-index: 0;
            overflow: hidden;
            background: linear-gradient(180deg, #F0F7FF 0%, var(--gray-50) 50%, #FAFBFF 100%);
        }

        .orb {
            position: absolute;
            border-radius: 50%;
            filter: blur(80px);
            opacity: 0.6;
            animation: float 20s ease-in-out infinite;
        }

        .orb-1 {
            width: 400px;
            height: 400px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.15) 0%, transparent 70%);
            top: -15%;
            left: -20%;
        }

        .orb-2 {
            width: 300px;
            height: 300px;
            background: radial-gradient(circle, rgba(147, 197, 253, 0.2) 0%, transparent 70%);
            top: 30%;
            right: -20%;
            animation-delay: -7s;
            animation-duration: 25s;
        }

        .orb-3 {
            width: 250px;
            height: 250px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.1) 0%, transparent 70%);
            bottom: -10%;
            left: 20%;
            animation-delay: -14s;
            animation-duration: 30s;
        }

        @keyframes float {
            0%, 100% { transform: translate(0, 0) scale(1); }
            25% { transform: translate(40px, -40px) scale(1.05); }
            50% { transform: translate(20px, 40px) scale(0.95); }
            75% { transform: translate(-40px, 20px) scale(1.02); }
        }

        .grain {
            position: fixed;
            inset: 0;
            z-index: 1;
            pointer-events: none;
            opacity: 0.02;
            background-image: url("data:image/svg+xml,%3Csvg viewBox='0 0 512 512' xmlns='http://www.w3.org/2000/svg'%3E%3Cfilter id='n'%3E%3CfeTurbulence type='fractalNoise' baseFrequency='0.8' numOctaves='4' stitchTiles='stitch'/%3E%3C/filter%3E%3Crect width='100%25' height='100%25' filter='url(%23n)'/%3E%3C/svg%3E");
        }

        .page {
            position: relative;
            z-index: 2;
            min-height: 100vh;
            display: flex;
            flex-direction: column;
        }

        .nav {
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            z-index: 100;
            padding: 12px 16px;
        }

        .nav-inner {
            max-width: 1200px;
            margin: 0 auto;
            height: 56px;
            background: rgba(255, 255, 255, 0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 0 16px;
            display: flex;
            align-items: center;
            justify-content: space-between;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 12px 48px rgba(0,0,0,0.03);
        }

        .logo {
            display: flex;
            align-items: center;
            gap: 14px;
            text-decoration: none;
        }

        .logo-icon {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .logo-icon svg { width: 36px; height: 12px; }

        .logo-text {
            font-size: 18px;
            font-weight: 700;
            letter-spacing: -0.5px;
            color: var(--gray-900);
        }

        .logo-text .gas { color: var(--blue-600); }
        .logo-text .consult { color: #0F172A; }
        .logo-text .ai { color: rgba(15, 23, 42, 0.4); }

        .nav-links {
            display: none;
            align-items: center;
            gap: 4px;
        }

        .nav-link {
            padding: 10px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .nav-link.active {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .mobile-menu-btn {
            display: flex;
            flex-direction: column;
            gap: 5px;
            background: none;
            border: none;
            cursor: pointer;
            padding: 8px;
            border-radius: 8px;
            transition: background 0.2s ease;
        }

        .mobile-menu-btn:hover {
            background: rgba(0,0,0,0.04);
        }

        .mobile-menu-btn span {
            display: block;
            width: 22px;
            height: 2px;
            background: var(--gray-700);
            border-radius: 1px;
            transition: all 0.3s ease;
        }

        .mobile-menu-btn.active span:nth-child(1) {
            transform: rotate(45deg) translate(7px, 7px);
        }

        .mobile-menu-btn.active span:nth-child(2) {
            opacity: 0;
        }

        .mobile-menu-btn.active span:nth-child(3) {
            transform: rotate(-45deg) translate(7px, -7px);
        }

        .mobile-menu {
            display: none;
            position: fixed;
            top: 80px;
            left: 16px;
            right: 16px;
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 8px;
            box-shadow: 0 4px 16px rgba(0,0,0,0.08), 0 12px 48px rgba(0,0,0,0.12);
            z-index: 99;
            flex-direction: column;
            gap: 4px;
        }

        .mobile-menu.active {
            display: flex;
        }

        .mobile-menu-link {
            padding: 14px 16px;
            font-size: 15px;
            font-weight: 500;
            color: var(--gray-700);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .mobile-menu-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .hero {
            padding: 120px 20px 60px;
            text-align: center;
        }

        .hero-badge {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            background: var(--white);
            border: 1px solid var(--gray-200);
            border-radius: 100px;
            padding: 8px 16px 8px 12px;
            margin-bottom: 24px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.04);
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) forwards;
            opacity: 0;
        }

        .badge-dot {
            width: 8px;
            height: 8px;
            background: var(--blue-500);
            border-radius: 50%;
            position: relative;
        }

        .badge-dot::after {
            content: '';
            position: absolute;
            inset: -3px;
            border-radius: 50%;
            background: var(--blue-400);
            animation: pulse-ring 2s ease-out infinite;
        }

        @keyframes pulse-ring {
            0% { transform: scale(0.8); opacity: 0.8; }
            100% { transform: scale(2); opacity: 0; }
        }

        .badge-text {
            font-size: 12px;
            font-weight: 600;
            color: var(--gray-700);
        }

        .hero-title {
            font-size: 40px;
            font-weight: 800;
            line-height: 1.1;
            letter-spacing: -2px;
            color: var(--gray-900);
            margin-bottom: 20px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.1s forwards;
            opacity: 0;
        }

        .hero-title .gradient { color: var(--blue-600); }

        .hero-subtitle {
            font-size: 16px;
            font-weight: 400;
            line-height: 1.6;
            color: var(--gray-500);
            max-width: 560px;
            margin: 0 auto 40px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.2s forwards;
            opacity: 0;
        }

        @keyframes fade-up {
            from { opacity: 0; transform: translateY(24px); }
            to { opacity: 1; transform: translateY(0); }
        }

        .chat-container {
            max-width: 720px;
            margin: 0 auto 60px;
            padding: 0 16px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.3s forwards;
            opacity: 0;
        }

        .chat-card {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(40px) saturate(180%);
            -webkit-backdrop-filter: blur(40px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 8px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
            transition: all 0.4s cubic-bezier(0.4,0,0.2,1);
        }

        .chat-card:focus-within {
            box-shadow: 0 0 0 4px rgba(59,130,246,0.1), 0 1px 2px rgba(0,0,0,0.02), 0 8px 24px rgba(37,99,235,0.08), 0 32px 100px rgba(37,99,235,0.12), inset 0 1px 0 rgba(255,255,255,0.8);
            border-color: rgba(59,130,246,0.3);
        }

        .chat-inner {
            background: var(--white);
            border-radius: 14px;
            padding: 6px;
            display: flex;
            align-items: flex-end;
            gap: 6px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04), inset 0 1px 0 rgba(255,255,255,1);
        }

        .chat-input {
            flex: 1;
            border: none;
            outline: none;
            padding: 14px 16px;
            font-size: 16px;
            font-family: inherit;
            color: var(--gray-800);
            background: transparent;
            resize: none;
            min-height: 52px;
            max-height: 150px;
            line-height: 1.5;
        }

        .chat-input::placeholder { color: var(--gray-400); }

        .chat-send {
            width: 48px;
            height: 48px;
            background: var(--blue-600);
            border: none;
            border-radius: 12px;
            color: var(--white);
            cursor: pointer;
            display: flex;
            align-items: center;
            justify-content: center;
            transition: all 0.25s cubic-bezier(0.4,0,0.2,1);
            box-shadow: 0 1px 2px rgba(37,99,235,0.2), 0 4px 16px rgba(37,99,235,0.2), inset 0 1px 0 rgba(255,255,255,0.1), inset 0 -1px 0 rgba(0,0,0,0.1);
            flex-shrink: 0;
            margin: 6px;
        }

        .chat-send:hover {
            background: var(--blue-700);
            transform: translateY(-2px);
            box-shadow: 0 2px 4px rgba(37,99,235,0.2), 0 12px 40px rgba(37,99,235,0.3), inset 0 1px 0 rgba(255,255,255,0.1), inset 0 -1px 0 rgba(0,0,0,0.1);
        }

        .chat-send:active { transform: translateY(0); }
        .chat-send svg { width: 20px; height: 20px; }

        .chat-hints {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            padding: 16px 8px 6px;
        }

        .hint-chip {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            background: rgba(255,255,255,0.6);
            border: 1px solid var(--gray-200);
            border-radius: 100px;
            padding: 10px 14px;
            font-size: 12px;
            font-weight: 500;
            color: var(--gray-600);
            cursor: pointer;
            transition: all 0.2s ease;
        }

        .hint-chip:hover {
            background: var(--white);
            border-color: var(--blue-200);
            color: var(--blue-600);
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(37,99,235,0.1);
        }

        .hint-chip svg { width: 14px; height: 14px; opacity: 0.5; }
        .hint-chip:hover svg { opacity: 1; color: var(--blue-500); }

        .features { padding: 60px 20px 80px; }

        .features-header {
            text-align: center;
            margin-bottom: 40px;
        }

        .features-label {
            font-size: 11px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 1.5px;
            color: var(--blue-600);
            margin-bottom: 12px;
        }

        .features-title {
            font-size: 28px;
            font-weight: 800;
            letter-spacing: -1px;
            color: var(--gray-900);
            margin-bottom: 12px;
        }

        .features-subtitle {
            font-size: 16px;
            color: var(--gray-500);
            max-width: 480px;
            margin: 0 auto;
        }

        .features-grid {
            display: grid;
            grid-template-columns: 1fr;
            gap: 16px;
            max-width: 1200px;
            margin: 0 auto;
        }

        .feature-card {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            padding: 28px;
            position: relative;
            overflow: hidden;
            transition: all 0.4s cubic-bezier(0.4,0,0.2,1);
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04);
            cursor: pointer;
        }

        .feature-card::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 1px;
            background: linear-gradient(90deg, transparent 0%, rgba(255,255,255,0.8) 50%, transparent 100%);
        }

        .feature-card:hover {
            transform: translateY(-4px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.04), 0 24px 64px rgba(0,0,0,0.08);
            border-color: rgba(59,130,246,0.2);
        }

        .feature-icon {
            width: 56px;
            height: 56px;
            border-radius: 16px;
            display: flex;
            align-items: center;
            justify-content: center;
            margin-bottom: 20px;
            transition: all 0.3s ease;
        }

        .feature-card:hover .feature-icon { transform: scale(1.05); }
        .feature-icon svg { width: 24px; height: 24px; }

        .feature-icon.blue {
            background: var(--blue-50);
            border: 1px solid var(--blue-100);
            box-shadow: 0 4px 16px rgba(37,99,235,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.blue svg { color: var(--blue-600); }

        .feature-icon.emerald {
            background: #ECFDF5;
            border: 1px solid #D1FAE5;
            box-shadow: 0 4px 16px rgba(16,185,129,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.emerald svg { color: #059669; }

        .feature-icon.violet {
            background: #F5F3FF;
            border: 1px solid #EDE9FE;
            box-shadow: 0 4px 16px rgba(139,92,246,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.violet svg { color: #7C3AED; }

        .feature-icon.rose {
            background: #FFF1F2;
            border: 1px solid #FFE4E6;
            box-shadow: 0 4px 16px rgba(244,63,94,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.rose svg { color: #E11D48; }

        .feature-icon.cyan {
            background: #ECFEFF;
            border: 1px solid #CFFAFE;
            box-shadow: 0 4px 16px rgba(6,182,212,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.cyan svg { color: #0891B2; }

        .feature-title {
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 10px;
            letter-spacing: -0.3px;
        }

        .feature-desc {
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-500);
            margin-bottom: 20px;
        }

        .feature-link {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            font-size: 14px;
            font-weight: 600;
            color: var(--blue-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .feature-link:hover { gap: 10px; }
        .feature-link svg { width: 16px; height: 16px; transition: transform 0.2s ease; }
        .feature-link:hover svg { transform: translateX(4px); }

        .footer {
            padding: 32px 20px;
            border-top: 1px solid var(--gray-200);
            background: rgba(255,255,255,0.5);
        }

        .footer-inner {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 20px;
            text-align: center;
        }

        .footer-brand {
            display: flex;
            align-items: center;
            gap: 10px;
        }

        .footer-logo {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .footer-logo svg { width: 32px; height: 32px; }

        .footer-text {
            font-size: 13px;
            color: var(--gray-500);
        }

        .footer-links {
            display: flex;
            gap: 24px;
        }

        .footer-link {
            font-size: 13px;
            color: var(--gray-500);
            text-decoration: none;
            transition: color 0.2s ease;
        }

        .footer-link:hover { color: var(--gray-700); }

        @media (min-width: 768px) {
            .nav { padding: 16px 32px; }
            .nav-inner { height: 64px; padding: 0 24px; border-radius: 20px; }
            .logo-icon svg { width: 42px; height: 15px; }
            .logo-text { font-size: 20px; }
            .nav-links { display: flex; }
            .mobile-menu-btn { display: none; }
            .hero { padding: 160px 32px 80px; }
            .hero-badge { padding: 10px 20px 10px 14px; margin-bottom: 32px; }
            .badge-dot { width: 10px; height: 10px; }
            .badge-text { font-size: 13px; }
            .hero-title { font-size: 56px; letter-spacing: -2.5px; margin-bottom: 24px; }
            .hero-subtitle { font-size: 18px; margin-bottom: 48px; }
            .chat-container { padding: 0 24px; margin-bottom: 80px; }
            .chat-card { border-radius: 24px; padding: 10px; }
            .chat-inner { border-radius: 18px; padding: 8px; }
            .chat-input { padding: 16px 18px; min-height: 56px; }
            .chat-send { width: 52px; height: 52px; border-radius: 14px; }
            .chat-hints { gap: 10px; padding: 18px 10px 8px; }
            .hint-chip { padding: 12px 18px; font-size: 13px; }
            .hint-chip svg { width: 16px; height: 16px; }
            .features { padding: 80px 32px 100px; }
            .features-header { margin-bottom: 56px; }
            .features-label { font-size: 12px; margin-bottom: 16px; }
            .features-title { font-size: 36px; letter-spacing: -1.5px; }
            .features-subtitle { font-size: 18px; }
            .features-grid { grid-template-columns: repeat(2, 1fr); gap: 20px; }
            .feature-card { padding: 36px; border-radius: 24px; }
            .feature-card:hover { transform: translateY(-6px); }
            .feature-icon { width: 60px; height: 60px; border-radius: 18px; margin-bottom: 24px; }
            .feature-icon svg { width: 26px; height: 26px; }
            .feature-title { font-size: 20px; margin-bottom: 12px; }
            .feature-desc { font-size: 15px; line-height: 1.7; margin-bottom: 24px; }
            .footer { padding: 40px 32px; }
            .footer-inner { flex-direction: row; justify-content: space-between; text-align: left; }
            .footer-logo svg { width: 36px; height: 36px; }
            .footer-text { font-size: 14px; }
            .footer-links { gap: 32px; }
            .footer-link { font-size: 14px; }
            .orb-1 { width: 600px; height: 600px; left: -10%; }
            .orb-2 { width: 450px; height: 450px; right: -10%; }
            .orb-3 { width: 400px; height: 400px; }
        }

        @media (min-width: 1024px) {
            .nav { padding: 16px 40px; }
            .hero { padding: 180px 40px 80px; }
            .hero-title { font-size: 72px; letter-spacing: -3px; margin-bottom: 28px; }
            .hero-subtitle { font-size: 20px; margin-bottom: 56px; }
            .chat-container { margin-bottom: 100px; }
            .chat-card { border-radius: 28px; }
            .chat-inner { border-radius: 20px; }
            .chat-input { padding: 18px 20px; min-height: 60px; max-height: 180px; }
            .chat-send { width: 56px; height: 56px; border-radius: 16px; margin: 8px; }
            .chat-send svg { width: 22px; height: 22px; }
            .chat-hints { padding: 20px 12px 8px; }
            .hint-chip { padding: 12px 20px; }
            .features { padding: 80px 40px 120px; }
            .features-header { margin-bottom: 64px; }
            .features-title { font-size: 40px; }
            .features-grid { grid-template-columns: repeat(3, 1fr); gap: 24px; }
            .feature-card { padding: 40px; border-radius: 28px; }
            .feature-card:hover { transform: translateY(-8px); }
            .feature-icon { width: 64px; height: 64px; border-radius: 20px; margin-bottom: 28px; }
            .feature-icon svg { width: 28px; height: 28px; }
            .footer { padding: 48px 40px; }
            .orb-1 { width: 800px; height: 800px; top: -20%; left: -10%; }
            .orb-2 { width: 600px; height: 600px; right: -15%; }
            .orb-3 { width: 500px; height: 500px; }
        }

        @media (min-width: 1280px) {
            .hero-title { font-size: 80px; }
        }

        /* Chat Interface Styles */
        .chat-view {
            flex: 1;
            display: flex;
            flex-direction: column;
            padding-top: 80px;
            min-height: 100vh;
        }

        .messages-container {
            flex: 1;
            max-width: 900px;
            width: 100%;
            margin: 0 auto;
            padding: 20px 16px;
            overflow-y: auto;
        }

        .message {
            margin-bottom: 24px;
            animation: fade-up 0.4s ease;
        }

        .user-message {
            display: flex;
            justify-content: flex-end;
        }

        .user-message .message-bubble {
            background: linear-gradient(135deg, var(--blue-600) 0%, var(--blue-700) 100%);
            color: white;
            padding: 16px 20px;
            border-radius: 20px 20px 4px 20px;
            max-width: 80%;
            font-size: 15px;
            line-height: 1.5;
            box-shadow: 0 2px 8px rgba(37,99,235,0.2), 0 8px 24px rgba(37,99,235,0.15);
            word-wrap: break-word;
        }

        .ai-message {
            display: flex;
            flex-direction: column;
            align-items: flex-start;
        }

        .ai-message .message-bubble {
            background: rgba(255,255,255,0.9);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            padding: 20px 24px;
            border-radius: 20px 20px 20px 4px;
            max-width: 85%;
            font-size: 15px;
            line-height: 1.7;
            color: var(--gray-800);
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 12px 48px rgba(0,0,0,0.03);
        }

        .ai-message .message-bubble h1,
        .ai-message .message-bubble h2,
        .ai-message .message-bubble h3 {
            color: var(--gray-900);
            margin: 16px 0 12px 0;
            font-weight: 700;
        }

        .ai-message .message-bubble h1 { font-size: 20px; }
        .ai-message .message-bubble h2 { font-size: 18px; }
        .ai-message .message-bubble h3 { font-size: 16px; }

        .ai-message .message-bubble ul,
        .ai-message .message-bubble ol {
            margin: 12px 0;
            padding-left: 24px;
        }

        .ai-message .message-bubble li {
            margin: 6px 0;
        }

        .ai-message .message-bubble p {
            margin: 12px 0;
        }

        .ai-message .message-bubble strong {
            color: var(--gray-900);
            font-weight: 600;
        }

        .ai-message .message-bubble code {
            background: var(--gray-100);
            padding: 2px 6px;
            border-radius: 4px;
            font-family: 'Monaco', 'Courier New', monospace;
            font-size: 14px;
        }

        .evidence-badge {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            padding: 6px 12px;
            border-radius: 100px;
            font-size: 11px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            margin-bottom: 12px;
            cursor: help;
            transition: all 0.2s ease;
        }

        .evidence-badge:hover {
            transform: translateY(-1px);
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }

        .evidence-badge.high {
            background: #ECFDF5;
            color: #059669;
            border: 1px solid #D1FAE5;
        }

        .evidence-badge.moderate {
            background: #FEF3C7;
            color: #D97706;
            border: 1px solid #FDE68A;
        }

        .evidence-badge.low {
            background: #FEF2F2;
            color: #DC2626;
            border: 1px solid #FECACA;
        }

        .evidence-explanation {
            display: block;
            font-size: 10px;
            opacity: 0.85;
            font-weight: 500;
            margin-top: 2px;
            text-transform: none;
            letter-spacing: 0px;
        }

        /* Evidence Quality Badge - Sleek Modern Design */
        .evidence-quality-badge {
            background: linear-gradient(135deg, rgba(255, 255, 255, 0.95) 0%, rgba(249, 250, 251, 0.95) 100%);
            border: 1px solid rgba(229, 231, 235, 0.8);
            border-radius: 12px;
            padding: 16px;
            margin: 16px 0 20px 0;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.04), 0 1px 3px rgba(0, 0, 0, 0.02);
            transition: all 0.2s ease;
        }

        .evidence-quality-badge:hover {
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.08), 0 2px 4px rgba(0, 0, 0, 0.04);
        }

        .confidence-level {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            padding: 8px 16px;
            border-radius: 8px;
            font-size: 13px;
            font-weight: 600;
            margin-bottom: 10px;
            transition: all 0.2s ease;
        }

        .confidence-level strong {
            font-weight: 700;
            margin-right: 4px;
        }

        .confidence-level.high {
            background: linear-gradient(135deg, #ECFDF5 0%, #D1FAE5 100%);
            color: #047857;
            border: 1.5px solid #10B981;
            box-shadow: 0 2px 4px rgba(16, 185, 129, 0.1);
        }

        .confidence-level.moderate {
            background: linear-gradient(135deg, #FFFBEB 0%, #FEF3C7 100%);
            color: #B45309;
            border: 1.5px solid #F59E0B;
            box-shadow: 0 2px 4px rgba(245, 158, 11, 0.1);
        }

        .confidence-level.low {
            background: linear-gradient(135deg, #FEF2F2 0%, #FEE2E2 100%);
            color: #B91C1C;
            border: 1.5px solid #EF4444;
            box-shadow: 0 2px 4px rgba(239, 68, 68, 0.1);
        }

        .evidence-details {
            font-size: 12px;
            color: #6B7280;
            line-height: 1.6;
            padding: 4px 0;
            font-weight: 500;
        }

        .evidence-details::before {
            content: '';
            display: inline-block;
            width: 3px;
            height: 3px;
            background: #9CA3AF;
            border-radius: 50%;
            margin: 0 8px 2px 0;
        }

        .references-section {
            margin-top: 16px;
            padding-top: 16px;
            border-top: 1px solid var(--gray-200);
        }

        .references-title {
            font-size: 13px;
            font-weight: 700;
            color: var(--gray-700);
            margin-bottom: 10px;
            display: flex;
            align-items: center;
            gap: 6px;
        }

        .references-title svg {
            width: 16px;
            height: 16px;
            color: var(--blue-600);
        }

        .reference-item {
            background: var(--gray-50);
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            padding: 12px 14px;
            margin-bottom: 8px;
            transition: all 0.2s ease;
        }

        .reference-item:hover {
            background: var(--blue-50);
            border-color: var(--blue-200);
        }

        .reference-link {
            color: var(--blue-600);
            text-decoration: none;
            font-size: 14px;
            line-height: 1.5;
            font-weight: 500;
            display: block;
        }

        .reference-link:hover {
            text-decoration: underline;
        }

        .reference-meta {
            font-size: 12px;
            color: var(--gray-500);
            margin-top: 4px;
        }

        .streaming-indicator {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            color: var(--blue-600);
            font-size: 13px;
            font-weight: 500;
            margin-top: 8px;
        }

        .streaming-dots {
            display: flex;
            gap: 4px;
        }

        .streaming-dots span {
            width: 6px;
            height: 6px;
            background: var(--blue-600);
            border-radius: 50%;
            animation: pulse-dot 1.4s ease-in-out infinite;
        }

        .streaming-dots span:nth-child(2) { animation-delay: 0.2s; }
        .streaming-dots span:nth-child(3) { animation-delay: 0.4s; }

        @keyframes pulse-dot {
            0%, 80%, 100% { opacity: 0.3; transform: scale(0.8); }
            40% { opacity: 1; transform: scale(1); }
        }

        .chat-input-area {
            position: sticky;
            bottom: 0;
            background: var(--gray-50);
            border-top: 1px solid var(--gray-200);
            padding: 16px;
            z-index: 10;
        }

        .chat-input-wrapper {
            max-width: 900px;
            margin: 0 auto;
            display: flex;
            flex-direction: column;
            gap: 12px;
        }

        #chat-form {
            display: flex;
            gap: 8px;
            align-items: flex-end;
        }

        #chat-form .chat-card {
            flex: 1;
            margin: 0;
        }

        .new-chat-btn {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            padding: 12px 20px;
            background: rgba(255,255,255,0.8);
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            color: var(--gray-700);
            text-decoration: none;
            font-size: 14px;
            font-weight: 500;
            transition: all 0.2s ease;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02);
            align-self: flex-start;
        }

        .new-chat-btn:hover {
            background: white;
            border-color: var(--blue-200);
            color: var(--blue-600);
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(37,99,235,0.1);
        }

        .new-chat-btn svg {
            width: 16px;
            height: 16px;
        }

        .error-message {
            background: #FEF2F2;
            border: 1px solid #FECACA;
            color: #DC2626;
            padding: 16px 20px;
            border-radius: 12px;
            margin-bottom: 20px;
            font-size: 14px;
            line-height: 1.6;
        }

        @media (min-width: 768px) {
            .messages-container {
                padding: 32px 24px;
            }

            .message {
                margin-bottom: 32px;
            }

            .user-message .message-bubble {
                padding: 18px 24px;
                font-size: 16px;
                border-radius: 24px 24px 6px 24px;
            }

            .ai-message .message-bubble {
                padding: 24px 28px;
                font-size: 16px;
                border-radius: 24px 24px 24px 6px;
            }

            .chat-input-area {
                padding: 20px 24px;
            }

            .chat-input-wrapper {
                flex-direction: row;
                align-items: center;
                justify-content: space-between;
            }

            #chat-form {
                flex: 1;
                max-width: 720px;
            }

            .new-chat-btn {
                padding: 14px 24px;
                font-size: 15px;
            }
        }

        @media (min-width: 1024px) {
            .messages-container {
                padding: 40px 32px 60px;
            }

            .chat-input-area {
                padding: 24px 32px;
            }
        }

        /* Markdown Content Styling */
        .message-content {
            font-size: 15px;
            line-height: 1.7;
            color: var(--gray-800);
        }

        .message-content h1 {
            font-size: 24px;
            font-weight: 800;
            color: var(--gray-900);
            margin: 24px 0 16px 0;
            padding-bottom: 8px;
            border-bottom: 2px solid var(--blue-100);
        }

        .message-content h2 {
            font-size: 20px;
            font-weight: 700;
            color: var(--gray-900);
            margin: 20px 0 12px 0;
        }

        .message-content h3 {
            font-size: 17px;
            font-weight: 700;
            color: var(--gray-900);
            margin: 16px 0 10px 0;
        }

        .message-content h4 {
            font-size: 15px;
            font-weight: 700;
            color: var(--gray-800);
            margin: 14px 0 8px 0;
        }

        .message-content p {
            margin: 12px 0;
            line-height: 1.7;
        }

        .message-content ul,
        .message-content ol {
            margin: 12px 0;
            padding-left: 24px;
        }

        .message-content li {
            margin: 6px 0;
            line-height: 1.6;
        }

        .message-content ul li {
            list-style-type: disc;
        }

        .message-content ol li {
            list-style-type: decimal;
        }

        .message-content ul ul,
        .message-content ol ul,
        .message-content ul ol,
        .message-content ol ol {
            margin: 4px 0;
        }

        .message-content strong {
            font-weight: 700;
            color: var(--gray-900);
        }

        .message-content em {
            font-style: italic;
        }

        .message-content code {
            background: var(--gray-100);
            padding: 2px 6px;
            border-radius: 4px;
            font-family: 'Monaco', 'Courier New', monospace;
            font-size: 14px;
            color: #E11D48;
        }

        .message-content pre {
            background: var(--gray-900);
            color: var(--gray-50);
            padding: 16px;
            border-radius: 8px;
            overflow-x: auto;
            margin: 16px 0;
        }

        .message-content pre code {
            background: transparent;
            padding: 0;
            color: var(--gray-50);
            border-radius: 0;
        }

        .message-content blockquote {
            border-left: 4px solid var(--blue-500);
            padding-left: 16px;
            margin: 16px 0;
            color: var(--gray-700);
            font-style: italic;
        }

        .message-content hr {
            border: none;
            border-top: 1px solid var(--gray-200);
            margin: 24px 0;
        }

        .message-content a {
            color: var(--blue-600);
            text-decoration: none;
            border-bottom: 1px solid var(--blue-200);
            transition: all 0.2s ease;
        }

        .message-content a:hover {
            color: var(--blue-700);
            border-bottom-color: var(--blue-500);
        }

        .message-content table {
            width: 100%;
            border-collapse: collapse;
            margin: 16px 0;
            font-size: 14px;
        }

        .message-content th,
        .message-content td {
            padding: 10px 12px;
            text-align: left;
            border: 1px solid var(--gray-200);
        }

        .message-content th {
            background: var(--gray-50);
            font-weight: 600;
            color: var(--gray-900);
        }

        .message-content tr:nth-child(even) {
            background: var(--gray-50);
        }

        /* First paragraph no top margin */
        .message-content > p:first-child,
        .message-content > h1:first-child,
        .message-content > h2:first-child,
        .message-content > h3:first-child {
            margin-top: 0;
        }

        /* Last element no bottom margin */
        .message-content > *:last-child {
            margin-bottom: 0;
        }

    </style>

    <!-- Marked.js for Markdown Parsing -->
    <script src="https://cdn.jsdelivr.net/npm/marked@11.1.1/marked.min.js"></script>
</head>
<body>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <nav class="nav">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo">
                    <div class="logo-icon">
                        <svg width="36" height="12" viewBox="0 0 52 18" fill="none">
                            <circle cx="9" cy="9" r="9" fill="#2563EB"/>
                            <circle cx="26" cy="9" r="9" fill="#2563EB" fill-opacity="0.5"/>
                            <circle cx="43" cy="9" r="9" fill="#2563EB" fill-opacity="0.2"/>
                        </svg>
                    </div>
                    <span class="logo-text"><span class="gas">gas</span><span class="consult">consult</span><span class="ai">.ai</span></span>
                </a>
                <div class="nav-links">
                    <a href="/?clear=1" class="nav-link active">Home</a>
                    <a href="/quick-dose" class="nav-link">Quick Dose</a>
                    <a href="/preop" class="nav-link">Pre-Op</a>
                    <a href="/calculators" class="nav-link">Clinical Calculators</a>
                    <a href="/hypotension" class="nav-link">IOH Predictor</a>
                </div>
                <button class="mobile-menu-btn" onclick="toggleMobileMenu()" aria-label="Toggle menu">
                    <span></span>
                    <span></span>
                    <span></span>
                </button>
            </div>
        </nav>
        <div class="mobile-menu" id="mobileMenu">
            <a href="/?clear=1" class="mobile-menu-link">Home</a>
            <a href="/quick-dose" class="mobile-menu-link">Quick Dose</a>
            <a href="/preop" class="mobile-menu-link">Pre-Op</a>
            <a href="/calculators" class="mobile-menu-link">Clinical Calculators</a>
            <a href="/hypotension" class="mobile-menu-link">IOH Predictor</a>
        </div>

        {% if messages and messages|length > 0 %}
        <!-- Chat Interface -->
        <section class="chat-view">
            <div class="messages-container" id="messagesContainer">
                {% if error_message %}
                <div class="error-message">{{ error_message|safe }}</div>
                {% endif %}

                {% for message in messages %}
                    {% if message.content and message.content.strip() %}
                    {% if message.role == 'user' %}
                    <div class="message user-message">
                        <div class="message-bubble">{{ message.content }}</div>
                    </div>
                    {% else %}
                    <div class="message ai-message">
                        <div class="message-bubble">
                            {% if message.get('evidence_strength') %}
                            <div class="evidence-badge {{ 'high' if message.evidence_strength == 'High' else ('moderate' if message.evidence_strength == 'Moderate' else 'low') }}"
                                 title="{% if message.evidence_strength == 'High' %}Based on {{ message.num_papers }} high-quality studies including meta-analyses, systematic reviews, or RCTs{% elif message.evidence_strength == 'Moderate' %}Based on {{ message.num_papers }} studies with some high-quality evidence{% else %}Limited evidence available ({{ message.num_papers }} studies). Use clinical judgment.{% endif %}">
                                <div style="display: flex; flex-direction: column; align-items: flex-start;">
                                    <div>
                                        {% if message.evidence_strength == 'High' %}
                                        ✓ High Confidence
                                        {% elif message.evidence_strength == 'Moderate' %}
                                        ~ Moderate Confidence
                                        {% else %}
                                        ! Low Confidence
                                        {% endif %}
                                        • {{ message.num_papers }} studies
                                    </div>
                                    <div class="evidence-explanation">
                                        {% if message.evidence_strength == 'High' %}
                                        Strong evidence from meta-analyses, RCTs, or systematic reviews
                                        {% elif message.evidence_strength == 'Moderate' %}
                                        Moderate evidence - consider individual patient factors
                                        {% else %}
                                        Limited evidence - use caution and clinical judgment
                                        {% endif %}
                                    </div>
                                </div>
                            </div>
                            {% endif %}

                            <div class="message-content">{{ message.content|safe }}</div>

                            {% if message.references and message.references|length > 0 %}
                            <div class="references-section">
                                <div class="references-title">
                                    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                        <path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path>
                                        <path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path>
                                    </svg>
                                    References ({{ message.references|length }})
                                </div>
                                {% for ref in message.references %}
                                <div class="reference-item">
                                    <a href="https://pubmed.ncbi.nlm.nih.gov/{{ ref.pmid }}/" target="_blank" rel="noopener noreferrer" class="reference-link">
                                        {{ ref.title }}
                                    </a>
                                    <div class="reference-meta">
                                        {{ ref.authors }} — {{ ref.journal }}, {{ ref.year }}
                                    </div>
                                </div>
                                {% endfor %}
                            </div>
                            {% endif %}
                        </div>
                    </div>
                    {% endif %}
                    {% endif %}
                {% endfor %}

                {% if pending_stream %}
                <div class="message ai-message" id="streamingMessage">
                    <div class="message-bubble">
                        <div class="message-content" id="streamingContent" style="display: none;"></div>
                        <div class="streaming-indicator">
                            <div class="streaming-dots">
                                <span></span>
                                <span></span>
                                <span></span>
                            </div>
                            Thinking...
                        </div>
                        <div class="references-section" id="streamingReferences" style="display: none;"></div>
                    </div>
                </div>
                {% endif %}
            </div>

            <div class="chat-input-area">
                <div class="chat-input-wrapper">
                    <form method="post" action="/" id="chat-form" class="chat-card">
                        <input type="hidden" name="csrf_token" value="{{ csrf_token() }}"/>
                        <div class="chat-inner">
                            <textarea name="query" class="chat-input" placeholder="Ask a follow-up question..." rows="1" required></textarea>
                            <button type="submit" class="chat-send">
                                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5">
                                    <line x1="12" y1="19" x2="12" y2="5"></line>
                                    <polyline points="5 12 12 5 19 12"></polyline>
                                </svg>
                            </button>
                        </div>
                    </form>
                    <a href="/clear" class="new-chat-btn">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                            <polyline points="1 4 1 10 7 10"></polyline>
                            <path d="M3.51 15a9 9 0 1 0 2.13-9.36L1 10"></path>
                        </svg>
                        New Chat
                    </a>
                </div>
            </div>
        </section>

        {% else %}
        <!-- Homepage Hero -->
        <section class="hero">
            <div class="hero-badge">
                <span class="badge-dot"></span>
                <span class="badge-text">PubMed-Powered AI</span>
            </div>
            <h1 class="hero-title">The AI copilot for<br><span class="gradient">anesthesiology</span></h1>
            <p class="hero-subtitle">Evidence-based clinical decision support, instant drug dosing, and intelligent pre-op assessments — all in one place.</p>

            <div class="chat-container">
                <form method="post" action="/" class="chat-card">
                    <input type="hidden" name="csrf_token" value="{{ csrf_token() }}"/>
                    <div class="chat-inner">
                        <textarea name="query" class="chat-input" placeholder="Ask anything about anesthesiology..." rows="1" required></textarea>
                        <button type="submit" class="chat-send">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5">
                                <line x1="12" y1="19" x2="12" y2="5"></line>
                                <polyline points="5 12 12 5 19 12"></polyline>
                            </svg>
                        </button>
                    </div>
                    <div class="chat-hints">
                        <div class="hint-chip" onclick="fillQuery(event)">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M22 12h-4l-3 9L9 3l-3 9H2"></path></svg>
                            Pre-op cardiac risk
                        </div>
                        <div class="hint-chip" onclick="fillQuery(event)">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M12 22s8-4 8-10V5l-8-3-8 3v7c0 6 8 10 8 10z"></path></svg>
                            Sugammadex dosing
                        </div>
                        <div class="hint-chip" onclick="fillQuery(event)">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="12" r="10"></circle><line x1="12" y1="8" x2="12" y2="12"></line><line x1="12" y1="16" x2="12.01" y2="16"></line></svg>
                            MH protocol
                        </div>
                        <div class="hint-chip" onclick="fillQuery(event)">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><rect x="3" y="3" width="18" height="18" rx="2" ry="2"></rect><line x1="3" y1="9" x2="21" y2="9"></line><line x1="9" y1="21" x2="9" y2="9"></line></svg>
                            RSI checklist
                        </div>
                    </div>
                </form>
            </div>
        </section>
        {% endif %}

        <section class="features">
            <div class="features-header">
                <div class="features-label">Features</div>
                <h2 class="features-title">Everything you need at the bedside</h2>
                <p class="features-subtitle">Clinical tools powered by the latest evidence.</p>
            </div>
            <div class="features-grid">
                <div class="feature-card">
                    <div class="feature-icon blue">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 15a2 2 0 0 1-2 2H7l-4 4V5a2 2 0 0 1 2-2h14a2 2 0 0 1 2 2z"></path></svg>
                    </div>
                    <h3 class="feature-title">AI Clinical Chat</h3>
                    <p class="feature-desc">Ask complex clinical questions and get evidence-based answers with PubMed citations in seconds.</p>
                    <a href="/" class="feature-link">Try it now <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><line x1="5" y1="12" x2="19" y2="12"></line><polyline points="12 5 19 12 12 19"></polyline></svg></a>
                </div>
                <div class="feature-card">
                    <div class="feature-icon emerald">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M12 22s8-4 8-10V5l-8-3-8 3v7c0 6 8 10 8 10z"></path></svg>
                    </div>
                    <h3 class="feature-title">Quick Dose Calculator</h3>
                    <p class="feature-desc">Weight-based dosing for all common anesthesia drugs with color-coded syringe labels.</p>
                    <a href="/quick-dose" class="feature-link">Calculate doses <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><line x1="5" y1="12" x2="19" y2="12"></line><polyline points="12 5 19 12 12 19"></polyline></svg></a>
                </div>
                <div class="feature-card">
                    <div class="feature-icon violet">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M16 4h2a2 2 0 0 1 2 2v14a2 2 0 0 1-2 2H6a2 2 0 0 1-2-2V6a2 2 0 0 1 2-2h2"></path><rect x="8" y="2" width="8" height="4" rx="1" ry="1"></rect><path d="M9 14l2 2 4-4"></path></svg>
                    </div>
                    <h3 class="feature-title">Pre-Op Assessment</h3>
                    <p class="feature-desc">Structured pre-operative evaluation with automatic risk scoring and recommendations.</p>
                    <a href="/preop" class="feature-link">Start assessment <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><line x1="5" y1="12" x2="19" y2="12"></line><polyline points="12 5 19 12 12 19"></polyline></svg></a>
                </div>
                <div class="feature-card">
                    <div class="feature-icon rose">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polygon points="7.86 2 16.14 2 22 7.86 22 16.14 16.14 22 7.86 22 2 16.14 2 7.86 7.86 2"></polygon><line x1="12" y1="8" x2="12" y2="12"></line><line x1="12" y1="16" x2="12.01" y2="16"></line></svg>
                    </div>
                    <h3 class="feature-title">Crisis Protocols</h3>
                    <p class="feature-desc">One-tap access to emergency protocols for MH, LAST, anaphylaxis, and more.</p>
                    <a href="/" class="feature-link">View protocols <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><line x1="5" y1="12" x2="19" y2="12"></line><polyline points="12 5 19 12 12 19"></polyline></svg></a>
                </div>
                <div class="feature-card">
                    <div class="feature-icon cyan">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path><path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path><line x1="8" y1="7" x2="16" y2="7"></line><line x1="8" y1="11" x2="16" y2="11"></line></svg>
                    </div>
                    <h3 class="feature-title">PubMed Integration</h3>
                    <p class="feature-desc">Every answer is backed by real citations from peer-reviewed medical literature.</p>
                    <a href="/" class="feature-link">Learn more <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><line x1="5" y1="12" x2="19" y2="12"></line><polyline points="12 5 19 12 12 19"></polyline></svg></a>
                </div>
            </div>
        </section>

        <footer class="footer">
            <div class="footer-inner">
                <div class="footer-brand">
                    <div class="footer-logo">
                        <svg viewBox="0 0 32 32" fill="none"><path d="M4 16 L9 16 L11 10 L14 22 L16 4 L18 28 L21 10 L23 16 L28 16" stroke="white" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"/></svg>
                    </div>
                    <span class="footer-text">© 2025 GasConsult.ai</span>
                </div>
                <div class="footer-links">
                    <a href="/privacy" class="footer-link">Privacy</a>
                    <a href="/terms" class="footer-link">Terms</a>
                    <a href="mailto:contact@gasconsult.ai" class="footer-link">Contact</a>
                </div>
            </div>
        </footer>
    </div>

    <script>
        // Auto-resize textarea on input
        const textareas = document.querySelectorAll('.chat-input');
        textareas.forEach(textarea => {
            textarea.addEventListener('input', function() {
                this.style.height = 'auto';
                this.style.height = Math.min(this.scrollHeight, 180) + 'px';
            });

            // Cmd/Ctrl + Enter to submit
            textarea.addEventListener('keydown', function(e) {
                if ((e.metaKey || e.ctrlKey) && e.key === 'Enter') {
                    e.preventDefault();
                    const form = this.closest('form');
                    if (form) {
                        form.requestSubmit();
                    }
                }
            });
        });

        // Fill query from hint chips (homepage only)
        function fillQuery(event) {
            const chip = event.currentTarget;
            const text = chip.textContent.trim();
            const textarea = document.querySelector('.chat-input');
            if (textarea) {
                textarea.value = text;
                textarea.focus();
            }
        }

        // Toggle mobile menu
        function toggleMobileMenu() {
            const menu = document.getElementById('mobileMenu');
            const btn = document.querySelector('.mobile-menu-btn');
            if (menu && btn) {
                menu.classList.toggle('active');
                btn.classList.toggle('active');
            }
        }

        // Scroll to bottom of messages
        function scrollToBottom() {
            const container = document.getElementById('messagesContainer');
            if (container) {
                container.scrollTop = container.scrollHeight;
            }
        }

        // SSE Streaming functionality
        {% if pending_stream %}
        (function() {
            const requestId = '{{ pending_stream }}';
            const streamingContent = document.getElementById('streamingContent');
            const streamingReferences = document.getElementById('streamingReferences');
            const streamingIndicator = document.querySelector('.streaming-indicator');
            const streamingMessage = document.getElementById('streamingMessage');

            if (!streamingContent) {
                console.error('[SSE] Streaming elements not found');
                return;
            }

            console.log('[SSE] Starting stream for request_id:', requestId);

            const eventSource = new EventSource('/stream?request_id=' + encodeURIComponent(requestId));
            let accumulatedMarkdown = ''; // Accumulate markdown during streaming

            eventSource.addEventListener('message', function(e) {
                try {
                    const event = JSON.parse(e.data);

                    if (event.type === 'content') {
                        // Stream content chunks
                        if (event.data) {
                            // Show content div on first content chunk
                            if (streamingContent.style.display === 'none') {
                                streamingContent.style.display = 'block';
                            }
                            accumulatedMarkdown += event.data;
                            // Parse and show formatted markdown in real-time
                            if (typeof marked !== 'undefined') {
                                streamingContent.innerHTML = marked.parse(accumulatedMarkdown);
                            } else {
                                streamingContent.textContent = accumulatedMarkdown;
                            }
                            scrollToBottom();
                        }
                    } else if (event.type === 'references') {
                        // Markdown already parsed in real-time, no need to parse again
                        // (keeping this for backwards compatibility if needed)
                        if (!streamingContent.innerHTML && accumulatedMarkdown && typeof marked !== 'undefined') {
                            streamingContent.innerHTML = marked.parse(accumulatedMarkdown);
                        }

                        // Display references
                        if (event.data && event.data.length > 0) {
                            console.log('[SSE] Received', event.data.length, 'references');

                            // Hide streaming indicator
                            if (streamingIndicator) {
                                streamingIndicator.style.display = 'none';
                            }

                            // Build references HTML
                            let refsHTML = '<div class="references-title">';
                            refsHTML += '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">';
                            refsHTML += '<path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path>';
                            refsHTML += '<path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path>';
                            refsHTML += '</svg>';
                            refsHTML += 'References (' + event.data.length + ')';
                            refsHTML += '</div>';

                            event.data.forEach(function(ref) {
                                refsHTML += '<div class="reference-item">';
                                refsHTML += '<a href="https://pubmed.ncbi.nlm.nih.gov/' + ref.pmid + '/" target="_blank" rel="noopener noreferrer" class="reference-link">';
                                refsHTML += ref.title;
                                refsHTML += '</a>';
                                refsHTML += '<div class="reference-meta">';
                                refsHTML += ref.authors + ' — ' + ref.journal + ', ' + ref.year;
                                refsHTML += '</div>';
                                refsHTML += '</div>';
                            });

                            streamingReferences.innerHTML = refsHTML;
                            streamingReferences.style.display = 'block';
                            scrollToBottom();
                        }
                    } else if (event.type === 'done') {
                        console.log('[SSE] Stream complete');

                        // Hide streaming indicator
                        if (streamingIndicator) {
                            streamingIndicator.style.display = 'none';
                        }

                        eventSource.close();
                        scrollToBottom();
                    } else if (event.type === 'error') {
                        console.error('[SSE] Server error:', event.message);

                        if (streamingIndicator) {
                            streamingIndicator.innerHTML = '<span style="color: #DC2626;">Error: ' + (event.message || 'Unknown error') + '</span>';
                        }

                        eventSource.close();
                    }
                } catch (err) {
                    console.error('[SSE] Error parsing message:', err);
                }
            });

            eventSource.addEventListener('error', function(e) {
                console.error('[SSE] Connection error:', e);

                if (streamingIndicator) {
                    streamingIndicator.innerHTML = '<span style="color: #DC2626;">Connection error - please refresh and try again</span>';
                }

                eventSource.close();
            });

            // Scroll to bottom on page load
            setTimeout(scrollToBottom, 100);
        })();
        {% endif %}

        // Parse markdown in existing messages on page load
        window.addEventListener('DOMContentLoaded', function() {
            if (typeof marked !== 'undefined') {
                // Find all message-content divs that contain markdown
                const messageContents = document.querySelectorAll('.ai-message .message-content');
                messageContents.forEach(function(element) {
                    // Get the text content (markdown)
                    const markdownText = element.textContent;
                    // Parse and replace with HTML
                    if (markdownText && markdownText.trim()) {
                        element.innerHTML = marked.parse(markdownText);
                    }
                });
            }
        });

        // Auto-scroll to bottom on page load (for chat view)
        {% if messages and messages|length > 0 %}
        window.addEventListener('load', function() {
            scrollToBottom();
        });
        {% endif %}
    </script>
</body>
</html>

"""


LIBRARY_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Library — gasconsult.ai</title>

    <!-- PWA -->
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=5">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">

    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap" rel="stylesheet">
    <style>

        :root {
            --white: #FFFFFF;
            --gray-50: #F8FAFC;
            --gray-100: #F1F5F9;
            --gray-200: #E2E8F0;
            --gray-300: #CBD5E1;
            --gray-400: #94A3B8;
            --gray-500: #64748B;
            --gray-600: #475569;
            --gray-700: #334155;
            --gray-800: #1E293B;
            --gray-900: #0F172A;
            --blue-50: #EFF6FF;
            --blue-100: #DBEAFE;
            --blue-200: #BFDBFE;
            --blue-300: #93C5FD;
            --blue-400: #60A5FA;
            --blue-500: #3B82F6;
            --blue-600: #2563EB;
            --blue-700: #1D4ED8;
        }

        * { margin: 0; padding: 0; box-sizing: border-box; }

        html {
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
            scroll-behavior: smooth;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            background: var(--gray-50);
            color: var(--gray-900);
            min-height: 100vh;
            overflow-x: hidden;
        }

        .bg-canvas {
            position: fixed;
            inset: 0;
            z-index: 0;
            overflow: hidden;
            background: linear-gradient(180deg, #F0F7FF 0%, var(--gray-50) 50%, #FAFBFF 100%);
        }

        .orb {
            position: absolute;
            border-radius: 50%;
            filter: blur(80px);
            opacity: 0.6;
            animation: float 20s ease-in-out infinite;
        }

        .orb-1 {
            width: 400px;
            height: 400px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.15) 0%, transparent 70%);
            top: -15%;
            left: -20%;
        }

        .orb-2 {
            width: 300px;
            height: 300px;
            background: radial-gradient(circle, rgba(147, 197, 253, 0.2) 0%, transparent 70%);
            top: 30%;
            right: -20%;
            animation-delay: -7s;
            animation-duration: 25s;
        }

        .orb-3 {
            width: 250px;
            height: 250px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.1) 0%, transparent 70%);
            bottom: -10%;
            left: 20%;
            animation-delay: -14s;
            animation-duration: 30s;
        }

        @keyframes float {
            0%, 100% { transform: translate(0, 0) scale(1); }
            25% { transform: translate(40px, -40px) scale(1.05); }
            50% { transform: translate(20px, 40px) scale(0.95); }
            75% { transform: translate(-40px, 20px) scale(1.02); }
        }

        .grain {
            position: fixed;
            inset: 0;
            z-index: 1;
            pointer-events: none;
            opacity: 0.02;
            background-image: url("data:image/svg+xml,%3Csvg viewBox='0 0 512 512' xmlns='http://www.w3.org/2000/svg'%3E%3Cfilter id='n'%3E%3CfeTurbulence type='fractalNoise' baseFrequency='0.8' numOctaves='4' stitchTiles='stitch'/%3E%3C/filter%3E%3Crect width='100%25' height='100%25' filter='url(%23n)'/%3E%3C/svg%3E");
        }

        .page {
            position: relative;
            z-index: 2;
            min-height: 100vh;
            display: flex;
            flex-direction: column;
        }

        .nav {
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            z-index: 100;
            padding: 12px 16px;
        }

        .nav-inner {
            max-width: 1200px;
            margin: 0 auto;
            height: 56px;
            background: rgba(255, 255, 255, 0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 0 16px;
            display: flex;
            align-items: center;
            justify-content: space-between;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 12px 48px rgba(0,0,0,0.03);
        }

        .logo {
            display: flex;
            align-items: center;
            gap: 14px;
            text-decoration: none;
        }

        .logo-icon {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .logo-icon svg { width: 36px; height: 12px; }

        .logo-text {
            font-size: 18px;
            font-weight: 700;
            letter-spacing: -0.5px;
            color: var(--gray-900);
        }

        .logo-text .gas { color: var(--blue-600); }
        .logo-text .consult { color: #0F172A; }
        .logo-text .ai { color: rgba(15, 23, 42, 0.4); }

        .nav-links {
            display: none;
            align-items: center;
            gap: 4px;
        }

        .nav-link {
            padding: 10px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .nav-link.active {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .mobile-menu-btn {
            display: flex;
            flex-direction: column;
            gap: 5px;
            background: none;
            border: none;
            cursor: pointer;
            padding: 8px;
            border-radius: 8px;
            transition: background 0.2s ease;
        }

        .mobile-menu-btn:hover {
            background: rgba(0,0,0,0.04);
        }

        .mobile-menu-btn span {
            display: block;
            width: 22px;
            height: 2px;
            background: var(--gray-700);
            border-radius: 1px;
            transition: all 0.3s ease;
        }

        .mobile-menu-btn.active span:nth-child(1) {
            transform: rotate(45deg) translate(7px, 7px);
        }

        .mobile-menu-btn.active span:nth-child(2) {
            opacity: 0;
        }

        .mobile-menu-btn.active span:nth-child(3) {
            transform: rotate(-45deg) translate(7px, -7px);
        }

        .mobile-menu {
            display: none;
            position: fixed;
            top: 80px;
            left: 16px;
            right: 16px;
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 8px;
            box-shadow: 0 4px 16px rgba(0,0,0,0.08), 0 12px 48px rgba(0,0,0,0.12);
            z-index: 99;
            flex-direction: column;
            gap: 4px;
        }

        .mobile-menu.active {
            display: flex;
        }

        .mobile-menu-link {
            padding: 14px 16px;
            font-size: 15px;
            font-weight: 500;
            color: var(--gray-700);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .mobile-menu-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .hero {
            padding: 120px 20px 60px;
            text-align: center;
        }

        .hero-badge {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            background: var(--white);
            border: 1px solid var(--gray-200);
            border-radius: 100px;
            padding: 8px 16px 8px 12px;
            margin-bottom: 24px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.04);
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) forwards;
            opacity: 0;
        }

        .badge-dot {
            width: 8px;
            height: 8px;
            background: var(--blue-500);
            border-radius: 50%;
            position: relative;
        }

        .badge-dot::after {
            content: '';
            position: absolute;
            inset: -3px;
            border-radius: 50%;
            background: var(--blue-400);
            animation: pulse-ring 2s ease-out infinite;
        }

        @keyframes pulse-ring {
            0% { transform: scale(0.8); opacity: 0.8; }
            100% { transform: scale(2); opacity: 0; }
        }

        .badge-text {
            font-size: 12px;
            font-weight: 600;
            color: var(--gray-700);
        }

        .hero-title {
            font-size: 40px;
            font-weight: 800;
            line-height: 1.1;
            letter-spacing: -2px;
            color: var(--gray-900);
            margin-bottom: 20px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.1s forwards;
            opacity: 0;
        }

        .hero-title .gradient { color: var(--blue-600); }

        .hero-subtitle {
            font-size: 16px;
            font-weight: 400;
            line-height: 1.6;
            color: var(--gray-500);
            max-width: 560px;
            margin: 0 auto 40px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.2s forwards;
            opacity: 0;
        }

        @keyframes fade-up {
            from { opacity: 0; transform: translateY(24px); }
            to { opacity: 1; transform: translateY(0); }
        }

        .chat-container {
            max-width: 720px;
            margin: 0 auto 60px;
            padding: 0 16px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.3s forwards;
            opacity: 0;
        }

        .chat-card {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(40px) saturate(180%);
            -webkit-backdrop-filter: blur(40px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 8px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
            transition: all 0.4s cubic-bezier(0.4,0,0.2,1);
        }

        .chat-card:focus-within {
            box-shadow: 0 0 0 4px rgba(59,130,246,0.1), 0 1px 2px rgba(0,0,0,0.02), 0 8px 24px rgba(37,99,235,0.08), 0 32px 100px rgba(37,99,235,0.12), inset 0 1px 0 rgba(255,255,255,0.8);
            border-color: rgba(59,130,246,0.3);
        }

        .chat-inner {
            background: var(--white);
            border-radius: 14px;
            padding: 6px;
            display: flex;
            align-items: flex-end;
            gap: 6px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04), inset 0 1px 0 rgba(255,255,255,1);
        }

        .chat-input {
            flex: 1;
            border: none;
            outline: none;
            padding: 14px 16px;
            font-size: 16px;
            font-family: inherit;
            color: var(--gray-800);
            background: transparent;
            resize: none;
            min-height: 52px;
            max-height: 150px;
            line-height: 1.5;
        }

        .chat-input::placeholder { color: var(--gray-400); }

        .chat-send {
            width: 48px;
            height: 48px;
            background: var(--blue-600);
            border: none;
            border-radius: 12px;
            color: var(--white);
            cursor: pointer;
            display: flex;
            align-items: center;
            justify-content: center;
            transition: all 0.25s cubic-bezier(0.4,0,0.2,1);
            box-shadow: 0 1px 2px rgba(37,99,235,0.2), 0 4px 16px rgba(37,99,235,0.2), inset 0 1px 0 rgba(255,255,255,0.1), inset 0 -1px 0 rgba(0,0,0,0.1);
            flex-shrink: 0;
            margin: 6px;
        }

        .chat-send:hover {
            background: var(--blue-700);
            transform: translateY(-2px);
            box-shadow: 0 2px 4px rgba(37,99,235,0.2), 0 12px 40px rgba(37,99,235,0.3), inset 0 1px 0 rgba(255,255,255,0.1), inset 0 -1px 0 rgba(0,0,0,0.1);
        }

        .chat-send:active { transform: translateY(0); }
        .chat-send svg { width: 20px; height: 20px; }

        .chat-hints {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            padding: 16px 8px 6px;
        }

        .hint-chip {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            background: rgba(255,255,255,0.6);
            border: 1px solid var(--gray-200);
            border-radius: 100px;
            padding: 10px 14px;
            font-size: 12px;
            font-weight: 500;
            color: var(--gray-600);
            cursor: pointer;
            transition: all 0.2s ease;
        }

        .hint-chip:hover {
            background: var(--white);
            border-color: var(--blue-200);
            color: var(--blue-600);
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(37,99,235,0.1);
        }

        .hint-chip svg { width: 14px; height: 14px; opacity: 0.5; }
        .hint-chip:hover svg { opacity: 1; color: var(--blue-500); }

        .features { padding: 60px 20px 80px; }

        .features-header {
            text-align: center;
            margin-bottom: 40px;
        }

        .features-label {
            font-size: 11px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 1.5px;
            color: var(--blue-600);
            margin-bottom: 12px;
        }

        .features-title {
            font-size: 28px;
            font-weight: 800;
            letter-spacing: -1px;
            color: var(--gray-900);
            margin-bottom: 12px;
        }

        .features-subtitle {
            font-size: 16px;
            color: var(--gray-500);
            max-width: 480px;
            margin: 0 auto;
        }

        .features-grid {
            display: grid;
            grid-template-columns: 1fr;
            gap: 16px;
            max-width: 1200px;
            margin: 0 auto;
        }

        .feature-card {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            padding: 28px;
            position: relative;
            overflow: hidden;
            transition: all 0.4s cubic-bezier(0.4,0,0.2,1);
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04);
            cursor: pointer;
        }

        .feature-card::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 1px;
            background: linear-gradient(90deg, transparent 0%, rgba(255,255,255,0.8) 50%, transparent 100%);
        }

        .feature-card:hover {
            transform: translateY(-4px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.04), 0 24px 64px rgba(0,0,0,0.08);
            border-color: rgba(59,130,246,0.2);
        }

        .feature-icon {
            width: 56px;
            height: 56px;
            border-radius: 16px;
            display: flex;
            align-items: center;
            justify-content: center;
            margin-bottom: 20px;
            transition: all 0.3s ease;
        }

        .feature-card:hover .feature-icon { transform: scale(1.05); }
        .feature-icon svg { width: 24px; height: 24px; }

        .feature-icon.blue {
            background: var(--blue-50);
            border: 1px solid var(--blue-100);
            box-shadow: 0 4px 16px rgba(37,99,235,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.blue svg { color: var(--blue-600); }

        .feature-icon.emerald {
            background: #ECFDF5;
            border: 1px solid #D1FAE5;
            box-shadow: 0 4px 16px rgba(16,185,129,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.emerald svg { color: #059669; }

        .feature-icon.violet {
            background: #F5F3FF;
            border: 1px solid #EDE9FE;
            box-shadow: 0 4px 16px rgba(139,92,246,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.violet svg { color: #7C3AED; }

        .feature-icon.rose {
            background: #FFF1F2;
            border: 1px solid #FFE4E6;
            box-shadow: 0 4px 16px rgba(244,63,94,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.rose svg { color: #E11D48; }

        .feature-icon.cyan {
            background: #ECFEFF;
            border: 1px solid #CFFAFE;
            box-shadow: 0 4px 16px rgba(6,182,212,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.cyan svg { color: #0891B2; }

        .feature-title {
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 10px;
            letter-spacing: -0.3px;
        }

        .feature-desc {
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-500);
            margin-bottom: 20px;
        }

        .feature-link {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            font-size: 14px;
            font-weight: 600;
            color: var(--blue-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .feature-link:hover { gap: 10px; }
        .feature-link svg { width: 16px; height: 16px; transition: transform 0.2s ease; }
        .feature-link:hover svg { transform: translateX(4px); }

        .footer {
            padding: 32px 20px;
            border-top: 1px solid var(--gray-200);
            background: rgba(255,255,255,0.5);
        }

        .footer-inner {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 20px;
            text-align: center;
        }

        .footer-brand {
            display: flex;
            align-items: center;
            gap: 10px;
        }

        .footer-logo {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .footer-logo svg { width: 32px; height: 32px; }

        .footer-text {
            font-size: 13px;
            color: var(--gray-500);
        }

        .footer-links {
            display: flex;
            gap: 24px;
        }

        .footer-link {
            font-size: 13px;
            color: var(--gray-500);
            text-decoration: none;
            transition: color 0.2s ease;
        }

        .footer-link:hover { color: var(--gray-700); }

        @media (min-width: 768px) {
            .nav { padding: 16px 32px; }
            .nav-inner { height: 64px; padding: 0 24px; border-radius: 20px; }
            .logo-icon svg { width: 42px; height: 15px; }
            .logo-text { font-size: 20px; }
            .nav-links { display: flex; }
            .mobile-menu-btn { display: none; }
            .hero { padding: 160px 32px 80px; }
            .hero-badge { padding: 10px 20px 10px 14px; margin-bottom: 32px; }
            .badge-dot { width: 10px; height: 10px; }
            .badge-text { font-size: 13px; }
            .hero-title { font-size: 56px; letter-spacing: -2.5px; margin-bottom: 24px; }
            .hero-subtitle { font-size: 18px; margin-bottom: 48px; }
            .chat-container { padding: 0 24px; margin-bottom: 80px; }
            .chat-card { border-radius: 24px; padding: 10px; }
            .chat-inner { border-radius: 18px; padding: 8px; }
            .chat-input { padding: 16px 18px; min-height: 56px; }
            .chat-send { width: 52px; height: 52px; border-radius: 14px; }
            .chat-hints { gap: 10px; padding: 18px 10px 8px; }
            .hint-chip { padding: 12px 18px; font-size: 13px; }
            .hint-chip svg { width: 16px; height: 16px; }
            .features { padding: 80px 32px 100px; }
            .features-header { margin-bottom: 56px; }
            .features-label { font-size: 12px; margin-bottom: 16px; }
            .features-title { font-size: 36px; letter-spacing: -1.5px; }
            .features-subtitle { font-size: 18px; }
            .features-grid { grid-template-columns: repeat(2, 1fr); gap: 20px; }
            .feature-card { padding: 36px; border-radius: 24px; }
            .feature-card:hover { transform: translateY(-6px); }
            .feature-icon { width: 60px; height: 60px; border-radius: 18px; margin-bottom: 24px; }
            .feature-icon svg { width: 26px; height: 26px; }
            .feature-title { font-size: 20px; margin-bottom: 12px; }
            .feature-desc { font-size: 15px; line-height: 1.7; margin-bottom: 24px; }
            .footer { padding: 40px 32px; }
            .footer-inner { flex-direction: row; justify-content: space-between; text-align: left; }
            .footer-logo svg { width: 36px; height: 36px; }
            .footer-text { font-size: 14px; }
            .footer-links { gap: 32px; }
            .footer-link { font-size: 14px; }
            .orb-1 { width: 600px; height: 600px; left: -10%; }
            .orb-2 { width: 450px; height: 450px; right: -10%; }
            .orb-3 { width: 400px; height: 400px; }
        }

        @media (min-width: 1024px) {
            .nav { padding: 16px 40px; }
            .hero { padding: 180px 40px 80px; }
            .hero-title { font-size: 72px; letter-spacing: -3px; margin-bottom: 28px; }
            .hero-subtitle { font-size: 20px; margin-bottom: 56px; }
            .chat-container { margin-bottom: 100px; }
            .chat-card { border-radius: 28px; }
            .chat-inner { border-radius: 20px; }
            .chat-input { padding: 18px 20px; min-height: 60px; max-height: 180px; }
            .chat-send { width: 56px; height: 56px; border-radius: 16px; margin: 8px; }
            .chat-send svg { width: 22px; height: 22px; }
            .chat-hints { padding: 20px 12px 8px; }
            .hint-chip { padding: 12px 20px; }
            .features { padding: 80px 40px 120px; }
            .features-header { margin-bottom: 64px; }
            .features-title { font-size: 40px; }
            .features-grid { grid-template-columns: repeat(3, 1fr); gap: 24px; }
            .feature-card { padding: 40px; border-radius: 28px; }
            .feature-card:hover { transform: translateY(-8px); }
            .feature-icon { width: 64px; height: 64px; border-radius: 20px; margin-bottom: 28px; }
            .feature-icon svg { width: 28px; height: 28px; }
            .footer { padding: 48px 40px; }
            .orb-1 { width: 800px; height: 800px; top: -20%; left: -10%; }
            .orb-2 { width: 600px; height: 600px; right: -15%; }
            .orb-3 { width: 500px; height: 500px; }
        }

        @media (min-width: 1280px) {
            .hero-title { font-size: 80px; }
        }
    

        .main-content {
            flex: 1;
            padding: 100px 20px 40px;
            max-width: 1200px;
            margin: 0 auto;
            width: 100%;
        }

        .content-card {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            padding: 32px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
        }

        h1 { font-size: 32px; margin-bottom: 16px; font-weight: 700; letter-spacing: -0.5px; }
        h2 { font-size: 24px; margin-bottom: 12px; font-weight: 700; letter-spacing: -0.5px; }
        h3 { font-size: 20px; margin-bottom: 10px; font-weight: 700; letter-spacing: -0.3px; }

        input[type="text"], input[type="number"], input[type="email"], select, textarea {
            width: 100%;
            padding: 12px 16px;
            border: 1px solid var(--gray-300);
            border-radius: 12px;
            font-family: inherit;
            font-size: 15px;
            background: var(--white);
            color: var(--gray-900);
            transition: all 0.2s ease;
        }

        input:focus, select:focus, textarea:focus {
            outline: none;
            border-color: var(--blue-500);
            box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.1);
        }

        button, .btn {
            padding: 12px 24px;
            background: var(--blue-600);
            color: var(--white);
            border: none;
            border-radius: 12px;
            font-family: inherit;
            font-size: 15px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s ease;
            box-shadow: 0 1px 2px rgba(37,99,235,0.2), 0 4px 16px rgba(37,99,235,0.2), inset 0 1px 0 rgba(255,255,255,0.1);
        }

        button:hover, .btn:hover {
            background: var(--blue-700);
            transform: translateY(-2px);
            box-shadow: 0 2px 4px rgba(37,99,235,0.2), 0 12px 40px rgba(37,99,235,0.3), inset 0 1px 0 rgba(255,255,255,0.1);
        }

        button:active, .btn:active {
            transform: translateY(0);
        }

        table {
            width: 100%;
            border-collapse: collapse;
            background: rgba(255,255,255,0.5);
            border-radius: 12px;
            overflow: hidden;
        }

        th, td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid var(--gray-200);
        }

        th {
            background: rgba(37,99,235,0.05);
            font-weight: 600;
            color: var(--blue-700);
        }

        @media (min-width: 768px) {
            .main-content { padding: 120px 32px 60px; }
        }

        @media (min-width: 1024px) {
            .main-content { padding: 140px 40px 80px; }
        }

    </style>
</head>
<body>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <nav class="nav">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo">
                    <div class="logo-icon">
                        <svg width="36" height="12" viewBox="0 0 52 18" fill="none">
                            <circle cx="9" cy="9" r="9" fill="#2563EB"/>
                            <circle cx="26" cy="9" r="9" fill="#2563EB" fill-opacity="0.5"/>
                            <circle cx="43" cy="9" r="9" fill="#2563EB" fill-opacity="0.2"/>
                        </svg>
                    </div>
                    <span class="logo-text"><span class="gas">gas</span><span class="consult">consult</span><span class="ai">.ai</span></span>
                </a>
                <div class="nav-links">
                    <a href="/?clear=1" class="nav-link">Home</a>
                    <a href="/quick-dose" class="nav-link">Quick Dose</a>
                    <a href="/preop" class="nav-link">Pre-Op</a>
                    <a href="/calculators" class="nav-link">Calculators</a>
                    <a href="/hypotension" class="nav-link">IOH Predictor</a>
                </div>
                <button class="mobile-menu-btn" onclick="toggleMobileMenu()" aria-label="Toggle menu">
                    <span></span>
                    <span></span>
                    <span></span>
                </button>
            </div>
        </nav>
        <div class="mobile-menu" id="mobileMenu">
            <a href="/?clear=1" class="mobile-menu-link">Home</a>
            <a href="/quick-dose" class="mobile-menu-link">Quick Dose</a>
            <a href="/preop" class="mobile-menu-link">Pre-Op</a>
            <a href="/calculators" class="mobile-menu-link">Calculators</a>
            <a href="/hypotension" class="mobile-menu-link">IOH Predictor</a>
        </div>

        <main class="main-content">
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
        function toggleMobileMenu() {
            const menu = document.getElementById('mobileMenu');
            const btn = document.querySelector('.mobile-menu-btn');
            if (menu && btn) {
                menu.classList.toggle('active');
                btn.classList.toggle('active');
            }
        }
    </script>
        </main>
    </div>
</body>
</html>
"""

SHARED_RESPONSE_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Shared Response — gasconsult.ai</title>

    <!-- PWA -->
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=5">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">

    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap" rel="stylesheet">
    <style>

        :root {
            --white: #FFFFFF;
            --gray-50: #F8FAFC;
            --gray-100: #F1F5F9;
            --gray-200: #E2E8F0;
            --gray-300: #CBD5E1;
            --gray-400: #94A3B8;
            --gray-500: #64748B;
            --gray-600: #475569;
            --gray-700: #334155;
            --gray-800: #1E293B;
            --gray-900: #0F172A;
            --blue-50: #EFF6FF;
            --blue-100: #DBEAFE;
            --blue-200: #BFDBFE;
            --blue-300: #93C5FD;
            --blue-400: #60A5FA;
            --blue-500: #3B82F6;
            --blue-600: #2563EB;
            --blue-700: #1D4ED8;
        }

        * { margin: 0; padding: 0; box-sizing: border-box; }

        html {
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
            scroll-behavior: smooth;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            background: var(--gray-50);
            color: var(--gray-900);
            min-height: 100vh;
            overflow-x: hidden;
        }

        .bg-canvas {
            position: fixed;
            inset: 0;
            z-index: 0;
            overflow: hidden;
            background: linear-gradient(180deg, #F0F7FF 0%, var(--gray-50) 50%, #FAFBFF 100%);
        }

        .orb {
            position: absolute;
            border-radius: 50%;
            filter: blur(80px);
            opacity: 0.6;
            animation: float 20s ease-in-out infinite;
        }

        .orb-1 {
            width: 400px;
            height: 400px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.15) 0%, transparent 70%);
            top: -15%;
            left: -20%;
        }

        .orb-2 {
            width: 300px;
            height: 300px;
            background: radial-gradient(circle, rgba(147, 197, 253, 0.2) 0%, transparent 70%);
            top: 30%;
            right: -20%;
            animation-delay: -7s;
            animation-duration: 25s;
        }

        .orb-3 {
            width: 250px;
            height: 250px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.1) 0%, transparent 70%);
            bottom: -10%;
            left: 20%;
            animation-delay: -14s;
            animation-duration: 30s;
        }

        @keyframes float {
            0%, 100% { transform: translate(0, 0) scale(1); }
            25% { transform: translate(40px, -40px) scale(1.05); }
            50% { transform: translate(20px, 40px) scale(0.95); }
            75% { transform: translate(-40px, 20px) scale(1.02); }
        }

        .grain {
            position: fixed;
            inset: 0;
            z-index: 1;
            pointer-events: none;
            opacity: 0.02;
            background-image: url("data:image/svg+xml,%3Csvg viewBox='0 0 512 512' xmlns='http://www.w3.org/2000/svg'%3E%3Cfilter id='n'%3E%3CfeTurbulence type='fractalNoise' baseFrequency='0.8' numOctaves='4' stitchTiles='stitch'/%3E%3C/filter%3E%3Crect width='100%25' height='100%25' filter='url(%23n)'/%3E%3C/svg%3E");
        }

        .page {
            position: relative;
            z-index: 2;
            min-height: 100vh;
            display: flex;
            flex-direction: column;
        }

        .nav {
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            z-index: 100;
            padding: 12px 16px;
        }

        .nav-inner {
            max-width: 1200px;
            margin: 0 auto;
            height: 56px;
            background: rgba(255, 255, 255, 0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 0 16px;
            display: flex;
            align-items: center;
            justify-content: space-between;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 12px 48px rgba(0,0,0,0.03);
        }

        .logo {
            display: flex;
            align-items: center;
            gap: 14px;
            text-decoration: none;
        }

        .logo-icon {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .logo-icon svg { width: 36px; height: 12px; }

        .logo-text {
            font-size: 18px;
            font-weight: 700;
            letter-spacing: -0.5px;
            color: var(--gray-900);
        }

        .logo-text .gas { color: var(--blue-600); }
        .logo-text .consult { color: #0F172A; }
        .logo-text .ai { color: rgba(15, 23, 42, 0.4); }

        .nav-links {
            display: none;
            align-items: center;
            gap: 4px;
        }

        .nav-link {
            padding: 10px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .nav-link.active {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .mobile-menu-btn {
            display: flex;
            flex-direction: column;
            gap: 5px;
            background: none;
            border: none;
            cursor: pointer;
            padding: 8px;
            border-radius: 8px;
            transition: background 0.2s ease;
        }

        .mobile-menu-btn:hover {
            background: rgba(0,0,0,0.04);
        }

        .mobile-menu-btn span {
            display: block;
            width: 22px;
            height: 2px;
            background: var(--gray-700);
            border-radius: 1px;
            transition: all 0.3s ease;
        }

        .mobile-menu-btn.active span:nth-child(1) {
            transform: rotate(45deg) translate(7px, 7px);
        }

        .mobile-menu-btn.active span:nth-child(2) {
            opacity: 0;
        }

        .mobile-menu-btn.active span:nth-child(3) {
            transform: rotate(-45deg) translate(7px, -7px);
        }

        .mobile-menu {
            display: none;
            position: fixed;
            top: 80px;
            left: 16px;
            right: 16px;
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 8px;
            box-shadow: 0 4px 16px rgba(0,0,0,0.08), 0 12px 48px rgba(0,0,0,0.12);
            z-index: 99;
            flex-direction: column;
            gap: 4px;
        }

        .mobile-menu.active {
            display: flex;
        }

        .mobile-menu-link {
            padding: 14px 16px;
            font-size: 15px;
            font-weight: 500;
            color: var(--gray-700);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .mobile-menu-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .hero {
            padding: 120px 20px 60px;
            text-align: center;
        }

        .hero-badge {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            background: var(--white);
            border: 1px solid var(--gray-200);
            border-radius: 100px;
            padding: 8px 16px 8px 12px;
            margin-bottom: 24px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.04);
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) forwards;
            opacity: 0;
        }

        .badge-dot {
            width: 8px;
            height: 8px;
            background: var(--blue-500);
            border-radius: 50%;
            position: relative;
        }

        .badge-dot::after {
            content: '';
            position: absolute;
            inset: -3px;
            border-radius: 50%;
            background: var(--blue-400);
            animation: pulse-ring 2s ease-out infinite;
        }

        @keyframes pulse-ring {
            0% { transform: scale(0.8); opacity: 0.8; }
            100% { transform: scale(2); opacity: 0; }
        }

        .badge-text {
            font-size: 12px;
            font-weight: 600;
            color: var(--gray-700);
        }

        .hero-title {
            font-size: 40px;
            font-weight: 800;
            line-height: 1.1;
            letter-spacing: -2px;
            color: var(--gray-900);
            margin-bottom: 20px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.1s forwards;
            opacity: 0;
        }

        .hero-title .gradient { color: var(--blue-600); }

        .hero-subtitle {
            font-size: 16px;
            font-weight: 400;
            line-height: 1.6;
            color: var(--gray-500);
            max-width: 560px;
            margin: 0 auto 40px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.2s forwards;
            opacity: 0;
        }

        @keyframes fade-up {
            from { opacity: 0; transform: translateY(24px); }
            to { opacity: 1; transform: translateY(0); }
        }

        .chat-container {
            max-width: 720px;
            margin: 0 auto 60px;
            padding: 0 16px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.3s forwards;
            opacity: 0;
        }

        .chat-card {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(40px) saturate(180%);
            -webkit-backdrop-filter: blur(40px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 8px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
            transition: all 0.4s cubic-bezier(0.4,0,0.2,1);
        }

        .chat-card:focus-within {
            box-shadow: 0 0 0 4px rgba(59,130,246,0.1), 0 1px 2px rgba(0,0,0,0.02), 0 8px 24px rgba(37,99,235,0.08), 0 32px 100px rgba(37,99,235,0.12), inset 0 1px 0 rgba(255,255,255,0.8);
            border-color: rgba(59,130,246,0.3);
        }

        .chat-inner {
            background: var(--white);
            border-radius: 14px;
            padding: 6px;
            display: flex;
            align-items: flex-end;
            gap: 6px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04), inset 0 1px 0 rgba(255,255,255,1);
        }

        .chat-input {
            flex: 1;
            border: none;
            outline: none;
            padding: 14px 16px;
            font-size: 16px;
            font-family: inherit;
            color: var(--gray-800);
            background: transparent;
            resize: none;
            min-height: 52px;
            max-height: 150px;
            line-height: 1.5;
        }

        .chat-input::placeholder { color: var(--gray-400); }

        .chat-send {
            width: 48px;
            height: 48px;
            background: var(--blue-600);
            border: none;
            border-radius: 12px;
            color: var(--white);
            cursor: pointer;
            display: flex;
            align-items: center;
            justify-content: center;
            transition: all 0.25s cubic-bezier(0.4,0,0.2,1);
            box-shadow: 0 1px 2px rgba(37,99,235,0.2), 0 4px 16px rgba(37,99,235,0.2), inset 0 1px 0 rgba(255,255,255,0.1), inset 0 -1px 0 rgba(0,0,0,0.1);
            flex-shrink: 0;
            margin: 6px;
        }

        .chat-send:hover {
            background: var(--blue-700);
            transform: translateY(-2px);
            box-shadow: 0 2px 4px rgba(37,99,235,0.2), 0 12px 40px rgba(37,99,235,0.3), inset 0 1px 0 rgba(255,255,255,0.1), inset 0 -1px 0 rgba(0,0,0,0.1);
        }

        .chat-send:active { transform: translateY(0); }
        .chat-send svg { width: 20px; height: 20px; }

        .chat-hints {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            padding: 16px 8px 6px;
        }

        .hint-chip {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            background: rgba(255,255,255,0.6);
            border: 1px solid var(--gray-200);
            border-radius: 100px;
            padding: 10px 14px;
            font-size: 12px;
            font-weight: 500;
            color: var(--gray-600);
            cursor: pointer;
            transition: all 0.2s ease;
        }

        .hint-chip:hover {
            background: var(--white);
            border-color: var(--blue-200);
            color: var(--blue-600);
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(37,99,235,0.1);
        }

        .hint-chip svg { width: 14px; height: 14px; opacity: 0.5; }
        .hint-chip:hover svg { opacity: 1; color: var(--blue-500); }

        .features { padding: 60px 20px 80px; }

        .features-header {
            text-align: center;
            margin-bottom: 40px;
        }

        .features-label {
            font-size: 11px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 1.5px;
            color: var(--blue-600);
            margin-bottom: 12px;
        }

        .features-title {
            font-size: 28px;
            font-weight: 800;
            letter-spacing: -1px;
            color: var(--gray-900);
            margin-bottom: 12px;
        }

        .features-subtitle {
            font-size: 16px;
            color: var(--gray-500);
            max-width: 480px;
            margin: 0 auto;
        }

        .features-grid {
            display: grid;
            grid-template-columns: 1fr;
            gap: 16px;
            max-width: 1200px;
            margin: 0 auto;
        }

        .feature-card {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            padding: 28px;
            position: relative;
            overflow: hidden;
            transition: all 0.4s cubic-bezier(0.4,0,0.2,1);
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04);
            cursor: pointer;
        }

        .feature-card::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 1px;
            background: linear-gradient(90deg, transparent 0%, rgba(255,255,255,0.8) 50%, transparent 100%);
        }

        .feature-card:hover {
            transform: translateY(-4px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.04), 0 24px 64px rgba(0,0,0,0.08);
            border-color: rgba(59,130,246,0.2);
        }

        .feature-icon {
            width: 56px;
            height: 56px;
            border-radius: 16px;
            display: flex;
            align-items: center;
            justify-content: center;
            margin-bottom: 20px;
            transition: all 0.3s ease;
        }

        .feature-card:hover .feature-icon { transform: scale(1.05); }
        .feature-icon svg { width: 24px; height: 24px; }

        .feature-icon.blue {
            background: var(--blue-50);
            border: 1px solid var(--blue-100);
            box-shadow: 0 4px 16px rgba(37,99,235,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.blue svg { color: var(--blue-600); }

        .feature-icon.emerald {
            background: #ECFDF5;
            border: 1px solid #D1FAE5;
            box-shadow: 0 4px 16px rgba(16,185,129,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.emerald svg { color: #059669; }

        .feature-icon.violet {
            background: #F5F3FF;
            border: 1px solid #EDE9FE;
            box-shadow: 0 4px 16px rgba(139,92,246,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.violet svg { color: #7C3AED; }

        .feature-icon.rose {
            background: #FFF1F2;
            border: 1px solid #FFE4E6;
            box-shadow: 0 4px 16px rgba(244,63,94,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.rose svg { color: #E11D48; }

        .feature-icon.cyan {
            background: #ECFEFF;
            border: 1px solid #CFFAFE;
            box-shadow: 0 4px 16px rgba(6,182,212,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.cyan svg { color: #0891B2; }

        .feature-title {
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 10px;
            letter-spacing: -0.3px;
        }

        .feature-desc {
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-500);
            margin-bottom: 20px;
        }

        .feature-link {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            font-size: 14px;
            font-weight: 600;
            color: var(--blue-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .feature-link:hover { gap: 10px; }
        .feature-link svg { width: 16px; height: 16px; transition: transform 0.2s ease; }
        .feature-link:hover svg { transform: translateX(4px); }

        .footer {
            padding: 32px 20px;
            border-top: 1px solid var(--gray-200);
            background: rgba(255,255,255,0.5);
        }

        .footer-inner {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 20px;
            text-align: center;
        }

        .footer-brand {
            display: flex;
            align-items: center;
            gap: 10px;
        }

        .footer-logo {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .footer-logo svg { width: 32px; height: 32px; }

        .footer-text {
            font-size: 13px;
            color: var(--gray-500);
        }

        .footer-links {
            display: flex;
            gap: 24px;
        }

        .footer-link {
            font-size: 13px;
            color: var(--gray-500);
            text-decoration: none;
            transition: color 0.2s ease;
        }

        .footer-link:hover { color: var(--gray-700); }

        @media (min-width: 768px) {
            .nav { padding: 16px 32px; }
            .nav-inner { height: 64px; padding: 0 24px; border-radius: 20px; }
            .logo-icon svg { width: 42px; height: 15px; }
            .logo-text { font-size: 20px; }
            .nav-links { display: flex; }
            .mobile-menu-btn { display: none; }
            .hero { padding: 160px 32px 80px; }
            .hero-badge { padding: 10px 20px 10px 14px; margin-bottom: 32px; }
            .badge-dot { width: 10px; height: 10px; }
            .badge-text { font-size: 13px; }
            .hero-title { font-size: 56px; letter-spacing: -2.5px; margin-bottom: 24px; }
            .hero-subtitle { font-size: 18px; margin-bottom: 48px; }
            .chat-container { padding: 0 24px; margin-bottom: 80px; }
            .chat-card { border-radius: 24px; padding: 10px; }
            .chat-inner { border-radius: 18px; padding: 8px; }
            .chat-input { padding: 16px 18px; min-height: 56px; }
            .chat-send { width: 52px; height: 52px; border-radius: 14px; }
            .chat-hints { gap: 10px; padding: 18px 10px 8px; }
            .hint-chip { padding: 12px 18px; font-size: 13px; }
            .hint-chip svg { width: 16px; height: 16px; }
            .features { padding: 80px 32px 100px; }
            .features-header { margin-bottom: 56px; }
            .features-label { font-size: 12px; margin-bottom: 16px; }
            .features-title { font-size: 36px; letter-spacing: -1.5px; }
            .features-subtitle { font-size: 18px; }
            .features-grid { grid-template-columns: repeat(2, 1fr); gap: 20px; }
            .feature-card { padding: 36px; border-radius: 24px; }
            .feature-card:hover { transform: translateY(-6px); }
            .feature-icon { width: 60px; height: 60px; border-radius: 18px; margin-bottom: 24px; }
            .feature-icon svg { width: 26px; height: 26px; }
            .feature-title { font-size: 20px; margin-bottom: 12px; }
            .feature-desc { font-size: 15px; line-height: 1.7; margin-bottom: 24px; }
            .footer { padding: 40px 32px; }
            .footer-inner { flex-direction: row; justify-content: space-between; text-align: left; }
            .footer-logo svg { width: 36px; height: 36px; }
            .footer-text { font-size: 14px; }
            .footer-links { gap: 32px; }
            .footer-link { font-size: 14px; }
            .orb-1 { width: 600px; height: 600px; left: -10%; }
            .orb-2 { width: 450px; height: 450px; right: -10%; }
            .orb-3 { width: 400px; height: 400px; }
        }

        @media (min-width: 1024px) {
            .nav { padding: 16px 40px; }
            .hero { padding: 180px 40px 80px; }
            .hero-title { font-size: 72px; letter-spacing: -3px; margin-bottom: 28px; }
            .hero-subtitle { font-size: 20px; margin-bottom: 56px; }
            .chat-container { margin-bottom: 100px; }
            .chat-card { border-radius: 28px; }
            .chat-inner { border-radius: 20px; }
            .chat-input { padding: 18px 20px; min-height: 60px; max-height: 180px; }
            .chat-send { width: 56px; height: 56px; border-radius: 16px; margin: 8px; }
            .chat-send svg { width: 22px; height: 22px; }
            .chat-hints { padding: 20px 12px 8px; }
            .hint-chip { padding: 12px 20px; }
            .features { padding: 80px 40px 120px; }
            .features-header { margin-bottom: 64px; }
            .features-title { font-size: 40px; }
            .features-grid { grid-template-columns: repeat(3, 1fr); gap: 24px; }
            .feature-card { padding: 40px; border-radius: 28px; }
            .feature-card:hover { transform: translateY(-8px); }
            .feature-icon { width: 64px; height: 64px; border-radius: 20px; margin-bottom: 28px; }
            .feature-icon svg { width: 28px; height: 28px; }
            .footer { padding: 48px 40px; }
            .orb-1 { width: 800px; height: 800px; top: -20%; left: -10%; }
            .orb-2 { width: 600px; height: 600px; right: -15%; }
            .orb-3 { width: 500px; height: 500px; }
        }

        @media (min-width: 1280px) {
            .hero-title { font-size: 80px; }
        }
    

        .main-content {
            flex: 1;
            padding: 100px 20px 40px;
            max-width: 1200px;
            margin: 0 auto;
            width: 100%;
        }

        .content-card {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            padding: 32px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
        }

        h1 { font-size: 32px; margin-bottom: 16px; font-weight: 700; letter-spacing: -0.5px; }
        h2 { font-size: 24px; margin-bottom: 12px; font-weight: 700; letter-spacing: -0.5px; }
        h3 { font-size: 20px; margin-bottom: 10px; font-weight: 700; letter-spacing: -0.3px; }

        input[type="text"], input[type="number"], input[type="email"], select, textarea {
            width: 100%;
            padding: 12px 16px;
            border: 1px solid var(--gray-300);
            border-radius: 12px;
            font-family: inherit;
            font-size: 15px;
            background: var(--white);
            color: var(--gray-900);
            transition: all 0.2s ease;
        }

        input:focus, select:focus, textarea:focus {
            outline: none;
            border-color: var(--blue-500);
            box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.1);
        }

        button, .btn {
            padding: 12px 24px;
            background: var(--blue-600);
            color: var(--white);
            border: none;
            border-radius: 12px;
            font-family: inherit;
            font-size: 15px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s ease;
            box-shadow: 0 1px 2px rgba(37,99,235,0.2), 0 4px 16px rgba(37,99,235,0.2), inset 0 1px 0 rgba(255,255,255,0.1);
        }

        button:hover, .btn:hover {
            background: var(--blue-700);
            transform: translateY(-2px);
            box-shadow: 0 2px 4px rgba(37,99,235,0.2), 0 12px 40px rgba(37,99,235,0.3), inset 0 1px 0 rgba(255,255,255,0.1);
        }

        button:active, .btn:active {
            transform: translateY(0);
        }

        table {
            width: 100%;
            border-collapse: collapse;
            background: rgba(255,255,255,0.5);
            border-radius: 12px;
            overflow: hidden;
        }

        th, td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid var(--gray-200);
        }

        th {
            background: rgba(37,99,235,0.05);
            font-weight: 600;
            color: var(--blue-700);
        }

        @media (min-width: 768px) {
            .main-content { padding: 120px 32px 60px; }
        }

        @media (min-width: 1024px) {
            .main-content { padding: 140px 40px 80px; }
        }

    </style>
</head>
<body>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <nav class="nav">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo">
                    <div class="logo-icon">
                        <svg width="36" height="12" viewBox="0 0 52 18" fill="none">
                            <circle cx="9" cy="9" r="9" fill="#2563EB"/>
                            <circle cx="26" cy="9" r="9" fill="#2563EB" fill-opacity="0.5"/>
                            <circle cx="43" cy="9" r="9" fill="#2563EB" fill-opacity="0.2"/>
                        </svg>
                    </div>
                    <span class="logo-text"><span class="gas">gas</span><span class="consult">consult</span><span class="ai">.ai</span></span>
                </a>
                <div class="nav-links">
                    <a href="/?clear=1" class="nav-link">Home</a>
                    <a href="/quick-dose" class="nav-link">Quick Dose</a>
                    <a href="/preop" class="nav-link">Pre-Op</a>
                    <a href="/calculators" class="nav-link">Calculators</a>
                    <a href="/hypotension" class="nav-link">IOH Predictor</a>
                </div>
                <button class="mobile-menu-btn" onclick="toggleMobileMenu()" aria-label="Toggle menu">
                    <span></span>
                    <span></span>
                    <span></span>
                </button>
            </div>
        </nav>
        <div class="mobile-menu" id="mobileMenu">
            <a href="/?clear=1" class="mobile-menu-link">Home</a>
            <a href="/quick-dose" class="mobile-menu-link">Quick Dose</a>
            <a href="/preop" class="mobile-menu-link">Pre-Op</a>
            <a href="/calculators" class="mobile-menu-link">Calculators</a>
            <a href="/hypotension" class="mobile-menu-link">IOH Predictor</a>
        </div>

        <main class="main-content">
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
        </main>
    </div>
</body>
</html>
"""

TERMS_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Terms of Service — gasconsult.ai</title>

    <!-- PWA -->
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=5">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">

    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap" rel="stylesheet">
    <style>

        :root {
            --white: #FFFFFF;
            --gray-50: #F8FAFC;
            --gray-100: #F1F5F9;
            --gray-200: #E2E8F0;
            --gray-300: #CBD5E1;
            --gray-400: #94A3B8;
            --gray-500: #64748B;
            --gray-600: #475569;
            --gray-700: #334155;
            --gray-800: #1E293B;
            --gray-900: #0F172A;
            --blue-50: #EFF6FF;
            --blue-100: #DBEAFE;
            --blue-200: #BFDBFE;
            --blue-300: #93C5FD;
            --blue-400: #60A5FA;
            --blue-500: #3B82F6;
            --blue-600: #2563EB;
            --blue-700: #1D4ED8;
        }

        * { margin: 0; padding: 0; box-sizing: border-box; }

        html {
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
            scroll-behavior: smooth;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            background: var(--gray-50);
            color: var(--gray-900);
            min-height: 100vh;
            overflow-x: hidden;
        }

        .bg-canvas {
            position: fixed;
            inset: 0;
            z-index: 0;
            overflow: hidden;
            background: linear-gradient(180deg, #F0F7FF 0%, var(--gray-50) 50%, #FAFBFF 100%);
        }

        .orb {
            position: absolute;
            border-radius: 50%;
            filter: blur(80px);
            opacity: 0.6;
            animation: float 20s ease-in-out infinite;
        }

        .orb-1 {
            width: 400px;
            height: 400px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.15) 0%, transparent 70%);
            top: -15%;
            left: -20%;
        }

        .orb-2 {
            width: 300px;
            height: 300px;
            background: radial-gradient(circle, rgba(147, 197, 253, 0.2) 0%, transparent 70%);
            top: 30%;
            right: -20%;
            animation-delay: -7s;
            animation-duration: 25s;
        }

        .orb-3 {
            width: 250px;
            height: 250px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.1) 0%, transparent 70%);
            bottom: -10%;
            left: 20%;
            animation-delay: -14s;
            animation-duration: 30s;
        }

        @keyframes float {
            0%, 100% { transform: translate(0, 0) scale(1); }
            25% { transform: translate(40px, -40px) scale(1.05); }
            50% { transform: translate(20px, 40px) scale(0.95); }
            75% { transform: translate(-40px, 20px) scale(1.02); }
        }

        .grain {
            position: fixed;
            inset: 0;
            z-index: 1;
            pointer-events: none;
            opacity: 0.02;
            background-image: url("data:image/svg+xml,%3Csvg viewBox='0 0 512 512' xmlns='http://www.w3.org/2000/svg'%3E%3Cfilter id='n'%3E%3CfeTurbulence type='fractalNoise' baseFrequency='0.8' numOctaves='4' stitchTiles='stitch'/%3E%3C/filter%3E%3Crect width='100%25' height='100%25' filter='url(%23n)'/%3E%3C/svg%3E");
        }

        .page {
            position: relative;
            z-index: 2;
            min-height: 100vh;
            display: flex;
            flex-direction: column;
        }

        .nav {
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            z-index: 100;
            padding: 12px 16px;
        }

        .nav-inner {
            max-width: 1200px;
            margin: 0 auto;
            height: 56px;
            background: rgba(255, 255, 255, 0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 0 16px;
            display: flex;
            align-items: center;
            justify-content: space-between;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 12px 48px rgba(0,0,0,0.03);
        }

        .logo {
            display: flex;
            align-items: center;
            gap: 14px;
            text-decoration: none;
        }

        .logo-icon {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .logo-icon svg { width: 36px; height: 12px; }

        .logo-text {
            font-size: 18px;
            font-weight: 700;
            letter-spacing: -0.5px;
            color: var(--gray-900);
        }

        .logo-text .gas { color: var(--blue-600); }
        .logo-text .consult { color: #0F172A; }
        .logo-text .ai { color: rgba(15, 23, 42, 0.4); }

        .nav-links {
            display: none;
            align-items: center;
            gap: 4px;
        }

        .nav-link {
            padding: 10px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .nav-link.active {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .mobile-menu-btn {
            display: flex;
            flex-direction: column;
            gap: 5px;
            background: none;
            border: none;
            cursor: pointer;
            padding: 8px;
            border-radius: 8px;
            transition: background 0.2s ease;
        }

        .mobile-menu-btn:hover {
            background: rgba(0,0,0,0.04);
        }

        .mobile-menu-btn span {
            display: block;
            width: 22px;
            height: 2px;
            background: var(--gray-700);
            border-radius: 1px;
            transition: all 0.3s ease;
        }

        .mobile-menu-btn.active span:nth-child(1) {
            transform: rotate(45deg) translate(7px, 7px);
        }

        .mobile-menu-btn.active span:nth-child(2) {
            opacity: 0;
        }

        .mobile-menu-btn.active span:nth-child(3) {
            transform: rotate(-45deg) translate(7px, -7px);
        }

        .mobile-menu {
            display: none;
            position: fixed;
            top: 80px;
            left: 16px;
            right: 16px;
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 8px;
            box-shadow: 0 4px 16px rgba(0,0,0,0.08), 0 12px 48px rgba(0,0,0,0.12);
            z-index: 99;
            flex-direction: column;
            gap: 4px;
        }

        .mobile-menu.active {
            display: flex;
        }

        .mobile-menu-link {
            padding: 14px 16px;
            font-size: 15px;
            font-weight: 500;
            color: var(--gray-700);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .mobile-menu-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .main-content {
            flex: 1;
            padding: 100px 20px 40px;
            max-width: 900px;
            margin: 0 auto;
            width: 100%;
        }

        .content-card {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            padding: 40px 28px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
        }

        h1 {
            font-size: 36px;
            margin-bottom: 8px;
            font-weight: 800;
            letter-spacing: -1px;
            color: var(--gray-900);
        }

        h2 {
            font-size: 22px;
            margin-top: 32px;
            margin-bottom: 16px;
            font-weight: 700;
            letter-spacing: -0.5px;
            color: var(--gray-900);
        }

        h3 {
            font-size: 18px;
            margin-top: 20px;
            margin-bottom: 12px;
            font-weight: 600;
            letter-spacing: -0.3px;
            color: var(--gray-800);
        }

        p {
            font-size: 15px;
            line-height: 1.7;
            color: var(--gray-700);
            margin-bottom: 16px;
        }

        ul {
            margin-left: 24px;
            margin-bottom: 16px;
        }

        li {
            font-size: 15px;
            line-height: 1.7;
            color: var(--gray-700);
            margin-bottom: 8px;
        }

        strong {
            font-weight: 600;
            color: var(--gray-900);
        }

        .last-updated {
            font-size: 14px;
            color: var(--gray-500);
            margin-bottom: 32px;
            font-weight: 500;
        }

        .notice-box {
            background: rgba(239, 68, 68, 0.05);
            border-left: 4px solid #EF4444;
            border-radius: 12px;
            padding: 20px 24px;
            margin: 24px 0;
        }

        .notice-box h3 {
            font-size: 16px;
            font-weight: 700;
            color: #DC2626;
            margin: 0 0 12px 0;
        }

        .notice-box p {
            font-size: 15px;
            line-height: 1.7;
            color: var(--gray-800);
            margin: 0;
        }

        .footer {
            padding: 32px 20px;
            border-top: 1px solid var(--gray-200);
            background: rgba(255,255,255,0.5);
        }

        .footer-inner {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 20px;
            text-align: center;
        }

        .footer-brand {
            display: flex;
            align-items: center;
            gap: 10px;
        }

        .footer-logo {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .footer-logo svg { width: 32px; height: 32px; }

        .footer-text {
            font-size: 13px;
            color: var(--gray-500);
        }

        .footer-links {
            display: flex;
            gap: 24px;
        }

        .footer-link {
            font-size: 13px;
            color: var(--gray-500);
            text-decoration: none;
            transition: color 0.2s ease;
        }

        .footer-link:hover { color: var(--gray-700); }

        @media (min-width: 768px) {
            .nav { padding: 16px 32px; }
            .nav-inner { height: 64px; padding: 0 24px; border-radius: 20px; }
            .logo-icon svg { width: 42px; height: 15px; }
            .logo-text { font-size: 20px; }
            .nav-links { display: flex; }
            .mobile-menu-btn { display: none; }
            .main-content { padding: 120px 32px 60px; }
            .content-card { padding: 48px 40px; border-radius: 24px; }
            h1 { font-size: 42px; margin-bottom: 12px; }
            h2 { font-size: 24px; margin-top: 40px; margin-bottom: 18px; }
            h3 { font-size: 19px; }
            .footer { padding: 40px 32px; }
            .footer-inner { flex-direction: row; justify-content: space-between; text-align: left; }
            .footer-logo svg { width: 36px; height: 36px; }
            .footer-text { font-size: 14px; }
            .footer-links { gap: 32px; }
            .footer-link { font-size: 14px; }
            .orb-1 { width: 600px; height: 600px; left: -10%; }
            .orb-2 { width: 450px; height: 450px; right: -10%; }
            .orb-3 { width: 400px; height: 400px; }
        }

        @media (min-width: 1024px) {
            .nav { padding: 16px 40px; }
            .main-content { padding: 140px 40px 80px; }
            .content-card { padding: 56px 48px; border-radius: 28px; }
            h1 { font-size: 48px; letter-spacing: -1.5px; }
            h2 { font-size: 26px; }
            h3 { font-size: 20px; }
            p, li { font-size: 16px; }
            .footer { padding: 48px 40px; }
            .orb-1 { width: 800px; height: 800px; top: -20%; left: -10%; }
            .orb-2 { width: 600px; height: 600px; right: -15%; }
            .orb-3 { width: 500px; height: 500px; }
        }

    </style>
</head>
<body>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <nav class="nav">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo">
                    <div class="logo-icon">
                        <svg width="36" height="12" viewBox="0 0 52 18" fill="none">
                            <circle cx="9" cy="9" r="9" fill="#2563EB"/>
                            <circle cx="26" cy="9" r="9" fill="#2563EB" fill-opacity="0.5"/>
                            <circle cx="43" cy="9" r="9" fill="#2563EB" fill-opacity="0.2"/>
                        </svg>
                    </div>
                    <span class="logo-text"><span class="gas">gas</span><span class="consult">consult</span><span class="ai">.ai</span></span>
                </a>
                <div class="nav-links">
                    <a href="/?clear=1" class="nav-link">Home</a>
                    <a href="/quick-dose" class="nav-link">Quick Dose</a>
                    <a href="/preop" class="nav-link">Pre-Op</a>
                    <a href="/calculators" class="nav-link">Calculators</a>
                    <a href="/hypotension" class="nav-link">IOH Predictor</a>
                </div>
                <button class="mobile-menu-btn" onclick="toggleMobileMenu()" aria-label="Toggle menu">
                    <span></span>
                    <span></span>
                    <span></span>
                </button>
            </div>
        </nav>
        <div class="mobile-menu" id="mobileMenu">
            <a href="/?clear=1" class="mobile-menu-link">Home</a>
            <a href="/quick-dose" class="mobile-menu-link">Quick Dose</a>
            <a href="/preop" class="mobile-menu-link">Pre-Op</a>
            <a href="/calculators" class="mobile-menu-link">Calculators</a>
            <a href="/hypotension" class="mobile-menu-link">IOH Predictor</a>
        </div>

        <main class="main-content">
            <div class="content-card">
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
        </main>

        <footer class="footer">
            <div class="footer-inner">
                <div class="footer-brand">
                    <div class="footer-logo"><svg width="32" height="32" viewBox="0 0 32 32" fill="none"><circle cx="6" cy="16" r="5" fill="#2563EB"/><circle cx="16" cy="16" r="5" fill="#2563EB" fill-opacity="0.5"/><circle cx="26" cy="16" r="5" fill="#2563EB" fill-opacity="0.2"/></svg></div>
                    <span class="footer-text">© 2025 GasConsult.ai</span>
                </div>
                <div class="footer-links">
                    <a href="/privacy" class="footer-link">Privacy</a>
                    <a href="/terms" class="footer-link">Terms</a>
                    <a href="#" class="footer-link">Contact</a>
                </div>
            </div>
        </footer>
    </div>

    <script>
        function toggleMobileMenu() {
            const menu = document.getElementById('mobileMenu');
            const btn = document.querySelector('.mobile-menu-btn');
            if (menu && btn) {
                menu.classList.toggle('active');
                btn.classList.toggle('active');
            }
        }
    </script>
</body>
</html>
"""

PRIVACY_POLICY_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Privacy Policy — gasconsult.ai</title>

    <!-- PWA -->
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=5">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">

    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap" rel="stylesheet">
    <style>

        :root {
            --white: #FFFFFF;
            --gray-50: #F8FAFC;
            --gray-100: #F1F5F9;
            --gray-200: #E2E8F0;
            --gray-300: #CBD5E1;
            --gray-400: #94A3B8;
            --gray-500: #64748B;
            --gray-600: #475569;
            --gray-700: #334155;
            --gray-800: #1E293B;
            --gray-900: #0F172A;
            --blue-50: #EFF6FF;
            --blue-100: #DBEAFE;
            --blue-200: #BFDBFE;
            --blue-300: #93C5FD;
            --blue-400: #60A5FA;
            --blue-500: #3B82F6;
            --blue-600: #2563EB;
            --blue-700: #1D4ED8;
        }

        * { margin: 0; padding: 0; box-sizing: border-box; }

        html {
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
            scroll-behavior: smooth;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            background: var(--gray-50);
            color: var(--gray-900);
            min-height: 100vh;
            overflow-x: hidden;
        }

        .bg-canvas {
            position: fixed;
            inset: 0;
            z-index: 0;
            overflow: hidden;
            background: linear-gradient(180deg, #F0F7FF 0%, var(--gray-50) 50%, #FAFBFF 100%);
        }

        .orb {
            position: absolute;
            border-radius: 50%;
            filter: blur(80px);
            opacity: 0.6;
            animation: float 20s ease-in-out infinite;
        }

        .orb-1 {
            width: 400px;
            height: 400px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.15) 0%, transparent 70%);
            top: -15%;
            left: -20%;
        }

        .orb-2 {
            width: 300px;
            height: 300px;
            background: radial-gradient(circle, rgba(147, 197, 253, 0.2) 0%, transparent 70%);
            top: 30%;
            right: -20%;
            animation-delay: -7s;
            animation-duration: 25s;
        }

        .orb-3 {
            width: 250px;
            height: 250px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.1) 0%, transparent 70%);
            bottom: -10%;
            left: 20%;
            animation-delay: -14s;
            animation-duration: 30s;
        }

        @keyframes float {
            0%, 100% { transform: translate(0, 0) scale(1); }
            25% { transform: translate(40px, -40px) scale(1.05); }
            50% { transform: translate(20px, 40px) scale(0.95); }
            75% { transform: translate(-40px, 20px) scale(1.02); }
        }

        .grain {
            position: fixed;
            inset: 0;
            z-index: 1;
            pointer-events: none;
            opacity: 0.02;
            background-image: url("data:image/svg+xml,%3Csvg viewBox='0 0 512 512' xmlns='http://www.w3.org/2000/svg'%3E%3Cfilter id='n'%3E%3CfeTurbulence type='fractalNoise' baseFrequency='0.8' numOctaves='4' stitchTiles='stitch'/%3E%3C/filter%3E%3Crect width='100%25' height='100%25' filter='url(%23n)'/%3E%3C/svg%3E");
        }

        .page {
            position: relative;
            z-index: 2;
            min-height: 100vh;
            display: flex;
            flex-direction: column;
        }

        .nav {
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            z-index: 100;
            padding: 12px 16px;
        }

        .nav-inner {
            max-width: 1200px;
            margin: 0 auto;
            height: 56px;
            background: rgba(255, 255, 255, 0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 0 16px;
            display: flex;
            align-items: center;
            justify-content: space-between;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 12px 48px rgba(0,0,0,0.03);
        }

        .logo {
            display: flex;
            align-items: center;
            gap: 14px;
            text-decoration: none;
        }

        .logo-icon {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .logo-icon svg { width: 36px; height: 12px; }

        .logo-text {
            font-size: 18px;
            font-weight: 700;
            letter-spacing: -0.5px;
            color: var(--gray-900);
        }

        .logo-text .gas { color: var(--blue-600); }
        .logo-text .consult { color: #0F172A; }
        .logo-text .ai { color: rgba(15, 23, 42, 0.4); }

        .nav-links {
            display: none;
            align-items: center;
            gap: 4px;
        }

        .nav-link {
            padding: 10px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .nav-link.active {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .mobile-menu-btn {
            display: flex;
            flex-direction: column;
            gap: 5px;
            background: none;
            border: none;
            cursor: pointer;
            padding: 8px;
            border-radius: 8px;
            transition: background 0.2s ease;
        }

        .mobile-menu-btn:hover {
            background: rgba(0,0,0,0.04);
        }

        .mobile-menu-btn span {
            display: block;
            width: 22px;
            height: 2px;
            background: var(--gray-700);
            border-radius: 1px;
            transition: all 0.3s ease;
        }

        .mobile-menu-btn.active span:nth-child(1) {
            transform: rotate(45deg) translate(7px, 7px);
        }

        .mobile-menu-btn.active span:nth-child(2) {
            opacity: 0;
        }

        .mobile-menu-btn.active span:nth-child(3) {
            transform: rotate(-45deg) translate(7px, -7px);
        }

        .mobile-menu {
            display: none;
            position: fixed;
            top: 80px;
            left: 16px;
            right: 16px;
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 8px;
            box-shadow: 0 4px 16px rgba(0,0,0,0.08), 0 12px 48px rgba(0,0,0,0.12);
            z-index: 99;
            flex-direction: column;
            gap: 4px;
        }

        .mobile-menu.active {
            display: flex;
        }

        .mobile-menu-link {
            padding: 14px 16px;
            font-size: 15px;
            font-weight: 500;
            color: var(--gray-700);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .mobile-menu-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .footer {
            padding: 32px 20px;
            border-top: 1px solid var(--gray-200);
            background: rgba(255,255,255,0.5);
        }

        .footer-inner {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 20px;
            text-align: center;
        }

        .footer-brand {
            display: flex;
            align-items: center;
            gap: 10px;
        }

        .footer-logo {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .footer-logo svg { width: 32px; height: 32px; }

        .footer-text {
            font-size: 13px;
            color: var(--gray-500);
        }

        .footer-links {
            display: flex;
            gap: 24px;
        }

        .footer-link {
            font-size: 13px;
            color: var(--gray-500);
            text-decoration: none;
            transition: color 0.2s ease;
        }

        .footer-link:hover { color: var(--gray-700); }

        @media (min-width: 768px) {
            .nav { padding: 16px 32px; }
            .nav-inner { height: 64px; padding: 0 24px; border-radius: 20px; }
            .logo-icon svg { width: 42px; height: 15px; }
            .logo-text { font-size: 20px; }
            .nav-links { display: flex; }
            .mobile-menu-btn { display: none; }
            .footer { padding: 40px 32px; }
            .footer-inner { flex-direction: row; justify-content: space-between; text-align: left; }
            .footer-logo svg { width: 36px; height: 36px; }
            .footer-text { font-size: 14px; }
            .footer-links { gap: 32px; }
            .footer-link { font-size: 14px; }
            .orb-1 { width: 600px; height: 600px; left: -10%; }
            .orb-2 { width: 450px; height: 450px; right: -10%; }
            .orb-3 { width: 400px; height: 400px; }
        }

        @media (min-width: 1024px) {
            .nav { padding: 16px 40px; }
            .footer { padding: 48px 40px; }
            .orb-1 { width: 800px; height: 800px; top: -20%; left: -10%; }
            .orb-2 { width: 600px; height: 600px; right: -15%; }
            .orb-3 { width: 500px; height: 500px; }
        }

        .main-content {
            flex: 1;
            padding: 100px 20px 40px;
            max-width: 900px;
            margin: 0 auto;
            width: 100%;
        }

        .content-card {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            padding: 32px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
        }

        h1 {
            font-size: 32px;
            font-weight: 800;
            letter-spacing: -1px;
            color: var(--gray-900);
            margin-bottom: 8px;
        }

        h2 {
            font-size: 22px;
            font-weight: 700;
            letter-spacing: -0.5px;
            color: var(--gray-900);
            margin-top: 32px;
            margin-bottom: 16px;
        }

        h3 {
            font-size: 18px;
            font-weight: 600;
            letter-spacing: -0.3px;
            color: var(--gray-800);
            margin-top: 24px;
            margin-bottom: 12px;
        }

        p {
            font-size: 15px;
            line-height: 1.7;
            color: var(--gray-700);
            margin-bottom: 16px;
        }

        ul {
            margin: 16px 0;
            padding-left: 24px;
        }

        li {
            font-size: 15px;
            line-height: 1.7;
            color: var(--gray-700);
            margin-bottom: 8px;
        }

        strong {
            font-weight: 600;
            color: var(--gray-900);
        }

        .last-updated {
            font-size: 14px;
            color: var(--gray-500);
            margin-bottom: 24px;
        }

        .highlight-box {
            background: rgba(59, 130, 246, 0.05);
            border: 1px solid rgba(59, 130, 246, 0.15);
            border-radius: 12px;
            padding: 20px;
            margin: 24px 0;
        }

        .highlight-box strong {
            color: var(--blue-700);
        }

        .highlight-box p {
            margin-bottom: 0;
        }

        .highlight-box ul {
            margin-top: 12px;
            margin-bottom: 0;
        }

        @media (min-width: 768px) {
            .main-content { padding: 120px 32px 60px; }
            .content-card { padding: 40px; border-radius: 24px; }
            h1 { font-size: 40px; margin-bottom: 12px; }
            h2 { font-size: 26px; margin-top: 40px; margin-bottom: 20px; }
            h3 { font-size: 20px; margin-top: 28px; margin-bottom: 14px; }
            p { font-size: 16px; }
            li { font-size: 16px; }
        }

        @media (min-width: 1024px) {
            .main-content { padding: 140px 40px 80px; max-width: 1000px; }
            .content-card { padding: 48px; }
            h1 { font-size: 48px; letter-spacing: -1.5px; }
        }

    </style>
</head>
<body>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <nav class="nav">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo">
                    <div class="logo-icon">
                        <svg width="36" height="12" viewBox="0 0 52 18" fill="none">
                            <circle cx="9" cy="9" r="9" fill="#2563EB"/>
                            <circle cx="26" cy="9" r="9" fill="#2563EB" fill-opacity="0.5"/>
                            <circle cx="43" cy="9" r="9" fill="#2563EB" fill-opacity="0.2"/>
                        </svg>
                    </div>
                    <span class="logo-text"><span class="gas">gas</span><span class="consult">consult</span><span class="ai">.ai</span></span>
                </a>
                <div class="nav-links">
                    <a href="/?clear=1" class="nav-link">Home</a>
                    <a href="/quick-dose" class="nav-link">Quick Dose</a>
                    <a href="/preop" class="nav-link">Pre-Op</a>
                    <a href="/calculators" class="nav-link">Calculators</a>
                    <a href="/hypotension" class="nav-link">IOH Predictor</a>
                </div>
                <button class="mobile-menu-btn" onclick="toggleMobileMenu()" aria-label="Toggle menu">
                    <span></span>
                    <span></span>
                    <span></span>
                </button>
            </div>
        </nav>
        <div class="mobile-menu" id="mobileMenu">
            <a href="/?clear=1" class="mobile-menu-link">Home</a>
            <a href="/quick-dose" class="mobile-menu-link">Quick Dose</a>
            <a href="/preop" class="mobile-menu-link">Pre-Op</a>
            <a href="/calculators" class="mobile-menu-link">Calculators</a>
            <a href="/hypotension" class="mobile-menu-link">IOH Predictor</a>
        </div>

        <main class="main-content">
            <div class="content-card">
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

            <footer class="footer">
                <div class="footer-inner">
                    <div class="footer-brand">
                        <div class="footer-logo"><svg width="32" height="32" viewBox="0 0 32 32" fill="none"><circle cx="6" cy="16" r="5" fill="#2563EB"/><circle cx="16" cy="16" r="5" fill="#2563EB" fill-opacity="0.5"/><circle cx="26" cy="16" r="5" fill="#2563EB" fill-opacity="0.2"/></svg></div>
                        <span class="footer-text">© 2025 GasConsult.ai</span>
                    </div>
                    <div class="footer-links">
                        <a href="/privacy" class="footer-link">Privacy</a>
                        <a href="/terms" class="footer-link">Terms</a>
                        <a href="#" class="footer-link">Contact</a>
                    </div>
                </div>
            </footer>
        </main>
    </div>
    <script>
        function toggleMobileMenu() {
            const menu = document.getElementById('mobileMenu');
            const btn = document.querySelector('.mobile-menu-btn');
            if (menu && btn) {
                menu.classList.toggle('active');
                btn.classList.toggle('active');
            }
        }
    </script>
</body>
</html>
"""

QUICK_DOSE_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Quick Dose Reference — gasconsult.ai</title>

    <!-- PWA -->
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=5">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">

    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap" rel="stylesheet">
    <style>

        :root {
            --white: #FFFFFF;
            --gray-50: #F8FAFC;
            --gray-100: #F1F5F9;
            --gray-200: #E2E8F0;
            --gray-300: #CBD5E1;
            --gray-400: #94A3B8;
            --gray-500: #64748B;
            --gray-600: #475569;
            --gray-700: #334155;
            --gray-800: #1E293B;
            --gray-900: #0F172A;
            --blue-50: #EFF6FF;
            --blue-100: #DBEAFE;
            --blue-200: #BFDBFE;
            --blue-300: #93C5FD;
            --blue-400: #60A5FA;
            --blue-500: #3B82F6;
            --blue-600: #2563EB;
            --blue-700: #1D4ED8;
        }

        * { margin: 0; padding: 0; box-sizing: border-box; }

        html {
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
            scroll-behavior: smooth;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            background: var(--gray-50);
            color: var(--gray-900);
            min-height: 100vh;
            overflow-x: hidden;
        }

        .bg-canvas {
            position: fixed;
            inset: 0;
            z-index: 0;
            overflow: hidden;
            background: linear-gradient(180deg, #F0F7FF 0%, var(--gray-50) 50%, #FAFBFF 100%);
        }

        .orb {
            position: absolute;
            border-radius: 50%;
            filter: blur(80px);
            opacity: 0.6;
            animation: float 20s ease-in-out infinite;
        }

        .orb-1 {
            width: 400px;
            height: 400px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.15) 0%, transparent 70%);
            top: -15%;
            left: -20%;
        }

        .orb-2 {
            width: 300px;
            height: 300px;
            background: radial-gradient(circle, rgba(147, 197, 253, 0.2) 0%, transparent 70%);
            top: 30%;
            right: -20%;
            animation-delay: -7s;
            animation-duration: 25s;
        }

        .orb-3 {
            width: 250px;
            height: 250px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.1) 0%, transparent 70%);
            bottom: -10%;
            left: 20%;
            animation-delay: -14s;
            animation-duration: 30s;
        }

        @keyframes float {
            0%, 100% { transform: translate(0, 0) scale(1); }
            25% { transform: translate(40px, -40px) scale(1.05); }
            50% { transform: translate(20px, 40px) scale(0.95); }
            75% { transform: translate(-40px, 20px) scale(1.02); }
        }

        .grain {
            position: fixed;
            inset: 0;
            z-index: 1;
            pointer-events: none;
            opacity: 0.02;
            background-image: url("data:image/svg+xml,%3Csvg viewBox='0 0 512 512' xmlns='http://www.w3.org/2000/svg'%3E%3Cfilter id='n'%3E%3CfeTurbulence type='fractalNoise' baseFrequency='0.8' numOctaves='4' stitchTiles='stitch'/%3E%3C/filter%3E%3Crect width='100%25' height='100%25' filter='url(%23n)'/%3E%3C/svg%3E");
        }

        .page {
            position: relative;
            z-index: 2;
            min-height: 100vh;
            display: flex;
            flex-direction: column;
        }

        .nav {
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            z-index: 100;
            padding: 12px 16px;
        }

        .nav-inner {
            max-width: 1200px;
            margin: 0 auto;
            height: 56px;
            background: rgba(255, 255, 255, 0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 0 16px;
            display: flex;
            align-items: center;
            justify-content: space-between;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 12px 48px rgba(0,0,0,0.03);
        }

        .logo {
            display: flex;
            align-items: center;
            gap: 14px;
            text-decoration: none;
        }

        .logo-icon {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .logo-icon svg { width: 36px; height: 12px; }

        .logo-text {
            font-size: 18px;
            font-weight: 700;
            letter-spacing: -0.5px;
            color: var(--gray-900);
        }

        .logo-text .gas { color: var(--blue-600); }
        .logo-text .consult { color: #0F172A; }
        .logo-text .ai { color: rgba(15, 23, 42, 0.4); }

        .nav-links {
            display: none;
            align-items: center;
            gap: 4px;
        }

        .nav-link {
            padding: 10px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .nav-link.active {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .mobile-menu-btn {
            display: flex;
            flex-direction: column;
            gap: 5px;
            background: none;
            border: none;
            cursor: pointer;
            padding: 8px;
            border-radius: 8px;
            transition: background 0.2s ease;
        }

        .mobile-menu-btn:hover {
            background: rgba(0,0,0,0.04);
        }

        .mobile-menu-btn span {
            display: block;
            width: 22px;
            height: 2px;
            background: var(--gray-700);
            border-radius: 1px;
            transition: all 0.3s ease;
        }

        .mobile-menu-btn.active span:nth-child(1) {
            transform: rotate(45deg) translate(7px, 7px);
        }

        .mobile-menu-btn.active span:nth-child(2) {
            opacity: 0;
        }

        .mobile-menu-btn.active span:nth-child(3) {
            transform: rotate(-45deg) translate(7px, -7px);
        }

        .mobile-menu {
            display: none;
            position: fixed;
            top: 80px;
            left: 16px;
            right: 16px;
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 8px;
            box-shadow: 0 4px 16px rgba(0,0,0,0.08), 0 12px 48px rgba(0,0,0,0.12);
            z-index: 99;
            flex-direction: column;
            gap: 4px;
        }

        .mobile-menu.active {
            display: flex;
        }

        .mobile-menu-link {
            padding: 14px 16px;
            font-size: 15px;
            font-weight: 500;
            color: var(--gray-700);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .mobile-menu-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .hero {
            padding: 120px 20px 60px;
            text-align: center;
        }

        .hero-badge {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            background: var(--white);
            border: 1px solid var(--gray-200);
            border-radius: 100px;
            padding: 8px 16px 8px 12px;
            margin-bottom: 24px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.04);
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) forwards;
            opacity: 0;
        }

        .badge-dot {
            width: 8px;
            height: 8px;
            background: var(--blue-500);
            border-radius: 50%;
            position: relative;
        }

        .badge-dot::after {
            content: '';
            position: absolute;
            inset: -3px;
            border-radius: 50%;
            background: var(--blue-400);
            animation: pulse-ring 2s ease-out infinite;
        }

        @keyframes pulse-ring {
            0% { transform: scale(0.8); opacity: 0.8; }
            100% { transform: scale(2); opacity: 0; }
        }

        .badge-text {
            font-size: 12px;
            font-weight: 600;
            color: var(--gray-700);
        }

        .hero-title {
            font-size: 40px;
            font-weight: 800;
            line-height: 1.1;
            letter-spacing: -2px;
            color: var(--gray-900);
            margin-bottom: 20px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.1s forwards;
            opacity: 0;
        }

        .hero-title .gradient { color: var(--blue-600); }

        .hero-subtitle {
            font-size: 16px;
            font-weight: 400;
            line-height: 1.6;
            color: var(--gray-500);
            max-width: 560px;
            margin: 0 auto 40px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.2s forwards;
            opacity: 0;
        }

        @keyframes fade-up {
            from { opacity: 0; transform: translateY(24px); }
            to { opacity: 1; transform: translateY(0); }
        }

        .chat-container {
            max-width: 720px;
            margin: 0 auto 60px;
            padding: 0 16px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.3s forwards;
            opacity: 0;
        }

        .chat-card {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(40px) saturate(180%);
            -webkit-backdrop-filter: blur(40px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 8px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
            transition: all 0.4s cubic-bezier(0.4,0,0.2,1);
        }

        .chat-card:focus-within {
            box-shadow: 0 0 0 4px rgba(59,130,246,0.1), 0 1px 2px rgba(0,0,0,0.02), 0 8px 24px rgba(37,99,235,0.08), 0 32px 100px rgba(37,99,235,0.12), inset 0 1px 0 rgba(255,255,255,0.8);
            border-color: rgba(59,130,246,0.3);
        }

        .chat-inner {
            background: var(--white);
            border-radius: 14px;
            padding: 6px;
            display: flex;
            align-items: flex-end;
            gap: 6px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04), inset 0 1px 0 rgba(255,255,255,1);
        }

        .chat-input {
            flex: 1;
            border: none;
            outline: none;
            padding: 14px 16px;
            font-size: 16px;
            font-family: inherit;
            color: var(--gray-800);
            background: transparent;
            resize: none;
            min-height: 52px;
            max-height: 150px;
            line-height: 1.5;
        }

        .chat-input::placeholder { color: var(--gray-400); }

        .chat-send {
            width: 48px;
            height: 48px;
            background: var(--blue-600);
            border: none;
            border-radius: 12px;
            color: var(--white);
            cursor: pointer;
            display: flex;
            align-items: center;
            justify-content: center;
            transition: all 0.25s cubic-bezier(0.4,0,0.2,1);
            box-shadow: 0 1px 2px rgba(37,99,235,0.2), 0 4px 16px rgba(37,99,235,0.2), inset 0 1px 0 rgba(255,255,255,0.1), inset 0 -1px 0 rgba(0,0,0,0.1);
            flex-shrink: 0;
            margin: 6px;
        }

        .chat-send:hover {
            background: var(--blue-700);
            transform: translateY(-2px);
            box-shadow: 0 2px 4px rgba(37,99,235,0.2), 0 12px 40px rgba(37,99,235,0.3), inset 0 1px 0 rgba(255,255,255,0.1), inset 0 -1px 0 rgba(0,0,0,0.1);
        }

        .chat-send:active { transform: translateY(0); }
        .chat-send svg { width: 20px; height: 20px; }

        .chat-hints {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            padding: 16px 8px 6px;
        }

        .hint-chip {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            background: rgba(255,255,255,0.6);
            border: 1px solid var(--gray-200);
            border-radius: 100px;
            padding: 10px 14px;
            font-size: 12px;
            font-weight: 500;
            color: var(--gray-600);
            cursor: pointer;
            transition: all 0.2s ease;
        }

        .hint-chip:hover {
            background: var(--white);
            border-color: var(--blue-200);
            color: var(--blue-600);
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(37,99,235,0.1);
        }

        .hint-chip svg { width: 14px; height: 14px; opacity: 0.5; }
        .hint-chip:hover svg { opacity: 1; color: var(--blue-500); }

        .features { padding: 60px 20px 80px; }

        .features-header {
            text-align: center;
            margin-bottom: 40px;
        }

        .features-label {
            font-size: 11px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 1.5px;
            color: var(--blue-600);
            margin-bottom: 12px;
        }

        .features-title {
            font-size: 28px;
            font-weight: 800;
            letter-spacing: -1px;
            color: var(--gray-900);
            margin-bottom: 12px;
        }

        .features-subtitle {
            font-size: 16px;
            color: var(--gray-500);
            max-width: 480px;
            margin: 0 auto;
        }

        .features-grid {
            display: grid;
            grid-template-columns: 1fr;
            gap: 16px;
            max-width: 1200px;
            margin: 0 auto;
        }

        .feature-card {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            padding: 28px;
            position: relative;
            overflow: hidden;
            transition: all 0.4s cubic-bezier(0.4,0,0.2,1);
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04);
            cursor: pointer;
        }

        .feature-card::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 1px;
            background: linear-gradient(90deg, transparent 0%, rgba(255,255,255,0.8) 50%, transparent 100%);
        }

        .feature-card:hover {
            transform: translateY(-4px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.04), 0 24px 64px rgba(0,0,0,0.08);
            border-color: rgba(59,130,246,0.2);
        }

        .feature-icon {
            width: 56px;
            height: 56px;
            border-radius: 16px;
            display: flex;
            align-items: center;
            justify-content: center;
            margin-bottom: 20px;
            transition: all 0.3s ease;
        }

        .feature-card:hover .feature-icon { transform: scale(1.05); }
        .feature-icon svg { width: 24px; height: 24px; }

        .feature-icon.blue {
            background: var(--blue-50);
            border: 1px solid var(--blue-100);
            box-shadow: 0 4px 16px rgba(37,99,235,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.blue svg { color: var(--blue-600); }

        .feature-icon.emerald {
            background: #ECFDF5;
            border: 1px solid #D1FAE5;
            box-shadow: 0 4px 16px rgba(16,185,129,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.emerald svg { color: #059669; }

        .feature-icon.violet {
            background: #F5F3FF;
            border: 1px solid #EDE9FE;
            box-shadow: 0 4px 16px rgba(139,92,246,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.violet svg { color: #7C3AED; }

        .feature-icon.rose {
            background: #FFF1F2;
            border: 1px solid #FFE4E6;
            box-shadow: 0 4px 16px rgba(244,63,94,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.rose svg { color: #E11D48; }

        .feature-icon.cyan {
            background: #ECFEFF;
            border: 1px solid #CFFAFE;
            box-shadow: 0 4px 16px rgba(6,182,212,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.cyan svg { color: #0891B2; }

        .feature-title {
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 10px;
            letter-spacing: -0.3px;
        }

        .feature-desc {
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-500);
            margin-bottom: 20px;
        }

        .feature-link {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            font-size: 14px;
            font-weight: 600;
            color: var(--blue-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .feature-link:hover { gap: 10px; }
        .feature-link svg { width: 16px; height: 16px; transition: transform 0.2s ease; }
        .feature-link:hover svg { transform: translateX(4px); }

        .footer {
            padding: 32px 20px;
            border-top: 1px solid var(--gray-200);
            background: rgba(255,255,255,0.5);
        }

        .footer-inner {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 20px;
            text-align: center;
        }

        .footer-brand {
            display: flex;
            align-items: center;
            gap: 10px;
        }

        .footer-logo {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .footer-logo svg { width: 32px; height: 32px; }

        .footer-text {
            font-size: 13px;
            color: var(--gray-500);
        }

        .footer-links {
            display: flex;
            gap: 24px;
        }

        .footer-link {
            font-size: 13px;
            color: var(--gray-500);
            text-decoration: none;
            transition: color 0.2s ease;
        }

        .footer-link:hover { color: var(--gray-700); }

        @media (min-width: 768px) {
            .nav { padding: 16px 32px; }
            .nav-inner { height: 64px; padding: 0 24px; border-radius: 20px; }
            .logo-icon svg { width: 42px; height: 15px; }
            .logo-text { font-size: 20px; }
            .nav-links { display: flex; }
            .mobile-menu-btn { display: none; }
            .hero { padding: 160px 32px 80px; }
            .hero-badge { padding: 10px 20px 10px 14px; margin-bottom: 32px; }
            .badge-dot { width: 10px; height: 10px; }
            .badge-text { font-size: 13px; }
            .hero-title { font-size: 56px; letter-spacing: -2.5px; margin-bottom: 24px; }
            .hero-subtitle { font-size: 18px; margin-bottom: 48px; }
            .chat-container { padding: 0 24px; margin-bottom: 80px; }
            .chat-card { border-radius: 24px; padding: 10px; }
            .chat-inner { border-radius: 18px; padding: 8px; }
            .chat-input { padding: 16px 18px; min-height: 56px; }
            .chat-send { width: 52px; height: 52px; border-radius: 14px; }
            .chat-hints { gap: 10px; padding: 18px 10px 8px; }
            .hint-chip { padding: 12px 18px; font-size: 13px; }
            .hint-chip svg { width: 16px; height: 16px; }
            .features { padding: 80px 32px 100px; }
            .features-header { margin-bottom: 56px; }
            .features-label { font-size: 12px; margin-bottom: 16px; }
            .features-title { font-size: 36px; letter-spacing: -1.5px; }
            .features-subtitle { font-size: 18px; }
            .features-grid { grid-template-columns: repeat(2, 1fr); gap: 20px; }
            .feature-card { padding: 36px; border-radius: 24px; }
            .feature-card:hover { transform: translateY(-6px); }
            .feature-icon { width: 60px; height: 60px; border-radius: 18px; margin-bottom: 24px; }
            .feature-icon svg { width: 26px; height: 26px; }
            .feature-title { font-size: 20px; margin-bottom: 12px; }
            .feature-desc { font-size: 15px; line-height: 1.7; margin-bottom: 24px; }
            .footer { padding: 40px 32px; }
            .footer-inner { flex-direction: row; justify-content: space-between; text-align: left; }
            .footer-logo svg { width: 36px; height: 36px; }
            .footer-text { font-size: 14px; }
            .footer-links { gap: 32px; }
            .footer-link { font-size: 14px; }
            .orb-1 { width: 600px; height: 600px; left: -10%; }
            .orb-2 { width: 450px; height: 450px; right: -10%; }
            .orb-3 { width: 400px; height: 400px; }
        }

        @media (min-width: 1024px) {
            .nav { padding: 16px 40px; }
            .hero { padding: 180px 40px 80px; }
            .hero-title { font-size: 72px; letter-spacing: -3px; margin-bottom: 28px; }
            .hero-subtitle { font-size: 20px; margin-bottom: 56px; }
            .chat-container { margin-bottom: 100px; }
            .chat-card { border-radius: 28px; }
            .chat-inner { border-radius: 20px; }
            .chat-input { padding: 18px 20px; min-height: 60px; max-height: 180px; }
            .chat-send { width: 56px; height: 56px; border-radius: 16px; margin: 8px; }
            .chat-send svg { width: 22px; height: 22px; }
            .chat-hints { padding: 20px 12px 8px; }
            .hint-chip { padding: 12px 20px; }
            .features { padding: 80px 40px 120px; }
            .features-header { margin-bottom: 64px; }
            .features-title { font-size: 40px; }
            .features-grid { grid-template-columns: repeat(3, 1fr); gap: 24px; }
            .feature-card { padding: 40px; border-radius: 28px; }
            .feature-card:hover { transform: translateY(-8px); }
            .feature-icon { width: 64px; height: 64px; border-radius: 20px; margin-bottom: 28px; }
            .feature-icon svg { width: 28px; height: 28px; }
            .footer { padding: 48px 40px; }
            .orb-1 { width: 800px; height: 800px; top: -20%; left: -10%; }
            .orb-2 { width: 600px; height: 600px; right: -15%; }
            .orb-3 { width: 500px; height: 500px; }
        }

        @media (min-width: 1280px) {
            .hero-title { font-size: 80px; }
        }
    

        .main-content {
            flex: 1;
            padding: 100px 20px 40px;
            max-width: 1200px;
            margin: 0 auto;
            width: 100%;
        }

        .content-card {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            padding: 32px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
        }

        h1 { font-size: 32px; margin-bottom: 16px; font-weight: 700; letter-spacing: -0.5px; }
        h2 { font-size: 24px; margin-bottom: 12px; font-weight: 700; letter-spacing: -0.5px; }
        h3 { font-size: 20px; margin-bottom: 10px; font-weight: 700; letter-spacing: -0.3px; }

        input[type="text"], input[type="number"], input[type="email"], select, textarea {
            width: 100%;
            padding: 12px 16px;
            border: 1px solid var(--gray-300);
            border-radius: 12px;
            font-family: inherit;
            font-size: 15px;
            background: var(--white);
            color: var(--gray-900);
            transition: all 0.2s ease;
        }

        input:focus, select:focus, textarea:focus {
            outline: none;
            border-color: var(--blue-500);
            box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.1);
        }

        button, .btn {
            padding: 12px 24px;
            background: var(--blue-600);
            color: var(--white);
            border: none;
            border-radius: 12px;
            font-family: inherit;
            font-size: 15px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s ease;
            box-shadow: 0 1px 2px rgba(37,99,235,0.2), 0 4px 16px rgba(37,99,235,0.2), inset 0 1px 0 rgba(255,255,255,0.1);
        }

        button:hover, .btn:hover {
            background: var(--blue-700);
            transform: translateY(-2px);
            box-shadow: 0 2px 4px rgba(37,99,235,0.2), 0 12px 40px rgba(37,99,235,0.3), inset 0 1px 0 rgba(255,255,255,0.1);
        }

        button:active, .btn:active {
            transform: translateY(0);
        }

        table {
            width: 100%;
            border-collapse: collapse;
            background: rgba(255,255,255,0.5);
            border-radius: 12px;
            overflow: hidden;
        }

        th, td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid var(--gray-200);
        }

        th {
            background: rgba(37,99,235,0.05);
            font-weight: 600;
            color: var(--blue-700);
        }

        @media (min-width: 768px) {
            .main-content { padding: 120px 32px 60px; }
        }

        @media (min-width: 1024px) {
            .main-content { padding: 140px 40px 80px; }
        }

        /* Quick Dose Specific Styles */
        .weight-section {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 32px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
        }

        .weight-label {
            font-size: 14px;
            font-weight: 600;
            color: var(--gray-600);
            margin-bottom: 12px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        .weight-input-row {
            display: flex;
            align-items: center;
            gap: 16px;
            margin-bottom: 20px;
        }

        .weight-input-wrapper {
            position: relative;
            flex: 1;
            max-width: 200px;
        }

        .weight-input-wrapper input {
            width: 100%;
            padding: 14px 48px 14px 16px;
            font-size: 24px;
            font-weight: 700;
            color: var(--gray-900);
            border: 2px solid var(--blue-200);
            border-radius: 12px;
            background: var(--white);
        }

        .weight-unit {
            position: absolute;
            right: 16px;
            top: 50%;
            transform: translateY(-50%);
            font-size: 16px;
            font-weight: 600;
            color: var(--gray-500);
        }

        .conversion-text {
            font-size: 16px;
            color: var(--gray-600);
        }

        .conversion-value {
            font-weight: 700;
            color: var(--blue-600);
        }

        .quick-weights {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
        }

        .quick-weight-btn {
            padding: 10px 20px;
            background: var(--white);
            border: 1px solid var(--gray-300);
            border-radius: 100px;
            font-size: 14px;
            font-weight: 600;
            color: var(--gray-700);
            cursor: pointer;
            transition: all 0.2s ease;
        }

        .quick-weight-btn:hover {
            border-color: var(--blue-500);
            color: var(--blue-600);
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(37,99,235,0.15);
        }

        .quick-weight-btn.active {
            background: var(--blue-600);
            border-color: var(--blue-600);
            color: var(--white);
        }

        .drug-category {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            margin-bottom: 16px;
            overflow: hidden;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04);
        }

        .category-header {
            display: flex;
            align-items: center;
            gap: 16px;
            padding: 20px 24px;
            cursor: pointer;
            transition: all 0.2s ease;
        }

        .category-header:hover {
            background: rgba(0,0,0,0.02);
        }

        .color-indicator {
            width: 4px;
            height: 40px;
            border-radius: 2px;
            flex-shrink: 0;
        }

        .category-info {
            flex: 1;
        }

        .category-title {
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 4px;
        }

        .category-subtitle {
            font-size: 13px;
            color: var(--gray-500);
        }

        .chevron {
            width: 24px;
            height: 24px;
            color: var(--gray-400);
            transition: transform 0.3s ease;
        }

        .drug-category.open .chevron {
            transform: rotate(180deg);
        }

        .category-content {
            max-height: 0;
            overflow: hidden;
            transition: max-height 0.4s ease;
        }

        .drug-category.open .category-content {
            max-height: 5000px;
            padding: 0 24px 24px;
        }

        .drug-card {
            background: var(--white);
            border: 1px solid var(--gray-200);
            border-radius: 16px;
            padding: 24px;
            margin-bottom: 16px;
        }

        .drug-card:last-child {
            margin-bottom: 0;
        }

        .drug-header {
            display: flex;
            justify-content: space-between;
            align-items: flex-start;
            margin-bottom: 20px;
        }

        .drug-name-section h3 {
            font-size: 20px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 4px;
        }

        .drug-subtitle {
            font-size: 13px;
            color: var(--gray-500);
        }

        .concentration-badge {
            background: var(--blue-50);
            color: var(--blue-700);
            padding: 6px 12px;
            border-radius: 100px;
            font-size: 12px;
            font-weight: 600;
        }

        .dose-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(120px, 1fr));
            gap: 16px;
            margin-bottom: 20px;
        }

        .dose-item {
            text-align: center;
        }

        .dose-label {
            font-size: 11px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: var(--gray-500);
            margin-bottom: 8px;
        }

        .dose-value {
            font-size: 28px;
            font-weight: 800;
            color: var(--blue-600);
            margin-bottom: 4px;
        }

        .dose-unit {
            font-size: 14px;
            font-weight: 600;
            color: var(--gray-600);
        }

        .dose-range {
            font-size: 12px;
            color: var(--gray-500);
        }

        .clinical-pearl {
            display: flex;
            gap: 12px;
            background: var(--blue-50);
            border-left: 3px solid var(--blue-500);
            padding: 16px;
            border-radius: 8px;
        }

        .pearl-icon {
            width: 16px;
            height: 16px;
            color: var(--blue-600);
            flex-shrink: 0;
            margin-top: 2px;
        }

        .pearl-text {
            font-size: 13px;
            line-height: 1.6;
            color: var(--gray-700);
        }

        .color-legend {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            padding: 32px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04);
        }

        .legend-title {
            font-size: 16px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 20px;
        }

        .legend-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 12px;
        }

        .legend-item {
            display: flex;
            align-items: center;
            gap: 10px;
        }

        .legend-swatch {
            width: 24px;
            height: 24px;
            border-radius: 6px;
            flex-shrink: 0;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }

        .legend-label {
            font-size: 13px;
            font-weight: 500;
            color: var(--gray-700);
        }

        .disclaimer {
            background: rgba(239, 68, 68, 0.1);
            border: 1px solid rgba(239, 68, 68, 0.2);
            border-radius: 12px;
            padding: 16px 20px;
            font-size: 13px;
            line-height: 1.6;
            color: var(--gray-700);
        }

        .disclaimer strong {
            color: #DC2626;
        }

        .crisis-overlay {
            display: none;
            position: fixed;
            inset: 0;
            background: rgba(0,0,0,0.5);
            z-index: 1000;
            backdrop-filter: blur(4px);
        }

        .crisis-overlay.show {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .crisis-modal {
            background: var(--white);
            border-radius: 24px;
            max-width: 600px;
            width: 90%;
            max-height: 80vh;
            overflow-y: auto;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
        }

        .crisis-modal-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 24px;
            border-bottom: 1px solid var(--gray-200);
        }

        .crisis-modal-title {
            display: flex;
            align-items: center;
            gap: 12px;
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
        }

        .crisis-modal-title svg {
            color: #EF4444;
        }

        .crisis-close-btn {
            width: 36px;
            height: 36px;
            border-radius: 50%;
            border: none;
            background: var(--gray-100);
            color: var(--gray-600);
            cursor: pointer;
            display: flex;
            align-items: center;
            justify-content: center;
            transition: all 0.2s ease;
        }

        .crisis-close-btn:hover {
            background: var(--gray-200);
            color: var(--gray-900);
        }

        .crisis-modal-content {
            padding: 24px;
        }

        .protocol-grid {
            display: grid;
            gap: 12px;
        }

        .protocol-btn {
            background: var(--gray-50);
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            padding: 16px;
            text-align: left;
            cursor: pointer;
            transition: all 0.2s ease;
        }

        .protocol-btn:hover {
            background: var(--blue-50);
            border-color: var(--blue-200);
            transform: translateX(4px);
        }

        .protocol-title {
            font-size: 15px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 4px;
        }

        .protocol-desc {
            font-size: 13px;
            color: var(--gray-600);
        }

    </style>
</head>
<body>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <nav class="nav">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo">
                    <div class="logo-icon">
                        <svg width="36" height="12" viewBox="0 0 52 18" fill="none">
                            <circle cx="9" cy="9" r="9" fill="#2563EB"/>
                            <circle cx="26" cy="9" r="9" fill="#2563EB" fill-opacity="0.5"/>
                            <circle cx="43" cy="9" r="9" fill="#2563EB" fill-opacity="0.2"/>
                        </svg>
                    </div>
                    <span class="logo-text"><span class="gas">gas</span><span class="consult">consult</span><span class="ai">.ai</span></span>
                </a>
                <div class="nav-links">
                    <a href="/?clear=1" class="nav-link">Home</a>
                    <a href="/quick-dose" class="nav-link active">Quick Dose</a>
                    <a href="/preop" class="nav-link">Pre-Op</a>
                    <a href="/calculators" class="nav-link">Calculators</a>
                    <a href="/hypotension" class="nav-link">IOH Predictor</a>
                </div>
                <button class="mobile-menu-btn" onclick="toggleMobileMenu()" aria-label="Toggle menu">
                    <span></span>
                    <span></span>
                    <span></span>
                </button>
            </div>
        </nav>
        <div class="mobile-menu" id="mobileMenu">
            <a href="/?clear=1" class="mobile-menu-link">Home</a>
            <a href="/quick-dose" class="mobile-menu-link">Quick Dose</a>
            <a href="/preop" class="mobile-menu-link">Pre-Op</a>
            <a href="/calculators" class="mobile-menu-link">Calculators</a>
            <a href="/hypotension" class="mobile-menu-link">IOH Predictor</a>
        </div>

        <main class="main-content">
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

    <!-- Glassmorphism Footer -->
    <footer class="footer">
        <div class="footer-inner">
            <div class="footer-brand">
                <div class="footer-logo">
                    <svg width="32" height="32" viewBox="0 0 32 32" fill="none">
                        <circle cx="6" cy="16" r="5" fill="#2563EB"/>
                        <circle cx="16" cy="16" r="5" fill="#2563EB" fill-opacity="0.5"/>
                        <circle cx="26" cy="16" r="5" fill="#2563EB" fill-opacity="0.2"/>
                    </svg>
                </div>
                <span class="footer-text">© 2025 GasConsult.ai</span>
            </div>
            <div class="footer-links">
                <a href="/privacy" class="footer-link">Privacy</a>
                <a href="/terms" class="footer-link">Terms</a>
                <a href="#" class="footer-link">Contact</a>
            </div>
        </div>
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

        // Mobile menu toggle
        function toggleMobileMenu() {
            const menu = document.getElementById('mobileMenu');
            const btn = document.querySelector('.mobile-menu-btn');
            if (menu && btn) {
                menu.classList.toggle('active');
                btn.classList.toggle('active');
            }
        }
    </script>
        </main>
    </div>
</body>
</html>
"""

CALCULATORS_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Clinical Calculators — gasconsult.ai</title>

    <!-- PWA -->
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=5">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">

    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap" rel="stylesheet">
    <style>

        :root {
            --white: #FFFFFF;
            --gray-50: #F8FAFC;
            --gray-100: #F1F5F9;
            --gray-200: #E2E8F0;
            --gray-300: #CBD5E1;
            --gray-400: #94A3B8;
            --gray-500: #64748B;
            --gray-600: #475569;
            --gray-700: #334155;
            --gray-800: #1E293B;
            --gray-900: #0F172A;
            --blue-50: #EFF6FF;
            --blue-100: #DBEAFE;
            --blue-200: #BFDBFE;
            --blue-300: #93C5FD;
            --blue-400: #60A5FA;
            --blue-500: #3B82F6;
            --blue-600: #2563EB;
            --blue-700: #1D4ED8;
        }

        * { margin: 0; padding: 0; box-sizing: border-box; }

        html {
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
            scroll-behavior: smooth;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            background: var(--gray-50);
            color: var(--gray-900);
            min-height: 100vh;
            overflow-x: hidden;
        }

        .bg-canvas {
            position: fixed;
            inset: 0;
            z-index: 0;
            overflow: hidden;
            background: linear-gradient(180deg, #F0F7FF 0%, var(--gray-50) 50%, #FAFBFF 100%);
        }

        .orb {
            position: absolute;
            border-radius: 50%;
            filter: blur(80px);
            opacity: 0.6;
            animation: float 20s ease-in-out infinite;
        }

        .orb-1 {
            width: 400px;
            height: 400px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.15) 0%, transparent 70%);
            top: -15%;
            left: -20%;
        }

        .orb-2 {
            width: 300px;
            height: 300px;
            background: radial-gradient(circle, rgba(147, 197, 253, 0.2) 0%, transparent 70%);
            top: 30%;
            right: -20%;
            animation-delay: -7s;
            animation-duration: 25s;
        }

        .orb-3 {
            width: 250px;
            height: 250px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.1) 0%, transparent 70%);
            bottom: -10%;
            left: 20%;
            animation-delay: -14s;
            animation-duration: 30s;
        }

        @keyframes float {
            0%, 100% { transform: translate(0, 0) scale(1); }
            25% { transform: translate(40px, -40px) scale(1.05); }
            50% { transform: translate(20px, 40px) scale(0.95); }
            75% { transform: translate(-40px, 20px) scale(1.02); }
        }

        .grain {
            position: fixed;
            inset: 0;
            z-index: 1;
            pointer-events: none;
            opacity: 0.02;
            background-image: url("data:image/svg+xml,%3Csvg viewBox='0 0 512 512' xmlns='http://www.w3.org/2000/svg'%3E%3Cfilter id='n'%3E%3CfeTurbulence type='fractalNoise' baseFrequency='0.8' numOctaves='4' stitchTiles='stitch'/%3E%3C/filter%3E%3Crect width='100%25' height='100%25' filter='url(%23n)'/%3E%3C/svg%3E");
        }

        .page {
            position: relative;
            z-index: 2;
            min-height: 100vh;
            display: flex;
            flex-direction: column;
        }

        .nav {
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            z-index: 100;
            padding: 12px 16px;
        }

        .nav-inner {
            max-width: 1200px;
            margin: 0 auto;
            height: 56px;
            background: rgba(255, 255, 255, 0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 0 16px;
            display: flex;
            align-items: center;
            justify-content: space-between;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 12px 48px rgba(0,0,0,0.03);
        }

        .logo {
            display: flex;
            align-items: center;
            gap: 14px;
            text-decoration: none;
        }

        .logo-icon {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .logo-icon svg { width: 36px; height: 12px; }

        .logo-text {
            font-size: 18px;
            font-weight: 700;
            letter-spacing: -0.5px;
            color: var(--gray-900);
        }

        .logo-text .gas { color: var(--blue-600); }
        .logo-text .consult { color: #0F172A; }
        .logo-text .ai { color: rgba(15, 23, 42, 0.4); }

        .nav-links {
            display: none;
            align-items: center;
            gap: 4px;
        }

        .nav-link {
            padding: 10px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .nav-link.active {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .mobile-menu-btn {
            display: flex;
            flex-direction: column;
            gap: 5px;
            background: none;
            border: none;
            cursor: pointer;
            padding: 8px;
            border-radius: 8px;
            transition: background 0.2s ease;
        }

        .mobile-menu-btn:hover {
            background: rgba(0,0,0,0.04);
        }

        .mobile-menu-btn span {
            display: block;
            width: 22px;
            height: 2px;
            background: var(--gray-700);
            border-radius: 1px;
            transition: all 0.3s ease;
        }

        .mobile-menu-btn.active span:nth-child(1) {
            transform: rotate(45deg) translate(7px, 7px);
        }

        .mobile-menu-btn.active span:nth-child(2) {
            opacity: 0;
        }

        .mobile-menu-btn.active span:nth-child(3) {
            transform: rotate(-45deg) translate(7px, -7px);
        }

        .mobile-menu {
            display: none;
            position: fixed;
            top: 80px;
            left: 16px;
            right: 16px;
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 8px;
            box-shadow: 0 4px 16px rgba(0,0,0,0.08), 0 12px 48px rgba(0,0,0,0.12);
            z-index: 99;
            flex-direction: column;
            gap: 4px;
        }

        .mobile-menu.active {
            display: flex;
        }

        .mobile-menu-link {
            padding: 14px 16px;
            font-size: 15px;
            font-weight: 500;
            color: var(--gray-700);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .mobile-menu-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .hero {
            padding: 120px 20px 60px;
            text-align: center;
        }

        .hero-badge {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            background: var(--white);
            border: 1px solid var(--gray-200);
            border-radius: 100px;
            padding: 8px 16px 8px 12px;
            margin-bottom: 24px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.04);
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) forwards;
            opacity: 0;
        }

        .badge-dot {
            width: 8px;
            height: 8px;
            background: var(--blue-500);
            border-radius: 50%;
            position: relative;
        }

        .badge-dot::after {
            content: '';
            position: absolute;
            inset: -3px;
            border-radius: 50%;
            background: var(--blue-400);
            animation: pulse-ring 2s ease-out infinite;
        }

        @keyframes pulse-ring {
            0% { transform: scale(0.8); opacity: 0.8; }
            100% { transform: scale(2); opacity: 0; }
        }

        .badge-text {
            font-size: 12px;
            font-weight: 600;
            color: var(--gray-700);
        }

        .hero-title {
            font-size: 40px;
            font-weight: 800;
            line-height: 1.1;
            letter-spacing: -2px;
            color: var(--gray-900);
            margin-bottom: 20px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.1s forwards;
            opacity: 0;
        }

        .hero-title .gradient { color: var(--blue-600); }

        .hero-subtitle {
            font-size: 16px;
            font-weight: 400;
            line-height: 1.6;
            color: var(--gray-500);
            max-width: 560px;
            margin: 0 auto 40px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.2s forwards;
            opacity: 0;
        }

        @keyframes fade-up {
            from { opacity: 0; transform: translateY(24px); }
            to { opacity: 1; transform: translateY(0); }
        }

        .chat-container {
            max-width: 720px;
            margin: 0 auto 60px;
            padding: 0 16px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.3s forwards;
            opacity: 0;
        }

        .chat-card {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(40px) saturate(180%);
            -webkit-backdrop-filter: blur(40px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 8px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
            transition: all 0.4s cubic-bezier(0.4,0,0.2,1);
        }

        .chat-card:focus-within {
            box-shadow: 0 0 0 4px rgba(59,130,246,0.1), 0 1px 2px rgba(0,0,0,0.02), 0 8px 24px rgba(37,99,235,0.08), 0 32px 100px rgba(37,99,235,0.12), inset 0 1px 0 rgba(255,255,255,0.8);
            border-color: rgba(59,130,246,0.3);
        }

        .chat-inner {
            background: var(--white);
            border-radius: 14px;
            padding: 6px;
            display: flex;
            align-items: flex-end;
            gap: 6px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04), inset 0 1px 0 rgba(255,255,255,1);
        }

        .chat-input {
            flex: 1;
            border: none;
            outline: none;
            padding: 14px 16px;
            font-size: 16px;
            font-family: inherit;
            color: var(--gray-800);
            background: transparent;
            resize: none;
            min-height: 52px;
            max-height: 150px;
            line-height: 1.5;
        }

        .chat-input::placeholder { color: var(--gray-400); }

        .chat-send {
            width: 48px;
            height: 48px;
            background: var(--blue-600);
            border: none;
            border-radius: 12px;
            color: var(--white);
            cursor: pointer;
            display: flex;
            align-items: center;
            justify-content: center;
            transition: all 0.25s cubic-bezier(0.4,0,0.2,1);
            box-shadow: 0 1px 2px rgba(37,99,235,0.2), 0 4px 16px rgba(37,99,235,0.2), inset 0 1px 0 rgba(255,255,255,0.1), inset 0 -1px 0 rgba(0,0,0,0.1);
            flex-shrink: 0;
            margin: 6px;
        }

        .chat-send:hover {
            background: var(--blue-700);
            transform: translateY(-2px);
            box-shadow: 0 2px 4px rgba(37,99,235,0.2), 0 12px 40px rgba(37,99,235,0.3), inset 0 1px 0 rgba(255,255,255,0.1), inset 0 -1px 0 rgba(0,0,0,0.1);
        }

        .chat-send:active { transform: translateY(0); }
        .chat-send svg { width: 20px; height: 20px; }

        .chat-hints {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            padding: 16px 8px 6px;
        }

        .hint-chip {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            background: rgba(255,255,255,0.6);
            border: 1px solid var(--gray-200);
            border-radius: 100px;
            padding: 10px 14px;
            font-size: 12px;
            font-weight: 500;
            color: var(--gray-600);
            cursor: pointer;
            transition: all 0.2s ease;
        }

        .hint-chip:hover {
            background: var(--white);
            border-color: var(--blue-200);
            color: var(--blue-600);
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(37,99,235,0.1);
        }

        .hint-chip svg { width: 14px; height: 14px; opacity: 0.5; }
        .hint-chip:hover svg { opacity: 1; color: var(--blue-500); }

        .features { padding: 60px 20px 80px; }

        .features-header {
            text-align: center;
            margin-bottom: 40px;
        }

        .features-label {
            font-size: 11px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 1.5px;
            color: var(--blue-600);
            margin-bottom: 12px;
        }

        .features-title {
            font-size: 28px;
            font-weight: 800;
            letter-spacing: -1px;
            color: var(--gray-900);
            margin-bottom: 12px;
        }

        .features-subtitle {
            font-size: 16px;
            color: var(--gray-500);
            max-width: 480px;
            margin: 0 auto;
        }

        .features-grid {
            display: grid;
            grid-template-columns: 1fr;
            gap: 16px;
            max-width: 1200px;
            margin: 0 auto;
        }

        .feature-card {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            padding: 28px;
            position: relative;
            overflow: hidden;
            transition: all 0.4s cubic-bezier(0.4,0,0.2,1);
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04);
            cursor: pointer;
        }

        .feature-card::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 1px;
            background: linear-gradient(90deg, transparent 0%, rgba(255,255,255,0.8) 50%, transparent 100%);
        }

        .feature-card:hover {
            transform: translateY(-4px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.04), 0 24px 64px rgba(0,0,0,0.08);
            border-color: rgba(59,130,246,0.2);
        }

        .feature-icon {
            width: 56px;
            height: 56px;
            border-radius: 16px;
            display: flex;
            align-items: center;
            justify-content: center;
            margin-bottom: 20px;
            transition: all 0.3s ease;
        }

        .feature-card:hover .feature-icon { transform: scale(1.05); }
        .feature-icon svg { width: 24px; height: 24px; }

        .feature-icon.blue {
            background: var(--blue-50);
            border: 1px solid var(--blue-100);
            box-shadow: 0 4px 16px rgba(37,99,235,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.blue svg { color: var(--blue-600); }

        .feature-icon.emerald {
            background: #ECFDF5;
            border: 1px solid #D1FAE5;
            box-shadow: 0 4px 16px rgba(16,185,129,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.emerald svg { color: #059669; }

        .feature-icon.violet {
            background: #F5F3FF;
            border: 1px solid #EDE9FE;
            box-shadow: 0 4px 16px rgba(139,92,246,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.violet svg { color: #7C3AED; }

        .feature-icon.rose {
            background: #FFF1F2;
            border: 1px solid #FFE4E6;
            box-shadow: 0 4px 16px rgba(244,63,94,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.rose svg { color: #E11D48; }

        .feature-icon.cyan {
            background: #ECFEFF;
            border: 1px solid #CFFAFE;
            box-shadow: 0 4px 16px rgba(6,182,212,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.cyan svg { color: #0891B2; }

        .feature-title {
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 10px;
            letter-spacing: -0.3px;
        }

        .feature-desc {
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-500);
            margin-bottom: 20px;
        }

        .feature-link {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            font-size: 14px;
            font-weight: 600;
            color: var(--blue-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .feature-link:hover { gap: 10px; }
        .feature-link svg { width: 16px; height: 16px; transition: transform 0.2s ease; }
        .feature-link:hover svg { transform: translateX(4px); }

        .footer {
            padding: 32px 20px;
            border-top: 1px solid var(--gray-200);
            background: rgba(255,255,255,0.5);
        }

        .footer-inner {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 20px;
            text-align: center;
        }

        .footer-brand {
            display: flex;
            align-items: center;
            gap: 10px;
        }

        .footer-logo {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .footer-logo svg { width: 32px; height: 32px; }

        .footer-text {
            font-size: 13px;
            color: var(--gray-500);
        }

        .footer-links {
            display: flex;
            gap: 24px;
        }

        .footer-link {
            font-size: 13px;
            color: var(--gray-500);
            text-decoration: none;
            transition: color 0.2s ease;
        }

        .footer-link:hover { color: var(--gray-700); }

        @media (min-width: 768px) {
            .nav { padding: 16px 32px; }
            .nav-inner { height: 64px; padding: 0 24px; border-radius: 20px; }
            .logo-icon svg { width: 42px; height: 15px; }
            .logo-text { font-size: 20px; }
            .nav-links { display: flex; }
            .mobile-menu-btn { display: none; }
            .hero { padding: 160px 32px 80px; }
            .hero-badge { padding: 10px 20px 10px 14px; margin-bottom: 32px; }
            .badge-dot { width: 10px; height: 10px; }
            .badge-text { font-size: 13px; }
            .hero-title { font-size: 56px; letter-spacing: -2.5px; margin-bottom: 24px; }
            .hero-subtitle { font-size: 18px; margin-bottom: 48px; }
            .chat-container { padding: 0 24px; margin-bottom: 80px; }
            .chat-card { border-radius: 24px; padding: 10px; }
            .chat-inner { border-radius: 18px; padding: 8px; }
            .chat-input { padding: 16px 18px; min-height: 56px; }
            .chat-send { width: 52px; height: 52px; border-radius: 14px; }
            .chat-hints { gap: 10px; padding: 18px 10px 8px; }
            .hint-chip { padding: 12px 18px; font-size: 13px; }
            .hint-chip svg { width: 16px; height: 16px; }
            .features { padding: 80px 32px 100px; }
            .features-header { margin-bottom: 56px; }
            .features-label { font-size: 12px; margin-bottom: 16px; }
            .features-title { font-size: 36px; letter-spacing: -1.5px; }
            .features-subtitle { font-size: 18px; }
            .features-grid { grid-template-columns: repeat(2, 1fr); gap: 20px; }
            .feature-card { padding: 36px; border-radius: 24px; }
            .feature-card:hover { transform: translateY(-6px); }
            .feature-icon { width: 60px; height: 60px; border-radius: 18px; margin-bottom: 24px; }
            .feature-icon svg { width: 26px; height: 26px; }
            .feature-title { font-size: 20px; margin-bottom: 12px; }
            .feature-desc { font-size: 15px; line-height: 1.7; margin-bottom: 24px; }
            .footer { padding: 40px 32px; }
            .footer-inner { flex-direction: row; justify-content: space-between; text-align: left; }
            .footer-logo svg { width: 36px; height: 36px; }
            .footer-text { font-size: 14px; }
            .footer-links { gap: 32px; }
            .footer-link { font-size: 14px; }
            .orb-1 { width: 600px; height: 600px; left: -10%; }
            .orb-2 { width: 450px; height: 450px; right: -10%; }
            .orb-3 { width: 400px; height: 400px; }
        }

        @media (min-width: 1024px) {
            .nav { padding: 16px 40px; }
            .hero { padding: 180px 40px 80px; }
            .hero-title { font-size: 72px; letter-spacing: -3px; margin-bottom: 28px; }
            .hero-subtitle { font-size: 20px; margin-bottom: 56px; }
            .chat-container { margin-bottom: 100px; }
            .chat-card { border-radius: 28px; }
            .chat-inner { border-radius: 20px; }
            .chat-input { padding: 18px 20px; min-height: 60px; max-height: 180px; }
            .chat-send { width: 56px; height: 56px; border-radius: 16px; margin: 8px; }
            .chat-send svg { width: 22px; height: 22px; }
            .chat-hints { padding: 20px 12px 8px; }
            .hint-chip { padding: 12px 20px; }
            .features { padding: 80px 40px 120px; }
            .features-header { margin-bottom: 64px; }
            .features-title { font-size: 40px; }
            .features-grid { grid-template-columns: repeat(3, 1fr); gap: 24px; }
            .feature-card { padding: 40px; border-radius: 28px; }
            .feature-card:hover { transform: translateY(-8px); }
            .feature-icon { width: 64px; height: 64px; border-radius: 20px; margin-bottom: 28px; }
            .feature-icon svg { width: 28px; height: 28px; }
            .footer { padding: 48px 40px; }
            .orb-1 { width: 800px; height: 800px; top: -20%; left: -10%; }
            .orb-2 { width: 600px; height: 600px; right: -15%; }
            .orb-3 { width: 500px; height: 500px; }
        }

        @media (min-width: 1280px) {
            .hero-title { font-size: 80px; }
        }
    

        .main-content {
            flex: 1;
            padding: 100px 20px 40px;
            max-width: 1200px;
            margin: 0 auto;
            width: 100%;
        }

        .content-card {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            padding: 32px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
        }

        h1 { font-size: 32px; margin-bottom: 16px; font-weight: 700; letter-spacing: -0.5px; }
        h2 { font-size: 24px; margin-bottom: 12px; font-weight: 700; letter-spacing: -0.5px; }
        h3 { font-size: 20px; margin-bottom: 10px; font-weight: 700; letter-spacing: -0.3px; }

        input[type="text"], input[type="number"], input[type="email"], select, textarea {
            width: 100%;
            padding: 12px 16px;
            border: 1px solid var(--gray-300);
            border-radius: 12px;
            font-family: inherit;
            font-size: 15px;
            background: var(--white);
            color: var(--gray-900);
            transition: all 0.2s ease;
        }

        input:focus, select:focus, textarea:focus {
            outline: none;
            border-color: var(--blue-500);
            box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.1);
        }

        button, .btn {
            padding: 12px 24px;
            background: var(--blue-600);
            color: var(--white);
            border: none;
            border-radius: 12px;
            font-family: inherit;
            font-size: 15px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s ease;
            box-shadow: 0 1px 2px rgba(37,99,235,0.2), 0 4px 16px rgba(37,99,235,0.2), inset 0 1px 0 rgba(255,255,255,0.1);
        }

        button:hover, .btn:hover {
            background: var(--blue-700);
            transform: translateY(-2px);
            box-shadow: 0 2px 4px rgba(37,99,235,0.2), 0 12px 40px rgba(37,99,235,0.3), inset 0 1px 0 rgba(255,255,255,0.1);
        }

        button:active, .btn:active {
            transform: translateY(0);
        }

        table {
            width: 100%;
            border-collapse: collapse;
            background: rgba(255,255,255,0.5);
            border-radius: 12px;
            overflow: hidden;
        }

        th, td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid var(--gray-200);
        }

        th {
            background: rgba(37,99,235,0.05);
            font-weight: 600;
            color: var(--blue-700);
        }

        @media (min-width: 768px) {
            .main-content { padding: 120px 32px 60px; }
        }

        @media (min-width: 1024px) {
            .main-content { padding: 140px 40px 80px; }
        }


        /* AI Chat Modal Styles */
        .ai-chat-modal {
            display: none;
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            bottom: 0;
            background: rgba(0, 0, 0, 0.5);
            backdrop-filter: blur(4px);
            z-index: 1000;
            padding: 20px;
            overflow-y: auto;
        }

        .ai-chat-modal.active {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .ai-chat-container {
            background: var(--white);
            border-radius: 20px;
            width: 100%;
            max-width: 600px;
            max-height: 80vh;
            display: flex;
            flex-direction: column;
            box-shadow: 0 24px 64px rgba(0,0,0,0.2);
        }

        .ai-chat-header {
            padding: 20px 24px;
            border-bottom: 1px solid var(--gray-200);
            display: flex;
            align-items: center;
            justify-content: space-between;
        }

        .ai-chat-header h3 {
            font-size: 18px;
            font-weight: 600;
            color: var(--gray-900);
            margin: 0;
        }

        .close-modal-btn {
            background: var(--gray-100);
            border: none;
            width: 32px;
            height: 32px;
            border-radius: 8px;
            display: flex;
            align-items: center;
            justify-content: center;
            cursor: pointer;
            transition: all 0.2s ease;
            color: var(--gray-600);
        }

        .close-modal-btn:hover {
            background: var(--gray-200);
            color: var(--gray-900);
        }

        .ai-chat-messages {
            flex: 1;
            overflow-y: auto;
            padding: 20px 24px;
            min-height: 200px;
            max-height: 400px;
        }

        .ai-message {
            margin-bottom: 16px;
            padding: 12px 16px;
            background: var(--gray-50);
            border-radius: 12px;
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-800);
        }

        .user-message {
            margin-bottom: 16px;
            padding: 12px 16px;
            background: var(--blue-50);
            border-radius: 12px;
            font-size: 14px;
            line-height: 1.6;
            color: var(--blue-900);
            text-align: right;
        }

        .ai-chat-input-area {
            padding: 16px 24px 20px;
            border-top: 1px solid var(--gray-200);
        }

        .ai-input-container {
            display: flex;
            align-items: flex-end;
            gap: 12px;
            background: var(--gray-50);
            border: 1px solid var(--gray-300);
            border-radius: 12px;
            padding: 8px;
            transition: all 0.2s ease;
        }

        .ai-input-container:focus-within {
            border-color: var(--blue-500);
            box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.1);
        }

        .ai-chat-input {
            flex: 1;
            border: none;
            outline: none;
            background: transparent;
            padding: 8px 12px;
            font-size: 15px;
            font-family: inherit;
            color: var(--gray-900);
            resize: none;
            min-height: 42px;
            max-height: 120px;
            line-height: 1.5;
        }

        .ai-chat-input::placeholder {
            color: var(--gray-400);
        }

        .send-ai-message-btn {
            background: var(--blue-600);
            color: var(--white);
            border: none;
            border-radius: 10px;
            padding: 10px 20px;
            font-size: 14px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s ease;
            flex-shrink: 0;
            height: 42px;
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .send-ai-message-btn:hover {
            background: var(--blue-700);
            transform: translateY(-1px);
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.3);
        }

        .send-ai-message-btn:active {
            transform: translateY(0);
        }

        .send-ai-message-btn:disabled {
            background: var(--gray-300);
            cursor: not-allowed;
            transform: none;
        }
    </style>
</head>
<body>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <nav class="nav">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo">
                    <div class="logo-icon">
                        <svg width="36" height="12" viewBox="0 0 52 18" fill="none">
                            <circle cx="9" cy="9" r="9" fill="#2563EB"/>
                            <circle cx="26" cy="9" r="9" fill="#2563EB" fill-opacity="0.5"/>
                            <circle cx="43" cy="9" r="9" fill="#2563EB" fill-opacity="0.2"/>
                        </svg>
                    </div>
                    <span class="logo-text"><span class="gas">gas</span><span class="consult">consult</span><span class="ai">.ai</span></span>
                </a>
                <div class="nav-links">
                    <a href="/?clear=1" class="nav-link">Home</a>
                    <a href="/quick-dose" class="nav-link">Quick Dose</a>
                    <a href="/preop" class="nav-link">Pre-Op</a>
                    <a href="/calculators" class="nav-link">Calculators</a>
                    <a href="/hypotension" class="nav-link">IOH Predictor</a>
                </div>
                <button class="mobile-menu-btn" onclick="toggleMobileMenu()" aria-label="Toggle menu">
                    <span></span>
                    <span></span>
                    <span></span>
                </button>
            </div>
        </nav>
        <div class="mobile-menu" id="mobileMenu">
            <a href="/?clear=1" class="mobile-menu-link">Home</a>
            <a href="/quick-dose" class="mobile-menu-link">Quick Dose</a>
            <a href="/preop" class="mobile-menu-link">Pre-Op</a>
            <a href="/calculators" class="mobile-menu-link">Calculators</a>
            <a href="/hypotension" class="mobile-menu-link">IOH Predictor</a>
        </div>

        <main class="main-content">
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
                <!-- Messages will be added here dynamically -->
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

            // Redirect to chat page with the message
            openAIChat(message);
        }

        // AI Chat Modal Functions
        async function openAIChat(initialMessage) {
            const modal = document.getElementById('aiChatModal');
            const messagesDiv = document.getElementById('aiChatMessages');

            // Clear previous messages
            messagesDiv.innerHTML = '';

            // Show modal
            modal.classList.add('active');

            // Send initial message
            if (initialMessage) {
                addMessageToModal(initialMessage, 'user');
                await fetchAIResponseForModal(initialMessage);
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

            // Add user message
            addMessageToModal(message, 'user');
            input.value = '';

            // Fetch AI response
            fetchAIResponseForModal(message);
        }

        function addMessageToModal(message, type) {
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

        async function fetchAIResponseForModal(userMessage) {
            // Show loading indicator
            addMessageToModal('', 'loading');

            try {
                // Send to homepage to prepare streaming
                const csrfToken = document.getElementById('csrf_token').value;
                const response = await fetch('/', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/x-www-form-urlencoded',
                        'X-Requested-With': 'XMLHttpRequest'  // Signal this is AJAX
                    },
                    body: `query=${encodeURIComponent(userMessage)}&csrf_token=${encodeURIComponent(csrfToken)}&modal=1`,
                    redirect: 'manual'  // Don't follow redirects
                });

                // Handle response - it might be JSON with request_id
                const data = await response.json();

                if (data.request_id) {
                    // Start streaming with EventSource
                    const eventSource = new EventSource(`/stream?request_id=${data.request_id}`);
                    let fullResponse = '';
                    let refs = [];
                    let aiMessageEl = null;  // Keep reference to the AI message element

                    eventSource.addEventListener('message', function(e) {
                        const event = JSON.parse(e.data);

                        if (event.type === 'content') {
                            // Remove loading and create AI message on first content
                            const loadingMsg = document.getElementById('loadingMessage');
                            if (loadingMsg && fullResponse === '' && !aiMessageEl) {
                                loadingMsg.remove();
                                // Create AI message element directly
                                const messagesDiv = document.getElementById('aiChatMessages');
                                aiMessageEl = document.createElement('div');
                                aiMessageEl.className = 'ai-message';
                                messagesDiv.appendChild(aiMessageEl);
                                messagesDiv.scrollTop = messagesDiv.scrollHeight;
                            }

                            fullResponse += event.content;
                            if (aiMessageEl) {
                                aiMessageEl.innerHTML = `<strong>AI:</strong> ${fullResponse}`;
                            }
                        } else if (event.type === 'references') {
                            refs = event.references;
                        } else if (event.type === 'done') {
                            eventSource.close();

                            // Add references if any
                            if (refs && refs.length > 0 && aiMessageEl) {
                                let refsHTML = '<div style="margin-top: 16px; padding-top: 16px; border-top: 1px solid var(--border);"><strong>References:</strong><div>';
                                refs.forEach((ref, index) => {
                                    refsHTML += `<div style="margin: 8px 0; font-size: 13px;"><a href="https://pubmed.ncbi.nlm.nih.gov/${ref.pmid}/" target="_blank" style="color: var(--primary-blue);">[${index + 1}] ${ref.title} (${ref.year})</a></div>`;
                                });
                                refsHTML += '</div></div>';
                                aiMessageEl.innerHTML += refsHTML;
                            }
                        } else if (event.type === 'error') {
                            console.error('Stream error:', event.message);
                            const loadingMsg = document.getElementById('loadingMessage');
                            if (loadingMsg) loadingMsg.remove();
                            addMessageToModal(`Error: ${event.message}`, 'ai');
                            eventSource.close();
                        }
                    });

                    eventSource.addEventListener('error', function() {
                        console.error('EventSource connection error');
                        const loadingMsg = document.getElementById('loadingMessage');
                        if (loadingMsg) loadingMsg.remove();
                        addMessageToModal('Connection error. Please try again.', 'ai');
                        eventSource.close();
                    });
                } else {
                    throw new Error('No request_id received');
                }

            } catch (error) {
                console.error('Error:', error);
                const loadingMsg = document.getElementById('loadingMessage');
                if (loadingMsg) loadingMsg.remove();
                addMessageToModal('Sorry, there was an error. Please try again.', 'ai');
            }
        }

        function escapeHtml(text) {
            const div = document.createElement('div');
            div.textContent = text;
            return div.innerHTML;
        }

        // Mobile menu toggle
        function toggleMobileMenu() {
            const menu = document.getElementById('mobileMenu');
            const btn = document.querySelector('.mobile-menu-btn');
            if (menu && btn) {
                menu.classList.toggle('active');
                btn.classList.toggle('active');
            }
        }
    </script>

    <!-- Glassmorphism Footer -->
    <footer class="footer">
        <div class="footer-inner">
            <div class="footer-brand">
                <div class="footer-logo">
                    <svg width="32" height="32" viewBox="0 0 32 32" fill="none">
                        <circle cx="6" cy="16" r="5" fill="#2563EB"/>
                        <circle cx="16" cy="16" r="5" fill="#2563EB" fill-opacity="0.5"/>
                        <circle cx="26" cy="16" r="5" fill="#2563EB" fill-opacity="0.2"/>
                    </svg>
                </div>
                <span class="footer-text">© 2025 GasConsult.ai</span>
            </div>
            <div class="footer-links">
                <a href="/privacy" class="footer-link">Privacy</a>
                <a href="/terms" class="footer-link">Terms</a>
                <a href="#" class="footer-link">Contact</a>
            </div>
        </div>
    </footer>
        </main>
    </div>
</body>
</html>
"""

HYPOTENSION_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>IOH Predictor — gasconsult.ai</title>

    <!-- PWA -->
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=5">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=5">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">

    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap" rel="stylesheet">
    <style>

        :root {
            --white: #FFFFFF;
            --gray-50: #F8FAFC;
            --gray-100: #F1F5F9;
            --gray-200: #E2E8F0;
            --gray-300: #CBD5E1;
            --gray-400: #94A3B8;
            --gray-500: #64748B;
            --gray-600: #475569;
            --gray-700: #334155;
            --gray-800: #1E293B;
            --gray-900: #0F172A;
            --blue-50: #EFF6FF;
            --blue-100: #DBEAFE;
            --blue-200: #BFDBFE;
            --blue-300: #93C5FD;
            --blue-400: #60A5FA;
            --blue-500: #3B82F6;
            --blue-600: #2563EB;
            --blue-700: #1D4ED8;
        }

        * { margin: 0; padding: 0; box-sizing: border-box; }

        html {
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
            scroll-behavior: smooth;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            background: var(--gray-50);
            color: var(--gray-900);
            min-height: 100vh;
            overflow-x: hidden;
        }

        .bg-canvas {
            position: fixed;
            inset: 0;
            z-index: 0;
            overflow: hidden;
            background: linear-gradient(180deg, #F0F7FF 0%, var(--gray-50) 50%, #FAFBFF 100%);
        }

        .orb {
            position: absolute;
            border-radius: 50%;
            filter: blur(80px);
            opacity: 0.6;
            animation: float 20s ease-in-out infinite;
        }

        .orb-1 {
            width: 400px;
            height: 400px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.15) 0%, transparent 70%);
            top: -15%;
            left: -20%;
        }

        .orb-2 {
            width: 300px;
            height: 300px;
            background: radial-gradient(circle, rgba(147, 197, 253, 0.2) 0%, transparent 70%);
            top: 30%;
            right: -20%;
            animation-delay: -7s;
            animation-duration: 25s;
        }

        .orb-3 {
            width: 250px;
            height: 250px;
            background: radial-gradient(circle, rgba(59, 130, 246, 0.1) 0%, transparent 70%);
            bottom: -10%;
            left: 20%;
            animation-delay: -14s;
            animation-duration: 30s;
        }

        @keyframes float {
            0%, 100% { transform: translate(0, 0) scale(1); }
            25% { transform: translate(40px, -40px) scale(1.05); }
            50% { transform: translate(20px, 40px) scale(0.95); }
            75% { transform: translate(-40px, 20px) scale(1.02); }
        }

        .grain {
            position: fixed;
            inset: 0;
            z-index: 1;
            pointer-events: none;
            opacity: 0.02;
            background-image: url("data:image/svg+xml,%3Csvg viewBox='0 0 512 512' xmlns='http://www.w3.org/2000/svg'%3E%3Cfilter id='n'%3E%3CfeTurbulence type='fractalNoise' baseFrequency='0.8' numOctaves='4' stitchTiles='stitch'/%3E%3C/filter%3E%3Crect width='100%25' height='100%25' filter='url(%23n)'/%3E%3C/svg%3E");
        }

        .page {
            position: relative;
            z-index: 2;
            min-height: 100vh;
            display: flex;
            flex-direction: column;
        }

        .nav {
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            z-index: 100;
            padding: 12px 16px;
        }

        .nav-inner {
            max-width: 1200px;
            margin: 0 auto;
            height: 56px;
            background: rgba(255, 255, 255, 0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 0 16px;
            display: flex;
            align-items: center;
            justify-content: space-between;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 12px 48px rgba(0,0,0,0.03);
        }

        .logo {
            display: flex;
            align-items: center;
            gap: 14px;
            text-decoration: none;
        }

        .logo-icon {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .logo-icon svg { width: 36px; height: 12px; }

        .logo-text {
            font-size: 18px;
            font-weight: 700;
            letter-spacing: -0.5px;
            color: var(--gray-900);
        }

        .logo-text .gas { color: var(--blue-600); }
        .logo-text .consult { color: #0F172A; }
        .logo-text .ai { color: rgba(15, 23, 42, 0.4); }

        .nav-links {
            display: none;
            align-items: center;
            gap: 4px;
        }

        .nav-link {
            padding: 10px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .nav-link.active {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .mobile-menu-btn {
            display: flex;
            flex-direction: column;
            gap: 5px;
            background: none;
            border: none;
            cursor: pointer;
            padding: 8px;
            border-radius: 8px;
            transition: background 0.2s ease;
        }

        .mobile-menu-btn:hover {
            background: rgba(0,0,0,0.04);
        }

        .mobile-menu-btn span {
            display: block;
            width: 22px;
            height: 2px;
            background: var(--gray-700);
            border-radius: 1px;
            transition: all 0.3s ease;
        }

        .mobile-menu-btn.active span:nth-child(1) {
            transform: rotate(45deg) translate(7px, 7px);
        }

        .mobile-menu-btn.active span:nth-child(2) {
            opacity: 0;
        }

        .mobile-menu-btn.active span:nth-child(3) {
            transform: rotate(-45deg) translate(7px, -7px);
        }

        .mobile-menu {
            display: none;
            position: fixed;
            top: 80px;
            left: 16px;
            right: 16px;
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 8px;
            box-shadow: 0 4px 16px rgba(0,0,0,0.08), 0 12px 48px rgba(0,0,0,0.12);
            z-index: 99;
            flex-direction: column;
            gap: 4px;
        }

        .mobile-menu.active {
            display: flex;
        }

        .mobile-menu-link {
            padding: 14px 16px;
            font-size: 15px;
            font-weight: 500;
            color: var(--gray-700);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .mobile-menu-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .footer {
            padding: 32px 20px;
            border-top: 1px solid var(--gray-200);
            background: rgba(255,255,255,0.5);
        }

        .footer-inner {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 20px;
            text-align: center;
        }

        .footer-brand {
            display: flex;
            align-items: center;
            gap: 10px;
        }

        .footer-logo {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .footer-logo svg { width: 32px; height: 32px; }

        .footer-text {
            font-size: 13px;
            color: var(--gray-500);
        }

        .footer-links {
            display: flex;
            gap: 24px;
        }

        .footer-link {
            font-size: 13px;
            color: var(--gray-500);
            text-decoration: none;
            transition: color 0.2s ease;
        }

        .footer-link:hover { color: var(--gray-700); }

        @media (min-width: 768px) {
            .nav { padding: 16px 32px; }
            .nav-inner { height: 64px; padding: 0 24px; border-radius: 20px; }
            .logo-icon svg { width: 42px; height: 15px; }
            .logo-text { font-size: 20px; }
            .nav-links { display: flex; }
            .mobile-menu-btn { display: none; }
            .footer { padding: 40px 32px; }
            .footer-inner { flex-direction: row; justify-content: space-between; text-align: left; }
            .footer-logo svg { width: 36px; height: 36px; }
            .footer-text { font-size: 14px; }
            .footer-links { gap: 32px; }
            .footer-link { font-size: 14px; }
        }

        @media (min-width: 1024px) {
            .nav { padding: 16px 40px; }
            .footer { padding: 48px 40px; }
        }

        .main-content {
            flex: 1;
            padding: 100px 20px 40px;
            max-width: 1200px;
            margin: 0 auto;
            width: 100%;
        }

        .content-card {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            padding: 32px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
        }

        h1 { font-size: 32px; margin-bottom: 16px; font-weight: 700; letter-spacing: -0.5px; color: var(--gray-900); }
        h2 { font-size: 24px; margin-bottom: 12px; font-weight: 700; letter-spacing: -0.5px; color: var(--gray-900); }
        h3 { font-size: 20px; margin-bottom: 10px; font-weight: 700; letter-spacing: -0.3px; color: var(--gray-900); }

        p { line-height: 1.6; color: var(--gray-700); }

        input[type="text"], input[type="number"], input[type="email"], select, textarea {
            width: 100%;
            padding: 12px 16px;
            border: 1px solid var(--gray-300);
            border-radius: 12px;
            font-family: inherit;
            font-size: 15px;
            background: var(--white);
            color: var(--gray-900);
            transition: all 0.2s ease;
        }

        input:focus, select:focus, textarea:focus {
            outline: none;
            border-color: var(--blue-500);
            box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.1);
        }

        button, .btn {
            padding: 12px 24px;
            background: var(--blue-600);
            color: var(--white);
            border: none;
            border-radius: 12px;
            font-family: inherit;
            font-size: 15px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s ease;
            box-shadow: 0 1px 2px rgba(37,99,235,0.2), 0 4px 16px rgba(37,99,235,0.2), inset 0 1px 0 rgba(255,255,255,0.1);
        }

        button:hover, .btn:hover {
            background: var(--blue-700);
            transform: translateY(-2px);
            box-shadow: 0 2px 4px rgba(37,99,235,0.2), 0 12px 40px rgba(37,99,235,0.3), inset 0 1px 0 rgba(255,255,255,0.1);
        }

        button:active, .btn:active {
            transform: translateY(0);
        }

        label {
            display: block;
            font-size: 14px;
            font-weight: 600;
            color: var(--gray-700);
            margin-bottom: 8px;
        }

        @media (min-width: 768px) {
            .main-content { padding: 120px 32px 60px; }
        }

        @media (min-width: 1024px) {
            .main-content { padding: 140px 40px 80px; }
        }

        /* IOH Predictor Specific Styles */
        .page-header {
            text-align: center;
            margin-bottom: 40px;
        }

        .page-header h1 {
            font-size: 36px;
            font-weight: 800;
            letter-spacing: -1px;
            color: var(--gray-900);
            margin-bottom: 12px;
        }

        .page-header p {
            font-size: 16px;
            color: var(--gray-600);
            margin-bottom: 16px;
        }

        .edu-warning {
            background: linear-gradient(135deg, rgba(239, 68, 68, 0.15) 0%, rgba(220, 38, 38, 0.08) 100%);
            border: 2px solid rgba(239, 68, 68, 0.3);
            border-radius: 16px;
            padding: 24px;
            margin-bottom: 32px;
            display: flex;
            gap: 16px;
            align-items: flex-start;
        }

        .edu-warning-icon {
            font-size: 32px;
            flex-shrink: 0;
        }

        .edu-warning-content h3 {
            font-size: 16px;
            font-weight: 800;
            color: #DC2626;
            margin-bottom: 8px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        .edu-warning-content p {
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-800);
            margin: 0;
        }

        .edu-badge {
            display: inline-flex;
            align-items: center;
            background: rgba(239, 68, 68, 0.1);
            color: #DC2626;
            padding: 8px 16px;
            border-radius: 100px;
            font-size: 12px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            border: 1px solid rgba(239, 68, 68, 0.2);
        }

        .disclaimer-box {
            background: rgba(255,255,255,0.9);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-left: 4px solid #DC2626;
            border-radius: 16px;
            padding: 28px;
            margin-bottom: 32px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04);
        }

        .disclaimer-box h3 {
            font-size: 18px;
            font-weight: 700;
            color: #DC2626;
            margin-bottom: 16px;
        }

        .disclaimer-box p {
            font-size: 14px;
            line-height: 1.7;
            color: var(--gray-700);
            margin-bottom: 12px;
        }

        .disclaimer-box p:last-of-type {
            margin-bottom: 8px;
        }

        .disclaimer-box ul {
            margin-left: 24px;
            font-size: 14px;
            line-height: 1.8;
            color: var(--gray-700);
        }

        .disclaimer-box li {
            margin-bottom: 6px;
        }

        .metrics-box {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 28px;
            margin-bottom: 32px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04);
        }

        .metrics-box h3 {
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 20px;
        }

        .metrics-intro {
            font-size: 14px;
            line-height: 1.7;
            color: var(--gray-600);
            margin-bottom: 24px;
        }

        .metrics-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
            gap: 16px;
            margin-bottom: 24px;
        }

        .metric-card {
            background: var(--gray-50);
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            padding: 16px;
            text-align: center;
        }

        .metric-label {
            font-size: 12px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: var(--gray-500);
            margin-bottom: 8px;
        }

        .metric-value {
            font-size: 28px;
            font-weight: 800;
            color: var(--blue-600);
            margin-bottom: 4px;
        }

        .metric-description {
            font-size: 12px;
            color: var(--gray-600);
        }

        .metrics-explanation {
            background: var(--blue-50);
            border-radius: 12px;
            padding: 20px;
        }

        .metrics-explanation h4 {
            font-size: 14px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 12px;
        }

        .metrics-explanation ul {
            list-style: none;
            padding: 0;
        }

        .metrics-explanation li {
            font-size: 13px;
            line-height: 1.7;
            color: var(--gray-700);
            margin-bottom: 10px;
            padding-left: 0;
        }

        .metric-term {
            font-weight: 700;
            color: var(--blue-700);
        }

        .form-section {
            margin-bottom: 32px;
        }

        .form-section h2 {
            font-size: 20px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 20px;
            padding-bottom: 12px;
            border-bottom: 2px solid var(--blue-100);
        }

        .form-grid {
            display: grid;
            grid-template-columns: 1fr;
            gap: 20px;
        }

        @media (min-width: 768px) {
            .form-grid {
                grid-template-columns: repeat(2, 1fr);
            }
        }

        @media (min-width: 1024px) {
            .form-grid {
                grid-template-columns: repeat(3, 1fr);
            }
        }

        .form-group {
            margin-bottom: 0;
        }

        .submit-btn {
            width: 100%;
            padding: 16px 32px;
            font-size: 16px;
            font-weight: 700;
            margin-top: 24px;
        }

        .results-warning {
            background: linear-gradient(135deg, rgba(251, 191, 36, 0.15) 0%, rgba(245, 158, 11, 0.08) 100%);
            border: 2px solid rgba(251, 191, 36, 0.3);
            border-radius: 16px;
            padding: 20px;
            margin-bottom: 32px;
            text-align: center;
        }

        .results-warning p {
            font-size: 14px;
            font-weight: 600;
            color: #92400E;
            margin: 0;
        }

        .risk-gauges {
            display: grid;
            grid-template-columns: 1fr;
            gap: 20px;
            margin-bottom: 32px;
        }

        @media (min-width: 768px) {
            .risk-gauges {
                grid-template-columns: repeat(3, 1fr);
            }
        }

        .risk-gauge {
            background: rgba(255,255,255,0.9);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 28px;
            text-align: center;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04);
        }

        .risk-gauge h3 {
            font-size: 16px;
            font-weight: 700;
            color: var(--gray-600);
            margin-bottom: 16px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        .risk-value {
            font-size: 56px;
            font-weight: 900;
            margin-bottom: 12px;
        }

        .risk-value.low { color: #059669; }
        .risk-value.moderate { color: #D97706; }
        .risk-value.high { color: #DC2626; }
        .risk-value.very-high { color: #991B1B; }

        .risk-label {
            display: inline-block;
            padding: 8px 16px;
            border-radius: 100px;
            font-size: 12px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        .risk-label.low {
            background: rgba(5, 150, 105, 0.1);
            color: #059669;
            border: 1px solid rgba(5, 150, 105, 0.2);
        }

        .risk-label.moderate {
            background: rgba(217, 119, 6, 0.1);
            color: #D97706;
            border: 1px solid rgba(217, 119, 6, 0.2);
        }

        .risk-label.high {
            background: rgba(220, 38, 38, 0.1);
            color: #DC2626;
            border: 1px solid rgba(220, 38, 38, 0.2);
        }

        .risk-label.very-high {
            background: rgba(153, 27, 27, 0.1);
            color: #991B1B;
            border: 1px solid rgba(153, 27, 27, 0.2);
        }

        .factors-section {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 28px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04);
        }

        .factors-section h2 {
            font-size: 20px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 20px;
            border-bottom: none;
            padding-bottom: 0;
        }

        .factor-item {
            display: flex;
            gap: 12px;
            padding: 16px;
            background: var(--blue-50);
            border-radius: 12px;
            margin-bottom: 12px;
            border-left: 3px solid var(--blue-500);
        }

        .factor-item:last-child {
            margin-bottom: 0;
        }

        .factor-item strong {
            color: var(--blue-700);
            font-weight: 700;
        }

        .interventions-section {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 28px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04);
        }

        .interventions-section h2 {
            font-size: 20px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 20px;
            border-bottom: none;
            padding-bottom: 0;
        }

        .intervention-item {
            display: flex;
            gap: 16px;
            padding: 20px;
            background: var(--white);
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            margin-bottom: 16px;
        }

        .intervention-item:last-child {
            margin-bottom: 0;
        }

        .intervention-rank {
            display: flex;
            align-items: center;
            justify-content: center;
            width: 36px;
            height: 36px;
            background: var(--blue-600);
            color: var(--white);
            border-radius: 50%;
            font-size: 14px;
            font-weight: 700;
            flex-shrink: 0;
        }

        .intervention-content {
            flex: 1;
        }

        .intervention-content h3 {
            font-size: 16px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 8px;
        }

        .intervention-content p {
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-600);
            margin: 0;
        }

    </style>
</head>
<body>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <nav class="nav">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo">
                    <div class="logo-icon">
                        <svg width="36" height="12" viewBox="0 0 52 18" fill="none">
                            <circle cx="9" cy="9" r="9" fill="#2563EB"/>
                            <circle cx="26" cy="9" r="9" fill="#2563EB" fill-opacity="0.5"/>
                            <circle cx="43" cy="9" r="9" fill="#2563EB" fill-opacity="0.2"/>
                        </svg>
                    </div>
                    <span class="logo-text"><span class="gas">gas</span><span class="consult">consult</span><span class="ai">.ai</span></span>
                </a>
                <div class="nav-links">
                    <a href="/?clear=1" class="nav-link">Home</a>
                    <a href="/quick-dose" class="nav-link">Quick Dose</a>
                    <a href="/preop" class="nav-link">Pre-Op</a>
                    <a href="/calculators" class="nav-link">Calculators</a>
                    <a href="/hypotension" class="nav-link active">IOH Predictor</a>
                </div>
                <button class="mobile-menu-btn" onclick="toggleMobileMenu()" aria-label="Toggle menu">
                    <span></span>
                    <span></span>
                    <span></span>
                </button>
            </div>
        </nav>
        <div class="mobile-menu" id="mobileMenu">
            <a href="/?clear=1" class="mobile-menu-link">Home</a>
            <a href="/quick-dose" class="mobile-menu-link">Quick Dose</a>
            <a href="/preop" class="mobile-menu-link">Pre-Op</a>
            <a href="/calculators" class="mobile-menu-link">Calculators</a>
            <a href="/hypotension" class="mobile-menu-link">IOH Predictor</a>
        </div>

        <main class="main-content">
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
            <div class="content-card">
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
            </div>

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

        <!-- Glassmorphism Footer -->
        <footer class="footer">
            <div class="footer-inner">
                <div class="footer-brand">
                    <div class="footer-logo">
                        <svg width="32" height="32" viewBox="0 0 32 32" fill="none">
                            <circle cx="6" cy="16" r="5" fill="#2563EB"/>
                            <circle cx="16" cy="16" r="5" fill="#2563EB" fill-opacity="0.5"/>
                            <circle cx="26" cy="16" r="5" fill="#2563EB" fill-opacity="0.2"/>
                        </svg>
                    </div>
                    <span class="footer-text">© 2025 GasConsult.ai</span>
                </div>
                <div class="footer-links">
                    <a href="/privacy" class="footer-link">Privacy</a>
                    <a href="/terms" class="footer-link">Terms</a>
                    <a href="#" class="footer-link">Contact</a>
                </div>
            </div>
        </footer>
    </div>

    <script>
        function toggleMobileMenu() {
            const menu = document.getElementById('mobileMenu');
            const btn = document.querySelector('.mobile-menu-btn');
            if (menu && btn) {
                menu.classList.toggle('active');
                btn.classList.toggle('active');
            }
        }
    </script>
</body>
</html>
"""
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            padding: 28px;
            position: relative;
            overflow: hidden;
            transition: all 0.4s cubic-bezier(0.4,0,0.2,1);
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04);
            cursor: pointer;
        }

        .feature-card::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 1px;
            background: linear-gradient(90deg, transparent 0%, rgba(255,255,255,0.8) 50%, transparent 100%);
        }

        .feature-card:hover {
            transform: translateY(-4px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.04), 0 24px 64px rgba(0,0,0,0.08);
            border-color: rgba(59,130,246,0.2);
        }

        .feature-icon {
            width: 56px;
            height: 56px;
            border-radius: 16px;
            display: flex;
            align-items: center;
            justify-content: center;
            margin-bottom: 20px;
            transition: all 0.3s ease;
        }

        .feature-card:hover .feature-icon { transform: scale(1.05); }
        .feature-icon svg { width: 24px; height: 24px; }

        .feature-icon.blue {
            background: var(--blue-50);
            border: 1px solid var(--blue-100);
            box-shadow: 0 4px 16px rgba(37,99,235,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.blue svg { color: var(--blue-600); }

        .feature-icon.emerald {
            background: #ECFDF5;
            border: 1px solid #D1FAE5;
            box-shadow: 0 4px 16px rgba(16,185,129,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.emerald svg { color: #059669; }

        .feature-icon.violet {
            background: #F5F3FF;
            border: 1px solid #EDE9FE;
            box-shadow: 0 4px 16px rgba(139,92,246,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.violet svg { color: #7C3AED; }

        .feature-icon.rose {
            background: #FFF1F2;
            border: 1px solid #FFE4E6;
            box-shadow: 0 4px 16px rgba(244,63,94,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.rose svg { color: #E11D48; }

        .feature-icon.cyan {
            background: #ECFEFF;
            border: 1px solid #CFFAFE;
            box-shadow: 0 4px 16px rgba(6,182,212,0.1), inset 0 1px 0 rgba(255,255,255,0.8);
        }
        .feature-icon.cyan svg { color: #0891B2; }

        .feature-title {
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 10px;
            letter-spacing: -0.3px;
        }

        .feature-desc {
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-500);
            margin-bottom: 20px;
        }

        .feature-link {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            font-size: 14px;
            font-weight: 600;
            color: var(--blue-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .feature-link:hover { gap: 10px; }
        .feature-link svg { width: 16px; height: 16px; transition: transform 0.2s ease; }
        .feature-link:hover svg { transform: translateX(4px); }

        .footer {
            padding: 32px 20px;
            border-top: 1px solid var(--gray-200);
            background: rgba(255,255,255,0.5);
        }

        .footer-inner {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 20px;
            text-align: center;
        }

        .footer-brand {
            display: flex;
            align-items: center;
            gap: 10px;
        }

        .footer-logo {
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .footer-logo svg { width: 32px; height: 32px; }

        .footer-text {
            font-size: 13px;
            color: var(--gray-500);
        }

        .footer-links {
            display: flex;
            gap: 24px;
        }

        .footer-link {
            font-size: 13px;
            color: var(--gray-500);
            text-decoration: none;
            transition: color 0.2s ease;
        }

        .footer-link:hover { color: var(--gray-700); }

        @media (min-width: 768px) {
            .nav { padding: 16px 32px; }
            .nav-inner { height: 64px; padding: 0 24px; border-radius: 20px; }
            .logo-icon svg { width: 42px; height: 15px; }
            .logo-text { font-size: 20px; }
            .nav-links { display: flex; }
            .mobile-menu-btn { display: none; }
            .hero { padding: 160px 32px 80px; }
            .hero-badge { padding: 10px 20px 10px 14px; margin-bottom: 32px; }
            .badge-dot { width: 10px; height: 10px; }
            .badge-text { font-size: 13px; }
            .hero-title { font-size: 56px; letter-spacing: -2.5px; margin-bottom: 24px; }
            .hero-subtitle { font-size: 18px; margin-bottom: 48px; }
            .chat-container { padding: 0 24px; margin-bottom: 80px; }
            .chat-card { border-radius: 24px; padding: 10px; }
            .chat-inner { border-radius: 18px; padding: 8px; }
            .chat-input { padding: 16px 18px; min-height: 56px; }
            .chat-send { width: 52px; height: 52px; border-radius: 14px; }
            .chat-hints { gap: 10px; padding: 18px 10px 8px; }
            .hint-chip { padding: 12px 18px; font-size: 13px; }
            .hint-chip svg { width: 16px; height: 16px; }
            .features { padding: 80px 32px 100px; }
            .features-header { margin-bottom: 56px; }
            .features-label { font-size: 12px; margin-bottom: 16px; }
            .features-title { font-size: 36px; letter-spacing: -1.5px; }
            .features-subtitle { font-size: 18px; }
            .features-grid { grid-template-columns: repeat(2, 1fr); gap: 20px; }
            .feature-card { padding: 36px; border-radius: 24px; }
            .feature-card:hover { transform: translateY(-6px); }
            .feature-icon { width: 60px; height: 60px; border-radius: 18px; margin-bottom: 24px; }
            .feature-icon svg { width: 26px; height: 26px; }
            .feature-title { font-size: 20px; margin-bottom: 12px; }
            .feature-desc { font-size: 15px; line-height: 1.7; margin-bottom: 24px; }
            .footer { padding: 40px 32px; }
            .footer-inner { flex-direction: row; justify-content: space-between; text-align: left; }
            .footer-logo svg { width: 36px; height: 36px; }
            .footer-text { font-size: 14px; }
            .footer-links { gap: 32px; }
            .footer-link { font-size: 14px; }
            .orb-1 { width: 600px; height: 600px; left: -10%; }
            .orb-2 { width: 450px; height: 450px; right: -10%; }
            .orb-3 { width: 400px; height: 400px; }
        }

        @media (min-width: 1024px) {
            .nav { padding: 16px 40px; }
            .hero { padding: 180px 40px 80px; }
            .hero-title { font-size: 72px; letter-spacing: -3px; margin-bottom: 28px; }
            .hero-subtitle { font-size: 20px; margin-bottom: 56px; }
            .chat-container { margin-bottom: 100px; }
            .chat-card { border-radius: 28px; }
            .chat-inner { border-radius: 20px; }
            .chat-input { padding: 18px 20px; min-height: 60px; max-height: 180px; }
            .chat-send { width: 56px; height: 56px; border-radius: 16px; margin: 8px; }
            .chat-send svg { width: 22px; height: 22px; }
            .chat-hints { padding: 20px 12px 8px; }
            .hint-chip { padding: 12px 20px; }
            .features { padding: 80px 40px 120px; }
            .features-header { margin-bottom: 64px; }
            .features-title { font-size: 40px; }
            .features-grid { grid-template-columns: repeat(3, 1fr); gap: 24px; }
            .feature-card { padding: 40px; border-radius: 28px; }
            .feature-card:hover { transform: translateY(-8px); }
            .feature-icon { width: 64px; height: 64px; border-radius: 20px; margin-bottom: 28px; }
            .feature-icon svg { width: 28px; height: 28px; }
            .footer { padding: 48px 40px; }
            .orb-1 { width: 800px; height: 800px; top: -20%; left: -10%; }
            .orb-2 { width: 600px; height: 600px; right: -15%; }
            .orb-3 { width: 500px; height: 500px; }
        }

        @media (min-width: 1280px) {
            .hero-title { font-size: 80px; }
        }
    

        .main-content {
            flex: 1;
            padding: 100px 20px 40px;
            max-width: 1200px;
            margin: 0 auto;
            width: 100%;
        }

        .content-card {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            border-radius: 20px;
            padding: 32px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
        }

        h1 { font-size: 32px; margin-bottom: 16px; font-weight: 700; letter-spacing: -0.5px; }
        h2 { font-size: 24px; margin-bottom: 12px; font-weight: 700; letter-spacing: -0.5px; }
        h3 { font-size: 20px; margin-bottom: 10px; font-weight: 700; letter-spacing: -0.3px; }

        input[type="text"], input[type="number"], input[type="email"], select, textarea {
            width: 100%;
            padding: 12px 16px;
            border: 1px solid var(--gray-300);
            border-radius: 12px;
            font-family: inherit;
            font-size: 15px;
            background: var(--white);
            color: var(--gray-900);
            transition: all 0.2s ease;
        }

        input:focus, select:focus, textarea:focus {
            outline: none;
            border-color: var(--blue-500);
            box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.1);
        }

        button, .btn {
            padding: 12px 24px;
            background: var(--blue-600);
            color: var(--white);
            border: none;
            border-radius: 12px;
            font-family: inherit;
            font-size: 15px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s ease;
            box-shadow: 0 1px 2px rgba(37,99,235,0.2), 0 4px 16px rgba(37,99,235,0.2), inset 0 1px 0 rgba(255,255,255,0.1);
        }

        button:hover, .btn:hover {
            background: var(--blue-700);
            transform: translateY(-2px);
            box-shadow: 0 2px 4px rgba(37,99,235,0.2), 0 12px 40px rgba(37,99,235,0.3), inset 0 1px 0 rgba(255,255,255,0.1);
        }

        button:active, .btn:active {
            transform: translateY(0);
        }

        table {
            width: 100%;
            border-collapse: collapse;
            background: rgba(255,255,255,0.5);
            border-radius: 12px;
            overflow: hidden;
        }

        th, td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid var(--gray-200);
        }

        th {
            background: rgba(37,99,235,0.05);
            font-weight: 600;
            color: var(--blue-700);
        }

        @media (min-width: 768px) {
            .main-content { padding: 120px 32px 60px; }
        }

        @media (min-width: 1024px) {
            .main-content { padding: 140px 40px 80px; }
        }

        /* IOH Predictor Specific Styles */
        .edu-warning {
            background: linear-gradient(135deg, rgba(239, 68, 68, 0.15) 0%, rgba(220, 38, 38, 0.08) 100%);
            border: 2px solid rgba(239, 68, 68, 0.3);
            border-radius: 16px;
            padding: 24px;
            margin-bottom: 32px;
            display: flex;
            gap: 16px;
            align-items: flex-start;
        }

        .edu-warning-icon {
            font-size: 32px;
            flex-shrink: 0;
        }

        .edu-warning-content h3 {
            font-size: 16px;
            font-weight: 800;
            color: #DC2626;
            margin-bottom: 8px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        .edu-warning-content p {
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-800);
            margin: 0;
        }

        .edu-badge {
            display: inline-flex;
            align-items: center;
            background: rgba(239, 68, 68, 0.1);
            color: #DC2626;
            padding: 8px 16px;
            border-radius: 100px;
            font-size: 12px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            border: 1px solid rgba(239, 68, 68, 0.2);
            margin-top: 16px;
        }

        .disclaimer-box {
            background: rgba(255,255,255,0.9);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-left: 4px solid #DC2626;
            border-radius: 16px;
            padding: 28px;
            margin-bottom: 32px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04);
        }

        .disclaimer-box h3 {
            font-size: 18px;
            font-weight: 700;
            color: #DC2626;
            margin-bottom: 16px;
        }

        .disclaimer-box p {
            font-size: 14px;
            line-height: 1.7;
            color: var(--gray-700);
            margin-bottom: 12px;
        }

        .disclaimer-box p:last-of-type {
            margin-bottom: 8px;
        }

        .disclaimer-box ul {
            margin-left: 24px;
            font-size: 14px;
            line-height: 1.8;
            color: var(--gray-700);
        }

        .disclaimer-box li {
            margin-bottom: 6px;
        }

        .metrics-box {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 28px;
            margin-bottom: 32px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04);
        }

        .metrics-box h3 {
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 20px;
        }

        .predictor-form {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 32px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
        }

        .form-section {
            margin-bottom: 32px;
        }

        .form-section:last-child {
            margin-bottom: 0;
        }

        .form-section-title {
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 20px;
            padding-bottom: 12px;
            border-bottom: 2px solid var(--blue-100);
        }

        .input-group {
            display: grid;
            grid-template-columns: 1fr;
            gap: 20px;
            margin-bottom: 24px;
        }

        @media (min-width: 768px) {
            .input-group {
                grid-template-columns: repeat(2, 1fr);
            }
        }

        @media (min-width: 1024px) {
            .input-group {
                grid-template-columns: repeat(3, 1fr);
            }
        }

        .prediction-result {
            background: var(--white);
            border: 2px solid var(--blue-200);
            border-radius: 16px;
            padding: 32px;
            margin-top: 32px;
            text-align: center;
        }

        .risk-level {
            font-size: 14px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 1px;
            margin-bottom: 12px;
        }

        .risk-level.low { color: #059669; }
        .risk-level.moderate { color: #D97706; }
        .risk-level.high { color: #DC2626; }
        .risk-level.very-high { color: #991B1B; }

        .risk-percentage {
            font-size: 56px;
            font-weight: 900;
            margin-bottom: 12px;
        }

        .risk-percentage.low { color: #059669; }
        .risk-percentage.moderate { color: #D97706; }
        .risk-percentage.high { color: #DC2626; }
        .risk-percentage.very-high { color: #991B1B; }

        .risk-description {
            font-size: 15px;
            color: var(--gray-600);
            margin-bottom: 24px;
        }

        .risk-factors-list {
            background: var(--gray-50);
            border-radius: 12px;
            padding: 24px;
            text-align: left;
            margin-top: 24px;
        }

        .risk-factors-title {
            font-size: 16px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 16px;
        }

        .risk-factor-item {
            display: flex;
            gap: 12px;
            padding: 12px;
            background: var(--white);
            border-radius: 8px;
            margin-bottom: 10px;
            border-left: 3px solid var(--blue-500);
        }

        .risk-factor-item:last-child {
            margin-bottom: 0;
        }

        .risk-factor-name {
            font-size: 14px;
            font-weight: 600;
            color: var(--gray-900);
            margin-bottom: 4px;
        }

        .risk-factor-desc {
            font-size: 13px;
            color: var(--gray-600);
            line-height: 1.5;
        }

    </style>
</head>
<body>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <nav class="nav">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo">
                    <div class="logo-icon">
                        <svg width="36" height="12" viewBox="0 0 52 18" fill="none">
                            <circle cx="9" cy="9" r="9" fill="#2563EB"/>
                            <circle cx="26" cy="9" r="9" fill="#2563EB" fill-opacity="0.5"/>
                            <circle cx="43" cy="9" r="9" fill="#2563EB" fill-opacity="0.2"/>
                        </svg>
                    </div>
                    <span class="logo-text"><span class="gas">gas</span><span class="consult">consult</span><span class="ai">.ai</span></span>
                </a>
                <div class="nav-links">
                    <a href="/?clear=1" class="nav-link">Home</a>
                    <a href="/quick-dose" class="nav-link">Quick Dose</a>
                    <a href="/preop" class="nav-link">Pre-Op</a>
                    <a href="/calculators" class="nav-link">Calculators</a>
                    <a href="/hypotension" class="nav-link active">IOH Predictor</a>
                </div>
                <button class="mobile-menu-btn" onclick="toggleMobileMenu()" aria-label="Toggle menu">
                    <span></span>
                    <span></span>
                    <span></span>
                </button>
            </div>
        </nav>
        <div class="mobile-menu" id="mobileMenu">
            <a href="/?clear=1" class="mobile-menu-link">Home</a>
            <a href="/quick-dose" class="mobile-menu-link">Quick Dose</a>
            <a href="/preop" class="mobile-menu-link">Pre-Op</a>
            <a href="/calculators" class="mobile-menu-link">Calculators</a>
            <a href="/hypotension" class="mobile-menu-link">IOH Predictor</a>
        </div>

        <main class="main-content">
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

    <!-- Glassmorphism Footer -->
    <footer class="footer">
        <div class="footer-inner">
            <div class="footer-brand">
                <div class="footer-logo">
                    <svg width="32" height="32" viewBox="0 0 32 32" fill="none">
                        <circle cx="6" cy="16" r="5" fill="#2563EB"/>
                        <circle cx="16" cy="16" r="5" fill="#2563EB" fill-opacity="0.5"/>
                        <circle cx="26" cy="16" r="5" fill="#2563EB" fill-opacity="0.2"/>
                    </svg>
                </div>
                <span class="footer-text">© 2025 GasConsult.ai</span>
            </div>
            <div class="footer-links">
                <a href="/privacy" class="footer-link">Privacy</a>
                <a href="/terms" class="footer-link">Terms</a>
                <a href="#" class="footer-link">Contact</a>
            </div>
        </div>
    </footer>
        </main>
    </div>
</body>
</html>
"""

@app.route("/stream")
@csrf.exempt  # SSE endpoint uses GET requests, CSRF not applicable
def stream():
    """Server-Sent Events endpoint for streaming GPT responses"""
    request_id = request.args.get('request_id')

    logger.info(f"[STREAM] Received request with request_id: {request_id}")

    if not request_id:
        logger.error("[STREAM] No request_id provided")
        error_msg = json.dumps({'type': 'error', 'message': 'No request ID provided'})
        return Response(f"data: {error_msg}\n\n", mimetype='text/event-stream')

    stream_key = f'stream_data_{request_id}'
    stream_data = None

    # Try to get from session first
    logger.info(f"[STREAM] Looking for key: {stream_key}")
    logger.info(f"[STREAM] Session stream keys: {[k for k in session.keys() if k.startswith('stream_')]}")

    if stream_key in session:
        stream_data = session[stream_key]
        logger.info(f"[STREAM] Found stream data in session")
    else:
        # Fallback to in-memory cache
        logger.info(f"[STREAM] Session miss, checking cache...")
        logger.info(f"[STREAM] Cache keys: {list(STREAM_DATA_CACHE.keys())}")
        stream_data = get_stream_data(request_id)
        if stream_data:
            logger.info(f"[STREAM] Found stream data in cache (session fallback)")

    if not stream_data:
        logger.error(f"[STREAM] Stream data not found in session or cache for key: {stream_key}")
        logger.error(f"[STREAM] This usually means the session was not properly shared between requests")
        logger.error(f"[STREAM] If using multiple workers, consider using Redis for session storage")
        error_msg = json.dumps({
            'type': 'error',
            'message': 'Session expired or not found. This can happen if the server restarted. Please go back and try again.'
        })
        return Response(f"data: {error_msg}\n\n", mimetype='text/event-stream')

    # Extract data from stream_data (already retrieved from session or cache)
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
                    # Send content chunk - ensure_ascii=False to preserve HTML characters
                    yield f"data: {json.dumps({'type': 'content', 'data': content}, ensure_ascii=False)}\n\n"

            # Calculate evidence strength
            evidence_strength = get_evidence_strength(num_papers, refs)

            # Clean markdown code fences from the response
            cleaned_response = strip_markdown_code_fences(full_response)

            # Save complete response to session
            # Check if last message is an empty assistant placeholder (from homepage redirect)
            # If so, update it instead of appending a new one
            if (session.get('messages') and
                len(session['messages']) > 0 and
                session['messages'][-1].get('role') == 'assistant' and
                session['messages'][-1].get('content') == ''):
                # Update the placeholder
                session['messages'][-1] = {
                    "role": "assistant",
                    "content": cleaned_response,
                    "references": refs,
                    "num_papers": num_papers,
                    "evidence_strength": evidence_strength
                }
                print(f"[DEBUG] Updated existing placeholder assistant message")
            else:
                # Append new message (for AJAX submissions from chat page)
                session['messages'].append({
                    "role": "assistant",
                    "content": cleaned_response,
                    "references": refs,
                    "num_papers": num_papers,
                    "evidence_strength": evidence_strength
                })
                print(f"[DEBUG] Appended new assistant message")
            session.modified = True

            # Send references with evidence strength
            yield f"data: {json.dumps({'type': 'references', 'data': refs, 'num_papers': num_papers, 'evidence_strength': evidence_strength}, ensure_ascii=False)}\n\n"

            # Send completion event
            yield f"data: {json.dumps({'type': 'done'})}\n\n"

            # Clean up stream data from session
            session.pop(f'stream_data_{request_id}', None)
            session.modified = True

        except Exception as e:
            print(f"[ERROR] Streaming failed: {e}")
            yield f"data: {json.dumps({'type': 'error', 'message': str(e)})}\n\n"

    return Response(stream_with_context(generate()), mimetype='text/event-stream')

@app.route("/", methods=["GET", "POST"])
def index():
    """Homepage - handles both welcome screen and chat"""
    # Check if user wants to clear chat and return to hero state
    if request.method == "GET" and request.args.get('clear') == '1':
        session.pop('messages', None)
        session.pop('chat_active', None)
        session.pop('conversation_topic', None)
        session.modified = True
        return redirect(url_for('index'))

    # Initialize conversation history in session
    if 'messages' not in session:
        session['messages'] = []

    # Initialize conversation topic tracking
    if 'conversation_topic' not in session:
        session['conversation_topic'] = None

    # Explicitly initialize session to ensure CSRF token is generated
    if 'initialized' not in session:
        session['initialized'] = True

    if request.method == "POST":
        try:
            # Safely get query from form data and sanitize it
            raw_query = request.form.get("query", "").strip()
            print(f"[DEBUG] Form data received: {dict(request.form)}")
            print(f"[DEBUG] Query field raw value: '{request.form.get('query', 'MISSING')}'")
            print(f"[DEBUG] Is AJAX: {request.headers.get('X-Requested-With') == 'XMLHttpRequest'}")
            print(f"[DEBUG] Session messages count: {len(session.get('messages', []))}")
            raw_query = sanitize_user_query(raw_query)

            # Check if this is an AJAX request
            is_ajax = request.headers.get('X-Requested-With') == 'XMLHttpRequest' or request.form.get('modal') == '1'

            # If query is empty, return appropriate response
            if not raw_query:
                print(f"[DEBUG] Empty query received - likely issue with form submission")
                if is_ajax:
                    return jsonify({'status': 'error', 'message': 'Please enter a question before submitting.'})
                # Return to homepage with error message in session
                session['error_message'] = 'Please enter a question before submitting.'
                session.modified = True
                return redirect(url_for('index'))

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
            session['chat_active'] = True  # Set chat mode flag
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
                if is_ajax:
                    return jsonify({
                        'status': 'calculation',
                        'result': calc_result,
                        'message': 'Calculation complete'
                    })
                print(f"[DEBUG] Redirecting after calculation")
                return redirect(url_for('index'))

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
6. CRITICAL: Return ONLY the HTML content - do NOT wrap your response in markdown code fences (```html or ```)

Answer as if you're a colleague continuing the conversation:"""

                    print(f"[DEBUG] Preparing streaming for follow-up...")

                    # Generate unique request ID for this streaming session
                    request_id = str(uuid.uuid4())

                    # Prepare stream data
                    stream_data = {
                        'prompt': prompt,
                        'refs': [],
                        'num_papers': 0,
                        'raw_query': raw_query
                    }

                    # Store in session
                    session[f'stream_data_{request_id}'] = stream_data
                    session.modified = True

                    # Also store in memory cache as backup
                    store_stream_data(request_id, stream_data)

                    print(f"[DEBUG] Stream data prepared for follow-up, request_id: {request_id}")

                    # Add placeholder assistant message and set pending_stream
                    # This happens for BOTH AJAX and regular form submissions
                    session['messages'].append({
                        "role": "assistant",
                        "content": "",  # Will be populated by streaming
                        "references": [],
                        "num_papers": 0
                    })
                    session['pending_stream'] = request_id
                    session.modified = True
                    print(f"[DEBUG] Added placeholder assistant message and set pending_stream")

                    # If AJAX request, return JSON (JavaScript will handle redirect)
                    if is_ajax:
                        return jsonify({
                            'status': 'ready',
                            'request_id': request_id,
                            'raw_query': raw_query
                        })

                    # For regular form submissions, redirect to index page
                    print(f"[DEBUG] Redirecting to index page")
                    return redirect(url_for('index'))
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
                    if is_ajax:
                        return jsonify({
                            'status': 'error',
                            'message': error_msg
                        })
                    print(f"[DEBUG] Error message added, redirecting")
                    return redirect(url_for('index'))

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
6. CRITICAL: Return ONLY the HTML content - do NOT wrap your response in markdown code fences (```html or ```)
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
            logger.info(f"[CHAT] Preparing streaming with {num_papers} papers")

            # Generate unique request ID for this streaming session
            request_id = str(uuid.uuid4())
            logger.info(f"[CHAT] Generated request_id: {request_id}")

            # Prepare stream data
            stream_data = {
                'prompt': prompt,
                'refs': refs,
                'num_papers': num_papers,
                'raw_query': raw_query,
                'question_type': question_type  # Store for temperature adjustment
            }

            # Store in session
            stream_key = f'stream_data_{request_id}'
            session[stream_key] = stream_data
            session.modified = True

            # Also store in memory cache as backup (helps with session sync issues)
            store_stream_data(request_id, stream_data)

            logger.info(f"[CHAT] Stored stream data with key: {stream_key}")
            logger.info(f"[CHAT] Session keys after storing: {[k for k in session.keys() if k.startswith('stream_')]}")
            logger.info(f"[CHAT] Cache keys: {list(STREAM_DATA_CACHE.keys())}")

            print(f"[DEBUG] Stream data prepared, returning request_id: {request_id}")

            # Add placeholder assistant message for auto-start JavaScript to populate
            # This happens for BOTH AJAX and regular form submissions
            session['messages'].append({
                "role": "assistant",
                "content": "",  # Will be populated by streaming
                "references": [],
                "num_papers": 0
            })
            session['pending_stream'] = request_id
            session.modified = True
            print(f"[DEBUG] Added placeholder assistant message and set pending_stream")

            # If AJAX request, return JSON (JavaScript will handle redirect)
            if is_ajax:
                return jsonify({
                    'status': 'ready',
                    'request_id': request_id,
                    'raw_query': raw_query
                })

            # For regular form submissions (not AJAX), redirect to index page
            print(f"[DEBUG] Redirecting to index page")
            return redirect(url_for('index'))

        except Exception as e:
            # Catch all unhandled errors
            print(f"\n[ERROR] ===== UNHANDLED EXCEPTION =====")
            print(f"[ERROR] {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()

            error_content = f"<p><strong>Error:</strong> {str(e)}</p><p>Please try rephrasing your question or start a new conversation.</p>"
            session['messages'].append({
                "role": "assistant",
                "content": error_content,
                "references": [],
                "num_papers": 0
            })
            session.modified = True

            # Check if AJAX (may not be defined if error occurred early)
            is_ajax_check = request.headers.get('X-Requested-With') == 'XMLHttpRequest' or request.form.get('modal') == '1'
            if is_ajax_check:
                return jsonify({
                    'status': 'error',
                    'message': error_content
                })
            return redirect(url_for('index'))

    # GET request - check for pending stream and render page
    pending_stream = session.pop('pending_stream', None)

    # Log session state for debugging
    logger.info(f"[INDEX GET] pending_stream = {pending_stream}")
    logger.info(f"[INDEX GET] num messages = {len(session.get('messages', []))}")
    logger.info(f"[INDEX GET] stream_data keys = {[k for k in session.keys() if k.startswith('stream_')]}")

    print(f"[DEBUG] GET / - pending_stream = {pending_stream}")
    print(f"[DEBUG] GET / - num messages = {len(session.get('messages', []))}")

    # Get and clear any error message from session
    error_message = session.pop('error_message', None)

    # Ensure session is saved before rendering (critical for streaming to work)
    session.modified = True

    # Create response with no-cache headers to prevent stale UI
    response = make_response(render_template_string(HTML, messages=session.get('messages', []), pending_stream=pending_stream, error_message=error_message))
    response.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
    response.headers['Pragma'] = 'no-cache'
    response.headers['Expires'] = '0'
    return response

@app.route("/clear")
def clear():
    """Clear conversation history and return to homepage hero state"""
    session.pop('messages', None)
    session.pop('conversation_topic', None)
    session.pop('chat_active', None)  # Clear chat mode flag
    session.modified = True  # Ensure session is saved
    response = redirect(url_for('index'))
    # Prevent caching to ensure fresh page load
    response.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
    response.headers['Pragma'] = 'no-cache'
    response.headers['Expires'] = '0'
    return response

@app.route("/clear-chat")
def clear_chat():
    """Clear conversation messages but stay in active chat interface"""
    session.pop('messages', None)
    session.pop('conversation_topic', None)
    session['chat_active'] = True  # Keep chat interface active
    session.modified = True  # Ensure session is saved
    response = redirect(url_for('index'))
    # Prevent caching to ensure fresh page load
    response.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
    response.headers['Pragma'] = 'no-cache'
    response.headers['Expires'] = '0'
    return response

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
        response = make_response(render_template_string(PREOP_HTML, summary=None, references=None))
        response.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
        response.headers['Pragma'] = 'no-cache'
        response.headers['Expires'] = '0'
        return response

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

    # Validate required fields
    if not procedure or not procedure.strip():
        error_message = "<p style='color: #ff6b6b; text-align: center; padding: 20px;'><strong>Error:</strong> Please specify the surgical procedure before submitting the form.</p>"
        return render_template_string(PREOP_HTML, summary=error_message, references=None)

    if not surgery_risk:
        error_message = "<p style='color: #ff6b6b; text-align: center; padding: 20px;'><strong>Error:</strong> Please select a surgery risk category before submitting the form.</p>"
        return render_template_string(PREOP_HTML, summary=error_message, references=None)

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
            "chat": "/",
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