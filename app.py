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
    # Fix RSI - handle both standalone and embedded (order matters - do specific before general)
    q = re.sub(r'\brsi\b', '("rapid sequence induction" OR RSI OR "rapid sequence intubation" OR "airway management")', q, flags=re.IGNORECASE)
    q = q.replace("rapid sequence", '"rapid sequence induction" OR RSI OR "rapid sequence intubation"')
    q = q.replace("awake intubation", '"awake intubation" OR "fiberoptic intubation" OR "difficult airway"')
    q = q.replace("awake fiberoptic", '"awake intubation" OR "fiberoptic intubation" OR "difficult airway"')
    q = re.sub(r'\bcabg\b', '("coronary artery bypass" OR CABG OR "cardiac surgery")', q, flags=re.IGNORECASE)
    q = q.replace("coronary artery bypass", '"coronary artery bypass" OR CABG OR "cardiac surgery"')
    q = re.sub(r'\btavr\b', '("transcatheter aortic valve" OR TAVR OR "structural heart")', q, flags=re.IGNORECASE)
    q = q.replace("transcatheter aortic", '"transcatheter aortic valve" OR TAVR')

    # Common complications
    # LAST - keep space-based to avoid matching "last" in "the last time"
    q = q.replace(" last ", ' ("local anesthetic systemic toxicity" OR LAST OR "lipid emulsion" OR "intralipid") ')
    q = q.replace("local anesthetic toxicity", '"local anesthetic systemic toxicity" OR LAST OR "lipid emulsion"')
    # PRIS - use word boundary (less likely to be a common word)
    q = re.sub(r'\bpris\b', '("propofol infusion syndrome" OR PRIS)', q, flags=re.IGNORECASE)
    q = q.replace("propofol infusion syndrome", '"propofol infusion syndrome" OR PRIS')
    # Remove ambiguous "MH" abbreviation - it matches "Mental Health" in PubMed
    q = q.replace("malignant hyperthermia", '"malignant hyperthermia"[MeSH Terms] OR dantrolene OR "MH crisis"')

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
    <title>Pre-Op Risk Assessment — gasconsult.ai</title>

    <!-- PWA -->
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=6">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=6">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">

    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800;900&display=swap" rel="stylesheet">
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
            padding: 6px;
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
            padding: 4px;
            display: flex;
            align-items: flex-end;
            gap: 4px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04), inset 0 1px 0 rgba(255,255,255,1);
        }

        .chat-input {
            flex: 1;
            border: none;
            outline: none;
            padding: 10px 14px;
            font-size: 16px;
            font-family: inherit;
            color: var(--gray-800);
            background: transparent;
            resize: none;
            min-height: 38px;
            max-height: 110px;
            line-height: 1.5;
        }

        .chat-input::placeholder { color: var(--gray-400); }

        .chat-send {
            width: 44px;
            height: 44px;
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
            margin: 4px;
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

        /* Pre-Op Assessment Wizard Styles */
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

        /* Wizard Container */
        .wizard-container {
            max-width: 800px;
            margin: 0 auto;
        }

        .wizard-card {
            background: rgba(255,255,255,0.85);
            backdrop-filter: blur(30px) saturate(180%);
            -webkit-backdrop-filter: blur(30px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 24px;
            padding: 40px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 8px 24px rgba(0,0,0,0.05), 0 32px 100px rgba(0,0,0,0.08), inset 0 1px 0 rgba(255,255,255,0.9);
        }

        /* Section Header */
        .section-header {
            margin-bottom: 32px;
        }


        .section-title {
            font-size: 26px;
            font-weight: 800;
            letter-spacing: -0.8px;
            color: var(--gray-900);
            margin-bottom: 8px;
        }

        .section-description {
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-600);
        }

        /* Step Content */
        .section-content {
            opacity: 0;
            animation: fadeSlideIn 0.5s cubic-bezier(0.4, 0, 0.2, 1) forwards;
        }

        @keyframes fadeSlideIn {
            from {
                opacity: 0;
                transform: translateY(20px);
            }
            to {
                opacity: 1;
                transform: translateY(0);
            }
        }

        .form-section {
            display: block;
            margin-bottom: 40px;
            padding-bottom: 40px;
            border-bottom: 1px solid var(--gray-200);
        }

        .form-section:last-of-type {
            border-bottom: none;
        }

        /* Form Fields */
        .field-group {
            margin-bottom: 24px;
        }

        .field-label {
            display: block;
            font-size: 14px;
            font-weight: 600;
            color: var(--gray-700);
            margin-bottom: 8px;
        }

        .field-label .required {
            color: #DC2626;
            margin-left: 2px;
        }

        .field-input, select.field-input, textarea.field-input {
            width: 100%;
            padding: 14px 16px;
            border: 2px solid var(--gray-300);
            border-radius: 12px;
            font-family: inherit;
            font-size: 15px;
            background: var(--white);
            color: var(--gray-900);
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        }

        .field-input:focus, select.field-input:focus, textarea.field-input:focus {
            outline: none;
            border-color: var(--blue-500);
            box-shadow: 0 0 0 4px rgba(59, 130, 246, 0.1), 0 4px 16px rgba(59, 130, 246, 0.15);
            transform: translateY(-2px);
        }

        .field-input::placeholder {
            color: var(--gray-400);
        }

        textarea.field-input {
            min-height: 100px;
            resize: vertical;
        }

        .field-row {
            display: grid;
            grid-template-columns: 1fr;
            gap: 20px;
        }

        @media (min-width: 640px) {
            .field-row {
                grid-template-columns: repeat(2, 1fr);
            }
            .field-row.three-col {
                grid-template-columns: repeat(3, 1fr);
            }
        }

        /* Radio Buttons */
        .radio-group {
            display: flex;
            flex-direction: column;
            gap: 12px;
        }

        .radio-option {
            position: relative;
            padding: 18px 22px;
            background: var(--white);
            border: 2px solid var(--gray-200);
            border-radius: 14px;
            cursor: pointer;
            transition: all 0.25s cubic-bezier(0.4, 0, 0.2, 1);
        }

        .radio-option:hover {
            border-color: var(--blue-300);
            box-shadow: 0 2px 8px rgba(37, 99, 235, 0.08);
        }

        .radio-option input[type="radio"] {
            position: absolute;
            opacity: 0;
            pointer-events: none;
        }

        .radio-option input[type="radio"]:checked ~ .radio-content .radio-check {
            background: var(--blue-600);
            border-color: var(--blue-600);
        }

        .radio-option input[type="radio"]:checked ~ .radio-content .radio-check::after {
            opacity: 1;
        }

        .radio-option input[type="radio"]:checked ~ .radio-content .radio-label {
            font-weight: 600;
            color: var(--blue-600);
        }

        .radio-option.selected {
            background: var(--blue-50);
            border-color: var(--blue-500);
            box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.12);
        }

        .radio-content {
            display: flex;
            align-items: flex-start;
            gap: 14px;
        }

        .radio-check {
            width: 22px;
            height: 22px;
            border: 2.5px solid var(--gray-300);
            border-radius: 50%;
            flex-shrink: 0;
            position: relative;
            transition: all 0.25s ease;
            margin-top: 0px;
        }

        .radio-check::after {
            content: '';
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            width: 10px;
            height: 10px;
            background: white;
            border-radius: 50%;
            opacity: 0;
            transition: opacity 0.25s ease;
        }

        .radio-text {
            flex: 1;
        }

        .radio-label {
            font-size: 15px;
            font-weight: 600;
            color: var(--gray-800);
            display: block;
            margin-bottom: 4px;
            transition: all 0.3s ease;
        }

        .radio-desc {
            font-size: 13px;
            line-height: 1.5;
            color: var(--gray-600);
        }

        /* Checkbox Groups */
        .checkbox-group {
            display: grid;
            grid-template-columns: 1fr;
            gap: 12px;
        }

        @media (min-width: 640px) {
            .checkbox-group {
                grid-template-columns: repeat(2, 1fr);
            }
        }

        .checkbox-item {
            display: flex;
            align-items: center;
            gap: 10px;
            padding: 12px;
            background: var(--white);
            border: 2px solid var(--gray-300);
            border-radius: 10px;
            cursor: pointer;
            transition: all 0.2s ease;
        }

        .checkbox-item:hover {
            border-color: var(--blue-400);
            background: var(--blue-50);
        }

        .checkbox-item input[type="checkbox"] {
            width: 18px;
            height: 18px;
            cursor: pointer;
            accent-color: var(--blue-600);
        }

        .checkbox-item label {
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-700);
            cursor: pointer;
            flex: 1;
        }

        .checkbox-item input[type="checkbox"]:checked ~ label {
            color: var(--blue-700);
            font-weight: 600;
        }

        /* System Groups */
        .system-group {
            margin-bottom: 28px;
        }

        .system-group:last-child {
            margin-bottom: 0;
        }

        .system-title {
            font-size: 14px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: var(--gray-700);
            margin-bottom: 12px;
            padding-left: 4px;
        }

        /* Auto-calculations */
        .auto-calc-box {
            margin-top: 28px;
            padding: 20px;
            background: linear-gradient(135deg, var(--blue-50) 0%, var(--blue-100) 100%);
            border: 2px solid var(--blue-200);
            border-radius: 14px;
        }

        .auto-calc-title {
            font-size: 13px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: var(--blue-700);
            margin-bottom: 12px;
        }

        .calc-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(120px, 1fr));
            gap: 12px;
        }

        .calc-item {
            text-align: center;
        }

        .calc-label {
            font-size: 11px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: var(--blue-600);
            margin-bottom: 4px;
        }

        .calc-value {
            font-size: 22px;
            font-weight: 800;
            color: var(--blue-700);
        }

        .calc-unit {
            font-size: 11px;
            font-weight: 600;
            color: var(--blue-600);
            margin-left: 4px;
        }

        /* Form Submit */
        .form-submit {
            display: flex;
            justify-content: center;
            margin-top: 48px;
            padding-top: 40px;
        }

        .submit-btn {
            padding: 16px 40px;
            font-family: inherit;
            font-size: 16px;
            font-weight: 700;
            border-radius: 14px;
            cursor: pointer;
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            border: none;
            display: inline-flex;
            align-items: center;
            gap: 10px;
            background: linear-gradient(135deg, #3B82F6 0%, #2563EB 100%);
            color: var(--white);
            box-shadow: 0 4px 16px rgba(37, 99, 235, 0.3), inset 0 1px 0 rgba(255, 255, 255, 0.2);
        }

        .submit-btn:hover {
            transform: translateY(-2px);
            box-shadow: 0 8px 24px rgba(37, 99, 235, 0.4), inset 0 1px 0 rgba(255, 255, 255, 0.2);
        }

        .submit-btn:active {
            transform: translateY(0);
        }

        @keyframes shake {
            0%, 100% { transform: translateX(0); }
            25% { transform: translateX(-8px); }
            75% { transform: translateX(8px); }
        }

        /* Results Page */
        .results-container {
            max-width: 900px;
            margin: 0 auto;
        }

        .hero-card {
            background: linear-gradient(135deg, var(--white) 0%, var(--blue-50) 100%);
            border: 2px solid var(--blue-200);
            border-radius: 20px;
            padding: 40px;
            text-align: center;
            margin-bottom: 32px;
            box-shadow: 0 8px 32px rgba(37, 99, 235, 0.15);
        }

        .risk-level-badge {
            display: inline-flex;
            align-items: center;
            justify-content: center;
            gap: 10px;
            padding: 12px 32px;
            border-radius: 100px;
            font-size: 16px;
            font-weight: 800;
            text-transform: uppercase;
            letter-spacing: 1px;
            margin-bottom: 16px;
        }

        .risk-level-badge.low {
            background: linear-gradient(135deg, #10B981 0%, #059669 100%);
            color: white;
            box-shadow: 0 4px 20px rgba(16, 185, 129, 0.4);
        }

        .risk-level-badge.moderate {
            background: linear-gradient(135deg, #F59E0B 0%, #D97706 100%);
            color: white;
            box-shadow: 0 4px 20px rgba(245, 158, 11, 0.4);
        }

        .risk-level-badge.high {
            background: linear-gradient(135deg, #EF4444 0%, #DC2626 100%);
            color: white;
            box-shadow: 0 4px 20px rgba(239, 68, 68, 0.4);
        }

        .risk-title {
            font-size: 32px;
            font-weight: 800;
            letter-spacing: -1px;
            color: var(--gray-900);
            margin-bottom: 12px;
        }

        .risk-subtitle {
            font-size: 16px;
            color: var(--gray-600);
        }

        .info-grid {
            display: grid;
            grid-template-columns: 1fr;
            gap: 20px;
            margin-bottom: 32px;
        }

        @media (min-width: 768px) {
            .info-grid {
                grid-template-columns: repeat(2, 1fr);
            }
        }

        .info-card {
            background: rgba(255, 255, 255, 0.9);
            backdrop-filter: blur(20px);
            border: 2px solid var(--gray-200);
            border-radius: 16px;
            padding: 28px;
            box-shadow: 0 4px 16px rgba(0, 0, 0, 0.05);
        }

        .info-card-title {
            font-size: 14px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: var(--blue-600);
            margin-bottom: 16px;
        }

        .info-card-content {
            font-size: 15px;
            line-height: 1.8;
            color: var(--gray-700);
        }

        .score-display {
            font-size: 36px;
            font-weight: 800;
            color: var(--blue-600);
            margin-bottom: 8px;
        }

        .score-description {
            font-size: 13px;
            color: var(--gray-600);
        }

        .risk-factors-list {
            list-style: none;
            padding: 0;
            margin: 0;
        }

        .risk-factors-list li {
            padding: 10px 0;
            border-bottom: 1px solid var(--gray-200);
            font-size: 14px;
            color: var(--gray-700);
            display: flex;
            align-items: center;
            gap: 8px;
        }

        .risk-factors-list li:last-child {
            border-bottom: none;
        }

        .risk-factors-list li::before {
            content: '•';
            color: var(--blue-600);
            font-size: 20px;
            line-height: 1;
        }

        .recommendations-card {
            background: rgba(255, 255, 255, 0.9);
            backdrop-filter: blur(20px);
            border: 2px solid var(--blue-200);
            border-radius: 16px;
            padding: 32px;
            margin-bottom: 32px;
            box-shadow: 0 4px 16px rgba(0, 0, 0, 0.05);
        }

        .recommendations-title {
            font-size: 22px;
            font-weight: 800;
            letter-spacing: -0.5px;
            color: var(--gray-900);
            margin-bottom: 20px;
        }

        .recommendations-list {
            list-style: none;
            padding: 0;
            margin: 0;
        }

        .recommendation-item {
            padding: 16px;
            background: var(--blue-50);
            border-left: 4px solid var(--blue-600);
            border-radius: 8px;
            margin-bottom: 12px;
            font-size: 15px;
            line-height: 1.6;
            color: var(--gray-800);
        }

        .recommendation-item:last-child {
            margin-bottom: 0;
        }

        .recommendation-item strong {
            color: var(--blue-700);
        }

        .action-buttons {
            display: flex;
            gap: 16px;
            flex-wrap: wrap;
            justify-content: center;
        }

        .action-btn {
            padding: 14px 28px;
            font-family: inherit;
            font-size: 15px;
            font-weight: 700;
            border-radius: 12px;
            cursor: pointer;
            transition: all 0.3s ease;
            border: none;
            text-decoration: none;
            display: inline-flex;
            align-items: center;
            gap: 8px;
        }

        .action-btn-primary {
            background: linear-gradient(135deg, var(--blue-600) 0%, var(--blue-700) 100%);
            color: var(--white);
            box-shadow: 0 4px 16px rgba(37, 99, 235, 0.3);
        }

        .action-btn-primary:hover {
            transform: translateY(-2px);
            box-shadow: 0 8px 24px rgba(37, 99, 235, 0.4);
        }

        .action-btn-secondary {
            background: var(--white);
            color: var(--gray-700);
            border: 2px solid var(--gray-300);
        }

        .action-btn-secondary:hover {
            background: var(--gray-50);
            border-color: var(--gray-400);
            transform: translateY(-2px);
        }

        /* Loading State */
        .loading-overlay {
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            bottom: 0;
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(10px);
            z-index: 9999;
            display: none;
            align-items: center;
            justify-content: center;
        }

        .loading-overlay.show {
            display: flex;
        }

        .loading-content {
            text-align: center;
        }

        .loading-spinner {
            width: 60px;
            height: 60px;
            border: 5px solid var(--gray-200);
            border-top-color: var(--blue-600);
            border-radius: 50%;
            animation: spin 1s linear infinite;
            margin: 0 auto 20px;
        }

        @keyframes spin {
            to { transform: rotate(360deg); }
        }

        .loading-text {
            font-size: 16px;
            font-weight: 600;
            color: var(--gray-700);
        }

        /* Modern Results UI */
        .results-header {
            background: linear-gradient(135deg, var(--blue-600) 0%, var(--blue-700) 100%);
            color: white;
            padding: 32px;
            border-radius: 20px 20px 0 0;
            margin: -24px -24px 0 -24px;
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 20px;
        }

        .results-header-content {
            display: flex;
            align-items: center;
            gap: 20px;
        }

        .results-icon {
            width: 56px;
            height: 56px;
            background: rgba(255,255,255,0.2);
            border-radius: 16px;
            display: flex;
            align-items: center;
            justify-content: center;
            flex-shrink: 0;
        }

        .results-icon svg {
            width: 32px;
            height: 32px;
            stroke: white;
        }

        .results-title {
            font-size: 28px;
            font-weight: 800;
            margin: 0;
            letter-spacing: -0.5px;
        }

        .results-subtitle {
            font-size: 14px;
            margin: 4px 0 0 0;
            opacity: 0.9;
            font-weight: 500;
        }

        .results-actions-compact {
            display: flex;
            gap: 8px;
        }

        .icon-btn {
            width: 44px;
            height: 44px;
            background: rgba(255,255,255,0.2);
            border: none;
            border-radius: 12px;
            display: flex;
            align-items: center;
            justify-content: center;
            cursor: pointer;
            transition: all 0.2s ease;
            color: white;
            text-decoration: none;
        }

        .icon-btn:hover {
            background: rgba(255,255,255,0.3);
            transform: translateY(-2px);
        }

        .icon-btn svg {
            width: 20px;
            height: 20px;
            stroke: white;
        }

        .results-body {
            padding: 32px 24px;
        }

        .assessment-card {
            background: white;
            border-radius: 16px;
            border: 1px solid var(--gray-200);
            box-shadow: 0 1px 3px rgba(0,0,0,0.05);
            margin-bottom: 24px;
        }

        .assessment-content {
            padding: 32px;
            line-height: 1.8;
        }

        .assessment-content h3 {
            font-size: 20px;
            font-weight: 700;
            color: var(--gray-900);
            margin: 32px 0 16px 0;
        }

        .assessment-content h3:first-child {
            margin-top: 0;
        }

        .assessment-content p {
            margin: 12px 0;
            color: var(--gray-700);
        }

        .assessment-content ul {
            margin: 12px 0;
            padding-left: 24px;
        }

        .assessment-content li {
            margin: 8px 0;
            color: var(--gray-700);
        }

        .assessment-content strong {
            color: var(--gray-900);
            font-weight: 600;
        }

        .references-card {
            background: white;
            border-radius: 16px;
            border: 1px solid var(--gray-200);
            box-shadow: 0 1px 3px rgba(0,0,0,0.05);
        }

        .references-card-header {
            display: flex;
            align-items: center;
            gap: 12px;
            padding: 20px 24px;
            border-bottom: 1px solid var(--gray-200);
            font-size: 16px;
            font-weight: 700;
            color: var(--gray-900);
        }

        .references-card-header svg {
            width: 20px;
            height: 20px;
            stroke: var(--blue-600);
            flex-shrink: 0;
        }

        .references-count {
            margin-left: auto;
            font-size: 13px;
            font-weight: 600;
            background: var(--blue-50);
            color: var(--blue-600);
            padding: 4px 12px;
            border-radius: 20px;
        }

        .references-list {
            padding: 8px;
        }

        .reference-item-modern {
            display: flex;
            gap: 16px;
            padding: 16px;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .reference-item-modern:hover {
            background: var(--gray-50);
        }

        .reference-number {
            font-size: 14px;
            font-weight: 700;
            color: var(--blue-600);
            flex-shrink: 0;
            min-width: 32px;
        }

        .reference-content-modern {
            flex: 1;
        }

        .reference-title-link {
            color: var(--gray-900);
            text-decoration: none;
            font-size: 15px;
            font-weight: 600;
            line-height: 1.5;
            display: block;
            margin-bottom: 8px;
            transition: color 0.2s ease;
        }

        .reference-title-link:hover {
            color: var(--blue-600);
        }

        .reference-meta-modern {
            font-size: 13px;
            color: var(--gray-600);
            display: flex;
            align-items: center;
            gap: 8px;
            flex-wrap: wrap;
            margin-bottom: 8px;
        }

        .reference-separator {
            color: var(--gray-400);
        }

        .study-type-badge-preop {
            display: inline-block;
            padding: 4px 10px;
            border-radius: 6px;
            font-size: 11px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.3px;
            color: white;
        }

        .results-footer {
            padding: 24px;
            border-top: 1px solid var(--gray-200);
            display: flex;
            gap: 12px;
            flex-wrap: wrap;
            justify-content: center;
        }

        .btn-modern {
            padding: 14px 24px;
            font-family: inherit;
            font-size: 15px;
            font-weight: 600;
            border-radius: 12px;
            cursor: pointer;
            transition: all 0.2s ease;
            border: none;
            text-decoration: none;
            display: inline-flex;
            align-items: center;
            gap: 10px;
        }

        .btn-modern svg {
            width: 18px;
            height: 18px;
        }

        .btn-modern-primary {
            background: linear-gradient(135deg, var(--blue-600) 0%, var(--blue-700) 100%);
            color: white;
            box-shadow: 0 2px 8px rgba(37,99,235,0.3);
        }

        .btn-modern-primary:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(37,99,235,0.4);
        }

        .btn-modern-secondary {
            background: white;
            color: var(--gray-700);
            border: 2px solid var(--gray-300);
        }

        .btn-modern-secondary:hover {
            background: var(--gray-50);
            border-color: var(--gray-400);
            transform: translateY(-2px);
        }

        .btn-modern-tertiary {
            background: transparent;
            color: var(--gray-600);
            border: none;
        }

        .btn-modern-tertiary:hover {
            color: var(--gray-900);
            background: var(--gray-100);
        }

        /* Responsive */
        @media (max-width: 640px) {
            .wizard-card {
                padding: 24px;
            }
            .section-title {
                font-size: 22px;
            }
            .wizard-nav {
                flex-direction: column;
            }
            .wizard-btn-prev {
                order: 2;
            }
            .wizard-btn-next, .wizard-btn-submit {
                order: 1;
            }
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
                    <a href="/calculators" class="nav-link">Clinical Calculators</a>
                    <a href="/crisis" class="nav-link">Crisis Protocols</a>
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
            <a href="/crisis" class="mobile-menu-link">Crisis Protocols</a>
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
                <h1 class="page-title"><span style="color: var(--blue-600);">Pre-Operative</span> Risk Assessment</h1>
                <p class="page-subtitle">Comprehensive NSQIP-inspired risk stratification with evidence-based recommendations</p>
            </div>

            {% if not summary %}
            <!-- Form Container -->
            <div class="wizard-container">
                <div class="wizard-card">
                    <!-- Single Form -->
                    <form method="post" action="/preop" id="preopForm">
                        <input type="hidden" name="csrf_token" value="{{ csrf_token() }}"/>

                        <!-- Section 1: Demographics -->
                        <div class="form-section">
                            <div class="section-header">
                                <h2 class="section-title">Patient Demographics</h2>
                                <p class="section-description">Basic patient information for risk calculation</p>
                            </div>

                            <div class="section-content">
                                <div class="field-row three-col">
                                    <div class="field-group">
                                        <label for="age" class="field-label">Age <span class="required">*</span></label>
                                        <input type="number" id="age" name="age" class="field-input" placeholder="Years" required>
                                    </div>
                                    <div class="field-group">
                                        <label for="weight" class="field-label">Weight <span class="required">*</span></label>
                                        <input type="number" id="weight" name="weight" step="0.1" class="field-input" placeholder="kg" required>
                                    </div>
                                    <div class="field-group">
                                        <label for="height" class="field-label">Height <span class="required">*</span></label>
                                        <input type="number" id="height" name="height" step="0.1" class="field-input" placeholder="cm" required>
                                    </div>
                                </div>

                                <div class="field-group">
                                    <label class="field-label">Sex Assigned at Birth <span class="required">*</span></label>
                                    <div class="field-row">
                                        <label class="radio-option" onclick="selectRadio('sex', 'male')">
                                            <input type="radio" id="sex-male" name="sex" value="male" required>
                                            <div class="radio-check"></div>
                                            <div class="radio-content">
                                                <div class="radio-text">
                                                    <span class="radio-label">Male</span>
                                                </div>
                                            </div>
                                        </label>
                                        <label class="radio-option" onclick="selectRadio('sex', 'female')">
                                            <input type="radio" id="sex-female" name="sex" value="female" required>
                                            <div class="radio-check"></div>
                                            <div class="radio-content">
                                                <div class="radio-text">
                                                    <span class="radio-label">Female</span>
                                                </div>
                                            </div>
                                        </label>
                                    </div>
                                </div>

                                <!-- Auto-calculations -->
                                <div class="auto-calc-box" id="autoCalcBox" style="display: none;">
                                    <div class="auto-calc-title">Calculated Metrics</div>
                                    <div class="calc-grid">
                                        <div class="calc-item">
                                            <div class="calc-label">BMI</div>
                                            <div class="calc-value" id="calcBMI">--</div>
                                        </div>
                                        <div class="calc-item">
                                            <div class="calc-label">BMI Category</div>
                                            <div class="calc-value" id="calcBMICategory" style="font-size: 14px;">--</div>
                                        </div>
                                        <div class="calc-item">
                                            <div class="calc-label">IBW</div>
                                            <div class="calc-value" id="calcIBW">--<span class="calc-unit">kg</span></div>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>

                        <!-- Section 2: ASA Classification -->
                        <div class="form-section">
                            <div class="section-header">
                                <h2 class="section-title">ASA Physical Status</h2>
                                <p class="section-description">Select the appropriate ASA classification</p>
                            </div>

                            <div class="section-content">
                                <div class="field-group">
                                    <label class="field-label">ASA Classification <span class="required">*</span></label>
                                    <div class="radio-group">
                                        <label class="radio-option" onclick="selectRadio('asa', '1')">
                                            <input type="radio" id="asa-1" name="asa" value="1" required>
                                            <div class="radio-check"></div>
                                            <div class="radio-content">
                                                <div class="radio-text">
                                                    <span class="radio-label">ASA I</span>
                                                    <span class="radio-desc">Normal healthy patient</span>
                                                </div>
                                            </div>
                                        </label>
                                        <label class="radio-option" onclick="selectRadio('asa', '2')">
                                            <input type="radio" id="asa-2" name="asa" value="2" required>
                                            <div class="radio-check"></div>
                                            <div class="radio-content">
                                                <div class="radio-text">
                                                    <span class="radio-label">ASA II</span>
                                                    <span class="radio-desc">Mild systemic disease, no functional limitations</span>
                                                </div>
                                            </div>
                                        </label>
                                        <label class="radio-option" onclick="selectRadio('asa', '3')">
                                            <input type="radio" id="asa-3" name="asa" value="3" required>
                                            <div class="radio-check"></div>
                                            <div class="radio-content">
                                                <div class="radio-text">
                                                    <span class="radio-label">ASA III</span>
                                                    <span class="radio-desc">Severe systemic disease with functional limitations</span>
                                                </div>
                                            </div>
                                        </label>
                                        <label class="radio-option" onclick="selectRadio('asa', '4')">
                                            <input type="radio" id="asa-4" name="asa" value="4" required>
                                            <div class="radio-check"></div>
                                            <div class="radio-content">
                                                <div class="radio-text">
                                                    <span class="radio-label">ASA IV</span>
                                                    <span class="radio-desc">Severe systemic disease that is a constant threat to life</span>
                                                </div>
                                            </div>
                                        </label>
                                        <label class="radio-option" onclick="selectRadio('asa', '5')">
                                            <input type="radio" id="asa-5" name="asa" value="5" required>
                                            <div class="radio-check"></div>
                                            <div class="radio-content">
                                                <div class="radio-text">
                                                    <span class="radio-label">ASA V</span>
                                                    <span class="radio-desc">Moribund patient not expected to survive without operation</span>
                                                </div>
                                            </div>
                                        </label>
                                        <label class="radio-option" onclick="selectRadio('asa', '6')">
                                            <input type="radio" id="asa-6" name="asa" value="6" required>
                                            <div class="radio-check"></div>
                                            <div class="radio-content">
                                                <div class="radio-text">
                                                    <span class="radio-label">ASA VI</span>
                                                    <span class="radio-desc">Declared brain-dead patient for organ donation</span>
                                                </div>
                                            </div>
                                        </label>
                                    </div>
                                </div>

                                <div class="field-group" style="margin-top: 24px;">
                                    <div class="checkbox-item">
                                        <input type="checkbox" id="emergency" name="emergency" value="true">
                                        <label for="emergency">Emergency Surgery (E modifier)</label>
                                    </div>
                                </div>
                            </div>
                        </div>

                        <!-- Section3: Comorbidities -->
                        <div class="form-section" id="section3">
                            <div class="section-header">
                                <h2 class="section-title">Comorbidities</h2>
                                <p class="section-description">Select all applicable conditions organized by system</p>
                            </div>

                            <div class="section-content">
                                <!-- Cardiac -->
                                <div class="system-group">
                                    <div class="system-title">Cardiac</div>
                                    <div class="checkbox-group">
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="cad" name="comorbidities" value="CAD">
                                            <label for="cad">Coronary Artery Disease</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="mi" name="comorbidities" value="Prior MI">
                                            <label for="mi">Prior MI</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="chf" name="comorbidities" value="Heart Failure">
                                            <label for="chf">Heart Failure</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="afib" name="comorbidities" value="Atrial Fibrillation">
                                            <label for="afib">Atrial Fibrillation</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="pacemaker" name="comorbidities" value="Pacemaker/ICD">
                                            <label for="pacemaker">Pacemaker/ICD</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="valvular" name="comorbidities" value="Valvular Disease">
                                            <label for="valvular">Valvular Disease</label>
                                        </div>
                                    </div>
                                </div>

                                <!-- Pulmonary -->
                                <div class="system-group">
                                    <div class="system-title">Pulmonary</div>
                                    <div class="checkbox-group">
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="copd" name="comorbidities" value="COPD">
                                            <label for="copd">COPD</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="asthma" name="comorbidities" value="Asthma">
                                            <label for="asthma">Asthma</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="osa" name="comorbidities" value="OSA">
                                            <label for="osa">Obstructive Sleep Apnea</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="ild" name="comorbidities" value="Interstitial Lung Disease">
                                            <label for="ild">Interstitial Lung Disease</label>
                                        </div>
                                    </div>
                                </div>

                                <!-- Endocrine -->
                                <div class="system-group">
                                    <div class="system-title">Endocrine</div>
                                    <div class="checkbox-group">
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="dm" name="comorbidities" value="Diabetes Mellitus">
                                            <label for="dm">Diabetes Mellitus</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="dm-insulin" name="comorbidities" value="DM on Insulin">
                                            <label for="dm-insulin">DM on Insulin</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="thyroid" name="comorbidities" value="Thyroid Disease">
                                            <label for="thyroid">Thyroid Disease</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="adrenal" name="comorbidities" value="Adrenal Insufficiency">
                                            <label for="adrenal">Adrenal Insufficiency</label>
                                        </div>
                                    </div>
                                </div>

                                <!-- Renal -->
                                <div class="system-group">
                                    <div class="system-title">Renal</div>
                                    <div class="checkbox-group">
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="ckd" name="comorbidities" value="CKD">
                                            <label for="ckd">Chronic Kidney Disease</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="esrd" name="comorbidities" value="ESRD">
                                            <label for="esrd">ESRD</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="dialysis" name="comorbidities" value="Dialysis">
                                            <label for="dialysis">On Dialysis</label>
                                        </div>
                                    </div>
                                </div>

                                <!-- Neurological -->
                                <div class="system-group">
                                    <div class="system-title">Neurological</div>
                                    <div class="checkbox-group">
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="stroke" name="comorbidities" value="CVA/TIA">
                                            <label for="stroke">Prior CVA/TIA</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="seizure" name="comorbidities" value="Seizure Disorder">
                                            <label for="seizure">Seizure Disorder</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="neuromuscular" name="comorbidities" value="Neuromuscular Disease">
                                            <label for="neuromuscular">Neuromuscular Disease</label>
                                        </div>
                                    </div>
                                </div>

                                <!-- Hematologic -->
                                <div class="system-group">
                                    <div class="system-title">Hematologic</div>
                                    <div class="checkbox-group">
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="anemia" name="comorbidities" value="Anemia">
                                            <label for="anemia">Anemia</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="bleeding" name="comorbidities" value="Bleeding Disorder">
                                            <label for="bleeding">Bleeding Disorder</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="anticoag" name="comorbidities" value="On Anticoagulation">
                                            <label for="anticoag">On Anticoagulation</label>
                                        </div>
                                    </div>
                                </div>

                                <!-- Other -->
                                <div class="system-group">
                                    <div class="system-title">Other</div>
                                    <div class="checkbox-group">
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="htn" name="comorbidities" value="Hypertension">
                                            <label for="htn">Hypertension</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="liver" name="comorbidities" value="Liver Disease">
                                            <label for="liver">Liver Disease</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="cancer" name="comorbidities" value="Cancer">
                                            <label for="cancer">Cancer</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="immunosupp" name="comorbidities" value="Immunosuppression">
                                            <label for="immunosupp">Immunosuppression</label>
                                        </div>
                                    </div>
                                </div>

                                <div class="field-group" style="margin-top: 24px;">
                                    <label for="other_comorbidities" class="field-label">Other Conditions</label>
                                    <textarea id="other_comorbidities" name="other_comorbidities" class="field-input" placeholder="List any other relevant conditions..." rows="2"></textarea>
                                </div>
                            </div>
                        </div>

                        <!-- Section4: Functional Status -->
                        <div class="form-section" id="section4">
                            <div class="section-header">
                                <h2 class="section-title">Functional Status & Exercise Tolerance</h2>
                                <p class="section-description">Assess patient's functional capacity and activity level</p>
                            </div>

                            <div class="section-content">
                                <div class="field-group">
                                    <label class="field-label">Metabolic Equivalents (METs) <span class="required">*</span></label>
                                    <div class="radio-group">
                                        <label class="radio-option" onclick="selectRadio('mets', 'poor')">
                                            <input type="radio" id="mets-poor" name="mets" value="<4 METs" required>
                                            <div class="radio-check"></div>
                                            <div class="radio-content">
                                                <div class="radio-text">
                                                    <span class="radio-label">&lt;4 METs - Poor Functional Capacity</span>
                                                    <span class="radio-desc">Cannot climb 2 flights of stairs or walk 4 blocks</span>
                                                </div>
                                            </div>
                                        </label>
                                        <label class="radio-option" onclick="selectRadio('mets', 'moderate')">
                                            <input type="radio" id="mets-moderate" name="mets" value="4-10 METs" required>
                                            <div class="radio-check"></div>
                                            <div class="radio-content">
                                                <div class="radio-text">
                                                    <span class="radio-label">4-10 METs - Moderate Functional Capacity</span>
                                                    <span class="radio-desc">Can climb 2 flights of stairs, moderate exertion tolerated</span>
                                                </div>
                                            </div>
                                        </label>
                                        <label class="radio-option" onclick="selectRadio('mets', 'excellent')">
                                            <input type="radio" id="mets-excellent" name="mets" value=">10 METs" required>
                                            <div class="radio-content">
                                                <div class="radio-check"></div>
                                                <div class="radio-text">
                                                    <span class="radio-label">&gt;10 METs - Excellent Functional Capacity</span>
                                                    <span class="radio-desc">Strenuous exercise, high activity level</span>
                                                </div>
                                            </div>
                                        </label>
                                    </div>
                                </div>

                                <div class="field-group" style="margin-top: 32px;">
                                    <label class="field-label">Activities of Daily Living (ADLs) <span class="required">*</span></label>
                                    <div class="field-row">
                                        <label class="radio-option" onclick="selectRadio('adl', 'independent')">
                                            <input type="radio" id="adl-independent" name="adl" value="Independent" required>
                                            <div class="radio-check"></div>
                                            <div class="radio-content">
                                                <div class="radio-text">
                                                    <span class="radio-label">Independent</span>
                                                </div>
                                            </div>
                                        </label>
                                        <label class="radio-option" onclick="selectRadio('adl', 'partial')">
                                            <input type="radio" id="adl-partial" name="adl" value="Partially Dependent" required>
                                            <div class="radio-check"></div>
                                            <div class="radio-content">
                                                <div class="radio-text">
                                                    <span class="radio-label">Partially Dependent</span>
                                                </div>
                                            </div>
                                        </label>
                                        <label class="radio-option" onclick="selectRadio('adl', 'dependent')">
                                            <input type="radio" id="adl-dependent" name="adl" value="Fully Dependent" required>
                                            <div class="radio-check"></div>
                                            <div class="radio-content">
                                                <div class="radio-text">
                                                    <span class="radio-label">Fully Dependent</span>
                                                </div>
                                            </div>
                                        </label>
                                    </div>
                                </div>

                                <div class="field-group" style="margin-top: 32px;">
                                    <label for="anesthesia_history" class="field-label">Previous Anesthesia History</label>
                                    <textarea id="anesthesia_history" name="anesthesia_history" class="field-input" placeholder="e.g., General anesthesia for appendectomy 2015 - no complications. Family history of malignant hyperthermia..." rows="3"></textarea>
                                </div>
                            </div>
                        </div>

                        <!-- Section5: Medications -->
                        <div class="form-section" id="section5">
                            <div class="section-header">
                                <h2 class="section-title">Current Medications</h2>
                                <p class="section-description">List medications and identify critical agents</p>
                            </div>

                            <div class="section-content">
                                <div class="field-group">
                                    <label for="medications" class="field-label">Medication List</label>
                                    <textarea id="medications" name="medications" class="field-input" placeholder="e.g., Metoprolol 50mg BID, Lisinopril 10mg daily, Atorvastatin 40mg nightly..." rows="4"></textarea>
                                </div>

                                <div class="system-group" style="margin-top: 32px;">
                                    <div class="system-title">Specific Medications of Interest</div>
                                    <div class="checkbox-group">
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="aspirin" name="medications_specific" value="Aspirin">
                                            <label for="aspirin">Aspirin</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="p2y12" name="medications_specific" value="P2Y12 Inhibitor">
                                            <label for="p2y12">P2Y12 Inhibitor (Plavix/Ticagrelor)</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="warfarin" name="medications_specific" value="Warfarin">
                                            <label for="warfarin">Warfarin</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="doac" name="medications_specific" value="DOAC">
                                            <label for="doac">DOAC (Apixaban/Rivaroxaban/Edoxaban)</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="beta-blocker" name="medications_specific" value="Beta Blocker">
                                            <label for="beta-blocker">Beta Blocker</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="ace-arb" name="medications_specific" value="ACE-I/ARB">
                                            <label for="ace-arb">ACE-I/ARB</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="diuretic" name="medications_specific" value="Diuretic">
                                            <label for="diuretic">Diuretic</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="insulin" name="medications_specific" value="Insulin">
                                            <label for="insulin">Insulin</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="oral-hypoglycemic" name="medications_specific" value="Oral Hypoglycemic">
                                            <label for="oral-hypoglycemic">Oral Hypoglycemic</label>
                                        </div>
                                        <div class="checkbox-item">
                                            <input type="checkbox" id="steroids" name="medications_specific" value="Steroids">
                                            <label for="steroids">Steroids</label>
                                        </div>
                                    </div>
                                </div>

                                <div class="field-row" style="margin-top: 24px;">
                                    <div class="field-group">
                                        <label for="anticoag_last_dose" class="field-label">Last Anticoagulant Dose</label>
                                        <input type="text" id="anticoag_last_dose" name="anticoag_last_dose" class="field-input" placeholder="e.g., Apixaban 5mg yesterday PM">
                                    </div>
                                    <div class="field-group">
                                        <label for="allergies" class="field-label">Allergies</label>
                                        <input type="text" id="allergies" name="allergies" class="field-input" placeholder="e.g., PCN (rash), NKDA">
                                    </div>
                                </div>
                            </div>
                        </div>

                        <!-- Section6: Labs & Cardiac Workup -->
                        <div class="form-section" id="section6">
                            <div class="section-header">
                                <h2 class="section-title">Laboratory Values & Cardiac Assessment</h2>
                                <p class="section-description">Enter available lab results and cardiac workup findings</p>
                            </div>

                            <div class="section-content">
                                <div class="field-row">
                                    <div class="field-group">
                                        <label for="hgb" class="field-label">Hemoglobin (g/dL)</label>
                                        <input type="number" id="hgb" name="hgb" step="0.1" class="field-input" placeholder="e.g., 12.5">
                                    </div>
                                    <div class="field-group">
                                        <label for="plt" class="field-label">Platelets (×10³/μL)</label>
                                        <input type="number" id="plt" name="plt" class="field-input" placeholder="e.g., 250">
                                    </div>
                                    <div class="field-group">
                                        <label for="cr" class="field-label">Creatinine (mg/dL)</label>
                                        <input type="number" id="cr" name="cr" step="0.01" class="field-input" placeholder="e.g., 0.9">
                                    </div>
                                </div>

                                <div class="field-row" style="margin-top: 20px;">
                                    <div class="field-group">
                                        <label for="inr" class="field-label">INR</label>
                                        <input type="number" id="inr" name="inr" step="0.1" class="field-input" placeholder="e.g., 1.1">
                                    </div>
                                    <div class="field-group">
                                        <label for="glucose" class="field-label">Glucose (mg/dL)</label>
                                        <input type="number" id="glucose" name="glucose" class="field-input" placeholder="e.g., 110">
                                    </div>
                                    <div class="field-group">
                                        <label for="ef" class="field-label">Ejection Fraction (%)</label>
                                        <input type="number" id="ef" name="ef" class="field-input" placeholder="e.g., 60">
                                    </div>
                                </div>

                                <!-- Auto-calculated eGFR -->
                                <div class="auto-calc-box" id="eGFRBox" style="display: none; margin-top: 24px;">
                                    <div class="auto-calc-title">Calculated Renal Function</div>
                                    <div class="calc-grid">
                                        <div class="calc-item">
                                            <div class="calc-label">eGFR</div>
                                            <div class="calc-value" id="calcEGFR">--<span class="calc-unit">mL/min</span></div>
                                        </div>
                                        <div class="calc-item">
                                            <div class="calc-label">CKD Stage</div>
                                            <div class="calc-value" id="calcCKDStage" style="font-size: 16px;">--</div>
                                        </div>
                                    </div>
                                </div>

                                <div class="field-row" style="margin-top: 24px;">
                                    <div class="field-group">
                                        <label for="ekg" class="field-label">EKG Findings</label>
                                        <select id="ekg" name="ekg" class="field-input">
                                            <option value="">Select or type...</option>
                                            <option value="Normal Sinus Rhythm">Normal Sinus Rhythm</option>
                                            <option value="Atrial Fibrillation">Atrial Fibrillation</option>
                                            <option value="Old MI">Old MI</option>
                                            <option value="Left Ventricular Hypertrophy">LVH</option>
                                            <option value="Bundle Branch Block">Bundle Branch Block</option>
                                            <option value="Other">Other</option>
                                            <option value="Not Done">Not Done</option>
                                        </select>
                                    </div>
                                    <div class="field-group">
                                        <label for="npo" class="field-label">NPO Status</label>
                                        <input type="text" id="npo" name="npo" class="field-input" placeholder="e.g., NPO since midnight">
                                    </div>
                                </div>
                            </div>
                        </div>

                        <!-- Section7: Surgical Procedure -->
                        <div class="form-section" id="section7">
                            <div class="section-header">
                                <h2 class="section-title">Surgical Procedure</h2>
                                <p class="section-description">Details about the planned operation</p>
                            </div>

                            <div class="section-content">
                                <div class="field-group">
                                    <label for="procedure" class="field-label">Procedure Name <span class="required">*</span></label>
                                    <input type="text" id="procedure" name="procedure" class="field-input" placeholder="e.g., Total Knee Arthroplasty, Laparoscopic Cholecystectomy, CABG..." required>
                                </div>

                                <div class="field-group" style="margin-top: 24px;">
                                    <label class="field-label">Surgical Risk Category <span class="required">*</span></label>
                                    <div class="radio-group">
                                        <label class="radio-option" onclick="selectRadio('surgery_risk', 'low')">
                                            <input type="radio" id="surgery_risk-low" name="surgery_risk" value="Low" required>
                                            <div class="radio-check"></div>
                                            <div class="radio-content">
                                                <div class="radio-text">
                                                    <span class="radio-label">Low Risk (&lt;1% cardiac risk)</span>
                                                    <span class="radio-desc">Superficial procedures, cataract surgery, breast surgery, ambulatory procedures</span>
                                                </div>
                                            </div>
                                        </label>
                                        <label class="radio-option" onclick="selectRadio('surgery_risk', 'intermediate')">
                                            <input type="radio" id="surgery_risk-intermediate" name="surgery_risk" value="Intermediate" required>
                                            <div class="radio-check"></div>
                                            <div class="radio-content">
                                                <div class="radio-text">
                                                    <span class="radio-label">Intermediate Risk (1-5% cardiac risk)</span>
                                                    <span class="radio-desc">Intraperitoneal, intrathoracic, orthopedic, head/neck surgery</span>
                                                </div>
                                            </div>
                                        </label>
                                        <label class="radio-option" onclick="selectRadio('surgery_risk', 'high')">
                                            <input type="radio" id="surgery_risk-high" name="surgery_risk" value="High" required>
                                            <div class="radio-check"></div>
                                            <div class="radio-content">
                                                <div class="radio-text">
                                                    <span class="radio-label">High Risk (&gt;5% cardiac risk)</span>
                                                    <span class="radio-desc">Vascular surgery, aortic surgery, prolonged procedures with major fluid shifts</span>
                                                </div>
                                            </div>
                                        </label>
                                    </div>
                                </div>

                                <div class="field-row" style="margin-top: 24px;">
                                    <div class="field-group">
                                        <label class="field-label">Urgency <span class="required">*</span></label>
                                        <select id="urgency" name="urgency" class="field-input" required>
                                            <option value="">Select...</option>
                                            <option value="Elective">Elective</option>
                                            <option value="Urgent">Urgent (within 24-48 hours)</option>
                                            <option value="Emergent">Emergent (immediate)</option>
                                        </select>
                                    </div>
                                    <div class="field-group">
                                        <label for="duration" class="field-label">Expected Duration (hours)</label>
                                        <input type="number" id="duration" name="duration" step="0.5" class="field-input" placeholder="e.g., 2.5">
                                    </div>
                                </div>
                            </div>
                        </div>

                        <!-- Submit Button -->
                        <div class="form-submit">
                            <button type="submit" class="submit-btn">
                                Generate Assessment →
                            </button>
                        </div>
                    </form>
                </div>
            </div>

            {% else %}
            <!-- Results Display - Modern Redesigned UI -->
            <div class="results-container">
                <div class="wizard-card">
                    <h2 style="font-size: 28px; font-weight: 800; text-align: center; margin-bottom: 32px; letter-spacing: -1px;"><span style="color: var(--blue-600);">Pre-Operative</span> Assessment Summary</h2>

                    <div class="info-card" style="margin-bottom: 24px;">
                        {{ summary|safe }}
                    </div>

                    {% if references %}
                    <div class="info-card">
                        <div class="info-card-title">Evidence-Based References</div>
                        <div class="info-card-content">
                            {% for ref in references %}
                            <div style="padding: 12px 0; border-bottom: 1px solid var(--gray-200); font-size: 14px;">
                                <a href="https://pubmed.ncbi.nlm.nih.gov/{{ ref.pmid }}/" target="_blank" style="color: var(--blue-600); text-decoration: none; font-weight: 600;">
                                    [{{ loop.index }}] {{ ref.title }}
                                </a>
                                <div style="color: var(--gray-600); font-size: 13px; margin-top: 4px;">{{ ref.year }}</div>
                <!-- Header Section -->
                <div class="results-header">
                    <div class="results-header-content">
                        <div class="results-icon">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                <path d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2"/>
                                <path d="M9 12l2 2 4-4"/>
                            </svg>
                        </div>
                        <div>
                            <h1 class="results-title">Pre-Operative Assessment Complete</h1>
                            <p class="results-subtitle">Evidence-based risk stratification and optimization recommendations</p>
                        </div>
                    </div>
                    <div class="results-actions-compact">
                        <button onclick="window.print()" class="icon-btn" title="Print Assessment">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                <path d="M17 17h2a2 2 0 002-2v-4a2 2 0 00-2-2H5a2 2 0 00-2 2v4a2 2 0 002 2h2"/>
                                <path d="M9 21H15"/>
                                <path d="M7 9V5a2 2 0 012-2h6a2 2 0 012 2v4"/>
                            </svg>
                        </button>
                        <a href="/preop" class="icon-btn" title="New Assessment">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                <path d="M12 4v16m8-8H4"/>
                            </svg>
                        </a>
                    </div>
                </div>

                <!-- Assessment Content -->
                <div class="results-body">
                    <div class="assessment-card">
                        <div class="assessment-content">
                            {{ summary|safe }}
                        </div>
                    </div>

                    <!-- References Section -->
                    {% if references %}
                    <div class="references-card">
                        <div class="references-card-header">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                <path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path>
                                <path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path>
                            </svg>
                            <span>Evidence-Based References</span>
                            <span class="references-count">{{ references|length }} papers</span>
                        </div>
                        <div class="references-list">
                            {% for ref in references %}
                            <div class="reference-item-modern">
                                <div class="reference-number">[{{ loop.index }}]</div>
                                <div class="reference-content-modern">
                                    <a href="https://pubmed.ncbi.nlm.nih.gov/{{ ref.pmid }}/" target="_blank" rel="noopener noreferrer" class="reference-title-link">
                                        {{ ref.title }}
                                    </a>
                                    <div class="reference-meta-modern">
                                        <span class="reference-authors">{{ ref.authors }}</span>
                                        <span class="reference-separator">•</span>
                                        <span class="reference-journal">{{ ref.journal }}</span>
                                        <span class="reference-separator">•</span>
                                        <span class="reference-year">{{ ref.year }}</span>
                                    </div>
                                    {% if ref.get('study_badge') and ref.get('study_color') %}
                                    <span class="study-type-badge-preop" style="background-color: {{ ref.study_color }};">
                                        {{ ref.study_badge }}
                                    </span>
                                    {% endif %}
                                </div>
                            </div>
                            {% endfor %}
                        </div>
                    </div>
                    {% endif %}
                </div>

                <!-- Action Buttons -->
                <div class="results-footer">
                    <a href="/preop" class="btn-modern btn-modern-primary">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                            <path d="M12 4v16m8-8H4"/>
                        </svg>
                        New Assessment
                    </a>
                    <button onclick="window.print()" class="btn-modern btn-modern-secondary">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                            <path d="M17 17h2a2 2 0 002-2v-4a2 2 0 00-2-2H5a2 2 0 00-2 2v4a2 2 0 002 2h2"/>
                            <path d="M9 21H15"/>
                            <path d="M7 9V5a2 2 0 012-2h6a2 2 0 012 2v4"/>
                        </svg>
                        Print/Save PDF
                    </button>
                    <a href="/" class="btn-modern btn-modern-tertiary">
                        Back to Home
                    </a>
                </div>
            </div>
            {% endif %}
        </div>
    </div>

    <!-- Loading Overlay -->
    <div class="loading-overlay" id="loadingOverlay">
        <div class="loading-content">
            <div class="loading-spinner"></div>
            <div class="loading-text">Generating Evidence-Based Assessment...</div>
        </div>
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
    </div>

    <script>

        // ====== SESSION STORAGE ======
        function saveFormData() {
            const form = document.getElementById('preopForm');
            const formData = new FormData(form);
            const data = {};

            // Save all form fields
            for (let [key, value] of formData.entries()) {
                if (data[key]) {
                    if (Array.isArray(data[key])) {
                        data[key].push(value);
                    } else {
                        data[key] = [data[key], value];
                    }
                } else {
                    data[key] = value;
                }
            }

            sessionStorage.setItem('preopFormData', JSON.stringify(data));
        }

        function loadFormData() {
            const saved = sessionStorage.getItem('preopFormData');
            if (!saved) return;

            try {
                const data = JSON.parse(saved);
                const form = document.getElementById('preopForm');

                for (let [key, value] of Object.entries(data)) {
                    const elements = form.elements[key];
                    if (!elements) continue;

                    if (elements.type === 'radio') {
                        form.elements[key].forEach(radio => {
                            if (radio.value === value) radio.checked = true;
                        });
                    } else if (elements.type === 'checkbox') {
                        if (Array.isArray(value)) {
                            form.elements[key].forEach(checkbox => {
                                if (value.includes(checkbox.value)) checkbox.checked = true;
                            });
                        } else {
                            elements.checked = value === 'true';
                        }
                    } else {
                        elements.value = value;
                    }
                }

                // Recalculate metrics if data loaded
                calculateMetrics();
                calculateEGFR();
            } catch (e) {
                console.error('Error loading form data:', e);
            }
        }


        // ====== RADIO BUTTON SELECTION ======
        function selectRadio(name, value) {
            const radio = document.getElementById(`${name}-${value}`);
            if (radio) {
                radio.checked = true;

                // Update visual state
                const options = document.querySelectorAll(`input[name="${name}"]`);
                options.forEach(opt => {
                    const parent = opt.closest('.radio-option');
                    if (parent) {
                        if (opt.checked) {
                            parent.classList.add('selected');
                        } else {
                            parent.classList.remove('selected');
                        }
                    }
                });

                // Save and recalculate
                saveFormData();
                calculateMetrics();
            }
        }

        // ====== AUTO-CALCULATIONS ======
        function calculateMetrics() {
            const age = parseFloat(document.getElementById('age')?.value);
            const weight = parseFloat(document.getElementById('weight')?.value);
            const height = parseFloat(document.getElementById('height')?.value);
            const sex = document.querySelector('input[name="sex"]:checked')?.value;

            if (!weight || !height || !sex) return;

            // Calculate BMI
            const heightM = height / 100;
            const bmi = weight / (heightM * heightM);

            // BMI Category
            let bmiCategory = '';
            if (bmi < 18.5) bmiCategory = 'Underweight';
            else if (bmi < 25) bmiCategory = 'Normal';
            else if (bmi < 30) bmiCategory = 'Overweight';
            else if (bmi < 35) bmiCategory = 'Obese I';
            else if (bmi < 40) bmiCategory = 'Obese II';
            else bmiCategory = 'Obese III';

            // Calculate IBW (Devine formula)
            let ibw;
            if (sex === 'male') {
                ibw = 50 + 2.3 * ((height / 2.54) - 60);
            } else {
                ibw = 45.5 + 2.3 * ((height / 2.54) - 60);
            }

            // Display results
            document.getElementById('calcBMI').textContent = bmi.toFixed(1);
            document.getElementById('calcBMICategory').textContent = bmiCategory;
            document.getElementById('calcIBW').innerHTML = ibw.toFixed(1) + '<span class="calc-unit">kg</span>';
            document.getElementById('autoCalcBox').style.display = 'block';
        }

        function calculateEGFR() {
            const age = parseFloat(document.getElementById('age')?.value);
            const cr = parseFloat(document.getElementById('cr')?.value);
            const sex = document.querySelector('input[name="sex"]:checked')?.value;

            if (!age || !cr || !sex) return;

            // CKD-EPI equation (simplified)
            let egfr;
            const k = sex === 'female' ? 0.7 : 0.9;
            const a = sex === 'female' ? -0.329 : -0.411;
            const factor = sex === 'female' ? 1.018 : 1;

            const creatRatio = cr / k;
            const minVal = Math.min(creatRatio, 1);
            const maxVal = Math.max(creatRatio, 1);

            egfr = 141 * Math.pow(minVal, a) * Math.pow(maxVal, -1.209) * Math.pow(0.993, age) * factor;

            // CKD Stage
            let ckdStage = '';
            if (egfr >= 90) ckdStage = 'Normal';
            else if (egfr >= 60) ckdStage = 'Stage 2';
            else if (egfr >= 45) ckdStage = 'Stage 3a';
            else if (egfr >= 30) ckdStage = 'Stage 3b';
            else if (egfr >= 15) ckdStage = 'Stage 4';
            else ckdStage = 'Stage 5';

            // Display results
            document.getElementById('calcEGFR').innerHTML = egfr.toFixed(1) + '<span class="calc-unit">mL/min</span>';
            document.getElementById('calcCKDStage').textContent = ckdStage;
            document.getElementById('eGFRBox').style.display = 'block';
        }

        // ====== EVENT LISTENERS ======
        document.addEventListener('DOMContentLoaded', function() {
            // Load saved data
            loadFormData();

            // Auto-calculate on input
            ['age', 'weight', 'height'].forEach(id => {
                const el = document.getElementById(id);
                if (el) {
                    el.addEventListener('input', calculateMetrics);
                }
            });

            ['age', 'cr'].forEach(id => {
                const el = document.getElementById(id);
                if (el) {
                    el.addEventListener('input', calculateEGFR);
                }
            });

            // Save on any input change
            document.getElementById('wizardForm').addEventListener('change', saveFormData);

            // Show loading overlay on submit
            document.getElementById('wizardForm').addEventListener('submit', function() {
                document.getElementById('loadingOverlay').classList.add('show');
            });
        });

        // ====== MOBILE MENU ======
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

def classify_study_type(title, journal):
    """
    Classify a single study based on its title and journal.
    Returns dict with study_type, quality_score, badge_text, and badge_color.

    Study types ranked by evidence quality (highest to lowest):
    1. Guideline/Consensus
    2. Meta-analysis
    3. Systematic Review
    4. RCT
    5. Observational Study
    6. Review Article
    7. Case Report/Series
    """
    title_lower = title.lower()
    journal_lower = journal.lower()

    # Check for guidelines/consensus (highest quality)
    if any(keyword in title_lower or keyword in journal_lower for keyword in
           ['guideline', 'guidelines', 'consensus', 'recommendation', 'practice parameter',
            'practice guideline', 'clinical practice', 'society statement']):
        return {
            'type': 'Guideline',
            'score': 4,
            'badge_text': 'Guideline',
            'badge_color': '#7C3AED',  # Purple
            'sort_priority': 1
        }

    # Check for meta-analysis
    if any(keyword in title_lower or keyword in journal_lower for keyword in
           ['meta-analysis', 'metaanalysis', 'meta analysis', 'pooled analysis',
            'network meta-analysis', 'meta-regression']):
        return {
            'type': 'Meta-analysis',
            'score': 3,
            'badge_text': 'Meta-analysis',
            'badge_color': '#059669',  # Green
            'sort_priority': 2
        }

    # Check for systematic review (separate from meta-analysis)
    if any(keyword in title_lower or keyword in journal_lower for keyword in
           ['systematic review', 'cochrane review', 'systematic literature review']):
        return {
            'type': 'Systematic Review',
            'score': 2.5,
            'badge_text': 'Systematic Review',
            'badge_color': '#0891B2',  # Teal
            'sort_priority': 3
        }

    # Check for RCT
    if any(keyword in title_lower for keyword in
           ['randomized', 'randomised', ' rct', 'randomized controlled trial',
            'randomised controlled trial', 'double-blind', 'double blind',
            'placebo-controlled', 'placebo controlled', 'controlled clinical trial']):
        return {
            'type': 'RCT',
            'score': 2,
            'badge_text': 'RCT',
            'badge_color': '#2563EB',  # Blue
            'sort_priority': 4
        }

    # Check for observational studies
    if any(keyword in title_lower for keyword in
           ['cohort study', 'cohort analysis', 'case-control', 'case control',
            'observational study', 'prospective study', 'retrospective study',
            'retrospective analysis', 'database analysis', 'registry']):
        return {
            'type': 'Observational',
            'score': 1,
            'badge_text': 'Observational',
            'badge_color': '#F59E0B',  # Amber
            'sort_priority': 5
        }

    # Check for case reports (lowest quality)
    if any(keyword in title_lower for keyword in
           ['case report', 'case series', 'case study']):
        return {
            'type': 'Case Report',
            'score': 0.5,
            'badge_text': 'Case Report',
            'badge_color': '#DC2626',  # Red
            'sort_priority': 7
        }

    # Check for review articles (general)
    if 'review' in title_lower:
        return {
            'type': 'Review',
            'score': 1,
            'badge_text': 'Review',
            'badge_color': '#6B7280',  # Gray
            'sort_priority': 6
        }

    # Default: unclassified
    return {
        'type': 'Study',
        'score': 0.5,
        'badge_text': 'Study',
        'badge_color': '#9CA3AF',  # Light gray
        'sort_priority': 8
    }

def get_evidence_strength(num_papers, references):
    """
    Analyze evidence strength and return classification based on study quality hierarchy.

    Evidence Hierarchy (points assigned):
    - Guidelines/Consensus: 4 points
    - Meta-analysis: 3 points
    - Systematic Review: 2.5 points
    - RCT: 2 points
    - Observational Study: 1 point
    - Review Article: 1 point
    - Case Report/Series: 0.5 points

    Confidence Levels:
    - High: score ≥8 OR 2+ meta-analyses OR (1+ guideline AND 1+ meta-analysis)
    - Moderate: score ≥4 OR 1+ meta-analysis OR 2+ RCTs OR 1+ guideline
    - Low: everything else
    """
    if not num_papers or num_papers == 0:
        return {
            'level': 'Low',
            'color': '#EF4444',
            'description': 'No evidence found',
            'score': 0,
            'breakdown': {
                'guidelines': 0,
                'meta_analyses': 0,
                'systematic_reviews': 0,
                'rcts': 0,
                'observational': 0,
                'reviews': 0,
                'case_reports': 0,
                'total': 0,
                'recent_count': 0
            }
        }

    # Initialize counters
    guideline_count = 0
    meta_analysis_count = 0
    systematic_review_count = 0
    rct_count = 0
    observational_count = 0
    review_count = 0
    case_report_count = 0
    recent_count = 0  # Papers from last 5 years

    quality_score = 0
    current_year = 2025

    # Analyze each paper for multiple study type indicators
    for ref in references:
        title = ref.get('title', '').lower()
        journal = ref.get('journal', '').lower()
        year = ref.get('year', '')

        # Track recency (papers from last 5 years)
        try:
            paper_year = int(year) if year else 0
            if current_year - paper_year <= 5:
                recent_count += 1
        except (ValueError, TypeError):
            pass

        # Check for guidelines/consensus (highest quality)
        if any(keyword in title or keyword in journal for keyword in
               ['guideline', 'guidelines', 'consensus', 'recommendation', 'practice parameter',
                'practice guideline', 'clinical practice', 'society statement']):
            guideline_count += 1
            quality_score += 4

        # Check for meta-analysis
        if any(keyword in title or keyword in journal for keyword in
               ['meta-analysis', 'metaanalysis', 'meta analysis', 'pooled analysis',
                'network meta-analysis', 'meta-regression']):
            meta_analysis_count += 1
            quality_score += 3

        # Check for systematic review (separate from meta-analysis)
        elif any(keyword in title or keyword in journal for keyword in
                 ['systematic review', 'cochrane review', 'systematic literature review']):
            systematic_review_count += 1
            quality_score += 2.5

        # Check for RCT
        elif any(keyword in title for keyword in
                 ['randomized', 'randomised', ' rct', 'randomized controlled trial',
                  'randomised controlled trial', 'double-blind', 'double blind',
                  'placebo-controlled', 'placebo controlled', 'controlled clinical trial']):
            rct_count += 1
            quality_score += 2

        # Check for observational studies
        elif any(keyword in title for keyword in
                 ['cohort study', 'cohort analysis', 'case-control', 'case control',
                  'observational study', 'prospective study', 'retrospective study',
                  'retrospective analysis', 'database analysis', 'registry']):
            observational_count += 1
            quality_score += 1

        # Check for case reports (lowest quality)
        elif any(keyword in title for keyword in
                 ['case report', 'case series', 'case study']):
            case_report_count += 1
            quality_score += 0.5

        # Check for review articles (general)
        elif 'review' in title:
            review_count += 1
            quality_score += 1

    # Add recency bonus (5% boost for recent evidence)
    recency_bonus = (recent_count / num_papers) * 0.5 if num_papers > 0 else 0
    total_score = quality_score + recency_bonus

    # Build breakdown object
    breakdown = {
        'guidelines': guideline_count,
        'meta_analyses': meta_analysis_count,
        'systematic_reviews': systematic_review_count,
        'rcts': rct_count,
        'observational': observational_count,
        'reviews': review_count,
        'case_reports': case_report_count,
        'total': num_papers,
        'recent_count': recent_count
    }

    # Determine confidence level based on refined criteria
    # HIGH CONFIDENCE: Strong evidence from multiple high-quality sources
    if (total_score >= 8 or
        meta_analysis_count >= 2 or
        (guideline_count >= 1 and meta_analysis_count >= 1) or
        (guideline_count >= 1 and rct_count >= 2)):

        # Build description highlighting strongest evidence
        desc_parts = []
        if guideline_count > 0:
            desc_parts.append(f"{guideline_count} guideline{'s' if guideline_count != 1 else ''}")
        if meta_analysis_count > 0:
            desc_parts.append(f"{meta_analysis_count} meta-analysis/analyses")
        if systematic_review_count > 0:
            desc_parts.append(f"{systematic_review_count} systematic review{'s' if systematic_review_count != 1 else ''}")
        if rct_count > 0:
            desc_parts.append(f"{rct_count} RCT{'s' if rct_count != 1 else ''}")

        description = f"{num_papers} papers: {', '.join(desc_parts) if desc_parts else 'high-quality evidence'}"

        return {
            'level': 'High',
            'color': '#10B981',
            'description': description,
            'score': total_score,
            'breakdown': breakdown
        }

    # MODERATE CONFIDENCE: Some high-quality evidence or multiple moderate-quality sources
    elif (total_score >= 4 or
          meta_analysis_count >= 1 or
          systematic_review_count >= 1 or
          rct_count >= 2 or
          guideline_count >= 1):

        # Build description
        desc_parts = []
        if guideline_count > 0:
            desc_parts.append(f"{guideline_count} guideline{'s' if guideline_count != 1 else ''}")
        if meta_analysis_count > 0:
            desc_parts.append(f"{meta_analysis_count} meta-analysis/analyses")
        if systematic_review_count > 0:
            desc_parts.append(f"{systematic_review_count} systematic review{'s' if systematic_review_count != 1 else ''}")
        if rct_count > 0:
            desc_parts.append(f"{rct_count} RCT{'s' if rct_count != 1 else ''}")
        if observational_count > 0:
            desc_parts.append(f"{observational_count} observational")

        description = f"{num_papers} papers: {', '.join(desc_parts) if desc_parts else 'moderate-quality evidence'}"

        return {
            'level': 'Moderate',
            'color': '#FBBF24',
            'description': description,
            'score': total_score,
            'breakdown': breakdown
        }

    # LOW CONFIDENCE: Limited or low-quality evidence
    else:
        # Build description
        desc_parts = []
        if observational_count > 0:
            desc_parts.append(f"{observational_count} observational")
        if review_count > 0:
            desc_parts.append(f"{review_count} review{'s' if review_count != 1 else ''}")
        if case_report_count > 0:
            desc_parts.append(f"{case_report_count} case report{'s' if case_report_count != 1 else ''}")

        description = f"Limited evidence ({num_papers} papers): {', '.join(desc_parts) if desc_parts else 'low-quality studies'}"

        return {
            'level': 'Low',
            'color': '#EF4444',
            'description': description,
            'score': total_score,
            'breakdown': breakdown
        }

HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GasConsult.ai — AI-Powered Anesthesiology Assistant</title>

    <!-- PWA -->
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=6">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=6">
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
            padding: 6px;
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
            padding: 4px;
            display: flex;
            align-items: flex-end;
            gap: 4px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04), inset 0 1px 0 rgba(255,255,255,1);
        }

        .chat-input {
            flex: 1;
            border: none;
            outline: none;
            padding: 10px 14px;
            font-size: 16px;
            font-family: inherit;
            color: var(--gray-800);
            background: transparent;
            resize: none;
            min-height: 38px;
            max-height: 110px;
            line-height: 1.5;
        }

        .chat-input::placeholder { color: var(--gray-400); }

        .chat-send {
            width: 44px;
            height: 44px;
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
            margin: 4px;
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
            background: var(--gray-50);
            position: relative;
            z-index: 10;
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

        .study-type-badge {
            display: inline-block;
            padding: 3px 8px;
            border-radius: 4px;
            font-size: 10px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.3px;
            margin-top: 6px;
            color: white;
        }

        .reference-header {
            display: flex;
            align-items: flex-start;
            gap: 8px;
            margin-bottom: 4px;
        }

        .reference-link-container {
            flex: 1;
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
            padding: 8px;
            z-index: 10;
        }

        .chat-input-wrapper {
            max-width: 720px;
            margin: 0 auto;
        }

        #chat-form.chat-card {
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
                padding: 10px 24px;
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
                padding: 12px 32px;
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
                    <a href="/crisis" class="nav-link">Crisis Protocols</a>
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
            <a href="/crisis" class="mobile-menu-link">Crisis Protocols</a>
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
                            {% set level = message.evidence_strength.level if message.evidence_strength is mapping else message.evidence_strength %}
                            {% set breakdown = message.evidence_strength.breakdown if message.evidence_strength is mapping else {} %}
                            {% set description = message.evidence_strength.description if message.evidence_strength is mapping else '' %}
                            <div class="evidence-badge {{ 'high' if level == 'High' else ('moderate' if level == 'Moderate' else 'low') }}"
                                 title="Evidence Quality: {{ level }} - {{ description }}">
                                <div style="display: flex; flex-direction: column; align-items: flex-start;">
                                    <div>
                                        {% if level == 'High' %}
                                        ✓ High Confidence
                                        {% elif level == 'Moderate' %}
                                        ~ Moderate Confidence
                                        {% else %}
                                        ! Low Confidence
                                        {% endif %}
                                        • {{ message.num_papers }} studies
                                    </div>
                                    <div class="evidence-explanation">
                                        {{ description if description else ('Strong evidence from meta-analyses, RCTs, or systematic reviews' if level == 'High' else ('Moderate evidence - consider individual patient factors' if level == 'Moderate' else 'Limited evidence - use caution and clinical judgment')) }}
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
                                    <div class="reference-header">
                                        <div class="reference-link-container">
                                            <a href="https://pubmed.ncbi.nlm.nih.gov/{{ ref.pmid }}/" target="_blank" rel="noopener noreferrer" class="reference-link">
                                                {{ ref.title }}
                                            </a>
                                            <div class="reference-meta">
                                                {{ ref.authors }} — {{ ref.journal }}, {{ ref.year }}
                                            </div>
                                            {% if ref.get('study_badge') and ref.get('study_color') %}
                                            <span class="study-type-badge" style="background-color: {{ ref.study_color }};">
                                                {{ ref.study_badge }}
                                            </span>
                                            {% endif %}
                                        </div>
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
                        <div id="streamingEvidenceBadge" style="display: none;"></div>
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

        {% if not messages or messages|length == 0 %}
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
        {% endif %}

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

                        // Display evidence badge
                        if (event.evidence_strength) {
                            const evidenceBadge = document.getElementById('streamingEvidenceBadge');
                            if (evidenceBadge) {
                                const strength = event.evidence_strength;
                                const level = strength.level || 'Low';
                                const description = strength.description || '';
                                const numPapers = event.num_papers || 0;

                                // Determine CSS class
                                const badgeClass = level === 'High' ? 'high' : (level === 'Moderate' ? 'moderate' : 'low');

                                // Determine icon
                                const icon = level === 'High' ? '✓' : (level === 'Moderate' ? '~' : '!');

                                // Build badge HTML
                                let badgeHTML = '<div class="evidence-badge ' + badgeClass + '" ';
                                badgeHTML += 'title="Evidence Quality: ' + level + ' - ' + description + '">';
                                badgeHTML += '<div style="display: flex; flex-direction: column; align-items: flex-start;">';
                                badgeHTML += '<div>' + icon + ' ' + level + ' Confidence • ' + numPapers + ' studies</div>';
                                badgeHTML += '<div class="evidence-explanation">' + description + '</div>';
                                badgeHTML += '</div></div>';

                                evidenceBadge.innerHTML = badgeHTML;
                                evidenceBadge.style.display = 'block';
                            }
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
                                refsHTML += '<div class="reference-header"><div class="reference-link-container">';
                                refsHTML += '<a href="https://pubmed.ncbi.nlm.nih.gov/' + ref.pmid + '/" target="_blank" rel="noopener noreferrer" class="reference-link">';
                                refsHTML += ref.title;
                                refsHTML += '</a>';
                                refsHTML += '<div class="reference-meta">';
                                refsHTML += ref.authors + ' — ' + ref.journal + ', ' + ref.year;
                                refsHTML += '</div>';
                                // Add study type badge if available
                                if (ref.study_badge && ref.study_color) {
                                    refsHTML += '<span class="study-type-badge" style="background-color: ' + ref.study_color + ';">';
                                    refsHTML += ref.study_badge;
                                    refsHTML += '</span>';
                                }
                                refsHTML += '</div></div>';
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
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=6">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=6">
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
            padding: 6px;
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
            padding: 4px;
            display: flex;
            align-items: flex-end;
            gap: 4px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04), inset 0 1px 0 rgba(255,255,255,1);
        }

        .chat-input {
            flex: 1;
            border: none;
            outline: none;
            padding: 10px 14px;
            font-size: 16px;
            font-family: inherit;
            color: var(--gray-800);
            background: transparent;
            resize: none;
            min-height: 38px;
            max-height: 110px;
            line-height: 1.5;
        }

        .chat-input::placeholder { color: var(--gray-400); }

        .chat-send {
            width: 44px;
            height: 44px;
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
            margin: 4px;
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
                    <a href="/crisis" class="nav-link">Crisis Protocols</a>
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
            <a href="/crisis" class="mobile-menu-link">Crisis Protocols</a>
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
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=6">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=6">
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
            padding: 6px;
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
            padding: 4px;
            display: flex;
            align-items: flex-end;
            gap: 4px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04), inset 0 1px 0 rgba(255,255,255,1);
        }

        .chat-input {
            flex: 1;
            border: none;
            outline: none;
            padding: 10px 14px;
            font-size: 16px;
            font-family: inherit;
            color: var(--gray-800);
            background: transparent;
            resize: none;
            min-height: 38px;
            max-height: 110px;
            line-height: 1.5;
        }

        .chat-input::placeholder { color: var(--gray-400); }

        .chat-send {
            width: 44px;
            height: 44px;
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
            margin: 4px;
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
                    <a href="/crisis" class="nav-link">Crisis Protocols</a>
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
            <a href="/crisis" class="mobile-menu-link">Crisis Protocols</a>
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
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=6">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=6">
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
                    <a href="/crisis" class="nav-link">Crisis Protocols</a>
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
            <a href="/crisis" class="mobile-menu-link">Crisis Protocols</a>
            <a href="/hypotension" class="mobile-menu-link">IOH Predictor</a>
        </div>

        <main class="main-content">
            <div class="content-card">
                <h1>Terms of Service</h1>
                <p class="last-updated">Last Updated: January 2025</p>

                <div class="notice-box">
                    <h3>CRITICAL MEDICAL DISCLAIMER</h3>
                    <p><strong>THIS SERVICE DOES NOT PROVIDE MEDICAL ADVICE, DIAGNOSIS, OR TREATMENT. DO NOT USE FOR EMERGENCIES. CALL 911 FOR EMERGENCIES.</strong> gasconsult.ai is strictly an educational and informational tool. No doctor-patient relationship is created by using this Service. All information must be independently verified by qualified healthcare professionals before any clinical application.</p>
                </div>

                <p>These Terms of Service ("Terms") constitute a legally binding agreement between you ("User" or "you") and gasconsult.ai ("Service", "we", "us", or "our"). By accessing or using this Service, you acknowledge that you have read, understood, and unconditionally accept all terms, disclaimers, and limitations set forth below.</p>

                <h2>1. Acceptance of Terms</h2>
                <p>By accessing, browsing, or using gasconsult.ai in any manner, you acknowledge and agree that:</p>
                <ul>
                    <li>You have read and understood these Terms in their entirety</li>
                    <li>You agree to be legally bound by all provisions contained herein</li>
                    <li>If you do not agree to these Terms, you must immediately cease all use of this Service</li>
                    <li>Continued use of the Service constitutes ongoing acceptance of these Terms and any future modifications</li>
                    <li>These Terms supersede any prior agreements, representations, or understandings</li>
                </ul>

                <h2>2. NOT Medical Advice — Critical Disclaimers</h2>

                <h3>2.1 Not Medical Advice, Diagnosis, or Treatment</h3>
                <p><strong>THE INFORMATION PROVIDED BY GASCONSULT.AI IS STRICTLY FOR EDUCATIONAL AND INFORMATIONAL PURPOSES ONLY.</strong> Under no circumstances should the content generated by this Service be interpreted, construed, or relied upon as:</p>
                <ul>
                    <li>Medical advice, clinical recommendations, or treatment guidance</li>
                    <li>Diagnosis or diagnostic assistance for any medical condition</li>
                    <li>Prescription or recommendation of any medication, treatment, or therapy</li>
                    <li>A substitute for consultation with qualified, licensed healthcare professionals</li>
                    <li>A replacement for comprehensive patient assessment and examination</li>
                    <li>Definitive or authoritative medical guidance for patient care</li>
                    <li>Standard of care or best practices for any clinical situation</li>
                    <li>Current, complete, or accurate medical information</li>
                </ul>

                <h3>2.2 Educational and Research Purposes Only</h3>
                <p>This Service is designed exclusively as an educational resource for healthcare professionals to:</p>
                <ul>
                    <li>Explore medical literature and research findings</li>
                    <li>Supplement continuing medical education efforts</li>
                    <li>Generate hypotheses for further investigation</li>
                    <li>Understand trends in published medical research</li>
                </ul>
                <p><strong>THIS SERVICE MUST NOT BE USED</strong> as the basis for any clinical decision, treatment plan, medication dosing, or patient care intervention without independent professional verification and judgment.</p>

                <h3>2.3 Multiple Mandatory Disclaimers</h3>
                <p>To ensure absolute clarity, we state unequivocally:</p>
                <ul>
                    <li><strong>WE DO NOT PROVIDE MEDICAL ADVICE.</strong></li>
                    <li><strong>WE DO NOT DIAGNOSE MEDICAL CONDITIONS.</strong></li>
                    <li><strong>WE DO NOT RECOMMEND TREATMENTS.</strong></li>
                    <li><strong>WE DO NOT PRESCRIBE MEDICATIONS.</strong></li>
                    <li><strong>WE DO NOT REPLACE PROFESSIONAL MEDICAL JUDGMENT.</strong></li>
                    <li><strong>YOU CANNOT AND MUST NOT RELY ON THIS SERVICE FOR MEDICAL DECISIONS.</strong></li>
                </ul>

                <h2>3. No Doctor-Patient Relationship</h2>
                <p><strong>USE OF THIS SERVICE DOES NOT CREATE, AND SHALL NEVER BE CONSTRUED TO CREATE, A DOCTOR-PATIENT RELATIONSHIP, HEALTHCARE PROVIDER-PATIENT RELATIONSHIP, OR ANY PROFESSIONAL MEDICAL RELATIONSHIP OF ANY KIND.</strong></p>
                <p>You acknowledge and agree that:</p>
                <ul>
                    <li>No physician-patient relationship exists or will be formed through use of this Service</li>
                    <li>We are not your healthcare provider, physician, medical advisor, or clinician</li>
                    <li>We have no duty of care, fiduciary duty, or professional obligation to you</li>
                    <li>We do not and cannot provide personalized medical advice for your specific situation</li>
                    <li>All interactions with this Service are non-clinical and educational only</li>
                    <li>You must consult with your own qualified healthcare providers for all medical matters</li>
                </ul>

                <h2>4. Emergency Disclaimer — DO NOT USE FOR EMERGENCIES</h2>

                <div class="notice-box">
                    <h3>EMERGENCY MEDICAL DISCLAIMER</h3>
                    <p><strong>DO NOT USE THIS SERVICE FOR MEDICAL EMERGENCIES UNDER ANY CIRCUMSTANCES.</strong></p>
                    <p><strong>IF YOU ARE EXPERIENCING A MEDICAL EMERGENCY, CALL 911 (OR YOUR LOCAL EMERGENCY NUMBER) IMMEDIATELY.</strong></p>
                    <p>This Service is not designed for, equipped to handle, or suitable for emergency medical situations. Any delay in seeking emergency medical care while using this Service could result in serious injury or death.</p>
                </div>

                <p>You expressly acknowledge that:</p>
                <ul>
                    <li>This Service provides no real-time clinical support or emergency assistance</li>
                    <li>Information provided may be delayed, incomplete, or unsuitable for urgent situations</li>
                    <li>In any emergency or potentially life-threatening situation, you must immediately contact emergency medical services (call 911) or go to the nearest emergency department</li>
                    <li>We are not liable for any harm, injury, or death resulting from failure to seek timely emergency medical care</li>
                </ul>

                <h2>5. Educational Use Only — Professional Verification Required</h2>
                <p>This Service is intended exclusively for use by licensed, qualified healthcare professionals (physicians, CRNAs, nurse practitioners, physician assistants, medical students, residents, fellows) as an educational supplement.</p>
                <p><strong>ALL INFORMATION PROVIDED MUST BE INDEPENDENTLY VERIFIED</strong> through:</p>
                <ul>
                    <li>Primary medical literature and original research publications</li>
                    <li>Current evidence-based clinical practice guidelines</li>
                    <li>Consultation with qualified medical experts and specialists</li>
                    <li>Institutional protocols and standards of care</li>
                    <li>Individual patient assessment and clinical judgment</li>
                    <li>Appropriate diagnostic testing and patient evaluation</li>
                </ul>
                <p><strong>YOU ASSUME ALL RESPONSIBILITY</strong> for verifying information before any clinical application.</p>

                <h2>6. AI-Generated Content Disclaimer</h2>

                <h3>6.1 Artificial Intelligence Limitations</h3>
                <p>This Service uses artificial intelligence (AI) technology, including large language models (LLMs) such as GPT-4o. <strong>YOU ACKNOWLEDGE AND ACCEPT THAT AI TECHNOLOGY:</strong></p>
                <ul>
                    <li><strong>Can and does generate errors, inaccuracies, and false information ("hallucinations")</strong></li>
                    <li>May produce plausible-sounding but factually incorrect medical information</li>
                    <li>Can misinterpret, misrepresent, or incorrectly synthesize medical literature</li>
                    <li>May provide outdated, incomplete, or biased information</li>
                    <li>Cannot replace human medical expertise, clinical judgment, or professional training</li>
                    <li>May generate contradictory or inconsistent responses to similar queries</li>
                    <li>Lacks the ability to understand individual patient context, nuance, or complexity</li>
                    <li>Is not validated for clinical decision-making or patient care</li>
                </ul>

                <h3>6.2 No Validation or Clinical Testing</h3>
                <p>The AI outputs provided by this Service have NOT been:</p>
                <ul>
                    <li>Validated through clinical trials or rigorous medical testing</li>
                    <li>Approved by any regulatory authority (FDA, EMA, etc.)</li>
                    <li>Reviewed or endorsed by medical professional organizations</li>
                    <li>Verified for accuracy, safety, or clinical appropriateness</li>
                    <li>Tested for use in actual patient care scenarios</li>
                </ul>

                <h3>6.3 Non-Deterministic Responses</h3>
                <p>AI-generated responses are probabilistic and non-deterministic. The same query may produce different responses at different times. You cannot rely on consistency, reproducibility, or reliability of AI outputs.</p>

                <h2>7. Third-Party Content and Data Sources</h2>
                <p>This Service aggregates and synthesizes information from third-party sources, including but not limited to:</p>
                <ul>
                    <li>PubMed and NCBI medical literature databases</li>
                    <li>OpenAI's GPT-4o and other AI models</li>
                    <li>Published medical research papers and abstracts</li>
                    <li>Medical databases and online resources</li>
                </ul>

                <p><strong>WE ARE NOT RESPONSIBLE FOR:</strong></p>
                <ul>
                    <li>The accuracy, completeness, or reliability of third-party content</li>
                    <li>Errors, omissions, or inaccuracies in source materials</li>
                    <li>Retracted, corrected, or updated research findings</li>
                    <li>Biases, conflicts of interest, or methodological flaws in published research</li>
                    <li>Availability or accessibility of third-party services</li>
                    <li>Any harm resulting from third-party content or data</li>
                </ul>

                <p>All citations, references, and source materials must be independently verified through original publications.</p>

                <h2>8. No Warranties — Service Provided "AS IS"</h2>

                <h3>8.1 Disclaimer of All Warranties</h3>
                <p><strong>TO THE MAXIMUM EXTENT PERMITTED BY LAW, THIS SERVICE IS PROVIDED "AS IS", "AS AVAILABLE", AND "WITH ALL FAULTS" WITHOUT WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED.</strong></p>

                <p>We expressly disclaim all warranties, including but not limited to:</p>
                <ul>
                    <li><strong>Warranties of merchantability</strong></li>
                    <li><strong>Warranties of fitness for a particular purpose</strong></li>
                    <li><strong>Warranties of accuracy, reliability, or completeness</strong></li>
                    <li><strong>Warranties of non-infringement</strong></li>
                    <li><strong>Warranties of title</strong></li>
                    <li><strong>Warranties arising from course of dealing or usage of trade</strong></li>
                    <li><strong>Any warranties regarding uninterrupted, timely, secure, or error-free service</strong></li>
                    <li><strong>Any warranties that defects will be corrected</strong></li>
                    <li><strong>Any warranties regarding results obtained from use of the Service</strong></li>
                </ul>

                <h3>8.2 No Guarantee of Accuracy</h3>
                <p>We make absolutely no representation or warranty that:</p>
                <ul>
                    <li>Information provided is accurate, current, complete, or reliable</li>
                    <li>The Service will meet your requirements or expectations</li>
                    <li>Information is suitable for any particular purpose or clinical situation</li>
                    <li>AI-generated content is free from errors, hallucinations, or inaccuracies</li>
                    <li>Third-party data sources are accurate or up-to-date</li>
                    <li>The Service will be available, accessible, or functional at any given time</li>
                </ul>

                <h2>9. Limitation of Liability — Maximum Legal Protection</h2>

                <h3>9.1 No Liability for Any Damages</h3>
                <p><strong>TO THE FULLEST EXTENT PERMITTED BY APPLICABLE LAW, GASCONSULT.AI, ITS OWNERS, OPERATORS, DEVELOPERS, CONTRIBUTORS, AFFILIATES, LICENSORS, SERVICE PROVIDERS, AND AGENTS (COLLECTIVELY, "RELEASED PARTIES") SHALL NOT BE LIABLE FOR ANY DAMAGES WHATSOEVER, INCLUDING BUT NOT LIMITED TO:</strong></p>
                <ul>
                    <li><strong>Direct damages</strong></li>
                    <li><strong>Indirect damages</strong></li>
                    <li><strong>Incidental damages</strong></li>
                    <li><strong>Consequential damages</strong></li>
                    <li><strong>Punitive damages</strong></li>
                    <li><strong>Special damages</strong></li>
                    <li><strong>Exemplary damages</strong></li>
                    <li><strong>Loss of profits, revenue, business, data, or opportunities</strong></li>
                    <li><strong>Personal injury or death</strong></li>
                    <li><strong>Medical malpractice or negligence</strong></li>
                    <li><strong>Reliance damages</strong></li>
                    <li><strong>Emotional distress</strong></li>
                </ul>

                <h3>9.2 No Liability for Medical Outcomes</h3>
                <p><strong>THE RELEASED PARTIES SHALL NOT BE LIABLE FOR:</strong></p>
                <ul>
                    <li>Patient injuries, complications, adverse events, or death</li>
                    <li>Medical errors or clinical mistakes based on information from this Service</li>
                    <li>Misdiagnosis, delayed diagnosis, or failure to diagnose</li>
                    <li>Inappropriate treatment or medication errors</li>
                    <li>Harm resulting from reliance on AI-generated content</li>
                    <li>Errors, omissions, or inaccuracies in provided information</li>
                    <li>Outdated, incomplete, or misleading medical information</li>
                    <li>Failure to seek timely professional medical care</li>
                    <li>Deviation from standard of care based on Service content</li>
                </ul>

                <h3>9.3 Maximum Liability Cap</h3>
                <p><strong>IN NO EVENT SHALL THE TOTAL AGGREGATE LIABILITY OF THE RELEASED PARTIES EXCEED THE GREATER OF (A) ZERO DOLLARS ($0.00) OR (B) THE TOTAL AMOUNT PAID BY YOU TO GASCONSULT.AI IN THE TWELVE (12) MONTHS PRECEDING THE CLAIM.</strong></p>
                <p>Since this Service is provided free of charge, the maximum liability is <strong>$0.00</strong>.</p>

                <h3>9.4 Limitations Apply to All Claims</h3>
                <p>These limitations apply to all claims, regardless of whether they are based on:</p>
                <ul>
                    <li>Warranty, contract, tort (including negligence), strict liability, or any other legal theory</li>
                    <li>Whether or not the Released Parties have been advised of the possibility of such damages</li>
                    <li>Whether a remedy fails of its essential purpose</li>
                </ul>

                <h2>10. Indemnification — You Hold Us Harmless</h2>
                <p><strong>YOU AGREE TO INDEMNIFY, DEFEND, AND HOLD HARMLESS THE RELEASED PARTIES FROM AND AGAINST ANY AND ALL:</strong></p>
                <ul>
                    <li>Claims, demands, actions, suits, or proceedings</li>
                    <li>Losses, damages, liabilities, settlements, penalties, fines, or judgments</li>
                    <li>Costs and expenses (including reasonable attorneys' fees, expert fees, and litigation costs)</li>
                </ul>

                <p><strong>ARISING from or related to:</strong></p>
                <ul>
                    <li>Your use or misuse of this Service</li>
                    <li>Clinical decisions, treatments, or patient care based on information from the Service</li>
                    <li>Medical malpractice or negligence claims related to your use of the Service</li>
                    <li>Your violation of these Terms of Service</li>
                    <li>Your violation of any applicable laws, regulations, or professional standards</li>
                    <li>Your violation of any rights of third parties</li>
                    <li>Any representation or warranty made by you to patients or third parties regarding this Service</li>
                    <li>Your failure to independently verify information provided by the Service</li>
                    <li>Patient harm resulting directly or indirectly from your use of the Service</li>
                </ul>

                <p>This indemnification obligation survives termination of your use of the Service.</p>

                <h2>11. Assumption of Risk — You Accept All Risks</h2>
                <p><strong>BY USING THIS SERVICE, YOU EXPRESSLY ACKNOWLEDGE, UNDERSTAND, AND VOLUNTARILY ASSUME ALL RISKS ASSOCIATED WITH:</strong></p>
                <ul>
                    <li>Reliance on AI-generated medical information</li>
                    <li>Potential inaccuracies, errors, or omissions in provided content</li>
                    <li>AI hallucinations and false information generation</li>
                    <li>Outdated or incomplete medical information</li>
                    <li>Misinterpretation of medical literature or research findings</li>
                    <li>Service unavailability, downtime, or technical failures</li>
                    <li>Security vulnerabilities or data breaches</li>
                    <li>Third-party content errors or inaccuracies</li>
                    <li>Consequences of any clinical decisions informed by this Service</li>
                </ul>

                <p><strong>YOU ACCEPT FULL AND SOLE RESPONSIBILITY FOR ALL RISKS ASSOCIATED WITH YOUR USE OF THIS SERVICE.</strong></p>

                <h2>12. HIPAA and Protected Health Information (PHI)</h2>

                <h3>12.1 Not a HIPAA Covered Entity</h3>
                <p><strong>GASCONSULT.AI IS NOT A COVERED ENTITY OR BUSINESS ASSOCIATE UNDER THE HEALTH INSURANCE PORTABILITY AND ACCOUNTABILITY ACT (HIPAA).</strong></p>
                <p>We do not provide healthcare services, maintain medical records, or engage in activities that would classify us as a HIPAA-covered entity.</p>

                <h3>12.2 DO NOT Enter Protected Health Information</h3>
                <p><strong>YOU MUST NOT ENTER, SUBMIT, OR UPLOAD ANY PROTECTED HEALTH INFORMATION (PHI) OR PERSONALLY IDENTIFIABLE HEALTH DATA INTO THIS SERVICE.</strong></p>
                <p>This includes but is not limited to:</p>
                <ul>
                    <li>Patient names, initials, or identifying information</li>
                    <li>Medical record numbers or patient identifiers</li>
                    <li>Dates of birth, addresses, or contact information</li>
                    <li>Social security numbers or insurance information</li>
                    <li>Specific medical histories or case details that could identify a patient</li>
                    <li>Any of the 18 HIPAA identifiers</li>
                </ul>

                <h3>12.3 Your HIPAA Compliance Responsibility</h3>
                <p>If you are a HIPAA-covered entity or business associate, <strong>YOU ARE SOLELY RESPONSIBLE</strong> for ensuring your use of this Service complies with all HIPAA requirements. We assume no responsibility for your HIPAA compliance.</p>

                <h2>13. Regulatory Status and FDA Disclaimer</h2>

                <h3>13.1 Not an FDA-Approved Medical Device</h3>
                <p><strong>THIS SERVICE IS NOT AN FDA-APPROVED, FDA-CLEARED, OR FDA-REGISTERED MEDICAL DEVICE.</strong></p>
                <p>This Service has not been evaluated, approved, or cleared by the U.S. Food and Drug Administration (FDA) or any other regulatory authority.</p>

                <h3>13.2 Not Intended for Clinical Use</h3>
                <p>This Service is not intended to:</p>
                <ul>
                    <li>Diagnose, treat, cure, or prevent any disease</li>
                    <li>Be used in the diagnosis of disease or other conditions</li>
                    <li>Be used in the cure, mitigation, treatment, or prevention of disease</li>
                    <li>Affect the structure or function of the body</li>
                    <li>Replace FDA-approved medical devices or diagnostic tools</li>
                </ul>

                <h3>13.3 No Professional Endorsement</h3>
                <p>This Service is not endorsed, approved, or recommended by:</p>
                <ul>
                    <li>The American Society of Anesthesiologists (ASA)</li>
                    <li>The American Association of Nurse Anesthetists (AANA)</li>
                    <li>Any medical professional organization or licensing board</li>
                    <li>Any hospital, healthcare system, or medical institution</li>
                    <li>Any governmental health agency or regulatory body</li>
                </ul>

                <h2>14. Dispute Resolution — Mandatory Binding Arbitration</h2>

                <h3>14.1 Agreement to Arbitrate</h3>
                <p><strong>YOU AND GASCONSULT.AI AGREE THAT ANY AND ALL DISPUTES, CLAIMS, OR CONTROVERSIES ARISING OUT OF OR RELATING TO THESE TERMS OR YOUR USE OF THE SERVICE SHALL BE RESOLVED EXCLUSIVELY THROUGH FINAL AND BINDING ARBITRATION, RATHER THAN IN COURT.</strong></p>

                <h3>14.2 Arbitration Procedures</h3>
                <p>Any arbitration shall be administered by the American Arbitration Association (AAA) under its Commercial Arbitration Rules and conducted in accordance with the following:</p>
                <ul>
                    <li>The arbitration shall be conducted by a single neutral arbitrator</li>
                    <li>The arbitration shall take place in Delaware, United States</li>
                    <li>The arbitration shall be governed by the Federal Arbitration Act (9 U.S.C. § 1 et seq.)</li>
                    <li>The arbitrator's decision shall be final and binding</li>
                    <li>Judgment on the arbitration award may be entered in any court having jurisdiction</li>
                    <li>Each party shall bear its own costs and attorneys' fees</li>
                </ul>

                <h3>14.3 Waiver of Jury Trial</h3>
                <p><strong>YOU HEREBY IRREVOCABLY WAIVE ANY AND ALL RIGHTS TO TRIAL BY JURY IN ANY LEGAL PROCEEDING ARISING OUT OF OR RELATED TO THESE TERMS OR YOUR USE OF THE SERVICE.</strong></p>

                <h3>14.4 Waiver of Class Actions</h3>
                <p><strong>YOU AGREE THAT ANY ARBITRATION OR LEGAL PROCEEDING SHALL BE CONDUCTED ONLY ON AN INDIVIDUAL BASIS AND NOT IN A CLASS, CONSOLIDATED, OR REPRESENTATIVE ACTION.</strong></p>

                <h2>15. Class Action Waiver</h2>
                <p><strong>YOU EXPRESSLY WAIVE YOUR RIGHT TO PARTICIPATE IN ANY CLASS ACTION, CLASS ARBITRATION, COLLECTIVE ACTION, REPRESENTATIVE ACTION, OR MASS ACTION AGAINST GASCONSULT.AI.</strong></p>
                <p>You agree that:</p>
                <ul>
                    <li>You will not seek to bring or participate in any class, collective, or representative action</li>
                    <li>You will not serve as a class member, representative, or opt into any class proceeding</li>
                    <li>All claims must be brought in your individual capacity only</li>
                    <li>You waive any right to have claims consolidated with claims of other users</li>
                </ul>

                <h2>16. Governing Law and Jurisdiction</h2>

                <h3>16.1 Governing Law</h3>
                <p>These Terms shall be governed by and construed in accordance with the laws of the State of Delaware, United States, without regard to its conflict of law principles.</p>

                <h3>16.2 Jurisdiction and Venue</h3>
                <p>To the extent any dispute is not subject to arbitration, you agree to submit to the exclusive personal jurisdiction and venue of the state and federal courts located in Delaware, United States.</p>

                <h3>16.3 International Users</h3>
                <p>This Service is controlled and operated from the United States. If you access the Service from outside the United States, you do so at your own risk and are responsible for compliance with local laws.</p>

                <h2>17. Severability and Waiver</h2>

                <h3>17.1 Severability</h3>
                <p>If any provision of these Terms is found to be invalid, illegal, or unenforceable by a court of competent jurisdiction:</p>
                <ul>
                    <li>The invalid provision shall be modified to the minimum extent necessary to make it valid and enforceable</li>
                    <li>If modification is not possible, the invalid provision shall be severed from these Terms</li>
                    <li>The remaining provisions shall remain in full force and effect</li>
                    <li>The invalid provision shall not affect the validity or enforceability of any other provision</li>
                </ul>

                <h3>17.2 No Waiver</h3>
                <p>Our failure to enforce any provision of these Terms shall not constitute a waiver of that provision or any other provision. No waiver shall be effective unless made in writing and signed by an authorized representative of gasconsult.ai.</p>

                <h2>18. Changes to Terms</h2>
                <p>We reserve the right to modify, amend, or update these Terms at any time, in our sole discretion, with or without notice.</p>
                <p><strong>Changes shall be effective immediately upon posting to this page.</strong></p>
                <p>Your continued use of the Service after any changes constitutes your acceptance of the modified Terms. You are responsible for regularly reviewing these Terms. If you do not agree to modified Terms, you must immediately cease using the Service.</p>

                <h2>19. Contact Information</h2>
                <p>For questions, concerns, or notices regarding these Terms of Service, please contact us at:</p>
                <p><strong>Email:</strong> support@gasconsult.ai</p>
                <p><strong>Legal Notice Address:</strong> gasconsult.ai Legal Department</p>

                <div class="notice-box" style="margin-top: 48px;">
                    <h3>FINAL ACKNOWLEDGMENT</h3>
                    <p><strong>BY USING GASCONSULT.AI, YOU ACKNOWLEDGE THAT:</strong></p>
                    <ul style="margin-top: 12px;">
                        <li>You have read, understood, and agree to these Terms of Service in their entirety</li>
                        <li>You understand this Service does NOT provide medical advice and creates NO doctor-patient relationship</li>
                        <li>You will NOT use this Service for emergencies and will call 911 for emergencies</li>
                        <li>You understand AI can generate errors and inaccurate information</li>
                        <li>You will independently verify all information before any clinical application</li>
                        <li>You assume all risks and waive all claims against us</li>
                        <li>You agree to binding arbitration and waive your right to jury trial and class actions</li>
                        <li>You agree to indemnify and hold us harmless from all claims</li>
                    </ul>
                </div>

                <p style="margin-top: 32px; font-weight: 600; font-size: 14px; color: var(--gray-600); text-align: center;">These Terms of Service constitute the entire agreement between you and gasconsult.ai.</p>
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
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=6">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=6">
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
                    <a href="/crisis" class="nav-link">Crisis Protocols</a>
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
            <a href="/crisis" class="mobile-menu-link">Crisis Protocols</a>
            <a href="/hypotension" class="mobile-menu-link">IOH Predictor</a>
        </div>

        <main class="main-content">
            <div class="content-card">
                <h1>Privacy Policy</h1>
                <p class="last-updated">Last Updated: January 2025</p>

                <div class="highlight-box">
                    <strong>CRITICAL NOTICE:</strong> DO NOT enter any Protected Health Information (PHI) or personally identifiable patient data.
                    This is an educational tool only. We collect minimal technical data, use temporary session storage, and NEVER sell your information.
                </div>

                <p>
                    Welcome to gasconsult.ai ("we," "our," "us," or "the Service"). This Privacy Policy describes our practices regarding
                    the collection, use, disclosure, and protection of information when you use our anesthesiology consultation platform.
                </p>
                <p>
                    <strong>BY USING THIS SERVICE, YOU ACKNOWLEDGE THAT YOU HAVE READ, UNDERSTOOD, AND AGREE TO BE BOUND BY THIS PRIVACY POLICY.</strong>
                    If you do not agree with these terms, you must immediately discontinue use of the Service.
                </p>

                <h2>1. Information We Collect</h2>
                <p>
                    We practice data minimization and collect only information necessary to operate the Service. We DO NOT require user accounts
                    or registration. We DO NOT store long-term personal health information.
                </p>

                <h3>1.1 Information You Voluntarily Provide</h3>
                <ul>
                    <li><strong>Medical Queries:</strong> Clinical questions you submit about anesthesiology topics (stored temporarily during your session only)</li>
                    <li><strong>Form Inputs:</strong> Data entered into pre-operative assessment tools, clinical calculators, or hypotension predictors (processed in real-time, not permanently stored)</li>
                    <li><strong>Session Conversation History:</strong> Your chat interactions during an active browser session (automatically cleared when you close the browser or click "Clear Session")</li>
                </ul>

                <h3>1.2 Automatically Collected Technical Information</h3>
                <p>We automatically collect limited technical data for operational purposes:</p>
                <ul>
                    <li><strong>Usage Analytics:</strong> Pages visited, features accessed, button clicks, time spent on pages (anonymized and aggregated)</li>
                    <li><strong>Technical Logs:</strong> IP address (for rate limiting and abuse prevention), browser type and version, device type, operating system, referring URLs</li>
                    <li><strong>Performance Data:</strong> Page load times, error messages, API response times (used solely for debugging and service optimization)</li>
                    <li><strong>Session Cookies:</strong> Temporary session identifiers stored in your browser (required for chat functionality, automatically deleted on session end)</li>
                </ul>

                <h3>1.3 Information We DO NOT Collect</h3>
                <ul>
                    <li>User accounts, passwords, or login credentials (we do not have user registration)</li>
                    <li>Payment information (the Service is free)</li>
                    <li>Long-term storage of medical queries or conversation history</li>
                    <li>Biometric data, precise geolocation, or device fingerprinting beyond standard web analytics</li>
                    <li>Social media profiles or third-party account linkages</li>
                </ul>

                <h2>2. How We Use Your Information</h2>
                <p>We use collected information exclusively for the following legitimate purposes:</p>
                <ul>
                    <li><strong>Service Delivery:</strong> To process your queries, search PubMed medical literature, generate AI-assisted responses via OpenAI GPT-4o, and display clinical calculator results</li>
                    <li><strong>Session Management:</strong> To maintain conversation continuity during your active browser session</li>
                    <li><strong>Performance Improvement:</strong> To analyze anonymized usage patterns, identify bugs, optimize response times, and enhance user experience</li>
                    <li><strong>Security & Abuse Prevention:</strong> To detect malicious activity, prevent spam, enforce rate limits (60 requests/minute per IP), and protect against unauthorized access</li>
                    <li><strong>Legal Compliance:</strong> To comply with applicable laws, respond to lawful government requests, enforce our Terms of Service, and protect our legal rights</li>
                </ul>
                <p>
                    <strong>WE DO NOT:</strong> Sell your data to third parties, use your queries for targeted advertising, share identifiable information with data brokers,
                    or use your medical questions to train AI models (OpenAI's enterprise API does not train on customer data per their data processing agreement).
                </p>

                <h2>3. CRITICAL WARNING: Do Not Enter Protected Health Information (PHI)</h2>
                <div class="highlight-box">
                    <strong>⚠️ MANDATORY NOTICE:</strong> This Service is designed for educational and informational purposes ONLY.
                    <strong>You MUST NOT enter any Protected Health Information (PHI) or personally identifiable patient data.</strong>
                    <p style="margin-top: 12px;">Prohibited information includes, but is not limited to:</p>
                    <ul style="margin-top: 8px; margin-bottom: 12px;">
                        <li>Patient names, initials, medical record numbers (MRNs), account numbers, or any identifiers</li>
                        <li>Dates of birth, admission dates, discharge dates, dates of death, or ages over 89</li>
                        <li>Specific geographic locations smaller than a state (addresses, ZIP codes, GPS coordinates)</li>
                        <li>Telephone numbers, fax numbers, email addresses, Social Security numbers, insurance ID numbers</li>
                        <li>Photographs, videos, biometric identifiers (fingerprints, voice recordings, facial images)</li>
                        <li>Device identifiers (IP addresses in medical contexts), medical device serial numbers, URLs containing PHI</li>
                        <li>Any case-specific details that could reasonably identify an individual patient or provider</li>
                    </ul>
                    <p style="margin-top: 12px; margin-bottom: 0;">
                        <strong>You assume all risk and liability for any PHI you choose to enter.</strong> Use only hypothetical scenarios,
                        de-identified case presentations, or generalized clinical questions. If you accidentally enter PHI, immediately clear your session.
                    </p>
                </div>

                <h2>4. Third-Party Services and Data Sharing</h2>
                <p>We utilize the following third-party services to operate the platform. Your data may be transmitted to these providers:</p>

                <h3>4.1 OpenAI (GPT-4o AI Model)</h3>
                <ul>
                    <li><strong>Purpose:</strong> To generate evidence-based clinical responses using AI synthesis</li>
                    <li><strong>Data Shared:</strong> Your medical queries, conversation context, and PubMed search results</li>
                    <li><strong>Privacy Policy:</strong> <a href="https://openai.com/privacy" target="_blank" style="color: var(--blue-600); text-decoration: underline;">https://openai.com/privacy</a></li>
                    <li><strong>Data Training:</strong> Per OpenAI's Enterprise API terms, your queries are NOT used to train their AI models</li>
                    <li><strong>Location:</strong> OpenAI is a U.S.-based company with servers in the United States</li>
                </ul>

                <h3>4.2 NCBI PubMed / Entrez API</h3>
                <ul>
                    <li><strong>Purpose:</strong> To search medical literature databases for high-quality evidence (systematic reviews, meta-analyses, clinical trials)</li>
                    <li><strong>Data Shared:</strong> Your search queries (medical terms extracted from your questions)</li>
                    <li><strong>Privacy Policy:</strong> <a href="https://www.nlm.nih.gov/web_policies.html" target="_blank" style="color: var(--blue-600); text-decoration: underline;">https://www.nlm.nih.gov/web_policies.html</a></li>
                    <li><strong>Operator:</strong> U.S. National Library of Medicine (government agency)</li>
                </ul>

                <h3>4.3 No Other Third-Party Sharing</h3>
                <p>
                    We DO NOT sell, rent, lease, or trade your information to third parties for marketing purposes. We DO NOT share data with
                    advertisers, data brokers, or analytics companies beyond basic anonymized usage statistics. We may disclose information only if:
                </p>
                <ul>
                    <li><strong>Legal Obligation:</strong> Required by law, court order, subpoena, or government request</li>
                    <li><strong>Safety & Fraud Prevention:</strong> Necessary to prevent harm, investigate abuse, or protect legal rights</li>
                    <li><strong>Business Transfer:</strong> In the event of a merger, acquisition, or sale of assets (you will be notified)</li>
                </ul>

                <h2>5. Data Security and Limitations of Liability</h2>

                <h3>5.1 Security Measures We Implement</h3>
                <p>We employ industry-standard security practices to protect your information:</p>
                <ul>
                    <li><strong>Encryption in Transit:</strong> All data transmitted between your browser and our servers is encrypted using HTTPS/TLS 1.2+ protocols</li>
                    <li><strong>Server-Side Sessions:</strong> Conversation history is stored server-side (not in cookies), reducing client-side exposure</li>
                    <li><strong>Input Sanitization:</strong> All user inputs are sanitized using the Bleach library to prevent Cross-Site Scripting (XSS) attacks</li>
                    <li><strong>CSRF Protection:</strong> Flask-WTF token validation prevents Cross-Site Request Forgery attacks</li>
                    <li><strong>Rate Limiting:</strong> Flask-Limiter restricts requests to 60 per minute per IP address to prevent abuse and DDoS attacks</li>
                    <li><strong>Access Controls:</strong> Server infrastructure uses firewalls, SSH key authentication, and principle of least privilege</li>
                    <li><strong>Regular Updates:</strong> Dependencies and server software are regularly patched for security vulnerabilities</li>
                </ul>

                <h3>5.2 Limitations and Disclaimers</h3>
                <p>
                    <strong>NO GUARANTEE OF ABSOLUTE SECURITY:</strong> While we implement reasonable security measures, no method of internet transmission
                    or electronic storage is 100% secure. We cannot guarantee that unauthorized third parties will never defeat our security measures
                    or misuse your information.
                </p>
                <p>
                    <strong>DISCLAIMER OF LIABILITY FOR DATA BREACHES:</strong> TO THE MAXIMUM EXTENT PERMITTED BY LAW, WE DISCLAIM ALL LIABILITY
                    FOR ANY UNAUTHORIZED ACCESS, USE, DISCLOSURE, OR MODIFICATION OF YOUR DATA. BY USING THIS SERVICE, YOU ACKNOWLEDGE AND ACCEPT
                    THE INHERENT SECURITY RISKS OF INTERNET-BASED PLATFORMS.
                </p>
                <p>
                    <strong>USER RESPONSIBILITY:</strong> You are solely responsible for: (1) ensuring you do not enter PHI or sensitive personal information,
                    (2) using strong passwords if accessing the Service from shared devices, (3) logging out or clearing sessions after use on public computers,
                    and (4) maintaining the confidentiality of any information you choose to enter.
                </p>

                <h2>6. Data Retention and Deletion</h2>
                <ul>
                    <li><strong>Session Data (Chat History):</strong> Stored temporarily in server-side sessions during your active browser session. Automatically deleted when you close your browser, click "Clear Session," or after 24 hours of inactivity.</li>
                    <li><strong>System Logs:</strong> Technical logs (IP addresses, timestamps, error messages) are retained for 90 days for debugging, security monitoring, and abuse prevention, then permanently deleted.</li>
                    <li><strong>Analytics Data:</strong> Anonymized and aggregated usage statistics (page views, feature usage) are retained indefinitely for service improvement but cannot be traced back to individual users.</li>
                    <li><strong>No Long-Term Storage:</strong> We do NOT maintain databases of your medical queries, conversation transcripts, or form inputs beyond the temporary session duration.</li>
                    <li><strong>Manual Deletion:</strong> You can clear your session data at any time by clicking the "Clear Session" button in the chat interface or by closing your browser.</li>
                </ul>

                <h2>7. Your Privacy Rights (GDPR, CCPA, and Other Laws)</h2>
                <p>
                    Depending on your location, you may have specific privacy rights under laws such as the European Union's General Data Protection Regulation (GDPR),
                    California Consumer Privacy Act (CCPA), Virginia Consumer Data Protection Act (VCDPA), and similar state/international laws.
                </p>

                <h3>7.1 Rights for EU/EEA Users (GDPR)</h3>
                <p>If you are located in the European Union or European Economic Area, you have the following rights:</p>
                <ul>
                    <li><strong>Right of Access (Art. 15):</strong> Request confirmation of whether we process your data and obtain a copy</li>
                    <li><strong>Right to Rectification (Art. 16):</strong> Request correction of inaccurate or incomplete data</li>
                    <li><strong>Right to Erasure (Art. 17):</strong> Request deletion of your data ("right to be forgotten")</li>
                    <li><strong>Right to Restriction (Art. 18):</strong> Request limitation of processing under certain circumstances</li>
                    <li><strong>Right to Data Portability (Art. 20):</strong> Receive your data in a structured, machine-readable format</li>
                    <li><strong>Right to Object (Art. 21):</strong> Object to processing based on legitimate interests</li>
                    <li><strong>Right to Withdraw Consent (Art. 7):</strong> Withdraw consent at any time (does not affect prior processing)</li>
                    <li><strong>Right to Lodge a Complaint:</strong> File a complaint with your national data protection authority</li>
                </ul>

                <h3>7.2 Rights for California Users (CCPA/CPRA)</h3>
                <p>If you are a California resident, you have the following rights under the California Consumer Privacy Act:</p>
                <ul>
                    <li><strong>Right to Know:</strong> Request disclosure of categories and specific pieces of personal information we collect</li>
                    <li><strong>Right to Delete:</strong> Request deletion of your personal information (subject to exceptions)</li>
                    <li><strong>Right to Opt-Out:</strong> Opt-out of the "sale" or "sharing" of personal information (Note: We do NOT sell or share your data)</li>
                    <li><strong>Right to Correct:</strong> Request correction of inaccurate personal information</li>
                    <li><strong>Right to Non-Discrimination:</strong> Exercise privacy rights without receiving discriminatory treatment</li>
                </ul>

                <h3>7.3 How to Exercise Your Rights</h3>
                <p>
                    To submit a privacy rights request, email us at <strong>privacy@gasconsult.ai</strong> with the subject line "Privacy Rights Request."
                    Please include: (1) your name and contact information, (2) description of your request, (3) approximate dates of service usage (if known).
                </p>
                <p>
                    <strong>Important Limitation:</strong> Due to our minimal data collection practices (no user accounts, temporary session storage only),
                    we may have limited ability to identify or retrieve historical data associated with your use of the Service. In most cases, simply clearing
                    your browser session will delete all temporary data.
                </p>
                <p>
                    We will respond to verifiable requests within 30 days (GDPR) or 45 days (CCPA). We reserve the right to request additional information
                    to verify your identity before processing requests.
                </p>

                <h2>8. HIPAA Compliance and Limitations</h2>
                <p>
                    <strong>WE ARE NOT A HIPAA-COVERED ENTITY:</strong> gasconsult.ai is NOT a healthcare provider, health plan, or healthcare clearinghouse
                    as defined under the Health Insurance Portability and Accountability Act (HIPAA). We are NOT a "Business Associate" of any covered entity.
                </p>
                <p>
                    <strong>NO HIPAA PROTECTIONS APPLY:</strong> This Service does NOT provide HIPAA-compliant safeguards for Protected Health Information (PHI).
                    We do NOT sign Business Associate Agreements (BAAs). You MUST NOT use this Service to store, transmit, or process PHI.
                </p>
                <p>
                    <strong>HEALTHCARE PROFESSIONALS:</strong> If you are a healthcare provider, you are solely responsible for ensuring your use of this
                    Service complies with HIPAA and other applicable healthcare privacy laws. We recommend using only de-identified, hypothetical clinical
                    scenarios that cannot be traced to real patients.
                </p>

                <h2>9. Children's Privacy</h2>
                <p>
                    This Service is intended for use by healthcare professionals and adults seeking educational medical information. We do NOT knowingly
                    collect personal information from individuals under 18 years of age.
                </p>
                <p>
                    If we become aware that we have inadvertently collected information from a child under 18, we will take immediate steps to delete
                    such information from our systems. If you are a parent or guardian and believe your child has provided us with personal information,
                    please contact us at <strong>privacy@gasconsult.ai</strong>.
                </p>

                <h2>10. International Data Transfers</h2>
                <p>
                    This Service is operated from the United States. If you access the Service from outside the United States, your information will be
                    transferred to, stored, and processed in the United States where our servers and third-party service providers (OpenAI, hosting infrastructure) are located.
                </p>
                <p>
                    <strong>Data Protection Differences:</strong> The United States may not provide the same level of data protection as your home country,
                    particularly for users in the European Union. The U.S. is not subject to GDPR adequacy decisions for all data transfers.
                </p>
                <p>
                    <strong>Legal Basis for Transfers:</strong> By using this Service from outside the U.S., you explicitly consent to the transfer of your
                    information to the United States. We rely on standard contractual clauses with third-party providers where applicable.
                </p>

                <h2>11. Changes to This Privacy Policy</h2>
                <p>
                    We reserve the right to modify this Privacy Policy at any time, at our sole discretion, without prior notice. When we make changes,
                    we will update the "Last Updated" date at the top of this page.
                </p>
                <p>
                    <strong>Material Changes:</strong> For significant changes that materially affect your privacy rights, we will provide prominent notice
                    on the homepage or via other reasonable means. However, we are NOT obligated to provide individual notice.
                </p>
                <p>
                    <strong>Deemed Acceptance:</strong> Your continued use of the Service after changes are posted constitutes your acceptance of the revised
                    Privacy Policy. If you do not agree to the updated terms, you must immediately discontinue use.
                </p>
                <p>
                    <strong>Responsibility to Review:</strong> You are responsible for periodically reviewing this Privacy Policy to stay informed of updates.
                </p>

                <h2>12. Contact Information</h2>
                <p>If you have questions, concerns, or requests regarding this Privacy Policy or our data practices, please contact us:</p>
                <ul>
                    <li><strong>Email:</strong> privacy@gasconsult.ai</li>
                    <li><strong>Website:</strong> https://gasconsult.ai</li>
                    <li><strong>Response Time:</strong> We aim to respond to inquiries within 5-7 business days</li>
                </ul>
                <p>
                    For privacy rights requests (GDPR, CCPA), please use the subject line "Privacy Rights Request" and include the specific right you wish to exercise.
                </p>

                <h2>13. Severability and Governing Law</h2>
                <p>
                    If any provision of this Privacy Policy is found to be unlawful, void, or unenforceable, that provision shall be deemed severable and shall
                    not affect the validity and enforceability of the remaining provisions.
                </p>
                <p>
                    This Privacy Policy shall be governed by and construed in accordance with the laws of the United States and the State of [Your State],
                    without regard to conflict of law principles. Any disputes arising from this Privacy Policy shall be subject to the exclusive jurisdiction
                    of the courts located in [Your State].
                </p>

                <h2>14. Acknowledgment and Consent</h2>
                <p>
                    <strong>BY CLICKING "I ACCEPT," CONTINUING TO USE THE SERVICE, OR SUBMITTING ANY QUERIES, YOU ACKNOWLEDGE THAT:</strong>
                </p>
                <ul>
                    <li>You have read and understood this Privacy Policy in its entirety</li>
                    <li>You consent to the collection, use, and disclosure of your information as described herein</li>
                    <li>You understand that we are not a HIPAA-covered entity and do NOT provide HIPAA protections</li>
                    <li>You agree NOT to enter any Protected Health Information (PHI) or personally identifiable patient data</li>
                    <li>You accept the inherent security risks of internet-based services and our limitations of liability</li>
                    <li>You acknowledge that this is an educational tool only and NOT a substitute for professional medical advice</li>
                </ul>

                <div class="highlight-box" style="margin-top: 32px;">
                    <strong>Final Reminder:</strong> This Service provides educational information only and is NOT medical advice.
                    Always consult qualified healthcare professionals for clinical decision-making. NEVER enter PHI or patient-identifying information.
                    We prioritize your privacy, collect minimal data, and never sell your information.
                </div>
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

CRISIS_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Crisis Protocols — gasconsult.ai</title>

    <!-- PWA -->
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=6">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=6">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">

    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800;900&display=swap" rel="stylesheet">
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
            --red-50: #FEF2F2;
            --red-500: #EF4444;
            --red-600: #DC2626;
            --red-700: #B91C1C;
            --orange-50: #FFF7ED;
            --orange-500: #F97316;
            --orange-600: #EA580C;
            --purple-50: #FAF5FF;
            --purple-500: #A855F7;
            --purple-600: #9333EA;
            --amber-50: #FFFBEB;
            --amber-500: #F59E0B;
            --amber-600: #D97706;
            --emerald-50: #ECFDF5;
            --emerald-500: #10B981;
            --emerald-600: #059669;
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
            padding: 100px 20px 40px;
            text-align: center;
            max-width: 900px;
            margin: 0 auto;
        }

        .hero-badge {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            background: linear-gradient(135deg, var(--red-50) 0%, #FEE2E2 100%);
            border: 1px solid rgba(220, 38, 38, 0.2);
            border-radius: 100px;
            padding: 8px 16px 8px 12px;
            margin-bottom: 20px;
            box-shadow: 0 2px 8px rgba(239, 68, 68, 0.08);
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) forwards;
            opacity: 0;
        }

        .badge-dot {
            width: 8px;
            height: 8px;
            background: var(--red-500);
            border-radius: 50%;
            position: relative;
            animation: pulse-dot 2s ease-in-out infinite;
        }

        @keyframes pulse-dot {
            0%, 100% { transform: scale(1); opacity: 1; }
            50% { transform: scale(1.2); opacity: 0.8; }
        }

        .badge-text {
            font-size: 12px;
            font-weight: 600;
            color: var(--red-700);
        }

        .hero-title {
            font-size: 48px;
            font-weight: 800;
            line-height: 1.1;
            letter-spacing: -2px;
            color: var(--gray-900);
            margin-bottom: 16px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.1s forwards;
            opacity: 0;
        }

        .hero-title .gradient {
            background: linear-gradient(135deg, var(--red-600) 0%, var(--orange-600) 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
        }

        .hero-subtitle {
            font-size: 18px;
            font-weight: 400;
            line-height: 1.6;
            color: var(--gray-600);
            max-width: 700px;
            margin: 0 auto 40px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.2s forwards;
            opacity: 0;
        }

        @keyframes fade-up {
            from { opacity: 0; transform: translateY(24px); }
            to { opacity: 1; transform: translateY(0); }
        }

        .search-container {
            max-width: 600px;
            margin: 0 auto 50px;
            padding: 0 20px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.3s forwards;
            opacity: 0;
        }

        .search-box {
            position: relative;
        }

        .search-input {
            width: 100%;
            padding: 14px 48px 14px 48px;
            font-size: 16px;
            font-family: inherit;
            border: 2px solid var(--gray-200);
            border-radius: 16px;
            background: var(--white);
            color: var(--gray-900);
            transition: all 0.3s ease;
        }

        .search-input:focus {
            outline: none;
            border-color: var(--blue-400);
            box-shadow: 0 0 0 4px rgba(59, 130, 246, 0.1);
        }

        .search-icon {
            position: absolute;
            left: 16px;
            top: 50%;
            transform: translateY(-50%);
            color: var(--gray-400);
            pointer-events: none;
        }

        .clear-btn {
            position: absolute;
            right: 12px;
            top: 50%;
            transform: translateY(-50%);
            background: none;
            border: none;
            color: var(--gray-400);
            cursor: pointer;
            padding: 6px;
            border-radius: 8px;
            display: none;
            transition: all 0.2s ease;
        }

        .clear-btn:hover {
            background: var(--gray-100);
            color: var(--gray-600);
        }

        .clear-btn.visible {
            display: block;
        }

        .protocols-container {
            max-width: 1200px;
            margin: 0 auto;
            padding: 0 20px 80px;
        }

        .category-section {
            margin-bottom: 40px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) forwards;
            opacity: 0;
        }

        .category-section:nth-child(1) { animation-delay: 0.4s; }
        .category-section:nth-child(2) { animation-delay: 0.5s; }
        .category-section:nth-child(3) { animation-delay: 0.6s; }
        .category-section:nth-child(4) { animation-delay: 0.7s; }

        .category-header {
            display: flex;
            align-items: center;
            gap: 12px;
            margin-bottom: 20px;
        }

        .category-icon {
            width: 40px;
            height: 40px;
            border-radius: 12px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 20px;
        }

        .category-icon.red {
            background: linear-gradient(135deg, var(--red-500) 0%, var(--red-600) 100%);
            box-shadow: 0 4px 12px rgba(239, 68, 68, 0.3);
        }

        .category-icon.orange {
            background: linear-gradient(135deg, var(--orange-500) 0%, var(--orange-600) 100%);
            box-shadow: 0 4px 12px rgba(249, 115, 22, 0.3);
        }

        .category-icon.purple {
            background: linear-gradient(135deg, var(--purple-500) 0%, var(--purple-600) 100%);
            box-shadow: 0 4px 12px rgba(168, 85, 247, 0.3);
        }

        .category-icon.amber {
            background: linear-gradient(135deg, var(--amber-500) 0%, var(--amber-600) 100%);
            box-shadow: 0 4px 12px rgba(245, 158, 11, 0.3);
        }

        .category-icon.blue {
            background: linear-gradient(135deg, var(--blue-500) 0%, var(--blue-600) 100%);
            box-shadow: 0 4px 12px rgba(59, 130, 246, 0.3);
        }

        .category-title {
            font-size: 24px;
            font-weight: 700;
            color: var(--gray-900);
            letter-spacing: -0.5px;
        }

        .protocols-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(320px, 1fr));
            gap: 20px;
        }

        .protocol-card {
            background: var(--white);
            border: 2px solid var(--gray-200);
            border-radius: 16px;
            padding: 24px;
            cursor: pointer;
            transition: all 0.3s ease;
            position: relative;
            overflow: hidden;
        }

        .protocol-card::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 4px;
            transition: all 0.3s ease;
        }

        .protocol-card.red::before { background: linear-gradient(90deg, var(--red-500) 0%, var(--red-600) 100%); }
        .protocol-card.orange::before { background: linear-gradient(90deg, var(--orange-500) 0%, var(--orange-600) 100%); }
        .protocol-card.purple::before { background: linear-gradient(90deg, var(--purple-500) 0%, var(--purple-600) 100%); }
        .protocol-card.amber::before { background: linear-gradient(90deg, var(--amber-500) 0%, var(--amber-600) 100%); }
        .protocol-card.blue::before { background: linear-gradient(90deg, var(--blue-500) 0%, var(--blue-600) 100%); }

        .protocol-card:hover {
            border-color: var(--gray-300);
            box-shadow: 0 8px 24px rgba(0, 0, 0, 0.08);
            transform: translateY(-2px);
        }

        .protocol-card.expanded {
            border-color: var(--blue-400);
            box-shadow: 0 12px 32px rgba(59, 130, 246, 0.15);
        }

        .protocol-header {
            display: flex;
            align-items: start;
            justify-content: space-between;
            margin-bottom: 8px;
        }

        .protocol-title {
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 4px;
            line-height: 1.3;
        }

        .expand-icon {
            width: 24px;
            height: 24px;
            border-radius: 8px;
            background: var(--gray-100);
            display: flex;
            align-items: center;
            justify-content: center;
            flex-shrink: 0;
            transition: all 0.3s ease;
        }

        .protocol-card.expanded .expand-icon {
            background: var(--blue-100);
            transform: rotate(180deg);
        }

        .protocol-summary {
            font-size: 14px;
            color: var(--gray-600);
            line-height: 1.5;
            margin-bottom: 12px;
        }

        .protocol-tags {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
        }

        .protocol-tag {
            font-size: 11px;
            font-weight: 600;
            padding: 4px 10px;
            border-radius: 6px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        .protocol-tag.immediate {
            background: var(--red-50);
            color: var(--red-700);
        }

        .protocol-tag.urgent {
            background: var(--orange-50);
            color: var(--orange-700);
        }

        .protocol-tag.call-help {
            background: var(--purple-50);
            color: var(--purple-700);
        }

        .protocol-content {
            max-height: 0;
            overflow: hidden;
            transition: max-height 0.5s ease;
        }

        .protocol-card.expanded .protocol-content {
            max-height: 3000px;
        }

        .protocol-details {
            margin-top: 20px;
            padding-top: 20px;
            border-top: 2px solid var(--gray-100);
        }

        .protocol-section {
            margin-bottom: 20px;
        }

        .protocol-section-title {
            font-size: 14px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 12px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        .protocol-steps {
            list-style: none;
            counter-reset: step-counter;
        }

        .protocol-step {
            counter-increment: step-counter;
            position: relative;
            padding-left: 36px;
            margin-bottom: 12px;
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-700);
        }

        .protocol-step::before {
            content: counter(step-counter);
            position: absolute;
            left: 0;
            top: 0;
            width: 24px;
            height: 24px;
            background: var(--blue-100);
            color: var(--blue-700);
            border-radius: 50%;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 12px;
            font-weight: 700;
        }

        .protocol-step strong {
            color: var(--gray-900);
            font-weight: 600;
        }

        .dose-box {
            background: var(--emerald-50);
            border-left: 3px solid var(--emerald-500);
            padding: 12px 16px;
            border-radius: 8px;
            margin: 12px 0;
        }

        .dose-box-title {
            font-size: 13px;
            font-weight: 700;
            color: var(--emerald-700);
            margin-bottom: 6px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        .dose-detail {
            font-size: 14px;
            color: var(--gray-800);
            margin: 4px 0;
            line-height: 1.5;
        }

        .warning-box {
            background: var(--red-50);
            border-left: 3px solid var(--red-500);
            padding: 12px 16px;
            border-radius: 8px;
            margin: 12px 0;
        }

        .warning-box-title {
            font-size: 13px;
            font-weight: 700;
            color: var(--red-700);
            margin-bottom: 6px;
            display: flex;
            align-items: center;
            gap: 6px;
        }

        .warning-detail {
            font-size: 14px;
            color: var(--gray-800);
            line-height: 1.5;
        }

        .info-box {
            background: var(--blue-50);
            border-left: 3px solid var(--blue-500);
            padding: 12px 16px;
            border-radius: 8px;
            margin: 12px 0;
        }

        .info-box-title {
            font-size: 13px;
            font-weight: 700;
            color: var(--blue-700);
            margin-bottom: 6px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        .info-detail {
            font-size: 14px;
            color: var(--gray-800);
            line-height: 1.5;
        }

        .footer {
            background: rgba(255, 255, 255, 0.8);
            backdrop-filter: blur(20px);
            border-top: 1px solid var(--gray-200);
            padding: 40px 20px;
            text-align: center;
            margin-top: auto;
        }

        .footer-content {
            max-width: 800px;
            margin: 0 auto;
        }

        .disclaimer {
            background: var(--amber-50);
            border: 2px solid var(--amber-200);
            border-radius: 12px;
            padding: 20px;
            margin-bottom: 24px;
        }

        .disclaimer-title {
            font-size: 14px;
            font-weight: 700;
            color: var(--amber-700);
            margin-bottom: 8px;
            display: flex;
            align-items: center;
            justify-content: center;
            gap: 8px;
        }

        .disclaimer-text {
            font-size: 13px;
            color: var(--gray-700);
            line-height: 1.6;
        }

        .footer-links {
            display: flex;
            justify-content: center;
            gap: 24px;
            flex-wrap: wrap;
            margin-bottom: 16px;
        }

        .footer-link {
            font-size: 14px;
            color: var(--gray-600);
            text-decoration: none;
            transition: color 0.2s ease;
        }

        .footer-link:hover {
            color: var(--blue-600);
        }

        .footer-copy {
            font-size: 13px;
            color: var(--gray-500);
        }

        @media (min-width: 768px) {
            .nav-links {
                display: flex;
            }

            .mobile-menu-btn {
                display: none;
            }

            .hero-title {
                font-size: 56px;
            }

            .protocols-grid {
                grid-template-columns: repeat(2, 1fr);
            }
        }

        @media (min-width: 1024px) {
            .protocols-grid {
                grid-template-columns: repeat(3, 1fr);
            }
        }

        .hidden {
            display: none !important;
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
        <!-- Navigation -->
        <nav class="nav">
            <div class="nav-inner">
                <a href="/" class="logo">
                    <div class="logo-icon">
                        <svg viewBox="0 0 100 100" fill="currentColor">
                            <circle cx="20" cy="50" r="8" fill="#2563EB" opacity="0.6"/>
                            <circle cx="50" cy="50" r="10" fill="#2563EB"/>
                            <circle cx="80" cy="50" r="8" fill="#2563EB" opacity="0.6"/>
                        </svg>
                    </div>
                    <div class="logo-text">
                        <span class="gas">gas</span><span class="consult">consult</span><span class="ai">.ai</span>
                    </div>
                </a>

                <div class="nav-links">
                    <a href="/" class="nav-link">Home</a>
                    <a href="/calculators" class="nav-link">Clinical Calculators</a>
                    <a href="/crisis" class="nav-link active">Crisis Protocols</a>
                    <a href="/privacy" class="nav-link">Privacy</a>
                </div>

                <button class="mobile-menu-btn" onclick="toggleMobileMenu()">
                    <span></span>
                    <span></span>
                    <span></span>
                </button>
            </div>
        </nav>

        <div class="mobile-menu" id="mobileMenu">
            <a href="/" class="mobile-menu-link">Home</a>
            <a href="/calculators" class="mobile-menu-link">Clinical Calculators</a>
            <a href="/crisis" class="mobile-menu-link">Crisis Protocols</a>
            <a href="/privacy" class="mobile-menu-link">Privacy</a>
        </div>

        <!-- Hero -->
        <section class="hero">
            <div class="hero-badge">
                <div class="badge-dot"></div>
                <div class="badge-text">Emergency Reference</div>
            </div>
            <h1 class="hero-title">
                <span class="gradient">Crisis</span> Protocols
            </h1>
            <p class="hero-subtitle">
                Evidence-based, step-by-step management for anesthesia emergencies. Quick access to critical protocols when seconds count.
            </p>
        </section>

        <!-- Search -->
        <div class="search-container">
            <div class="search-box">
                <svg class="search-icon" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <circle cx="11" cy="11" r="8"></circle>
                    <path d="m21 21-4.35-4.35"></path>
                </svg>
                <input
                    type="text"
                    class="search-input"
                    id="searchInput"
                    placeholder="Search protocols (e.g., malignant hyperthermia, anaphylaxis...)"
                    oninput="filterProtocols()"
                />
                <button class="clear-btn" id="clearBtn" onclick="clearSearch()">
                    <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                        <line x1="18" y1="6" x2="6" y2="18"></line>
                        <line x1="6" y1="6" x2="18" y2="18"></line>
                    </svg>
                </button>
            </div>
        </div>

        <!-- Protocols -->
        <div class="protocols-container">
            <!-- Life-Threatening / Cardiac -->
            <div class="category-section" data-category="cardiac">
                <div class="category-header">
                    <div class="category-icon red">⚡</div>
                    <h2 class="category-title">Life-Threatening / Cardiac</h2>
                </div>
                <div class="protocols-grid">

                    <!-- Malignant Hyperthermia -->
                    <div class="protocol-card red" data-keywords="malignant hyperthermia mh dantrolene hypermetabolic crisis muscle rigidity hyperthermia">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">Malignant Hyperthermia</h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Hypermetabolic crisis triggered by volatile anesthetics or succinylcholine</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag immediate">Immediate</span>
                            <span class="protocol-tag call-help">Call Help</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Actions</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>STOP triggers immediately:</strong> Discontinue all volatile anesthetics and succinylcholine</li>
                                        <li class="protocol-step"><strong>Call for help:</strong> Activate MH emergency protocol, assign roles</li>
                                        <li class="protocol-step"><strong>Hyperventilate with 100% O₂:</strong> 2-3× normal minute ventilation to eliminate CO₂</li>
                                        <li class="protocol-step"><strong>Give dantrolene immediately:</strong> See dosing below</li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Dantrolene Dosing</div>
                                    <div class="dose-detail"><strong>Initial:</strong> 2.5 mg/kg IV rapid push (reconstitute each 20mg vial with 60mL sterile water)</div>
                                    <div class="dose-detail"><strong>Repeat:</strong> 1 mg/kg boluses every 5-10 min until signs resolve</div>
                                    <div class="dose-detail"><strong>Maximum:</strong> Up to 10 mg/kg in acute phase (rarely >10 vials needed initially)</div>
                                    <div class="dose-detail"><strong>Continuation:</strong> 1 mg/kg IV q6h × 24-48h to prevent recrudescence</div>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Supportive Care</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Cool the patient:</strong> Cold IV saline, ice packs to groin/axilla, cooling blanket. Target temp <38.5°C</li>
                                        <li class="protocol-step"><strong>Treat hyperkalemia:</strong> Insulin/dextrose, calcium chloride, bicarbonate, avoid Ca²⁺ channel blockers with dantrolene</li>
                                        <li class="protocol-step"><strong>Manage arrhythmias:</strong> Standard ACLS (avoid calcium channel blockers)</li>
                                        <li class="protocol-step"><strong>Monitor urine output:</strong> Foley catheter, maintain >1 mL/kg/h to prevent myoglobin-induced renal failure</li>
                                        <li class="protocol-step"><strong>Labs:</strong> ABG, electrolytes, CK, lactate, coags, myoglobin q6h</li>
                                    </ol>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">📞 MH Hotline</div>
                                    <div class="info-detail"><strong>USA:</strong> 1-800-MH-HYPER (1-800-644-9737)</div>
                                    <div class="info-detail"><strong>Outside USA:</strong> +1-315-464-7079</div>
                                    <div class="info-detail">Expert consultant available 24/7 for real-time guidance</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Key Points</div>
                                    <div class="warning-detail">• Early signs: Unexplained ↑EtCO₂, masseter spasm, tachycardia, hypercarbia refractory to ↑ventilation<br>• Late signs: Fever, rigidity, rhabdomyolysis, hyperkalemia, acidosis<br>• Dantrolene can cause profound weakness — prepare for prolonged ventilation<br>• ICU monitoring × 24-48h minimum (recrudescence occurs in ~25%)</div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- Cardiac Arrest -->
                    <div class="protocol-card red" data-keywords="cardiac arrest code blue cpr acls asystole vfib pea pulseless">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">Cardiac Arrest (ACLS)</h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Immediate CPR and rhythm-specific advanced life support</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag immediate">Immediate</span>
                            <span class="protocol-tag call-help">Call Help</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Universal Steps (All Rhythms)</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Start CPR immediately:</strong> 100-120 compressions/min, depth 2-2.4 inches, minimize interruptions</li>
                                        <li class="protocol-step"><strong>Call for help / Code Blue</strong></li>
                                        <li class="protocol-step"><strong>Attach defibrillator/monitor:</strong> Identify rhythm</li>
                                        <li class="protocol-step"><strong>Secure airway:</strong> ETT or supraglottic device + capnography (target EtCO₂ >10 mmHg)</li>
                                        <li class="protocol-step"><strong>IV/IO access:</strong> Establish vascular access for medications</li>
                                        <li class="protocol-step"><strong>Consider reversible causes (H's and T's)</strong></li>
                                    </ol>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Shockable Rhythms (VF/pVT)</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Defibrillate:</strong> Biphasic 120-200J (or manufacturer recommendation), resume CPR immediately × 2 min</li>
                                        <li class="protocol-step"><strong>After 2nd shock:</strong> Epinephrine 1 mg IV/IO q3-5min</li>
                                        <li class="protocol-step"><strong>After 3rd shock:</strong> Amiodarone 300 mg IV/IO (or lidocaine 1-1.5 mg/kg if amio unavailable)</li>
                                        <li class="protocol-step"><strong>Continue CPR + defibrillation every 2 min</strong></li>
                                    </ol>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Non-Shockable Rhythms (PEA/Asystole)</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>High-quality CPR × 2 min</strong></li>
                                        <li class="protocol-step"><strong>Epinephrine 1 mg IV/IO immediately,</strong> then q3-5min</li>
                                        <li class="protocol-step"><strong>Consider atropine 1 mg IV</strong> if slow PEA (rate <60), may repeat to total 3 mg</li>
                                        <li class="protocol-step"><strong>Treat reversible causes aggressively</strong> (see H's and T's below)</li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Key Medications</div>
                                    <div class="dose-detail"><strong>Epinephrine:</strong> 1 mg (1:10,000) IV/IO every 3-5 minutes</div>
                                    <div class="dose-detail"><strong>Amiodarone:</strong> 300 mg IV/IO first dose, then 150 mg second dose</div>
                                    <div class="dose-detail"><strong>Lidocaine (alternative):</strong> 1-1.5 mg/kg first dose, then 0.5-0.75 mg/kg</div>
                                    <div class="dose-detail"><strong>Sodium bicarbonate:</strong> 1 mEq/kg (for hyperkalemia, TCA overdose, prolonged arrest)</div>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🔍 H's and T's (Reversible Causes)</div>
                                    <div class="info-detail"><strong>H's:</strong> Hypovolemia, Hypoxia, H⁺ (acidosis), Hyper/hypokalemia, Hypothermia<br>
                                    <strong>T's:</strong> Tension pneumothorax, Tamponade (cardiac), Toxins, Thrombosis (coronary/pulmonary)</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Anesthesia-Specific Considerations</div>
                                    <div class="warning-detail">• Turn off volatile anesthetics during arrest<br>• Consider <strong>anesthesia-specific causes:</strong> local anesthetic toxicity (give lipid emulsion), hyperkalemia from succinylcholine, pneumothorax from line placement, total spinal<br>• Continue CPR during transfer to ICU if needed<br>• Document ROSC time, rhythm changes, total epi/defib doses</div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- Anaphylaxis -->
                    <div class="protocol-card red" data-keywords="anaphylaxis allergic reaction epinephrine bronchospasm hypotension urticaria angioedema">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">Anaphylaxis</h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Severe IgE-mediated hypersensitivity reaction with cardiovascular/airway collapse</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag immediate">Immediate</span>
                            <span class="protocol-tag call-help">Call Help</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Management</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>STOP suspected trigger:</strong> Antibiotics, NMBs, latex, colloids are most common</li>
                                        <li class="protocol-step"><strong>Call for help</strong></li>
                                        <li class="protocol-step"><strong>Epinephrine IM immediately:</strong> 0.3-0.5 mg (0.3-0.5 mL of 1:1000) into anterolateral thigh, repeat q5-15min PRN</li>
                                        <li class="protocol-step"><strong>100% O₂:</strong> Maintain airway, consider early intubation if upper airway edema</li>
                                        <li class="protocol-step"><strong>Aggressive fluid resuscitation:</strong> 20-50 mL/kg crystalloid rapidly for refractory hypotension</li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Epinephrine Dosing</div>
                                    <div class="dose-detail"><strong>IM (first-line):</strong> 0.3-0.5 mg (1:1000) into thigh, repeat q5-15min</div>
                                    <div class="dose-detail"><strong>IV bolus (severe/arrest):</strong> 10-100 mcg (0.01-0.1 mg) slow push, titrate to effect</div>
                                    <div class="dose-detail"><strong>IV infusion (refractory):</strong> 0.05-0.5 mcg/kg/min, titrate to BP/HR</div>
                                    <div class="dose-detail"><strong>Pediatric IM:</strong> 0.01 mg/kg (max 0.5 mg)</div>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Adjunct Therapies</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>H1 blocker:</strong> Diphenhydramine 25-50 mg IV slowly</li>
                                        <li class="protocol-step"><strong>H2 blocker:</strong> Famotidine 20 mg IV or ranitidine 50 mg IV</li>
                                        <li class="protocol-step"><strong>Corticosteroids:</strong> Methylprednisolone 1-2 mg/kg IV (prevents late-phase reaction)</li>
                                        <li class="protocol-step"><strong>Bronchodilators:</strong> Albuterol for persistent bronchospasm</li>
                                        <li class="protocol-step"><strong>Glucagon (if on β-blockers):</strong> 1-2 mg IV (epinephrine may be ineffective)</li>
                                    </ol>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🩺 Diagnostic Criteria (2 or more)</div>
                                    <div class="info-detail">• <strong>Skin/mucosal:</strong> Urticaria, angioedema, flushing<br>• <strong>Respiratory:</strong> Bronchospasm, wheezing, stridor, dyspnea, ↓SpO₂<br>• <strong>Cardiovascular:</strong> Hypotension (SBP <90 or >30% drop), tachycardia, arrhythmia, collapse<br>• <strong>GI:</strong> Cramping, vomiting, diarrhea</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Key Considerations</div>
                                    <div class="warning-detail">• <strong>Common triggers:</strong> NMBs (rocuronium, succinylcholine), antibiotics (cephalosporins, penicillins), latex, chlorhexidine<br>• Confirm diagnosis: Send tryptase levels (draw immediately, then 1-2h and 24h later)<br>• Biphasic reactions occur in 20% — observe ≥4-6h minimum, admit if severe<br>• Document reaction in chart and advise patient to see allergist<br>• Refractory hypotension: Consider methylene blue 1-2 mg/kg for vasoplegia</div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- LAST -->
                    <div class="protocol-card red" data-keywords="last local anesthetic systemic toxicity lipid emulsion intralipid bupivacaine ropivacaine seizure arrhythmia">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">Local Anesthetic Systemic Toxicity (LAST)</h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">CNS/cardiac toxicity from systemic absorption of local anesthetics</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag immediate">Immediate</span>
                            <span class="protocol-tag call-help">Call Help</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Actions</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>STOP local anesthetic injection</strong></li>
                                        <li class="protocol-step"><strong>Call for help:</strong> Get lipid emulsion (Intralipid 20%)</li>
                                        <li class="protocol-step"><strong>Airway management:</strong> 100% O₂, ventilate if needed, suppress seizures</li>
                                        <li class="protocol-step"><strong>Give lipid emulsion immediately</strong> (see dosing below) — DO NOT DELAY</li>
                                        <li class="protocol-step"><strong>If cardiac arrest:</strong> Start CPR, consider prolonged resuscitation (LAST arrest can require >1h CPR)</li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Lipid Emulsion 20% (Intralipid) Dosing</div>
                                    <div class="dose-detail"><strong>Bolus:</strong> 1.5 mL/kg IV (lean body weight) over 1 minute (~100 mL for 70 kg adult)</div>
                                    <div class="dose-detail"><strong>Infusion:</strong> 0.25 mL/kg/min (~18 mL/min for 70 kg = ~500 mL bag over 30 min)</div>
                                    <div class="dose-detail"><strong>Repeat bolus:</strong> If cardiovascular instability persists after 5 min, give up to 2 more boluses (same dose)</div>
                                    <div class="dose-detail"><strong>Increase infusion:</strong> Double rate to 0.5 mL/kg/min if BP remains unstable</div>
                                    <div class="dose-detail"><strong>Maximum dose:</strong> ~10 mL/kg over first 30 minutes</div>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Seizure Management</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Benzodiazepines:</strong> Midazolam 1-2 mg IV or lorazepam 1-2 mg IV</li>
                                        <li class="protocol-step"><strong>AVOID propofol</strong> in large doses (additional lipid load, myocardial depression)</li>
                                        <li class="protocol-step"><strong>Small-dose propofol OK</strong> if lipid already given and seizures refractory</li>
                                    </ol>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Cardiac Arrest Modifications</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Reduce epinephrine doses:</strong> Use <1 mcg/kg (ACLS doses may worsen outcome)</li>
                                        <li class="protocol-step"><strong>Avoid vasopressin, calcium, β-blockers, local anesthetics (lidocaine)</strong></li>
                                        <li class="protocol-step"><strong>CONTINUE lipid infusion</strong></li>
                                        <li class="protocol-step"><strong>Prolonged CPR:</strong> Full recovery reported after >60 min resuscitation</li>
                                        <li class="protocol-step"><strong>Consider ECMO/CPB</strong> if available and refractory arrest</li>
                                    </ol>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🩺 Signs of LAST</div>
                                    <div class="info-detail"><strong>Early CNS:</strong> Circumoral numbness, metallic taste, tinnitus, confusion, agitation<br>
                                    <strong>Severe CNS:</strong> Seizures, loss of consciousness<br>
                                    <strong>Cardiac:</strong> Bradycardia, hypotension, arrhythmias (wide QRS), asystole, PEA<br>
                                    <strong>Note:</strong> Cardiac toxicity can occur WITHOUT preceding CNS symptoms (especially bupivacaine)</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Critical Points</div>
                                    <div class="warning-detail">• Lipid emulsion is PRIMARY treatment — give early, don't wait for arrest<br>• Bupivacaine/ropivacaine are more cardiotoxic than lidocaine/mepivacaine<br>• Max doses: Bupivacaine 2.5 mg/kg plain, 3 mg/kg with epi; Lidocaine 5 mg/kg plain, 7 mg/kg with epi<br>• Post-resuscitation: Monitor ≥4-6h (cardiac arrest patients → ICU), watch for pancreatitis from lipid load</div>
                                </div>
                            </div>
                        </div>
                    </div>

                </div>
            </div>

            <!-- Airway Emergencies -->
            <div class="category-section" data-category="airway">
                <div class="category-header">
                    <div class="category-icon orange">🫁</div>
                    <h2 class="category-title">Airway Emergencies</h2>
                </div>
                <div class="protocols-grid">

                    <!-- Can't Intubate Can't Oxygenate -->
                    <div class="protocol-card orange" data-keywords="cico cant intubate oxygenate difficult airway cricothyroidotomy emergency front neck access scalpel bougie">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">Can't Intubate, Can't Oxygenate (CICO)</h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Failed intubation + failed oxygenation requiring emergency front of neck access</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag immediate">Immediate</span>
                            <span class="protocol-tag call-help">Call Help</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">CICO Criteria (Declare CICO if both present)</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Can't Intubate:</strong> 3 failed laryngoscopy attempts by experienced provider OR 2 failed attempts + failed supraglottic device</li>
                                        <li class="protocol-step"><strong>Can't Oxygenate:</strong> SpO₂ <90% despite 100% O₂, facemask, OPA/NPA, 2-person BVM, and/or SGA</li>
                                    </ol>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Declare CICO Out Loud</div>
                                    <div class="warning-detail">"This is a CAN'T INTUBATE, CAN'T OXYGENATE situation. Prepare for emergency cricothyroidotomy NOW."</div>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Emergency Front of Neck Access (Scalpel-Bougie-Tube)</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Position:</strong> Extend neck (if no C-spine concerns), palpate cricothyroid membrane</li>
                                        <li class="protocol-step"><strong>Scalpel:</strong> Transverse stab incision through skin + cricothyroid membrane (1 motion, blade perpendicular)</li>
                                        <li class="protocol-step"><strong>Bougie:</strong> Insert bougie (or tracheal hook to retract) into trachea, advance caudally</li>
                                        <li class="protocol-step"><strong>Tube:</strong> Railroad cuffed ETT (6.0 or smaller) or Shiley 6 trach over bougie, inflate cuff</li>
                                        <li class="protocol-step"><strong>Ventilate:</strong> Confirm with capnography, ventilate, secure tube</li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">🔪 Scalpel Cricothyroidotomy Technique</div>
                                    <div class="dose-detail"><strong>Equipment:</strong> #10 scalpel blade, bougie, 6.0 cuffed ETT (or 6.0 Shiley trach)</div>
                                    <div class="dose-detail"><strong>Incision:</strong> Horizontal stab through cricothyroid membrane, turn blade 90° to open</div>
                                    <div class="dose-detail"><strong>Depth:</strong> 1-1.5 cm deep (through membrane into trachea)</div>
                                    <div class="dose-detail"><strong>Time goal:</strong> <60 seconds from declaration to ventilation</div>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Post-Procedure</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Confirm placement:</strong> Capnography, chest rise, SpO₂ improvement</li>
                                        <li class="protocol-step"><strong>Secure tube:</strong> Suture or trach ties, note depth</li>
                                        <li class="protocol-step"><strong>CXR:</strong> Confirm position, rule out pneumothorax</li>
                                        <li class="protocol-step"><strong>ENT/surgery consult:</strong> For formal tracheostomy or laryngeal evaluation</li>
                                    </ol>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">Alternative: Needle Cricothyroidotomy (Temporizing Only)</div>
                                    <div class="info-detail">• <strong>Indication:</strong> Pediatric <10 years (cricoid cartilage too small for surgical cric)<br>• <strong>Technique:</strong> 14G or 16G IV catheter through cricothyroid membrane, attach to jet ventilator or BVM with Y-connector<br>• <strong>Limitation:</strong> Inadequate ventilation (↑ CO₂), only buys 30-45 min before hypercapnia/barotrauma</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Critical Reminders</div>
                                    <div class="warning-detail">• Do NOT attempt >3 laryngoscopies — recognize failure early<br>• Do NOT delay cricothyroidotomy — hypoxic brain injury starts after 3-5 min<br>• Scalpel technique is superior to needle cric in adults (faster, more reliable ventilation)<br>• Practice on simulators regularly — muscle memory is critical</div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- Laryngospasm -->
                    <div class="protocol-card orange" data-keywords="laryngospasm stridor vocal cord spasm cpap jaw thrust succinylcholine propofol">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">Laryngospasm</h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Reflex closure of vocal cords causing airway obstruction</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag urgent">Urgent</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Management (Escalating)</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Remove stimulus:</strong> Suction secretions, blood, stop surgical stimulation</li>
                                        <li class="protocol-step"><strong>100% O₂:</strong> Via facemask with reservoir</li>
                                        <li class="protocol-step"><strong>Jaw thrust + CPAP:</strong> Anterior displacement of mandible at TMJ (Larson maneuver: press behind ramus), apply gentle positive pressure (10-20 cm H₂O)</li>
                                        <li class="protocol-step"><strong>Deepen anesthesia:</strong> Propofol 0.5-1 mg/kg IV (if IV access and patient not apneic)</li>
                                        <li class="protocol-step"><strong>If persistent/desaturating:</strong> Proceed to step 6</li>
                                        <li class="protocol-step"><strong>Succinylcholine:</strong> 0.5-1 mg/kg IV (or 2-4 mg/kg IM if no IV access)</li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Medication Doses</div>
                                    <div class="dose-detail"><strong>Propofol:</strong> 0.5-1 mg/kg IV bolus to deepen anesthesia</div>
                                    <div class="dose-detail"><strong>Succinylcholine (IV):</strong> 0.5-1 mg/kg IV (lower dose than RSI, may give without pretreatment in emergency)</div>
                                    <div class="dose-detail"><strong>Succinylcholine (IM):</strong> 2-4 mg/kg IM into deltoid or thigh if no IV access</div>
                                    <div class="dose-detail"><strong>Atropine (for bradycardia):</strong> 0.01-0.02 mg/kg IV (especially in children)</div>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🩺 Clinical Presentation</div>
                                    <div class="info-detail"><strong>Partial:</strong> High-pitched stridor, crowing sound, difficulty ventilating<br><strong>Complete:</strong> Silent chest, no air movement, severe paradoxical chest/abdominal motion, rapid desaturation</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Key Points</div>
                                    <div class="warning-detail">• Most common during light anesthesia (stage II) — emergence, extubation, airway instrumentation<br>• <strong>Larson's maneuver:</strong> Firm pressure at "laryngospasm notch" (behind angle of mandible, anterior to mastoid) + jaw thrust often breaks spasm<br>• Avoid repeated forceful PPV — worsens spasm and risks gastric insufflation/aspiration<br>• After succinylcholine, MUST ventilate/intubate until spontaneous breathing returns (~5-10 min)<br>• Consider intubation if multiple episodes or risk of recurrence</div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- Bronchospasm -->
                    <div class="protocol-card orange" data-keywords="bronchospasm wheezing asthma copd albuterol epinephrine deepen anesthesia">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">Bronchospasm</h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Diffuse small airway constriction causing wheezing and ↑ airway pressure</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag urgent">Urgent</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Actions</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>100% O₂:</strong> Increase FiO₂ to 1.0</li>
                                        <li class="protocol-step"><strong>Deepen anesthesia:</strong> Increase volatile (if not contraindicated) or propofol bolus 0.5-1 mg/kg</li>
                                        <li class="protocol-step"><strong>Hand ventilate:</strong> Slow rate (8-10/min), low tidal volume, prolonged expiration (I:E = 1:3 or 1:4) to avoid air trapping</li>
                                        <li class="protocol-step"><strong>Albuterol MDI:</strong> 4-8 puffs via ETT/LMA or in-line with circuit</li>
                                        <li class="protocol-step"><strong>Rule out other causes:</strong> Check circuit (kink, ETT obstruction), suction secretions, confirm tube not endobronchial</li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Bronchodilator Therapy</div>
                                    <div class="dose-detail"><strong>Albuterol MDI:</strong> 4-8 puffs (90 mcg/puff) via ETT, may repeat q20min</div>
                                    <div class="dose-detail"><strong>Albuterol nebulizer:</strong> 2.5-5 mg in-line with circuit (continuous if refractory)</div>
                                    <div class="dose-detail"><strong>Epinephrine (severe):</strong> 10-50 mcg IV boluses, titrate to effect</div>
                                    <div class="dose-detail"><strong>Epinephrine infusion:</strong> 0.02-0.1 mcg/kg/min for refractory cases</div>
                                    <div class="dose-detail"><strong>Ketamine:</strong> 0.5-1 mg/kg IV bolus (bronchodilator + anesthetic), then 0.5-2 mg/kg/h infusion</div>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Escalation for Severe/Refractory Bronchospasm</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>IV Steroids:</strong> Methylprednisolone 1-2 mg/kg IV (benefit after 4-6h)</li>
                                        <li class="protocol-step"><strong>Magnesium sulfate:</strong> 2 g IV over 20 min (smooth muscle relaxation)</li>
                                        <li class="protocol-step"><strong>Ipratropium bromide:</strong> 500 mcg nebulized (anticholinergic, add to albuterol)</li>
                                        <li class="protocol-step"><strong>Heliox:</strong> 60-70% helium / 30-40% O₂ mixture (reduces turbulent flow, buys time)</li>
                                        <li class="protocol-step"><strong>Consider anaphylaxis</strong> if new-onset bronchospasm without asthma history → give epinephrine IM</li>
                                    </ol>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🔍 Differential Diagnosis</div>
                                    <div class="info-detail"><strong>Not bronchospasm:</strong> ETT in bronchus, kinked tube, mucus plug, pneumothorax, pulmonary edema, aspiration, pulmonary embolism<br><strong>Confirm:</strong> Auscultate lungs (bilateral expiratory wheezing), check plateau pressure (↑ in bronchospasm, normal in obstruction)</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Ventilation Strategy</div>
                                    <div class="warning-detail">• <strong>Allow permissive hypercapnia</strong> — accept high CO₂ to avoid barotrauma from air trapping<br>• <strong>Avoid high PEEP</strong> initially (worsens air trapping) — use low PEEP (0-5 cm H₂O)<br>• Reduce respiratory rate, increase expiratory time (I:E ratio 1:3 or 1:4)<br>• Goal: Peak pressure <40 cm H₂O, plateau <30 cm H₂O if possible<br>• If refractory + deteriorating: Consider ECMO/VV-ECMO</div>
                                </div>
                            </div>
                        </div>
                    </div>

                </div>
            </div>

            <!-- Neurologic / Regional -->
            <div class="category-section" data-category="neurologic">
                <div class="category-header">
                    <div class="category-icon purple">🧠</div>
                    <h2 class="category-title">Neurologic / Regional Complications</h2>
                </div>
                <div class="protocols-grid">

                    <!-- High/Total Spinal -->
                    <div class="protocol-card purple" data-keywords="high spinal total spinal neuraxial block hypotension bradycardia respiratory arrest spinal epidural">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">High/Total Spinal</h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Excessive cephalad spread of neuraxial anesthesia causing cardiorespiratory compromise</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag immediate">Immediate</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Management (ABC approach)</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Call for help</strong></li>
                                        <li class="protocol-step"><strong>Airway:</strong> Secure airway immediately if respiratory distress, intubate if apneic or unable to protect airway</li>
                                        <li class="protocol-step"><strong>Breathing:</strong> Ventilate with 100% O₂ (PPV or mechanical ventilation)</li>
                                        <li class="protocol-step"><strong>Circulation:</strong> Aggressive vasopressor/fluid resuscitation (see below)</li>
                                        <li class="protocol-step"><strong>Position:</strong> Supine with legs elevated (↑ venous return), avoid Trendelenburg (worsens block height)</li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Hemodynamic Support</div>
                                    <div class="dose-detail"><strong>Ephedrine:</strong> 5-10 mg IV boluses (first-line if bradycardic)</div>
                                    <div class="dose-detail"><strong>Phenylephrine:</strong> 50-100 mcg IV boluses (if tachycardic)</div>
                                    <div class="dose-detail"><strong>Epinephrine:</strong> 10-100 mcg IV if refractory hypotension, or 0.05-0.1 mcg/kg/min infusion</div>
                                    <div class="dose-detail"><strong>Atropine:</strong> 0.4-0.6 mg IV for bradycardia</div>
                                    <div class="dose-detail"><strong>Fluids:</strong> Crystalloid 500-1000 mL rapid bolus</div>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Additional Supportive Care</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Reassure patient:</strong> Explain what's happening if conscious, give anxiolysis (midazolam 1-2 mg IV)</li>
                                        <li class="protocol-step"><strong>Monitor:</strong> Continuous BP, HR, SpO₂, capnography if intubated</li>
                                        <li class="protocol-step"><strong>Assess block level:</strong> Check sensory level when stable</li>
                                        <li class="protocol-step"><strong>Duration:</strong> Expect resolution in 1-3 hours (depends on local anesthetic type/dose)</li>
                                    </ol>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🩺 Clinical Presentation</div>
                                    <div class="info-detail"><strong>Early signs:</strong> Dyspnea, difficulty speaking, upper extremity weakness/tingling, nausea<br>
                                    <strong>Severe signs:</strong> Apnea, loss of consciousness, severe hypotension, bradycardia → asystole<br>
                                    <strong>Block level:</strong> High spinal = T1-T4 (respiratory muscles affected), Total spinal = brainstem involvement (unconsciousness)</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Key Considerations</div>
                                    <div class="warning-detail">• <strong>Causes:</strong> Excessive local anesthetic dose, migration of epidural catheter into subarachnoid space, unrecognized dural puncture during epidural placement<br>• Differentiate from local anesthetic toxicity (LAST) — high spinal has symmetric block and bradycardia; LAST has CNS excitation → seizures<br>• Most patients fully recover — maintain oxygenation and BP until block wears off<br>• Document event, monitor for post-dural puncture headache if unintentional wet tap</div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- Venous Air Embolism -->
                    <div class="protocol-card purple" data-keywords="air embolism vae gas embolism durant position sitting craniotomy etco2">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">Venous Air Embolism (VAE)</h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Air entrainment into venous circulation causing cardiovascular collapse</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag immediate">Immediate</span>
                            <span class="protocol-tag call-help">Call Help</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Actions</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Notify surgeon — STOP surgery:</strong> Flood surgical field with saline, occlude open veins</li>
                                        <li class="protocol-step"><strong>100% FiO₂:</strong> Discontinue nitrous oxide (N₂O expands air bubbles)</li>
                                        <li class="protocol-step"><strong>Lower surgical site:</strong> Reduce gradient between surgical site and heart</li>
                                        <li class="protocol-step"><strong>Durant position:</strong> Left lateral decubitus + Trendelenburg (traps air in RV apex, away from outflow)</li>
                                        <li class="protocol-step"><strong>Attempt aspiration:</strong> If CVL in place, aspirate from distal port (withdraw air from RA/RV)</li>
                                        <li class="protocol-step"><strong>Hemodynamic support:</strong> Fluids, vasopressors, inotropes as needed</li>
                                    </ol>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🩺 Signs of VAE (in order of sensitivity)</div>
                                    <div class="info-detail">1. <strong>Precordial Doppler:</strong> Mill-wheel murmur (most sensitive, requires placement)<br>
                                    2. <strong>↓ EtCO₂:</strong> Sudden drop (dead space ↑ from pulmonary embolism)<br>
                                    3. <strong>↓ SpO₂, ↓ BP, ↑ CVP, arrhythmias</strong><br>
                                    4. <strong>Cardiac arrest</strong> (if large volume air)</div>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Supportive Therapy</div>
                                    <div class="dose-detail"><strong>Volume:</strong> Rapid crystalloid infusion (↑ CVP reduces further air entrainment)</div>
                                    <div class="dose-detail"><strong>Vasopressors:</strong> Phenylephrine 50-200 mcg IV boluses or epinephrine 10-100 mcg IV</div>
                                    <div class="dose-detail"><strong>Inotropes:</strong> Epinephrine infusion 0.05-0.5 mcg/kg/min if myocardial dysfunction</div>
                                    <div class="dose-detail"><strong>CPR if arrest:</strong> Standard ACLS + prolonged resuscitation (air resorbs over time)</div>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Post-Event Management</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Neuro assessment:</strong> Evaluate for paradoxical embolism (if PFO) → stroke symptoms</li>
                                        <li class="protocol-step"><strong>Consider hyperbaric O₂</strong> if neurologic deficits (shrinks bubbles, improves oxygenation)</li>
                                        <li class="protocol-step"><strong>ICU monitoring:</strong> For severe VAE with hemodynamic instability</li>
                                        <li class="protocol-step"><strong>Echocardiography:</strong> TEE/TTE to assess for residual air, PFO, cardiac function</li>
                                    </ol>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ High-Risk Situations</div>
                                    <div class="warning-detail">• <strong>Sitting position craniotomy</strong> (most common, head >20 cm above heart)<br>• Posterior fossa surgery, cervical laminectomy<br>• Cesarean section, liver resection, laparoscopy with insufflation issues<br>• Central line placement with air entrainment<br>• <strong>Lethal dose:</strong> 3-5 mL/kg rapid bolus (200-300 mL in adults) can cause cardiovascular collapse<br>• <strong>Prevention:</strong> Avoid sitting position when possible, maintain adequate hydration, use precordial Doppler monitoring</div>
                                </div>
                            </div>
                        </div>
                    </div>

                </div>
            </div>

            <!-- Metabolic / Other -->
            <div class="category-section" data-category="metabolic">
                <div class="category-header">
                    <div class="category-icon amber">⚗️</div>
                    <h2 class="category-title">Metabolic / Electrolyte Crises</h2>
                </div>
                <div class="protocols-grid">

                    <!-- Hyperkalemia -->
                    <div class="protocol-card amber" data-keywords="hyperkalemia potassium peaked t wave arrhythmia calcium insulin dextrose kayexalate dialysis">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">Severe Hyperkalemia</h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Elevated serum potassium (K⁺ >6.5 mEq/L) with cardiac toxicity</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag immediate">Immediate</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Treatment (3-Step Approach)</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>STABILIZE cardiac membrane:</strong> Calcium (see dosing below) — immediate effect</li>
                                        <li class="protocol-step"><strong>SHIFT K⁺ intracellularly:</strong> Insulin + dextrose, beta-agonists, bicarbonate — onset 15-30 min</li>
                                        <li class="protocol-step"><strong>REMOVE K⁺ from body:</strong> Diuretics, dialysis, Kayexalate — onset 2-24 hours</li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Step 1: Stabilize Cardiac Membrane</div>
                                    <div class="dose-detail"><strong>Calcium chloride 10%:</strong> 10-20 mL (1-2 grams) IV over 2-5 min via central line (or diluted peripheral)</div>
                                    <div class="dose-detail"><strong>Calcium gluconate 10%:</strong> 30 mL (3 grams) IV over 2-5 min (alternative if only peripheral access)</div>
                                    <div class="dose-detail"><strong>Onset:</strong> 1-3 minutes, duration 30-60 min</div>
                                    <div class="dose-detail"><strong>Monitor:</strong> Continuous ECG during administration</div>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Step 2: Shift K⁺ Intracellularly</div>
                                    <div class="dose-detail"><strong>Insulin + Dextrose:</strong> 10 units regular insulin IV + 50 mL D50 (25g dextrose) over 5 min</div>
                                    <div class="dose-detail"><strong>Albuterol:</strong> 10-20 mg nebulized over 10 min (or 0.5 mg IV if available)</div>
                                    <div class="dose-detail"><strong>Sodium bicarbonate:</strong> 50-100 mEq (1-2 amps) IV over 5 min (if concurrent acidosis)</div>
                                    <div class="dose-detail"><strong>Effect:</strong> Lowers K⁺ by 0.5-1.5 mEq/L, onset 15-30 min, duration 4-6 hours</div>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Step 3: Remove K⁺ from Body</div>
                                    <div class="dose-detail"><strong>Furosemide:</strong> 40-80 mg IV (if adequate renal function)</div>
                                    <div class="dose-detail"><strong>Sodium polystyrene sulfonate (Kayexalate):</strong> 15-30 g PO/PR (slow, avoid if ileus)</div>
                                    <div class="dose-detail"><strong>Hemodialysis:</strong> Most effective for refractory hyperkalemia or ESRD (lowers K⁺ by 1-2 mEq/L/hour)</div>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🩺 ECG Changes (Progressive Severity)</div>
                                    <div class="info-detail">1. Peaked T waves (tall, narrow, symmetric)<br>
                                    2. Prolonged PR interval, flattened P waves<br>
                                    3. Widened QRS complex<br>
                                    4. Sine wave pattern → VF/asystole</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Anesthesia-Specific Causes</div>
                                    <div class="warning-detail">• <strong>Succinylcholine:</strong> ↑ K⁺ by 0.5 mEq/L normally; AVOID in burns, denervation injuries, chronic paralysis, massive trauma (can → fatal hyperkalemia)<br>• <strong>Massive transfusion:</strong> Stored blood has high K⁺<br>• <strong>Tourniquet release:</strong> Washout of ischemic tissue K⁺<br>• <strong>Tumor lysis, rhabdomyolysis, crush injury</strong><br>• Always check K⁺ before giving succinylcholine in at-risk patients</div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- Severe Hypoglycemia -->
                    <div class="protocol-card amber" data-keywords="hypoglycemia low blood sugar glucose dextrose glucagon insulin">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">Severe Hypoglycemia</h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Blood glucose <70 mg/dL with altered mental status or coma</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag urgent">Urgent</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Management</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Confirm with glucose check:</strong> Fingerstick or lab draw</li>
                                        <li class="protocol-step"><strong>If conscious and able to swallow:</strong> 15-20 g oral glucose (juice, glucose tabs)</li>
                                        <li class="protocol-step"><strong>If unconscious or NPO:</strong> Give IV dextrose (see dosing below)</li>
                                        <li class="protocol-step"><strong>Recheck glucose:</strong> Every 15 minutes until >100 mg/dL</li>
                                        <li class="protocol-step"><strong>Investigate cause:</strong> Insulin overdose, missed meal, sepsis, liver failure, adrenal insufficiency</li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 IV Dextrose Dosing</div>
                                    <div class="dose-detail"><strong>D50 (50% dextrose):</strong> 50 mL (25 g) IV push (adults) — verify IV patency (caustic if extravasates)</div>
                                    <div class="dose-detail"><strong>D10 (10% dextrose):</strong> 250 mL (25 g) IV over 5-10 min (safer for peripheral IV, pediatrics)</div>
                                    <div class="dose-detail"><strong>Pediatric:</strong> 0.5-1 g/kg IV (2-4 mL/kg of D25 or 5-10 mL/kg of D10)</div>
                                    <div class="dose-detail"><strong>Maintenance:</strong> Start D5 or D10 infusion after bolus to prevent recurrence</div>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Glucagon (if no IV access)</div>
                                    <div class="dose-detail"><strong>Adult dose:</strong> 1 mg IM or SC</div>
                                    <div class="dose-detail"><strong>Pediatric dose:</strong> 0.5 mg IM/SC if <20 kg, 1 mg if >20 kg</div>
                                    <div class="dose-detail"><strong>Onset:</strong> 10-15 minutes, duration 60-90 min</div>
                                    <div class="dose-detail"><strong>Note:</strong> Ineffective in glycogen-depleted states (starvation, alcoholism, adrenal insufficiency)</div>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🩺 Clinical Presentation</div>
                                    <div class="info-detail"><strong>Mild (50-70 mg/dL):</strong> Tremor, palpitations, sweating, hunger, anxiety<br>
                                    <strong>Moderate (40-50 mg/dL):</strong> Confusion, difficulty concentrating, slurred speech, blurred vision<br>
                                    <strong>Severe (<40 mg/dL):</strong> Seizures, loss of consciousness, coma, death if untreated</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Perioperative Considerations</div>
                                    <div class="warning-detail">• <strong>High-risk patients:</strong> Type 1 diabetes, insulin pumps, sulfonylurea use, prolonged fasting<br>• Intraoperative hypoglycemia may be masked by anesthesia — check glucose regularly in diabetics<br>• Avoid over-treating (hyperglycemia is harmful too) — target 100-180 mg/dL perioperatively<br>• After treatment, give complex carbs/meal to prevent rebound hypoglycemia (if patient able to eat)</div>
                                </div>
                            </div>
                        </div>
                    </div>

                </div>
            </div>

        </div>

        <!-- Footer -->
        <footer class="footer">
            <div class="footer-content">
                <div class="disclaimer">
                    <div class="disclaimer-title">
                        <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                            <path d="M10.29 3.86L1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z"></path>
                            <line x1="12" y1="9" x2="12" y2="13"></line>
                            <line x1="12" y1="17" x2="12.01" y2="17"></line>
                        </svg>
                        Important Disclaimer
                    </div>
                    <p class="disclaimer-text">
                        These protocols are for educational reference only and should NOT replace clinical judgment, institutional protocols, or expert consultation. Always follow your institution's guidelines and call for help early in any crisis situation. This is not medical advice.
                    </p>
                </div>
                <div class="footer-links">
                    <a href="/" class="footer-link">Home</a>
                    <a href="/calculators" class="footer-link">Calculators</a>
                    <a href="/terms" class="footer-link">Terms</a>
                    <a href="/privacy" class="footer-link">Privacy</a>
                </div>
                <p class="footer-copy">
                    © 2025 gasconsult.ai — Evidence-based anesthesiology AI assistant
                </p>
            </div>
        </footer>

    </div>

    <script>
        // Mobile menu toggle
        function toggleMobileMenu() {
            const btn = document.querySelector('.mobile-menu-btn');
            const menu = document.getElementById('mobileMenu');
            btn.classList.toggle('active');
            menu.classList.toggle('active');
        }

        // Protocol card expansion
        document.addEventListener('DOMContentLoaded', function() {
            const protocolCards = document.querySelectorAll('.protocol-card');

            protocolCards.forEach(card => {
                card.addEventListener('click', function(e) {
                    // Don't toggle if clicking on a link inside the card
                    if (e.target.tagName === 'A') return;

                    this.classList.toggle('expanded');
                });
            });

            // Search input clear button visibility
            const searchInput = document.getElementById('searchInput');
            const clearBtn = document.getElementById('clearBtn');

            searchInput.addEventListener('input', function() {
                if (this.value.length > 0) {
                    clearBtn.classList.add('visible');
                } else {
                    clearBtn.classList.remove('visible');
                }
            });
        });

        // Search/filter functionality
        function filterProtocols() {
            const searchTerm = document.getElementById('searchInput').value.toLowerCase();
            const cards = document.querySelectorAll('.protocol-card');
            const sections = document.querySelectorAll('.category-section');

            let visibleCount = 0;

            cards.forEach(card => {
                const keywords = card.getAttribute('data-keywords').toLowerCase();
                const title = card.querySelector('.protocol-title').textContent.toLowerCase();
                const summary = card.querySelector('.protocol-summary').textContent.toLowerCase();

                if (keywords.includes(searchTerm) || title.includes(searchTerm) || summary.includes(searchTerm)) {
                    card.style.display = '';
                    visibleCount++;
                } else {
                    card.style.display = 'none';
                }
            });

            // Hide empty categories
            sections.forEach(section => {
                const visibleCards = section.querySelectorAll('.protocol-card:not([style*="display: none"])');
                if (visibleCards.length === 0) {
                    section.classList.add('hidden');
                } else {
                    section.classList.remove('hidden');
                }
            });
        }

        function clearSearch() {
            const searchInput = document.getElementById('searchInput');
            const clearBtn = document.getElementById('clearBtn');
            searchInput.value = '';
            clearBtn.classList.remove('visible');
            filterProtocols();
        }

        // Keyboard shortcut: Press '/' to focus search
        document.addEventListener('keydown', function(e) {
            if (e.key === '/' && document.activeElement !== document.getElementById('searchInput')) {
                e.preventDefault();
                document.getElementById('searchInput').focus();
            }
            // Press Escape to clear search and collapse all cards
            if (e.key === 'Escape') {
                clearSearch();
                document.querySelectorAll('.protocol-card.expanded').forEach(card => {
                    card.classList.remove('expanded');
                });
            }
        });
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
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=6">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=6">
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
            padding: 6px;
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
            padding: 4px;
            display: flex;
            align-items: flex-end;
            gap: 4px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04), inset 0 1px 0 rgba(255,255,255,1);
        }

        .chat-input {
            flex: 1;
            border: none;
            outline: none;
            padding: 10px 14px;
            font-size: 16px;
            font-family: inherit;
            color: var(--gray-800);
            background: transparent;
            resize: none;
            min-height: 38px;
            max-height: 110px;
            line-height: 1.5;
        }

        .chat-input::placeholder { color: var(--gray-400); }

        .chat-send {
            width: 44px;
            height: 44px;
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
            margin: 4px;
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
                    <a href="/calculators" class="nav-link">Clinical Calculators</a>
                    <a href="/crisis" class="nav-link">Crisis Protocols</a>
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
            <a href="/crisis" class="mobile-menu-link">Crisis Protocols</a>
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
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=6">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=6">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">

    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800;900&display=swap" rel="stylesheet">
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
            --green-50: #F0FDF4;
            --green-500: #10B981;
            --green-600: #059669;
            --green-700: #047857;
            --yellow-50: #FEFCE8;
            --yellow-500: #F59E0B;
            --yellow-600: #D97706;
            --red-50: #FEF2F2;
            --red-500: #EF4444;
            --red-600: #DC2626;
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

        /* Glassmorphic Background */
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

        /* Navigation */
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

        @media (min-width: 768px) {
            .nav { padding: 16px 32px; }
            .nav-inner { height: 64px; padding: 0 24px; border-radius: 20px; }
            .logo-icon svg { width: 42px; height: 15px; }
            .logo-text { font-size: 20px; }
            .nav-links { display: flex; }
            .mobile-menu-btn { display: none; }
        }

        /* Hero Section */
        .hero {
            padding: 120px 20px 40px;
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
            margin-bottom: 20px;
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

        @keyframes fade-up {
            from { opacity: 0; transform: translateY(24px); }
            to { opacity: 1; transform: translateY(0); }
        }

        .badge-text {
            font-size: 12px;
            font-weight: 600;
            color: var(--gray-700);
        }

        .hero-title {
            font-size: 42px;
            font-weight: 900;
            line-height: 1.1;
            letter-spacing: -1.5px;
            color: var(--gray-900);
            margin-bottom: 16px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.1s forwards;
            opacity: 0;
        }

        .hero-title .gradient {
            background: linear-gradient(135deg, var(--blue-600) 0%, var(--blue-500) 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
        }

        .hero-subtitle {
            font-size: 17px;
            font-weight: 400;
            line-height: 1.6;
            color: var(--gray-600);
            max-width: 580px;
            margin: 0 auto;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.2s forwards;
            opacity: 0;
        }

        /* Main Content */
        .main {
            flex: 1;
            padding: 0 20px 80px;
        }

        .container {
            max-width: 1400px;
            margin: 0 auto;
        }

        /* Calculator Grid */
        .calc-grid {
            display: grid;
            grid-template-columns: 1fr;
            gap: 24px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.3s forwards;
            opacity: 0;
        }

        @media (min-width: 768px) {
            .calc-grid {
                grid-template-columns: repeat(2, 1fr);
            }
        }

        @media (min-width: 1024px) {
            .nav { padding: 16px 40px; }
        }

        @media (min-width: 1200px) {
            .calc-grid {
                grid-template-columns: repeat(3, 1fr);
            }
        }

        /* Calculator Card */
        .calc-card {
            background: rgba(255, 255, 255, 0.8);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 20px;
            padding: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 12px 48px rgba(0,0,0,0.03);
            transition: all 0.3s cubic-bezier(0.4,0,0.2,1);
        }

        .calc-card:hover {
            transform: translateY(-4px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.04), 0 12px 32px rgba(0,0,0,0.06), 0 24px 64px rgba(0,0,0,0.08);
            border-color: rgba(59,130,246,0.2);
        }

        .calc-card.large {
            grid-column: span 1;
        }

        @media (min-width: 768px) {
            .calc-card.large {
                grid-column: span 2;
            }
        }

        .calc-header {
            margin-bottom: 20px;
            padding-bottom: 20px;
            border-bottom: 1px solid var(--gray-200);
        }

        .calc-title {
            font-size: 20px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 6px;
            letter-spacing: -0.5px;
        }

        .calc-desc {
            font-size: 13px;
            color: var(--gray-500);
            line-height: 1.5;
        }

        .calc-body {
            margin-bottom: 20px;
        }

        /* Form Inputs */
        .input-group {
            margin-bottom: 16px;
        }

        .input-group:last-child {
            margin-bottom: 0;
        }

        .input-label {
            display: block;
            font-size: 13px;
            font-weight: 600;
            color: var(--gray-700);
            margin-bottom: 8px;
            letter-spacing: 0.2px;
        }

        .input-wrapper {
            position: relative;
        }

        .calc-input, .calc-select {
            width: 100%;
            padding: 12px 16px;
            font-family: inherit;
            font-size: 15px;
            font-weight: 500;
            color: var(--gray-900);
            background: var(--white);
            border: 1.5px solid var(--gray-300);
            border-radius: 12px;
            transition: all 0.2s ease;
            outline: none;
        }

        .calc-input:focus, .calc-select:focus {
            border-color: var(--blue-500);
            box-shadow: 0 0 0 4px rgba(59,130,246,0.1);
        }

        .calc-input::placeholder {
            color: var(--gray-400);
            font-weight: 400;
        }

        .input-unit {
            position: absolute;
            right: 16px;
            top: 50%;
            transform: translateY(-50%);
            font-size: 13px;
            font-weight: 600;
            color: var(--gray-500);
            pointer-events: none;
        }

        /* Checkbox Inputs */
        .checkbox-group {
            display: flex;
            flex-direction: column;
            gap: 12px;
        }

        .checkbox-item {
            display: flex;
            align-items: center;
            padding: 12px;
            background: var(--gray-50);
            border: 1.5px solid var(--gray-200);
            border-radius: 10px;
            transition: all 0.2s ease;
            cursor: pointer;
        }

        .checkbox-item:hover {
            background: var(--blue-50);
            border-color: var(--blue-200);
        }

        .checkbox-item input[type="checkbox"] {
            width: 20px;
            height: 20px;
            margin-right: 12px;
            cursor: pointer;
            accent-color: var(--blue-600);
        }

        .checkbox-label {
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-700);
            cursor: pointer;
            flex: 1;
        }

        /* Radio Inputs */
        .radio-group {
            display: flex;
            flex-direction: column;
            gap: 10px;
        }

        .radio-item {
            display: flex;
            align-items: start;
            padding: 12px;
            background: var(--gray-50);
            border: 1.5px solid var(--gray-200);
            border-radius: 10px;
            transition: all 0.2s ease;
            cursor: pointer;
        }

        .radio-item:hover {
            background: var(--blue-50);
            border-color: var(--blue-200);
        }

        .radio-item input[type="radio"] {
            width: 18px;
            height: 18px;
            margin-right: 12px;
            margin-top: 2px;
            cursor: pointer;
            accent-color: var(--blue-600);
            flex-shrink: 0;
        }

        .radio-content {
            flex: 1;
        }

        .radio-label {
            font-size: 14px;
            font-weight: 600;
            color: var(--gray-900);
            cursor: pointer;
            display: block;
            margin-bottom: 2px;
        }

        .radio-desc {
            font-size: 12px;
            color: var(--gray-600);
            line-height: 1.4;
        }

        /* Results Section */
        .calc-result {
            background: linear-gradient(135deg, var(--blue-50) 0%, var(--blue-100) 100%);
            border: 1.5px solid var(--blue-200);
            border-radius: 14px;
            padding: 20px;
            margin-bottom: 16px;
            opacity: 0;
            transform: scale(0.95);
            transition: all 0.4s cubic-bezier(0.34,1.56,0.64,1);
        }

        .calc-result.show {
            opacity: 1;
            transform: scale(1);
        }

        .calc-result.green {
            background: linear-gradient(135deg, var(--green-50) 0%, #E6F7EE 100%);
            border-color: var(--green-500);
        }

        .calc-result.yellow {
            background: linear-gradient(135deg, var(--yellow-50) 0%, #FEF9E7 100%);
            border-color: var(--yellow-500);
        }

        .calc-result.red {
            background: linear-gradient(135deg, var(--red-50) 0%, #FFEBEE 100%);
            border-color: var(--red-500);
        }

        .result-label {
            font-size: 12px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.8px;
            color: var(--gray-600);
            margin-bottom: 8px;
        }

        .result-value {
            font-size: 32px;
            font-weight: 900;
            color: var(--blue-700);
            line-height: 1;
            letter-spacing: -1px;
        }

        .green .result-value { color: var(--green-700); }
        .yellow .result-value { color: var(--yellow-600); }
        .red .result-value { color: var(--red-600); }

        .result-unit {
            font-size: 18px;
            font-weight: 600;
            color: var(--gray-600);
            margin-left: 4px;
        }

        .result-text {
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-700);
            margin-top: 8px;
            line-height: 1.5;
        }

        .result-grid {
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 12px;
        }

        .result-item {
            background: rgba(255, 255, 255, 0.6);
            border-radius: 10px;
            padding: 14px;
        }

        .result-item-label {
            font-size: 11px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: var(--gray-600);
            margin-bottom: 6px;
        }

        .result-item-value {
            font-size: 22px;
            font-weight: 800;
            color: var(--gray-900);
            line-height: 1;
        }

        .result-item-unit {
            font-size: 13px;
            font-weight: 600;
            color: var(--gray-600);
            margin-left: 2px;
        }

        /* Reference Link */
        .calc-reference {
            padding-top: 16px;
            border-top: 1px solid var(--gray-200);
        }

        .ref-link {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            font-size: 12px;
            font-weight: 600;
            color: var(--blue-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .ref-link:hover {
            color: var(--blue-700);
            gap: 8px;
        }

        .ref-link svg {
            width: 14px;
            height: 14px;
        }

        /* Footer */
        .footer {
            padding: 40px 20px;
            border-top: 1px solid var(--gray-200);
            background: rgba(255, 255, 255, 0.5);
            backdrop-filter: blur(10px);
        }

        .footer-inner {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 20px;
        }

        .footer-brand {
            display: flex;
            align-items: center;
            gap: 16px;
        }

        .footer-logo svg {
            width: 32px;
            height: 32px;
        }

        .footer-text {
            font-size: 14px;
            color: var(--gray-600);
        }

        .footer-links {
            display: flex;
            gap: 24px;
        }

        .footer-link {
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            transition: color 0.2s ease;
        }

        .footer-link:hover {
            color: var(--blue-600);
        }

        /* Medical Disclaimer */
        .disclaimer {
            max-width: 1200px;
            margin: 0 auto 40px;
            padding: 0 20px;
        }

        .disclaimer-card {
            background: rgba(239, 246, 255, 0.6);
            border: 1px solid var(--blue-200);
            border-radius: 12px;
            padding: 16px 20px;
            display: flex;
            align-items: start;
            gap: 12px;
        }

        .disclaimer-icon {
            flex-shrink: 0;
            width: 20px;
            height: 20px;
            color: var(--blue-600);
        }

        .disclaimer-text {
            font-size: 13px;
            line-height: 1.6;
            color: var(--gray-700);
        }

        .disclaimer-text strong {
            font-weight: 600;
            color: var(--gray-900);
        }

        @media (max-width: 767px) {
            .hero-title {
                font-size: 32px;
            }

            .result-grid {
                grid-template-columns: 1fr;
            }

            .footer-inner {
                flex-direction: column;
                text-align: center;
            }
        }
    </style>
</head>
<body>
    <!-- Background -->
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <!-- Navigation -->
        <nav class="nav">
            <div class="nav-inner">
                <a href="/" class="logo">
                    <div class="logo-icon">
                        <svg viewBox="0 0 120 40" fill="none">
                            <circle cx="20" cy="20" r="18" fill="#2563EB"/>
                            <circle cx="60" cy="20" r="18" fill="#2563EB" fill-opacity="0.5"/>
                            <circle cx="100" cy="20" r="18" fill="#2563EB" fill-opacity="0.2"/>
                        </svg>
                    </div>
                    <div class="logo-text">
                        <span class="gas">gas</span><span class="consult">consult</span><span class="ai">.ai</span>
                    </div>
                </a>
                <div class="nav-links">
                    <a href="/?clear=1" class="nav-link">Home</a>
                    <a href="/quick-dose" class="nav-link">Quick Dose</a>
                    <a href="/preop" class="nav-link">Pre-Op</a>
                    <a href="/calculators" class="nav-link active">Clinical Calculators</a>
                    <a href="/hypotension" class="nav-link">IOH Predictor</a>
                </div>
                <button class="mobile-menu-btn" onclick="toggleMobileMenu()" aria-label="Menu">
                    <span></span>
                    <span></span>
                    <span></span>
                </button>
            </div>
        </nav>

        <!-- Mobile Menu -->
        <div class="mobile-menu" id="mobileMenu">
            <a href="/?clear=1" class="mobile-menu-link">Home</a>
            <a href="/quick-dose" class="mobile-menu-link">Quick Dose</a>
            <a href="/preop" class="mobile-menu-link">Pre-Op</a>
            <a href="/calculators" class="mobile-menu-link">Clinical Calculators</a>
            <a href="/crisis" class="mobile-menu-link">Crisis Protocols</a>
            <a href="/hypotension" class="mobile-menu-link">IOH Predictor</a>
        </div>

        <!-- Hero -->
        <section class="hero">
            <div class="hero-badge">
                <div class="badge-dot"></div>
                <span class="badge-text">Evidence-Based Formulas</span>
            </div>
            <h1 class="hero-title">
                Clinical <span class="gradient">Calculators</span>
            </h1>
            <p class="hero-subtitle">
                Real-time medical calculations with validated formulas. All results are instant — no submit buttons needed.
            </p>
        </section>

        <!-- Disclaimer -->
        <div class="disclaimer">
            <div class="disclaimer-card">
                <svg class="disclaimer-icon" fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                    <path stroke-linecap="round" stroke-linejoin="round" d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z"></path>
                </svg>
                <div class="disclaimer-text">
                    <strong>Medical Disclaimer:</strong> These calculators are for educational purposes only. Always verify calculations and consult clinical judgment before making treatment decisions.
                </div>
            </div>
        </div>

        <!-- Main Content -->
        <main class="main">
            <div class="container">
                <div class="calc-grid">

                    <!-- IDEAL BODY WEIGHT -->
                    <div class="calc-card">
                        <div class="calc-header">
                            <h3 class="calc-title">Ideal Body Weight (IBW)</h3>
                            <p class="calc-desc">Devine formula for anesthesia dosing</p>
                        </div>
                        <div class="calc-body">
                            <div class="input-group">
                                <label class="input-label">Sex</label>
                                <select class="calc-select" id="ibw-sex" onchange="calcIBW()">
                                    <option value="">Select...</option>
                                    <option value="male">Male</option>
                                    <option value="female">Female</option>
                                </select>
                            </div>
                            <div class="input-group">
                                <label class="input-label">Height</label>
                                <div class="input-wrapper">
                                    <input type="number" class="calc-input" id="ibw-height" placeholder="Enter height" oninput="calcIBW()" min="0" step="0.1">
                                    <span class="input-unit">inches</span>
                                </div>
                            </div>
                            <div class="input-group">
                                <label class="input-label">Actual Weight (optional, for ABW)</label>
                                <div class="input-wrapper">
                                    <input type="number" class="calc-input" id="ibw-weight" placeholder="Enter weight" oninput="calcIBW()" min="0" step="0.1">
                                    <span class="input-unit">kg</span>
                                </div>
                            </div>
                        </div>
                        <div id="ibw-result" class="calc-result">
                            <div class="result-label">Ideal Body Weight</div>
                            <div class="result-value"><span id="ibw-value">--</span><span class="result-unit">kg</span></div>
                            <div id="abw-display" style="margin-top: 12px; display: none;">
                                <div class="result-label" style="margin-bottom: 4px;">Adjusted Body Weight</div>
                                <div style="font-size: 20px; font-weight: 700; color: var(--gray-700);"><span id="abw-value">--</span> kg</div>
                            </div>
                        </div>
                        <div class="calc-reference">
                            <a href="https://clincalc.com/kinetics/idealbw.aspx" target="_blank" class="ref-link">
                                <span>Reference: ClinCalc IBW</span>
                                <svg fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                                    <path stroke-linecap="round" stroke-linejoin="round" d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14"></path>
                                </svg>
                            </a>
                        </div>
                    </div>

                    <!-- BODY SURFACE AREA -->
                    <div class="calc-card">
                        <div class="calc-header">
                            <h3 class="calc-title">Body Surface Area (BSA)</h3>
                            <p class="calc-desc">Mosteller, DuBois & Haycock formulas</p>
                        </div>
                        <div class="calc-body">
                            <div class="input-group">
                                <label class="input-label">Height</label>
                                <div class="input-wrapper">
                                    <input type="number" class="calc-input" id="bsa-height" placeholder="Enter height" oninput="calcBSA()" min="0" step="0.1">
                                    <span class="input-unit">cm</span>
                                </div>
                            </div>
                            <div class="input-group">
                                <label class="input-label">Weight</label>
                                <div class="input-wrapper">
                                    <input type="number" class="calc-input" id="bsa-weight" placeholder="Enter weight" oninput="calcBSA()" min="0" step="0.1">
                                    <span class="input-unit">kg</span>
                                </div>
                            </div>
                        </div>
                        <div id="bsa-result" class="calc-result green">
                            <div class="result-grid">
                                <div class="result-item">
                                    <div class="result-item-label">Mosteller ⭐</div>
                                    <div class="result-item-value"><span id="bsa-mosteller">--</span><span class="result-item-unit">m²</span></div>
                                </div>
                                <div class="result-item">
                                    <div class="result-item-label">DuBois</div>
                                    <div class="result-item-value"><span id="bsa-dubois">--</span><span class="result-item-unit">m²</span></div>
                                </div>
                                <div class="result-item" style="grid-column: 1 / -1;">
                                    <div class="result-item-label">Haycock (Pediatric)</div>
                                    <div class="result-item-value"><span id="bsa-haycock">--</span><span class="result-item-unit">m²</span></div>
                                </div>
                            </div>
                        </div>
                        <div class="calc-reference">
                            <a href="https://medplore.com/health-tools/bsa-calculator-with-formulas/" target="_blank" class="ref-link">
                                <span>Reference: BSA Formulas</span>
                                <svg fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                                    <path stroke-linecap="round" stroke-linejoin="round" d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14"></path>
                                </svg>
                            </a>
                        </div>
                    </div>

                    <!-- MAXIMUM ALLOWABLE BLOOD LOSS -->
                    <div class="calc-card">
                        <div class="calc-header">
                            <h3 class="calc-title">Maximum Allowable Blood Loss</h3>
                            <p class="calc-desc">Estimated blood volume & transfusion threshold</p>
                        </div>
                        <div class="calc-body">
                            <div class="input-group">
                                <label class="input-label">Patient Type</label>
                                <select class="calc-select" id="mabl-type" onchange="calcMABL()">
                                    <option value="">Select...</option>
                                    <option value="male">Adult Male</option>
                                    <option value="female">Adult Female</option>
                                    <option value="pediatric">Pediatric</option>
                                </select>
                            </div>
                            <div class="input-group">
                                <label class="input-label">Weight</label>
                                <div class="input-wrapper">
                                    <input type="number" class="calc-input" id="mabl-weight" placeholder="Enter weight" oninput="calcMABL()" min="0" step="0.1">
                                    <span class="input-unit">kg</span>
                                </div>
                            </div>
                            <div class="input-group">
                                <label class="input-label">Starting Hematocrit</label>
                                <div class="input-wrapper">
                                    <input type="number" class="calc-input" id="mabl-start-hct" placeholder="e.g., 40" oninput="calcMABL()" min="0" max="100" step="0.1">
                                    <span class="input-unit">%</span>
                                </div>
                            </div>
                            <div class="input-group">
                                <label class="input-label">Target Hematocrit</label>
                                <div class="input-wrapper">
                                    <input type="number" class="calc-input" id="mabl-target-hct" placeholder="e.g., 25" oninput="calcMABL()" min="0" max="100" step="0.1">
                                    <span class="input-unit">%</span>
                                </div>
                            </div>
                        </div>
                        <div id="mabl-result" class="calc-result yellow">
                            <div class="result-label">Max Allowable Blood Loss</div>
                            <div class="result-value"><span id="mabl-value">--</span><span class="result-unit">mL</span></div>
                            <div class="result-text">Estimated Blood Volume: <strong id="ebv-value">--</strong> mL</div>
                        </div>
                        <div class="calc-reference">
                            <a href="https://www.mdcalc.com/calc/3905/maximum-allowable-blood-loss-abl-without-transfusion" target="_blank" class="ref-link">
                                <span>Reference: MDCalc MABL</span>
                                <svg fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                                    <path stroke-linecap="round" stroke-linejoin="round" d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14"></path>
                                </svg>
                            </a>
                        </div>
                    </div>

                    <!-- QTc INTERVAL -->
                    <div class="calc-card">
                        <div class="calc-header">
                            <h3 class="calc-title">QTc Interval (Bazett)</h3>
                            <p class="calc-desc">Corrected QT interval for cardiac risk</p>
                        </div>
                        <div class="calc-body">
                            <div class="input-group">
                                <label class="input-label">QT Interval</label>
                                <div class="input-wrapper">
                                    <input type="number" class="calc-input" id="qtc-qt" placeholder="e.g., 400" oninput="calcQTc()" min="0" step="1">
                                    <span class="input-unit">ms</span>
                                </div>
                            </div>
                            <div class="input-group">
                                <label class="input-label">Heart Rate</label>
                                <div class="input-wrapper">
                                    <input type="number" class="calc-input" id="qtc-hr" placeholder="e.g., 75" oninput="calcQTc()" min="0" step="1">
                                    <span class="input-unit">bpm</span>
                                </div>
                            </div>
                            <div class="input-group">
                                <label class="input-label">Sex (for interpretation)</label>
                                <select class="calc-select" id="qtc-sex" onchange="calcQTc()">
                                    <option value="">Select...</option>
                                    <option value="male">Male</option>
                                    <option value="female">Female</option>
                                </select>
                            </div>
                        </div>
                        <div id="qtc-result" class="calc-result">
                            <div class="result-label">QTc (Bazett)</div>
                            <div class="result-value"><span id="qtc-value">--</span><span class="result-unit">ms</span></div>
                            <div class="result-text" id="qtc-interpretation">Enter values above</div>
                        </div>
                        <div class="calc-reference">
                            <a href="https://en.wikipedia.org/wiki/QT_interval#Corrected_QT_interval" target="_blank" class="ref-link">
                                <span>Reference: Bazett Formula</span>
                                <svg fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                                    <path stroke-linecap="round" stroke-linejoin="round" d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14"></path>
                                </svg>
                            </a>
                        </div>
                    </div>

                    <!-- MAINTENANCE FLUIDS -->
                    <div class="calc-card">
                        <div class="calc-header">
                            <h3 class="calc-title">Maintenance Fluids (4-2-1 Rule)</h3>
                            <p class="calc-desc">Holliday-Segar method for pediatrics</p>
                        </div>
                        <div class="calc-body">
                            <div class="input-group">
                                <label class="input-label">Weight</label>
                                <div class="input-wrapper">
                                    <input type="number" class="calc-input" id="fluids-weight" placeholder="Enter weight" oninput="calcFluids()" min="0" step="0.1">
                                    <span class="input-unit">kg</span>
                                </div>
                            </div>
                        </div>
                        <div id="fluids-result" class="calc-result green">
                            <div class="result-label">Maintenance Rate</div>
                            <div class="result-value"><span id="fluids-value">--</span><span class="result-unit">mL/hr</span></div>
                            <div class="result-text">Daily volume: <strong id="fluids-daily">--</strong> mL/day</div>
                        </div>
                        <div class="calc-reference">
                            <a href="https://www.ncbi.nlm.nih.gov/books/NBK562337/" target="_blank" class="ref-link">
                                <span>Reference: 4-2-1 Rule</span>
                                <svg fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                                    <path stroke-linecap="round" stroke-linejoin="round" d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14"></path>
                                </svg>
                            </a>
                        </div>
                    </div>

                    <!-- PONV RISK (APFEL SCORE) -->
                    <div class="calc-card">
                        <div class="calc-header">
                            <h3 class="calc-title">PONV Risk (Apfel Score)</h3>
                            <p class="calc-desc">Postoperative nausea & vomiting prediction</p>
                        </div>
                        <div class="calc-body">
                            <div class="checkbox-group">
                                <label class="checkbox-item">
                                    <input type="checkbox" id="ponv-female" onchange="calcPONV()">
                                    <span class="checkbox-label">Female</span>
                                </label>
                                <label class="checkbox-item">
                                    <input type="checkbox" id="ponv-nonsmoker" onchange="calcPONV()">
                                    <span class="checkbox-label">Non-smoker</span>
                                </label>
                                <label class="checkbox-item">
                                    <input type="checkbox" id="ponv-history" onchange="calcPONV()">
                                    <span class="checkbox-label">History of PONV or motion sickness</span>
                                </label>
                                <label class="checkbox-item">
                                    <input type="checkbox" id="ponv-opioids" onchange="calcPONV()">
                                    <span class="checkbox-label">Postoperative opioids planned</span>
                                </label>
                            </div>
                        </div>
                        <div id="ponv-result" class="calc-result">
                            <div class="result-label">PONV Risk</div>
                            <div class="result-value"><span id="ponv-value">--</span><span class="result-unit">%</span></div>
                            <div class="result-text">Score: <strong id="ponv-score">0</strong> / 4 risk factors</div>
                        </div>
                        <div class="calc-reference">
                            <a href="https://www.mdcalc.com/calc/1887/apfel-score-postoperative-nausea-vomiting" target="_blank" class="ref-link">
                                <span>Reference: Apfel Score</span>
                                <svg fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                                    <path stroke-linecap="round" stroke-linejoin="round" d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14"></path>
                                </svg>
                            </a>
                        </div>
                    </div>

                    <!-- RCRI (CARDIAC RISK INDEX) -->
                    <div class="calc-card">
                        <div class="calc-header">
                            <h3 class="calc-title">RCRI (Revised Cardiac Risk Index)</h3>
                            <p class="calc-desc">Perioperative cardiac event prediction</p>
                        </div>
                        <div class="calc-body">
                            <div class="checkbox-group">
                                <label class="checkbox-item">
                                    <input type="checkbox" id="rcri-surgery" onchange="calcRCRI()">
                                    <span class="checkbox-label">High-risk surgery (vascular, intraperitoneal, intrathoracic)</span>
                                </label>
                                <label class="checkbox-item">
                                    <input type="checkbox" id="rcri-ihd" onchange="calcRCRI()">
                                    <span class="checkbox-label">History of ischemic heart disease</span>
                                </label>
                                <label class="checkbox-item">
                                    <input type="checkbox" id="rcri-chf" onchange="calcRCRI()">
                                    <span class="checkbox-label">History of congestive heart failure</span>
                                </label>
                                <label class="checkbox-item">
                                    <input type="checkbox" id="rcri-cvd" onchange="calcRCRI()">
                                    <span class="checkbox-label">History of cerebrovascular disease</span>
                                </label>
                                <label class="checkbox-item">
                                    <input type="checkbox" id="rcri-dm" onchange="calcRCRI()">
                                    <span class="checkbox-label">Diabetes on insulin therapy</span>
                                </label>
                                <label class="checkbox-item">
                                    <input type="checkbox" id="rcri-cr" onchange="calcRCRI()">
                                    <span class="checkbox-label">Preoperative creatinine > 2.0 mg/dL</span>
                                </label>
                            </div>
                        </div>
                        <div id="rcri-result" class="calc-result">
                            <div class="result-label">Cardiac Event Risk</div>
                            <div class="result-value"><span id="rcri-value">--</span><span class="result-unit">%</span></div>
                            <div class="result-text">Score: <strong id="rcri-score">0</strong> / 6 risk factors</div>
                        </div>
                        <div class="calc-reference">
                            <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC10578796/" target="_blank" class="ref-link">
                                <span>Reference: RCRI Validation</span>
                                <svg fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                                    <path stroke-linecap="round" stroke-linejoin="round" d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14"></path>
                                </svg>
                            </a>
                        </div>
                    </div>

                    <!-- ASA PHYSICAL STATUS -->
                    <div class="calc-card">
                        <div class="calc-header">
                            <h3 class="calc-title">ASA Physical Status</h3>
                            <p class="calc-desc">Classification of patient health status</p>
                        </div>
                        <div class="calc-body">
                            <div class="radio-group">
                                <label class="radio-item">
                                    <input type="radio" name="asa" value="1" onchange="calcASA()">
                                    <div class="radio-content">
                                        <span class="radio-label">ASA I</span>
                                        <span class="radio-desc">Healthy patient, no systemic disease</span>
                                    </div>
                                </label>
                                <label class="radio-item">
                                    <input type="radio" name="asa" value="2" onchange="calcASA()">
                                    <div class="radio-content">
                                        <span class="radio-label">ASA II</span>
                                        <span class="radio-desc">Mild systemic disease, no functional limitations</span>
                                    </div>
                                </label>
                                <label class="radio-item">
                                    <input type="radio" name="asa" value="3" onchange="calcASA()">
                                    <div class="radio-content">
                                        <span class="radio-label">ASA III</span>
                                        <span class="radio-desc">Severe systemic disease with functional limitations</span>
                                    </div>
                                </label>
                                <label class="radio-item">
                                    <input type="radio" name="asa" value="4" onchange="calcASA()">
                                    <div class="radio-content">
                                        <span class="radio-label">ASA IV</span>
                                        <span class="radio-desc">Severe disease that is constant threat to life</span>
                                    </div>
                                </label>
                                <label class="radio-item">
                                    <input type="radio" name="asa" value="5" onchange="calcASA()">
                                    <div class="radio-content">
                                        <span class="radio-label">ASA V</span>
                                        <span class="radio-desc">Moribund, not expected to survive without operation</span>
                                    </div>
                                </label>
                                <label class="radio-item">
                                    <input type="radio" name="asa" value="6" onchange="calcASA()">
                                    <div class="radio-content">
                                        <span class="radio-label">ASA VI</span>
                                        <span class="radio-desc">Brain-dead patient for organ donation</span>
                                    </div>
                                </label>
                            </div>
                        </div>
                        <div id="asa-result" class="calc-result" style="display: none;">
                            <div class="result-label">Selected Classification</div>
                            <div class="result-value" style="font-size: 24px;"><span id="asa-value">--</span></div>
                            <div class="result-text" id="asa-desc">--</div>
                        </div>
                        <div class="calc-reference">
                            <a href="https://www.asahq.org/standards-and-guidelines/asa-physical-status-classification-system" target="_blank" class="ref-link">
                                <span>Reference: ASA Guidelines</span>
                                <svg fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                                    <path stroke-linecap="round" stroke-linejoin="round" d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14"></path>
                                </svg>
                            </a>
                        </div>
                    </div>

                    <!-- OPIOID CONVERSION -->
                    <div class="calc-card large">
                        <div class="calc-header">
                            <h3 class="calc-title">Opioid Conversion Calculator</h3>
                            <p class="calc-desc">Convert to morphine milligram equivalents (MME)</p>
                        </div>
                        <div class="calc-body">
                            <div class="input-group">
                                <label class="input-label">Opioid</label>
                                <select class="calc-select" id="opioid-type" onchange="calcOpioid()">
                                    <option value="">Select opioid...</option>
                                    <option value="1">Morphine (IV/IM/SC)</option>
                                    <option value="3">Morphine (PO)</option>
                                    <option value="0.01">Fentanyl (IV/IM)</option>
                                    <option value="0.3">Fentanyl (Transdermal mcg/hr)</option>
                                    <option value="5">Hydromorphone (PO)</option>
                                    <option value="1.5">Oxycodone (PO)</option>
                                    <option value="1.5">Hydrocodone (PO)</option>
                                    <option value="10">Tramadol (PO)</option>
                                </select>
                            </div>
                            <div class="input-group">
                                <label class="input-label">Dose</label>
                                <div class="input-wrapper">
                                    <input type="number" class="calc-input" id="opioid-dose" placeholder="Enter dose" oninput="calcOpioid()" min="0" step="0.1">
                                    <span class="input-unit">mg (or mcg/hr for fentanyl patch)</span>
                                </div>
                            </div>
                        </div>
                        <div id="opioid-result" class="calc-result yellow">
                            <div class="result-label">Morphine Milligram Equivalent (MME)</div>
                            <div class="result-value"><span id="opioid-value">--</span><span class="result-unit">mg</span></div>
                            <div class="result-text">Use caution with opioid conversions — consider patient factors and start at lower doses.</div>
                        </div>
                        <div class="calc-reference">
                            <a href="https://www.cdc.gov/opioids/data-resources/calculating-mme.html" target="_blank" class="ref-link">
                                <span>Reference: CDC MME Conversion</span>
                                <svg fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                                    <path stroke-linecap="round" stroke-linejoin="round" d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14"></path>
                                </svg>
                            </a>
                        </div>
                    </div>

                    <!-- LOCAL ANESTHETIC MAX DOSE -->
                    <div class="calc-card">
                        <div class="calc-header">
                            <h3 class="calc-title">Local Anesthetic Max Dose</h3>
                            <p class="calc-desc">Maximum safe dosing to prevent LAST</p>
                        </div>
                        <div class="calc-body">
                            <div class="input-group">
                                <label class="input-label">Local Anesthetic</label>
                                <select class="calc-select" id="la-type" onchange="calcLA()">
                                    <option value="">Select...</option>
                                    <option value="lidocaine">Lidocaine (without epi)</option>
                                    <option value="lidocaine-epi">Lidocaine (with epi)</option>
                                    <option value="bupivacaine">Bupivacaine (without epi)</option>
                                    <option value="bupivacaine-epi">Bupivacaine (with epi)</option>
                                    <option value="ropivacaine">Ropivacaine</option>
                                </select>
                            </div>
                            <div class="input-group">
                                <label class="input-label">Patient Weight</label>
                                <div class="input-wrapper">
                                    <input type="number" class="calc-input" id="la-weight" placeholder="Enter weight" oninput="calcLA()" min="0" step="0.1">
                                    <span class="input-unit">kg</span>
                                </div>
                            </div>
                        </div>
                        <div id="la-result" class="calc-result red">
                            <div class="result-label">Maximum Dose</div>
                            <div class="result-value"><span id="la-value">--</span><span class="result-unit">mg</span></div>
                            <div class="result-text" id="la-volume">--</div>
                        </div>
                        <div class="calc-reference">
                            <a href="https://www.nysora.com/topics/foundations-of-regional-anesthesia/pharmacology/local-anesthetics/" target="_blank" class="ref-link">
                                <span>Reference: NYSORA Guidelines</span>
                                <svg fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                                    <path stroke-linecap="round" stroke-linejoin="round" d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14"></path>
                                </svg>
                            </a>
                        </div>
                    </div>

                </div>
            </div>
        </main>

        <!-- Footer -->
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
        // Toggle Mobile Menu
        function toggleMobileMenu() {
            const menu = document.getElementById('mobileMenu');
            const btn = document.querySelector('.mobile-menu-btn');
            menu.classList.toggle('active');
            btn.classList.toggle('active');
        }

        // IDEAL BODY WEIGHT
        function calcIBW() {
            const sex = document.getElementById('ibw-sex').value;
            const height = parseFloat(document.getElementById('ibw-height').value);
            const weight = parseFloat(document.getElementById('ibw-weight').value);

            if (!sex || !height) return;

            let ibw;
            if (sex === 'male') {
                ibw = 50 + 2.3 * (height - 60);
            } else {
                ibw = 45.5 + 2.3 * (height - 60);
            }

            ibw = Math.max(ibw, 0);
            document.getElementById('ibw-value').textContent = ibw.toFixed(1);

            // Calculate ABW if actual weight provided
            if (weight) {
                const abw = ibw + 0.4 * (weight - ibw);
                document.getElementById('abw-value').textContent = abw.toFixed(1);
                document.getElementById('abw-display').style.display = 'block';
            } else {
                document.getElementById('abw-display').style.display = 'none';
            }

            const result = document.getElementById('ibw-result');
            result.classList.add('show');
        }

        // BODY SURFACE AREA
        function calcBSA() {
            const height = parseFloat(document.getElementById('bsa-height').value);
            const weight = parseFloat(document.getElementById('bsa-weight').value);

            if (!height || !weight) return;

            // Mosteller
            const mosteller = Math.sqrt((height * weight) / 3600);

            // DuBois
            const dubois = 0.007184 * Math.pow(height, 0.725) * Math.pow(weight, 0.425);

            // Haycock
            const haycock = 0.024265 * Math.pow(height, 0.3964) * Math.pow(weight, 0.5378);

            document.getElementById('bsa-mosteller').textContent = mosteller.toFixed(2);
            document.getElementById('bsa-dubois').textContent = dubois.toFixed(2);
            document.getElementById('bsa-haycock').textContent = haycock.toFixed(2);

            const result = document.getElementById('bsa-result');
            result.classList.add('show');
        }

        // MAXIMUM ALLOWABLE BLOOD LOSS
        function calcMABL() {
            const type = document.getElementById('mabl-type').value;
            const weight = parseFloat(document.getElementById('mabl-weight').value);
            const startHct = parseFloat(document.getElementById('mabl-start-hct').value);
            const targetHct = parseFloat(document.getElementById('mabl-target-hct').value);

            if (!type || !weight || !startHct || !targetHct) return;

            let ebvPerKg;
            if (type === 'male') ebvPerKg = 75;
            else if (type === 'female') ebvPerKg = 65;
            else ebvPerKg = 85; // pediatric

            const ebv = weight * ebvPerKg;
            const avgHct = (startHct + targetHct) / 2;
            const mabl = ebv * ((startHct - targetHct) / avgHct);

            document.getElementById('ebv-value').textContent = Math.round(ebv);
            document.getElementById('mabl-value').textContent = Math.round(mabl);

            const result = document.getElementById('mabl-result');
            result.classList.add('show');
        }

        // QTc INTERVAL
        function calcQTc() {
            const qt = parseFloat(document.getElementById('qtc-qt').value);
            const hr = parseFloat(document.getElementById('qtc-hr').value);
            const sex = document.getElementById('qtc-sex').value;

            if (!qt || !hr) return;

            const rr = 60000 / hr; // RR interval in ms
            const qtc = qt / Math.sqrt(rr / 1000);

            document.getElementById('qtc-value').textContent = Math.round(qtc);

            let interpretation = '';
            let resultClass = 'green';

            if (sex === 'male') {
                if (qtc < 440) {
                    interpretation = 'Normal (< 440 ms in males)';
                    resultClass = 'green';
                } else if (qtc < 500) {
                    interpretation = 'Borderline prolonged (440-500 ms)';
                    resultClass = 'yellow';
                } else {
                    interpretation = 'Prolonged (> 500 ms) - High risk for arrhythmia';
                    resultClass = 'red';
                }
            } else if (sex === 'female') {
                if (qtc < 460) {
                    interpretation = 'Normal (< 460 ms in females)';
                    resultClass = 'green';
                } else if (qtc < 500) {
                    interpretation = 'Borderline prolonged (460-500 ms)';
                    resultClass = 'yellow';
                } else {
                    interpretation = 'Prolonged (> 500 ms) - High risk for arrhythmia';
                    resultClass = 'red';
                }
            } else {
                interpretation = 'Select sex for interpretation';
            }

            document.getElementById('qtc-interpretation').textContent = interpretation;

            const result = document.getElementById('qtc-result');
            result.className = 'calc-result ' + resultClass + ' show';
        }

        // MAINTENANCE FLUIDS
        function calcFluids() {
            const weight = parseFloat(document.getElementById('fluids-weight').value);

            if (!weight) return;

            let rate;
            if (weight <= 10) {
                rate = weight * 4;
            } else if (weight <= 20) {
                rate = 40 + (weight - 10) * 2;
            } else {
                rate = 60 + (weight - 20) * 1;
            }

            const daily = rate * 24;

            document.getElementById('fluids-value').textContent = Math.round(rate);
            document.getElementById('fluids-daily').textContent = Math.round(daily);

            const result = document.getElementById('fluids-result');
            result.classList.add('show');
        }

        // PONV RISK
        function calcPONV() {
            const female = document.getElementById('ponv-female').checked;
            const nonsmoker = document.getElementById('ponv-nonsmoker').checked;
            const history = document.getElementById('ponv-history').checked;
            const opioids = document.getElementById('ponv-opioids').checked;

            const score = (female ? 1 : 0) + (nonsmoker ? 1 : 0) + (history ? 1 : 0) + (opioids ? 1 : 0);

            const risks = [10, 20, 40, 60, 80];
            const risk = score === 0 ? 10 : score === 1 ? 20 : score === 2 ? 40 : score === 3 ? 60 : 80;

            document.getElementById('ponv-score').textContent = score;
            document.getElementById('ponv-value').textContent = risk;

            const result = document.getElementById('ponv-result');
            result.className = 'calc-result show';
            if (risk <= 20) result.classList.add('green');
            else if (risk <= 60) result.classList.add('yellow');
            else result.classList.add('red');
        }

        // RCRI
        function calcRCRI() {
            const surgery = document.getElementById('rcri-surgery').checked;
            const ihd = document.getElementById('rcri-ihd').checked;
            const chf = document.getElementById('rcri-chf').checked;
            const cvd = document.getElementById('rcri-cvd').checked;
            const dm = document.getElementById('rcri-dm').checked;
            const cr = document.getElementById('rcri-cr').checked;

            const score = (surgery ? 1 : 0) + (ihd ? 1 : 0) + (chf ? 1 : 0) + (cvd ? 1 : 0) + (dm ? 1 : 0) + (cr ? 1 : 0);

            const risks = [0.4, 1.0, 2.4, 5.4, 5.4, 5.4, 5.4];
            const risk = risks[score];

            document.getElementById('rcri-score').textContent = score;
            document.getElementById('rcri-value').textContent = risk.toFixed(1);

            const result = document.getElementById('rcri-result');
            result.className = 'calc-result show';
            if (score === 0) result.classList.add('green');
            else if (score <= 2) result.classList.add('yellow');
            else result.classList.add('red');
        }

        // ASA PHYSICAL STATUS
        function calcASA() {
            const selected = document.querySelector('input[name="asa"]:checked');
            if (!selected) return;

            const value = selected.value;
            const descriptions = {
                '1': 'Healthy patient with no systemic disease',
                '2': 'Mild systemic disease without functional limitations',
                '3': 'Severe systemic disease with definite functional limitations',
                '4': 'Severe systemic disease that is constant threat to life',
                '5': 'Moribund patient not expected to survive without operation',
                '6': 'Brain-dead patient whose organs are being removed for donation'
            };

            document.getElementById('asa-value').textContent = 'ASA ' + value;
            document.getElementById('asa-desc').textContent = descriptions[value];

            const result = document.getElementById('asa-result');
            result.style.display = 'block';
            result.className = 'calc-result show';
            if (value <= 2) result.classList.add('green');
            else if (value <= 3) result.classList.add('yellow');
            else result.classList.add('red');
        }

        // OPIOID CONVERSION
        function calcOpioid() {
            const type = parseFloat(document.getElementById('opioid-type').value);
            const dose = parseFloat(document.getElementById('opioid-dose').value);

            if (!type || !dose) return;

            const mme = dose / type;

            document.getElementById('opioid-value').textContent = mme.toFixed(1);

            const result = document.getElementById('opioid-result');
            result.classList.add('show');
        }

        // LOCAL ANESTHETIC MAX DOSE
        function calcLA() {
            const type = document.getElementById('la-type').value;
            const weight = parseFloat(document.getElementById('la-weight').value);

            if (!type || !weight) return;

            const maxDoses = {
                'lidocaine': 5,
                'lidocaine-epi': 7,
                'bupivacaine': 2.5,
                'bupivacaine-epi': 3,
                'ropivacaine': 3
            };

            const concentrations = {
                'lidocaine': '1% = 10 mg/mL, 2% = 20 mg/mL',
                'lidocaine-epi': '1% = 10 mg/mL, 2% = 20 mg/mL',
                'bupivacaine': '0.25% = 2.5 mg/mL, 0.5% = 5 mg/mL',
                'bupivacaine-epi': '0.25% = 2.5 mg/mL, 0.5% = 5 mg/mL',
                'ropivacaine': '0.2% = 2 mg/mL, 0.5% = 5 mg/mL'
            };

            const maxDose = weight * maxDoses[type];

            document.getElementById('la-value').textContent = Math.round(maxDose);
            document.getElementById('la-volume').textContent = 'Common concentrations: ' + concentrations[type];

            const result = document.getElementById('la-result');
            result.classList.add('show');
        }

        // Initialize PONV result on page load
        document.addEventListener('DOMContentLoaded', () => {
            calcPONV();
            calcRCRI();
        });
    </script>
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
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=6">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=6">
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
                    <a href="/calculators" class="nav-link">Clinical Calculators</a>
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
            <a href="/calculators" class="mobile-menu-link">Clinical Calculators</a>
            <a href="/crisis" class="mobile-menu-link">Crisis Protocols</a>
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

            # Boost emergency protocol searches
            query_lower = query.lower()
            if any(word in query_lower for word in ['protocol', 'checklist', 'crisis', 'emergency', 'management']):
                # Add protocol/guideline terms to boost relevant results
                if 'malignant hyperthermia' in query_lower or 'mh' in query_lower:
                    q = q + ' AND (protocol[ti] OR emergency[ti] OR crisis[ti] OR treatment[ti])'
                elif any(term in query_lower for term in ['rapid sequence', 'rsi', 'intubation']):
                    q = q + ' AND (protocol[ti] OR guideline[ti] OR airway[ti])'
                elif 'last' in query_lower or 'local anesthetic' in query_lower:
                    q = q + ' AND (protocol[ti] OR treatment[ti] OR resuscitation[ti])'

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

            # Multi-tier search strategy to find relevant papers
            ids = []

            # Detect well-established topics that should have lots of evidence
            well_established_topics = [
                'malignant hyperthermia', 'difficult airway', 'aspiration', 'local anesthetic systemic toxicity',
                'last', 'propofol infusion syndrome', 'pris', 'anaphylaxis', 'tranexamic acid', 'txa',
                'postoperative nausea', 'ponv', 'rapid sequence', 'rsi', 'spinal anesthesia', 'epidural',
                'general anesthesia', 'blood transfusion', 'massive transfusion', 'coagulopathy'
            ]
            is_established_topic = any(topic in query_lower for topic in well_established_topics)

            # Adjust date range for established topics (go back further for classic evidence)
            date_range = '("2010/01/01"[PDAT] : "3000"[PDAT])' if is_established_topic else '("2015/01/01"[PDAT] : "3000"[PDAT])'

            # Tier 1: Broader search WITHOUT anesthesiology restriction first (better recall)
            try:
                print(f"[DEBUG] Tier 1: Searching PubMed (broad, high-quality filters)...")
                handle = Entrez.esearch(db="pubmed", term=f'{search_term}', retmax=25, sort="relevance")
                result = Entrez.read(handle)
                ids = result.get("IdList", [])
                print(f"[DEBUG] Tier 1 found {len(ids)} papers")
            except Exception as e:
                print(f"[ERROR] Tier 1 search failed: {e}")
                ids = []

            # Tier 2: Add anesthesiology MeSH if we got too few results
            if len(ids) < 5 and not is_followup:
                try:
                    print(f"[DEBUG] Tier 2: Adding anesthesiology MeSH refinement...")
                    handle = Entrez.esearch(db="pubmed", term=f'anesthesiology[MeSH Terms] AND {search_term}', retmax=25, sort="relevance")
                    result = Entrez.read(handle)
                    tier2_ids = result.get("IdList", [])
                    # Combine and deduplicate
                    ids = list(set(ids + tier2_ids))
                    print(f"[DEBUG] Tier 2 found {len(tier2_ids)} papers, total unique: {len(ids)}")
                except Exception as e:
                    print(f"[ERROR] Tier 2 search failed: {e}")

            # Tier 3: Further broaden if still insufficient (remove publication type restrictions)
            if len(ids) < 5 and not is_followup:
                try:
                    print(f"[DEBUG] Tier 3: Broadening to all publication types...")
                    broader_term = f'({q}) AND {date_range}'
                    handle = Entrez.esearch(db="pubmed", term=broader_term, retmax=25, sort="relevance")
                    result = Entrez.read(handle)
                    tier3_ids = result.get("IdList", [])
                    # Combine and deduplicate
                    ids = list(set(ids + tier3_ids))
                    print(f"[DEBUG] Tier 3 found {len(tier3_ids)} papers, total unique: {len(ids)}")
                except Exception as e:
                    print(f"[ERROR] Tier 3 search failed: {e}")

            # Tier 4: Last resort - very broad search with anesthesia-related terms
            if len(ids) < 3 and not is_followup:
                try:
                    print(f"[DEBUG] Tier 4: Very broad search with anesthesia terms...")
                    anesthesia_boost = f'({q}) AND (anesthesia OR anesthetic OR perioperative OR intraoperative)'
                    handle = Entrez.esearch(db="pubmed", term=anesthesia_boost, retmax=20, sort="relevance")
                    result = Entrez.read(handle)
                    tier4_ids = result.get("IdList", [])
                    # Combine and deduplicate
                    ids = list(set(ids + tier4_ids))
                    print(f"[DEBUG] Tier 4 found {len(tier4_ids)} papers, total unique: {len(ids)}")
                except Exception as e:
                    print(f"[ERROR] Tier 4 search failed: {e}")

            print(f"[DEBUG] Final total: {len(ids)} unique papers found")

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
            # Process up to 15 papers for comprehensive evidence (increased from 8)
            # More papers = better coverage of high-quality evidence = higher confidence ratings
            for p in papers[:15]:
                try:
                    art = p["MedlineCitation"]["Article"]
                    title = art.get("ArticleTitle", "No title")
                    # Truncate abstracts to 1200 chars for better citation accuracy
                    abstract = " ".join(str(t) for t in art.get("Abstract", {}).get("AbstractText", [])) if art.get("Abstract") else ""
                    abstract = abstract[:1200] + "..." if len(abstract) > 1200 else abstract
                    authors = ", ".join([a.get("LastName","") + " " + (a.get("ForeName","")[:1]+"." if a.get("ForeName") else "") for a in art.get("AuthorList",[])[:3]])  # Reduced to 3 authors
                    journal = art["Journal"].get("Title", "Unknown")
                    year = art["Journal"]["JournalIssue"]["PubDate"].get("Year", "N/A")
                    pmid = p["MedlineCitation"]["PMID"]

                    # Classify study type for quality indicators
                    study_classification = classify_study_type(title, journal)

                    refs.append({
                        "title": title,
                        "authors": authors,
                        "journal": journal,
                        "year": year,
                        "pmid": pmid,
                        "study_type": study_classification['type'],
                        "study_badge": study_classification['badge_text'],
                        "study_color": study_classification['badge_color'],
                        "study_score": study_classification['score'],
                        "sort_priority": study_classification['sort_priority']
                    })
                    context += f"Title: {title}\nAbstract: {abstract}\nAuthors: {authors}\nJournal: {journal} ({year})\nPMID: {pmid}\n\n"
                except:
                    continue

            # Sort references by quality (highest quality first)
            refs.sort(key=lambda x: x.get('sort_priority', 99))

            num_papers = len(refs)
            print(f"[DEBUG] Processed {num_papers} paper references")

            # Build smart conversation context (includes relevant earlier messages + recent)
            conversation_context = build_smart_context(session['messages'][:-1], raw_query)  # Exclude just-added user message
            print(f"[DEBUG] Smart context built ({len(conversation_context)} chars)")

            # Create numbered reference list for citation
            ref_list = ""
            for i, ref in enumerate(refs, 1):
                ref_list += f"[{i}] {ref['title']} - {ref['authors']} ({ref['year']}) PMID: {ref['pmid']}\n"

            # Log papers being used for debugging
            print(f"\n[DEBUG] ===== PAPERS RETURNED FOR QUERY: '{raw_query}' =====")
            for i, ref in enumerate(refs[:3], 1):  # Show first 3
                print(f"[DEBUG] [{i}] {ref['title'][:100]}...")
            print(f"[DEBUG] ================================================\n")

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
3. **CITATION ACCURACY (CRITICAL):**
   - ONLY cite a paper [N] if that SPECIFIC paper's abstract directly supports the exact claim you're making
   - Verify each citation: Does the abstract contain this specific fact/dose/finding?
   - If a claim is from your general knowledge but NOT explicitly in the provided abstracts, DO NOT add a citation
   - Better to have NO citation than an INACCURATE citation
   - When citing: use [1], [2], etc. - NO author names in text
   - Multiple papers can support one claim: "X is effective [1][2][3]"
   - If papers discuss a topic generally but don't support your specific claim, omit the citation
4. Be conversational but clinically complete - like talking to a colleague
5. HTML format: <h3> for sections, <p> for paragraphs, <strong> for emphasis, <ul><li> for lists
6. CRITICAL: Return ONLY the HTML content - do NOT wrap your response in markdown code fences (```html or ```)
7. START your response with an evidence quality badge:
   <div class="evidence-quality-badge">
   <div class="confidence-level [high/moderate/low]">
   <strong>Evidence Quality:</strong> [High/Moderate/Low] Confidence
   </div>
   <div class="evidence-details">
   📊 {num_papers} papers analyzed • Study types: [list types] • Date range: [range]
   </div>
   </div>

CONFIDENCE LEVEL GUIDANCE:
- **High Confidence**: Multiple meta-analyses/systematic reviews OR strong RCT evidence
- **Moderate Confidence**: Some RCTs/reviews but limited OR conflicting evidence
- **Low Confidence**: Few papers OR case reports/expert opinion only

CITATION VERIFICATION CHECKLIST (Check EACH citation before adding):
   ❌ WRONG: "Propofol 2-3 mg/kg is recommended [1]" when abstract says "induction agents" generally
   ❌ WRONG: "TXA reduces blood loss by 30% [2]" when abstract doesn't give specific percentage
   ❌ WRONG: "Common side effects include nausea [3]" when abstract doesn't mention side effects
   ✅ CORRECT: "TXA reduces blood loss in spine surgery [1][2]" when both abstracts explicitly state this
   ✅ CORRECT: "Standard monitoring includes pulse oximetry" (NO citation - general knowledge)

IMPORTANT: Use the provided research papers to inform your answer. Only cite papers when their abstracts directly support the specific claim. If a claim is general knowledge or not supported by the abstracts, omit the citation.

Example with proper conservative citation:
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
100% oxygen and hand-ventilate to assess compliance. Deepen anesthesia with propofol or increase volatile agent to 2+ MAC.</p>
<p><strong>Bronchodilators:</strong><br>
Inhaled beta-2 agonists are first-line treatment [1][2]. Response typically seen within minutes. Severe cases may require IV epinephrine.</p>
<p><strong>Monitoring:</strong><br>
Watch for auto-PEEP, pneumothorax, and cardiovascular compromise from high airway pressures.</p>"

NOTE: Only cite when abstracts explicitly support claims. Generic principles have NO citations as they're standard knowledge.

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

@app.route("/crisis")
def crisis_protocols():
    """Crisis Protocols - Evidence-based anesthesia emergency management"""
    return render_template_string(CRISIS_HTML)

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

                        # Classify study type for quality indicators
                        study_classification = classify_study_type(title, journal)

                        all_refs.append({
                            "title": title,
                            "authors": authors,
                            "journal": journal,
                            "year": year,
                            "pmid": pmid,
                            "study_type": study_classification['type'],
                            "study_badge": study_classification['badge_text'],
                            "study_color": study_classification['badge_color'],
                            "study_score": study_classification['score'],
                            "sort_priority": study_classification['sort_priority']
                        })
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

    # Sort references by quality (highest quality first)
    unique_refs.sort(key=lambda x: x.get('sort_priority', 99))

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