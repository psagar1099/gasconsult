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
import uuid
import json
import bleach
import redis
import hashlib
import datetime
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

app = Flask(__name__)
app.secret_key = os.getenv('FLASK_SECRET_KEY', secrets.token_hex(32))

# Configure server-side sessions with Redis backend (production-grade, multi-worker safe)
# Falls back to localhost Redis if REDIS_URL not set (for local development)
redis_url = os.getenv('REDIS_URL', 'redis://localhost:6379')
try:
    # Initialize Redis connection
    redis_client = redis.from_url(
        redis_url,
        decode_responses=False,  # Keep binary for Flask-Session compatibility
        socket_connect_timeout=5,
        socket_keepalive=True,
        health_check_interval=30
    )
    # Test connection
    redis_client.ping()
    print(f"[REDIS] Successfully connected to Redis at {redis_url.split('@')[-1]}")  # Hide credentials in logs

    app.config['SESSION_TYPE'] = 'redis'
    app.config['SESSION_REDIS'] = redis_client
except (redis.ConnectionError, redis.TimeoutError) as e:
    print(f"[WARNING] Redis connection failed: {e}")
    print(f"[WARNING] Falling back to filesystem sessions (not recommended for production)")
    # Fallback to filesystem sessions for development
    app.config['SESSION_TYPE'] = 'filesystem'
    app.config['SESSION_FILE_DIR'] = tempfile.gettempdir()

app.config['SESSION_PERMANENT'] = True
app.config['PERMANENT_SESSION_LIFETIME'] = 3600  # 1 hour
app.config['SESSION_USE_SIGNER'] = True
app.config['SESSION_COOKIE_SAMESITE'] = 'Lax'
app.config['SESSION_COOKIE_HTTPONLY'] = True
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
    # MH - Expand abbreviation to malignant hyperthermia AND exclude psychiatric papers
    # Use word boundary to avoid matching "smith" or other words containing "mh"
    # SIMPLIFIED: Use MeSH term only with exclusions - don't add synonyms that complicate the query
    q = re.sub(r'\bmh\b', '("malignant hyperthermia"[MeSH Terms] NOT (mental[ti] OR psychiatric[ti] OR psychology[ti] OR depression[ti]))', q, flags=re.IGNORECASE)
    # Expand full term "malignant hyperthermia" to use MeSH for better indexing
    # Use negative lookbehind/lookahead to avoid double-quoting already expanded terms
    q = re.sub(r'(?<!")malignant hyperthermia(?!")', '"malignant hyperthermia"[MeSH Terms]', q, flags=re.IGNORECASE)

    # Common scenarios
    q = q.replace("full stomach", '"aspiration"[MeSH Terms] OR "rapid sequence" OR RSI OR "aspiration risk"')
    q = q.replace("difficult airway", '"difficult airway"[MeSH Terms] OR "airway management" OR intubation OR "difficult intubation"')
    q = q.replace("anticipated difficult", '"difficult airway"[MeSH Terms] OR "airway management" OR "awake intubation"')

    # Regional anesthesia
    q = q.replace("nerve block", '"nerve block"[MeSH Terms] OR "regional anesthesia" OR "peripheral nerve block"')
    q = q.replace("epidural", '"epidural"[MeSH Terms] OR "epidural anesthesia" OR "neuraxial"')
    q = q.replace("spinal", '"spinal anesthesia"[MeSH Terms] OR "subarachnoid" OR "neuraxial"')

    # Other common anesthesia abbreviations that might be ambiguous
    q = re.sub(r'\bosa\b', '("obstructive sleep apnea" OR OSA OR "sleep apnea")', q, flags=re.IGNORECASE)
    q = re.sub(r'\bcopd\b', '("chronic obstructive pulmonary disease" OR COPD OR "obstructive lung disease")', q, flags=re.IGNORECASE)
    q = re.sub(r'\bchf\b', '("congestive heart failure" OR CHF OR "heart failure")', q, flags=re.IGNORECASE)
    q = re.sub(r'\bcad\b', '("coronary artery disease" OR CAD OR "ischemic heart disease")', q, flags=re.IGNORECASE)
    q = re.sub(r'\bckd\b', '("chronic kidney disease" OR CKD OR "renal insufficiency")', q, flags=re.IGNORECASE)
    # PE - exclude physical exam papers
    q = re.sub(r'\bpe\b', '(("pulmonary embolism"[MeSH Terms] OR "pulmonary embolus" OR thromboembolic) NOT "physical exam"[ti])', q, flags=re.IGNORECASE)
    # GA - general anesthesia (exclude genetic analysis)
    q = re.sub(r'\bga\b', '(("general anesthesia"[MeSH Terms] OR "general anaesthesia") NOT genetic[ti])', q, flags=re.IGNORECASE)
    # MAC - anesthesia context (both meanings are anesthesia-related, so just expand)
    q = re.sub(r'\bmac\b', '("monitored anesthesia care" OR "minimum alveolar concentration" OR MAC)', q, flags=re.IGNORECASE)

    # NOTE: DO NOT expand common words like "protocol", "crisis", "management"
    # These generic expansions create overly broad queries that return wrong papers
    # Let them remain as-is so PubMed can properly combine them with the medical terms

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
    <!-- Google Analytics 4 -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=G-01NZYD1DPP"></script>
    <script>
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-01NZYD1DPP');
    </script>

    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pre-Operative Assessment - gasconsult.ai</title>

    <!-- SEO Meta Tags -->
    <meta name="description" content="AI-powered pre-operative assessment tool combining evidence-based protocols with patient-specific factors for comprehensive anesthesia planning.">

    <!-- Open Graph / Facebook -->
    <meta property="og:type" content="website">
    <meta property="og:url" content="https://gasconsult.ai/preop">
    <meta property="og:title" content="Pre-Operative Assessment - gasconsult.ai">
    <meta property="og:description" content="AI-powered pre-operative assessment tool for comprehensive anesthesia planning.">
    <meta property="og:image" content="https://gasconsult.ai/static/logo.png">

    <!-- Twitter -->
    <meta property="twitter:card" content="summary_large_image">
    <meta property="twitter:url" content="https://gasconsult.ai/preop">
    <meta property="twitter:title" content="Pre-Operative Assessment - gasconsult.ai">
    <meta property="twitter:description" content="AI-powered pre-operative assessment tool for comprehensive anesthesia planning.">
    <meta property="twitter:image" content="https://gasconsult.ai/static/logo.png">

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

        /* Background Canvas */
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

        .nav-dropdown:has(.nav-dropdown-link.active) .nav-dropdown-toggle {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .nav-dropdown {
            position: relative;
            display: inline-block;
        }

        .nav-dropdown-toggle {
            cursor: pointer;
            background: none;
            border: none;
            font-family: inherit;
        }

        .nav-dropdown-menu {
            display: none;
            position: absolute;
            top: 100%;
            right: 0;
            background: white;
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            min-width: 200px;
            margin-top: 4px;
            z-index: 1000;
            overflow: hidden;
        }

        .nav-dropdown-menu.show {
            display: block;
        }

        .nav-dropdown-link {
            display: block;
            padding: 12px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .nav-dropdown-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
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
        }

        .mobile-menu-btn span {
            display: block;
            width: 22px;
            height: 2px;
            background: var(--gray-700);
            border-radius: 1px;
            transition: all 0.3s ease;
        }

        .mobile-menu {
            display: none;
            position: fixed;
            top: 80px;
            left: 16px;
            right: 16px;
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 8px;
            box-shadow: 0 4px 16px rgba(0,0,0,0.08), 0 12px 48px rgba(0,0,0,0.12);
            z-index: 99;
            flex-direction: column;
            gap: 4px;
        }

        .mobile-menu.active { display: flex; }

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

        /* Main Content */
        .main-content {
            flex: 1;
            padding: 88px 16px 32px;
            max-width: 800px;
            margin: 0 auto;
            width: 100%;
        }

        /* Header */
        .header {
            text-align: center;
            margin-bottom: 32px;
            animation: fade-up 0.6s ease forwards;
        }

        @keyframes fade-up {
            from { opacity: 0; transform: translateY(20px); }
            to { opacity: 1; transform: translateY(0); }
        }

        .header-title {
            font-size: 28px;
            font-weight: 800;
            letter-spacing: -1px;
            color: var(--gray-900);
            margin-bottom: 8px;
        }

        .header-title .blue { color: var(--blue-600); }

        .header-subtitle {
            font-size: 15px;
            color: var(--gray-500);
            line-height: 1.6;
        }

        /* Form Card */
        .form-card {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(40px) saturate(180%);
            -webkit-backdrop-filter: blur(40px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 24px;
            padding: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06);
            animation: fade-up 0.6s 0.1s ease forwards;
            opacity: 0;
        }

        /* Form Groups */
        .form-group {
            margin-bottom: 24px;
        }

        .form-label {
            display: block;
            font-size: 14px;
            font-weight: 600;
            color: var(--gray-700);
            margin-bottom: 8px;
        }

        .form-label .required {
            color: #EF4444;
            margin-left: 2px;
        }

        /* Input Styles */
        .form-input {
            width: 100%;
            padding: 12px 16px;
            font-size: 15px;
            font-family: inherit;
            color: var(--gray-800);
            background: var(--white);
            border: 2px solid var(--gray-200);
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .form-input:focus {
            outline: none;
            border-color: var(--blue-500);
            box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.1);
        }

        .form-input::placeholder {
            color: var(--gray-400);
        }

        /* Input Row for Multiple Fields */
        .input-row {
            display: grid;
            gap: 12px;
            grid-template-columns: 1fr;
        }

        /* Radio Buttons */
        .radio-group {
            display: flex;
            flex-direction: column;
            gap: 10px;
        }

        .radio-option {
            position: relative;
            display: flex;
            align-items: center;
            padding: 14px 16px;
            background: var(--white);
            border: 2px solid var(--gray-200);
            border-radius: 12px;
            cursor: pointer;
            transition: all 0.2s ease;
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

        .radio-option input[type="radio"]:checked + .radio-visual {
            background: var(--blue-600);
            border-color: var(--blue-600);
        }

        .radio-option input[type="radio"]:checked + .radio-visual::after {
            opacity: 1;
        }

        .radio-option input[type="radio"]:checked ~ .radio-label {
            font-weight: 600;
            color: var(--blue-600);
        }

        .radio-option.selected {
            background: var(--blue-50);
            border-color: var(--blue-500);
            box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.12);
        }

        .radio-visual {
            width: 20px;
            height: 20px;
            border: 2px solid var(--gray-300);
            border-radius: 50%;
            flex-shrink: 0;
            position: relative;
            transition: all 0.2s ease;
            margin-right: 12px;
        }

        .radio-visual::after {
            content: '';
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            width: 8px;
            height: 8px;
            border-radius: 50%;
            background: white;
            opacity: 0;
            transition: opacity 0.2s ease;
        }

        .radio-label {
            font-size: 14px;
            color: var(--gray-700);
            flex: 1;
            transition: all 0.2s ease;
        }

        /* Checkboxes */
        .checkbox-group {
            display: grid;
            grid-template-columns: 1fr;
            gap: 10px;
        }

        .checkbox-option {
            display: flex;
            align-items: center;
            padding: 12px 14px;
            background: var(--white);
            border: 2px solid var(--gray-200);
            border-radius: 10px;
            cursor: pointer;
            transition: all 0.2s ease;
        }

        .checkbox-option:hover {
            border-color: var(--blue-300);
            box-shadow: 0 2px 6px rgba(37, 99, 235, 0.06);
        }

        .checkbox-option input[type="checkbox"] {
            width: 18px;
            height: 18px;
            margin-right: 10px;
            cursor: pointer;
            accent-color: var(--blue-600);
        }

        .checkbox-label {
            font-size: 14px;
            color: var(--gray-700);
            cursor: pointer;
            flex: 1;
        }

        /* Textarea */
        textarea.form-input {
            resize: vertical;
            min-height: 80px;
        }

        /* Section Divider */
        .section-divider {
            height: 1px;
            background: linear-gradient(90deg, transparent 0%, var(--gray-200) 50%, transparent 100%);
            margin: 32px 0;
        }

        .section-title {
            font-size: 16px;
            font-weight: 700;
            color: var(--gray-800);
            margin-bottom: 16px;
        }

        /* Submit Button */
        .submit-btn {
            width: 100%;
            padding: 16px;
            font-size: 16px;
            font-weight: 600;
            color: white;
            background: linear-gradient(135deg, var(--blue-600) 0%, var(--blue-700) 100%);
            border: none;
            border-radius: 14px;
            cursor: pointer;
            transition: all 0.3s ease;
            box-shadow: 0 2px 8px rgba(37,99,235,0.3), 0 8px 24px rgba(37,99,235,0.2);
        }

        .submit-btn:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(37,99,235,0.4), 0 12px 32px rgba(37,99,235,0.25);
        }

        .submit-btn:active {
            transform: translateY(0);
        }

        /* Results State */
        .results-container {
            animation: fade-up 0.6s ease forwards;
        }

        .results-header {
            text-align: center;
            padding: 40px 32px;
            background: linear-gradient(135deg, var(--blue-50) 0%, var(--blue-100) 100%);
            border-radius: 24px;
            margin-bottom: 32px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 12px 48px rgba(37,99,235,0.08);
        }

        .results-icon {
            width: 56px;
            height: 56px;
            margin: 0 auto 16px;
            background: var(--blue-600);
            border-radius: 50%;
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .results-icon svg {
            width: 28px;
            height: 28px;
            stroke: white;
        }

        .results-title {
            font-size: 24px;
            font-weight: 800;
            color: var(--gray-900);
            margin-bottom: 8px;
            letter-spacing: -0.5px;
        }

        .results-subtitle {
            font-size: 14px;
            color: var(--gray-600);
        }

        .results-card {
            background: rgba(255,255,255,0.85);
            backdrop-filter: blur(40px) saturate(180%);
            -webkit-backdrop-filter: blur(40px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 24px;
            padding: 32px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06);
        }

        .results-content {
            font-size: 15px;
            line-height: 1.8;
            color: var(--gray-700);
        }

        /* Assessment Content Typography */
        .results-content h3 {
            font-size: 20px;
            font-weight: 800;
            color: var(--gray-900);
            margin-top: 32px;
            margin-bottom: 16px;
            padding-bottom: 12px;
            border-bottom: 2px solid var(--blue-100);
            letter-spacing: -0.5px;
            display: flex;
            align-items: center;
            gap: 10px;
        }

        .results-content h3:first-child {
            margin-top: 0;
        }

        .results-content h3::before {
            content: '';
            width: 4px;
            height: 24px;
            background: linear-gradient(135deg, var(--blue-500) 0%, var(--blue-600) 100%);
            border-radius: 2px;
            flex-shrink: 0;
        }

        .results-content h4 {
            font-size: 17px;
            font-weight: 700;
            color: var(--gray-800);
            margin-top: 24px;
            margin-bottom: 12px;
            letter-spacing: -0.3px;
        }

        .results-content p {
            margin-bottom: 16px;
            line-height: 1.8;
        }

        .results-content p:last-child {
            margin-bottom: 0;
        }

        .results-content strong {
            font-weight: 700;
            color: var(--gray-900);
        }

        .results-content ul,
        .results-content ol {
            margin: 16px 0;
            padding-left: 0;
            list-style: none;
        }

        .results-content ul li,
        .results-content ol li {
            position: relative;
            padding-left: 32px;
            margin-bottom: 12px;
            line-height: 1.8;
        }

        .results-content ul li::before {
            content: '';
            position: absolute;
            left: 8px;
            top: 11px;
            width: 6px;
            height: 6px;
            background: var(--blue-500);
            border-radius: 50%;
        }

        .results-content ol {
            counter-reset: item;
        }

        .results-content ol li {
            counter-increment: item;
        }

        .results-content ol li::before {
            content: counter(item);
            position: absolute;
            left: 0;
            top: 0;
            width: 24px;
            height: 24px;
            background: var(--blue-50);
            color: var(--blue-600);
            border-radius: 6px;
            font-size: 12px;
            font-weight: 700;
            display: flex;
            align-items: center;
            justify-content: center;
        }

        /* Info boxes within assessment */
        .results-content blockquote {
            margin: 20px 0;
            padding: 16px 20px;
            background: linear-gradient(135deg, var(--blue-50) 0%, rgba(239, 246, 255, 0.5) 100%);
            border-left: 4px solid var(--blue-500);
            border-radius: 12px;
            font-style: normal;
        }

        /* Citation styling */
        .results-content a {
            color: var(--blue-600);
            text-decoration: none;
            font-weight: 600;
            transition: all 0.2s ease;
        }

        .results-content a:hover {
            color: var(--blue-700);
            text-decoration: underline;
        }

        /* Code/monospace elements */
        .results-content code {
            background: var(--gray-100);
            padding: 2px 6px;
            border-radius: 4px;
            font-family: 'Monaco', 'Courier New', monospace;
            font-size: 14px;
            color: var(--gray-800);
        }

        /* References */
        .references-card {
            background: rgba(255,255,255,0.85);
            backdrop-filter: blur(40px) saturate(180%);
            -webkit-backdrop-filter: blur(40px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 24px;
            padding: 32px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06);
        }

        .references-title {
            display: flex;
            align-items: center;
            gap: 10px;
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 20px;
        }

        .references-title svg {
            width: 20px;
            height: 20px;
            stroke: var(--blue-600);
        }

        .reference-item {
            padding: 16px 0;
            border-bottom: 1px solid var(--gray-200);
        }

        .reference-item:last-child {
            border-bottom: none;
        }

        .reference-number {
            display: inline-block;
            width: 32px;
            height: 32px;
            background: var(--blue-50);
            color: var(--blue-600);
            border-radius: 8px;
            text-align: center;
            line-height: 32px;
            font-size: 13px;
            font-weight: 700;
            margin-right: 12px;
            flex-shrink: 0;
        }

        .reference-link {
            color: var(--blue-600);
            text-decoration: none;
            font-weight: 600;
            font-size: 14px;
            line-height: 1.5;
        }

        .reference-link:hover {
            text-decoration: underline;
        }

        .reference-meta {
            font-size: 13px;
            color: var(--gray-600);
            margin-top: 6px;
        }

        /* Action Buttons */
        .action-buttons {
            display: flex;
            flex-direction: column;
            gap: 12px;
        }

        .btn {
            padding: 14px 20px;
            font-size: 15px;
            font-weight: 600;
            text-align: center;
            text-decoration: none;
            border-radius: 12px;
            cursor: pointer;
            transition: all 0.2s ease;
            display: inline-flex;
            align-items: center;
            justify-content: center;
            gap: 8px;
        }

        .btn-primary {
            background: linear-gradient(135deg, var(--blue-600) 0%, var(--blue-700) 100%);
            color: white;
            border: none;
            box-shadow: 0 2px 8px rgba(37,99,235,0.3);
        }

        .btn-primary:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(37,99,235,0.4);
        }

        .btn-secondary {
            background: white;
            color: var(--gray-700);
            border: 2px solid var(--gray-300);
        }

        .btn-secondary:hover {
            background: var(--gray-50);
            border-color: var(--gray-400);
        }

        /* Footer */
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

        .footer-link:hover {
            color: var(--gray-700);
        }

        /* Tablet & Desktop */
        @media (min-width: 640px) {
            .input-row {
                grid-template-columns: repeat(2, 1fr);
            }

            .checkbox-group {
                grid-template-columns: repeat(2, 1fr);
            }

            .radio-group {
                flex-direction: row;
                flex-wrap: wrap;
            }

            .radio-option {
                flex: 1;
                min-width: calc(50% - 5px);
            }
        }

        @media (min-width: 768px) {
            .nav { padding: 16px 32px; }
            .nav-inner { height: 64px; padding: 0 24px; border-radius: 20px; }
            .logo-icon svg { width: 42px; height: 15px; }
            .logo-text { font-size: 20px; }
            .nav-links { display: flex; }
            .mobile-menu-btn { display: none; }

            .main-content {
                padding: 108px 32px 48px;
            }

            .header-title {
                font-size: 36px;
            }

            .form-card {
                padding: 32px;
            }

            .action-buttons {
                flex-direction: row;
            }

            .footer { padding: 40px 32px; }
            .footer-inner { flex-direction: row; justify-content: space-between; text-align: left; }
            .footer-text { font-size: 14px; }
            .footer-links { gap: 32px; }
            .footer-link { font-size: 14px; }
        }

        @media (min-width: 1024px) {
            .nav { padding: 16px 40px; }

            .main-content {
                max-width: 900px;
            }

            .header-title {
                font-size: 42px;
            }

            .form-card,
            .results-card,
            .references-card {
                padding: 40px;
            }

            .results-header {
                padding: 48px 40px;
            }
        }

        /* Mobile-specific adjustments */
        @media (max-width: 639px) {
            .results-card,
            .references-card {
                padding: 20px;
                border-radius: 20px;
            }

            .results-header {
                padding: 24px 20px;
                border-radius: 20px;
            }

            .results-content h3 {
                font-size: 18px;
            }

            .results-content h4 {
                font-size: 16px;
            }

            .results-title {
                font-size: 20px;
            }
        }
    </style>
</head>
<body>
    <a href="#main-content" class="skip-to-content">Skip to main content</a>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <!-- Navigation -->
        <nav class="nav" role="navigation" aria-label="Main navigation">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo" aria-label="GasConsult.ai home">
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
                    <div class="nav-dropdown">
                        <button class="nav-link nav-dropdown-toggle" onclick="toggleNavDropdown(event)">More ▼</button>
                        <div class="nav-dropdown-menu">
                            <a href="/hypotension" class="nav-dropdown-link">IOH Predictor</a>
                            <a href="/difficult-airway" class="nav-dropdown-link">Difficult Airway</a>
                            <a href="/informed-consent" class="nav-dropdown-link">Informed Consent</a>
                        </div>
                    </div>
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
            <a href="/difficult-airway" class="mobile-menu-link">Difficult Airway</a>
            <a href="/informed-consent" class="mobile-menu-link">Informed Consent</a>
        </div>

        <!-- Main Content -->
        <main class="main-content">
            {% if not summary %}
            <!-- Form State -->
            <div class="header">
                <h1 class="header-title"><span class="blue">Pre-Operative</span> Assessment</h1>
                <p class="header-subtitle">Evidence-based risk stratification with personalized recommendations</p>
            </div>

            <div class="form-card">
                <form method="POST" action="/preop">
                    <input type="hidden" name="csrf_token" value="{{ csrf_token() }}"/>

                    <!-- Patient Demographics -->
                    <div class="section-title">Patient Information</div>
                    
                    <div class="input-row">
                        <div class="form-group">
                            <label class="form-label">Age<span class="required">*</span></label>
                            <input type="number" name="age" class="form-input" placeholder="65" required min="0" max="120">
                        </div>
                        <div class="form-group">
                            <label class="form-label">Weight (kg)<span class="required">*</span></label>
                            <input type="number" name="weight" class="form-input" placeholder="70" required min="1" max="300" step="0.1">
                        </div>
                    </div>

                    <div class="input-row">
                        <div class="form-group">
                            <label class="form-label">Height (cm)<span class="required">*</span></label>
                            <input type="number" name="height" class="form-input" placeholder="170" required min="50" max="250">
                        </div>
                        <div class="form-group">
                            <label class="form-label">Sex<span class="required">*</span></label>
                            <select name="sex" class="form-input" required>
                                <option value="">Select...</option>
                                <option value="male">Male</option>
                                <option value="female">Female</option>
                            </select>
                        </div>
                    </div>

                    <div class="section-divider"></div>

                    <!-- Functional Capacity -->
                    <div class="section-title">Functional Capacity</div>
                    <div class="form-group">
                        <label class="form-label">METs Capacity<span class="required">*</span></label>
                        <div class="radio-group">
                            <label class="radio-option">
                                <input type="radio" name="mets" value="<4 METs" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">&lt;4 METs - Poor</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="mets" value="4-10 METs" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">4-10 METs - Moderate</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="mets" value=">10 METs" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">&gt;10 METs - Excellent</span>
                            </label>
                        </div>
                    </div>

                    <div class="section-divider"></div>

                    <!-- Comorbidities -->
                    <div class="section-title">Comorbidities</div>
                    <div class="form-group">
                        <label class="form-label">Select all that apply</label>
                        <div class="checkbox-group">
                            <label class="checkbox-option">
                                <input type="checkbox" name="comorbidities" value="Hypertension">
                                <span class="checkbox-label">Hypertension</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="comorbidities" value="Diabetes Mellitus">
                                <span class="checkbox-label">Diabetes</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="comorbidities" value="Coronary Artery Disease">
                                <span class="checkbox-label">CAD</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="comorbidities" value="Heart Failure">
                                <span class="checkbox-label">Heart Failure</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="comorbidities" value="Chronic Kidney Disease">
                                <span class="checkbox-label">CKD</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="comorbidities" value="Obstructive Sleep Apnea">
                                <span class="checkbox-label">OSA</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="comorbidities" value="COPD">
                                <span class="checkbox-label">COPD</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="comorbidities" value="Prior Stroke">
                                <span class="checkbox-label">Prior Stroke</span>
                            </label>
                        </div>
                    </div>

                    <div class="form-group">
                        <label class="form-label">Other Comorbidities</label>
                        <input type="text" name="other_comorbidities" class="form-input" placeholder="Enter any other conditions...">
                    </div>

                    <div class="section-divider"></div>

                    <!-- Medications -->
                    <div class="section-title">Medications & Allergies</div>
                    <div class="form-group">
                        <label class="form-label">Current Medications</label>
                        <textarea name="medications" class="form-input" placeholder="List medications (e.g., aspirin, metoprolol, lisinopril...)"></textarea>
                    </div>

                    <div class="form-group">
                        <label class="form-label">Allergies</label>
                        <input type="text" name="allergies" class="form-input" placeholder="Drug allergies and reactions...">
                    </div>

                    <div class="section-divider"></div>

                    <!-- Procedure Details -->
                    <div class="section-title">Surgical Procedure</div>
                    <div class="form-group">
                        <label class="form-label">Planned Procedure<span class="required">*</span></label>
                        <input type="text" name="procedure" class="form-input" placeholder="e.g., Total hip arthroplasty" required>
                    </div>

                    <div class="form-group">
                        <label class="form-label">Surgery Risk Category<span class="required">*</span></label>
                        <div class="radio-group">
                            <label class="radio-option">
                                <input type="radio" name="surgery_risk" value="low" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">Low Risk</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="surgery_risk" value="intermediate" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">Intermediate</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="surgery_risk" value="high" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">High Risk</span>
                            </label>
                        </div>
                    </div>

                    <div class="section-divider"></div>

                    <!-- Additional Details -->
                    <div class="section-title">Additional Information (Optional)</div>
                    
                    <div class="input-row">
                        <div class="form-group">
                            <label class="form-label">Hemoglobin (g/dL)</label>
                            <input type="text" name="hgb" class="form-input" placeholder="12.5">
                        </div>
                        <div class="form-group">
                            <label class="form-label">Creatinine (mg/dL)</label>
                            <input type="text" name="cr" class="form-input" placeholder="1.0">
                        </div>
                    </div>

                    <div class="input-row">
                        <div class="form-group">
                            <label class="form-label">Platelets (×10³/µL)</label>
                            <input type="text" name="plt" class="form-input" placeholder="250">
                        </div>
                        <div class="form-group">
                            <label class="form-label">INR</label>
                            <input type="text" name="inr" class="form-input" placeholder="1.0">
                        </div>
                    </div>

                    <div class="form-group">
                        <label class="form-label">Ejection Fraction (%)</label>
                        <input type="text" name="ef" class="form-input" placeholder="55-60%">
                    </div>

                    <div class="form-group">
                        <label class="form-label">EKG Findings</label>
                        <input type="text" name="ekg" class="form-input" placeholder="Normal sinus rhythm">
                    </div>

                    <div class="form-group">
                        <label class="form-label">NPO Status</label>
                        <input type="text" name="npo" class="form-input" placeholder="NPO after midnight">
                    </div>

                    <div class="form-group">
                        <label class="form-label">Previous Anesthesia Issues</label>
                        <textarea name="previous_anesthesia" class="form-input" placeholder="Any prior anesthesia complications, difficult airway, PONV history..."></textarea>
                    </div>

                    <button type="submit" class="submit-btn">Generate Assessment</button>
                </form>
            </div>

            {% else %}
            <!-- Results State -->
            <div class="results-container">
                <div class="results-header">
                    <div class="results-icon">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <path d="M22 11.08V12a10 10 0 1 1-5.93-9.14"></path>
                            <polyline points="22 4 12 14.01 9 11.01"></polyline>
                        </svg>
                    </div>
                    <h1 class="results-title"><span style="color: var(--blue-600);">Pre-Operative</span> Assessment Complete</h1>
                    <p class="results-subtitle">Evidence-based risk stratification and optimization recommendations</p>
                </div>

                <div class="results-card">
                    <div class="results-content">
                        {{ summary|safe }}
                    </div>
                </div>

                {% if references %}
                <div class="references-card">
                    <div class="references-title">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path>
                            <path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path>
                        </svg>
                        <span>Evidence-Based References</span>
                    </div>
                    {% for ref in references %}
                    <div class="reference-item">
                        <span class="reference-number">[{{ loop.index }}]</span>
                        <a href="https://pubmed.ncbi.nlm.nih.gov/{{ ref.pmid }}/" target="_blank" rel="noopener noreferrer" class="reference-link">
                            {{ ref.title }}
                        </a>
                        <div class="reference-meta">{{ ref.authors }} - {{ ref.journal }}, {{ ref.year }}</div>
                    </div>
                    {% endfor %}
                </div>
                {% endif %}

                <div class="action-buttons">
                    <a href="/preop" class="btn btn-primary">
                        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <line x1="12" y1="5" x2="12" y2="19"></line>
                            <line x1="5" y1="12" x2="19" y2="12"></line>
                        </svg>
                        New Assessment
                    </a>
                    <button onclick="window.print()" class="btn btn-secondary">
                        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <polyline points="6 9 6 2 18 2 18 9"></polyline>
                            <path d="M6 18H4a2 2 0 0 1-2-2v-5a2 2 0 0 1 2-2h16a2 2 0 0 1 2 2v5a2 2 0 0 1-2 2h-2"></path>
                            <rect x="6" y="14" width="12" height="8"></rect>
                        </svg>
                        Print/Save PDF
                    </button>
                </div>
            </div>
            {% endif %}
        </main>

        <!-- Footer -->
        <footer class="footer">
            <div class="footer-inner">
                <span class="footer-text">© 2025 GasConsult.ai</span>
                <div class="footer-links">
                    <a href="/privacy" class="footer-link">Privacy</a>
                    <a href="/terms" class="footer-link">Terms</a>
                    <a href="mailto:contact@gasconsult.ai" class="footer-link">Contact</a>
                </div>
            </div>
        </footer>
    </div>

    <script>
        // Mobile menu toggle
        function toggleMobileMenu() {
            const menu = document.getElementById('mobileMenu');
            const btn = document.querySelector('.mobile-menu-btn');
            if (menu && btn) {
                menu.classList.toggle('active');
                btn.classList.toggle('active');
            }
        }

        function toggleNavDropdown(e) {
            e.preventDefault();
            e.stopPropagation();
            const menu = e.target.nextElementSibling;
            if (menu) {
                menu.classList.toggle('show');
            }
        }

        // Close dropdown when clicking outside
        document.addEventListener('click', function() {
            document.querySelectorAll('.nav-dropdown-menu').forEach(m => m.classList.remove('show'));
        });

        // Radio button selection handler
        document.querySelectorAll('.radio-option').forEach(option => {
            option.addEventListener('click', function() {
                const radio = this.querySelector('input[type="radio"]');
                if (radio) {
                    radio.checked = true;
                    // Remove selected class from siblings
                    const name = radio.getAttribute('name');
                    document.querySelectorAll(`input[name="${name}"]`).forEach(r => {
                        r.closest('.radio-option').classList.remove('selected');
                    });
                    // Add selected class to clicked option
                    this.classList.add('selected');
                }
            });
        });

        // Form submission loading state
        document.querySelector('form')?.addEventListener('submit', function(e) {
            const submitBtn = this.querySelector('.submit-btn');
            if (submitBtn) {
                submitBtn.textContent = 'Generating Assessment...';
                submitBtn.style.opacity = '0.7';
                submitBtn.disabled = true;
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
                'practice guideline', 'clinical practice', 'society statement', 'expert consensus',
                'position statement', 'standards of care', 'best practice', 'clinical pathway']):
            guideline_count += 1
            quality_score += 4

        # Check for meta-analysis
        if any(keyword in title or keyword in journal for keyword in
               ['meta-analysis', 'metaanalysis', 'meta analysis', 'pooled analysis',
                'network meta-analysis', 'meta-regression', 'systematic review and meta-analysis',
                'quantitative synthesis', 'evidence synthesis']):
            meta_analysis_count += 1
            quality_score += 3

        # Check for systematic review (separate from meta-analysis)
        elif any(keyword in title or keyword in journal for keyword in
                 ['systematic review', 'cochrane review', 'systematic literature review',
                  'evidence-based review', 'integrative review', 'scoping review',
                  'umbrella review', 'critical review']):
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

        # Check for review articles (general) - increased value since they're clinically useful
        elif any(keyword in title for keyword in ['review', 'update', 'overview', 'state of the art',
                                                   'current concepts', 'advances in', 'recent advances']):
            review_count += 1
            quality_score += 1.5  # Increased from 1.0 since reviews are still valuable

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
    # Adjusted thresholds to be more realistic while maintaining quality standards
    if (total_score >= 6 or  # Lowered from 8 (e.g., 1 guideline + 1 systematic review)
        meta_analysis_count >= 2 or
        (guideline_count >= 1 and meta_analysis_count >= 1) or
        (guideline_count >= 1 and rct_count >= 2) or
        (systematic_review_count >= 2 and num_papers >= 5) or  # 2+ systematic reviews with good paper count
        (meta_analysis_count >= 1 and systematic_review_count >= 1) or  # 1 meta + 1 systematic
        (num_papers >= 8 and total_score >= 4)):  # Many papers with decent quality

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
    # More lenient thresholds to avoid "Low Confidence" for good answers
    elif (total_score >= 2.5 or  # Lowered from 4 (e.g., 1 systematic review OR 1 RCT + 1 observational)
          meta_analysis_count >= 1 or
          systematic_review_count >= 1 or
          rct_count >= 2 or
          rct_count >= 1 and num_papers >= 4 or  # 1 RCT + multiple supporting papers
          guideline_count >= 1 or
          num_papers >= 5 and total_score >= 1.5):  # 5+ papers with some quality

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
    <!-- Google Analytics 4 -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=G-01NZYD1DPP"></script>
    <script>
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-01NZYD1DPP');
    </script>

    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GasConsult.ai - AI-Powered Anesthesiology Assistant</title>

    <!-- SEO Meta Tags -->
    <meta name="description" content="Evidence-based anesthesiology AI consultant combining PubMed research with GPT-4o for hallucination-free, citation-backed clinical answers.">
    <meta name="keywords" content="anesthesiology, anesthesia, medical AI, PubMed, clinical guidelines, evidence-based medicine, GPT-4o, systematic reviews">

    <!-- Open Graph / Facebook -->
    <meta property="og:type" content="website">
    <meta property="og:url" content="https://gasconsult.ai/">
    <meta property="og:title" content="GasConsult.ai - Evidence-Based Anesthesiology AI">
    <meta property="og:description" content="Get hallucination-free anesthesiology answers backed by PubMed research. GPT-4o + systematic reviews, RCTs, and clinical guidelines.">
    <meta property="og:image" content="https://gasconsult.ai/static/logo.png">

    <!-- Twitter -->
    <meta property="twitter:card" content="summary_large_image">
    <meta property="twitter:url" content="https://gasconsult.ai/">
    <meta property="twitter:title" content="GasConsult.ai - Evidence-Based Anesthesiology AI">
    <meta property="twitter:description" content="Get hallucination-free anesthesiology answers backed by PubMed research.">
    <meta property="twitter:image" content="https://gasconsult.ai/static/logo.png">

    <!-- PWA -->
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=6">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=6">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">

    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap" rel="stylesheet">

    <!-- Structured Data for SEO -->
    <script type="application/ld+json">
    {
      "@context": "https://schema.org",
      "@type": "MedicalWebPage",
      "name": "GasConsult.ai",
      "description": "Evidence-based anesthesiology AI consultant combining PubMed research with GPT-4o",
      "specialty": "Anesthesiology",
      "audience": {
        "@type": "MedicalAudience",
        "audienceType": "Healthcare Professional"
      }
    }
    </script>

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

        /* Skip to Content Link for Accessibility */
        .skip-to-content {
            position: absolute;
            top: -40px;
            left: 0;
            background: var(--blue-600);
            color: white;
            padding: 8px 16px;
            text-decoration: none;
            border-radius: 0 0 4px 0;
            z-index: 1000;
            font-weight: 600;
            transition: top 0.2s;
        }

        .skip-to-content:focus {
            top: 0;
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

        .nav-dropdown:has(.nav-dropdown-link.active) .nav-dropdown-toggle {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .nav-dropdown {
            position: relative;
            display: inline-block;
        }

        .nav-dropdown-toggle {
            cursor: pointer;
            background: none;
            border: none;
            font-family: inherit;
        }

        .nav-dropdown-menu {
            display: none;
            position: absolute;
            top: 100%;
            right: 0;
            background: white;
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            min-width: 200px;
            margin-top: 4px;
            z-index: 1000;
            overflow: hidden;
        }

        .nav-dropdown-menu.show {
            display: block;
        }

        .nav-dropdown-link {
            display: block;
            padding: 12px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .nav-dropdown-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
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
            padding: 3px;
            display: flex;
            align-items: flex-end;
            gap: 3px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04), inset 0 1px 0 rgba(255,255,255,1);
        }

        .chat-input {
            flex: 1;
            border: none;
            outline: none;
            padding: 8px 12px;
            font-size: 16px;
            font-family: inherit;
            color: var(--gray-800);
            background: transparent;
            resize: none;
            min-height: 32px;
            max-height: 110px;
            line-height: 1.5;
        }

        .chat-input::placeholder { color: var(--gray-400); }

        .chat-send {
            width: 40px;
            height: 40px;
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
            .chat-inner { border-radius: 16px; padding: 4px; }
            .chat-input { padding: 10px 14px; min-height: 40px; }
            .chat-send { width: 44px; height: 44px; border-radius: 12px; }
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
            .chat-inner { border-radius: 18px; }
            .chat-input { padding: 12px 16px; min-height: 44px; max-height: 160px; }
            .chat-send { width: 48px; height: 48px; border-radius: 14px; margin: 6px; }
            .chat-send svg { width: 20px; height: 20px; }
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
            margin-bottom: 32px;
            animation: fade-up 0.4s ease;
        }

        /* User Message - Clean, No Bubble Design */
        .user-message {
            display: flex;
            flex-direction: column;
            align-items: flex-end;
            width: 100%;
        }

        .user-message .message-bubble {
            background: transparent;
            color: var(--gray-900);
            padding: 0 0 20px 0;
            border-radius: 0;
            max-width: 100%;
            font-size: 15px;
            line-height: 1.6;
            box-shadow: none;
            position: relative;
            width: 100%;
            text-align: right;
        }

        .user-message .message-bubble::before {
            content: 'You';
            display: block;
            font-size: 11px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: var(--blue-600);
            margin-bottom: 8px;
        }

        .user-message .message-bubble::after {
            content: '';
            display: block;
            width: 100%;
            height: 1px;
            background: linear-gradient(90deg, transparent 0%, var(--gray-200) 50%, transparent 100%);
            margin-top: 16px;
            position: absolute;
            bottom: 0;
            left: 0;
        }

        /* AI Message - Full Width Glass Card, No Bubble */
        .ai-message {
            display: flex;
            flex-direction: column;
            align-items: stretch;
            width: 100%;
        }

        .ai-message .message-bubble {
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.8);
            padding: 24px;
            border-radius: 12px;
            max-width: 100%;
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

        .reference-number {
            font-size: 13px;
            font-weight: 700;
            color: var(--blue-600);
            min-width: 32px;
            flex-shrink: 0;
            padding-top: 1px;
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

        /* Skeleton Loader */
        .skeleton-loader {
            padding: 24px;
            background: rgba(255,255,255,0.7);
            backdrop-filter: blur(20px) saturate(180%);
            border-radius: 12px;
            margin-bottom: 12px;
        }

        .skeleton-line {
            height: 14px;
            background: linear-gradient(90deg, #e2e8f0 25%, #f1f5f9 50%, #e2e8f0 75%);
            background-size: 200% 100%;
            border-radius: 4px;
            animation: skeleton-pulse 1.5s ease-in-out infinite;
            margin-bottom: 10px;
        }

        .skeleton-line:last-child {
            margin-bottom: 0;
            width: 70%;
        }

        @keyframes skeleton-pulse {
            0% { background-position: 200% 0; }
            100% { background-position: -200% 0; }
        }

        .progress-indicator {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            font-size: 12px;
            color: var(--gray-600);
            margin-top: 12px;
            padding: 8px 12px;
            background: rgba(255, 255, 255, 0.5);
            border-radius: 8px;
            backdrop-filter: blur(10px);
        }

        .progress-indicator .progress-icon {
            width: 14px;
            height: 14px;
            border: 2px solid var(--blue-600);
            border-top-color: transparent;
            border-radius: 50%;
            animation: spin 0.8s linear infinite;
        }

        @keyframes spin {
            to { transform: rotate(360deg); }
        }

        .chat-input-area {
            position: sticky;
            bottom: 0;
            background: var(--gray-50);
            border-top: 1px solid var(--gray-200);
            padding: 6px;
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
            gap: 6px;
            padding: 8px 14px;
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
                margin-bottom: 40px;
            }

            .user-message .message-bubble {
                padding: 0 0 24px 0;
                font-size: 16px;
            }

            .ai-message .message-bubble {
                padding: 28px 32px;
                font-size: 16px;
                border-radius: 14px;
            }

            .chat-input-area {
                padding: 8px 20px;
            }

            .new-chat-btn {
                padding: 10px 18px;
                font-size: 14px;
            }
        }

        @media (min-width: 1024px) {
            .messages-container {
                padding: 40px 32px 60px;
            }

            .message {
                margin-bottom: 48px;
            }

            .user-message .message-bubble {
                font-size: 17px;
            }

            .ai-message .message-bubble {
                padding: 32px 36px;
                font-size: 17px;
                border-radius: 16px;
            }

            .chat-input-area {
                padding: 10px 28px;
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
    <a href="#main-content" class="skip-to-content">Skip to main content</a>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <nav class="nav" role="navigation" aria-label="Main navigation">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo" aria-label="GasConsult.ai home">
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
                    <div class="nav-dropdown">
                        <button class="nav-link nav-dropdown-toggle" onclick="toggleNavDropdown(event)">More ▼</button>
                        <div class="nav-dropdown-menu">
                            <a href="/hypotension" class="nav-dropdown-link">IOH Predictor</a>
                            <a href="/difficult-airway" class="nav-dropdown-link">Difficult Airway</a>
                            <a href="/informed-consent" class="nav-dropdown-link">Informed Consent</a>
                        </div>
                    </div>
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
            <a href="/difficult-airway" class="mobile-menu-link">Difficult Airway</a>
            <a href="/informed-consent" class="mobile-menu-link">Informed Consent</a>
        </div>

        {% if messages and messages|length > 0 %}
        <!-- Chat Interface -->
        <section class="chat-view" id="main-content">
            <div class="messages-container" id="messagesContainer" role="log" aria-live="polite" aria-label="Conversation history">
                {% if error_message %}
                <div class="error-message" role="alert">{{ error_message|safe }}</div>
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
                                        <div class="reference-number">[{{ loop.index }}]</div>
                                        <div class="reference-link-container">
                                            <a href="https://pubmed.ncbi.nlm.nih.gov/{{ ref.pmid }}/" target="_blank" rel="noopener noreferrer" class="reference-link">
                                                {{ ref.title }}
                                            </a>
                                            <div class="reference-meta">
                                                {{ ref.authors }} - {{ ref.journal }}, {{ ref.year }}
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
                        <div id="skeletonLoader" class="skeleton-loader">
                            <div class="skeleton-line"></div>
                            <div class="skeleton-line"></div>
                            <div class="skeleton-line"></div>
                            <div class="skeleton-line"></div>
                        </div>
                        <div class="message-content" id="streamingContent" style="display: none;"></div>
                        <div class="streaming-indicator">
                            <div class="streaming-dots">
                                <span></span>
                                <span></span>
                                <span></span>
                            </div>
                            <span id="progressText">Searching PubMed...</span>
                        </div>
                        <div class="references-section" id="streamingReferences" style="display: none;"></div>
                    </div>
                </div>
                {% endif %}
            </div>

            <div class="chat-input-area">
                <div class="chat-input-wrapper">
                    <form method="post" action="/" id="chat-form" class="chat-card" role="form" aria-label="Follow-up question form">
                        <input type="hidden" name="csrf_token" value="{{ csrf_token() }}"/>
                        <div class="chat-inner">
                            <textarea name="query"
                                      id="chat-query-input"
                                      class="chat-input"
                                      placeholder="Ask a follow-up question..."
                                      rows="1"
                                      required
                                      aria-label="Ask a follow-up question"></textarea>
                            <button type="submit" class="chat-send" aria-label="Submit follow-up question">
                                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" aria-hidden="true">
                                    <line x1="12" y1="19" x2="12" y2="5"></line>
                                    <polyline points="5 12 12 5 19 12"></polyline>
                                </svg>
                            </button>
                        </div>
                    </form>
                    <a href="/clear" class="new-chat-btn" aria-label="Start a new chat conversation">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" aria-hidden="true">
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
        <section class="hero" id="main-content">
            <div class="hero-badge">
                <span class="badge-dot"></span>
                <span class="badge-text">PubMed-Powered AI</span>
            </div>
            <h1 class="hero-title">The AI copilot for<br><span class="gradient">anesthesiology</span></h1>
            <p class="hero-subtitle">Evidence-based clinical decision support, instant drug dosing, and intelligent pre-op assessments - all in one place.</p>

            <div class="chat-container">
                <form method="post" action="/" class="chat-card" role="search" aria-label="Ask anesthesiology question">
                    <input type="hidden" name="csrf_token" value="{{ csrf_token() }}"/>
                    <div class="chat-inner">
                        <textarea name="query"
                                  id="query-input"
                                  class="chat-input"
                                  placeholder="Ask anything about anesthesiology..."
                                  rows="1"
                                  required
                                  aria-label="Ask your anesthesiology question"
                                  aria-describedby="search-hint"></textarea>
                        <button type="submit" class="chat-send" aria-label="Submit question">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" aria-hidden="true">
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
                    <a href="/crisis" class="feature-link">View protocols <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><line x1="5" y1="12" x2="19" y2="12"></line><polyline points="12 5 19 12 12 19"></polyline></svg></a>
                </div>
                <div class="feature-card">
                    <div class="feature-icon cyan">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path><path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path><line x1="8" y1="7" x2="16" y2="7"></line><line x1="8" y1="11" x2="16" y2="11"></line></svg>
                    </div>
                    <h3 class="feature-title">PubMed Integration</h3>
                    <p class="feature-desc">Every answer is backed by real citations from peer-reviewed medical literature.</p>
                    <a href="/evidence" class="feature-link">Learn more <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><line x1="5" y1="12" x2="19" y2="12"></line><polyline points="12 5 19 12 12 19"></polyline></svg></a>
                </div>
            </div>
        </section>
        {% endif %}

        <footer class="footer">
            <div class="footer-inner">
                <span class="footer-text">© 2025 GasConsult.ai</span>
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

        function toggleNavDropdown(e) {
            e.preventDefault();
            e.stopPropagation();
            const menu = e.target.nextElementSibling;
            if (menu) {
                menu.classList.toggle('show');
            }
        }

        // Close dropdown when clicking outside
        document.addEventListener('click', function() {
            document.querySelectorAll('.nav-dropdown-menu').forEach(m => m.classList.remove('show'));
        });

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

            // CRITICAL: Scroll to bottom immediately when streaming starts
            // This ensures user sees the "Thinking..." indicator right away
            scrollToBottom();

            const eventSource = new EventSource('/stream?request_id=' + encodeURIComponent(requestId));
            let accumulatedMarkdown = ''; // Accumulate markdown during streaming

            // Aggressive auto-scroll: Keep scrolling during streaming to ensure user sees new content
            // This runs every 100ms to aggressively follow the streaming content
            const scrollInterval = setInterval(function() {
                scrollToBottom();
            }, 100);

            eventSource.addEventListener('message', function(e) {
                try {
                    const event = JSON.parse(e.data);

                    if (event.type === 'connected') {
                        // Initial connection established - scroll to show we're ready
                        const progressText = document.getElementById('progressText');
                        if (progressText) progressText.textContent = 'Analyzing evidence...';
                        scrollToBottom();
                    } else if (event.type === 'content') {
                        // Stream content chunks
                        if (event.data) {
                            // Hide skeleton loader and show content on first chunk
                            const skeletonLoader = document.getElementById('skeletonLoader');
                            if (skeletonLoader && skeletonLoader.style.display !== 'none') {
                                skeletonLoader.style.display = 'none';
                                const progressText = document.getElementById('progressText');
                                if (progressText) progressText.textContent = 'Streaming answer...';
                            }
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
                            // Scroll after each content chunk to keep new text visible
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

                            event.data.forEach(function(ref, index) {
                                refsHTML += '<div class="reference-item">';
                                refsHTML += '<div class="reference-header">';
                                refsHTML += '<div class="reference-number">[' + (index + 1) + ']</div>';
                                refsHTML += '<div class="reference-link-container">';
                                refsHTML += '<a href="https://pubmed.ncbi.nlm.nih.gov/' + ref.pmid + '/" target="_blank" rel="noopener noreferrer" class="reference-link">';
                                refsHTML += ref.title;
                                refsHTML += '</a>';
                                refsHTML += '<div class="reference-meta">';
                                refsHTML += ref.authors + ' - ' + ref.journal + ', ' + ref.year;
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

                        // Stop the aggressive auto-scroll interval
                        clearInterval(scrollInterval);

                        eventSource.close();
                        scrollToBottom();
                    } else if (event.type === 'error') {
                        console.error('[SSE] Server error:', event.message);

                        if (streamingIndicator) {
                            streamingIndicator.innerHTML = '<span style="color: #DC2626;">Error: ' + (event.message || 'Unknown error') + '</span>';
                        }

                        // Stop the aggressive auto-scroll interval
                        clearInterval(scrollInterval);

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

                // Stop the aggressive auto-scroll interval
                clearInterval(scrollInterval);

                eventSource.close();
            });

            // Scroll to bottom on page load
            setTimeout(scrollToBottom, 100);
        })();
        {% endif %}

        // Parse markdown in existing messages on page load
        // IMPORTANT: Only parse if content is plain text/markdown, not already HTML
        window.addEventListener('DOMContentLoaded', function() {
            if (typeof marked !== 'undefined') {
                // Find all message-content divs that contain markdown
                const messageContents = document.querySelectorAll('.ai-message .message-content');
                messageContents.forEach(function(element) {
                    // Check if content is already HTML (from session) or plain text/markdown
                    const html = element.innerHTML.trim();

                    // Skip if already contains HTML tags (loaded from session)
                    // Messages from session already have HTML, don't re-parse!
                    if (html.includes('<h') || html.includes('<p') || html.includes('<ul') || html.includes('<strong')) {
                        return;  // Already HTML, skip parsing
                    }

                    // Get the text content (markdown) for plain text messages
                    const markdownText = element.textContent;
                    if (markdownText && markdownText.trim()) {
                        element.innerHTML = marked.parse(markdownText);
                    }
                });
            }
        });

        // Auto-scroll to bottom on page load (for chat view)
        {% if messages and messages|length > 0 %}
        // Scroll immediately when DOM is ready
        window.addEventListener('DOMContentLoaded', function() {
            scrollToBottom();
        });

        // Scroll again after everything loads (images, etc)
        window.addEventListener('load', function() {
            scrollToBottom();
        });

        // Scroll one more time after a brief delay to ensure pending stream is visible
        {% if pending_stream %}
        setTimeout(function() {
            scrollToBottom();
        }, 100);
        {% endif %}
        {% endif %}
    </script>
</body>
</html>

"""


LIBRARY_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <!-- Google Analytics 4 -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=G-01NZYD1DPP"></script>
    <script>
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-01NZYD1DPP');
    </script>

    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Library - gasconsult.ai</title>

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

        /* Skip to Content Link for Accessibility */
        .skip-to-content {
            position: absolute;
            top: -40px;
            left: 0;
            background: var(--blue-600);
            color: white;
            padding: 8px 16px;
            text-decoration: none;
            border-radius: 0 0 4px 0;
            z-index: 1000;
            font-weight: 600;
            transition: top 0.2s;
        }

        .skip-to-content:focus {
            top: 0;
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

        .nav-dropdown:has(.nav-dropdown-link.active) .nav-dropdown-toggle {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .nav-dropdown {
            position: relative;
            display: inline-block;
        }

        .nav-dropdown-toggle {
            cursor: pointer;
            background: none;
            border: none;
            font-family: inherit;
        }

        .nav-dropdown-menu {
            display: none;
            position: absolute;
            top: 100%;
            right: 0;
            background: white;
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            min-width: 200px;
            margin-top: 4px;
            z-index: 1000;
            overflow: hidden;
        }

        .nav-dropdown-menu.show {
            display: block;
        }

        .nav-dropdown-link {
            display: block;
            padding: 12px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .nav-dropdown-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
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
            padding: 3px;
            display: flex;
            align-items: flex-end;
            gap: 3px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04), inset 0 1px 0 rgba(255,255,255,1);
        }

        .chat-input {
            flex: 1;
            border: none;
            outline: none;
            padding: 8px 12px;
            font-size: 16px;
            font-family: inherit;
            color: var(--gray-800);
            background: transparent;
            resize: none;
            min-height: 32px;
            max-height: 110px;
            line-height: 1.5;
        }

        .chat-input::placeholder { color: var(--gray-400); }

        .chat-send {
            width: 40px;
            height: 40px;
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
            .chat-inner { border-radius: 16px; padding: 4px; }
            .chat-input { padding: 10px 14px; min-height: 40px; }
            .chat-send { width: 44px; height: 44px; border-radius: 12px; }
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
            .chat-inner { border-radius: 18px; }
            .chat-input { padding: 12px 16px; min-height: 44px; max-height: 160px; }
            .chat-send { width: 48px; height: 48px; border-radius: 14px; margin: 6px; }
            .chat-send svg { width: 20px; height: 20px; }
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
    <a href="#main-content" class="skip-to-content">Skip to main content</a>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <nav class="nav" role="navigation" aria-label="Main navigation">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo" aria-label="GasConsult.ai home">
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
                    <a href="/crisis" class="nav-link">Crisis Protocols</a>
                    <div class="nav-dropdown">
                        <button class="nav-link nav-dropdown-toggle" onclick="toggleNavDropdown(event)">More ▼</button>
                        <div class="nav-dropdown-menu">
                            <a href="/hypotension" class="nav-dropdown-link">IOH Predictor</a>
                            <a href="/difficult-airway" class="nav-dropdown-link">Difficult Airway</a>
                            <a href="/informed-consent" class="nav-dropdown-link">Informed Consent</a>
                        </div>
                    </div>
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
            <a href="/difficult-airway" class="mobile-menu-link">Difficult Airway</a>
            <a href="/informed-consent" class="mobile-menu-link">Informed Consent</a>
        </div>

        <main class="main-content">
<main>
        <h1><svg style="width:14px;height:14px;margin-right:2px;" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path><path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path><path d="M8 7h8"></path><path d="M8 11h8"></path></svg> My Library</h1>

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

        function toggleNavDropdown(e) {
            e.preventDefault();
            e.stopPropagation();
            const menu = e.target.nextElementSibling;
            if (menu) {
                menu.classList.toggle('show');
            }
        }

        // Close dropdown when clicking outside
        document.addEventListener('click', function() {
            document.querySelectorAll('.nav-dropdown-menu').forEach(m => m.classList.remove('show'));
        });
    </script>
        </main>
    </div>
</body>
</html>
"""

SHARED_RESPONSE_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <!-- Google Analytics 4 -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=G-01NZYD1DPP"></script>
    <script>
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-01NZYD1DPP');
    </script>

    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Shared Response - gasconsult.ai</title>

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

        /* Skip to Content Link for Accessibility */
        .skip-to-content {
            position: absolute;
            top: -40px;
            left: 0;
            background: var(--blue-600);
            color: white;
            padding: 8px 16px;
            text-decoration: none;
            border-radius: 0 0 4px 0;
            z-index: 1000;
            font-weight: 600;
            transition: top 0.2s;
        }

        .skip-to-content:focus {
            top: 0;
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

        .nav-dropdown:has(.nav-dropdown-link.active) .nav-dropdown-toggle {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .nav-dropdown {
            position: relative;
            display: inline-block;
        }

        .nav-dropdown-toggle {
            cursor: pointer;
            background: none;
            border: none;
            font-family: inherit;
        }

        .nav-dropdown-menu {
            display: none;
            position: absolute;
            top: 100%;
            right: 0;
            background: white;
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            min-width: 200px;
            margin-top: 4px;
            z-index: 1000;
            overflow: hidden;
        }

        .nav-dropdown-menu.show {
            display: block;
        }

        .nav-dropdown-link {
            display: block;
            padding: 12px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .nav-dropdown-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
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
            padding: 3px;
            display: flex;
            align-items: flex-end;
            gap: 3px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04), inset 0 1px 0 rgba(255,255,255,1);
        }

        .chat-input {
            flex: 1;
            border: none;
            outline: none;
            padding: 8px 12px;
            font-size: 16px;
            font-family: inherit;
            color: var(--gray-800);
            background: transparent;
            resize: none;
            min-height: 32px;
            max-height: 110px;
            line-height: 1.5;
        }

        .chat-input::placeholder { color: var(--gray-400); }

        .chat-send {
            width: 40px;
            height: 40px;
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
            .chat-inner { border-radius: 16px; padding: 4px; }
            .chat-input { padding: 10px 14px; min-height: 40px; }
            .chat-send { width: 44px; height: 44px; border-radius: 12px; }
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
            .chat-inner { border-radius: 18px; }
            .chat-input { padding: 12px 16px; min-height: 44px; max-height: 160px; }
            .chat-send { width: 48px; height: 48px; border-radius: 14px; margin: 6px; }
            .chat-send svg { width: 20px; height: 20px; }
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
    <a href="#main-content" class="skip-to-content">Skip to main content</a>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <nav class="nav" role="navigation" aria-label="Main navigation">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo" aria-label="GasConsult.ai home">
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
                    <a href="/crisis" class="nav-link">Crisis Protocols</a>
                    <div class="nav-dropdown">
                        <button class="nav-link nav-dropdown-toggle" onclick="toggleNavDropdown(event)">More ▼</button>
                        <div class="nav-dropdown-menu">
                            <a href="/hypotension" class="nav-dropdown-link">IOH Predictor</a>
                            <a href="/difficult-airway" class="nav-dropdown-link">Difficult Airway</a>
                            <a href="/informed-consent" class="nav-dropdown-link">Informed Consent</a>
                        </div>
                    </div>
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
            <a href="/difficult-airway" class="mobile-menu-link">Difficult Airway</a>
            <a href="/informed-consent" class="mobile-menu-link">Informed Consent</a>
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
    <!-- Google Analytics 4 -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=G-01NZYD1DPP"></script>
    <script>
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-01NZYD1DPP');
    </script>

    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Terms of Service - gasconsult.ai</title>

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

        /* Skip to Content Link for Accessibility */
        .skip-to-content {
            position: absolute;
            top: -40px;
            left: 0;
            background: var(--blue-600);
            color: white;
            padding: 8px 16px;
            text-decoration: none;
            border-radius: 0 0 4px 0;
            z-index: 1000;
            font-weight: 600;
            transition: top 0.2s;
        }

        .skip-to-content:focus {
            top: 0;
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

        .nav-dropdown:has(.nav-dropdown-link.active) .nav-dropdown-toggle {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .nav-dropdown {
            position: relative;
            display: inline-block;
        }

        .nav-dropdown-toggle {
            cursor: pointer;
            background: none;
            border: none;
            font-family: inherit;
        }

        .nav-dropdown-menu {
            display: none;
            position: absolute;
            top: 100%;
            right: 0;
            background: white;
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            min-width: 200px;
            margin-top: 4px;
            z-index: 1000;
            overflow: hidden;
        }

        .nav-dropdown-menu.show {
            display: block;
        }

        .nav-dropdown-link {
            display: block;
            padding: 12px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .nav-dropdown-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
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
    <a href="#main-content" class="skip-to-content">Skip to main content</a>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <nav class="nav" role="navigation" aria-label="Main navigation">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo" aria-label="GasConsult.ai home">
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
                    <a href="/crisis" class="nav-link">Crisis Protocols</a>
                    <div class="nav-dropdown">
                        <button class="nav-link nav-dropdown-toggle" onclick="toggleNavDropdown(event)">More ▼</button>
                        <div class="nav-dropdown-menu">
                            <a href="/hypotension" class="nav-dropdown-link">IOH Predictor</a>
                            <a href="/difficult-airway" class="nav-dropdown-link">Difficult Airway</a>
                            <a href="/informed-consent" class="nav-dropdown-link">Informed Consent</a>
                        </div>
                    </div>
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
            <a href="/difficult-airway" class="mobile-menu-link">Difficult Airway</a>
            <a href="/informed-consent" class="mobile-menu-link">Informed Consent</a>
        </div>

        <main class="main-content">
            <div class="content-card">
                <h1>Terms of Service</h1>
                <p class="last-updated">Last Updated: December 2025</p>

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

                <h2>2. NOT Medical Advice - Critical Disclaimers</h2>

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

                <h2>4. Emergency Disclaimer - DO NOT USE FOR EMERGENCIES</h2>

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

                <h2>5. Educational Use Only - Professional Verification Required</h2>
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

                <h2>8. No Warranties - Service Provided "AS IS"</h2>

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

                <h2>9. Limitation of Liability - Maximum Legal Protection</h2>

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

                <h2>10. Indemnification - You Hold Us Harmless</h2>
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

                <h2>11. Assumption of Risk - You Accept All Risks</h2>
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

                <h2>14. Dispute Resolution - Mandatory Binding Arbitration</h2>

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
                <span class="footer-text">© 2025 GasConsult.ai</span>
                <div class="footer-links">
                    <a href="/privacy" class="footer-link">Privacy</a>
                    <a href="/terms" class="footer-link">Terms</a>
                    <a href="mailto:contact@gasconsult.ai" class="footer-link">Contact</a>
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

        function toggleNavDropdown(e) {
            e.preventDefault();
            e.stopPropagation();
            const menu = e.target.nextElementSibling;
            if (menu) {
                menu.classList.toggle('show');
            }
        }

        // Close dropdown when clicking outside
        document.addEventListener('click', function() {
            document.querySelectorAll('.nav-dropdown-menu').forEach(m => m.classList.remove('show'));
        });
    </script>
</body>
</html>
"""

PRIVACY_POLICY_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <!-- Google Analytics 4 -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=G-01NZYD1DPP"></script>
    <script>
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-01NZYD1DPP');
    </script>

    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Privacy Policy - gasconsult.ai</title>

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

        /* Skip to Content Link for Accessibility */
        .skip-to-content {
            position: absolute;
            top: -40px;
            left: 0;
            background: var(--blue-600);
            color: white;
            padding: 8px 16px;
            text-decoration: none;
            border-radius: 0 0 4px 0;
            z-index: 1000;
            font-weight: 600;
            transition: top 0.2s;
        }

        .skip-to-content:focus {
            top: 0;
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

        .nav-dropdown:has(.nav-dropdown-link.active) .nav-dropdown-toggle {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .nav-dropdown {
            position: relative;
            display: inline-block;
        }

        .nav-dropdown-toggle {
            cursor: pointer;
            background: none;
            border: none;
            font-family: inherit;
        }

        .nav-dropdown-menu {
            display: none;
            position: absolute;
            top: 100%;
            right: 0;
            background: white;
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            min-width: 200px;
            margin-top: 4px;
            z-index: 1000;
            overflow: hidden;
        }

        .nav-dropdown-menu.show {
            display: block;
        }

        .nav-dropdown-link {
            display: block;
            padding: 12px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .nav-dropdown-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
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
            .orb-1 { width: 600px; height: 600px; left: -10%; }
            .orb-2 { width: 450px; height: 450px; right: -10%; }
            .orb-3 { width: 400px; height: 400px; }
        }

        @media (min-width: 1024px) {
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
    <a href="#main-content" class="skip-to-content">Skip to main content</a>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <nav class="nav" role="navigation" aria-label="Main navigation">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo" aria-label="GasConsult.ai home">
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
                    <a href="/crisis" class="nav-link">Crisis Protocols</a>
                    <div class="nav-dropdown">
                        <button class="nav-link nav-dropdown-toggle" onclick="toggleNavDropdown(event)">More ▼</button>
                        <div class="nav-dropdown-menu">
                            <a href="/hypotension" class="nav-dropdown-link">IOH Predictor</a>
                            <a href="/difficult-airway" class="nav-dropdown-link">Difficult Airway</a>
                            <a href="/informed-consent" class="nav-dropdown-link">Informed Consent</a>
                        </div>
                    </div>
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
            <a href="/difficult-airway" class="mobile-menu-link">Difficult Airway</a>
            <a href="/informed-consent" class="mobile-menu-link">Informed Consent</a>
        </div>

        <main class="main-content">
            <div class="content-card">
                <h1>Privacy Policy</h1>
                <p class="last-updated">Last Updated: December 2025</p>

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
                <span class="footer-text">© 2025 GasConsult.ai</span>
                <div class="footer-links">
                    <a href="/privacy" class="footer-link">Privacy</a>
                    <a href="/terms" class="footer-link">Terms</a>
                    <a href="mailto:contact@gasconsult.ai" class="footer-link">Contact</a>
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

        function toggleNavDropdown(e) {
            e.preventDefault();
            e.stopPropagation();
            const menu = e.target.nextElementSibling;
            if (menu) {
                menu.classList.toggle('show');
            }
        }

        // Close dropdown when clicking outside
        document.addEventListener('click', function() {
            document.querySelectorAll('.nav-dropdown-menu').forEach(m => m.classList.remove('show'));
        });
    </script>
</body>
</html>
"""

EVIDENCE_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <!-- Google Analytics 4 -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=G-01NZYD1DPP"></script>
    <script>
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-01NZYD1DPP');
    </script>

    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Evidence-Based Methodology - gasconsult.ai</title>

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
            --green-50: #ECFDF5;
            --green-500: #10B981;
            --green-600: #059669;
            --amber-50: #FFFBEB;
            --amber-500: #F59E0B;
            --amber-600: #D97706;
            --purple-50: #FAF5FF;
            --purple-500: #A855F7;
            --purple-600: #9333EA;
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

        /* Skip to Content Link for Accessibility */
        .skip-to-content {
            position: absolute;
            top: -40px;
            left: 0;
            background: var(--blue-600);
            color: white;
            padding: 8px 16px;
            text-decoration: none;
            border-radius: 0 0 4px 0;
            z-index: 1000;
            font-weight: 600;
            transition: top 0.2s;
        }

        .skip-to-content:focus {
            top: 0;
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

        .nav-dropdown:has(.nav-dropdown-link.active) .nav-dropdown-toggle {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .nav-dropdown {
            position: relative;
            display: inline-block;
        }

        .nav-dropdown-toggle {
            cursor: pointer;
            background: none;
            border: none;
            font-family: inherit;
        }

        .nav-dropdown-menu {
            display: none;
            position: absolute;
            top: 100%;
            right: 0;
            background: white;
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            min-width: 200px;
            margin-top: 4px;
            z-index: 1000;
            overflow: hidden;
        }

        .nav-dropdown-menu.show {
            display: block;
        }

        .nav-dropdown-link {
            display: block;
            padding: 12px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .nav-dropdown-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
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

        /* Main Content */
        .main-content {
            flex: 1;
            padding: 100px 16px 40px;
            max-width: 1000px;
            margin: 0 auto;
            width: 100%;
        }

        /* Hero Section */
        .hero {
            text-align: center;
            margin-bottom: 48px;
        }

        .hero-badge {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            padding: 8px 16px;
            background: rgba(37, 99, 235, 0.1);
            border: 1px solid rgba(37, 99, 235, 0.2);
            border-radius: 100px;
            font-size: 12px;
            font-weight: 600;
            color: var(--blue-600);
            text-transform: uppercase;
            letter-spacing: 0.5px;
            margin-bottom: 20px;
        }

        .hero-badge svg {
            width: 16px;
            height: 16px;
        }

        .hero-title {
            font-size: 32px;
            font-weight: 800;
            letter-spacing: -1.5px;
            color: var(--gray-900);
            margin-bottom: 12px;
            line-height: 1.2;
        }

        .hero-title .gradient {
            background: linear-gradient(135deg, var(--blue-600) 0%, var(--blue-700) 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
        }

        .hero-subtitle {
            font-size: 15px;
            color: var(--gray-500);
            max-width: 700px;
            margin: 0 auto 28px;
            line-height: 1.6;
        }

        /* Content Cards */
        .content-card {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(40px) saturate(180%);
            -webkit-backdrop-filter: blur(40px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 24px;
            padding: 32px 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
            margin-bottom: 24px;
        }

        .content-card h2 {
            font-size: 22px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 16px;
            letter-spacing: -0.5px;
        }

        .content-card p {
            font-size: 15px;
            color: var(--gray-600);
            line-height: 1.7;
            margin-bottom: 16px;
        }

        .content-card p:last-child {
            margin-bottom: 0;
        }

        /* Study Hierarchy */
        .hierarchy-grid {
            display: grid;
            grid-template-columns: 1fr;
            gap: 16px;
            margin-top: 24px;
        }

        .hierarchy-item {
            background: rgba(255,255,255,0.6);
            border: 1px solid var(--gray-200);
            border-radius: 16px;
            padding: 20px;
            display: flex;
            align-items: flex-start;
            gap: 16px;
            transition: all 0.2s ease;
        }

        .hierarchy-item:hover {
            background: rgba(255,255,255,0.9);
            border-color: var(--gray-300);
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(0,0,0,0.06);
        }

        .hierarchy-rank {
            width: 40px;
            height: 40px;
            border-radius: 12px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 16px;
            font-weight: 700;
            flex-shrink: 0;
        }

        .hierarchy-rank.gold {
            background: linear-gradient(135deg, #F59E0B 0%, #D97706 100%);
            color: white;
            box-shadow: 0 4px 12px rgba(245, 158, 11, 0.3);
        }

        .hierarchy-rank.silver {
            background: linear-gradient(135deg, #6B7280 0%, #4B5563 100%);
            color: white;
            box-shadow: 0 4px 12px rgba(107, 114, 128, 0.3);
        }

        .hierarchy-rank.bronze {
            background: linear-gradient(135deg, #CD7F32 0%, #8B5A2B 100%);
            color: white;
            box-shadow: 0 4px 12px rgba(205, 127, 50, 0.3);
        }

        .hierarchy-rank.standard {
            background: linear-gradient(135deg, var(--blue-500) 0%, var(--blue-600) 100%);
            color: white;
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.3);
        }

        .hierarchy-content h3 {
            font-size: 16px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 6px;
        }

        .hierarchy-content p {
            font-size: 14px;
            color: var(--gray-600);
            line-height: 1.6;
            margin: 0;
        }

        /* Features Grid */
        .features-grid {
            display: grid;
            grid-template-columns: 1fr;
            gap: 16px;
            margin-top: 24px;
        }

        .feature-box {
            background: rgba(255,255,255,0.6);
            border: 1px solid var(--gray-200);
            border-radius: 16px;
            padding: 24px;
            display: flex;
            align-items: flex-start;
            gap: 16px;
        }

        .feature-icon-box {
            width: 48px;
            height: 48px;
            border-radius: 12px;
            display: flex;
            align-items: center;
            justify-content: center;
            flex-shrink: 0;
        }

        .feature-icon-box.blue {
            background: linear-gradient(135deg, var(--blue-500) 0%, var(--blue-600) 100%);
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.3);
        }

        .feature-icon-box.green {
            background: linear-gradient(135deg, var(--green-500) 0%, var(--green-600) 100%);
            box-shadow: 0 4px 12px rgba(16, 185, 129, 0.3);
        }

        .feature-icon-box.purple {
            background: linear-gradient(135deg, var(--purple-500) 0%, var(--purple-600) 100%);
            box-shadow: 0 4px 12px rgba(168, 85, 247, 0.3);
        }

        .feature-icon-box.amber {
            background: linear-gradient(135deg, var(--amber-500) 0%, var(--amber-600) 100%);
            box-shadow: 0 4px 12px rgba(245, 158, 11, 0.3);
        }

        .feature-icon-box svg {
            width: 24px;
            height: 24px;
            stroke: white;
        }

        .feature-box-content h3 {
            font-size: 16px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 8px;
        }

        .feature-box-content p {
            font-size: 14px;
            color: var(--gray-600);
            line-height: 1.6;
            margin: 0;
        }

        /* Code Block */
        .code-block {
            background: var(--gray-900);
            border-radius: 12px;
            padding: 20px;
            margin: 24px 0;
            overflow-x: auto;
        }

        .code-block code {
            font-family: 'Monaco', 'Menlo', monospace;
            font-size: 13px;
            color: #93C5FD;
            line-height: 1.6;
        }

        /* Footer */
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

        .footer-link:hover {
            color: var(--gray-700);
        }

        /* Responsive */
        @media (min-width: 640px) {
            .hierarchy-grid {
                grid-template-columns: repeat(2, 1fr);
            }

            .features-grid {
                grid-template-columns: repeat(2, 1fr);
            }
        }

        @media (min-width: 768px) {
            .nav { padding: 16px 32px; }
            .nav-inner { height: 64px; padding: 0 24px; border-radius: 20px; }
            .logo-icon svg { width: 42px; height: 15px; }
            .logo-text { font-size: 20px; }
            .nav-links { display: flex; }
            .mobile-menu-btn { display: none; }

            .main-content { padding: 120px 32px 60px; }

            .hero-title { font-size: 48px; }
            .hero-subtitle { font-size: 17px; }

            .content-card { padding: 40px; }

            .footer { padding: 40px 32px; }
            .footer-inner { flex-direction: row; justify-content: space-between; text-align: left; }
            .footer-text { font-size: 14px; }
            .footer-links { gap: 32px; }
            .footer-link { font-size: 14px; }
        }

        @media (min-width: 1024px) {
            .nav { padding: 16px 40px; }
            .main-content { padding: 140px 40px 80px; }
        }
    </style>
</head>
<body>
    <a href="#main-content" class="skip-to-content">Skip to main content</a>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <!-- Navigation -->
        <nav class="nav" role="navigation" aria-label="Main navigation">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo" aria-label="GasConsult.ai home">
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
                    <a href="/crisis" class="nav-link">Crisis Protocols</a>
                    <div class="nav-dropdown">
                        <button class="nav-link nav-dropdown-toggle" onclick="toggleNavDropdown(event)">More ▼</button>
                        <div class="nav-dropdown-menu">
                            <a href="/hypotension" class="nav-dropdown-link">IOH Predictor</a>
                            <a href="/difficult-airway" class="nav-dropdown-link">Difficult Airway</a>
                            <a href="/informed-consent" class="nav-dropdown-link">Informed Consent</a>
                        </div>
                    </div>
                </div>
                <button class="mobile-menu-btn" onclick="toggleMobileMenu()">
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
            <a href="/difficult-airway" class="mobile-menu-link">Difficult Airway</a>
            <a href="/informed-consent" class="mobile-menu-link">Informed Consent</a>
        </div>

        <!-- Main Content -->
        <main class="main-content">
            <!-- Hero -->
            <div class="hero">
                <div class="hero-badge">
                    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                        <path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path>
                        <path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path>
                        <line x1="8" y1="7" x2="16" y2="7"></line>
                        <line x1="8" y1="11" x2="16" y2="11"></line>
                    </svg>
                    <span>Evidence-Based</span>
                </div>
                <h1 class="hero-title">
                    Our <span class="gradient">Methodology</span>
                </h1>
                <p class="hero-subtitle">
                    Every answer on GasConsult.ai is backed by peer-reviewed medical literature from PubMed, stratified by study quality and research strength
                </p>
            </div>

            <!-- PubMed Integration -->
            <div class="content-card">
                <h2>PubMed Integration: How It Works</h2>
                <p>
                    When you ask a clinical question, GasConsult.ai doesn't rely on pre-trained knowledge or outdated information. Instead, we perform real-time searches of PubMed's database of over 35 million biomedical citations to find the most relevant, high-quality evidence.
                </p>
                <p>
                    Our intelligent search algorithm automatically expands medical abbreviations, applies anesthesia-specific filters, and prioritizes the highest levels of evidence to ensure you receive answers grounded in current, peer-reviewed research.
                </p>
            </div>

            <!-- Study Quality Hierarchy -->
            <div class="content-card">
                <h2>Evidence Hierarchy: How We Rank Studies</h2>
                <p>
                    Not all medical evidence is created equal. We stratify search results using the established evidence pyramid, prioritizing study designs with the highest methodological rigor and least susceptibility to bias.
                </p>

                <div class="hierarchy-grid">
                    <div class="hierarchy-item">
                        <div class="hierarchy-rank gold">1</div>
                        <div class="hierarchy-content">
                            <h3>Systematic Reviews & Meta-Analyses</h3>
                            <p>Comprehensive synthesis of multiple studies with statistical pooling. Highest level of evidence for clinical questions.</p>
                        </div>
                    </div>

                    <div class="hierarchy-item">
                        <div class="hierarchy-rank gold">2</div>
                        <div class="hierarchy-content">
                            <h3>Cochrane Reviews</h3>
                            <p>Gold-standard systematic reviews with rigorous methodology and regular updates to reflect new evidence.</p>
                        </div>
                    </div>

                    <div class="hierarchy-item">
                        <div class="hierarchy-rank silver">3</div>
                        <div class="hierarchy-content">
                            <h3>Randomized Controlled Trials (RCTs)</h3>
                            <p>Prospective studies with randomization and control groups. Minimizes bias and establishes causality.</p>
                        </div>
                    </div>

                    <div class="hierarchy-item">
                        <div class="hierarchy-rank bronze">4</div>
                        <div class="hierarchy-content">
                            <h3>Clinical Practice Guidelines</h3>
                            <p>Evidence-based recommendations from professional societies (ASA, ESA, etc.) synthesizing best practices.</p>
                        </div>
                    </div>

                    <div class="hierarchy-item">
                        <div class="hierarchy-rank standard">5</div>
                        <div class="hierarchy-content">
                            <h3>Cohort & Case-Control Studies</h3>
                            <p>Observational studies that analyze associations and outcomes in patient populations.</p>
                        </div>
                    </div>

                    <div class="hierarchy-item">
                        <div class="hierarchy-rank standard">6</div>
                        <div class="hierarchy-content">
                            <h3>Case Series & Expert Reviews</h3>
                            <p>Descriptive studies and expert opinion. Useful when higher-quality evidence is unavailable.</p>
                        </div>
                    </div>
                </div>
            </div>

            <!-- Search Strategy -->
            <div class="content-card">
                <h2>Intelligent Search Strategy</h2>
                <p>
                    Our PubMed search employs sophisticated filtering to maximize relevance and evidence quality:
                </p>

                <div class="features-grid">
                    <div class="feature-box">
                        <div class="feature-icon-box blue">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <circle cx="11" cy="11" r="8"></circle>
                                <path d="m21 21-4.35-4.35"></path>
                            </svg>
                        </div>
                        <div class="feature-box-content">
                            <h3>Anesthesia-First Filtering</h3>
                            <p>Searches prioritize anesthesiology-specific MeSH terms and journals before expanding to general medicine.</p>
                        </div>
                    </div>

                    <div class="feature-box">
                        <div class="feature-icon-box green">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <polyline points="20 6 9 17 4 12"></polyline>
                            </svg>
                        </div>
                        <div class="feature-box-content">
                            <h3>Study Type Prioritization</h3>
                            <p>Filters for systematic reviews, meta-analyses, RCTs, and Cochrane reviews before other study types.</p>
                        </div>
                    </div>

                    <div class="feature-box">
                        <div class="feature-icon-box purple">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <rect x="3" y="4" width="18" height="18" rx="2" ry="2"></rect>
                                <line x1="16" y1="2" x2="16" y2="6"></line>
                                <line x1="8" y1="2" x2="8" y2="6"></line>
                                <line x1="3" y1="10" x2="21" y2="10"></line>
                            </svg>
                        </div>
                        <div class="feature-box-content">
                            <h3>Recency Weighting (2015+)</h3>
                            <p>Focuses on studies published since 2015 to ensure recommendations reflect current medical practice.</p>
                        </div>
                    </div>

                    <div class="feature-box">
                        <div class="feature-icon-box amber">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <path d="M21 16V8a2 2 0 0 0-1-1.73l-7-4a2 2 0 0 0-2 0l-7 4A2 2 0 0 0 3 8v8a2 2 0 0 0 1 1.73l7 4a2 2 0 0 0 2 0l7-4A2 2 0 0 0 21 16z"></path>
                                <polyline points="7.5 4.21 12 6.81 16.5 4.21"></polyline>
                                <polyline points="7.5 19.79 7.5 14.6 3 12"></polyline>
                                <polyline points="21 12 16.5 14.6 16.5 19.79"></polyline>
                                <polyline points="3.27 6.96 12 12.01 20.73 6.96"></polyline>
                                <line x1="12" y1="22.08" x2="12" y2="12"></line>
                            </svg>
                        </div>
                        <div class="feature-box-content">
                            <h3>Synonym Expansion</h3>
                            <p>Automatically expands medical abbreviations (e.g., "roc" → rocuronium, "TXA" → tranexamic acid) for comprehensive results.</p>
                        </div>
                    </div>
                </div>
            </div>

            <!-- Citation Quality -->
            <div class="content-card">
                <h2>Citation Transparency & Quality Indicators</h2>
                <p>
                    Every answer includes full citations with direct PubMed links, allowing you to verify sources instantly. We display confidence badges based on the quantity and quality of retrieved evidence:
                </p>

                <div class="features-grid">
                    <div class="feature-box">
                        <div class="feature-icon-box green">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <path d="M12 22c5.523 0 10-4.477 10-10S17.523 2 12 2 2 6.477 2 12s4.477 10 10 10z"></path>
                                <path d="m9 12 2 2 4-4"></path>
                            </svg>
                        </div>
                        <div class="feature-box-content">
                            <h3>High Confidence</h3>
                            <p>10+ high-quality papers including systematic reviews or meta-analyses. Evidence is robust and well-established.</p>
                        </div>
                    </div>

                    <div class="feature-box">
                        <div class="feature-icon-box amber">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <circle cx="12" cy="12" r="10"></circle>
                                <path d="M12 8v4"></path>
                                <path d="M12 16h.01"></path>
                            </svg>
                        </div>
                        <div class="feature-box-content">
                            <h3>Moderate Confidence</h3>
                            <p>5-10 papers with mixed study types. Evidence exists but may be less conclusive or require clinical judgment.</p>
                        </div>
                    </div>

                    <div class="feature-box">
                        <div class="feature-icon-box blue">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <path d="M10.29 3.86 1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z"></path>
                                <line x1="12" y1="9" x2="12" y2="13"></line>
                                <line x1="12" y1="17" x2="12.01" y2="17"></line>
                            </svg>
                        </div>
                        <div class="feature-box-content">
                            <h3>Low Confidence</h3>
                            <p>Fewer than 5 papers or limited study types. Evidence is emerging or question is poorly studied. Use with caution.</p>
                        </div>
                    </div>

                    <div class="feature-box">
                        <div class="feature-icon-box purple">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <path d="M21 15a2 2 0 0 1-2 2H7l-4 4V5a2 2 0 0 1 2-2h14a2 2 0 0 1 2 2z"></path>
                            </svg>
                        </div>
                        <div class="feature-box-content">
                            <h3>Full Reference List</h3>
                            <p>Each answer displays complete citations with title, authors, journal, year, and clickable PubMed ID for verification.</p>
                        </div>
                    </div>
                </div>
            </div>

            <!-- Example Search -->
            <div class="content-card">
                <h2>Example: Real-Time PubMed Query</h2>
                <p>
                    When you ask "What's the evidence for TXA in spine surgery?", here's what happens:
                </p>
                <div class="code-block">
                    <code>
1. Query Expansion: "TXA" → "tranexamic acid"<br>
2. PubMed Search Filters:<br>
   - Study types: Systematic Review, Meta-Analysis, RCT<br>
   - Field: anesthesiology, spine surgery<br>
   - Date range: 2015-2025<br>
3. Retrieve Top 20 Results<br>
4. Extract Metadata: Title, Authors, PMID, Year, Abstract<br>
5. GPT-4o Synthesis: Analyze papers, synthesize answer<br>
6. Output: Evidence-based answer + Full citations<br>
7. Confidence Badge: High (15 RCTs + 3 meta-analyses found)
                    </code>
                </div>
            </div>

            <!-- AI Self-Verification Protocol -->
            <div class="content-card">
                <h2>AI Self-Verification & Quality Assurance</h2>
                <p>
                    To minimize errors and ensure clinical accuracy, every AI-generated response undergoes a rigorous internal self-verification protocol before being presented. This multi-layered quality control system operates transparently within each answer generation process.
                </p>

                <div class="features-grid">
                    <div class="feature-box">
                        <div class="feature-icon-box green">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"></path>
                                <polyline points="14 2 14 8 20 8"></polyline>
                                <line x1="12" y1="18" x2="12" y2="12"></line>
                                <line x1="9" y1="15" x2="15" y2="15"></line>
                            </svg>
                        </div>
                        <div class="feature-box-content">
                            <h3>1. Dosing Accuracy Verification</h3>
                            <p>All drug doses are cross-checked against ASA guidelines, FDA package inserts, and major anesthesiology textbooks (Miller's, Barash, Stoelting's). Out-of-range doses trigger automatic correction.</p>
                        </div>
                    </div>

                    <div class="feature-box">
                        <div class="feature-icon-box blue">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <path d="M10.29 3.86 1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z"></path>
                                <line x1="12" y1="9" x2="12" y2="13"></line>
                                <line x1="12" y1="17" x2="12.01" y2="17"></line>
                            </svg>
                        </div>
                        <div class="feature-box-content">
                            <h3>2. Contraindication & Safety Check</h3>
                            <p>Before finalizing answers, the AI explicitly verifies that absolute contraindications are mentioned and critical safety warnings are not omitted. Missing warnings trigger revision.</p>
                        </div>
                    </div>

                    <div class="feature-box">
                        <div class="feature-icon-box purple">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <path d="M21 15a2 2 0 0 1-2 2H7l-4 4V5a2 2 0 0 1 2-2h14a2 2 0 0 1 2 2z"></path>
                                <line x1="9" y1="10" x2="15" y2="10"></line>
                                <line x1="9" y1="14" x2="15" y2="14"></line>
                            </svg>
                        </div>
                        <div class="feature-box-content">
                            <h3>3. Citation Verification Protocol</h3>
                            <p>Each cited paper's abstract is verified to ACTUALLY support the specific claim made. Citations that don't directly support statements are automatically removed to prevent misleading references.</p>
                        </div>
                    </div>

                    <div class="feature-box">
                        <div class="feature-icon-box amber">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <polyline points="20 6 9 17 4 12"></polyline>
                            </svg>
                        </div>
                        <div class="feature-box-content">
                            <h3>4. Guideline Consistency Check</h3>
                            <p>Answers are cross-referenced against current ASA, ESA, and specialty society guidelines. Conflicts with established practice guidelines trigger clarification or revision.</p>
                        </div>
                    </div>

                    <div class="feature-box">
                        <div class="feature-icon-box green">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <circle cx="12" cy="12" r="10"></circle>
                                <path d="M9.09 9a3 3 0 0 1 5.83 1c0 2-3 3-3 3"></path>
                                <line x1="12" y1="17" x2="12.01" y2="17"></line>
                            </svg>
                        </div>
                        <div class="feature-box-content">
                            <h3>5. Self-Questioning Protocol</h3>
                            <p>The AI asks itself "How do I know this is correct?" before each answer. Uncertain claims are either verified against provided papers or explicitly acknowledged as uncertain.</p>
                        </div>
                    </div>

                    <div class="feature-box">
                        <div class="feature-icon-box blue">
                            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <rect x="3" y="3" width="18" height="18" rx="2" ry="2"></rect>
                                <path d="m9 11 3 3L22 4"></path>
                            </svg>
                        </div>
                        <div class="feature-box-content">
                            <h3>6. Completeness Validation</h3>
                            <p>For dosing questions: verified inclusion of route, typical range, and maximum doses. For safety questions: confirmed coverage of both common and serious risks.</p>
                        </div>
                    </div>
                </div>

                <p style="margin-top: 1.5rem;">
                    <strong>Dynamic Precision Control:</strong> The AI's response temperature (creativity vs. precision) is automatically adjusted based on query type. Dosing questions use ultra-low temperature (0.05) for exact accuracy, while safety questions use 0.1 for factual rigor. This ensures clinical precision where it matters most.
                </p>
            </div>

            <!-- Limitations -->
            <div class="content-card">
                <h2>Limitations & Responsible Use</h2>
                <p>
                    While our methodology maximizes evidence quality, users should be aware of important limitations:
                </p>
                <p>
                    <strong>• Not a substitute for clinical judgment:</strong> Evidence-based answers inform, but never replace, individualized patient assessment and decision-making.
                </p>
                <p>
                    <strong>• Publication bias:</strong> Published literature may overrepresent positive findings and underrepresent negative or null results.
                </p>
                <p>
                    <strong>• AI synthesis limitations:</strong> GPT-4o may occasionally misinterpret nuanced findings or fail to capture context that experts would recognize.
                </p>
                <p>
                    <strong>• Rapidly evolving field:</strong> New evidence emerges constantly. Always verify critical decisions with the latest guidelines and institutional protocols.
                </p>
                <p>
                    <strong>• Educational tool only:</strong> GasConsult.ai is designed for learning and clinical decision support, not as a diagnostic or treatment tool.
                </p>
            </div>
        </main>

        <!-- Footer -->
        <footer class="footer">
            <div class="footer-inner">
                <span class="footer-text">© 2025 GasConsult.ai</span>
                <div class="footer-links">
                    <a href="/privacy" class="footer-link">Privacy</a>
                    <a href="/terms" class="footer-link">Terms</a>
                    <a href="mailto:contact@gasconsult.ai" class="footer-link">Contact</a>
                </div>
            </div>
        </footer>
    </div>

    <script>
        function toggleMobileMenu() {
            const menu = document.getElementById('mobileMenu');
            const btn = document.querySelector('.mobile-menu-btn');
            menu.classList.toggle('active');
            btn.classList.toggle('active');
        }
    </script>
</body>
</html>
"""

CRISIS_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <!-- Google Analytics 4 -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=G-01NZYD1DPP"></script>
    <script>
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-01NZYD1DPP');
    </script>

    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Crisis Protocols - gasconsult.ai</title>

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

        /* Skip to Content Link for Accessibility */
        .skip-to-content {
            position: absolute;
            top: -40px;
            left: 0;
            background: var(--blue-600);
            color: white;
            padding: 8px 16px;
            text-decoration: none;
            border-radius: 0 0 4px 0;
            z-index: 1000;
            font-weight: 600;
            transition: top 0.2s;
        }

        .skip-to-content:focus {
            top: 0;
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

        .nav-dropdown:has(.nav-dropdown-link.active) .nav-dropdown-toggle {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .nav-dropdown {
            position: relative;
            display: inline-block;
        }

        .nav-dropdown-toggle {
            cursor: pointer;
            background: none;
            border: none;
            font-family: inherit;
        }

        .nav-dropdown-menu {
            display: none;
            position: absolute;
            top: 100%;
            right: 0;
            background: white;
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            min-width: 200px;
            margin-top: 4px;
            z-index: 1000;
            overflow: hidden;
        }

        .nav-dropdown-menu.show {
            display: block;
        }

        .nav-dropdown-link {
            display: block;
            padding: 12px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .nav-dropdown-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
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

        .demo-note {
            background: linear-gradient(135deg, var(--blue-50) 0%, var(--purple-50) 100%);
            border: 2px solid var(--blue-200);
            border-radius: 16px;
            padding: 20px 24px;
            margin-bottom: 40px;
            text-align: center;
        }

        .demo-note-title {
            font-size: 16px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 8px;
        }

        .demo-note-text {
            font-size: 14px;
            color: var(--gray-600);
            line-height: 1.6;
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

        .category-icon svg {
            width: 22px;
            height: 22px;
            stroke: white;
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
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
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
            box-shadow: 0 12px 32px rgba(0, 0, 0, 0.1);
            transform: translateY(-4px);
        }

        .protocol-card.expanded {
            border-color: var(--blue-400);
            box-shadow: 0 16px 48px rgba(59, 130, 246, 0.2);
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

        .protocol-ref-count {
            display: inline-flex;
            align-items: center;
            gap: 4px;
            background: var(--blue-50);
            color: var(--blue-700);
            font-size: 11px;
            font-weight: 600;
            padding: 4px 8px;
            border-radius: 6px;
            margin-left: 8px;
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
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
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
            transition: max-height 0.5s cubic-bezier(0.4, 0, 0.2, 1);
        }

        .protocol-card.expanded .protocol-content {
            max-height: 6000px;
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

        .ref-num {
            display: inline-flex;
            align-items: center;
            justify-content: center;
            background: var(--blue-100);
            color: var(--blue-700);
            font-size: 10px;
            font-weight: 700;
            min-width: 18px;
            height: 18px;
            border-radius: 4px;
            margin-left: 2px;
            cursor: help;
            vertical-align: super;
            line-height: 1;
            padding: 2px 4px;
            transition: all 0.2s ease;
        }

        .ref-num:hover {
            background: var(--blue-600);
            color: white;
            transform: scale(1.15);
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

        .protocol-references {
            margin-top: 24px;
            padding-top: 16px;
            border-top: 1px solid var(--gray-200);
        }

        .references-toggle {
            display: flex;
            align-items: center;
            gap: 8px;
            background: var(--gray-50);
            padding: 10px 14px;
            border-radius: 10px;
            cursor: pointer;
            border: none;
            width: 100%;
            font-family: inherit;
            font-size: 13px;
            font-weight: 600;
            color: var(--gray-700);
            transition: all 0.2s ease;
        }

        .references-toggle:hover {
            background: var(--gray-100);
        }

        .references-toggle svg {
            width: 16px;
            height: 16px;
            transition: transform 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        }

        .references-toggle.expanded svg {
            transform: rotate(180deg);
        }

        .references-list {
            max-height: 0;
            overflow: hidden;
            transition: max-height 0.4s cubic-bezier(0.4, 0, 0.2, 1);
            padding-top: 0;
        }

        .references-list.show {
            max-height: 1200px;
            padding-top: 12px;
        }

        .reference-item {
            padding: 12px 14px;
            font-size: 13px;
            line-height: 1.6;
            color: var(--gray-700);
            border-left: 3px solid transparent;
            border-radius: 6px;
            margin-bottom: 4px;
            transition: all 0.2s ease;
        }

        .reference-item:hover {
            background: var(--gray-50);
            border-left-color: var(--blue-500);
        }

        .reference-num {
            display: inline-flex;
            align-items: center;
            justify-content: center;
            background: var(--blue-100);
            color: var(--blue-700);
            font-size: 11px;
            font-weight: 700;
            min-width: 22px;
            height: 22px;
            border-radius: 5px;
            margin-right: 10px;
        }

        .reference-citation {
            color: var(--gray-600);
        }

        .reference-citation strong {
            color: var(--gray-800);
        }

        .reference-citation a {
            color: var(--blue-600);
            text-decoration: none;
        }

        .reference-citation a:hover {
            text-decoration: underline;
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

        .footer-link:hover {
            color: var(--gray-700);
        }

        @media (min-width: 768px) {
            .nav { padding: 16px 32px; }
            .nav-inner { height: 64px; padding: 0 24px; border-radius: 20px; }
            .logo-icon svg { width: 42px; height: 15px; }
            .logo-text { font-size: 20px; }

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

            .footer { padding: 40px 32px; }
            .footer-inner { flex-direction: row; justify-content: space-between; text-align: left; }
            .footer-text { font-size: 14px; }
            .footer-links { gap: 32px; }
            .footer-link { font-size: 14px; }
        }

        @media (min-width: 1024px) {
            .nav { padding: 16px 40px; }

            .protocols-grid {
                grid-template-columns: repeat(3, 1fr);
            }

            .footer { padding: 48px 40px; }
        }

        .hidden {
            display: none !important;
        }

    </style>
</head>
<body>
    <a href="#main-content" class="skip-to-content">Skip to main content</a>
    <div class="bg-canvas">
        <div class="orb orb-1"></div>
        <div class="orb orb-2"></div>
        <div class="orb orb-3"></div>
    </div>
    <div class="grain"></div>

    <div class="page">
        <!-- Navigation -->
        <nav class="nav" role="navigation" aria-label="Main navigation">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo" aria-label="GasConsult.ai home">
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
                    <a href="/crisis" class="nav-link active">Crisis Protocols</a>
                    <div class="nav-dropdown">
                        <button class="nav-link nav-dropdown-toggle" onclick="toggleNavDropdown(event)">More ▼</button>
                        <div class="nav-dropdown-menu">
                            <a href="/hypotension" class="nav-dropdown-link">IOH Predictor</a>
                            <a href="/difficult-airway" class="nav-dropdown-link">Difficult Airway</a>
                            <a href="/informed-consent" class="nav-dropdown-link">Informed Consent</a>
                        </div>
                    </div>
                </div>

                <button class="mobile-menu-btn" onclick="toggleMobileMenu()">
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
            <a href="/difficult-airway" class="mobile-menu-link">Difficult Airway</a>
            <a href="/informed-consent" class="mobile-menu-link">Informed Consent</a>
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
                    <div class="category-icon red">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <path d="M22 12h-4l-3 9L9 3l-3 9H2"/>
                        </svg>
                    </div>
                    <h2 class="category-title">Life-Threatening / Cardiac</h2>
                </div>
                <div class="protocols-grid">

                    <!-- Malignant Hyperthermia (ENHANCED) -->
                    <div class="protocol-card red" data-keywords="malignant hyperthermia mh dantrolene hypermetabolic crisis muscle rigidity hyperthermia" onclick="toggleProtocol(this)">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">
                                    Malignant Hyperthermia
                                    <span class="protocol-ref-count"><svg style="width:14px;height:14px;margin-right:2px;" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path><path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path><path d="M8 7h8"></path><path d="M8 11h8"></path></svg> 5 refs</span>
                                </h3>
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
                                        <li class="protocol-step"><strong>STOP triggers immediately:</strong> Discontinue all volatile anesthetics and succinylcholine<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Call for help:</strong> Activate MH emergency protocol, assign roles<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Hyperventilate with 100% O₂:</strong> 2-3× normal minute ventilation to eliminate CO₂<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Give dantrolene immediately:</strong> See dosing below<sup class="ref-num">2</sup></li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Dantrolene Dosing<sup class="ref-num">2</sup></div>
                                    <div class="dose-detail"><strong>Initial:</strong> 2.5 mg/kg IV rapid push (reconstitute each 20mg vial with 60mL sterile water)</div>
                                    <div class="dose-detail"><strong>Repeat:</strong> 1 mg/kg boluses every 5-10 min until signs resolve</div>
                                    <div class="dose-detail"><strong>Maximum:</strong> Up to 10 mg/kg in acute phase (rarely >10 vials needed initially)</div>
                                    <div class="dose-detail"><strong>Continuation:</strong> 1 mg/kg IV q6h × 24-48h to prevent recrudescence</div>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Supportive Care</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Cool the patient:</strong> Cold IV saline, ice packs to groin/axilla, cooling blanket. Target temp <38.5°C<sup class="ref-num">1,3</sup></li>
                                        <li class="protocol-step"><strong>Treat hyperkalemia:</strong> Insulin/dextrose, calcium chloride, bicarbonate, avoid Ca²⁺ channel blockers with dantrolene<sup class="ref-num">3</sup></li>
                                        <li class="protocol-step"><strong>Manage arrhythmias:</strong> Standard ACLS (avoid calcium channel blockers)<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Monitor urine output:</strong> Foley catheter, maintain >1 mL/kg/h to prevent myoglobin-induced renal failure<sup class="ref-num">3</sup></li>
                                        <li class="protocol-step"><strong>Labs:</strong> ABG, electrolytes, CK, lactate, coags, myoglobin q6h<sup class="ref-num">1</sup></li>
                                    </ol>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">📞 MH Hotline</div>
                                    <div class="info-detail"><strong>USA:</strong> 1-800-MH-HYPER (1-800-644-9737)</div>
                                    <div class="info-detail"><strong>Outside USA:</strong> +1-315-464-7079</div>
                                    <div class="info-detail">Expert consultant available 24/7 for real-time guidance<sup class="ref-num">4</sup></div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Key Points</div>
                                    <div class="warning-detail">• Early signs: Unexplained ↑EtCO₂, masseter spasm, tachycardia, hypercarbia refractory to ↑ventilation<br>• Late signs: Fever, rigidity, rhabdomyolysis, hyperkalemia, acidosis<br>• Dantrolene can cause profound weakness - prepare for prolonged ventilation<br>• ICU monitoring × 24-48h minimum (recrudescence occurs in ~25%)<sup class="ref-num">5</sup></div>
                                </div>

                                <!-- Inline References Section -->
                                <div class="protocol-references">
                                    <button class="references-toggle" onclick="event.stopPropagation(); toggleReferences(this)">
                                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <polyline points="6 9 12 15 18 9"></polyline>
                                        </svg>
                                        <span>View References (5)</span>
                                    </button>
                                    <div class="references-list">
                                        <div class="reference-item">
                                            <span class="reference-num">1</span>
                                            <span class="reference-citation">Litman RS, Griggs SM, Dowling JJ, et al. <strong>Malignant Hyperthermia Susceptibility and Related Diseases.</strong> <em>Anesthesiology</em>. 2018;128(1):159-167. PMID: 29200006</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">2</span>
                                            <span class="reference-citation">MHAUS Emergency Therapy for MH (2022 Update). <a href="https://www.mhaus.org/healthcare-professionals/be-prepared/managing-a-crisis/" target="_blank">mhaus.org/managing-a-crisis</a></span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">3</span>
                                            <span class="reference-citation">Rosenberg H, Pollock N, Schiemann A, et al. <strong>Malignant hyperthermia: a review.</strong> <em>Orphanet J Rare Dis</em>. 2015;10:93. PMID: 26238698</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">4</span>
                                            <span class="reference-citation">Malignant Hyperthermia Association of the United States (MHAUS). <strong>24/7 Emergency Hotline.</strong> <a href="https://www.mhaus.org" target="_blank">mhaus.org</a></span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">5</span>
                                            <span class="reference-citation">Larach MG, Brandom BW, Allen GC, et al. <strong>Malignant Hyperthermia Deaths Related to Inadequate Temperature Monitoring.</strong> <em>Anesthesiology</em>. 2019;130(1):40-51. PMID: 30475256</span>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- Cardiac Arrest (ENHANCED) -->
                    <div class="protocol-card red" data-keywords="cardiac arrest code blue cpr acls asystole vfib pea pulseless" onclick="toggleProtocol(this)">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">
                                    Cardiac Arrest (ACLS)
                                    <span class="protocol-ref-count"><svg style="width:14px;height:14px;margin-right:2px;" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path><path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path><path d="M8 7h8"></path><path d="M8 11h8"></path></svg> 4 refs</span>
                                </h3>
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
                                        <li class="protocol-step"><strong>Start CPR immediately:</strong> 100-120 compressions/min, depth 2-2.4 inches, minimize interruptions<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Call for help / Code Blue</strong></li>
                                        <li class="protocol-step"><strong>Attach defibrillator/monitor:</strong> Identify rhythm</li>
                                        <li class="protocol-step"><strong>Secure airway:</strong> ETT or supraglottic device + capnography (target EtCO₂ >10 mmHg)<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>IV/IO access:</strong> Establish vascular access for medications</li>
                                        <li class="protocol-step"><strong>Consider reversible causes (H's and T's)</strong><sup class="ref-num">2</sup></li>
                                    </ol>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Shockable Rhythms (VF/pVT)</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Defibrillate:</strong> Biphasic 120-200J (or manufacturer recommendation), resume CPR immediately × 2 min<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>After 2nd shock:</strong> Epinephrine 1 mg IV/IO q3-5min<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>After 3rd shock:</strong> Amiodarone 300 mg IV/IO (or lidocaine 1-1.5 mg/kg if amio unavailable)<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Continue CPR + defibrillation every 2 min</strong></li>
                                    </ol>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Non-Shockable Rhythms (PEA/Asystole)</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>High-quality CPR × 2 min</strong><sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Epinephrine 1 mg IV/IO immediately,</strong> then q3-5min<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Consider atropine 1 mg IV</strong> if slow PEA (rate <60), may repeat to total 3 mg</li>
                                        <li class="protocol-step"><strong>Treat reversible causes aggressively</strong> (see H's and T's below)<sup class="ref-num">2</sup></li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Key Medications<sup class="ref-num">1</sup></div>
                                    <div class="dose-detail"><strong>Epinephrine:</strong> 1 mg (1:10,000) IV/IO every 3-5 minutes</div>
                                    <div class="dose-detail"><strong>Amiodarone:</strong> 300 mg IV/IO first dose, then 150 mg second dose</div>
                                    <div class="dose-detail"><strong>Lidocaine (alternative):</strong> 1-1.5 mg/kg first dose, then 0.5-0.75 mg/kg</div>
                                    <div class="dose-detail"><strong>Sodium bicarbonate:</strong> 1 mEq/kg (for hyperkalemia, TCA overdose, prolonged arrest)</div>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🔍 H's and T's (Reversible Causes)<sup class="ref-num">2</sup></div>
                                    <div class="info-detail"><strong>H's:</strong> Hypovolemia, Hypoxia, H⁺ (acidosis), Hyper/hypokalemia, Hypothermia<br>
                                    <strong>T's:</strong> Tension pneumothorax, Tamponade (cardiac), Toxins, Thrombosis (coronary/pulmonary)</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Anesthesia-Specific Considerations</div>
                                    <div class="warning-detail">• Turn off volatile anesthetics during arrest<br>• Consider <strong>anesthesia-specific causes:</strong> local anesthetic toxicity (give lipid emulsion), hyperkalemia from succinylcholine, pneumothorax from line placement, total spinal<sup class="ref-num">3</sup><br>• Continue CPR during transfer to ICU if needed<br>• Document ROSC time, rhythm changes, total epi/defib doses<sup class="ref-num">4</sup></div>
                                </div>

                                <!-- Inline References Section -->
                                <div class="protocol-references">
                                    <button class="references-toggle" onclick="event.stopPropagation(); toggleReferences(this)">
                                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <polyline points="6 9 12 15 18 9"></polyline>
                                        </svg>
                                        <span>View References (4)</span>
                                    </button>
                                    <div class="references-list">
                                        <div class="reference-item">
                                            <span class="reference-num">1</span>
                                            <span class="reference-citation">Panchal AR, et al. <strong>Part 3: Adult Basic and Advanced Life Support. 2020 AHA Guidelines.</strong> <em>Circulation</em>. 2020;142(16_suppl_2):S366-S468. PMID: 33081529</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">2</span>
                                            <span class="reference-citation">Link MS, et al. <strong>Part 7: Adult Advanced Cardiovascular Life Support.</strong> <em>Circulation</em>. 2015;132(18 Suppl 2):S444-64. PMID: 26472995</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">3</span>
                                            <span class="reference-citation">Neal JM, et al. <strong>Anesthesia-Related Cardiac Arrest: A Registry Analysis.</strong> <em>Anesthesiology</em>. 2014;120(4):829-838. PMID: 24694841</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">4</span>
                                            <span class="reference-citation">Soar J, et al. <strong>European Resuscitation Council Guidelines.</strong> <em>Resuscitation</em>. 2021;161:98-114. PMID: 33773831</span>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- Anaphylaxis (ENHANCED) -->
                    <div class="protocol-card red" data-keywords="anaphylaxis allergic reaction epinephrine bronchospasm hypotension urticaria angioedema" onclick="toggleProtocol(this)">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">
                                    Anaphylaxis
                                    <span class="protocol-ref-count"><svg style="width:14px;height:14px;margin-right:2px;" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path><path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path><path d="M8 7h8"></path><path d="M8 11h8"></path></svg> 4 refs</span>
                                </h3>
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
                                        <li class="protocol-step"><strong>STOP suspected trigger:</strong> Antibiotics, NMBs, latex, colloids are most common<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Call for help</strong></li>
                                        <li class="protocol-step"><strong>Epinephrine IM immediately:</strong> 0.3-0.5 mg (0.3-0.5 mL of 1:1000) into anterolateral thigh, repeat q5-15min PRN<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>100% O₂:</strong> Maintain airway, consider early intubation if upper airway edema<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Aggressive fluid resuscitation:</strong> 20-50 mL/kg crystalloid rapidly for refractory hypotension<sup class="ref-num">2</sup></li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Epinephrine Dosing<sup class="ref-num">2</sup></div>
                                    <div class="dose-detail"><strong>IM (first-line):</strong> 0.3-0.5 mg (1:1000) into thigh, repeat q5-15min</div>
                                    <div class="dose-detail"><strong>IV bolus (severe/arrest):</strong> 10-100 mcg (0.01-0.1 mg) slow push, titrate to effect</div>
                                    <div class="dose-detail"><strong>IV infusion (refractory):</strong> 0.05-0.5 mcg/kg/min, titrate to BP/HR</div>
                                    <div class="dose-detail"><strong>Pediatric IM:</strong> 0.01 mg/kg (max 0.5 mg)</div>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Adjunct Therapies</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>H1 blocker:</strong> Diphenhydramine 25-50 mg IV slowly<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>H2 blocker:</strong> Famotidine 20 mg IV or ranitidine 50 mg IV<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Corticosteroids:</strong> Methylprednisolone 1-2 mg/kg IV (prevents late-phase reaction)<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Bronchodilators:</strong> Albuterol for persistent bronchospasm</li>
                                        <li class="protocol-step"><strong>Glucagon (if on β-blockers):</strong> 1-2 mg IV (epinephrine may be ineffective)<sup class="ref-num">2</sup></li>
                                    </ol>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🩺 Diagnostic Criteria (2 or more)<sup class="ref-num">2</sup></div>
                                    <div class="info-detail">• <strong>Skin/mucosal:</strong> Urticaria, angioedema, flushing<br>• <strong>Respiratory:</strong> Bronchospasm, wheezing, stridor, dyspnea, ↓SpO₂<br>• <strong>Cardiovascular:</strong> Hypotension (SBP <90 or >30% drop), tachycardia, arrhythmia, collapse<br>• <strong>GI:</strong> Cramping, vomiting, diarrhea</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Key Considerations</div>
                                    <div class="warning-detail">• <strong>Common triggers:</strong> NMBs (rocuronium, succinylcholine), antibiotics (cephalosporins, penicillins), latex, chlorhexidine<sup class="ref-num">3</sup><br>• Confirm diagnosis: Send tryptase levels (draw immediately, then 1-2h and 24h later)<sup class="ref-num">4</sup><br>• Biphasic reactions occur in 20% - observe ≥4-6h minimum, admit if severe<sup class="ref-num">2</sup><br>• Document reaction in chart and advise patient to see allergist<br>• Refractory hypotension: Consider methylene blue 1-2 mg/kg for vasoplegia</div>
                                </div>

                                <!-- Inline References Section -->
                                <div class="protocol-references">
                                    <button class="references-toggle" onclick="event.stopPropagation(); toggleReferences(this)">
                                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <polyline points="6 9 12 15 18 9"></polyline>
                                        </svg>
                                        <span>View References (4)</span>
                                    </button>
                                    <div class="references-list">
                                        <div class="reference-item">
                                            <span class="reference-num">1</span>
                                            <span class="reference-citation">Dewachter P, et al. <strong>Perioperative Anaphylaxis.</strong> <em>Anesthesiology</em>. 2009;111(5):1141-1150. PMID: 19858877</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">2</span>
                                            <span class="reference-citation">Shaker MS, et al. <strong>Anaphylaxis: A 2020 Practice Parameter Update.</strong> <em>Ann Allergy Asthma Immunol</em>. 2020;125(4):346-373. PMID: 32846301</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">3</span>
                                            <span class="reference-citation">Mertes PM, et al. <strong>Reducing the Risk of Anaphylaxis During Anesthesia.</strong> <em>J Allergy Clin Immunol Pract</em>. 2020;8(8):2544-2555. PMID: 32505781</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">4</span>
                                            <span class="reference-citation">Fisher MM, Baldo BA. <strong>Mast Cell Tryptase in Anaesthetic Anaphylactoid Reactions.</strong> <em>Br J Anaesth</em>. 1998;80(1):26-29. PMID: 9505773</span>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- LAST (ENHANCED) -->
                    <div class="protocol-card red" data-keywords="last local anesthetic systemic toxicity lipid emulsion intralipid bupivacaine ropivacaine seizure arrhythmia" onclick="toggleProtocol(this)">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">
                                    Local Anesthetic Systemic Toxicity (LAST)
                                    <span class="protocol-ref-count"><svg style="width:14px;height:14px;margin-right:2px;" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path><path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path><path d="M8 7h8"></path><path d="M8 11h8"></path></svg> 3 refs</span>
                                </h3>
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
                                        <li class="protocol-step"><strong>Call for help:</strong> Get lipid emulsion (Intralipid 20%)<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Airway management:</strong> 100% O₂, ventilate if needed, suppress seizures</li>
                                        <li class="protocol-step"><strong>Give lipid emulsion immediately</strong> (see dosing below) - DO NOT DELAY<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>If cardiac arrest:</strong> Start CPR, consider prolonged resuscitation (LAST arrest can require >1h CPR)<sup class="ref-num">2</sup></li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Lipid Emulsion 20% (Intralipid) Dosing<sup class="ref-num">1</sup></div>
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

                                <div class="info-box">
                                    <div class="info-box-title">🩺 Signs of LAST</div>
                                    <div class="info-detail"><strong>Early CNS:</strong> Circumoral numbness, metallic taste, tinnitus, confusion, agitation<br>
                                    <strong>Severe CNS:</strong> Seizures, loss of consciousness<br>
                                    <strong>Cardiac:</strong> Bradycardia, hypotension, arrhythmias (wide QRS), asystole, PEA<sup class="ref-num">1</sup><br>
                                    <strong>Note:</strong> Cardiac toxicity can occur WITHOUT preceding CNS symptoms (especially bupivacaine)</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Critical Points</div>
                                    <div class="warning-detail">• Lipid emulsion is PRIMARY treatment - give early, don't wait for arrest<sup class="ref-num">1</sup><br>• Bupivacaine/ropivacaine are more cardiotoxic than lidocaine/mepivacaine<br>• Max doses: Bupivacaine 2.5 mg/kg plain, 3 mg/kg with epi; Lidocaine 5 mg/kg plain, 7 mg/kg with epi<sup class="ref-num">3</sup><br>• Post-resuscitation: Monitor ≥4-6h (cardiac arrest patients → ICU), watch for pancreatitis from lipid load</div>
                                </div>

                                <!-- Inline References Section -->
                                <div class="protocol-references">
                                    <button class="references-toggle" onclick="event.stopPropagation(); toggleReferences(this)">
                                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <polyline points="6 9 12 15 18 9"></polyline>
                                        </svg>
                                        <span>View References (3)</span>
                                    </button>
                                    <div class="references-list">
                                        <div class="reference-item">
                                            <span class="reference-num">1</span>
                                            <span class="reference-citation">Neal JM, et al. <strong>ASRA Practice Advisory on Local Anesthetic Systemic Toxicity.</strong> <em>Reg Anesth Pain Med</em>. 2018;43(2):113-123. PMID: 29356773</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">2</span>
                                            <span class="reference-citation">Gitman M, Barrington MJ. <strong>Local Anesthetic Systemic Toxicity: A Review of Recent Case Reports.</strong> <em>Reg Anesth Pain Med</em>. 2018;43(2):124-130. PMID: 29095244</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">3</span>
                                            <span class="reference-citation">American Society of Regional Anesthesia and Pain Medicine. <strong>Checklist for Treatment of Local Anesthetic Systemic Toxicity.</strong> <em>Reg Anesth Pain Med</em>. 2012;37(1):16-18</span>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- High/Total Spinal (ENHANCED) -->
                    <div class="protocol-card red" data-keywords="high spinal total spinal epidural anesthesia neuraxial hypotension bradycardia respiratory arrest paralysis" onclick="toggleProtocol(this)">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">
                                    High/Total Spinal Anesthesia
                                    <span class="protocol-ref-count"><svg style="width:14px;height:14px;margin-right:2px;" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path><path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path><path d="M8 7h8"></path><path d="M8 11h8"></path></svg> 3 refs</span>
                                </h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Excessive cephalad spread of neuraxial anesthetic causing cardiorespiratory compromise</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag immediate">Immediate</span>
                            <span class="protocol-tag call-help">Call Help</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Actions - ABC Approach</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Airway:</strong> Secure airway immediately if respiratory distress. Intubate if apneic or unable to protect airway<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Breathing:</strong> Positive pressure ventilation with 100% O₂ (bag-mask or intubation)<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Circulation:</strong> Treat hypotension and bradycardia (see below)<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Position:</strong> Supine or slight Trendelenburg (controversy: may worsen vs improve symptoms)<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Reassure patient:</strong> If awake, explain that this is temporary and reversible</li>
                                    </ol>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Hemodynamic Management</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Hypotension:</strong> Phenylephrine 50-200 mcg IV boluses or ephedrine 5-10 mg IV boluses<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Severe hypotension:</strong> Epinephrine 10-100 mcg IV boluses (or start infusion)<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Bradycardia:</strong> Atropine 0.4-1 mg IV or glycopyrrolate 0.2-0.4 mg IV<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Severe bradycardia/arrest:</strong> Epinephrine 1 mg IV, start CPR if pulseless<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Volume:</strong> Rapid IV fluid bolus 500-1000 mL crystalloid<sup class="ref-num">1</sup></li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Key Medications<sup class="ref-num">2</sup></div>
                                    <div class="dose-detail"><strong>Phenylephrine:</strong> 50-200 mcg IV boluses (or 0.5-1 mcg/kg/min infusion)</div>
                                    <div class="dose-detail"><strong>Ephedrine:</strong> 5-10 mg IV boluses</div>
                                    <div class="dose-detail"><strong>Epinephrine (severe):</strong> 10-100 mcg IV boluses or 0.01-0.1 mcg/kg/min infusion</div>
                                    <div class="dose-detail"><strong>Atropine:</strong> 0.4-1 mg IV (for bradycardia)</div>
                                    <div class="dose-detail"><strong>Glycopyrrolate:</strong> 0.2-0.4 mg IV (alternative for bradycardia)</div>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🔍 Clinical Features<sup class="ref-num">1</sup></div>
                                    <div class="info-detail"><strong>Sensory level:</strong> T1-T4 (high spinal) vs C3-C8 (total spinal)<br>
                                    <strong>Motor:</strong> Upper extremity weakness, difficulty breathing, inability to speak<br>
                                    <strong>Cardiovascular:</strong> Hypotension, bradycardia (Bezold-Jarisch reflex)<br>
                                    <strong>Respiratory:</strong> Dyspnea, respiratory arrest (phrenic nerve paralysis C3-C5)<br>
                                    <strong>Neurologic:</strong> Unconsciousness (if brainstem affected), dilated pupils, nystagmus</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Key Points</div>
                                    <div class="warning-detail">• <strong>Self-limited:</strong> Block will spontaneously regress over 1-3 hours<sup class="ref-num">3</sup><br>• <strong>Early intubation:</strong> Don't wait for complete respiratory arrest - intubate early if distress<br>• <strong>Epinephrine early:</strong> Don't hesitate to use epi if severe hypotension/bradycardia<br>• <strong>Avoid sedation:</strong> Patient may already be unconscious from brainstem anesthesia<br>• <strong>Document sensory level:</strong> Check bilateral sensation to assess block height</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Prevention & Risk Factors</div>
                                    <div class="warning-detail"><strong>Risk factors for high spinal:</strong><br>• Excessive local anesthetic dose or volume<br>• Rapid injection or patient positioning (head-down)<br>• Inadvertent subdural or subarachnoid injection during epidural<br>• Short patient height, pregnancy (reduced CSF volume)<br>• Barbotage technique<sup class="ref-num">1</sup><br><br><strong>Prevention:</strong> Use appropriate dose for patient height, test dose, incremental dosing for epidurals, careful patient positioning</div>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🔍 Differential Diagnosis</div>
                                    <div class="info-detail">• <strong>Local anesthetic systemic toxicity (LAST):</strong> CNS excitation (seizures) before cardiac arrest<br>• <strong>Vasovagal syncope:</strong> Bradycardia + hypotension but normal respirations<br>• <strong>Anaphylaxis:</strong> Bronchospasm, urticaria, angioedema<br>• <strong>Myocardial infarction:</strong> ECG changes, chest pain<br>• <strong>Pulmonary embolism:</strong> Sudden hypoxia, tachycardia</div>
                                </div>

                                <!-- Inline References Section -->
                                <div class="protocol-references">
                                    <button class="references-toggle" onclick="event.stopPropagation(); toggleReferences(this)">
                                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <polyline points="6 9 12 15 18 9"></polyline>
                                        </svg>
                                        <span>View References (3)</span>
                                    </button>
                                    <div class="references-list">
                                        <div class="reference-item">
                                            <span class="reference-num">1</span>
                                            <span class="reference-citation">Reina MA, et al. <strong>Clinical implications of epidural fat in the spinal canal: a scanning electron microscopic study.</strong> <em>Acta Anaesthesiol Scand</em>. 2009;53(5):641-647. PMID: 19419359</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">2</span>
                                            <span class="reference-citation">Pollard JB. <strong>Cardiac arrest during spinal anesthesia: common mechanisms and strategies for prevention.</strong> <em>Anesth Analg</em>. 2001;92(1):252-256. PMID: 11133637</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">3</span>
                                            <span class="reference-citation">Auroux P, et al. <strong>Total spinal anesthesia after epidural test dose.</strong> <em>Anesthesiology</em>. 2000;92(5):1514-1516. PMID: 10781306</span>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- VAE - Venous Air Embolism (ENHANCED) -->
                    <div class="protocol-card red" data-keywords="vae venous air embolism gas embolism neurosurgery sitting position hypotension cardiovascular collapse mill wheel murmur" onclick="toggleProtocol(this)">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">
                                    Venous Air Embolism (VAE)
                                    <span class="protocol-ref-count"><svg style="width:14px;height:14px;margin-right:2px;" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path><path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path><path d="M8 7h8"></path><path d="M8 11h8"></path></svg> 4 refs</span>
                                </h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Air entrainment into venous system causing cardiovascular and/or neurologic compromise</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag immediate">Immediate</span>
                            <span class="protocol-tag call-help">Call Help</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Actions</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Notify surgeon:</strong> STOP surgery, flood field with saline, apply bone wax to exposed bone<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>100% O₂:</strong> Discontinue N₂O immediately (if used), switch to FiO₂ 1.0<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Position:</strong> Lower surgical site below heart level if possible (reverse Trendelenburg → supine or head down)<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Attempt aspiration:</strong> If central line in place, aspirate from distal port to remove air<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Support hemodynamics:</strong> IV fluids, vasopressors/inotropes as needed (see below)<sup class="ref-num">2</sup></li>
                                    </ol>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Advanced Management</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Hemodynamic support:</strong> Phenylephrine, ephedrine, or epinephrine as needed for hypotension<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>CPR if cardiac arrest:</strong> Chest compressions may help break up air lock and disperse air<sup class="ref-num">3</sup></li>
                                        <li class="protocol-step"><strong>Consider hyperbaric oxygen:</strong> If neurologic deficits persist (paradoxical embolism)<sup class="ref-num">4</sup></li>
                                        <li class="protocol-step"><strong>Durant maneuver (controversial):</strong> Left lateral decubitus + Trendelenburg to trap air in right atrium<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>TEE/TTE:</strong> If available, can visualize air in heart chambers and assess severity<sup class="ref-num">1</sup></li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Hemodynamic Support<sup class="ref-num">2</sup></div>
                                    <div class="dose-detail"><strong>Phenylephrine:</strong> 50-200 mcg IV boluses or infusion 0.5-1 mcg/kg/min</div>
                                    <div class="dose-detail"><strong>Ephedrine:</strong> 5-10 mg IV boluses</div>
                                    <div class="dose-detail"><strong>Epinephrine:</strong> 10-100 mcg IV boluses (or 0.01-0.1 mcg/kg/min infusion) for severe hypotension</div>
                                    <div class="dose-detail"><strong>Atropine:</strong> 0.4-1 mg IV if bradycardic</div>
                                    <div class="dose-detail"><strong>IV fluids:</strong> Rapid bolus 500-1000 mL crystalloid</div>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🔍 Clinical Signs (by Monitor Sensitivity)<sup class="ref-num">1</sup></div>
                                    <div class="info-detail"><strong>Most sensitive (earliest):</strong><br>• Precordial Doppler: "Mill wheel" murmur (detects 0.25 mL air)<br>• TEE: Visualize air in heart chambers<br>• Sudden ↓EtCO₂ (dead space ventilation)<br><br><strong>Moderate sensitivity:</strong><br>• Hypotension, tachycardia<br>• Hypoxia (↓SpO₂)<br>• Arrhythmias<br><br><strong>Late/severe:</strong><br>• Cardiovascular collapse<br>• Cardiac arrest (air lock in RV outflow)</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ High-Risk Procedures & Prevention</div>
                                    <div class="warning-detail"><strong>High-risk surgeries:</strong><br>• Sitting position craniotomy (most common)<br>• Posterior fossa surgery<br>• Neurosurgery with head elevated >15°<br>• Spine surgery (especially cervical)<br>• Laparoscopy, hepatic resection, cesarean section<sup class="ref-num">1</sup><br><br><strong>Prevention:</strong> Avoid N₂O, optimize patient positioning (minimize head elevation), adequate hydration, ensure good communication with surgeon, consider central venous catheter for aspiration, use precordial Doppler or TEE monitoring</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Paradoxical Air Embolism (PAE)</div>
                                    <div class="warning-detail">• Occurs when air crosses from right → left heart via PFO (present in ~25% of population)<sup class="ref-num">4</sup><br>• Results in <strong>stroke, MI, or organ ischemia</strong><br>• Higher risk with sitting position (negative intrathoracic pressure)<br>• Signs: Sudden neurologic deficit, ST-segment changes on ECG<br>• Management: Immediate 100% O₂, hemodynamic support, consider hyperbaric oxygen therapy<br>• Screen high-risk patients for PFO with TEE or bubble study (controversial)</div>
                                </div>

                                <!-- Inline References Section -->
                                <div class="protocol-references">
                                    <button class="references-toggle" onclick="event.stopPropagation(); toggleReferences(this)">
                                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <polyline points="6 9 12 15 18 9"></polyline>
                                        </svg>
                                        <span>View References (4)</span>
                                    </button>
                                    <div class="references-list">
                                        <div class="reference-item">
                                            <span class="reference-num">1</span>
                                            <span class="reference-citation">Mirski MA, et al. <strong>Diagnosis and treatment of vascular air embolism.</strong> <em>Anesthesiology</em>. 2007;106(1):164-177. PMID: 17197859</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">2</span>
                                            <span class="reference-citation">Shaikh N, Ummunisa F. <strong>Acute management of vascular air embolism.</strong> <em>J Emerg Trauma Shock</em>. 2009;2(3):180-185. PMID: 20009308</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">3</span>
                                            <span class="reference-citation">Vesely TM, et al. <strong>Air embolism during insertion of central venous catheters.</strong> <em>J Vasc Interv Radiol</em>. 2001;12(11):1291-1295. PMID: 11698626</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">4</span>
                                            <span class="reference-citation">Muth CM, Shank ES. <strong>Gas embolism.</strong> <em>N Engl J Med</em>. 2000;342(7):476-482. PMID: 10675429</span>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- Hyperkalemia (ENHANCED) -->
                    <div class="protocol-card red" data-keywords="hyperkalemia potassium cardiac arrest peaked t waves arrhythmia succinylcholine renal failure" onclick="toggleProtocol(this)">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">
                                    Hyperkalemia
                                    <span class="protocol-ref-count"><svg style="width:14px;height:14px;margin-right:2px;" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path><path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path><path d="M8 7h8"></path><path d="M8 11h8"></path></svg> 3 refs</span>
                                </h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Life-threatening elevation in serum potassium causing cardiac arrhythmias and arrest</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag immediate">Immediate</span>
                            <span class="protocol-tag call-help">Call Help</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Actions (3-Step Approach)</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>STEP 1 - Membrane stabilization (fastest):</strong> Calcium chloride 10% 10-20 mL IV over 2-5 min OR calcium gluconate 10% 30-60 mL IV<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>STEP 2 - Shift K⁺ intracellularly:</strong> Insulin 10 units IV + dextrose 25 g (D50W 50 mL) IV push<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>STEP 3 - Remove K⁺ from body:</strong> Diuretics (furosemide 40-80 mg IV) if renal function intact<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Additional shift therapy:</strong> Albuterol 10-20 mg nebulized (lowers K⁺ by 0.5-1 mEq/L)<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>If cardiac arrest:</strong> Start CPR, repeat calcium, consider emergency dialysis<sup class="ref-num">3</sup></li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Treatment Protocol (K⁺ >6.5 mEq/L or ECG changes)<sup class="ref-num">1</sup></div>
                                    <div class="dose-detail"><strong>1. Calcium chloride 10%:</strong> 10-20 mL (1-2 g) IV over 2-5 min (onset 1-3 min, duration 30-60 min) - OR -</div>
                                    <div class="dose-detail"><strong>1. Calcium gluconate 10%:</strong> 30-60 mL (3-6 g) IV over 2-5 min (less tissue necrosis if extravasates)</div>
                                    <div class="dose-detail"><strong>2. Regular insulin:</strong> 10 units IV + <strong>D50W 50 mL</strong> (25 g dextrose) IV push (onset 15-30 min, duration 4-6 h)</div>
                                    <div class="dose-detail"><strong>3. Albuterol:</strong> 10-20 mg (2-4 mL of 0.5% solution) nebulized over 10 min (onset 30 min, lowers K⁺ 0.5-1 mEq/L)</div>
                                    <div class="dose-detail"><strong>4. Sodium bicarbonate:</strong> 50-100 mEq IV (controversial, mainly for acidosis)</div>
                                    <div class="dose-detail"><strong>5. Furosemide:</strong> 40-80 mg IV (if renal function intact)</div>
                                    <div class="dose-detail"><strong>6. Kayexalate/Patiromer:</strong> NOT for acute management (slow onset, hours to days)</div>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🔍 ECG Changes by Severity<sup class="ref-num">1</sup></div>
                                    <div class="info-detail"><strong>Mild (K⁺ 5.5-6.5):</strong> Peaked, narrow T waves<br>
                                    <strong>Moderate (K⁺ 6.5-8.0):</strong> PR prolongation, P wave flattening/loss, QRS widening<br>
                                    <strong>Severe (K⁺ >8.0):</strong> Sine wave pattern, ventricular fibrillation, asystole</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Anesthesia-Specific Causes</div>
                                    <div class="warning-detail">• <strong>Succinylcholine:</strong> Especially in burns, crush injuries, denervation, prolonged immobility, neuromuscular disease<sup class="ref-num">2</sup><br>• <strong>Massive transfusion:</strong> Stored blood has high K⁺ (up to 50 mEq/L in old units)<br>• <strong>Tourniquet release:</strong> Sudden K⁺ release from ischemic limb<br>• <strong>Tumor lysis syndrome:</strong> Chemotherapy, large tumor burden<br>• <strong>Medications:</strong> ACE-I, ARBs, K⁺-sparing diuretics, NSAIDs<br>• <strong>Renal failure:</strong> Most common chronic cause</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Key Points</div>
                                    <div class="warning-detail">• <strong>Calcium first:</strong> Most important immediate treatment - stabilizes cardiac membrane<br>• <strong>Don't mix calcium with bicarb:</strong> Forms precipitate (give via separate IV)<br>• <strong>Monitor glucose:</strong> After insulin/dextrose therapy (risk of hypoglycemia 4-6h later)<br>• <strong>Repeat labs:</strong> Check K⁺ q2-4h until normalized<br>• <strong>Emergency dialysis:</strong> Consider if K⁺ >7.5 mEq/L, refractory, or cardiac arrest<sup class="ref-num">3</sup><br>• <strong>Avoid succinylcholine:</strong> If hyperkalemia suspected or patient at risk</div>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🔍 Pseudohyperkalemia (Lab Error)</div>
                                    <div class="info-detail">• <strong>Hemolysis:</strong> Most common cause of falsely elevated K⁺ (check sample for pink/red tint)<br>• <strong>Fist clenching:</strong> During phlebotomy<br>• <strong>Prolonged tourniquet time</strong><br>• <strong>Thrombocytosis or leukocytosis:</strong> Cell lysis in sample<br>• If suspected: Repeat with non-hemolyzed sample, correlate with ECG changes</div>
                                </div>

                                <!-- Inline References Section -->
                                <div class="protocol-references">
                                    <button class="references-toggle" onclick="event.stopPropagation(); toggleReferences(this)">
                                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <polyline points="6 9 12 15 18 9"></polyline>
                                        </svg>
                                        <span>View References (3)</span>
                                    </button>
                                    <div class="references-list">
                                        <div class="reference-item">
                                            <span class="reference-num">1</span>
                                            <span class="reference-citation">Hollander-Rodriguez JC, Calvert JF Jr. <strong>Hyperkalemia.</strong> <em>Am Fam Physician</em>. 2006;73(2):283-290. PMID: 16445274</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">2</span>
                                            <span class="reference-citation">Thapa S, Brull SJ. <strong>Succinylcholine-Induced Hyperkalemia in Patients with Renal Failure.</strong> <em>Anesth Analg</em>. 2000;91(5):1140-1145. PMID: 11049897</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">3</span>
                                            <span class="reference-citation">Alfonzo AV, et al. <strong>Potassium disorders-clinical spectrum and emergency management.</strong> <em>Resuscitation</em>. 2006;70(1):10-25. PMID: 16600469</span>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- PRIS - Propofol Infusion Syndrome (ENHANCED) -->
                    <div class="protocol-card red" data-keywords="pris propofol infusion syndrome rhabdomyolysis metabolic acidosis cardiac failure propofol icu sedation" onclick="toggleProtocol(this)">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">
                                    PRIS (Propofol Infusion Syndrome)
                                    <span class="protocol-ref-count"><svg style="width:14px;height:14px;margin-right:2px;" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path><path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path><path d="M8 7h8"></path><path d="M8 11h8"></path></svg> 4 refs</span>
                                </h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Rare but life-threatening syndrome from prolonged high-dose propofol infusion</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag immediate">Immediate</span>
                            <span class="protocol-tag call-help">Call Help</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Actions</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>STOP propofol immediately:</strong> Discontinue infusion - this is the MOST important intervention<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Switch sedation:</strong> Use alternative sedative (midazolam, dexmedetomidine, or ketamine)<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Supportive care:</strong> Hemodynamic support with vasopressors/inotropes as needed<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Labs:</strong> ABG, lactate, CK, troponin, lipase, electrolytes (especially K⁺), triglycerides, liver function<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Consider hemodialysis or ECMO:</strong> If refractory metabolic acidosis or cardiac failure<sup class="ref-num">3</sup></li>
                                    </ol>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Specific Management</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Metabolic acidosis:</strong> Bicarbonate therapy (controversial), treat underlying cause<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Hyperkalemia:</strong> Standard hyperkalemia protocol (calcium, insulin/dextrose, albuterol)<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Cardiac dysfunction:</strong> Inotropes (dobutamine, milrinone), avoid further myocardial stress<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Rhabdomyolysis:</strong> Aggressive fluid resuscitation, maintain urine output >1 mL/kg/h<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Renal replacement therapy:</strong> If severe metabolic derangement or renal failure<sup class="ref-num">3</sup></li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Alternative Sedation Options<sup class="ref-num">1</sup></div>
                                    <div class="dose-detail"><strong>Midazolam:</strong> Load 0.05-0.1 mg/kg IV, then 0.05-0.2 mg/kg/h infusion</div>
                                    <div class="dose-detail"><strong>Dexmedetomidine:</strong> Load 1 mcg/kg over 10 min, then 0.2-0.7 mcg/kg/h infusion</div>
                                    <div class="dose-detail"><strong>Ketamine:</strong> 0.5-1 mg/kg bolus, then 0.5-1 mg/kg/h infusion</div>
                                    <div class="dose-detail"><strong>Volatile anesthetics:</strong> Consider if in OR setting (sevoflurane, isoflurane)</div>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🔍 Classic Pentad of PRIS<sup class="ref-num">1</sup></div>
                                    <div class="info-detail">1. <strong>Metabolic acidosis:</strong> Severe, refractory (pH <7.2, lactate >5 mmol/L)<br>
                                    2. <strong>Rhabdomyolysis:</strong> Elevated CK (often >5000 U/L), myoglobinuria<br>
                                    3. <strong>Cardiac failure:</strong> Bradycardia → heart block → asystole, dilated cardiomyopathy<br>
                                    4. <strong>Renal failure:</strong> Acute kidney injury, often requiring dialysis<br>
                                    5. <strong>Hypertriglyceridemia/hepatomegaly:</strong> Lipemic serum, fatty liver</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Risk Factors<sup class="ref-num">1,4</sup></div>
                                    <div class="warning-detail"><strong>Classic triad:</strong><br>• <strong>High dose:</strong> >4 mg/kg/h (67 mcg/kg/min) for >48-72 hours<br>• <strong>Young age:</strong> Children > adults (but can occur in adults)<br>• <strong>Critical illness:</strong> Sepsis, traumatic brain injury, status epilepticus<br><br><strong>Other risk factors:</strong><br>• Catecholamine or steroid co-administration<br>• Inadequate carbohydrate intake<br>• Mitochondrial disease or inborn errors of metabolism</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Early Warning Signs (Monitor Daily)</div>
                                    <div class="warning-detail">• <strong>Unexplained metabolic acidosis</strong> (base deficit >−10, lactate >2 mmol/L)<br>• <strong>↑ Triglycerides</strong> (>400 mg/dL) - lipemic serum<br>• <strong>↑ CK</strong> (>1000 U/L) or myoglobinuria<br>• <strong>New ECG changes:</strong> Bradycardia, Brugada-like pattern, QRS widening, ST elevation<br>• <strong>↑ Troponin</strong> without MI<br>• <strong>Hepatomegaly</strong> or ↑ liver enzymes<br><br><strong>If ANY early signs → strongly consider switching sedation</strong><sup class="ref-num">4</sup></div>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🔍 Prevention Strategies</div>
                                    <div class="info-detail">• <strong>Limit dose:</strong> Keep propofol <4 mg/kg/h (<67 mcg/kg/min) if prolonged infusion needed<sup class="ref-num">4</sup><br>• <strong>Limit duration:</strong> Avoid continuous infusion >48-72 hours<br>• <strong>Monitor:</strong> Daily CK, triglycerides, lactate if high-dose or prolonged infusion<br>• <strong>Alternative sedation:</strong> Use multimodal approach (combine with dexmedetomidine, ketamine, or benzodiazepines to reduce propofol dose)<br>• <strong>Adequate nutrition:</strong> Ensure carbohydrate intake (failure of fat oxidation is part of pathophysiology)</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Prognosis</div>
                                    <div class="warning-detail">• <strong>Mortality:</strong> ~30-60% once fully developed (cardiac arrest is often terminal)<sup class="ref-num">3</sup><br>• <strong>Key to survival:</strong> Early recognition and immediate cessation of propofol<br>• <strong>Recovery:</strong> If caught early (before cardiac failure), most metabolic abnormalities resolve within 24-72 hours after stopping propofol</div>
                                </div>

                                <!-- Inline References Section -->
                                <div class="protocol-references">
                                    <button class="references-toggle" onclick="event.stopPropagation(); toggleReferences(this)">
                                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <polyline points="6 9 12 15 18 9"></polyline>
                                        </svg>
                                        <span>View References (4)</span>
                                    </button>
                                    <div class="references-list">
                                        <div class="reference-item">
                                            <span class="reference-num">1</span>
                                            <span class="reference-citation">Krajčová A, et al. <strong>Propofol infusion syndrome: a structured review of experimental studies and 153 published case reports.</strong> <em>Crit Care</em>. 2015;19:398. PMID: 26563768</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">2</span>
                                            <span class="reference-citation">Otterspoor LC, et al. <strong>Update on the Propofol Infusion Syndrome in ICU Management of Patients with Head Injury.</strong> <em>Curr Opin Anaesthesiol</em>. 2008;21(5):544-551. PMID: 18784477</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">3</span>
                                            <span class="reference-citation">Mirrakhimov AE, et al. <strong>Propofol Infusion Syndrome in Adults: A Clinical Update.</strong> <em>Crit Care Res Pract</em>. 2015;2015:260385. PMID: 26078890</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">4</span>
                                            <span class="reference-citation">Roberts RJ, et al. <strong>Propofol concentration and the risk of infusion syndrome.</strong> <em>Anesth Analg</em>. 2009;109(4):1058-1062. PMID: 19762733</span>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- Massive Transfusion Protocol (ENHANCED) -->
                    <div class="protocol-card red" data-keywords="massive transfusion mtp hemorrhage shock trauma coagulopathy prbc ffp platelets blood products" onclick="toggleProtocol(this)">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">
                                    Massive Transfusion Protocol
                                    <span class="protocol-ref-count"><svg style="width:14px;height:14px;margin-right:2px;" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path><path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path><path d="M8 7h8"></path><path d="M8 11h8"></path></svg> 5 refs</span>
                                </h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Life-threatening hemorrhage requiring massive blood product replacement</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag immediate">Immediate</span>
                            <span class="protocol-tag call-help">Call Help</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Actions & Activation</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Activate MTP:</strong> Call blood bank immediately, activate institutional massive transfusion protocol<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Large-bore IV access:</strong> Two 14-16G peripheral IVs or central access (consider rapid infusion device)<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Control hemorrhage:</strong> Direct pressure, surgical hemostasis, consider damage control surgery<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Balanced resuscitation (1:1:1 ratio):</strong> PRBC : FFP : Platelets in equal proportions<sup class="ref-num">3</sup></li>
                                        <li class="protocol-step"><strong>Tranexamic acid (TXA):</strong> Give EARLY if traumatic hemorrhage (within 3 hours of injury)<sup class="ref-num">4</sup></li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Balanced Resuscitation (1:1:1 Ratio)<sup class="ref-num">3</sup></div>
                                    <div class="dose-detail"><strong>PRBC (packed RBCs):</strong> 1 unit (~350 mL) raises Hgb ~1 g/dL</div>
                                    <div class="dose-detail"><strong>FFP (fresh frozen plasma):</strong> 1 unit (~250 mL) per 1 unit PRBC</div>
                                    <div class="dose-detail"><strong>Platelets:</strong> 1 apheresis unit (or 6-pack pooled) per 6 units PRBC</div>
                                    <div class="dose-detail"><strong>Example:</strong> 6 PRBC : 6 FFP : 1 platelet apheresis unit</div>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Adjunct Medications<sup class="ref-num">4,5</sup></div>
                                    <div class="dose-detail"><strong>Tranexamic acid (TXA):</strong> 1 g IV over 10 min, then 1 g over 8 hours (trauma/obstetric hemorrhage)</div>
                                    <div class="dose-detail"><strong>Calcium chloride:</strong> 1-2 g IV q4-6 units citrated blood (prevent citrate toxicity)</div>
                                    <div class="dose-detail"><strong>Fibrinogen concentrate:</strong> 3-4 g IV if fibrinogen <150 mg/dL (alternative to cryoprecipitate)</div>
                                    <div class="dose-detail"><strong>Cryoprecipitate:</strong> 10 units (1 pool) if fibrinogen <100 mg/dL or ongoing bleeding</div>
                                    <div class="dose-detail"><strong>Prothrombin complex concentrate (PCC):</strong> 25-50 units/kg for warfarin reversal</div>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Monitoring & Labs</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Labs q30-60min:</strong> CBC, PT/INR, PTT, fibrinogen, ionized calcium, ABG/lactate, TEG/ROTEM if available<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Goals:</strong> Hgb >7 g/dL, platelets >50K (>100K if ongoing bleeding/CNS injury), INR <1.5, fibrinogen >150 mg/dL<sup class="ref-num">3</sup></li>
                                        <li class="protocol-step"><strong>Temperature:</strong> Maintain normothermia >36°C (use warmer, forced-air warming)<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Acidosis:</strong> Correct pH >7.2 (ventilation, bicarbonate if severe metabolic acidosis)<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Hypocalcemia:</strong> Maintain ionized Ca²⁺ >1.0 mmol/L (give calcium with citrated blood products)<sup class="ref-num">2</sup></li>
                                    </ol>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🔍 MTP Activation Criteria<sup class="ref-num">1</sup></div>
                                    <div class="info-detail"><strong>Clinical triggers:</strong><br>• Hemodynamic instability despite initial resuscitation<br>• Anticipated transfusion >10 units PRBC in 24h<br>• >4 units PRBC in 1 hour with ongoing bleeding<br>• Replacement of ≥50% blood volume in 3 hours<br><br><strong>Specific scenarios:</strong><br>• Trauma with shock (SBP <90, HR >120)<br>• Ruptured AAA, massive GI bleed, postpartum hemorrhage<br>• Intraoperative hemorrhage (liver trauma, aortic surgery, placenta accreta)</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ The "Lethal Triad" (Trauma Death Spiral)</div>
                                    <div class="warning-detail">1. <strong>Hypothermia</strong> (<36°C) - worsens coagulopathy, impairs platelet function<br>2. <strong>Acidosis</strong> (pH <7.2) - worsens coagulopathy, impairs hemostasis<br>3. <strong>Coagulopathy</strong> (INR >1.5, platelets <50K) - perpetuates bleeding<sup class="ref-num">2</sup><br><br><strong>Prevention is critical:</strong> Aggressive warming, early balanced transfusion, damage control surgery, permissive hypotension (SBP 80-90 until hemorrhage control)</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Complications of Massive Transfusion</div>
                                    <div class="warning-detail">• <strong>TACO (Transfusion-Associated Circulatory Overload):</strong> Pulmonary edema, hypoxia - slow transfusion, diuretics<br>• <strong>TRALI (Transfusion-Related Acute Lung Injury):</strong> Noncardiogenic pulmonary edema within 6h - supportive care<br>• <strong>Citrate toxicity:</strong> Hypocalcemia (tremor, arrhythmias) - give calcium chloride<br>• <strong>Hyperkalemia:</strong> From stored blood (especially old units) - monitor K⁺, treat if >6 mEq/L<br>• <strong>Hypothermia:</strong> Use blood warmers, forced-air warming<br>• <strong>Dilutional coagulopathy:</strong> Prevented by balanced 1:1:1 resuscitation<sup class="ref-num">5</sup></div>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🔍 Special Populations</div>
                                    <div class="info-detail"><strong>Obstetric hemorrhage (PPH):</strong><br>• TXA 1 g IV within 3h of delivery<br>• Consider recombinant factor VIIa if refractory (off-label)<br>• Uterotonics: oxytocin, methylergonovine, carboprost, misoprostol<br><br><strong>Anticoagulated patients:</strong><br>• Warfarin → PCC 25-50 units/kg + vitamin K 10 mg IV<br>• Dabigatran → idarucizumab 5 g IV<br>• Xa inhibitors (rivaroxaban, apixaban) → andexanet alfa or PCC<sup class="ref-num">5</sup></div>
                                </div>

                                <!-- Inline References Section -->
                                <div class="protocol-references">
                                    <button class="references-toggle" onclick="event.stopPropagation(); toggleReferences(this)">
                                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <polyline points="6 9 12 15 18 9"></polyline>
                                        </svg>
                                        <span>View References (5)</span>
                                    </button>
                                    <div class="references-list">
                                        <div class="reference-item">
                                            <span class="reference-num">1</span>
                                            <span class="reference-citation">Holcomb JB, et al. <strong>The Prospective, Observational, Multicenter, Major Trauma Transfusion (PROMMTT) Study.</strong> <em>JAMA Surg</em>. 2013;148(2):127-136. PMID: 23560283</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">2</span>
                                            <span class="reference-citation">Spahn DR, et al. <strong>The European guideline on management of major bleeding and coagulopathy following trauma.</strong> <em>Crit Care</em>. 2019;23(1):98. PMID: 30917843</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">3</span>
                                            <span class="reference-citation">Holcomb JB, et al. <strong>Transfusion of plasma, platelets, and red blood cells in a 1:1:1 vs 1:1:2 ratio (PROPPR trial).</strong> <em>JAMA</em>. 2015;313(5):471-482. PMID: 25647203</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">4</span>
                                            <span class="reference-citation">CRASH-2 Trial Collaborators. <strong>Effects of tranexamic acid on death in trauma patients (CRASH-2).</strong> <em>Lancet</em>. 2010;376(9734):23-32. PMID: 20554319</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">5</span>
                                            <span class="reference-citation">Ghadimi K, et al. <strong>Perioperative management of the bleeding patient.</strong> <em>Br J Anaesth</em>. 2016;117(suppl 3):iii18-iii30. PMID: 27940452</span>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>

                </div>
            </div>

            <!-- Airway Emergencies -->
            <div class="category-section" data-category="airway">
                <div class="category-header">
                    <div class="category-icon orange">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <path d="M9 11a3 3 0 1 0 6 0a3 3 0 0 0 -6 0"/>
                            <path d="M12 12v6"/>
                            <path d="M6 15c-2.213 -1.246 -3.5 -3.154 -3.5 -5.294 0 -3.314 2.686 -6 6 -6h7c3.314 0 6 2.686 6 6 0 2.14 -1.287 4.048 -3.5 5.294"/>
                        </svg>
                    </div>
                    <h2 class="category-title">Airway Emergencies</h2>
                </div>
                <div class="protocols-grid">

                    <!-- CICO (ENHANCED) -->
                    <div class="protocol-card orange" data-keywords="cico cant intubate oxygenate difficult airway cricothyroidotomy emergency front neck access scalpel bougie" onclick="toggleProtocol(this)">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">
                                    Can't Intubate, Can't Oxygenate (CICO)
                                    <span class="protocol-ref-count"><svg style="width:14px;height:14px;margin-right:2px;" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path><path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path><path d="M8 7h8"></path><path d="M8 11h8"></path></svg> 2 refs</span>
                                </h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Life-threatening failure to intubate AND ventilate requiring emergency surgical airway</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag immediate">Immediate</span>
                            <span class="protocol-tag call-help">Call Help</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Actions - Scalpel Cricothyroidotomy</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Declare CICO:</strong> Call for help, assign roles<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Position:</strong> Extend neck, palpate cricothyroid membrane</li>
                                        <li class="protocol-step"><strong>Stab incision:</strong> Horizontal skin incision through cricothyroid membrane with scalpel (blade #10 or #20)<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Bougie:</strong> Insert bougie through membrane into trachea (feel "clicks" of tracheal rings)</li>
                                        <li class="protocol-step"><strong>Tube:</strong> Railroad 6.0 cuffed ETT or tracheostomy tube over bougie into trachea</li>
                                        <li class="protocol-step"><strong>Confirm:</strong> Inflate cuff, ventilate, confirm placement with EtCO₂</li>
                                    </ol>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Key Points</div>
                                    <div class="warning-detail">• Scalpel technique preferred over needle cricothyroidotomy (Seldinger kits have high failure rate)<sup class="ref-num">2</sup><br>• Don't delay - permanent brain damage occurs after 3-5 minutes of hypoxia<br>• If anatomy unclear: Make vertical skin incision, then palpate membrane and make horizontal membrane incision<br>• Post-procedure: Secure tube, get ENT/surgery consult, chest X-ray</div>
                                </div>

                                <!-- Inline References Section -->
                                <div class="protocol-references">
                                    <button class="references-toggle" onclick="event.stopPropagation(); toggleReferences(this)">
                                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <polyline points="6 9 12 15 18 9"></polyline>
                                        </svg>
                                        <span>View References (2)</span>
                                    </button>
                                    <div class="references-list">
                                        <div class="reference-item">
                                            <span class="reference-num">1</span>
                                            <span class="reference-citation">Apfelbaum JL, et al. <strong>2022 ASA Practice Guidelines for Management of the Difficult Airway.</strong> <em>Anesthesiology</em>. 2022;136(1):31-81. PMID: 34762729</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">2</span>
                                            <span class="reference-citation">Frerk C, et al. <strong>Difficult Airway Society Guidelines for Emergency Front-of-Neck Access.</strong> <em>Br J Anaesth</em>. 2015;115(6):827-848. PMID: 26556848</span>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- Laryngospasm (ENHANCED) -->
                    <div class="protocol-card orange" data-keywords="laryngospasm larynx spasm stridor negative pressure pulmonary edema npppe vocal cord glottic closure" onclick="toggleProtocol(this)">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">
                                    Laryngospasm
                                    <span class="protocol-ref-count"><svg style="width:14px;height:14px;margin-right:2px;" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path><path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path><path d="M8 7h8"></path><path d="M8 11h8"></path></svg> 4 refs</span>
                                </h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Reflex glottic closure causing complete or partial airway obstruction</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag immediate">Immediate</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Actions</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Remove stimulus:</strong> Stop surgery, suction oropharynx of blood/secretions<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>100% O₂ with positive pressure:</strong> Gentle jaw thrust + CPAP (try 5-10 cm H₂O first)<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Deepen anesthesia:</strong> If inadequate depth, give propofol 0.5-1 mg/kg IV<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Larson's maneuver:</strong> Firm pressure on "laryngospasm notch" (posterior to mandible angle, anterior to mastoid) while applying jaw thrust<sup class="ref-num">3</sup></li>
                                    </ol>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ If Laryngospasm Persists (>30 seconds)</div>
                                    <div class="warning-detail"><strong>Give succinylcholine:</strong> 0.1-0.5 mg/kg IV (or 2-4 mg/kg IM if no IV access)<sup class="ref-num">2</sup><br>• Prepare to ventilate and intubate if needed<br>• Monitor for bradycardia (especially in children) - have atropine ready</div>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Key Dosing<sup class="ref-num">2</sup></div>
                                    <div class="dose-detail"><strong>Propofol (deepen):</strong> 0.5-1 mg/kg IV bolus</div>
                                    <div class="dose-detail"><strong>Succinylcholine (if refractory):</strong> 0.1-0.5 mg/kg IV or 2-4 mg/kg IM</div>
                                    <div class="dose-detail"><strong>Atropine (if bradycardia):</strong> 0.01-0.02 mg/kg IV (minimum 0.1 mg)</div>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🔍 Post-Laryngospasm Management</div>
                                    <div class="info-detail">• <strong>Monitor for negative-pressure pulmonary edema (NPPE):</strong> Occurs in ~0.1% of laryngospasm cases<sup class="ref-num">4</sup><br>• Signs: Pink frothy sputum, decreased SpO₂, crackles on auscultation<br>• Treatment: Supplemental O₂, PEEP/CPAP, diuretics if needed, rarely intubation<br>• Monitor oxygenation for 1-2 hours post-event</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Prevention Strategies</div>
                                    <div class="warning-detail">• Extubate <strong>deep</strong> (under anesthesia) or <strong>awake</strong> - avoid "light" stage<sup class="ref-num">1</sup><br>• Suction oropharynx thoroughly before emergence (avoid pharyngeal stimulation during light anesthesia)<br>• Lidocaine 1-1.5 mg/kg IV 2 min before extubation may reduce risk<sup class="ref-num">2</sup><br>• Higher risk: Pediatrics, airway surgery, GERD, recent URI, reactive airway disease</div>
                                </div>

                                <!-- Inline References Section -->
                                <div class="protocol-references">
                                    <button class="references-toggle" onclick="event.stopPropagation(); toggleReferences(this)">
                                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <polyline points="6 9 12 15 18 9"></polyline>
                                        </svg>
                                        <span>View References (4)</span>
                                    </button>
                                    <div class="references-list">
                                        <div class="reference-item">
                                            <span class="reference-num">1</span>
                                            <span class="reference-citation">Orliaguet GA, et al. <strong>Management of laryngospasm in children.</strong> <em>Paediatr Anaesth</em>. 2019;29(7):774-780. PMID: 31025445</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">2</span>
                                            <span class="reference-citation">Visvanathan T, et al. <strong>Laryngospasm in anaesthesia.</strong> <em>Contin Educ Anaesth Crit Care Pain</em>. 2015;15(3):136-141.</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">3</span>
                                            <span class="reference-citation">Larson CP Jr. <strong>Laryngospasm--the best treatment.</strong> <em>Anesthesiology</em>. 1998;89(5):1293-1294. PMID: 9822034</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">4</span>
                                            <span class="reference-citation">Bhattacharya M, et al. <strong>Negative Pressure Pulmonary Edema.</strong> <em>Chest</em>. 2016;150(4):927-933. PMID: 27167224</span>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- Bronchospasm (ENHANCED) -->
                    <div class="protocol-card orange" data-keywords="bronchospasm wheezing asthma reactive airway bronchoconstriction albuterol beta agonist" onclick="toggleProtocol(this)">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">
                                    Bronchospasm
                                    <span class="protocol-ref-count"><svg style="width:14px;height:14px;margin-right:2px;" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path><path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path><path d="M8 7h8"></path><path d="M8 11h8"></path></svg> 4 refs</span>
                                </h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Acute airway obstruction from smooth muscle constriction and inflammation</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag immediate">Act Fast</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Actions</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>100% O₂:</strong> Increase FiO₂ to 1.0, ensure adequate oxygenation<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Deepen anesthesia:</strong> Increase volatile anesthetic (sevoflurane/isoflurane have bronchodilator effects)<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Rule out mechanical causes:</strong> Check for kinked ETT, mucus plug, endobronchial intubation, pneumothorax<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Beta-2 agonist (first-line):</strong> Albuterol 4-8 puffs via MDI with spacer into circuit (or 2.5-5 mg nebulized)<sup class="ref-num">1</sup></li>
                                    </ol>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Additional Therapies (if refractory)</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Epinephrine:</strong> 10-50 mcg IV boluses (or 0.3 mg IM if severe)<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Ketamine:</strong> 0.5-1 mg/kg IV bolus (bronchodilator via sympathomimetic effects)<sup class="ref-num">3</sup></li>
                                        <li class="protocol-step"><strong>Magnesium sulfate:</strong> 2 g IV over 20 min (smooth muscle relaxation)<sup class="ref-num">3</sup></li>
                                        <li class="protocol-step"><strong>Corticosteroids:</strong> Methylprednisolone 1-2 mg/kg IV or hydrocortisone 2-4 mg/kg IV (delayed onset ~6 hours)<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Anticholinergic:</strong> Ipratropium 0.5 mg nebulized (adjunct to beta-agonist)<sup class="ref-num">1</sup></li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Key Medications<sup class="ref-num">1,2,3</sup></div>
                                    <div class="dose-detail"><strong>Albuterol (MDI):</strong> 4-8 puffs into circuit with spacer, repeat q20min PRN</div>
                                    <div class="dose-detail"><strong>Albuterol (nebulized):</strong> 2.5-5 mg in 3 mL NS</div>
                                    <div class="dose-detail"><strong>Epinephrine:</strong> 10-50 mcg IV boluses (titrate) or 0.3 mg IM if severe</div>
                                    <div class="dose-detail"><strong>Ketamine:</strong> 0.5-1 mg/kg IV bolus, then 0.5-1 mg/kg/h infusion</div>
                                    <div class="dose-detail"><strong>Magnesium sulfate:</strong> 2 g (40 mg/kg peds) IV over 20 minutes</div>
                                    <div class="dose-detail"><strong>Methylprednisolone:</strong> 1-2 mg/kg IV (or hydrocortisone 2-4 mg/kg)</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Differential Diagnosis</div>
                                    <div class="warning-detail"><strong>Rule out these mechanical causes first:</strong><br>• Kinked or obstructed ETT<br>• Endobronchial intubation (check bilateral breath sounds)<br>• Mucus plug (consider bronchoscopy/suction)<br>• Tension pneumothorax<br>• Pulmonary edema or aspiration<br>• Anaphylaxis (check for hypotension, urticaria)<sup class="ref-num">2</sup></div>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🔍 Ventilation Strategies</div>
                                    <div class="info-detail">• <strong>Permissive hypercapnia:</strong> Accept higher PaCO₂ to avoid barotrauma<sup class="ref-num">4</sup><br>• <strong>Prolonged expiratory time:</strong> Decrease RR, increase I:E ratio to 1:3 or 1:4<br>• <strong>Avoid high peak pressures:</strong> Use pressure-control if needed, minimize auto-PEEP<br>• <strong>Manual ventilation:</strong> Gives better "feel" for airway resistance and allows slower rates</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Prevention & Risk Factors</div>
                                    <div class="warning-detail">• <strong>High-risk patients:</strong> Asthma, COPD, smokers, recent URI, reactive airway disease<br>• <strong>Preoperative optimization:</strong> Optimize asthma control, consider preop bronchodilators<br>• <strong>Avoid triggers:</strong> Deep extubation, avoid histamine-releasing drugs (morphine, atracurium), minimize airway instrumentation<br>• <strong>LMA vs ETT:</strong> Consider LMA for low-risk surgery to reduce airway stimulation<sup class="ref-num">1</sup></div>
                                </div>

                                <!-- Inline References Section -->
                                <div class="protocol-references">
                                    <button class="references-toggle" onclick="event.stopPropagation(); toggleReferences(this)">
                                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <polyline points="6 9 12 15 18 9"></polyline>
                                        </svg>
                                        <span>View References (4)</span>
                                    </button>
                                    <div class="references-list">
                                        <div class="reference-item">
                                            <span class="reference-num">1</span>
                                            <span class="reference-citation">Woods BD, Sladen RN. <strong>Perioperative considerations for the patient with asthma and bronchospasm.</strong> <em>Br J Anaesth</em>. 2009;103 Suppl 1:i57-65. PMID: 20007991</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">2</span>
                                            <span class="reference-citation">Mitsuhata H, et al. <strong>Mechanisms and management of intraoperative bronchospasm.</strong> <em>Curr Opin Anaesthesiol</em>. 1996;9(3):238-242.</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">3</span>
                                            <span class="reference-citation">Rodrigo GJ, et al. <strong>Acute asthma in adults: a review.</strong> <em>Chest</em>. 2004;125(3):1081-1102. PMID: 15006973</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">4</span>
                                            <span class="reference-citation">Scalese MJ, et al. <strong>Severe Refractory Status Asthmaticus: A Review.</strong> <em>J Intensive Care Med</em>. 2020;35(10):977-988. PMID: 30987529</span>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- Airway Fire (ENHANCED) -->
                    <div class="protocol-card orange" data-keywords="airway fire operating room fire laser surgery oxygen combustion ett burning thermal injury" onclick="toggleProtocol(this)">
                        <div class="protocol-header">
                            <div>
                                <h3 class="protocol-title">
                                    Airway Fire
                                    <span class="protocol-ref-count"><svg style="width:14px;height:14px;margin-right:2px;" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path><path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path><path d="M8 7h8"></path><path d="M8 11h8"></path></svg> 3 refs</span>
                                </h3>
                            </div>
                            <div class="expand-icon">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="6 9 12 15 18 9"></polyline>
                                </svg>
                            </div>
                        </div>
                        <p class="protocol-summary">Combustion in airway from heat source + oxygen-enriched environment + flammable material</p>
                        <div class="protocol-tags">
                            <span class="protocol-tag immediate">Immediate</span>
                            <span class="protocol-tag call-help">Call Help</span>
                        </div>
                        <div class="protocol-content">
                            <div class="protocol-details">
                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Immediate Actions (Remember: STOP-DROP-ROLL)</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>STOP gas flow:</strong> Immediately disconnect O₂ source and stop ventilation<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>DROP the ETT:</strong> Remove burning endotracheal tube from airway<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>ROLL patient to side:</strong> Pour saline into airway and oropharynx to extinguish fire<sup class="ref-num">1</sup></li>
                                        <li class="protocol-step"><strong>Mask ventilate with air or 21% O₂:</strong> Resume ventilation with lowest FiO₂ possible<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Reintubate:</strong> Use new ETT (smaller size if edema present), assess damage with laryngoscopy/bronchoscopy<sup class="ref-num">1</sup></li>
                                    </ol>
                                </div>

                                <div class="protocol-section">
                                    <h4 class="protocol-section-title">Secondary Assessment & Management</h4>
                                    <ol class="protocol-steps">
                                        <li class="protocol-step"><strong>Bronchoscopy:</strong> Assess airway injury severity (mucosal burns, soot, debris)<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Saline lavage:</strong> Irrigate airway to remove debris and carbonized material<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>Steroids:</strong> Consider dexamethasone 0.5 mg/kg IV (controversial, may reduce edema)<sup class="ref-num">2</sup></li>
                                        <li class="protocol-step"><strong>ICU admission:</strong> Monitor for delayed airway edema (peaks 12-24h), ARDS, pneumonia<sup class="ref-num">3</sup></li>
                                        <li class="protocol-step"><strong>No extubation:</strong> Keep intubated until edema resolves (typically 3-7 days)<sup class="ref-num">3</sup></li>
                                    </ol>
                                </div>

                                <div class="dose-box">
                                    <div class="dose-box-title">💊 Post-Fire Management</div>
                                    <div class="dose-detail"><strong>Dexamethasone (controversial):</strong> 0.5 mg/kg IV (may reduce edema)</div>
                                    <div class="dose-detail"><strong>Bronchodilators:</strong> Albuterol if bronchospasm develops</div>
                                    <div class="dose-detail"><strong>Antibiotics:</strong> NOT routinely indicated (no benefit for prophylaxis)</div>
                                    <div class="dose-detail"><strong>Ventilation strategy:</strong> Lung-protective ventilation, lowest FiO₂ maintaining SpO₂ >90%</div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Fire Triangle (All 3 Required for Fire)</div>
                                    <div class="warning-detail">1. <strong>Oxidizer:</strong> O₂ (>21%), N₂O<br>2. <strong>Fuel:</strong> ETT, surgical drapes, gauze, alcohol prep solutions<br>3. <strong>Ignition source:</strong> Electrocautery, laser, fiberoptic light<br><br><strong>Break one element to prevent fire</strong> - Lower FiO₂ to ≤30% during cautery near airway, use laser-safe ETT, keep flammables wet<sup class="ref-num">1</sup></div>
                                </div>

                                <div class="info-box">
                                    <div class="info-box-title">🔍 High-Risk Procedures</div>
                                    <div class="info-detail">• <strong>Head & neck surgery:</strong> Laser laryngoscopy, tonsillectomy, tracheostomy<br>• <strong>Oral/facial surgery:</strong> Cautery near airway<br>• <strong>ENT procedures:</strong> Any procedure with electrocautery + high FiO₂<br>• <strong>Prevention:</strong> Use lowest safe FiO₂ (≤30%), laser-resistant ETTs, wet sponges around surgical field, pause O₂ during cautery<sup class="ref-num">1</sup></div>
                                </div>

                                <div class="warning-box">
                                    <div class="warning-box-title">⚠️ Post-Fire Complications</div>
                                    <div class="warning-detail">• <strong>Acute:</strong> Airway edema (delayed, peaks 12-24h), ARDS, pneumonia<br>• <strong>Delayed:</strong> Tracheal stenosis, granulation tissue formation<br>• <strong>Monitoring:</strong> Serial bronchoscopies, prolonged intubation (3-7 days typical)<br>• <strong>Tracheostomy:</strong> May be needed if severe injury or prolonged intubation anticipated<sup class="ref-num">3</sup></div>
                                </div>

                                <!-- Inline References Section -->
                                <div class="protocol-references">
                                    <button class="references-toggle" onclick="event.stopPropagation(); toggleReferences(this)">
                                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                            <polyline points="6 9 12 15 18 9"></polyline>
                                        </svg>
                                        <span>View References (3)</span>
                                    </button>
                                    <div class="references-list">
                                        <div class="reference-item">
                                            <span class="reference-num">1</span>
                                            <span class="reference-citation">ASA Task Force on Operating Room Fires. <strong>Practice Advisory for the Prevention and Management of Operating Room Fires.</strong> <em>Anesthesiology</em>. 2013;118(2):271-290. PMID: 23287706</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">2</span>
                                            <span class="reference-citation">Pruitt BA, Cioffi WG. <strong>Management of burns in the airway and face.</strong> <em>Clin Plast Surg</em>. 2009;36(4):555-567. PMID: 19793552</span>
                                        </div>
                                        <div class="reference-item">
                                            <span class="reference-num">3</span>
                                            <span class="reference-citation">Worley SL. <strong>Fire safety in the operating room.</strong> <em>AORN J</em>. 2012;95(5):606-618. PMID: 22541771</span>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>


                </div>
            </div>

                </div>
            </div>
        </div>

        <!-- Footer -->
        <footer class="footer">
            <div class="footer-inner">
                <span class="footer-text">© 2025 GasConsult.ai</span>
                <div class="footer-links">
                    <a href="/privacy" class="footer-link">Privacy</a>
                    <a href="/terms" class="footer-link">Terms</a>
                    <a href="mailto:contact@gasconsult.ai" class="footer-link">Contact</a>
                </div>
            </div>
        </footer>
    </div>

    <script>
        function toggleMobileMenu() {
            const menu = document.getElementById('mobileMenu');
            const btn = document.querySelector('.mobile-menu-btn');
            menu.classList.toggle('active');
            btn.classList.toggle('active');
        }

        function toggleNavDropdown(event) {
            event.stopPropagation();
            const menu = event.target.nextElementSibling;
            menu.classList.toggle('show');
        }

        document.addEventListener('click', function(event) {
            if (!event.target.matches('.nav-dropdown-toggle')) {
                document.querySelectorAll('.nav-dropdown-menu').forEach(menu => {
                    menu.classList.remove('show');
                });
            }
        });

        function toggleProtocol(card) {
            const wasExpanded = card.classList.contains('expanded');

            // Close ALL protocol cards first
            document.querySelectorAll('.protocol-card').forEach(c => {
                c.classList.remove('expanded');
            });

            // If this card wasn't already open, open it
            if (!wasExpanded) {
                card.classList.add('expanded');
            }
        }

        function toggleReferences(button) {
            button.classList.toggle('expanded');
            const refList = button.nextElementSibling;
            refList.classList.toggle('show');
        }

        function filterProtocols() {
            const input = document.getElementById('searchInput');
            const filter = input.value.toLowerCase();
            const clearBtn = document.getElementById('clearBtn');
            const cards = document.querySelectorAll('.protocol-card');
            const categories = document.querySelectorAll('.category-section');

            // Show/hide clear button
            if (filter.length > 0) {
                clearBtn.classList.add('visible');
            } else {
                clearBtn.classList.remove('visible');
            }

            // Filter protocol cards
            let hasVisibleCards = new Set();

            cards.forEach(card => {
                const keywords = card.getAttribute('data-keywords') || '';
                const title = card.querySelector('.protocol-title').textContent;
                const summary = card.querySelector('.protocol-summary').textContent;

                const searchableText = (keywords + ' ' + title + ' ' + summary).toLowerCase();

                if (searchableText.includes(filter) || filter === '') {
                    card.style.display = '';
                    hasVisibleCards.add(card.closest('.category-section'));
                } else {
                    card.style.display = 'none';
                }
            });

            // Show/hide category sections
            categories.forEach(category => {
                if (hasVisibleCards.has(category) || filter === '') {
                    category.style.display = '';
                } else {
                    category.style.display = 'none';
                }
            });
        }

        function clearSearch() {
            document.getElementById('searchInput').value = '';
            filterProtocols();
        }
    </script>

</body>
</html>
"""


QUICK_DOSE_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <!-- Google Analytics 4 -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=G-01NZYD1DPP"></script>
    <script>
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-01NZYD1DPP');
    </script>

    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Quick Dose - GasConsult.ai</title>

    <!-- SEO Meta Tags -->
    <meta name="description" content="Fast drug dosing reference for anesthesia: induction agents, opioids, neuromuscular blockers, vasopressors, and emergency medications with evidence-based dosing.">

    <!-- Open Graph / Facebook -->
    <meta property="og:type" content="website">
    <meta property="og:url" content="https://gasconsult.ai/quick-dose">
    <meta property="og:title" content="Quick Dose - GasConsult.ai">
    <meta property="og:description" content="Fast drug dosing reference for anesthesia with evidence-based dosing guidelines.">
    <meta property="og:image" content="https://gasconsult.ai/static/logo.png">

    <!-- Twitter -->
    <meta property="twitter:card" content="summary_large_image">
    <meta property="twitter:url" content="https://gasconsult.ai/quick-dose">
    <meta property="twitter:title" content="Quick Dose - GasConsult.ai">
    <meta property="twitter:description" content="Fast drug dosing reference for anesthesia with evidence-based dosing guidelines.">
    <meta property="twitter:image" content="https://gasconsult.ai/static/logo.png">

    <!-- PWA -->
    <link rel="icon" type="image/svg+xml" href="/static/favicon.svg?v=6">
    <link rel="apple-touch-icon" href="/static/favicon.svg?v=6">
    <link rel="manifest" href="/static/manifest.json">
    <meta name="theme-color" content="#2563EB">

    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap" rel="stylesheet">
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

            /* ISO 26825 Syringe Colors */
            --color-induction: #FFFFFF;
            --color-induction-border: #e5e7eb;
            --color-induction-text: #374151;

            --color-opioid: #89CFF0;
            --color-opioid-dark: #0EA5E9;
            --color-opioid-text: #0369a1;

            --color-paralytic: #FF6B6B;
            --color-paralytic-dark: #DC2626;
            --color-paralytic-text: #b91c1c;

            --color-reversal: #E5E7EB;
            --color-reversal-border: #d1d5db;
            --color-reversal-text: #4b5563;

            --color-vasopressor: #DDA0DD;
            --color-vasopressor-dark: #A855F7;
            --color-vasopressor-text: #7c3aed;

            --color-hypnotic: #FFA500;
            --color-hypnotic-dark: #EA580C;
            --color-hypnotic-text: #c2410c;

            --color-benzo: #FFA500;
            --color-benzo-dark: #EA580C;
            --color-benzo-text: #c2410c;

            --color-local: #808080;
            --color-local-dark: #525252;
            --color-local-text: #404040;

            --color-antiemetic: #E2CFBE;
            --color-antiemetic-dark: #A78B71;
            --color-antiemetic-text: #78716c;

            --color-other: #f3f4f6;
            --color-other-border: #e5e7eb;
            --color-other-text: #6b7280;
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

        /* Background Canvas */
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

        @keyframes fade-up {
            from { opacity: 0; transform: translateY(20px); }
            to { opacity: 1; transform: translateY(0); }
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

        /* Navigation - Homepage Style */
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

        .nav-dropdown:has(.nav-dropdown-link.active) .nav-dropdown-toggle {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .nav-dropdown {
            position: relative;
            display: inline-block;
        }

        .nav-dropdown-toggle {
            cursor: pointer;
            background: none;
            border: none;
            font-family: inherit;
        }

        .nav-dropdown-menu {
            display: none;
            position: absolute;
            top: 100%;
            right: 0;
            background: white;
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            min-width: 200px;
            margin-top: 4px;
            z-index: 1000;
            overflow: hidden;
        }

        .nav-dropdown-menu.show {
            display: block;
        }

        .nav-dropdown-link {
            display: block;
            padding: 10px 14px;
            color: var(--gray-600);
            text-decoration: none;
            border-radius: 8px;
            font-size: 14px;
            font-weight: 500;
            transition: all 0.2s ease;
        }

        .nav-dropdown-link:hover {
            background: var(--gray-50);
            color: var(--gray-900);
        }

        .mobile-menu-btn {
            display: flex;
            flex-direction: column;
            gap: 5px;
            background: none;
            border: none;
            cursor: pointer;
            padding: 8px;
        }

        .mobile-menu-btn span {
            display: block;
            width: 22px;
            height: 2px;
            background: var(--gray-700);
            transition: all 0.3s ease;
        }

        .mobile-menu-btn.active span:nth-child(1) {
            transform: rotate(45deg) translate(6px, 6px);
        }

        .mobile-menu-btn.active span:nth-child(2) {
            opacity: 0;
        }

        .mobile-menu-btn.active span:nth-child(3) {
            transform: rotate(-45deg) translate(6px, -6px);
        }

        .mobile-menu {
            position: fixed;
            top: 80px;
            left: 16px;
            right: 16px;
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(20px);
            border: 1px solid var(--gray-200);
            border-radius: 16px;
            padding: 16px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.1);
            display: none;
            flex-direction: column;
            gap: 4px;
            z-index: 99;
        }

        .mobile-menu.active {
            display: flex;
        }

        .mobile-menu-link {
            padding: 12px 16px;
            color: var(--gray-700);
            text-decoration: none;
            border-radius: 12px;
            font-size: 15px;
            font-weight: 500;
            transition: all 0.2s ease;
        }

        .mobile-menu-link:hover {
            background: var(--gray-100);
            color: var(--gray-900);
        }

        /* Content Area - With top padding for fixed nav */
        .container {
            max-width: 700px;
            margin: 0 auto;
            padding: 130px 16px 80px;
        }

        /* Weight Section - Mobile First */
        .weight-section {
            text-align: center;
            margin-bottom: 32px;
            animation: fade-up 0.6s ease forwards;
            opacity: 0;
        }

        .weight-label {
            font-size: 13px;
            font-weight: 600;
            color: #94a3b8;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            margin-bottom: 12px;
        }

        .weight-input-container {
            display: inline-flex;
            align-items: baseline;
            gap: 6px;
            margin-bottom: 20px;
        }

        /* Mobile: Smaller weight input for easier typing */
        .weight-input {
            background: transparent;
            border: none;
            border-bottom: 3px solid #2563eb;
            color: #0f172a;
            font-size: 48px;
            font-weight: 800;
            width: 110px;
            text-align: center;
            outline: none;
            letter-spacing: -2px;
            transition: border-color 0.2s;
            min-height: 44px;
        }

        .weight-input:focus {
            border-color: #1d4ed8;
        }

        .weight-unit {
            font-size: 20px;
            font-weight: 600;
            color: #94a3b8;
        }

        .quick-weights {
            display: flex;
            justify-content: center;
            gap: 8px;
            flex-wrap: wrap;
        }

        /* Mobile: Larger touch targets (44px minimum) */
        .quick-weight-btn {
            background: white;
            border: 1px solid #e2e8f0;
            color: #64748b;
            padding: 12px 16px;
            border-radius: 100px;
            font-size: 14px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s ease;
            min-height: 44px;
            min-width: 60px;
        }

        .quick-weight-btn:hover {
            border-color: #cbd5e1;
            color: #475569;
        }

        .quick-weight-btn:active {
            transform: scale(0.95);
        }

        .quick-weight-btn.active {
            background: #2563eb;
            border-color: #2563eb;
            color: white;
        }

        /* Color Legend - Mobile First */
        .color-legend {
            display: flex;
            flex-wrap: wrap;
            justify-content: center;
            gap: 8px;
            margin-bottom: 24px;
            padding: 12px;
            background: white;
            border-radius: 12px;
            border: 1px solid #e5e7eb;
        }

        .legend-item {
            display: flex;
            align-items: center;
            gap: 6px;
            font-size: 10px;
            color: #64748b;
        }

        .legend-swatch {
            width: 16px;
            height: 16px;
            border-radius: 4px;
            border: 1px solid rgba(0,0,0,0.1);
        }

        .legend-swatch.opioid { background: var(--color-opioid); }
        .legend-swatch.paralytic { background: var(--color-paralytic); }
        .legend-swatch.hypnotic { background: var(--color-hypnotic); }
        .legend-swatch.vasopressor { background: var(--color-vasopressor); }
        .legend-swatch.reversal { background: var(--color-reversal); border-color: #d1d5db; }
        .legend-swatch.local { background: var(--color-local); }

        /* Drug Category Section */
        .category-group {
            margin-bottom: 32px;
            animation: fade-up 0.6s ease forwards;
            opacity: 0;
        }

        .category-header {
            display: flex;
            align-items: center;
            gap: 10px;
            margin-bottom: 12px;
            padding: 0 4px;
        }

        .category-swatch {
            width: 14px;
            height: 14px;
            border-radius: 4px;
            border: 1px solid rgba(0,0,0,0.1);
        }

        .category-swatch.induction { background: var(--color-hypnotic); }
        .category-swatch.opioid { background: var(--color-opioid); }
        .category-swatch.paralytic { background: var(--color-paralytic); }
        .category-swatch.reversal { background: var(--color-reversal); border-color: #ccc; }
        .category-swatch.vasopressor { background: var(--color-vasopressor); }
        .category-swatch.other { background: var(--color-other); border-color: #ddd; }

        .category-title {
            font-size: 13px;
            font-weight: 700;
            color: #64748b;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        /* Drug List */
        .drug-list {
            display: flex;
            flex-direction: column;
            gap: 8px;
        }

        /* Drug Item */
        .drug-item {
            background: white;
            border-radius: 14px;
            overflow: hidden;
            transition: all 0.25s ease;
            border: 2px solid #e5e7eb;
            animation: fade-up 0.6s ease forwards;
            opacity: 0;
        }

        .drug-item.induction { border-left: 4px solid var(--color-hypnotic); }
        .drug-item.opioid { border-left: 4px solid var(--color-opioid-dark); }
        .drug-item.paralytic { border-left: 4px solid var(--color-paralytic-dark); }
        .drug-item.reversal { border-left: 4px solid #9ca3af; }
        .drug-item.vasopressor { border-left: 4px solid var(--color-vasopressor-dark); }
        .drug-item.other { border-left: 4px solid #9ca3af; }

        .drug-item:hover {
            border-color: #d1d5db;
            box-shadow: 0 2px 12px rgba(0, 0, 0, 0.04);
        }

        .drug-item.expanded {
            box-shadow: 0 4px 24px rgba(0, 0, 0, 0.08);
        }

        /* Mobile: Smaller padding, better touch targets */
        .drug-header {
            display: flex;
            align-items: center;
            justify-content: space-between;
            padding: 14px 16px;
            cursor: pointer;
            user-select: none;
            min-height: 56px;
        }

        .drug-name {
            font-size: 16px;
            font-weight: 600;
            color: #1e293b;
        }

        .drug-quick-dose {
            display: flex;
            align-items: center;
            gap: 8px;
        }

        .drug-standard-dose {
            display: flex;
            align-items: baseline;
            gap: 3px;
            background: #f1f5f9;
            padding: 6px 12px;
            border-radius: 8px;
        }

        .drug-quick-value {
            font-size: 17px;
            font-weight: 700;
            color: #0f172a;
        }

        .drug-quick-unit {
            font-size: 12px;
            color: #64748b;
            font-weight: 500;
        }

        .drug-chevron {
            color: #94a3b8;
            transition: transform 0.3s ease;
        }

        .drug-item.expanded .drug-chevron {
            transform: rotate(180deg);
        }

        /* Drug Details - Expandable */
        .drug-details {
            max-height: 0;
            overflow: hidden;
            transition: max-height 0.4s cubic-bezier(0.4, 0, 0.2, 1);
            background: #fafafa;
        }

        .drug-item.expanded .drug-details {
            max-height: 800px;
        }

        .drug-details-inner {
            padding: 20px;
            border-top: 1px solid #f1f5f9;
        }

        /* Dose Range Cards - Mobile First (stacked) */
        .dose-range-grid {
            display: grid;
            grid-template-columns: 1fr;
            gap: 8px;
            margin-bottom: 16px;
        }

        /* Mobile: Horizontal layout with left-aligned labels */
        .dose-card {
            background: white;
            border: 1px solid #e5e7eb;
            border-radius: 12px;
            padding: 12px 16px;
            display: flex;
            justify-content: space-between;
            align-items: center;
            text-align: left;
            transition: all 0.2s;
        }

        .dose-card:hover {
            border-color: #d1d5db;
        }

        .dose-card.low {
            border-left: 3px solid #22c55e;
        }

        .dose-card.standard {
            border-left: 3px solid #2563eb;
            background: #f8fafc;
        }

        .dose-card.high {
            border-left: 3px solid #f59e0b;
        }

        .dose-card-label {
            font-size: 10px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            margin-bottom: 0;
        }

        .dose-card.low .dose-card-label { color: #16a34a; }
        .dose-card.standard .dose-card-label { color: #2563eb; }
        .dose-card.high .dose-card-label { color: #d97706; }

        .dose-card-value {
            font-size: 20px;
            font-weight: 800;
            color: #0f172a;
            line-height: 1;
            margin-bottom: 0;
        }

        .dose-card-unit {
            font-size: 12px;
            color: #64748b;
        }

        .dose-card-right {
            display: flex;
            align-items: baseline;
            gap: 4px;
        }

        /* Hide per-kg on mobile to save space */
        .dose-card-perkg {
            display: none;
        }

        /* Additional indications */
        .additional-doses {
            margin-bottom: 16px;
        }

        .additional-dose-row {
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 12px 14px;
            background: white;
            border: 1px solid #e5e7eb;
            border-radius: 10px;
            margin-bottom: 8px;
        }

        .additional-dose-label {
            font-size: 13px;
            font-weight: 500;
            color: #475569;
        }

        .additional-dose-value {
            font-size: 16px;
            font-weight: 700;
            color: #2563eb;
        }

        .additional-dose-unit {
            font-size: 12px;
            color: #64748b;
            margin-left: 3px;
        }

        .drug-meta {
            display: flex;
            gap: 24px;
            margin-bottom: 14px;
            padding-bottom: 14px;
            border-bottom: 1px solid #e5e7eb;
        }

        .meta-item {
            display: flex;
            flex-direction: column;
            gap: 2px;
        }

        .meta-label {
            font-size: 10px;
            text-transform: uppercase;
            letter-spacing: 0.3px;
            color: #94a3b8;
            font-weight: 600;
        }

        .meta-value {
            font-size: 13px;
            color: #475569;
            font-weight: 500;
        }

        .drug-notes {
            font-size: 13px;
            color: #64748b;
            line-height: 1.6;
            margin-bottom: 12px;
        }

        .drug-reference {
            font-size: 11px;
            color: #94a3b8;
            display: flex;
            align-items: center;
            gap: 6px;
        }

        /* Warning */
        .drug-warning {
            background: #fef2f2;
            border: 1px solid #fecaca;
            border-radius: 10px;
            padding: 12px 14px;
            margin-top: 12px;
            display: flex;
            align-items: flex-start;
            gap: 10px;
        }

        .drug-warning svg {
            color: #dc2626;
            flex-shrink: 0;
            margin-top: 1px;
        }

        .drug-warning span {
            font-size: 12px;
            color: #b91c1c;
            line-height: 1.5;
        }

        /* Disclaimer */
        .disclaimer {
            text-align: center;
            margin-top: 48px;
            padding: 20px;
        }

        .disclaimer-text {
            font-size: 12px;
            color: #94a3b8;
            line-height: 1.6;
        }

        /* Footer - Homepage Style */
        .footer {
            padding: 32px 20px;
            border-top: 1px solid var(--gray-200);
            background: rgba(255,255,255,0.5);
            margin-top: auto;
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

        /* Desktop Styles - Scale up from mobile */
        @media (min-width: 768px) {
            /* Desktop nav */
            .nav {
                padding: 16px 32px;
            }

            .nav-inner {
                height: 64px;
                padding: 0 24px;
                border-radius: 20px;
            }

            .logo-icon svg {
                width: 42px;
                height: 15px;
            }

            .logo-text {
                font-size: 20px;
            }

            .nav-links {
                display: flex;
            }

            .mobile-menu-btn {
                display: none;
            }

            /* Desktop container - wider layout */
            .container {
                max-width: 900px;
                padding: 140px 32px 100px;
            }

            /* Larger weight input on desktop */
            .weight-section {
                margin-bottom: 56px;
            }

            .weight-input {
                font-size: 72px;
                width: 160px;
            }

            .weight-unit {
                font-size: 28px;
            }

            /* Desktop: Smaller weight buttons */
            .quick-weight-btn {
                padding: 10px 20px;
            }

            /* Desktop: More spacing */
            .color-legend {
                gap: 16px;
                padding: 20px 24px;
                margin-bottom: 40px;
            }

            .legend-item {
                font-size: 12px;
            }

            /* Desktop: More padding */
            .drug-header {
                padding: 18px 24px;
            }

            /* Desktop: 3-column grid with more spacing */
            .dose-range-grid {
                grid-template-columns: repeat(3, 1fr);
                gap: 12px;
                margin-bottom: 24px;
            }

            /* Desktop: Centered card layout with more breathing room */
            .dose-card {
                display: block;
                text-align: center;
                padding: 18px 16px;
                border-left: 1px solid #e5e7eb;
            }

            .dose-card.low {
                border-left: 1px solid #e5e7eb;
                border-top: 4px solid #22c55e;
            }

            .dose-card.standard {
                border-left: 1px solid #e5e7eb;
                border-top: 4px solid #2563eb;
            }

            .dose-card.high {
                border-left: 1px solid #e5e7eb;
                border-top: 4px solid #f59e0b;
            }

            .dose-card-label {
                margin-bottom: 8px;
                font-size: 11px;
            }

            .dose-card-value {
                font-size: 28px;
                margin-bottom: 6px;
            }

            .dose-card-unit {
                font-size: 14px;
            }

            /* Show per-kg on desktop */
            .dose-card-perkg {
                display: block;
                font-size: 11px;
                color: #94a3b8;
                margin-top: 8px;
            }

            .dose-card-right {
                display: block;
            }

            /* Desktop: Larger drug details padding */
            .drug-details-inner {
                padding: 24px;
            }

            .additional-dose-row {
                padding: 14px 18px;
            }

            .drug-meta {
                gap: 32px;
                margin-bottom: 18px;
                padding-bottom: 18px;
            }

            .meta-label {
                font-size: 11px;
            }

            .meta-value {
                font-size: 14px;
            }

            .drug-notes {
                font-size: 14px;
                line-height: 1.7;
            }

            .drug-reference {
                font-size: 12px;
            }

            /* Desktop footer */
            .footer {
                padding: 40px 32px;
            }

            .footer-inner {
                flex-direction: row;
                justify-content: space-between;
                text-align: left;
            }

            .footer-text {
                font-size: 14px;
            }

            .footer-links {
                gap: 32px;
            }

            .footer-link {
                font-size: 14px;
            }
        }

        @media (min-width: 1024px) {
            .nav { padding: 16px 40px; }
            .footer { padding: 48px 40px; }
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
        <nav class="nav" role="navigation" aria-label="Main navigation">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo" aria-label="GasConsult.ai home">
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
                    <div class="nav-dropdown">
                        <button class="nav-link nav-dropdown-toggle" onclick="toggleNavDropdown(event)">More ▼</button>
                        <div class="nav-dropdown-menu">
                            <a href="/hypotension" class="nav-dropdown-link">IOH Predictor</a>
                            <a href="/difficult-airway" class="nav-dropdown-link">Difficult Airway</a>
                            <a href="/informed-consent" class="nav-dropdown-link">Informed Consent</a>
                        </div>
                    </div>
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
            <a href="/difficult-airway" class="mobile-menu-link">Difficult Airway</a>
            <a href="/informed-consent" class="mobile-menu-link">Informed Consent</a>
        </div>

    <div class="container">
        <!-- Weight Input -->
        <div class="weight-section">
            <div class="weight-label">Patient Weight</div>
            <div class="weight-input-container">
                <input type="number" id="weight" class="weight-input" placeholder="70" value="70" min="1" max="500">
                <span class="weight-unit">kg</span>
            </div>
            <div class="quick-weights">
                <button class="quick-weight-btn" onclick="setWeight(50)">50 kg</button>
                <button class="quick-weight-btn active" onclick="setWeight(70)">70 kg</button>
                <button class="quick-weight-btn" onclick="setWeight(80)">80 kg</button>
                <button class="quick-weight-btn" onclick="setWeight(100)">100 kg</button>
                <button class="quick-weight-btn" onclick="setWeight(120)">120 kg</button>
            </div>
        </div>

        <!-- Color Legend -->
        <div class="color-legend">
            <div class="legend-item">
                <div class="legend-swatch hypnotic" style="background: var(--color-hypnotic);"></div>
                <span>Induction/Hypnotic</span>
            </div>
            <div class="legend-item">
                <div class="legend-swatch opioid"></div>
                <span>Opioid</span>
            </div>
            <div class="legend-item">
                <div class="legend-swatch paralytic"></div>
                <span>Paralytic</span>
            </div>
            <div class="legend-item">
                <div class="legend-swatch reversal"></div>
                <span>Reversal</span>
            </div>
            <div class="legend-item">
                <div class="legend-swatch vasopressor"></div>
                <span>Vasopressor</span>
            </div>
        </div>

        <!-- Drug Lists by Category -->
        <div id="drug-categories"></div>

        <!-- Disclaimer -->
        <div class="disclaimer">
            <div class="disclaimer-text">
                For reference only. Verify all doses and adjust for patient-specific factors. Color coding follows ISO 26825 standard.
            </div>
        </div>
    </div>

    <script>
        const drugData = {
            induction: {
                name: "Induction Agents",
                color: "hypnotic",
                drugs: [
                    {
                        name: "Propofol",
                        low: 1.5,
                        standard: 2,
                        high: 2.5,
                        unit: "mg",
                        additionalDoses: [
                            { label: "Elderly/Compromised", perKg: 1, unit: "mg" }
                        ],
                        onset: "30-45 sec",
                        duration: "5-10 min",
                        notes: "Reduce dose 30-50% in elderly, hypovolemic, or hemodynamically unstable. Lidocaine 20-40mg can reduce injection pain.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    },
                    {
                        name: "Etomidate",
                        low: 0.2,
                        standard: 0.3,
                        high: 0.4,
                        unit: "mg",
                        onset: "30-60 sec",
                        duration: "3-5 min",
                        notes: "Hemodynamically stable. Single-dose adrenal suppression is not clinically significant. Myoclonus common.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    },
                    {
                        name: "Ketamine",
                        low: 1,
                        standard: 1.5,
                        high: 2,
                        unit: "mg",
                        additionalDoses: [
                            { label: "IM Induction", perKg: 5, unit: "mg", range: "4-6 mg/kg" },
                            { label: "Analgesic (sub-dissociative)", perKg: 0.25, unit: "mg" }
                        ],
                        onset: "IV: 30-60 sec | IM: 3-4 min",
                        duration: "10-20 min",
                        notes: "Maintains airway reflexes and hemodynamics. Bronchodilator. Consider midazolam for emergence phenomena.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    },
                    {
                        name: "Midazolam",
                        low: 0.02,
                        standard: 0.03,
                        high: 0.05,
                        unit: "mg",
                        additionalDoses: [
                            { label: "Anxiolysis", fixed: 1, unit: "mg", range: "1-2 mg" }
                        ],
                        onset: "1-3 min",
                        duration: "30-60 min",
                        notes: "Marked synergy with opioids - reduce both when combined. Anterograde amnesia. Reduce dose in elderly.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    }
                ]
            },
            opioid: {
                name: "Opioids",
                color: "opioid",
                drugs: [
                    {
                        name: "Fentanyl",
                        low: 1,
                        standard: 2,
                        high: 3,
                        unit: "mcg",
                        additionalDoses: [
                            { label: "High-dose (cardiac)", perKg: 10, unit: "mcg", range: "5-15 mcg/kg" }
                        ],
                        onset: "1-2 min",
                        duration: "30-60 min",
                        notes: "Chest wall rigidity at high doses (>5 mcg/kg bolus). Context-sensitive half-time increases with prolonged infusion.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    },
                    {
                        name: "Sufentanil",
                        low: 0.1,
                        standard: 0.3,
                        high: 0.5,
                        unit: "mcg",
                        additionalDoses: [
                            { label: "Cardiac anesthesia", perKg: 5, unit: "mcg", range: "2-8 mcg/kg" }
                        ],
                        onset: "1-3 min",
                        duration: "20-45 min",
                        notes: "~10× potency of fentanyl. Excellent hemodynamic stability for cardiac surgery.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    },
                    {
                        name: "Remifentanil",
                        low: 0.5,
                        standard: 1,
                        high: 1.5,
                        unit: "mcg",
                        additionalDoses: [
                            { label: "Infusion", perKg: 0.15, unit: "mcg/kg/min", range: "0.05-0.2" }
                        ],
                        onset: "1-3 min",
                        duration: "3-5 min",
                        notes: "Ester hydrolysis - organ-independent. Rapid offset requires transition to longer-acting opioid before emergence.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    },
                    {
                        name: "Hydromorphone",
                        low: 0.01,
                        standard: 0.015,
                        high: 0.02,
                        unit: "mg",
                        onset: "5 min",
                        duration: "3-4 hours",
                        notes: "5-7× morphine potency. Less histamine release. No active metabolites - preferred in renal impairment.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    },
                    {
                        name: "Morphine",
                        low: 0.05,
                        standard: 0.1,
                        high: 0.15,
                        unit: "mg",
                        onset: "5-10 min",
                        duration: "3-4 hours",
                        notes: "Histamine release - hypotension, bronchospasm possible. Active metabolite M6G accumulates in renal failure.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    }
                ]
            },
            paralytic: {
                name: "Neuromuscular Blockers",
                color: "paralytic",
                drugs: [
                    {
                        name: "Succinylcholine",
                        low: 1,
                        standard: 1.5,
                        high: 2,
                        unit: "mg",
                        onset: "30-60 sec",
                        duration: "5-10 min",
                        notes: "Only depolarizing NMBA. Fasciculations precede paralysis. ~0.5 mEq/L K+ rise in normal patients.",
                        reference: "Barash Clinical Anesthesia, 8th Ed",
                        warning: "Contraindicated: hyperkalemia risk (burns >24h, crush injury, denervation, prolonged immobility), MH history, myopathies"
                    },
                    {
                        name: "Rocuronium",
                        low: 0.6,
                        standard: 0.9,
                        high: 1.2,
                        unit: "mg",
                        additionalDoses: [
                            { label: "Maintenance", perKg: 0.15, unit: "mg", range: "0.1-0.2 mg/kg" }
                        ],
                        onset: "60-90 sec (high dose: 45-60 sec)",
                        duration: "30-60 min (dose-dependent)",
                        notes: "Aminosteroid NMBA. Fully reversible with sugammadex at any depth. RSI alternative to succinylcholine.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    },
                    {
                        name: "Vecuronium",
                        low: 0.08,
                        standard: 0.1,
                        high: 0.15,
                        unit: "mg",
                        additionalDoses: [
                            { label: "Maintenance", perKg: 0.02, unit: "mg" }
                        ],
                        onset: "2-3 min",
                        duration: "25-40 min",
                        notes: "Aminosteroid. Hepatic elimination - prolonged in liver failure. Reversible with sugammadex.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    },
                    {
                        name: "Cisatracurium",
                        low: 0.1,
                        standard: 0.15,
                        high: 0.2,
                        unit: "mg",
                        additionalDoses: [
                            { label: "Maintenance", perKg: 0.03, unit: "mg" }
                        ],
                        onset: "2-3 min",
                        duration: "40-60 min",
                        notes: "Hoffman elimination (organ-independent). Ideal for hepatic/renal failure. No histamine release.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    }
                ]
            },
            reversal: {
                name: "Reversal Agents",
                color: "reversal",
                drugs: [
                    {
                        name: "Sugammadex",
                        low: 2,
                        standard: 4,
                        high: 16,
                        unit: "mg",
                        lowLabel: "Moderate Block",
                        standardLabel: "Deep Block",
                        highLabel: "Immediate",
                        onset: "1-3 min",
                        duration: "N/A",
                        notes: "Selective relaxant binding agent for rocuronium/vecuronium. May reduce hormonal contraceptive efficacy for 7 days.",
                        reference: "FDA Package Insert"
                    },
                    {
                        name: "Neostigmine",
                        low: 0.03,
                        standard: 0.05,
                        high: 0.07,
                        unit: "mg",
                        max: 5,
                        additionalDoses: [
                            { label: "With glycopyrrolate", note: "0.2 mg glyco per 1 mg neo" }
                        ],
                        onset: "5-10 min",
                        duration: "40-60 min",
                        notes: "Anticholinesterase - MUST co-administer anticholinergic. Only effective at TOF count ≥2.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    },
                    {
                        name: "Glycopyrrolate",
                        low: 0.006,
                        standard: 0.01,
                        high: 0.014,
                        unit: "mg",
                        additionalDoses: [
                            { label: "Antisialagogue", fixed: 0.2, unit: "mg" }
                        ],
                        onset: "1 min",
                        duration: "2-4 hours",
                        notes: "Quaternary ammonium - does not cross BBB. Less tachycardia than atropine.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    },
                    {
                        name: "Naloxone",
                        isFixed: true,
                        low: 0.02,
                        standard: 0.04,
                        high: 0.4,
                        unit: "mg",
                        lowLabel: "Micro-dose",
                        standardLabel: "Titrated",
                        highLabel: "Full Reversal",
                        onset: "1-2 min",
                        duration: "30-45 min",
                        notes: "Opioid antagonist. Duration shorter than most opioids - resedation possible. Titrate to avoid withdrawal.",
                        reference: "ACLS Guidelines 2020"
                    },
                    {
                        name: "Flumazenil",
                        isFixed: true,
                        low: 0.1,
                        standard: 0.2,
                        high: 0.5,
                        unit: "mg",
                        onset: "1-2 min",
                        duration: "45-90 min",
                        notes: "Benzodiazepine antagonist. May precipitate seizures in chronic benzo users. Max total 3-5 mg.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    }
                ]
            },
            vasopressor: {
                name: "Vasopressors & Inotropes",
                color: "vasopressor",
                drugs: [
                    {
                        name: "Phenylephrine",
                        isFixed: true,
                        low: 50,
                        standard: 100,
                        high: 200,
                        unit: "mcg",
                        onset: "Immediate",
                        duration: "5-10 min",
                        notes: "Pure α₁-agonist. Increases SVR without inotropy. Reflex bradycardia - useful when HR elevated.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    },
                    {
                        name: "Ephedrine",
                        isFixed: true,
                        low: 5,
                        standard: 10,
                        high: 25,
                        unit: "mg",
                        onset: "1-2 min",
                        duration: "10-15 min",
                        notes: "Mixed α/β agonist + indirect NE release. Tachyphylaxis with repeated doses.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    },
                    {
                        name: "Epinephrine",
                        isFixed: true,
                        low: 5,
                        standard: 10,
                        high: 20,
                        unit: "mcg",
                        additionalDoses: [
                            { label: "Cardiac Arrest", fixed: 1, unit: "mg", range: "1 mg q3-5 min" },
                            { label: "Anaphylaxis IM", fixed: 0.3, unit: "mg", range: "0.3-0.5 mg" }
                        ],
                        onset: "Immediate",
                        duration: "5-10 min",
                        notes: "α + β agonist. First-line for anaphylaxis and cardiac arrest. Low-dose: β effects; high-dose: α effects.",
                        reference: "ACLS Guidelines 2020"
                    },
                    {
                        name: "Vasopressin",
                        isFixed: true,
                        low: 1,
                        standard: 2,
                        high: 4,
                        unit: "U",
                        additionalDoses: [
                            { label: "Infusion (shock)", note: "0.01-0.04 U/min" },
                            { label: "Cardiac Arrest", fixed: 40, unit: "U" }
                        ],
                        onset: "Immediate",
                        duration: "10-20 min",
                        notes: "V₁ receptor agonist. Catecholamine-sparing. No tachycardia. May cause splanchnic ischemia.",
                        reference: "Surviving Sepsis Guidelines 2021"
                    },
                    {
                        name: "Atropine",
                        isFixed: true,
                        low: 0.5,
                        standard: 1,
                        high: 1.5,
                        unit: "mg",
                        onset: "1-2 min",
                        duration: "30-60 min",
                        notes: "Anticholinergic. Crosses BBB. Paradoxical bradycardia possible with doses <0.5 mg. Max 3 mg.",
                        reference: "ACLS Guidelines 2020"
                    }
                ]
            },
            other: {
                name: "Other Agents",
                color: "other",
                drugs: [
                    {
                        name: "Ondansetron",
                        isFixed: true,
                        low: 4,
                        standard: 4,
                        high: 8,
                        unit: "mg",
                        onset: "5-10 min",
                        duration: "4-8 hours",
                        notes: "5-HT₃ antagonist. Give at end of surgery. QTc prolongation - caution with other QT-prolonging drugs.",
                        reference: "Consensus Guidelines for PONV, 2020"
                    },
                    {
                        name: "Dexamethasone",
                        isFixed: true,
                        low: 4,
                        standard: 8,
                        high: 10,
                        unit: "mg",
                        onset: "1-2 hours (peak)",
                        duration: "24-72 hours",
                        notes: "Give at induction for PONV (needs time for effect). Transient hyperglycemia. Perineal burning if rapid IV.",
                        reference: "Consensus Guidelines for PONV, 2020"
                    },
                    {
                        name: "Tranexamic Acid",
                        low: 10,
                        standard: 20,
                        high: 30,
                        unit: "mg",
                        additionalDoses: [
                            { label: "Maintenance infusion", perKg: 1, unit: "mg/kg/hr", range: "1-5 mg/kg/hr" }
                        ],
                        onset: "5-15 min",
                        duration: "~3 hours",
                        notes: "Antifibrinolytic. Contraindicated in active thromboembolic disease. Seizure risk at very high doses.",
                        reference: "CRASH-2 Trial, Lancet 2010"
                    },
                    {
                        name: "Dexmedetomidine",
                        low: 0.5,
                        standard: 1,
                        high: 1,
                        unit: "mcg",
                        additionalDoses: [
                            { label: "Infusion", perKg: 0.5, unit: "mcg/kg/hr", range: "0.2-0.7" }
                        ],
                        onset: "5-10 min",
                        duration: "Infusion-dependent",
                        notes: "α₂-agonist. Sedation without respiratory depression. Bradycardia and hypotension common - load slowly.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    },
                    {
                        name: "Lidocaine IV",
                        low: 1,
                        standard: 1.5,
                        high: 2,
                        unit: "mg",
                        onset: "1-2 min",
                        duration: "10-20 min",
                        notes: "Blunts airway reflexes, reduces MAC. Max 4.5 mg/kg. CNS toxicity precedes cardiac.",
                        reference: "Barash Clinical Anesthesia, 8th Ed"
                    }
                ]
            }
        };

        let currentWeight = 70;

        document.addEventListener('DOMContentLoaded', () => {
            renderAllCategories();

            document.getElementById('weight').addEventListener('input', (e) => {
                const val = parseFloat(e.target.value);
                if (val && val > 0) {
                    currentWeight = val;
                    updateQuickWeightButtons();
                    updateAllDoses();
                }
            });
        });

        function setWeight(w) {
            currentWeight = w;
            document.getElementById('weight').value = w;
            updateQuickWeightButtons();
            updateAllDoses();
        }

        function updateQuickWeightButtons() {
            document.querySelectorAll('.quick-weight-btn').forEach(btn => {
                btn.classList.remove('active');
                if (parseInt(btn.textContent) === currentWeight) {
                    btn.classList.add('active');
                }
            });
        }

        function formatDose(value) {
            if (value >= 1000) return (value / 1000).toFixed(1);
            if (value < 0.1) return value.toFixed(3);
            if (value < 1) return value.toFixed(2);
            if (value < 10) return value.toFixed(1);
            return Math.round(value);
        }

        function formatUnit(value, unit) {
            if (unit === 'mg' && value >= 1000) return 'g';
            if (unit === 'mcg' && value >= 1000) return 'mg';
            return unit;
        }

        function toggleDrug(el) {
            el.closest('.drug-item').classList.toggle('expanded');
        }

        function updateAllDoses() {
            Object.entries(drugData).forEach(([catKey, category]) => {
                category.drugs.forEach(drug => {
                    const item = document.querySelector(`.drug-item[data-drug="${drug.name}"]`);
                    if (!item) return;

                    // Update header quick dose
                    const standardVal = drug.isFixed ? drug.standard : drug.standard * currentWeight;
                    item.querySelector('.drug-quick-value').textContent = formatDose(standardVal);
                    item.querySelector('.drug-quick-unit').textContent = formatUnit(standardVal, drug.unit);

                    // Update dose cards
                    const lowVal = drug.isFixed ? drug.low : drug.low * currentWeight;
                    const highVal = drug.isFixed ? drug.high : drug.high * currentWeight;

                    const lowCard = item.querySelector('.dose-card.low .dose-card-value');
                    const stdCard = item.querySelector('.dose-card.standard .dose-card-value');
                    const highCard = item.querySelector('.dose-card.high .dose-card-value');

                    if (lowCard) lowCard.textContent = formatDose(lowVal);
                    if (stdCard) stdCard.textContent = formatDose(standardVal);
                    if (highCard) highCard.textContent = formatDose(highVal);

                    // Update units
                    item.querySelectorAll('.dose-card.low .dose-card-unit').forEach(el => el.textContent = formatUnit(lowVal, drug.unit));
                    item.querySelectorAll('.dose-card.standard .dose-card-unit').forEach(el => el.textContent = formatUnit(standardVal, drug.unit));
                    item.querySelectorAll('.dose-card.high .dose-card-unit').forEach(el => el.textContent = formatUnit(highVal, drug.unit));

                    // Update additional doses
                    if (drug.additionalDoses) {
                        drug.additionalDoses.forEach((dose, i) => {
                            const row = item.querySelectorAll('.additional-dose-row')[i];
                            if (row && dose.perKg) {
                                const val = dose.perKg * currentWeight;
                                row.querySelector('.additional-dose-value').textContent = formatDose(val);
                                row.querySelector('.additional-dose-unit').textContent = formatUnit(val, dose.unit);
                            }
                        });
                    }
                });
            });
        }

        function renderAllCategories() {
            const container = document.getElementById('drug-categories');

            container.innerHTML = Object.entries(drugData).map(([catKey, category]) => `
                <div class="category-group">
                    <div class="category-header">
                        <div class="category-swatch ${catKey}"></div>
                        <div class="category-title">${category.name}</div>
                    </div>
                    <div class="drug-list">
                        ${category.drugs.map(drug => {
                            const standardVal = drug.isFixed ? drug.standard : drug.standard * currentWeight;
                            const lowVal = drug.isFixed ? drug.low : drug.low * currentWeight;
                            const highVal = drug.isFixed ? drug.high : drug.high * currentWeight;

                            return `
                                <div class="drug-item ${catKey}" data-drug="${drug.name}" data-category="${catKey}">
                                    <div class="drug-header" onclick="toggleDrug(this)">
                                        <span class="drug-name">${drug.name}</span>
                                        <div class="drug-quick-dose">
                                            <div class="drug-standard-dose">
                                                <span class="drug-quick-value">${formatDose(standardVal)}</span>
                                                <span class="drug-quick-unit">${formatUnit(standardVal, drug.unit)}</span>
                                            </div>
                                            <svg class="drug-chevron" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                                <polyline points="6 9 12 15 18 9"></polyline>
                                            </svg>
                                        </div>
                                    </div>
                                    <div class="drug-details">
                                        <div class="drug-details-inner">
                                            <div class="dose-range-grid">
                                                <div class="dose-card low">
                                                    <div class="dose-card-label">${drug.lowLabel || 'Low'}</div>
                                                    <div class="dose-card-value">${formatDose(lowVal)}</div>
                                                    <div class="dose-card-unit">${formatUnit(lowVal, drug.unit)}</div>
                                                    ${!drug.isFixed ? `<div class="dose-card-perkg">${drug.low} ${drug.unit}/kg</div>` : ''}
                                                </div>
                                                <div class="dose-card standard">
                                                    <div class="dose-card-label">${drug.standardLabel || 'Standard'}</div>
                                                    <div class="dose-card-value">${formatDose(standardVal)}</div>
                                                    <div class="dose-card-unit">${formatUnit(standardVal, drug.unit)}</div>
                                                    ${!drug.isFixed ? `<div class="dose-card-perkg">${drug.standard} ${drug.unit}/kg</div>` : ''}
                                                </div>
                                                <div class="dose-card high">
                                                    <div class="dose-card-label">${drug.highLabel || 'High'}</div>
                                                    <div class="dose-card-value">${formatDose(highVal)}</div>
                                                    <div class="dose-card-unit">${formatUnit(highVal, drug.unit)}</div>
                                                    ${!drug.isFixed ? `<div class="dose-card-perkg">${drug.high} ${drug.unit}/kg</div>` : ''}
                                                </div>
                                            </div>

                                            ${drug.additionalDoses ? `
                                                <div class="additional-doses">
                                                    ${drug.additionalDoses.map(dose => {
                                                        if (dose.note) {
                                                            return `
                                                                <div class="additional-dose-row">
                                                                    <span class="additional-dose-label">${dose.label}</span>
                                                                    <span style="font-size: 13px; color: #64748b;">${dose.note}</span>
                                                                </div>
                                                            `;
                                                        }
                                                        const val = dose.fixed || dose.perKg * currentWeight;
                                                        return `
                                                            <div class="additional-dose-row">
                                                                <span class="additional-dose-label">${dose.label}${dose.range ? ` (${dose.range})` : ''}</span>
                                                                <div>
                                                                    <span class="additional-dose-value">${formatDose(val)}</span>
                                                                    <span class="additional-dose-unit">${formatUnit(val, dose.unit)}</span>
                                                                </div>
                                                            </div>
                                                        `;
                                                    }).join('')}
                                                </div>
                                            ` : ''}

                                            <div class="drug-meta">
                                                <div class="meta-item">
                                                    <span class="meta-label">Onset</span>
                                                    <span class="meta-value">${drug.onset}</span>
                                                </div>
                                                <div class="meta-item">
                                                    <span class="meta-label">Duration</span>
                                                    <span class="meta-value">${drug.duration}</span>
                                                </div>
                                            </div>

                                            <div class="drug-notes">${drug.notes}</div>

                                            <div class="drug-reference">
                                                <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                                    <path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path>
                                                    <path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path>
                                                </svg>
                                                ${drug.reference}
                                            </div>

                                            ${drug.warning ? `
                                                <div class="drug-warning">
                                                    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5">
                                                        <path d="M10.29 3.86L1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z"/>
                                                        <line x1="12" y1="9" x2="12" y2="13"/>
                                                        <line x1="12" y1="17" x2="12.01" y2="17"/>
                                                    </svg>
                                                    <span>${drug.warning}</span>
                                                </div>
                                            ` : ''}
                                        </div>
                                    </div>
                                </div>
                            `;
                        }).join('')}
                    </div>
                </div>
            `).join('');
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

        // Toggle nav dropdown
        function toggleNavDropdown(e) {
            e.preventDefault();
            e.stopPropagation();
            const menu = e.target.nextElementSibling;
            if (menu) {
                menu.classList.toggle('show');
            }
        }

        // Close dropdown when clicking outside
        document.addEventListener('click', function() {
            document.querySelectorAll('.nav-dropdown-menu').forEach(m => m.classList.remove('show'));
        });
    </script>

    <footer class="footer">
        <div class="footer-inner">
            <span class="footer-text">© 2025 GasConsult.ai</span>
            <div class="footer-links">
                <a href="/privacy" class="footer-link">Privacy</a>
                <a href="/terms" class="footer-link">Terms</a>
                <a href="mailto:contact@gasconsult.ai" class="footer-link">Contact</a>
            </div>
        </div>
    </footer>
    </div>
</body>
</html>
"""

CALCULATORS_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <!-- Google Analytics 4 -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=G-01NZYD1DPP"></script>
    <script>
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-01NZYD1DPP');
    </script>

    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Clinical Calculators - GasConsult.ai</title>
    <meta name="description" content="Evidence-based anesthesiology clinical calculators: IBW, BSA, MABL, QTc, PONV, RCRI, and more with validated formulas and references.">

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

        .nav-dropdown:has(.nav-dropdown-link.active) .nav-dropdown-toggle {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .nav-dropdown {
            position: relative;
            display: inline-block;
        }

        .nav-dropdown-toggle {
            cursor: pointer;
            background: none;
            border: none;
            font-family: inherit;
        }

        .nav-dropdown-menu {
            display: none;
            position: absolute;
            top: 100%;
            right: 0;
            background: white;
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            min-width: 200px;
            margin-top: 4px;
            z-index: 1000;
            overflow: hidden;
        }

        .nav-dropdown-menu.show {
            display: block;
        }

        .nav-dropdown-link {
            display: block;
            padding: 12px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .nav-dropdown-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
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

        /* Category Filters */
        .category-filters {
            display: flex;
            flex-wrap: wrap;
            gap: 12px;
            margin-bottom: 40px;
            padding: 0 4px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.25s forwards;
            opacity: 0;
        }

        .filter-btn {
            display: flex;
            align-items: center;
            gap: 10px;
            padding: 14px 20px;
            background: rgba(255, 255, 255, 0.7);
            backdrop-filter: blur(20px) saturate(180%);
            -webkit-backdrop-filter: blur(20px) saturate(180%);
            border: 1.5px solid var(--gray-200);
            border-radius: 14px;
            font-family: inherit;
            font-size: 14px;
            font-weight: 600;
            color: var(--gray-700);
            cursor: pointer;
            transition: all 0.3s cubic-bezier(0.4,0,0.2,1);
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 2px 8px rgba(0,0,0,0.03);
        }

        .filter-btn:hover {
            transform: translateY(-2px);
            box-shadow: 0 2px 4px rgba(0,0,0,0.04), 0 8px 16px rgba(0,0,0,0.06);
            border-color: var(--blue-300);
            background: rgba(255, 255, 255, 0.85);
        }

        .filter-btn.active {
            background: linear-gradient(135deg, var(--blue-600) 0%, var(--blue-500) 100%);
            border-color: var(--blue-600);
            color: white;
            box-shadow: 0 2px 4px rgba(37, 99, 235, 0.2), 0 8px 16px rgba(37, 99, 235, 0.15);
        }

        .filter-btn svg {
            flex-shrink: 0;
            opacity: 0.8;
        }

        .filter-btn.active svg {
            opacity: 1;
        }

        .filter-btn span:first-of-type {
            white-space: nowrap;
        }

        .filter-count {
            display: inline-flex;
            align-items: center;
            justify-content: center;
            min-width: 24px;
            height: 24px;
            padding: 0 8px;
            background: rgba(0, 0, 0, 0.08);
            border-radius: 12px;
            font-size: 12px;
            font-weight: 700;
            margin-left: auto;
        }

        .filter-btn.active .filter-count {
            background: rgba(255, 255, 255, 0.25);
        }

        @media (max-width: 767px) {
            .filter-btn {
                flex: 1 1 calc(50% - 6px);
                min-width: 0;
                padding: 12px 16px;
                font-size: 13px;
            }

            .filter-btn span:first-of-type {
                overflow: hidden;
                text-overflow: ellipsis;
            }

            .category-filters {
                gap: 8px;
            }
        }

        @media (min-width: 768px) and (max-width: 1199px) {
            .filter-btn {
                flex: 1 1 calc(33.333% - 8px);
            }
        }

        /* Calculator Grid */
        .calc-grid {
            display: grid;
            grid-template-columns: 1fr;
            gap: 24px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.35s forwards;
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

        /* Hide filtered calculators */
        .calc-card.hidden {
            display: none;
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

            .footer { padding: 40px 32px; }
            .footer-inner { flex-direction: row; justify-content: space-between; text-align: left; }
            .footer-text { font-size: 14px; }
            .footer-links { gap: 32px; }
            .footer-link { font-size: 14px; }
        }

        @media (min-width: 1024px) {
            .footer { padding: 48px 40px; }
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

        .footer-link:hover {
            color: var(--gray-700);
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
                    <a href="/crisis" class="nav-link">Crisis Protocols</a>
                    <div class="nav-dropdown">
                        <button class="nav-link nav-dropdown-toggle" onclick="toggleNavDropdown(event)">More ▼</button>
                        <div class="nav-dropdown-menu">
                            <a href="/hypotension" class="nav-dropdown-link">IOH Predictor</a>
                            <a href="/difficult-airway" class="nav-dropdown-link">Difficult Airway</a>
                            <a href="/informed-consent" class="nav-dropdown-link">Informed Consent</a>
                        </div>
                    </div>
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
            <a href="/difficult-airway" class="mobile-menu-link">Difficult Airway</a>
            <a href="/informed-consent" class="mobile-menu-link">Informed Consent</a>
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
                Evidence-based medical calculations organized by specialty. All results are instant - no submit buttons needed.
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
                <!-- Category Filters -->
                <div class="category-filters">
                    <button class="filter-btn active" onclick="filterCalculators('all')">
                        <svg width="18" height="18" fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                            <path stroke-linecap="round" stroke-linejoin="round" d="M4 6h16M4 12h16M4 18h16"></path>
                        </svg>
                        <span>All Calculators</span>
                        <span class="filter-count">10</span>
                    </button>
                    <button class="filter-btn" onclick="filterCalculators('dosing')">
                        <svg width="18" height="18" fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                            <path stroke-linecap="round" stroke-linejoin="round" d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 008 10.172V5L7 4z"></path>
                        </svg>
                        <span>Dosing & Pharmacology</span>
                        <span class="filter-count">3</span>
                    </button>
                    <button class="filter-btn" onclick="filterCalculators('risk')">
                        <svg width="18" height="18" fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                            <path stroke-linecap="round" stroke-linejoin="round" d="M9 12l2 2 4-4m5.618-4.016A11.955 11.955 0 0112 2.944a11.955 11.955 0 01-8.618 3.04A12.02 12.02 0 003 9c0 5.591 3.824 10.29 9 11.622 5.176-1.332 9-6.03 9-11.622 0-1.042-.133-2.052-.382-3.016z"></path>
                        </svg>
                        <span>Perioperative Risk</span>
                        <span class="filter-count">3</span>
                    </button>
                    <button class="filter-btn" onclick="filterCalculators('fluids')">
                        <svg width="18" height="18" fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                            <path stroke-linecap="round" stroke-linejoin="round" d="M7 16a4 4 0 01-.88-7.903A5 5 0 1115.9 6L16 6a5 5 0 011 9.9M9 19l3 3m0 0l3-3m-3 3V10"></path>
                        </svg>
                        <span>Hemodynamics & Fluids</span>
                        <span class="filter-count">3</span>
                    </button>
                    <button class="filter-btn" onclick="filterCalculators('cardiac')">
                        <svg width="18" height="18" fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                            <path stroke-linecap="round" stroke-linejoin="round" d="M4.318 6.318a4.5 4.5 0 000 6.364L12 20.364l7.682-7.682a4.5 4.5 0 00-6.364-6.364L12 7.636l-1.318-1.318a4.5 4.5 0 00-6.364 0z"></path>
                        </svg>
                        <span>Cardiac & Monitoring</span>
                        <span class="filter-count">1</span>
                    </button>
                </div>

                <!-- Calculator Grid -->
                <div class="calc-grid">

                    <!-- IDEAL BODY WEIGHT -->
                    <div class="calc-card" data-category="dosing">
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
                    <div class="calc-card" data-category="fluids">
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
                    <div class="calc-card" data-category="fluids">
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
                    <div class="calc-card" data-category="cardiac">
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
                    <div class="calc-card" data-category="fluids">
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
                    <div class="calc-card" data-category="risk">
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
                    <div class="calc-card" data-category="risk">
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
                    <div class="calc-card" data-category="risk">
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
                    <div class="calc-card large" data-category="dosing">
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
                            <div class="result-text">Use caution with opioid conversions - consider patient factors and start at lower doses.</div>
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
                    <div class="calc-card" data-category="dosing">
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
                <span class="footer-text">© 2025 GasConsult.ai</span>
                <div class="footer-links">
                    <a href="/privacy" class="footer-link">Privacy</a>
                    <a href="/terms" class="footer-link">Terms</a>
                    <a href="mailto:contact@gasconsult.ai" class="footer-link">Contact</a>
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

        // Toggle Navigation Dropdown
        function toggleNavDropdown(event) {
            event.stopPropagation();
            const dropdown = event.target.closest('.nav-dropdown');
            const menu = dropdown.querySelector('.nav-dropdown-menu');
            menu.classList.toggle('show');
        }

        // Close dropdown when clicking outside
        document.addEventListener('click', () => {
            document.querySelectorAll('.nav-dropdown-menu.show').forEach(menu => {
                menu.classList.remove('show');
            });
        });

        // Filter Calculators by Category
        function filterCalculators(category) {
            const cards = document.querySelectorAll('.calc-card');
            const buttons = document.querySelectorAll('.filter-btn');

            // Update button active states
            buttons.forEach(btn => btn.classList.remove('active'));
            event.target.closest('.filter-btn').classList.add('active');

            // Filter cards with smooth transition
            cards.forEach(card => {
                if (category === 'all' || card.dataset.category === category) {
                    card.classList.remove('hidden');
                    // Smooth fade-in
                    card.style.opacity = '0';
                    card.style.transform = 'translateY(10px)';
                    setTimeout(() => {
                        card.style.transition = 'opacity 0.3s ease, transform 0.3s ease';
                        card.style.opacity = '1';
                        card.style.transform = 'translateY(0)';
                    }, 10);
                } else {
                    card.classList.add('hidden');
                }
            });
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
    <!-- Google Analytics 4 -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=G-01NZYD1DPP"></script>
    <script>
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-01NZYD1DPP');
    </script>

    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>IOH Predictor - GasConsult.ai</title>

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
            --orange-50: #FFF7ED;
            --orange-500: #F97316;
            --red-50: #FEF2F2;
            --red-500: #EF4444;
            --red-600: #DC2626;
            --green-50: #ECFDF5;
            --green-500: #10B981;
            --amber-50: #FFFBEB;
            --amber-500: #F59E0B;
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

        /* Glass Background */
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

        @keyframes fade-up {
            from { opacity: 0; transform: translateY(20px); }
            to { opacity: 1; transform: translateY(0); }
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

        .nav-dropdown {
            position: relative;
        }

        .nav-dropdown-toggle {
            cursor: pointer;
            background: none;
            border: none;
            font-family: inherit;
        }

        .nav-dropdown-menu {
            display: none;
            position: absolute;
            top: 100%;
            right: 0;
            background: white;
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            min-width: 200px;
            margin-top: 4px;
            z-index: 1000;
            overflow: hidden;
        }

        .nav-dropdown-menu.show { display: block; }

        .nav-dropdown-link {
            display: block;
            padding: 12px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .nav-dropdown-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
        }

        .nav-dropdown-link.active {
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
        }

        .mobile-menu-btn span {
            display: block;
            width: 22px;
            height: 2px;
            background: var(--gray-700);
            border-radius: 1px;
            transition: all 0.3s ease;
        }

        .mobile-menu {
            display: none;
            position: fixed;
            top: 80px;
            left: 16px;
            right: 16px;
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 16px;
            padding: 8px;
            box-shadow: 0 4px 16px rgba(0,0,0,0.08);
            z-index: 99;
            flex-direction: column;
            gap: 4px;
        }

        .mobile-menu.active { display: flex; }

        .mobile-menu-link {
            padding: 14px 16px;
            font-size: 15px;
            font-weight: 500;
            color: var(--gray-700);
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        /* Main Content */
        .main {
            flex: 1;
            padding: 100px 20px 60px;
            max-width: 1000px;
            margin: 0 auto;
            width: 100%;
        }

        /* Hero */
        .hero {
            text-align: center;
            margin-bottom: 48px;
        }

        .hero-badge {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            background: linear-gradient(135deg, var(--orange-50) 0%, #FEE2E2 100%);
            border: 1px solid rgba(249, 115, 22, 0.2);
            border-radius: 100px;
            padding: 8px 16px 8px 12px;
            margin-bottom: 20px;
            box-shadow: 0 2px 8px rgba(249, 115, 22, 0.08);
        }

        .badge-dot {
            width: 8px;
            height: 8px;
            background: var(--orange-500);
            border-radius: 50%;
        }

        .badge-text {
            font-size: 12px;
            font-weight: 600;
            color: #EA580C;
        }

        .hero-title {
            font-size: 42px;
            font-weight: 800;
            line-height: 1.1;
            letter-spacing: -1.5px;
            color: var(--gray-900);
            margin-bottom: 16px;
        }

        .hero-title .gradient {
            background: linear-gradient(135deg, var(--orange-500) 0%, var(--red-600) 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
        }

        .hero-subtitle {
            font-size: 17px;
            font-weight: 400;
            line-height: 1.6;
            color: var(--gray-600);
            max-width: 680px;
            margin: 0 auto;
        }

        /* Glass Cards */
        .glass-card {
            background: rgba(255, 255, 255, 0.8);
            backdrop-filter: blur(20px) saturate(180%);
            border: 1px solid rgba(255, 255, 255, 0.8);
            border-radius: 20px;
            padding: 32px;
            margin-bottom: 24px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 12px 48px rgba(0,0,0,0.03);
            animation: fade-up 0.6s ease forwards;
            opacity: 0;
        }

        .section-title {
            font-size: 24px;
            font-weight: 700;
            letter-spacing: -0.5px;
            color: var(--gray-900);
            margin-bottom: 12px;
        }

        .section-subtitle {
            font-size: 15px;
            line-height: 1.7;
            color: var(--gray-600);
            margin-bottom: 28px;
        }

        /* Form Styles */
        .form-group {
            margin-bottom: 20px;
        }

        .form-label {
            display: block;
            font-size: 14px;
            font-weight: 600;
            color: var(--gray-700);
            margin-bottom: 8px;
        }

        .required {
            color: var(--red-600);
            margin-left: 2px;
        }

        .form-input {
            width: 100%;
            padding: 12px 16px;
            border: 1px solid var(--gray-300);
            border-radius: 12px;
            font-size: 14px;
            font-family: inherit;
            transition: all 0.2s ease;
            background: white;
        }

        .form-input:focus {
            outline: none;
            border-color: var(--blue-500);
            box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.1);
        }

        .input-row {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 16px;
        }

        .section-divider {
            height: 1px;
            background: var(--gray-200);
            margin: 24px 0;
        }

        .submit-btn {
            width: 100%;
            padding: 16px;
            background: linear-gradient(135deg, var(--blue-600) 0%, var(--blue-700) 100%);
            color: white;
            font-size: 16px;
            font-weight: 600;
            border: none;
            border-radius: 12px;
            cursor: pointer;
            transition: all 0.3s ease;
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.2);
        }

        .submit-btn:hover {
            transform: translateY(-2px);
            box-shadow: 0 8px 20px rgba(37, 99, 235, 0.3);
        }

        .submit-btn:active {
            transform: translateY(0);
        }

        /* Risk Cards */
        .risk-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 16px;
            margin: 24px 0;
        }

        .risk-card {
            border-radius: 16px;
            padding: 24px;
            text-align: center;
            border: 2px solid;
            transition: all 0.3s ease;
        }

        .risk-card:hover {
            transform: translateY(-4px);
            box-shadow: 0 12px 24px rgba(0,0,0,0.1);
        }

        .risk-card.low {
            background: var(--green-50);
            border-color: var(--green-500);
        }

        .risk-card.moderate {
            background: var(--amber-50);
            border-color: var(--amber-500);
        }

        .risk-card.high {
            background: var(--red-50);
            border-color: var(--red-500);
        }

        .risk-label {
            font-size: 13px;
            font-weight: 600;
            color: var(--gray-600);
            margin-bottom: 8px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        .risk-value {
            font-size: 40px;
            font-weight: 900;
            margin-bottom: 8px;
            line-height: 1;
        }

        .risk-card.low .risk-value { color: var(--green-500); }
        .risk-card.moderate .risk-value { color: var(--amber-500); }
        .risk-card.high .risk-value { color: var(--red-500); }

        .risk-text {
            font-size: 13px;
            font-weight: 600;
            color: var(--gray-700);
        }

        /* Factor/Intervention Cards */
        .factor-card {
            background: var(--blue-50);
            border-left: 4px solid var(--blue-600);
            padding: 20px;
            margin-bottom: 12px;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .factor-card:hover {
            transform: translateX(4px);
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.1);
        }

        .intervention-card {
            background: var(--green-50);
            border-left: 4px solid var(--green-500);
            padding: 20px;
            margin-bottom: 12px;
            border-radius: 12px;
            transition: all 0.2s ease;
        }

        .intervention-card:hover {
            transform: translateX(4px);
            box-shadow: 0 4px 12px rgba(16, 185, 129, 0.1);
        }

        .card-title {
            font-size: 15px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 6px;
        }

        .card-desc {
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-700);
        }

        /* Model Evidence Section */
        .evidence-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 20px;
            margin-top: 24px;
        }

        .evidence-card {
            background: var(--blue-50);
            border: 1px solid var(--blue-200);
            border-radius: 16px;
            padding: 24px;
            transition: all 0.3s ease;
        }

        .evidence-card:hover {
            transform: translateY(-4px);
            box-shadow: 0 8px 20px rgba(37, 99, 235, 0.15);
        }

        .evidence-number {
            width: 40px;
            height: 40px;
            background: var(--blue-600);
            color: white;
            border-radius: 10px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 18px;
            font-weight: 800;
            margin-bottom: 16px;
        }

        .evidence-title {
            font-size: 16px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 8px;
        }

        .evidence-text {
            font-size: 14px;
            line-height: 1.7;
            color: var(--gray-700);
        }

        .evidence-text strong {
            color: var(--blue-700);
            font-weight: 600;
        }

        /* Footer */
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
            gap: 16px;
            text-align: center;
        }

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

        .footer-link:hover {
            color: var(--gray-700);
        }

        .back-btn {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            padding: 12px 24px;
            background: var(--gray-200);
            color: var(--gray-900);
            font-size: 14px;
            font-weight: 600;
            text-decoration: none;
            border-radius: 12px;
            transition: all 0.2s ease;
            margin-top: 24px;
        }

        .back-btn:hover {
            background: var(--gray-300);
            transform: translateX(-2px);
        }

        /* Responsive */
        @media (min-width: 768px) {
            .nav { padding: 16px 32px; }
            .nav-inner { height: 64px; padding: 0 24px; border-radius: 20px; }
            .logo-icon svg { width: 42px; height: 15px; }
            .logo-text { font-size: 20px; }
            .nav-links { display: flex; }
            .mobile-menu-btn { display: none; }
            .main { padding: 120px 32px 80px; }
            .hero-title { font-size: 52px; }
            .footer { padding: 40px 32px; }
            .footer-inner { flex-direction: row; justify-content: space-between; text-align: left; }
        }

        @media (min-width: 1024px) {
            .nav { padding: 16px 40px; }
            .main { padding: 140px 40px 100px; }
            .footer { padding: 48px 40px; }
        }

        /* Disclaimer Banner */
        .disclaimer {
            background: linear-gradient(135deg, rgba(239, 68, 68, 0.08) 0%, rgba(220, 38, 38, 0.05) 100%);
            border: 2px solid rgba(239, 68, 68, 0.2);
            border-radius: 16px;
            padding: 20px 24px;
            margin-bottom: 32px;
            display: flex;
            align-items: start;
            gap: 16px;
        }

        .disclaimer-icon {
            flex-shrink: 0;
            width: 24px;
            height: 24px;
            color: var(--red-600);
        }

        .disclaimer-text {
            font-size: 14px;
            line-height: 1.7;
            color: var(--gray-800);
        }

        .disclaimer-text strong {
            font-weight: 700;
            color: var(--red-600);
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
                    <a href="/calculators" class="nav-link">Clinical Calculators</a>
                    <a href="/crisis" class="nav-link">Crisis Protocols</a>
                    <div class="nav-dropdown">
                        <button class="nav-link nav-dropdown-toggle" onclick="toggleNavDropdown(event)">More ▼</button>
                        <div class="nav-dropdown-menu">
                            <a href="/hypotension" class="nav-dropdown-link active">IOH Predictor</a>
                            <a href="/difficult-airway" class="nav-dropdown-link">Difficult Airway</a>
                            <a href="/informed-consent" class="nav-dropdown-link">Informed Consent</a>
                        </div>
                    </div>
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
            <a href="/difficult-airway" class="mobile-menu-link">Difficult Airway</a>
            <a href="/informed-consent" class="mobile-menu-link">Informed Consent</a>
        </div>

        <!-- Main Content -->
        <main class="main">
            <!-- Hero -->
            <div class="hero">
                <div class="hero-badge">
                    <div class="badge-dot"></div>
                    <span class="badge-text">Machine Learning</span>
                </div>
                <h1 class="hero-title">
                    Intraoperative <span class="gradient">Hypotension</span> Predictor
                </h1>
                <p class="hero-subtitle">
                    Evidence-based ML model that predicts IOH risk using 14 clinical features. Trained on 10,000 synthetic cases based on published research.
                </p>
            </div>

            <!-- Disclaimer -->
            <div class="disclaimer">
                <svg class="disclaimer-icon" fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                    <path stroke-linecap="round" stroke-linejoin="round" d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z"></path>
                </svg>
                <div class="disclaimer-text">
                    <strong>Educational Tool Only:</strong> This ML model is designed for learning about IOH risk prediction. It is NOT intended for clinical decision-making or actual patient care. Always consult evidence-based guidelines and clinical judgment for real patient management.
                </div>
            </div>

            {% if not prediction %}
            <!-- Prediction Form -->
            <div class="glass-card">
                <h2 class="section-title">Patient & Hemodynamic Parameters</h2>
                <p class="section-subtitle">
                    Enter current intraoperative data to estimate IOH risk using our RandomForest classifier.
                </p>

                <form method="POST" action="/hypotension">
                    <input type="hidden" name="csrf_token" value="{{ csrf_token() }}"/>

                    <!-- Patient Demographics -->
                    <h3 style="font-size: 18px; font-weight: 700; color: var(--gray-900); margin-bottom: 16px;">Patient Information</h3>

                    <div class="input-row">
                        <div class="form-group">
                            <label class="form-label">Age<span class="required">*</span></label>
                            <input type="number" name="age" class="form-input" placeholder="65" required min="18" max="120">
                        </div>
                        <div class="form-group">
                            <label class="form-label">Sex<span class="required">*</span></label>
                            <select name="sex" class="form-input" required>
                                <option value="">Select...</option>
                                <option value="male">Male</option>
                                <option value="female">Female</option>
                            </select>
                        </div>
                    </div>

                    <div class="input-row">
                        <div class="form-group">
                            <label class="form-label">Weight (kg)<span class="required">*</span></label>
                            <input type="number" name="weight" class="form-input" placeholder="70" required min="30" max="300" step="0.1">
                        </div>
                        <div class="form-group">
                            <label class="form-label">Height (cm)<span class="required">*</span></label>
                            <input type="number" name="height" class="form-input" placeholder="170" required min="100" max="250">
                        </div>
                    </div>

                    <div class="form-group">
                        <label class="form-label">ASA Class<span class="required">*</span></label>
                        <select name="asa" class="form-input" required>
                            <option value="">Select...</option>
                            <option value="1">ASA I - Healthy</option>
                            <option value="2">ASA II - Mild systemic disease</option>
                            <option value="3">ASA III - Severe systemic disease</option>
                            <option value="4">ASA IV - Life-threatening disease</option>
                        </select>
                    </div>

                    <div class="section-divider"></div>

                    <!-- Hemodynamic Parameters -->
                    <h3 style="font-size: 18px; font-weight: 700; color: var(--gray-900); margin-bottom: 16px;">Hemodynamic Parameters</h3>

                    <div class="input-row">
                        <div class="form-group">
                            <label class="form-label">Baseline MAP (mmHg)<span class="required">*</span></label>
                            <input type="number" name="baseline_map" class="form-input" placeholder="85" required min="40" max="150">
                        </div>
                        <div class="form-group">
                            <label class="form-label">Baseline HR (bpm)<span class="required">*</span></label>
                            <input type="number" name="baseline_hr" class="form-input" placeholder="75" required min="30" max="200">
                        </div>
                    </div>

                    <div class="input-row">
                        <div class="form-group">
                            <label class="form-label">Current MAP (mmHg)<span class="required">*</span></label>
                            <input type="number" name="current_map" class="form-input" placeholder="72" required min="40" max="150">
                        </div>
                        <div class="form-group">
                            <label class="form-label">MAP 5 min ago<span class="required">*</span></label>
                            <input type="number" name="map_5min" class="form-input" placeholder="76" required min="40" max="150">
                        </div>
                        <div class="form-group">
                            <label class="form-label">MAP 10 min ago<span class="required">*</span></label>
                            <input type="number" name="map_10min" class="form-input" placeholder="80" required min="40" max="150">
                        </div>
                    </div>

                    <div class="section-divider"></div>

                    <!-- Surgical & Anesthetic Factors -->
                    <h3 style="font-size: 18px; font-weight: 700; color: var(--gray-900); margin-bottom: 16px;">Surgical & Anesthetic Factors</h3>

                    <div class="input-row">
                        <div class="form-group">
                            <label class="form-label">Surgery Type<span class="required">*</span></label>
                            <select name="surgery_type" class="form-input" required>
                                <option value="">Select...</option>
                                <option value="minor">Minor Surgery</option>
                                <option value="moderate">Moderate Surgery</option>
                                <option value="major_abdominal">Major Abdominal</option>
                                <option value="cardiac">Cardiac Surgery</option>
                                <option value="vascular">Vascular Surgery</option>
                            </select>
                        </div>
                        <div class="form-group">
                            <label class="form-label">Surgery Duration (min)<span class="required">*</span></label>
                            <input type="number" name="surgery_duration" class="form-input" placeholder="120" required min="0" max="1000">
                        </div>
                    </div>

                    <div class="input-row">
                        <div class="form-group">
                            <label class="form-label">Induction Agent<span class="required">*</span></label>
                            <select name="induction_agent" class="form-input" required>
                                <option value="">Select...</option>
                                <option value="propofol">Propofol</option>
                                <option value="etomidate">Etomidate</option>
                                <option value="ketamine">Ketamine</option>
                            </select>
                        </div>
                        <div class="form-group">
                            <label class="form-label">Vasopressor in Use<span class="required">*</span></label>
                            <select name="vasopressor" class="form-input" required>
                                <option value="none">None</option>
                                <option value="phenylephrine">Phenylephrine</option>
                                <option value="ephedrine">Ephedrine</option>
                                <option value="norepinephrine">Norepinephrine</option>
                            </select>
                        </div>
                    </div>

                    <div class="form-group">
                        <label class="form-label">Emergency Surgery<span class="required">*</span></label>
                        <select name="emergency" class="form-input" required>
                            <option value="no">No</option>
                            <option value="yes">Yes</option>
                        </select>
                    </div>

                    <button type="submit" class="submit-btn">Calculate IOH Risk</button>
                </form>
            </div>

            {% else %}
            <!-- Prediction Results -->
            <div class="glass-card">
                <h2 class="section-title">Machine Learning Prediction Results</h2>
                <p class="section-subtitle">
                    RandomForest classifier analyzed 14 clinical features to estimate IOH probability.
                </p>

                <!-- Risk Probabilities -->
                <div class="risk-grid">
                    <div class="risk-card {{ prediction.risk_5min_class }}">
                        <div class="risk-label">5-Minute Risk</div>
                        <div class="risk-value">{{ prediction.prob_5min }}%</div>
                        <div class="risk-text">{{ prediction.risk_5min_text }}</div>
                    </div>
                    <div class="risk-card {{ prediction.risk_10min_class }}">
                        <div class="risk-label">10-Minute Risk</div>
                        <div class="risk-value">{{ prediction.prob_10min }}%</div>
                        <div class="risk-text">{{ prediction.risk_10min_text }}</div>
                    </div>
                    <div class="risk-card {{ prediction.risk_20min_class }}">
                        <div class="risk-label">20-Minute Risk</div>
                        <div class="risk-value">{{ prediction.prob_20min }}%</div>
                        <div class="risk-text">{{ prediction.risk_20min_text }}</div>
                    </div>
                </div>

                <!-- Risk Factors -->
                <h3 style="font-size: 20px; font-weight: 700; color: var(--gray-900); margin-top: 40px; margin-bottom: 20px;">Top Risk Factors Identified</h3>
                {% for factor in prediction.factors %}
                <div class="factor-card">
                    <div class="card-title">{{ factor.name }}</div>
                    <div class="card-desc">{{ factor.description }}</div>
                </div>
                {% endfor %}

                <!-- Interventions -->
                <h3 style="font-size: 20px; font-weight: 700; color: var(--gray-900); margin-top: 40px; margin-bottom: 20px;">Suggested Interventions</h3>
                {% for intervention in prediction.interventions %}
                <div class="intervention-card">
                    <div class="card-title">{{ intervention.name }}</div>
                    <div class="card-desc">{{ intervention.description }}</div>
                </div>
                {% endfor %}

                <a href="/hypotension" class="back-btn">
                    <svg width="16" height="16" fill="none" stroke="currentColor" stroke-width="2" viewBox="0 0 24 24">
                        <path stroke-linecap="round" stroke-linejoin="round" d="M10 19l-7-7m0 0l7-7m-7 7h18"></path>
                    </svg>
                    Calculate Another
                </a>
            </div>
            {% endif %}

            <!-- Model Evidence & Transparency -->
            <div class="glass-card">
                <h2 class="section-title">How the ML Model Works</h2>
                <p class="section-subtitle">
                    Our RandomForest classifier provides transparent, evidence-based IOH risk predictions. Here's how it works:
                </p>

                <div class="evidence-grid">
                    <div class="evidence-card">
                        <div class="evidence-number">1</div>
                        <h3 class="evidence-title">Training Dataset</h3>
                        <p class="evidence-text">
                            Trained on <strong>10,000 synthetic cases</strong> generated from published clinical research on IOH risk factors. Training data reflects real-world IOH prevalence and risk patterns from literature.
                        </p>
                    </div>

                    <div class="evidence-card">
                        <div class="evidence-number">2</div>
                        <h3 class="evidence-title">14 Clinical Features</h3>
                        <p class="evidence-text">
                            Model analyzes <strong>demographics</strong> (age, sex, BMI, ASA), <strong>hemodynamics</strong> (MAP trends, HR), and <strong>surgical factors</strong> (type, duration, induction agent, emergency status).
                        </p>
                    </div>

                    <div class="evidence-card">
                        <div class="evidence-number">3</div>
                        <h3 class="evidence-title">RandomForest Algorithm</h3>
                        <p class="evidence-text">
                            Uses <strong>200 decision trees</strong> with balanced class weights and feature scaling. Test accuracy: ~85%. Most important features: current MAP, MAP trends, age, and ASA class.
                        </p>
                    </div>

                    <div class="evidence-card">
                        <div class="evidence-number">4</div>
                        <h3 class="evidence-title">Evidence-Based</h3>
                        <p class="evidence-text">
                            Risk factors weighted based on <strong>published literature</strong>: Sessler 2015 (MAP thresholds), Hatib 2018 (ML prediction), Walsh 2013 (organ injury), Reich 2005 (anesthetic effects).
                        </p>
                    </div>

                    <div class="evidence-card">
                        <div class="evidence-number">5</div>
                        <h3 class="evidence-title">Interpretable Outputs</h3>
                        <p class="evidence-text">
                            Provides <strong>3 time-windowed predictions</strong> (5/10/20 min) with decay modeling. Identifies top risk factors and suggests evidence-based interventions for clinical context.
                        </p>
                    </div>

                    <div class="evidence-card">
                        <div class="evidence-number">6</div>
                        <h3 class="evidence-title">Continuous Improvement</h3>
                        <p class="evidence-text">
                            Model can be <strong>retrained and updated</strong> as new clinical evidence emerges. Currently optimized for general surgical cases. Not validated on real patient data.
                        </p>
                    </div>
                </div>
            </div>
        </main>

        <!-- Footer -->
        <footer class="footer">
            <div class="footer-inner">
                <span class="footer-text">© 2025 GasConsult.ai</span>
                <div class="footer-links">
                    <a href="/privacy" class="footer-link">Privacy</a>
                    <a href="/terms" class="footer-link">Terms</a>
                    <a href="mailto:contact@gasconsult.ai" class="footer-link">Contact</a>
                </div>
            </div>
        </footer>
    </div>

    <script>
        // Mobile menu toggle
        function toggleMobileMenu() {
            const menu = document.getElementById('mobileMenu');
            const btn = document.querySelector('.mobile-menu-btn');
            if (menu && btn) {
                menu.classList.toggle('active');
                btn.classList.toggle('active');
            }
        }

        // Dropdown toggle
        function toggleNavDropdown(event) {
            event.stopPropagation();
            const menu = event.target.nextElementSibling;
            menu.classList.toggle('show');

            // Close dropdown when clicking outside
            document.addEventListener('click', function closeDropdown(e) {
                if (!event.target.contains(e.target)) {
                    menu.classList.remove('show');
                    document.removeEventListener('click', closeDropdown);
                }
            });
        }
    </script>
</body>
</html>
"""

INFORMED_CONSENT_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <!-- Google Analytics 4 -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=G-01NZYD1DPP"></script>
    <script>
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-01NZYD1DPP');
    </script>

    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Informed Consent for Anesthesia - GasConsult.ai</title>
    <meta name="description" content="Evidence-based informed consent guide for anesthesia procedures with comprehensive risk discussion.">

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
            --green-500: #22C55E;
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

        .nav-dropdown {
            position: relative;
            display: inline-block;
        }

        .nav-dropdown:has(.nav-dropdown-link.active) .nav-dropdown-toggle {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .nav-dropdown-toggle {
            cursor: pointer;
            background: none;
            border: none;
            font-family: inherit;
        }

        .nav-dropdown-menu {
            display: none;
            position: absolute;
            top: 100%;
            right: 0;
            background: white;
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            min-width: 200px;
            margin-top: 4px;
            z-index: 1000;
            overflow: hidden;
        }

        .nav-dropdown-menu.show {
            display: block;
        }

        .nav-dropdown-link {
            display: block;
            padding: 12px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .nav-dropdown-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
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
            width: 100%;
            margin: 0 auto;
        }

        .header {
            text-align: center;
            margin-bottom: 48px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) forwards;
            opacity: 0;
        }

        .header-title {
            font-size: 36px;
            font-weight: 800;
            letter-spacing: -1.5px;
            color: var(--gray-900);
            margin-bottom: 12px;
        }

        .gradient {
            background: linear-gradient(135deg, var(--blue-600) 0%, var(--blue-500) 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
        }

        .header-subtitle {
            font-size: 16px;
            color: var(--gray-500);
            line-height: 1.6;
        }

        @keyframes fade-up {
            from { opacity: 0; transform: translateY(24px); }
            to { opacity: 1; transform: translateY(0); }
        }

        .content-card {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(40px) saturate(180%);
            -webkit-backdrop-filter: blur(40px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 40px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.1s forwards;
            opacity: 0;
            margin-bottom: 24px;
        }

        .section-title {
            font-size: 20px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 16px;
            padding-bottom: 12px;
            border-bottom: 2px solid var(--blue-100);
        }

        .section-content {
            font-size: 15px;
            line-height: 1.8;
            color: var(--gray-700);
            margin-bottom: 24px;
        }

        .section-content p {
            margin-bottom: 16px;
        }

        .section-content strong {
            color: var(--gray-900);
            font-weight: 600;
        }

        .risk-list {
            list-style: none;
            padding: 0;
            margin: 16px 0;
        }

        .risk-list li {
            padding: 12px 16px;
            margin-bottom: 8px;
            background: var(--gray-50);
            border-left: 3px solid var(--blue-500);
            border-radius: 8px;
            font-size: 14px;
            color: var(--gray-700);
        }

        .risk-list li strong {
            color: var(--gray-900);
            display: block;
            margin-bottom: 4px;
        }

        .info-box {
            background: var(--blue-50);
            border: 1px solid var(--blue-200);
            border-radius: 12px;
            padding: 20px;
            margin: 24px 0;
        }

        .info-box p {
            font-size: 14px;
            line-height: 1.7;
            color: var(--gray-700);
            margin: 0;
        }

        .references {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(40px) saturate(180%);
            -webkit-backdrop-filter: blur(40px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 32px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.2s forwards;
            opacity: 0;
        }

        .references-title {
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 20px;
            display: flex;
            align-items: center;
            gap: 8px;
        }

        .references-list {
            list-style: none;
            padding: 0;
            margin: 0;
        }

        .references-list li {
            padding: 16px 0;
            border-bottom: 1px solid var(--gray-100);
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-600);
        }

        .references-list li:last-child {
            border-bottom: none;
        }

        .references-list strong {
            color: var(--gray-900);
            font-weight: 600;
        }

        .references-list em {
            color: var(--gray-700);
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

        .footer-link:hover {
            color: var(--gray-700);
        }

        @media (min-width: 768px) {
            .nav-links { display: flex; }
            .mobile-menu-btn { display: none; }
            .header-title { font-size: 48px; }
            .content-card { padding: 48px; }
            .footer-inner { flex-direction: row; justify-content: space-between; text-align: left; }
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
        <nav class="nav" role="navigation" aria-label="Main navigation">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo" aria-label="GasConsult.ai home">
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
                    <a href="/crisis" class="nav-link">Crisis Protocols</a>
                    <div class="nav-dropdown">
                        <button class="nav-link nav-dropdown-toggle" onclick="toggleNavDropdown(event)">More ▼</button>
                        <div class="nav-dropdown-menu">
                            <a href="/hypotension" class="nav-dropdown-link">IOH Predictor</a>
                            <a href="/difficult-airway" class="nav-dropdown-link">Difficult Airway</a>
                            <a href="/informed-consent" class="nav-dropdown-link active">Informed Consent</a>
                        </div>
                    </div>
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
            <a href="/difficult-airway" class="mobile-menu-link">Difficult Airway</a>
            <a href="/informed-consent" class="mobile-menu-link">Informed Consent</a>
        </div>

        <!-- Main Content -->
        <main class="main-content" id="main-content">
            <div class="header">
                <h1 class="header-title"><span class="gradient">Informed Consent</span> Generator</h1>
                <p class="header-subtitle">AI-powered patient-specific consent discussion guide</p>
            </div>

            {% if not summary %}
            <!-- Form State -->
            <div class="form-card">
                <form method="POST" action="/informed-consent">
                    <input type="hidden" name="csrf_token" value="{{ csrf_token() }}"/>

                    <!-- Patient Demographics -->
                    <div class="section-title">Patient Information</div>
                    <div class="input-row">
                        <div class="form-group">
                            <label class="form-label">Age (years)<span class="required">*</span></label>
                            <input type="number" name="age" class="form-input" placeholder="50" required min="1" max="120">
                        </div>
                        <div class="form-group">
                            <label class="form-label">ASA Classification<span class="required">*</span></label>
                            <select name="asa" class="form-input" required>
                                <option value="">Select...</option>
                                <option value="I">ASA I - Healthy</option>
                                <option value="II">ASA II - Mild systemic disease</option>
                                <option value="III">ASA III - Severe systemic disease</option>
                                <option value="IV">ASA IV - Life-threatening disease</option>
                                <option value="V">ASA V - Moribund</option>
                            </select>
                        </div>
                    </div>

                    <div class="section-divider"></div>

                    <!-- Anesthesia Details -->
                    <div class="section-title">Anesthesia & Procedure</div>

                    <div class="form-group">
                        <label class="form-label">Anesthesia Type<span class="required">*</span></label>
                        <div class="radio-group">
                            <label class="radio-option">
                                <input type="radio" name="anesthesia_type" value="General" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">General</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="anesthesia_type" value="Regional/Neuraxial" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">Regional/Neuraxial</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="anesthesia_type" value="Sedation/MAC" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">Sedation/MAC</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="anesthesia_type" value="Local" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">Local</span>
                            </label>
                        </div>
                    </div>

                    <div class="form-group">
                        <label class="form-label">Surgical Procedure<span class="required">*</span></label>
                        <input type="text" name="procedure" class="form-input" placeholder="e.g., Total knee arthroplasty" required>
                    </div>

                    <div class="section-divider"></div>

                    <!-- Comorbidities -->
                    <div class="section-title">Medical History</div>
                    <div class="form-group">
                        <label class="form-label">Comorbidities (select all that apply)</label>
                        <div class="checkbox-group">
                            <label class="checkbox-option">
                                <input type="checkbox" name="comorbidities" value="Cardiac Disease">
                                <span class="checkbox-label">Cardiac Disease</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="comorbidities" value="Pulmonary Disease">
                                <span class="checkbox-label">Pulmonary Disease</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="comorbidities" value="Diabetes">
                                <span class="checkbox-label">Diabetes</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="comorbidities" value="Renal Disease">
                                <span class="checkbox-label">Renal Disease</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="comorbidities" value="Hepatic Disease">
                                <span class="checkbox-label">Hepatic Disease</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="comorbidities" value="OSA">
                                <span class="checkbox-label">OSA</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="comorbidities" value="Obesity">
                                <span class="checkbox-label">Obesity</span>
                            </label>
                        </div>
                    </div>

                    <div class="section-divider"></div>

                    <!-- Special Considerations -->
                    <div class="section-title">Special Considerations</div>
                    <div class="form-group">
                        <label class="form-label">Select all that apply</label>
                        <div class="checkbox-group">
                            <label class="checkbox-option">
                                <input type="checkbox" name="special_considerations" value="Pregnancy">
                                <span class="checkbox-label">Pregnancy</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="special_considerations" value="Malignant Hyperthermia History">
                                <span class="checkbox-label">Malignant Hyperthermia History</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="special_considerations" value="Difficult Airway">
                                <span class="checkbox-label">Difficult Airway</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="special_considerations" value="Latex Allergy">
                                <span class="checkbox-label">Latex Allergy</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="special_considerations" value="Chronic Pain">
                                <span class="checkbox-label">Chronic Pain</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="special_considerations" value="Anticoagulation">
                                <span class="checkbox-label">Anticoagulation</span>
                            </label>
                        </div>
                    </div>

                    <div class="form-group">
                        <label class="form-label">Additional Notes</label>
                        <textarea name="notes" class="form-input" placeholder="Any additional patient concerns, questions, or special circumstances..."></textarea>
                    </div>

                    <button type="submit" class="submit-btn">Generate Consent Discussion</button>
                </form>
            </div>

            {% else %}
            <!-- Results State -->
            <div class="results-container">
                <div class="results-header">
                    <div class="results-icon">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <path d="M22 11.08V12a10 10 0 1 1-5.93-9.14"></path>
                            <polyline points="22 4 12 14.01 9 11.01"></polyline>
                        </svg>
                    </div>
                    <h1 class="results-title"><span style="color: var(--blue-600);">Informed Consent</span> Discussion Guide</h1>
                    <p class="results-subtitle">Patient-specific consent discussion with evidence-based risk information</p>
                </div>

                <div class="results-card">
                    <div class="results-content">
                        {{ summary|safe }}
                    </div>
                </div>

                {% if references %}
                <div class="references-card">
                    <div class="references-title">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path>
                            <path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path>
                        </svg>
                        <span>Evidence-Based References</span>
                    </div>
                    {% for ref in references %}
                    <div class="reference-item">
                        <span class="reference-number">[{{ loop.index }}]</span>
                        <a href="https://pubmed.ncbi.nlm.nih.gov/{{ ref.pmid }}/" target="_blank" rel="noopener noreferrer" class="reference-link">
                            {{ ref.title }}
                        </a>
                        <div class="reference-meta">{{ ref.authors }} - {{ ref.journal }}, {{ ref.year }}</div>
                    </div>
                    {% endfor %}
                </div>
                {% endif %}

                <div class="action-buttons">
                    <a href="/informed-consent" class="btn btn-primary">
                        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <line x1="12" y1="5" x2="12" y2="19"></line>
                            <line x1="5" y1="12" x2="19" y2="12"></line>
                        </svg>
                        New Consent
                    </a>
                    <button onclick="window.print()" class="btn btn-secondary">
                        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <polyline points="6 9 6 2 18 2 18 9"></polyline>
                            <path d="M6 18H4a2 2 0 0 1-2-2v-5a2 2 0 0 1 2-2h16a2 2 0 0 1 2 2v5a2 2 0 0 1-2 2h-2"></path>
                            <rect x="6" y="14" width="12" height="8"></rect>
                        </svg>
                        Print/Save PDF
                    </button>
                </div>
            </div>
            {% endif %}
        </main>

        <footer class="footer">
            <div class="footer-inner">
                <span class="footer-text">© 2025 GasConsult.ai</span>
                <div class="footer-links">
                    <a href="/privacy" class="footer-link">Privacy</a>
                    <a href="/terms" class="footer-link">Terms</a>
                    <a href="mailto:contact@gasconsult.ai" class="footer-link">Contact</a>
                </div>
            </div>
        </footer>
    </div>

    <script>
        // Mobile menu toggle
        function toggleMobileMenu() {
            const menu = document.getElementById('mobileMenu');
            const btn = document.querySelector('.mobile-menu-btn');
            if (menu && btn) {
                menu.classList.toggle('active');
                btn.classList.toggle('active');
            }
        }

        function toggleNavDropdown(e) {
            e.preventDefault();
            e.stopPropagation();
            const menu = e.target.nextElementSibling;
            if (menu) {
                menu.classList.toggle('show');
            }
        }

        // Close dropdown when clicking outside
        document.addEventListener('click', function() {
            document.querySelectorAll('.nav-dropdown-menu').forEach(m => m.classList.remove('show'));
        });

        // Radio button selection handler
        document.querySelectorAll('.radio-option').forEach(option => {
            option.addEventListener('click', function() {
                const radio = this.querySelector('input[type="radio"]');
                if (radio) {
                    radio.checked = true;
                    // Remove selected class from siblings
                    const name = radio.getAttribute('name');
                    document.querySelectorAll(`input[name="${name}"]`).forEach(r => {
                        r.closest('.radio-option').classList.remove('selected');
                    });
                    // Add selected class to clicked option
                    this.classList.add('selected');
                }
            });
        });

        // Checkbox selection handler
        document.querySelectorAll('.checkbox-option').forEach(option => {
            option.addEventListener('click', function(e) {
                if (e.target.tagName !== 'INPUT') {
                    const checkbox = this.querySelector('input[type="checkbox"]');
                    if (checkbox) {
                        checkbox.checked = !checkbox.checked;
                    }
                }
            });
        });
    </script>
</body>
</html>
"""

DIFFICULT_AIRWAY_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <!-- Google Analytics 4 -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=G-01NZYD1DPP"></script>
    <script>
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-01NZYD1DPP');
    </script>

    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Difficult Airway Assessment - GasConsult.ai</title>
    <meta name="description" content="AI-powered difficult airway assessment tool combining validated predictors with evidence-based management strategies.">

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
            --amber-50: #FFFBEB;
            --amber-500: #F59E0B;
            --green-50: #F0FDF4;
            --green-500: #22C55E;
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

        .nav-dropdown {
            position: relative;
            display: inline-block;
        }

        .nav-dropdown:has(.nav-dropdown-link.active) .nav-dropdown-toggle {
            color: var(--blue-600);
            background: var(--blue-50);
        }

        .nav-dropdown-toggle {
            cursor: pointer;
            background: none;
            border: none;
            font-family: inherit;
        }

        .nav-dropdown-menu {
            display: none;
            position: absolute;
            top: 100%;
            right: 0;
            background: white;
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            min-width: 200px;
            margin-top: 4px;
            z-index: 1000;
            overflow: hidden;
        }

        .nav-dropdown-menu.show {
            display: block;
        }

        .nav-dropdown-link {
            display: block;
            padding: 12px 18px;
            font-size: 14px;
            font-weight: 500;
            color: var(--gray-600);
            text-decoration: none;
            transition: all 0.2s ease;
        }

        .nav-dropdown-link:hover {
            color: var(--gray-900);
            background: rgba(0,0,0,0.04);
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
            width: 100%;
            margin: 0 auto;
        }

        .header {
            text-align: center;
            margin-bottom: 48px;
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) forwards;
            opacity: 0;
        }

        .header-title {
            font-size: 36px;
            font-weight: 800;
            letter-spacing: -1.5px;
            color: var(--gray-900);
            margin-bottom: 12px;
        }

        .gradient {
            background: linear-gradient(135deg, var(--blue-600) 0%, var(--blue-500) 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
        }

        .header-subtitle {
            font-size: 16px;
            color: var(--gray-500);
            line-height: 1.6;
        }

        @keyframes fade-up {
            from { opacity: 0; transform: translateY(24px); }
            to { opacity: 1; transform: translateY(0); }
        }

        .content-card {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(40px) saturate(180%);
            -webkit-backdrop-filter: blur(40px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 40px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.1s forwards;
            opacity: 0;
            margin-bottom: 24px;
        }

        .section-title {
            font-size: 20px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 16px;
            padding-bottom: 12px;
            border-bottom: 2px solid var(--blue-100);
        }

        .section-content {
            font-size: 15px;
            line-height: 1.8;
            color: var(--gray-700);
            margin-bottom: 24px;
        }

        .section-content p {
            margin-bottom: 16px;
        }

        .section-content strong {
            color: var(--gray-900);
            font-weight: 600;
        }

        .predictor-list {
            list-style: none;
            padding: 0;
            margin: 16px 0;
        }

        .predictor-list li {
            padding: 12px 16px;
            margin-bottom: 8px;
            background: var(--gray-50);
            border-left: 3px solid var(--blue-500);
            border-radius: 8px;
            font-size: 14px;
            color: var(--gray-700);
        }

        .predictor-list li strong {
            color: var(--gray-900);
            display: block;
            margin-bottom: 4px;
        }

        .algorithm-box {
            background: var(--blue-50);
            border: 2px solid var(--blue-200);
            border-radius: 12px;
            padding: 24px;
            margin: 24px 0;
        }

        .algorithm-box h3 {
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 16px;
        }

        .algorithm-steps {
            list-style: none;
            padding: 0;
            margin: 0;
        }

        .algorithm-steps li {
            padding: 12px 16px 12px 40px;
            margin-bottom: 8px;
            background: white;
            border-radius: 8px;
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-700);
            position: relative;
        }

        .algorithm-steps li::before {
            content: attr(data-step);
            position: absolute;
            left: 12px;
            top: 12px;
            width: 20px;
            height: 20px;
            background: var(--blue-500);
            color: white;
            border-radius: 50%;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 11px;
            font-weight: 700;
        }

        .info-box {
            background: var(--amber-50);
            border: 1px solid var(--amber-500);
            border-radius: 12px;
            padding: 20px;
            margin: 24px 0;
        }

        .info-box p {
            font-size: 14px;
            line-height: 1.7;
            color: var(--gray-700);
            margin: 0;
        }

        .references {
            background: rgba(255,255,255,0.8);
            backdrop-filter: blur(40px) saturate(180%);
            -webkit-backdrop-filter: blur(40px) saturate(180%);
            border: 1px solid rgba(255,255,255,0.9);
            border-radius: 20px;
            padding: 32px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.02), 0 4px 16px rgba(0,0,0,0.04), 0 24px 80px rgba(0,0,0,0.06), inset 0 1px 0 rgba(255,255,255,0.8);
            animation: fade-up 0.8s cubic-bezier(0.16,1,0.3,1) 0.2s forwards;
            opacity: 0;
        }

        .references-title {
            font-size: 18px;
            font-weight: 700;
            color: var(--gray-900);
            margin-bottom: 20px;
            display: flex;
            align-items: center;
            gap: 8px;
        }

        .references-list {
            list-style: none;
            padding: 0;
            margin: 0;
        }

        .references-list li {
            padding: 16px 0;
            border-bottom: 1px solid var(--gray-100);
            font-size: 14px;
            line-height: 1.6;
            color: var(--gray-600);
        }

        .references-list li:last-child {
            border-bottom: none;
        }

        .references-list strong {
            color: var(--gray-900);
            font-weight: 600;
        }

        .references-list em {
            color: var(--gray-700);
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

        .footer-link:hover {
            color: var(--gray-700);
        }

        @media (min-width: 768px) {
            .nav-links { display: flex; }
            .mobile-menu-btn { display: none; }
            .header-title { font-size: 48px; }
            .content-card { padding: 48px; }
            .footer-inner { flex-direction: row; justify-content: space-between; text-align: left; }
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
        <nav class="nav" role="navigation" aria-label="Main navigation">
            <div class="nav-inner">
                <a href="/?clear=1" class="logo" aria-label="GasConsult.ai home">
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
                    <a href="/crisis" class="nav-link">Crisis Protocols</a>
                    <div class="nav-dropdown">
                        <button class="nav-link nav-dropdown-toggle" onclick="toggleNavDropdown(event)">More ▼</button>
                        <div class="nav-dropdown-menu">
                            <a href="/hypotension" class="nav-dropdown-link">IOH Predictor</a>
                            <a href="/difficult-airway" class="nav-dropdown-link active">Difficult Airway</a>
                            <a href="/informed-consent" class="nav-dropdown-link">Informed Consent</a>
                        </div>
                    </div>
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
            <a href="/difficult-airway" class="mobile-menu-link">Difficult Airway</a>
            <a href="/informed-consent" class="mobile-menu-link">Informed Consent</a>
        </div>

        <!-- Main Content -->
        <main class="main-content" id="main-content">
            <div class="header">
                <h1 class="header-title"><span class="gradient">Difficult Airway</span> Assessment</h1>
                <p class="header-subtitle">AI-powered risk stratification and evidence-based management planning</p>
            </div>

            {% if not summary %}
            <!-- Form State -->
            <div class="form-card">
                <form method="POST" action="/difficult-airway">
                    <input type="hidden" name="csrf_token" value="{{ csrf_token() }}"/>

                    <!-- Patient Demographics -->
                    <div class="section-title">Patient Demographics</div>
                    <div class="input-row">
                        <div class="form-group">
                            <label class="form-label">Age (years)<span class="required">*</span></label>
                            <input type="number" name="age" class="form-input" placeholder="50" required min="1" max="120">
                        </div>
                        <div class="form-group">
                            <label class="form-label">BMI (kg/m²)<span class="required">*</span></label>
                            <input type="number" name="bmi" class="form-input" placeholder="25" required min="10" max="100" step="0.1">
                        </div>
                    </div>

                    <div class="section-divider"></div>

                    <!-- Airway Assessment -->
                    <div class="section-title">Airway Examination</div>

                    <div class="form-group">
                        <label class="form-label">Mallampati Classification<span class="required">*</span></label>
                        <div class="radio-group">
                            <label class="radio-option">
                                <input type="radio" name="mallampati" value="I" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">Class I</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="mallampati" value="II" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">Class II</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="mallampati" value="III" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">Class III</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="mallampati" value="IV" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">Class IV</span>
                            </label>
                        </div>
                    </div>

                    <div class="form-group">
                        <label class="form-label">Thyromental Distance<span class="required">*</span></label>
                        <div class="radio-group">
                            <label class="radio-option">
                                <input type="radio" name="thyromental" value="<6cm" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">&lt;6 cm</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="thyromental" value="6-6.5cm" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">6-6.5 cm</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="thyromental" value=">6.5cm" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">&gt;6.5 cm</span>
                            </label>
                        </div>
                    </div>

                    <div class="form-group">
                        <label class="form-label">Mouth Opening<span class="required">*</span></label>
                        <div class="radio-group">
                            <label class="radio-option">
                                <input type="radio" name="mouth_opening" value="<3cm" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">&lt;3 cm</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="mouth_opening" value="3-4cm" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">3-4 cm</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="mouth_opening" value=">4cm" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">&gt;4 cm</span>
                            </label>
                        </div>
                    </div>

                    <div class="form-group">
                        <label class="form-label">Neck Extension<span class="required">*</span></label>
                        <div class="radio-group">
                            <label class="radio-option">
                                <input type="radio" name="neck_extension" value="Limited" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">Limited</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="neck_extension" value="Moderate" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">Moderate</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="neck_extension" value="Normal" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">Normal</span>
                            </label>
                        </div>
                    </div>

                    <div class="section-divider"></div>

                    <!-- Risk Factors -->
                    <div class="section-title">Risk Factors</div>
                    <div class="form-group">
                        <label class="form-label">Select all that apply</label>
                        <div class="checkbox-group">
                            <label class="checkbox-option">
                                <input type="checkbox" name="risk_factors" value="Previous Difficult Intubation">
                                <span class="checkbox-label">Previous Difficult Intubation</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="risk_factors" value="OSA">
                                <span class="checkbox-label">OSA</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="risk_factors" value="Beard">
                                <span class="checkbox-label">Beard</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="risk_factors" value="Prominent Incisors">
                                <span class="checkbox-label">Prominent Incisors</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="risk_factors" value="Short Neck">
                                <span class="checkbox-label">Short Neck</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="risk_factors" value="Large Tongue">
                                <span class="checkbox-label">Large Tongue</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="risk_factors" value="Facial Trauma">
                                <span class="checkbox-label">Facial Trauma</span>
                            </label>
                            <label class="checkbox-option">
                                <input type="checkbox" name="risk_factors" value="C-Spine Pathology">
                                <span class="checkbox-label">C-Spine Pathology</span>
                            </label>
                        </div>
                    </div>

                    <div class="section-divider"></div>

                    <!-- Procedure Details -->
                    <div class="section-title">Procedure Information</div>
                    <div class="form-group">
                        <label class="form-label">Planned Procedure<span class="required">*</span></label>
                        <input type="text" name="procedure" class="form-input" placeholder="e.g., Laparoscopic cholecystectomy" required>
                    </div>

                    <div class="form-group">
                        <label class="form-label">Case Type<span class="required">*</span></label>
                        <div class="radio-group">
                            <label class="radio-option">
                                <input type="radio" name="case_type" value="Elective" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">Elective</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="case_type" value="Urgent" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">Urgent</span>
                            </label>
                            <label class="radio-option">
                                <input type="radio" name="case_type" value="Emergency" required>
                                <div class="radio-visual"></div>
                                <span class="radio-label">Emergency</span>
                            </label>
                        </div>
                    </div>

                    <div class="form-group">
                        <label class="form-label">Additional Notes</label>
                        <textarea name="notes" class="form-input" placeholder="Any additional airway concerns, previous anesthesia complications, or special considerations..."></textarea>
                    </div>

                    <button type="submit" class="submit-btn">Generate Assessment</button>
                </form>
            </div>

            {% else %}
            <!-- Results State -->
            <div class="results-container">
                <div class="results-header">
                    <div class="results-icon">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <path d="M22 11.08V12a10 10 0 1 1-5.93-9.14"></path>
                            <polyline points="22 4 12 14.01 9 11.01"></polyline>
                        </svg>
                    </div>
                    <h1 class="results-title"><span style="color: var(--blue-600);">Difficult Airway</span> Assessment Complete</h1>
                    <p class="results-subtitle">Evidence-based risk stratification and management recommendations</p>
                </div>

                <div class="results-card">
                    <div class="results-content">
                        {{ summary|safe }}
                    </div>
                </div>

                {% if references %}
                <div class="references-card">
                    <div class="references-title">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path>
                            <path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path>
                        </svg>
                        <span>Evidence-Based References</span>
                    </div>
                    {% for ref in references %}
                    <div class="reference-item">
                        <span class="reference-number">[{{ loop.index }}]</span>
                        <a href="https://pubmed.ncbi.nlm.nih.gov/{{ ref.pmid }}/" target="_blank" rel="noopener noreferrer" class="reference-link">
                            {{ ref.title }}
                        </a>
                        <div class="reference-meta">{{ ref.authors }} - {{ ref.journal }}, {{ ref.year }}</div>
                    </div>
                    {% endfor %}
                </div>
                {% endif %}

                <div class="action-buttons">
                    <a href="/difficult-airway" class="btn btn-primary">
                        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <line x1="12" y1="5" x2="12" y2="19"></line>
                            <line x1="5" y1="12" x2="19" y2="12"></line>
                        </svg>
                        New Assessment
                    </a>
                    <button onclick="window.print()" class="btn btn-secondary">
                        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <polyline points="6 9 6 2 18 2 18 9"></polyline>
                            <path d="M6 18H4a2 2 0 0 1-2-2v-5a2 2 0 0 1 2-2h16a2 2 0 0 1 2 2v5a2 2 0 0 1-2 2h-2"></path>
                            <rect x="6" y="14" width="12" height="8"></rect>
                        </svg>
                        Print/Save PDF
                    </button>
                </div>
            </div>
            {% endif %}
        </main>

        <footer class="footer">
            <div class="footer-inner">
                <span class="footer-text">© 2025 GasConsult.ai</span>
                <div class="footer-links">
                    <a href="/privacy" class="footer-link">Privacy</a>
                    <a href="/terms" class="footer-link">Terms</a>
                    <a href="mailto:contact@gasconsult.ai" class="footer-link">Contact</a>
                </div>
            </div>
        </footer>
    </div>

    <script>
        // Mobile menu toggle
        function toggleMobileMenu() {
            const menu = document.getElementById('mobileMenu');
            const btn = document.querySelector('.mobile-menu-btn');
            if (menu && btn) {
                menu.classList.toggle('active');
                btn.classList.toggle('active');
            }
        }

        function toggleNavDropdown(e) {
            e.preventDefault();
            e.stopPropagation();
            const menu = e.target.nextElementSibling;
            if (menu) {
                menu.classList.toggle('show');
            }
        }

        // Close dropdown when clicking outside
        document.addEventListener('click', function() {
            document.querySelectorAll('.nav-dropdown-menu').forEach(m => m.classList.remove('show'));
        });

        // Radio button selection handler
        document.querySelectorAll('.radio-option').forEach(option => {
            option.addEventListener('click', function() {
                const radio = this.querySelector('input[type="radio"]');
                if (radio) {
                    radio.checked = true;
                    // Remove selected class from siblings
                    const name = radio.getAttribute('name');
                    document.querySelectorAll(`input[name="${name}"]`).forEach(r => {
                        r.closest('.radio-option').classList.remove('selected');
                    });
                    // Add selected class to clicked option
                    this.classList.add('selected');
                }
            });
        });

        // Checkbox selection handler
        document.querySelectorAll('.checkbox-option').forEach(option => {
            option.addEventListener('click', function(e) {
                if (e.target.tagName !== 'INPUT') {
                    const checkbox = this.querySelector('input[type="checkbox"]');
                    if (checkbox) {
                        checkbox.checked = !checkbox.checked;
                    }
                }
            });
        });
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
            print(f"[DEBUG] [STREAM] Before saving - session has {len(session.get('messages', []))} messages")
            for i, msg in enumerate(session.get('messages', [])):
                print(f"[DEBUG] [STREAM]   Message {i}: role={msg.get('role')}, content_length={len(msg.get('content', ''))}")

            # CRITICAL FIX: Get a copy of messages list to force Flask session change detection
            # Modifying nested list items directly (session['messages'][-1] = ...) doesn't
            # properly trigger Flask's filesystem session backend to persist changes
            messages = list(session.get('messages', []))

            # Check if last message is an empty assistant placeholder (from homepage redirect)
            # If so, update it instead of appending a new one
            if (messages and
                len(messages) > 0 and
                messages[-1].get('role') == 'assistant' and
                messages[-1].get('content') == ''):
                # Update the placeholder
                messages[-1] = {
                    "role": "assistant",
                    "content": cleaned_response,
                    "references": refs,
                    "num_papers": num_papers,
                    "evidence_strength": evidence_strength
                }
                print(f"[DEBUG] [STREAM] Updated existing placeholder assistant message")
            else:
                # Append new message (for AJAX submissions from chat page)
                messages.append({
                    "role": "assistant",
                    "content": cleaned_response,
                    "references": refs,
                    "num_papers": num_papers,
                    "evidence_strength": evidence_strength
                })
                print(f"[DEBUG] [STREAM] Appended new assistant message")

            # Reassign entire list to trigger Flask session change detection
            session['messages'] = messages
            session.modified = True

            # CRITICAL: Explicitly save session for SSE streaming responses
            # Even with Redis, generator functions need manual session save because
            # Flask only auto-saves at the end of the request, but generators
            # are still yielding data. We must force the save NOW before continuing.
            try:
                # Create a mock response object for save_session
                from werkzeug.wrappers import Response as WerkzeugResponse
                mock_response = WerkzeugResponse()
                app.session_interface.save_session(app, session, mock_response)
                print(f"[DEBUG] [STREAM] Explicitly saved session to Redis")
            except Exception as e:
                print(f"[ERROR] [STREAM] Failed to explicitly save session: {e}")

            print(f"[DEBUG] [STREAM] After saving - session has {len(session['messages'])} messages")
            print(f"[DEBUG] [STREAM] Verified last message content_length: {len(messages[-1].get('content', ''))}")

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

# ============================================================================
# PubMed Query Caching Functions
# ============================================================================

def generate_cache_key(search_term, retmax=20):
    """Generate a cache key for PubMed queries"""
    key_string = f"pubmed:{search_term}:{retmax}"
    return f"pubmed_query:{hashlib.md5(key_string.encode()).hexdigest()}"

def get_cached_pubmed_results(search_term, retmax=20):
    """Retrieve cached PubMed search results"""
    if not redis_client:
        return None

    try:
        cache_key = generate_cache_key(search_term, retmax)
        cached = redis_client.get(cache_key)
        if cached:
            print(f"[CACHE HIT] Found cached results for query: {search_term[:50]}")
            return json.loads(cached)
        print(f"[CACHE MISS] No cached results for query: {search_term[:50]}")
        return None
    except Exception as e:
        print(f"[CACHE ERROR] Failed to retrieve from cache: {e}")
        return None

def cache_pubmed_results(search_term, ids, retmax=20, ttl=86400):
    """Cache PubMed search results (default TTL: 24 hours)"""
    if not redis_client:
        return False

    try:
        cache_key = generate_cache_key(search_term, retmax)
        cache_data = {
            'ids': ids,
            'search_term': search_term,
            'timestamp': str(datetime.now()),
            'count': len(ids)
        }
        redis_client.setex(cache_key, ttl, json.dumps(cache_data))
        print(f"[CACHE SAVE] Cached {len(ids)} results for query: {search_term[:50]}")
        return True
    except Exception as e:
        print(f"[CACHE ERROR] Failed to save to cache: {e}")
        return False

def get_cached_paper_metadata(pmid):
    """Retrieve cached paper metadata"""
    if not redis_client:
        return None

    try:
        cache_key = f"pubmed_paper:{pmid}"
        cached = redis_client.get(cache_key)
        if cached:
            return json.loads(cached)
        return None
    except Exception as e:
        print(f"[CACHE ERROR] Failed to retrieve paper {pmid} from cache: {e}")
        return None

def cache_paper_metadata(pmid, metadata, ttl=604800):
    """Cache paper metadata (default TTL: 7 days)"""
    if not redis_client:
        return False

    try:
        cache_key = f"pubmed_paper:{pmid}"
        redis_client.setex(cache_key, ttl, json.dumps(metadata))
        return True
    except Exception as e:
        print(f"[CACHE ERROR] Failed to cache paper {pmid}: {e}")
        return False

# ============================================================================
# Main Route
# ============================================================================

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

            # Detect well-established topics that should have lots of evidence
            well_established_topics = [
                'malignant hyperthermia', 'difficult airway', 'aspiration', 'local anesthetic systemic toxicity',
                'last', 'propofol infusion syndrome', 'pris', 'anaphylaxis', 'tranexamic acid', 'txa',
                'postoperative nausea', 'ponv', 'rapid sequence', 'rsi', 'spinal anesthesia', 'epidural',
                'general anesthesia', 'blood transfusion', 'massive transfusion', 'coagulopathy'
            ]
            is_established_topic = any(topic in query_lower for topic in well_established_topics)

            # Adjust date range for established topics (go back further for classic evidence)
            date_range = '("2010/01/01"[PDAT] : "3000"[PDAT])' if is_established_topic else '("2018/01/01"[PDAT] : "3000"[PDAT])'

            # Customize search based on question type
            if is_followup:
                # Broader search for follow-ups
                search_term = f'({q}) AND ("2010/01/01"[PDAT] : "3000"[PDAT])'
            else:
                # Customize search filters based on question type with tighter date ranges
                if question_type == 'dosing':
                    # Prioritize guidelines and reviews for dosing questions
                    search_term = (
                        f'({q}) AND '
                        f'(guideline[pt] OR "practice guideline"[pt] OR review[pt] OR meta-analysis[pt]) AND '
                        f'{date_range}'
                    )
                elif question_type == 'safety':
                    # Include adverse effects and safety studies
                    search_term = (
                        f'({q}) AND '
                        f'(systematic review[pt] OR meta-analysis[pt] OR "randomized controlled trial"[pt] OR '
                        f'guideline[pt] OR "adverse effects"[sh] OR safety[ti]) AND '
                        f'{date_range}'
                    )
                elif question_type == 'comparison':
                    # Prioritize comparative studies
                    search_term = (
                        f'({q}) AND '
                        f'(systematic review[pt] OR meta-analysis[pt] OR "randomized controlled trial"[pt] OR '
                        f'"comparative study"[pt]) AND '
                        f'{date_range}'
                    )
                else:
                    # Default search (general, mechanism, management)
                    search_term = (
                        f'({q}) AND '
                        f'(systematic review[pt] OR meta-analysis[pt] OR "randomized controlled trial"[pt] OR '
                        f'"Cochrane Database Syst Rev"[ta] OR guideline[pt]) AND '
                        f'{date_range}'
                    )

            # Add negation modifier if detected
            if negation_search_modifier:
                search_term += negation_search_modifier

            print(f"[DEBUG] Search term: '{search_term[:150]}...'")

            # Multi-tier search strategy to find relevant papers
            ids = []

            # Tier 1: Broad perioperative/critical care context (inclusive of related specialties)
            # Anesthesia is closely tied to critical care, ICU, surgery, emergency medicine, and pain
            # Be INCLUSIVE of related fields while excluding obviously unrelated specialties
            periop_context = (
                '(anesthesia[MeSH Terms] OR anesthesiology[MeSH Terms] OR anesthetics[MeSH Terms] OR '
                'perioperative[tiab] OR intraoperative[tiab] OR postoperative[tiab] OR preoperative[tiab] OR '
                '"critical care"[MeSH Terms] OR "intensive care"[tiab] OR "ICU"[tiab] OR "critical illness"[tiab] OR '
                '"surgical anesthesia"[tiab] OR "anesthetic management"[tiab] OR '
                'surgery[MeSH Subheading] OR surgical[tiab] OR "operating room"[tiab] OR '
                '"emergency medicine"[MeSH Terms] OR resuscitation[MeSH Terms] OR "trauma"[tiab] OR '
                '"pain management"[tiab] OR "acute pain"[tiab] OR "pain medicine"[tiab] OR '
                'airway[tiab] OR "airway management"[tiab] OR ventilation[tiab] OR hemodynamic[tiab] OR '
                '"fluid management"[tiab] OR "blood transfusion"[tiab] OR shock[tiab])'
            )

            # Include papers from major related journals (anesthesia, critical care, surgery, emergency, general medical)
            # This is an OR boost, not a filter - papers in other journals are NOT excluded
            related_journals = (
                '("Anesthesiology"[Journal] OR "Anesthesia and Analgesia"[Journal] OR '
                '"British Journal of Anaesthesia"[Journal] OR "Anaesthesia"[Journal] OR '
                '"Critical Care Medicine"[Journal] OR "Intensive Care Medicine"[Journal] OR '
                '"Critical Care"[Journal] OR "Shock"[Journal] OR '
                '"JAMA Surgery"[Journal] OR "Annals of Surgery"[Journal] OR "Surgery"[Journal] OR '
                '"New England Journal of Medicine"[Journal] OR "JAMA"[Journal] OR '
                '"Lancet"[Journal] OR "Chest"[Journal] OR "Circulation"[Journal] OR '
                '"Emergency Medicine Journal"[Journal] OR "Academic Emergency Medicine"[Journal])'
            )

            tier1_search_term = f'{search_term} AND ({periop_context} OR {related_journals})'
            cached_result = get_cached_pubmed_results(tier1_search_term, retmax=20)

            if cached_result:
                ids = cached_result.get('ids', [])
                print(f"[DEBUG] Tier 1 (CACHED): Found {len(ids)} papers")
            else:
                try:
                    print(f"[DEBUG] Tier 1: Searching PubMed with anesthesiology context...")
                    handle = Entrez.esearch(db="pubmed", term=tier1_search_term, retmax=20, sort="relevance")
                    result = Entrez.read(handle)
                    ids = result.get("IdList", [])
                    print(f"[DEBUG] Tier 1 found {len(ids)} papers")
                    # Cache the results
                    cache_pubmed_results(tier1_search_term, ids, retmax=20)
                except Exception as e:
                    print(f"[ERROR] Tier 1 search failed: {e}")
                    ids = []

            # Tier 2: Even broader perioperative/critical care context (very inclusive)
            # Uses text words instead of MeSH for more flexibility
            # Includes critical care, ICU, surgical, emergency - all relevant to anesthesia practice
            if len(ids) < 5 and not is_followup:
                broader_periop = (
                    '(anesthesia[tiab] OR anesthesiology[tiab] OR anesthetic[tiab] OR '
                    'anaesthesia[tiab] OR anaesthetic[tiab] OR '
                    'perioperative[tiab] OR intraoperative[tiab] OR postoperative[tiab] OR '
                    'critical[tiab] OR ICU[tiab] OR intensive[tiab] OR '
                    'surgical[tiab] OR surgery[tiab] OR operative[tiab] OR '
                    'emergency[tiab] OR resuscitation[tiab] OR trauma[tiab] OR '
                    'airway[tiab] OR ventilat[tiab] OR hemodynamic[tiab] OR '
                    'pain[tiab] OR sedation[tiab])'
                )
                tier2_search_term = f'{search_term} AND {broader_periop}'

                # Check cache first
                cached_result = get_cached_pubmed_results(tier2_search_term, retmax=20)

                if cached_result:
                    tier2_ids = cached_result.get('ids', [])
                    ids = list(set(ids + tier2_ids))
                    print(f"[DEBUG] Tier 2 (CACHED): Found {len(tier2_ids)} papers, total unique: {len(ids)}")
                else:
                    try:
                        print(f"[DEBUG] Tier 2: Broader anesthesia context (still anesthesia-focused)...")
                        handle = Entrez.esearch(db="pubmed", term=tier2_search_term, retmax=20, sort="relevance")
                        result = Entrez.read(handle)
                        tier2_ids = result.get("IdList", [])
                        # Combine and deduplicate
                        ids = list(set(ids + tier2_ids))
                        print(f"[DEBUG] Tier 2 found {len(tier2_ids)} papers, total unique: {len(ids)}")
                        # Cache the results
                        cache_pubmed_results(tier2_search_term, tier2_ids, retmax=20)
                    except Exception as e:
                        print(f"[ERROR] Tier 2 search failed: {e}")

            # Tier 3: High-quality studies with broad perioperative context and extended date range
            if len(ids) < 5 and not is_followup:
                try:
                    print(f"[DEBUG] Tier 3: High-quality studies with broad perioperative context...")
                    # Maintain perioperative/critical care context - inclusive of related fields
                    tier3_periop = (
                        '(anesthesia[tiab] OR anesthetic[tiab] OR perioperative[tiab] OR '
                        'surgical[tiab] OR critical[tiab] OR ICU[tiab] OR intensive[tiab] OR '
                        'emergency[tiab] OR resuscitation[tiab])'
                    )
                    broader_term = (
                        f'({q}) AND {tier3_periop} AND '
                        f'(systematic review[pt] OR meta-analysis[pt] OR "randomized controlled trial"[pt] OR guideline[pt] OR review[pt]) AND '
                        f'("2010/01/01"[PDAT] : "3000"[PDAT])'
                    )
                    handle = Entrez.esearch(db="pubmed", term=broader_term, retmax=15, sort="relevance")
                    result = Entrez.read(handle)
                    tier3_ids = result.get("IdList", [])
                    # Combine and deduplicate
                    ids = list(set(ids + tier3_ids))
                    print(f"[DEBUG] Tier 3 found {len(tier3_ids)} papers, total unique: {len(ids)}")
                except Exception as e:
                    print(f"[ERROR] Tier 3 search failed: {e}")

            # Tier 4: Last resort - very broad perioperative context, extended date range
            if len(ids) < 3 and not is_followup:
                try:
                    print(f"[DEBUG] Tier 4: Very broad perioperative context, extended date range...")
                    # Most inclusive - captures anesthesia, critical care, surgical, emergency, ICU
                    tier4_periop = (
                        f'({q}) AND (anesthesia[tiab] OR anesthetic[tiab] OR perioperative[tiab] OR '
                        f'critical[tiab] OR ICU[tiab] OR surgical[tiab] OR emergency[tiab]) AND '
                        f'(review[pt] OR "randomized controlled trial"[pt] OR "clinical trial"[pt]) AND '
                        f'("2005/01/01"[PDAT] : "3000"[PDAT])'
                    )
                    handle = Entrez.esearch(db="pubmed", term=tier4_periop, retmax=15, sort="relevance")
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

CRITICAL: SELF-VERIFICATION PROTOCOL (Apply BEFORE finalizing your answer):
Before providing your answer, internally verify:
1. **Dosing Accuracy**: Are all drug doses within standard ranges? Cross-check against ASA guidelines and major textbooks
2. **Contraindication Check**: Have you mentioned absolute contraindications and safety warnings?
3. **Consistency Check**: Does your answer align with current ASA/ESA guidelines and standard practice?
4. **Self-Questioning**: Ask yourself "How do I know this is correct?" Acknowledge uncertainty where appropriate.
5. **Completeness**: For dosing questions, include route, typical range, and maximum doses. For safety questions, cover both common and serious risks.

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
                    # CRITICAL: Convert all Biopython StringElement objects to plain Python strings
                    # for Redis serialization compatibility
                    title = str(art.get("ArticleTitle", "No title"))
                    # Truncate abstracts to 1200 chars for better citation accuracy
                    abstract = " ".join(str(t) for t in art.get("Abstract", {}).get("AbstractText", [])) if art.get("Abstract") else ""
                    abstract = abstract[:1200] + "..." if len(abstract) > 1200 else abstract
                    authors = ", ".join([str(a.get("LastName","")) + " " + (str(a.get("ForeName",""))[:1]+"." if a.get("ForeName") else "") for a in art.get("AuthorList",[])[:3]])  # Reduced to 3 authors
                    journal = str(art["Journal"].get("Title", "Unknown"))
                    year = str(art["Journal"]["JournalIssue"]["PubDate"].get("Year", "N/A"))
                    pmid = str(p["MedlineCitation"]["PMID"])

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
                        "sort_priority": study_classification['sort_priority'],
                        "abstract": abstract  # Store abstract for numbered context
                    })
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

            # Build numbered context with full details (helps GPT understand which paper is which)
            context = ""
            for i, ref in enumerate(refs, 1):
                context += f"Paper [{i}]:\nTitle: {ref['title']}\nAbstract: {ref.get('abstract', 'No abstract available')}\nAuthors: {ref['authors']}\nJournal: {ref['journal']} ({ref['year']})\nPMID: {ref['pmid']}\n\n"

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

CRITICAL: SELF-VERIFICATION PROTOCOL (Apply BEFORE finalizing your answer):
Before providing your answer, internally verify:
1. **Dosing Accuracy**: Are all drug doses within standard ranges? Cross-check against ASA guidelines, package inserts, and major textbooks (Miller's, Barash, Stoelting's)
2. **Contraindication Check**: Have you mentioned absolute contraindications? Are there missing safety warnings?
3. **Citation Verification**: Does each cited paper's abstract ACTUALLY support the specific claim? Remove citations that don't directly support the statement.
4. **Consistency Check**: Does your answer align with current ASA/ESA guidelines and standard practice?
5. **Self-Questioning**: Ask yourself "How do I know this is correct?" If uncertain about any specific claim, either verify it against the provided papers or acknowledge uncertainty.
6. **Completeness**: For dosing questions, have you included route, typical range, and maximum doses? For safety questions, have you covered both common and serious risks?

INSTRUCTIONS:
1. Include specific dosages (mg/kg), contraindications, side effects, and monitoring when relevant
2. For acute situations, provide step-by-step protocols with drugs and doses
3. **IN-TEXT CITATIONS (CRITICAL):**
   - Use numbered citations [1], [2], [3] etc. corresponding to the numbered papers above
   - ONLY cite a paper [N] if that SPECIFIC paper's abstract directly supports the exact claim you're making
   - Verify each citation: Does Paper [N]'s abstract contain this specific fact/dose/finding?
   - If a claim is from your general knowledge but NOT explicitly in the provided abstracts, DO NOT add a citation
   - Better to have NO citation than an INACCURATE citation
   - When citing: use [1], [2], etc. - NO author names in text
   - Multiple papers can support one claim: "X is effective [1][2][3]"
   - Place citations immediately after the claim: "TXA reduces blood loss [1][2]."
   - If papers discuss a topic generally but don't support your specific claim, omit the citation
4. Be conversational but clinically complete - like talking to a colleague
5. HTML format: <h3> for sections, <p> for paragraphs, <strong> for emphasis, <ul><li> for lists
6. CRITICAL: Return ONLY the HTML content - do NOT wrap your response in markdown code fences (```html or ```)

CITATION VERIFICATION CHECKLIST (Check EACH citation before adding):
   ❌ WRONG: "Propofol 2-3 mg/kg is recommended [1]" when Paper [1]'s abstract says "induction agents" generally
   ❌ WRONG: "TXA reduces blood loss by 30% [2]" when Paper [2]'s abstract doesn't give specific percentage
   ❌ WRONG: "Common side effects include nausea [3]" when Paper [3]'s abstract doesn't mention side effects
   ✅ CORRECT: "TXA reduces blood loss in spine surgery [1][2]" when both Paper [1] and [2] abstracts explicitly state this
   ✅ CORRECT: "Standard monitoring includes pulse oximetry" (NO citation - general knowledge not from papers)

IMPORTANT: Use numbered citations [1], [2], [3] throughout your answer to reference the specific papers above. Only cite when the paper's abstract directly supports the specific claim. If a claim is general knowledge or not supported by the abstracts, omit the citation.

Example with proper numbered citations:
"<h3>Acute Bronchospasm Management</h3>
<p><strong>Immediate Actions:</strong><br>
100% oxygen and hand-ventilate to assess compliance. Deepen anesthesia with propofol or increase volatile agent to 2+ MAC.</p>
<p><strong>Bronchodilators:</strong><br>
Inhaled beta-2 agonists are first-line treatment for acute bronchospasm [1][2]. Response typically seen within minutes. Severe cases may require IV epinephrine [3].</p>
<p><strong>Monitoring:</strong><br>
Watch for auto-PEEP, pneumothorax, and cardiovascular compromise from high airway pressures.</p>"

NOTE: Numbers [1], [2], [3] must correspond to Paper [1], Paper [2], Paper [3] listed above. Only cite when abstracts explicitly support claims.

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
            print(f"[DEBUG] Session messages before redirect: {len(session.get('messages', []))} messages")
            for i, msg in enumerate(session.get('messages', [])):
                print(f"[DEBUG]   Message {i}: role={msg.get('role')}, content_length={len(msg.get('content', ''))}")
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
    print(f"[DEBUG] GET / - Session messages detail:")
    for i, msg in enumerate(session.get('messages', [])):
        print(f"[DEBUG]   Message {i}: role={msg.get('role')}, content_length={len(msg.get('content', ''))}, has_refs={len(msg.get('references', []))}")

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

@app.route("/evidence")
def evidence():
    """Evidence-Based Methodology - How PubMed integration works"""
    return render_template_string(EVIDENCE_HTML)

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
                        # CRITICAL: Convert all Biopython StringElement objects to plain Python strings
                        # for Redis serialization compatibility
                        title = str(art.get("ArticleTitle", "No title"))
                        abstract = " ".join(str(t) for t in art.get("Abstract", {}).get("AbstractText", [])) if art.get("Abstract") else ""
                        authors = ", ".join([str(a.get("LastName","")) + " " + (str(a.get("ForeName",""))[:1]+"." if a.get("ForeName") else "") for a in art.get("AuthorList",[])[:3]])
                        journal = str(art["Journal"].get("Title", "Unknown"))
                        year = str(art["Journal"]["JournalIssue"]["PubDate"].get("Year", "N/A"))
                        pmid = str(p["MedlineCitation"]["PMID"])

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

CRITICAL: PATIENT-SPECIFIC ASSESSMENT PROTOCOL
This is NOT a template - tailor EVERY recommendation to THIS specific patient. Before finalizing:
1. **Individualize Risk Assessment**: Don't use generic statements. Calculate actual RCRI score, cite specific comorbidities.
2. **Personalize Medication Recommendations**: Base timing on THIS patient's medications, procedures, and comorbidities - not general lists.
3. **Specific Risk Quantification**: Provide actual risk percentages for THIS patient (cardiac event, respiratory complication, AKI) using RCRI/NSQIP framework.
4. **Airway Individualization**: Address THIS patient's airway (BMI, OSA, age, Mallampati if mentioned, previous anesthesia).
5. **Cross-Check**: Are all recommendations consistent with THIS patient's labs, EF, creatinine, comorbidities?
6. **Avoid Generic Advice**: Don't say "consider monitoring" - specify WHICH monitors for THIS patient and WHY.

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

Use HTML formatting:
- <h3>Section Headers</h3>
- <p>Paragraphs</p>
- <strong>Bold for emphasis</strong>
- <br><br> for spacing

IMPORTANT: Use inline citations [1], [2], [3] throughout your assessment to reference the papers provided above. Do NOT create a separate "Evidence-Based Citations" or "References" section - the references will be displayed separately below your assessment.

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

@app.route("/informed-consent", methods=["GET", "POST"])
def informed_consent_generator():
    """Informed Consent Generator with AI-powered patient-specific consent discussion"""
    if request.method == "GET":
        response = make_response(render_template_string(INFORMED_CONSENT_HTML, summary=None, references=None))
        response.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
        response.headers['Pragma'] = 'no-cache'
        response.headers['Expires'] = '0'
        return response

    # Collect and sanitize form data
    age = int(request.form.get("age", 0))
    asa = sanitize_user_query(request.form.get("asa", ""))
    anesthesia_type = sanitize_user_query(request.form.get("anesthesia_type", ""))
    procedure = sanitize_user_query(request.form.get("procedure", ""))
    comorbidities = request.form.getlist("comorbidities")
    special_considerations = request.form.getlist("special_considerations")
    notes = sanitize_user_query(request.form.get("notes", ""))

    # Validate required fields
    if not procedure or not procedure.strip():
        error_message = "<p style='color: #ff6b6b; text-align: center; padding: 20px;'><strong>Error:</strong> Please specify the surgical procedure before submitting the form.</p>"
        return render_template_string(INFORMED_CONSENT_HTML, summary=error_message, references=None)

    if not anesthesia_type:
        error_message = "<p style='color: #ff6b6b; text-align: center; padding: 20px;'><strong>Error:</strong> Please select an anesthesia type before submitting the form.</p>"
        return render_template_string(INFORMED_CONSENT_HTML, summary=error_message, references=None)

    # Build targeted PubMed searches based on patient factors
    search_queries = []

    # Primary search: procedure-specific anesthesia risks
    search_queries.append(f"{procedure} anesthesia risks complications consent")

    # Anesthesia type-specific risks
    if anesthesia_type == "General":
        search_queries.append("general anesthesia complications informed consent")
    elif anesthesia_type == "Regional/Neuraxial":
        search_queries.append("neuraxial anesthesia spinal epidural complications consent")
    elif anesthesia_type == "Sedation/MAC":
        search_queries.append("monitored anesthesia care sedation complications")

    # Comorbidity-specific searches
    if "Cardiac Disease" in comorbidities:
        search_queries.append(f"{procedure} cardiac disease perioperative risk")
    if "Pulmonary Disease" in comorbidities:
        search_queries.append(f"{procedure} pulmonary complications COPD anesthesia")
    if "Diabetes" in comorbidities:
        search_queries.append("perioperative diabetes management anesthesia")
    if "Renal Disease" in comorbidities:
        search_queries.append("chronic kidney disease anesthesia perioperative")
    if "OSA" in comorbidities:
        search_queries.append("obstructive sleep apnea anesthesia perioperative risks")

    # Special considerations
    if "Pregnancy" in special_considerations:
        search_queries.append("anesthesia pregnancy safety obstetric")
    if "Malignant Hyperthermia History" in special_considerations:
        search_queries.append("malignant hyperthermia anesthesia management")
    if "Difficult Airway" in special_considerations:
        search_queries.append("difficult airway management anesthesia consent")
    if "Chronic Pain" in special_considerations:
        search_queries.append(f"{procedure} chronic pain postoperative analgesia")
    if "Anticoagulation" in special_considerations:
        search_queries.append("anticoagulation perioperative management bleeding risk")

    # Search PubMed for all queries and collect papers
    all_refs = []
    all_context = ""

    for query in search_queries[:4]:  # Limit to 4 searches
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
                        title = str(art.get("ArticleTitle", "No title"))
                        abstract = " ".join(str(t) for t in art.get("Abstract", {}).get("AbstractText", [])) if art.get("Abstract") else ""
                        authors = ", ".join([str(a.get("LastName","")) + " " + (str(a.get("ForeName",""))[:1]+"." if a.get("ForeName") else "") for a in art.get("AuthorList",[])[:3]])
                        journal = str(art["Journal"].get("Title", "Unknown"))
                        year = str(art["Journal"]["JournalIssue"]["PubDate"].get("Year", "N/A"))
                        pmid = str(p["MedlineCitation"]["PMID"])

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

    # Sort references by quality
    unique_refs.sort(key=lambda x: x.get('sort_priority', 99))

    # Create numbered reference list for GPT
    ref_list = ""
    for i, ref in enumerate(unique_refs, 1):
        ref_list += f"[{i}] {ref['title']} - {ref['authors']} ({ref['year']}) PMID: {ref['pmid']}\n"

    # Build patient summary for GPT
    all_comorbidities = ', '.join(comorbidities) if comorbidities else 'None'
    all_special = ', '.join(special_considerations) if special_considerations else 'None'

    patient_data = f"""
Patient Information:
- Age: {age} years
- ASA Classification: ASA {asa}

Planned Anesthesia & Procedure:
- Anesthesia Type: {anesthesia_type}
- Surgical Procedure: {procedure}

Medical History:
- Comorbidities: {all_comorbidities}
- Special Considerations: {all_special}

Additional Notes: {notes if notes else 'None'}
"""

    # Generate GPT summary
    prompt = f"""You are an expert anesthesiologist helping create a comprehensive, patient-specific informed consent discussion guide. Based on the patient's information and evidence-based literature, provide a detailed consent discussion that balances medical accuracy with patient-friendly language.

Patient Information:
{patient_data}

Available Evidence (use numbered citations [1], [2], etc.):
{ref_list}

Paper Details:
{all_context}

CRITICAL: PATIENT-SPECIFIC CONSENT PROTOCOL
This is NOT a generic template - personalize EVERY section to THIS specific patient. Before finalizing:
1. **Individualize Risks**: Focus on risks relevant to THIS patient's age ({age}), ASA {asa} status, {anesthesia_type} anesthesia, and comorbidities ({all_comorbidities}).
2. **Procedure-Specific Information**: Tailor common and serious risks to {procedure} specifically, not general anesthesia.
3. **Evidence-Based Incidence Rates**: Provide actual risk percentages from literature when available (e.g., "PONV occurs in X% of patients undergoing {procedure}").
4. **Special Population Considerations**: If patient is elderly, pediatric, pregnant, or has special considerations ({all_special}), address these specifically.
5. **Cross-Check Comorbidities**: Are all comorbidity-specific risks addressed? (cardiac, pulmonary, renal, etc.)

Generate a comprehensive informed consent discussion guide including:

1. **Procedure Overview & Anesthesia Type**:
   - Brief explanation of {anesthesia_type} anesthesia for {procedure}
   - What the patient will experience before, during, and after
   - Why this anesthesia type is appropriate for their procedure

2. **Common Side Effects** (with frequencies when available):
   - List 5-7 most common side effects specific to {anesthesia_type} and {procedure}
   - Include incidence rates (e.g., "Nausea occurs in 20-30% of patients")
   - Explain duration and management strategies
   - Use patient-friendly language while maintaining medical accuracy

3. **Serious but Rare Risks** (with incidence rates):
   - List serious complications relevant to THIS patient's ASA {asa} status and comorbidities
   - Include specific incidence rates (e.g., "1 in X,XXX patients")
   - Explain risk factors that increase likelihood (age, comorbidities, procedure type)
   - Address {anesthesia_type}-specific serious risks (e.g., nerve injury for regional, aspiration for general)

4. **Patient-Specific Risk Factors**:
   - How THIS patient's comorbidities ({all_comorbidities}) affect anesthesia risks
   - Special precautions being taken for their conditions
   - Why their ASA {asa} classification matters
   - Management of special considerations: {all_special}

5. **Alternative Options**:
   - Alternative anesthesia techniques available for {procedure}
   - Pros and cons of each option for THIS patient
   - Why {anesthesia_type} is recommended

6. **Questions Patients Should Ask**:
   - 5-7 important questions this patient should discuss with their anesthesiologist
   - Questions specific to their comorbidities and special considerations
   - Questions about postoperative pain management, recovery, etc.

7. **What to Expect & Preparation**:
   - Pre-operative instructions (NPO status, medications to continue/hold)
   - What happens in the holding area and OR
   - Recovery room expectations
   - When they can eat, drink, resume activities

Use HTML formatting:
- <h3>Section Headers</h3>
- <p>Paragraphs</p>
- <strong>Bold for emphasis and risk percentages</strong>
- <ul><li>Bulleted lists for side effects and risks</li></ul>
- <br><br> for spacing

IMPORTANT: Use inline citations [1], [2], [3] throughout your discussion to reference the papers provided above. Do NOT create a separate "References" section - the references will be displayed separately below.

Tone: Professional but accessible - explain medical concepts in plain English while maintaining accuracy. The goal is informed consent that patients can truly understand."""

    try:
        response = openai_client.chat.completions.create(
            model="gpt-4o",
            messages=[{"role": "user", "content": prompt}],
            temperature=0.1
        ).choices[0].message.content
    except Exception as e:
        response = f"<p>Error generating consent discussion: {str(e)}</p>"

    return render_template_string(INFORMED_CONSENT_HTML, summary=response, references=unique_refs)

@app.route("/difficult-airway", methods=["GET", "POST"])
def difficult_airway_assessment():
    """Difficult Airway Assessment with AI-powered risk stratification"""
    if request.method == "GET":
        response = make_response(render_template_string(DIFFICULT_AIRWAY_HTML, summary=None, references=None))
        response.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
        response.headers['Pragma'] = 'no-cache'
        response.headers['Expires'] = '0'
        return response

    # Collect and sanitize form data
    age = int(request.form.get("age", 0))
    bmi = float(request.form.get("bmi", 0))
    mallampati = sanitize_user_query(request.form.get("mallampati", ""))
    thyromental = sanitize_user_query(request.form.get("thyromental", ""))
    mouth_opening = sanitize_user_query(request.form.get("mouth_opening", ""))
    neck_extension = sanitize_user_query(request.form.get("neck_extension", ""))
    risk_factors = request.form.getlist("risk_factors")
    procedure = sanitize_user_query(request.form.get("procedure", ""))
    case_type = sanitize_user_query(request.form.get("case_type", ""))
    notes = sanitize_user_query(request.form.get("notes", ""))

    # Validate required fields
    if not procedure or not procedure.strip():
        error_message = "<p style='color: #ff6b6b; text-align: center; padding: 20px;'><strong>Error:</strong> Please specify the planned procedure before submitting the form.</p>"
        return render_template_string(DIFFICULT_AIRWAY_HTML, summary=error_message, references=None)

    # Calculate risk score based on validated predictors
    risk_score = 0
    risk_factors_list = []

    # Mallampati III-IV (1 point)
    if mallampati in ["III", "IV"]:
        risk_score += 1
        risk_factors_list.append(f"Mallampati Class {mallampati}")

    # Thyromental distance <6cm (1 point)
    if thyromental == "<6cm":
        risk_score += 1
        risk_factors_list.append("Thyromental distance <6 cm")

    # Mouth opening <3cm (1 point)
    if mouth_opening == "<3cm":
        risk_score += 1
        risk_factors_list.append("Mouth opening <3 cm")

    # Limited neck extension (1 point)
    if neck_extension == "Limited":
        risk_score += 1
        risk_factors_list.append("Limited neck extension")

    # High BMI >35 (1 point)
    if bmi > 35:
        risk_score += 1
        risk_factors_list.append(f"Obesity (BMI {bmi:.1f} kg/m²)")

    # Previous difficult intubation (2 points - strong predictor)
    if "Previous Difficult Intubation" in risk_factors:
        risk_score += 2
        risk_factors_list.append("Previous difficult intubation")

    # OSA (1 point)
    if "OSA" in risk_factors:
        risk_score += 1
        risk_factors_list.append("Obstructive sleep apnea")

    # Other anatomical factors (0.5 points each)
    for factor in ["Beard", "Prominent Incisors", "Short Neck", "Large Tongue", "Facial Trauma", "C-Spine Pathology"]:
        if factor in risk_factors:
            risk_score += 0.5
            risk_factors_list.append(factor.lower())

    # Determine risk category
    if risk_score >= 4:
        risk_category = "High"
    elif risk_score >= 2:
        risk_category = "Moderate"
    else:
        risk_category = "Low"

    # Build targeted PubMed searches
    search_queries = []

    # Primary airway management search
    search_queries.append("difficult airway management guidelines anesthesia")

    # If high risk, search for advanced techniques
    if risk_category == "High":
        search_queries.append("awake fiberoptic intubation indications technique")
        search_queries.append("video laryngoscopy difficult airway")

    # OSA-specific management
    if "OSA" in risk_factors:
        search_queries.append("obstructive sleep apnea difficult airway perioperative")

    # Obesity-specific positioning
    if bmi > 35:
        search_queries.append("obese patient airway management positioning HELP")

    # Emergency airway if emergency case
    if case_type == "Emergency":
        search_queries.append("emergency airway management rapid sequence intubation")

    # Previous difficult intubation
    if "Previous Difficult Intubation" in risk_factors:
        search_queries.append("previous difficult intubation management anesthesia")

    # Search PubMed for all queries and collect papers
    all_refs = []
    all_context = ""

    for query in search_queries[:4]:  # Limit to 4 searches
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
                        title = str(art.get("ArticleTitle", "No title"))
                        abstract = " ".join(str(t) for t in art.get("Abstract", {}).get("AbstractText", [])) if art.get("Abstract") else ""
                        authors = ", ".join([str(a.get("LastName","")) + " " + (str(a.get("ForeName",""))[:1]+"." if a.get("ForeName") else "") for a in art.get("AuthorList",[])[:3]])
                        journal = str(art["Journal"].get("Title", "Unknown"))
                        year = str(art["Journal"]["JournalIssue"]["PubDate"].get("Year", "N/A"))
                        pmid = str(p["MedlineCitation"]["PMID"])

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

    # Sort references by quality
    unique_refs.sort(key=lambda x: x.get('sort_priority', 99))

    # Create numbered reference list for GPT
    ref_list = ""
    for i, ref in enumerate(unique_refs, 1):
        ref_list += f"[{i}] {ref['title']} - {ref['authors']} ({ref['year']}) PMID: {ref['pmid']}\n"

    # Build patient summary for GPT
    all_risk_factors = ', '.join(risk_factors_list) if risk_factors_list else 'None identified'

    patient_data = f"""
Patient Demographics:
- Age: {age} years
- BMI: {bmi:.1f} kg/m²

Airway Assessment:
- Mallampati Classification: Class {mallampati}
- Thyromental Distance: {thyromental}
- Mouth Opening: {mouth_opening}
- Neck Extension: {neck_extension}

Risk Factors Present: {all_risk_factors}

Risk Score: {risk_score}/7 points
Risk Category: {risk_category}

Procedure: {procedure}
Case Type: {case_type}

Additional Notes: {notes if notes else 'None'}
"""

    # Generate GPT summary
    prompt = f"""You are an expert anesthesiologist performing a comprehensive difficult airway assessment. Based on validated clinical predictors and recent evidence, provide a detailed risk stratification and management plan.

Patient Information:
{patient_data}

Available Evidence (use numbered citations [1], [2], etc.):
{ref_list}

Paper Details:
{all_context}

CRITICAL: PATIENT-SPECIFIC AIRWAY ASSESSMENT PROTOCOL
This is NOT a template - tailor EVERY recommendation to THIS specific patient. Before finalizing:
1. **Individualize Risk Assessment**: Calculate the actual difficult airway risk score based on THIS patient's Mallampati, thyromental distance, mouth opening, neck extension, BMI, and history.
2. **Specific Management Strategy**: Base airway plan on THIS patient's risk category ({risk_category}), case urgency ({case_type}), and specific anatomical challenges.
3. **Evidence-Based Equipment Selection**: Recommend specific airway devices and techniques appropriate for THIS patient's risk factors.
4. **Personalized Positioning**: If patient is obese (BMI {bmi:.1f}), specify HELP position or ramping details.
5. **Cross-Check Case Urgency**: If {case_type} case, adjust awake vs asleep intubation recommendations accordingly.

Generate a comprehensive difficult airway assessment including:

1. **Risk Stratification Summary**:
   - Overall risk category ({risk_category}) with justification
   - Specific predictors identified in this patient
   - Estimated difficulty level for mask ventilation, supraglottic airway, and intubation

2. **Pre-Intubation Optimization**:
   - Patient positioning (HELP position if BMI >35, ramping if needed)
   - Pre-oxygenation strategy (target EtO₂, apneic oxygenation via nasal cannula)
   - Equipment preparation (specific devices: video laryngoscope type, LMA size, bougie, etc.)
   - Personnel needs (extra anesthesiologist, surgeon availability, etc.)

3. **Primary Airway Management Plan**:
   - Recommended primary technique (awake fiberoptic vs video laryngoscopy vs direct laryngoscopy)
   - Drug selection and dosing if asleep intubation planned
   - Maximum number of attempts before escalation
   - Specific technique modifications for this patient's anatomy

4. **Backup Plans (ASA Difficult Airway Algorithm 2022)**:
   - Plan A: Primary intubation approach
   - Plan B: Alternative intubation technique if Plan A fails
   - Plan C: Supraglottic airway rescue (specific device and size)
   - Plan D: Emergency front-of-neck access preparation (cricothyrotomy kit location)
   - Decision point: When to awaken patient vs proceed with emergency pathway

5. **Case-Specific Considerations**:
   - Special considerations for {case_type} case
   - Management of {', '.join(risk_factors) if risk_factors else 'patient-specific factors'}
   - Postoperative extubation planning
   - Communication plan with surgical team

Use HTML formatting:
- <h3>Section Headers</h3>
- <p>Paragraphs</p>
- <strong>Bold for emphasis</strong>
- <ul><li>Bulleted lists</li></ul>
- <br><br> for spacing

IMPORTANT: Use inline citations [1], [2], [3] throughout your assessment to reference the papers provided above. Do NOT create a separate "References" section - the references will be displayed separately below your assessment.

Provide maximum clinical utility with specific, actionable recommendations backed by evidence. This assessment should be directly usable for safe airway management planning."""

    try:
        response = openai_client.chat.completions.create(
            model="gpt-4o",
            messages=[{"role": "user", "content": prompt}],
            temperature=0.1
        ).choices[0].message.content
    except Exception as e:
        response = f"<p>Error generating assessment: {str(e)}</p>"

    return render_template_string(DIFFICULT_AIRWAY_HTML, summary=response, references=unique_refs)

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

    # Machine Learning-Based Prediction System
    # Load trained RandomForest model and make prediction

    try:
        import pickle
        import numpy as np

        # Load model and scaler
        with open('ioh_model.pkl', 'rb') as f:
            model = pickle.load(f)
        with open('ioh_scaler.pkl', 'rb') as f:
            scaler = pickle.load(f)

        # Encode categorical variables
        sex_encoded = 1 if sex == "male" else 0

        vasopressor_map = {"none": 0, "phenylephrine": 1, "ephedrine": 2, "norepinephrine": 3}
        vasopressor_encoded = vasopressor_map.get(vasopressor, 0)

        surgery_type_map = {"minor": 0, "moderate": 1, "major_abdominal": 2, "cardiac": 3, "vascular": 4}
        surgery_type_encoded = surgery_type_map.get(surgery_type, 0)

        induction_agent_map = {"propofol": 0, "etomidate": 1, "ketamine": 2}
        induction_agent_encoded = induction_agent_map.get(induction_agent, 0)

        emergency_encoded = 1 if emergency == "yes" else 0

        # Prepare features in same order as training
        # Features: age, sex, bmi, asa, baseline_map, baseline_hr, current_map, map_5min, map_10min,
        #           surgery_duration, vasopressor, surgery_type, induction_agent, emergency
        features = np.array([[
            age, sex_encoded, bmi, asa, baseline_map, baseline_hr,
            current_map, map_5min, map_10min, surgery_duration,
            vasopressor_encoded, surgery_type_encoded, induction_agent_encoded, emergency_encoded
        ]])

        # Scale features
        features_scaled = scaler.transform(features)

        # Get prediction probability
        ioh_prob = model.predict_proba(features_scaled)[0][1]  # Probability of IOH

        # Convert to percentage
        prob_5min = int(ioh_prob * 100)

        # Estimate probabilities for 10 and 20 minute windows (with decay)
        prob_10min = int(ioh_prob * 0.85 * 100)
        prob_20min = int(ioh_prob * 0.70 * 100)

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

        # Identify top risk factors based on feature analysis
        factors = []

        # MAP trend analysis
        map_trend = (current_map - map_5min) + (current_map - map_10min) / 2
        if map_trend < -10:
            factors.append({
                "name": "Declining MAP Trend",
                "description": f"MAP decreased by {abs(int(map_trend))} mmHg over last 10 minutes, indicating hemodynamic instability"
            })
        elif map_trend < -5:
            factors.append({
                "name": "Moderate MAP Decline",
                "description": f"MAP decreased by {abs(int(map_trend))} mmHg, suggesting early hemodynamic changes"
            })

        # Current MAP level
        if current_map < 70:
            factors.append({
                "name": "Low Current MAP",
                "description": f"Current MAP of {current_map} mmHg is approaching hypotension threshold (65 mmHg)"
            })

        # Baseline deviation
        map_drop_pct = ((baseline_map - current_map) / baseline_map) * 100 if baseline_map > 0 else 0
        if map_drop_pct > 20:
            factors.append({
                "name": "Significant MAP Drop from Baseline",
                "description": f"{int(map_drop_pct)}% decrease from baseline MAP, exceeding critical threshold"
            })

        # Age risk
        if age > 70:
            factors.append({
                "name": "Advanced Age",
                "description": f"Age {age} years associated with increased cardiovascular lability"
            })

        # ASA class
        if asa >= 3:
            factors.append({
                "name": "High ASA Classification",
                "description": f"ASA {asa} indicates significant comorbidities affecting hemodynamic stability"
            })

        # Surgery type
        high_risk_surgeries = ["cardiac", "major_abdominal", "vascular"]
        if surgery_type in high_risk_surgeries:
            factors.append({
                "name": "High-Risk Surgery Type",
                "description": f"{surgery_type.replace('_', ' ').title()} surgery associated with greater fluid shifts"
            })

        # Emergency
        if emergency == "yes":
            factors.append({
                "name": "Emergency Surgery",
                "description": "Emergency procedures have higher hypotension risk"
            })

        # Induction agent
        if induction_agent == "propofol":
            factors.append({
                "name": "Propofol Induction",
                "description": "Propofol associated with dose-dependent vasodilation"
            })

        # Vasopressor use
        if vasopressor != "none":
            factors.append({
                "name": "Current Vasopressor Requirement",
                "description": f"Ongoing {vasopressor} use indicates hemodynamic instability"
            })

        # Keep top 3 factors
        factors = factors[:3] if len(factors) > 3 else factors

        # If no factors, add default
        if not factors:
            factors.append({
                "name": "Stable Hemodynamics",
                "description": "No major risk factors identified based on current parameters"
            })

    except Exception as e:
        logger.error(f"ML model prediction error: {str(e)}")
        # Fallback to simple heuristic if model fails
        prob_5min = 30
        prob_10min = 25
        prob_20min = 20
        risk_5min_class, risk_5min_label, risk_5min_text = "moderate", "risk-moderate", "Moderate Risk"
        risk_10min_class, risk_10min_label, risk_10min_text = "low", "risk-low", "Low Risk"
        risk_20min_class, risk_20min_label, risk_20min_text = "low", "risk-low", "Low Risk"
        factors = [{
            "name": "Model Unavailable",
            "description": "ML model could not be loaded. Using fallback estimation."
        }]

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