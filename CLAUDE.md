# CLAUDE.md - AI Assistant Guide for gasconsult.ai

## Project Overview

**gasconsult.ai** is an evidence-based anesthesiology AI consultant that combines PubMed medical literature search with GPT-4o to provide hallucination-free, citation-backed clinical answers.

### Core Functionality
- Accepts clinical questions related to anesthesiology
- Searches PubMed for high-quality evidence (systematic reviews, meta-analyses, RCTs, Cochrane reviews, guidelines)
- Synthesizes answers using GPT-4o with real paper citations
- Returns structured response with answer + full reference list

## Technology Stack

| Component | Technology | Version |
|-----------|------------|---------|
| Backend | Flask | 3.0.3 |
| Medical Literature API | Biopython (Entrez) | 1.83 |
| AI/LLM | OpenAI (GPT-4o) | 1.53.0 |
| Production Server | Gunicorn | 22.0.0 |
| HTTP Client | httpx | 0.27.2 |
| Security | bleach, Flask-WTF, Flask-Limiter | 6.1.0, 1.2.1, 3.5.0 |
| Configuration | python-dotenv | 1.0.0 |
| Language | Python 3.x | - |

## Directory Structure

```
gasconsult/
├── app.py                  # Main Flask application (~8450 lines)
├── requirements.txt        # Python dependencies
├── .env.example           # Environment variable template
├── README.md              # Project description
├── CLAUDE.md              # This file - AI assistant guide
├── CHANGELOG.md           # Version history and changes
├── DEPLOYMENT_GUIDE.md    # Production deployment instructions
└── static/
    ├── favicon.svg        # Browser favicon
    ├── logo.png           # Site logo
    ├── manifest.json      # PWA manifest
    └── sw.js             # Service worker
```

## Architecture

### Single-File Application (`app.py`)

The entire application is contained in `app.py` with these components:

1. **Flask App Initialization** (lines 1-11)
   - Flask setup, Entrez email/API key, OpenAI client

2. **HTML Template** (lines 13-73)
   - Inline HTML with embedded CSS
   - Dark theme UI with Inter font
   - Form for query input, response/references display sections

3. **Main Route Handler** (`/`) (lines 75-154)
   - Query preprocessing with synonym expansion
   - PubMed search with fallback strategy
   - Paper metadata extraction
   - GPT-4o prompt construction and response generation

### Query Processing Flow

```
User Query → Synonym Expansion → PubMed Search (anesthesiology first, then general)
     → Paper Fetch → Metadata Extraction → GPT-4o Synthesis → Response + References
```

### Search Strategy

The app uses a two-tier search approach:
1. First attempts anesthesiology-specific search with MeSH terms
2. Falls back to general search if no results found
3. Filters for high-quality evidence types only (since 2015)

### Backend Intelligence Features (Added 2025-11-28)

The chat function includes 8 sophisticated query processing features that make it understand clinical questions better:

#### 1. **Expanded Medical Abbreviation Dictionary**
Automatically recognizes and expands 40+ medical abbreviations:
- **Drugs:** `roc` → rocuronium, `vec` → vecuronium, `sux` → succinylcholine, `dex` → dexmedetomidine, `neo` → phenylephrine, `epi` → epinephrine, `remi` → remifentanil, `versed` → midazolam
- **Procedures:** `RSI` → rapid sequence induction, `CABG` → coronary artery bypass, `TAVR` → transcatheter aortic valve
- **Complications:** `LAST` → local anesthetic systemic toxicity, `PRIS` → propofol infusion syndrome, `MH` → malignant hyperthermia
- **Regional:** `nerve block`, `epidural`, `spinal` → comprehensive neuraxial terms
- **Original terms:** `txa`, `ponv`, `peds`, `propofol`, `etomidate`, `ketamine`, `blood loss`, `spine surgery`

#### 2. **Pronoun & Reference Resolution**
Intelligently resolves pronouns and vague references in follow-up questions:
- "What about in pediatrics?" → Carries forward previous drug/topic
- "What are the contraindications for it?" → Replaces "it" with actual entity from context
- "How about in renal failure?" → Continues previous topic with new modifier
- Extracts medical entities from last 6 messages for context

#### 3. **Question Type Detection**
Classifies questions and customizes search strategy:
- **Dosing** → Prioritizes guidelines and review articles (e.g., "What's the dose of propofol?")
- **Safety** → Includes adverse effects and contraindication studies (e.g., "Is ketamine safe in...?")
- **Comparison** → Targets head-to-head trials (e.g., "Propofol vs etomidate")
- **Mechanism** → Focuses on review articles (e.g., "How does dexmedetomidine work?")
- **Management** → Searches for protocols and guidelines (e.g., "Difficult airway algorithm")

#### 4. **Negation Handling**
Detects contraindication/avoidance questions and modifies search:
- Detects keywords: "avoid", "not", "contraindication", "when to stop", "shouldn't use"
- Adds search filters: `contraindications[sh]`, `adverse effects[sh]`, `safety[ti]`
- Updates GPT prompt to focus on safety concerns and inappropriate use cases
- Example: "When NOT to use succinylcholine?" → Gets malignant hyperthermia, hyperkalemia risks

#### 5. **Dynamic Temperature Adjustment**
Adjusts GPT-4o temperature based on question criticality:
- **Dosing questions:** 0.05 (ultra-precise - dosing must be exact)
- **Safety questions:** 0.1 (factual accuracy critical)
- **Mechanism questions:** 0.15 (can be more explanatory)
- **Comparison/Management:** 0.1 (balanced)
- **General/No papers:** 0.2 (conversational)

#### 6. **Multi-Part Question Detection**
Identifies questions with multiple parts and ensures comprehensive answers:
- Detects patterns: "and what", "and how", multiple "?"
- Signals GPT to address each part thoroughly
- Example: "What's the dose of propofol and what are the side effects?" → Both parts answered

#### 7. **Smart Context Carryover**
Builds intelligent conversation context instead of just recent messages:
- Analyzes current query for key terms
- Finds earlier messages (beyond last 6) that share 2+ terms with current query
- Includes up to 2 relevant earlier Q&A pairs + last 6 messages
- Maintains continuity in long conversations across different topics

#### 8. **Conversation Topic Tracking**
Remembers the main topic of conversation:
- Extracts main topic from first message (keyword-based)
- Stores in session for potential query enhancement
- Available for expanding vague follow-up queries

### UX Features

#### Keyboard Shortcuts
- **Ctrl+Enter** (Windows/Linux) or **Cmd+Enter** (Mac) to submit queries
- Works on both homepage and chat interface
- Prevents accidental newlines while typing

#### Evidence Quality Badges
- Shows confidence level (High/Moderate/Low) based on paper quantity and types
- Displays number of papers, study types, and date range
- Color-coded: Green (high), Orange (moderate), Red (low)
- Transparent about evidence strength

#### Suggested Prompts
- Homepage shows 5 clickable example queries
- Auto-fills textarea when clicked
- Helps new users understand what to ask

### Synonym Expansion (Legacy - Integrated into Abbreviation Dictionary)

The original synonym expansion is now part of the comprehensive medical abbreviation dictionary (see Backend Intelligence Features #1 above).

## Environment Variables

| Variable | Required | Description |
|----------|----------|-------------|
| `OPENAI_API_KEY` | Yes | OpenAI API key for GPT-4o access |
| `ENTREZ_EMAIL` | Yes | Email for NCBI Entrez API |
| `ENTREZ_API_KEY` | Recommended | NCBI API key (increases rate limit) |
| `FLASK_SECRET_KEY` | Recommended | Secret key for session management |
| `FLASK_ENV` | Optional | `development` or `production` (default: production) |
| `LOG_LEVEL` | Optional | Logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL |
| `LOG_FILE` | Optional | Path to log file (empty = stdout only) |
| `RATE_LIMIT` | Optional | Rate limit (default: "60 per minute") |

**Setup:** Copy `.env.example` to `.env` and fill in your values.

## Development Commands

### Install Dependencies
```bash
pip install -r requirements.txt
```

### Run Development Server
```bash
python app.py
# Or with Flask CLI:
flask run --debug
```

### Run Production Server
```bash
gunicorn app:app -b 0.0.0.0:8000
```

## Code Conventions

### Style Guidelines
- Single-file Flask architecture (keep it simple)
- Inline HTML template with embedded CSS
- Use render_template_string for template rendering
- Minimal error handling with try/except blocks

### Naming Conventions
- Route handlers: lowercase function names (`index`)
- Variables: snake_case (`search_term`, `num_papers`)
- Constants: UPPERCASE (`HTML`)

### Error Handling Pattern
```python
try:
    # API call
except Exception as e:
    response = f"Error message: {str(e)}"
```

## Key Dependencies & APIs

### PubMed/Entrez API
- **Documentation:** https://www.ncbi.nlm.nih.gov/books/NBK25500/
- **Rate Limits:** 10 requests/second with API key
- **Methods Used:**
  - `Entrez.esearch()` - Search for paper IDs
  - `Entrez.efetch()` - Fetch paper metadata
  - `Entrez.read()` - Parse XML responses

### OpenAI API
- **Model:** gpt-4o
- **Temperature:** 0.1 (low for factual responses)
- **Usage:** Single chat completion call with structured prompt

## Important Considerations for AI Assistants

### When Modifying Code

1. **Preserve the single-file architecture** - This app is intentionally simple
2. **Keep inline HTML/CSS** - No need to extract to separate files
3. **Maintain synonym expansion** - Critical for search quality
4. **Preserve two-tier search strategy** - Ensures relevant results
5. **Keep low temperature (0.1)** - Essential for factual medical content

### Security Features (Production-Ready ✓)

- ✅ **All credentials in environment variables** - No hardcoded secrets
- ✅ **Input sanitization** - XSS protection using bleach library
- ✅ **CSRF protection** - All POST requests protected
- ✅ **Rate limiting** - 60 requests/minute per IP (configurable)
- ✅ **Server-side sessions** - No sensitive data in cookies
- ✅ **Comprehensive logging** - Request tracking and error monitoring
- ✅ **Health check endpoint** - `/health` for monitoring
- ✅ **HTTPS recommended** - Use reverse proxy with SSL

### Security Architecture

1. **Input Validation Layer**
   - `sanitize_user_query()` - Strips all HTML from user input
   - `sanitize_input()` - Allows safe HTML tags for GPT responses
   - Applied to all form inputs (chat, preop, hypotension)

2. **Rate Limiting**
   - Flask-Limiter with in-memory storage
   - Per-IP address tracking
   - Configurable limits via environment

3. **CSRF Protection**
   - Flask-WTF automatic token validation
   - Exempt endpoints: `/health`, `/api/status`

4. **Session Security**
   - Server-side filesystem storage
   - Signed sessions
   - Temporary conversation history only

### Testing Queries

Good test queries for the application:
- "TXA in spine surgery"
- "blood loss scoliosis"
- "propofol peds"
- "PONV prevention"

### Application Routes

| Route | Method | Description |
|-------|--------|-------------|
| `/` | GET | Homepage with search interface |
| `/chat` | GET/POST | Conversational AI chat interface |
| `/stream` | GET | Server-sent events for streaming responses |
| `/clear` | GET | Clear chat session |
| `/terms` | GET | Terms of Service |
| `/privacy` | GET | Privacy Policy |
| `/quick-dose` | GET | Drug dosing calculator |
| `/preop` | GET/POST | Pre-operative assessment tool |
| `/hypotension` | GET/POST | Intraoperative hypotension predictor |
| `/health` | GET | Health check endpoint (monitoring) |
| `/api/status` | GET | API status and endpoints list |

### Future Improvements

- ✅ ~~Move Entrez credentials to environment variables~~ (Completed)
- ✅ ~~Add input sanitization before rendering~~ (Completed)
- ✅ ~~Add logging for debugging~~ (Completed)
- ✅ ~~Consider rate limiting~~ (Completed)
- ⬜ Implement caching for repeated PubMed queries (Redis)
- ⬜ Add automated testing suite
- ⬜ Set up CI/CD pipeline
- ⬜ Enhanced medical disclaimers with acknowledgment
- ⬜ User accounts and saved queries (optional)
- ⬜ Mobile app (optional)

## Git Workflow

- Main development branch varies by feature
- Commit messages should be concise version tags (e.g., "v2.1")
- Recent commits indicate active development of search functionality

## Disclaimer

This application provides **educational content only** and is **not medical advice**. All outputs should be verified by qualified medical professionals before clinical use.
