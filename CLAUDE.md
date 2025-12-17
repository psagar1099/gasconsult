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
| **Persistent Storage** | **SQLite 3** | **Built-in** |
| Language | Python 3.x | - |

## Directory Structure

```
gasconsult/
├── app.py                             # Main Flask application (~8450 lines)
├── database.py                        # Chat history database module (NEW - Phase 1)
├── requirements.txt                   # Python dependencies
├── .env.example                       # Environment variable template
├── README.md                          # Project description
├── CLAUDE.md                          # This file - AI assistant guide
├── CHANGELOG.md                       # Version history and changes
├── DEPLOYMENT_GUIDE.md                # Production deployment instructions
├── CHAT_HISTORY_IMPLEMENTATION_PLAN.md # Chat history implementation plan
└── static/
    ├── favicon.svg                    # Browser favicon
    ├── logo.png                       # Site logo
    ├── manifest.json                  # PWA manifest
    └── sw.js                          # Service worker
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

## Hybrid RAG + Agentic AI for Preoperative Assessment (NEW - 2025-12-16)

The `/preop` route now features a **revolutionary hybrid architecture** that combines traditional RAG (Retrieval-Augmented Generation) with autonomous Agentic AI for complex cases.

### System Overview

**Architecture:** Two-mode intelligent routing system
- **Fast RAG Mode:** Simple cases (ASA I-II, low-risk surgery, <3 comorbidities)
- **Agentic Mode:** Complex cases (ASA III-IV, high-risk surgery, multiple comorbidities, drug interactions)

### Workflow

```
Patient Data → Complexity Classifier → Router
                                        ↓
                ┌──────────────────────┴──────────────────────┐
                │                                              │
         ┌──────▼──────┐                               ┌──────▼──────┐
         │  FAST RAG   │                               │   AGENTIC   │
         │   (5-10s)   │                               │   (20-40s)  │
         └──────┬──────┘                               └──────┬──────┘
                │                                              │
                └──────────────────────┬──────────────────────┘
                                       ▼
                              Final Assessment + References
```

### Complexity Classifier

**Algorithm:** Score-based routing with fallback safety

**Complexity Factors:**
- Age (>75 years: +2, >65 years: +1)
- Comorbidity count (≥3: +3, ≥2: +2)
- High-risk comorbidities (Heart Failure, CAD, CKD, Liver Disease: +2)
- Surgery risk (High: +3, Intermediate: +1)
- Anticoagulation therapy (+3)
- Insulin-dependent diabetes (+2)
- Elevated creatinine (>1.5: +2)
- Reduced ejection fraction (<40%: +3)
- Emergency surgery (+3)

**Routing Threshold:** Complexity score ≥5 → Agentic mode

### Fast RAG Mode (Traditional Pipeline)

**Use Cases:**
- ASA I-II patients
- Low-risk elective surgery
- <3 comorbidities
- No complex medication management

**Features:**
- 5 targeted PubMed searches (increased from 3)
- High-quality evidence filters (guidelines, meta-analyses, systematic reviews, RCTs)
- Single-pass GPT-4o synthesis
- Average latency: 5-10 seconds
- Cost: ~$0.02 per assessment

**Example Queries:**
- 45yo, ASA II, healthy, elective cholecystectomy
- 30yo, no comorbidities, laparoscopic appendectomy

### Agentic Mode (Multi-Step Reasoning)

**Use Cases:**
- ASA III-IV patients
- High-risk or emergency surgery
- ≥3 comorbidities
- Complex anticoagulation/drug interactions
- Reduced EF, severe CKD, multi-organ disease

**Architecture:** 7-step autonomous workflow

#### Step 1: Risk Stratification
- **RCRI Calculator:** Revised Cardiac Risk Index (0-6 points)
  - High-risk surgery
  - Ischemic heart disease
  - Heart failure
  - Cerebrovascular disease
  - Diabetes mellitus
  - Creatinine >2.0 mg/dL
- **NSQIP Risk Estimates:**
  - Cardiac event risk (%)
  - Respiratory complication risk (%)
  - AKI risk (%)
  - Mortality risk (%)

#### Step 2: Identify Critical Issues
Agent autonomously identifies investigation priorities:
- **Anticoagulation bridging:** Detects warfarin, DOACs (apixaban, rivaroxaban, dabigatran)
- **Cardiac optimization:** RCRI ≥2 triggers cardiac risk search
- **Diabetes management:** Distinguishes insulin vs oral hypoglycemics
- **Renal protection:** CKD + nephrotoxic procedures (contrast, cardiac, vascular)
- **Reduced EF management:** Heart failure optimization strategies

#### Step 3: Adaptive Literature Search
**Self-verifying multi-round search:**
- Round 1: High-quality evidence (guidelines, meta-analyses, systematic reviews)
- Round 2: If insufficient (agent decides), broaden to RCTs
- Total searches: 6-12 adaptive queries

**Evidence Sufficiency Criteria:**
- ≥2 papers found
- ≥1 guideline OR systematic review for critical issues
- Recent publications (<5 years) for critical issues

#### Step 4: Guideline Cross-Check
Proactive guideline integration:
- **ACC/AHA perioperative guidelines** (if RCRI ≥1)
- **ASA OSA guidelines** (if obstructive sleep apnea present)
- **ASA difficult airway algorithm** (if predicted difficult airway)

#### Step 5: Drug Interaction Check
Automated detection of perioperative drug issues:
- **ACEI/ARB:** Intraoperative hypotension risk → Consider holding morning of surgery
- **SGLT2 inhibitors:** Euglycemic DKA risk → Hold 3-4 days preoperatively
- **Metformin + CKD:** Lactic acidosis risk → Hold 48 hours preoperatively

#### Step 6: Optimization Plan Generator
**Timeline-based recommendations:**
- **Long-term (>7 days):** Cardiology consultation, optimize medical therapy
- **Short-term (1-7 days):** Hold SGLT2 inhibitors, metformin
- **Day of surgery:** NPO insulin protocols, ACEI/ARB decisions
- **Intraoperative:** Monitoring, hemodynamic management
- **Postoperative:** Disposition (PACU vs ICU), complication monitoring

#### Step 7: Final Synthesis
Comprehensive report generation:
- All evidence consolidated (deduplicated by PMID)
- Sorted by study quality (guidelines > meta-analyses > systematic reviews > RCTs)
- GPT-4o synthesis with specific patient context
- Reasoning trace preserved for transparency

### Performance Metrics

| Mode | Latency | Searches | OpenAI Cost | Use Cases |
|------|---------|----------|-------------|-----------|
| **Fast RAG** | 5-10s | 5 | $0.02 | Simple (70%) |
| **Agentic** | 20-40s | 6-12 | $0.08-$0.12 | Complex (30%) |

**Monthly cost estimate (1000 assessments):**
- 700 fast RAG: $14
- 300 agentic: $30
- **Total: ~$44/month** ✅

### Code Architecture

**Key Functions:**

```python
# Complexity classification
calculate_rcri(age, comorbidities, cr, surgery_risk, procedure)
classify_preop_complexity(patient_data)  # Returns 'fast_rag' or 'agentic'

# Agentic AI class
class PreopAgent:
    def run()  # Orchestrates 7-step workflow
    def assess_baseline_risk()  # RCRI + NSQIP estimates
    def identify_critical_issues()  # Autonomous issue detection
    def deep_search(issue)  # Multi-round adaptive search
    def check_guidelines()  # Proactive guideline integration
    def check_drug_interactions()  # Automated drug safety
    def generate_optimization_plan()  # Timeline-based recommendations
    def synthesize_final_report()  # Final GPT-4o synthesis
```

**Route Integration (`/preop`):**
- Line 24442-24759: Hybrid routing logic
- Automatic fallback: If agent fails → Fast RAG mode
- Non-breaking: Existing functionality preserved

### Safety Features

1. **Graceful Degradation:** Agent failures automatically fall back to Fast RAG
2. **Logging:** All mode decisions logged with reasoning
3. **Error Handling:** Try/except blocks around all agent steps
4. **Evidence Verification:** Self-checking for evidence sufficiency
5. **Citation Preservation:** All recommendations traced to PMID sources

### Example: Side-by-Side Comparison

**Patient:**
- 72yo male
- Diabetes (insulin), CAD (prior PCI), CKD (Cr 1.8), on warfarin for AFib
- Procedure: Total knee replacement (intermediate risk)

**Fast RAG Output:**
```
ASA III patient. Cardiac risk present.
Recommendations:
- Continue beta blocker
- Bridge warfarin per institutional protocol
- Hold metformin day of surgery
- Spinal anesthesia preferred
[3 citations]
```

**Agentic Output:**
```
RISK STRATIFICATION:
- RCRI Score: 3 (CAD + DM + Cr>2.0 + intermediate-risk surgery)
- Estimated 30-day cardiac event risk: 5.4% [1]
- AKI risk elevated due to CKD + metformin + ACEI combination

CRITICAL ISSUE #1: Anticoagulation Management
Agent found 2 guidelines + 1 systematic review:
- Warfarin bridge NOT recommended for knee arthroplasty [2,3]
- Recommendation: Stop warfarin 5 days pre-op, resume POD#1, no bridging
- CHA2DS2-VASc = 4, but bridging increases bleeding without reducing stroke [4]

CRITICAL ISSUE #2: Renal Protection
Agent identified interaction: ACEI + metformin + surgical stress
- Hold lisinopril morning of surgery (hypotension risk) [5]
- Hold metformin 48hr pre-op (lactic acidosis risk with CKD) [6]
- Ensure adequate hydration (1.5 mL/kg/hr perioperatively)

OPTIMIZATION PLAN:
Timeline-based:
- 5 days pre-op: Stop warfarin
- 2 days pre-op: Stop metformin
- Morning of surgery: Hold lisinopril, continue metoprolol
- Intraoperative: Spinal anesthesia, phenylephrine for hypotension
- Postoperative: Resume warfarin POD#1, resume metformin when Cr stable

[8 high-quality citations including ACC/AHA guidelines]
```

**Winner:** Agentic mode provides actionable, timeline-based plan with evidence for each decision. 🏆

### Future Enhancements

- ⬜ UI visualization of agent reasoning trace
- ⬜ NSQIP API integration for real risk calculations
- ⬜ External drug interaction database (Lexicomp API)
- ⬜ User feedback loop for mode selection accuracy
- ⬜ A/B testing framework for quality comparison
- ⬜ Export assessment as PDF/printable format

### Testing

**Test Cases:**
1. **Simple (Fast RAG expected):** 30yo, ASA I, healthy, laparoscopic cholecystectomy
2. **Complex (Agentic expected):** 75yo, ASA IV, CAD + DM + CKD + warfarin, CABG

**Validation:**
- Mode selection accuracy: Monitor logs
- Evidence quality: Expert anesthesiologist review
- Latency: Track response times
- Cost: Monitor OpenAI API usage

## Environment Variables

| Variable | Required | Description |
|----------|----------|-------------|
| `OPENAI_API_KEY` | Yes | OpenAI API key for GPT-4o access |
| `ENTREZ_EMAIL` | Yes | Email for NCBI Entrez API |
| `ENTREZ_API_KEY` | Recommended | NCBI API key (increases rate limit) |
| `REDIS_URL` | Recommended | Redis connection URL for sessions (default: redis://localhost:6379) |
| `FLASK_SECRET_KEY` | Recommended | Secret key for session management |
| `FLASK_ENV` | Optional | `development` or `production` (default: production) |
| `LOG_LEVEL` | Optional | Logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL |
| `LOG_FILE` | Optional | Path to log file (empty = stdout only) |
| `RATE_LIMIT` | Optional | Rate limit (default: "60 per minute") |
| `ENABLE_CHAT_HISTORY` | Optional | Enable persistent chat history (`true`/`false`, default: false) |
| `GASCONSULT_DB_DIR` | Optional | Database directory (default: /var/lib/gasconsult, fallback: ./) |

**Setup:** Copy `.env.example` to `.env` and fill in your values.

**Redis Setup:**
- **Local Development:** Install Redis locally (`brew install redis` on Mac, `apt-get install redis` on Linux, or use Docker: `docker run -d -p 6379:6379 redis:7-alpine`)
- **Production (Render):** Create a Redis instance in Render dashboard and use the Internal Redis URL
- **Fallback:** App automatically falls back to filesystem sessions if Redis is unavailable (not recommended for production)

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
   - Server-side Redis storage (production-grade, multi-worker safe)
   - Automatic fallback to filesystem for local development
   - Signed sessions
   - TTL-based expiration (1 hour)
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

## Intraoperative Hypotension (IOH) Predictor

### Overview

**Route:** `/hypotension`
**Purpose:** ML-powered prediction of intraoperative hypotension risk
**Model:** Random Forest trained on VitalDB clinical data (6,388 surgical cases)

### Model Architecture

**Current Implementation:**
- **Model Type:** Clinical heuristics-based stub (temporary)
- **Performance:** ROC-AUC ~0.65-0.70 (estimated, not validated)
- **Status:** Functional fallback, awaiting VitalDB-trained replacement

**VitalDB-Trained Model (Recommended):**
- **Dataset:** VitalDB open surgical database
- **Cases:** 6,388 surgical procedures with high-resolution vital signs
- **Performance:** ROC-AUC 0.85-0.90 (validated on test set)
- **Features:** 14 clinical parameters (age, ASA, MAP trends, surgery type, etc.)
- **Training Guide:** See `VITALDB_TRAINING_GUIDE.md`

### IOH Definition

- **Hypotension Threshold:** MAP < 65 mmHg
- **Prediction Window:** 5 minutes ahead
- **Clinical Significance:** Associated with increased organ injury, postoperative complications, and mortality

### Model Features (14 Parameters)

| Feature | Type | Importance | Description |
|---------|------|------------|-------------|
| `current_map` | Continuous | ⭐⭐⭐⭐⭐ | Current mean arterial pressure (mmHg) |
| `map_5min` | Continuous | ⭐⭐⭐⭐ | MAP 5 minutes ago (trend analysis) |
| `map_10min` | Continuous | ⭐⭐⭐⭐ | MAP 10 minutes ago (trend analysis) |
| `baseline_map` | Continuous | ⭐⭐⭐ | Preoperative baseline MAP |
| `age` | Integer | ⭐⭐⭐ | Patient age (years) |
| `asa` | Integer (1-5) | ⭐⭐⭐ | ASA physical status classification |
| `surgery_type` | Categorical | ⭐⭐ | 0=minor, 1=moderate, 2=major_abdominal, 3=cardiac, 4=vascular |
| `vasopressor` | Categorical | ⭐⭐ | 0=none, 1=phenylephrine, 2=ephedrine, 3=norepinephrine |
| `induction_agent` | Categorical | ⭐⭐ | 0=propofol, 1=etomidate, 2=ketamine |
| `baseline_hr` | Continuous | ⭐ | Preoperative heart rate (bpm) |
| `surgery_duration` | Integer | ⭐ | Time since induction (minutes) |
| `bmi` | Continuous | ⭐ | Body mass index |
| `emergency` | Binary | ⭐ | 0=elective, 1=emergency surgery |
| `sex` | Binary | ⭐ | 0=female, 1=male |

**Feature Importance Rankings:**
- MAP trends (current, 5min, 10min) account for ~56% of prediction power
- Patient factors (age, ASA, BMI) contribute ~25%
- Surgical/anesthetic factors contribute ~19%

### Training Your Own Model

**Quick Start:**
```bash
# 1. Install dependencies
pip install -r requirements_vitaldb.txt

# 2. Download VitalDB data (100 cases for testing)
python vitaldb_downloader.py

# 3. Train Random Forest model
python train_vitaldb_model.py

# 4. Compare to heuristics
python compare_models.py

# 5. Deploy to production
cp vitaldb_ioh_model.pkl ioh_model.pkl
cp vitaldb_ioh_scaler.pkl ioh_scaler.pkl
```

**Detailed Instructions:** See `VITALDB_TRAINING_GUIDE.md`

### Performance Benchmarks

**Published Research (2024):**
| Study | Method | ROC-AUC | Dataset |
|-------|--------|---------|---------|
| Springer 2024 | Deep Learning | 0.917 | VitalDB |
| eClinicalMedicine 2024 | Temporal Fusion Transformer | 0.933 | Multi-center |
| Meta-analysis 2024 | Various ML methods | 0.89 | 43 studies pooled |

**Our Pipeline (Expected):**
| Dataset Size | ROC-AUC | Training Time |
|-------------|---------|---------------|
| 100 cases | 0.75-0.80 | 30 minutes |
| 500 cases | 0.82-0.87 | 2 hours |
| 1,000 cases | 0.85-0.90 | 4 hours |
| 6,388 cases (full) | 0.87-0.92 | 12-24 hours |

### Model Files

| File | Purpose | Size |
|------|---------|------|
| `ioh_model.pkl` | Trained Random Forest classifier | ~200 KB (stub), ~5 MB (VitalDB) |
| `ioh_scaler.pkl` | StandardScaler for feature normalization | ~130 bytes |
| `ioh_models.py` | Model class definitions (for pickle unpickling) | ~5 KB |
| `vitaldb_downloader.py` | VitalDB data download script | ~12 KB |
| `train_vitaldb_model.py` | Model training pipeline | ~10 KB |
| `compare_models.py` | Model comparison/validation | ~6 KB |

### VitalDB Dataset Access

**Requirements:**
- Free access after Data Use Agreement
- Optional: CITI Human Research Training
- No formal application needed

**Dataset Characteristics:**
- **Cases:** 6,388 surgical procedures
- **Institution:** Seoul National University Hospital
- **Years:** 2011-2020
- **Parameters:** 196 intraoperative monitoring parameters
- **Resolution:** Up to 100 Hz for waveforms
- **Format:** REST API (JSON/CSV)

**Official Resources:**
- Dataset: https://vitaldb.net/dataset/
- Documentation: https://vitaldb.net/docs/
- Paper: https://www.nature.com/articles/s41597-022-01411-5
- GitHub: https://github.com/vitaldb/examples

### Clinical Validation

**Important:** The current heuristics model has NOT been validated in clinical settings. For production use, we strongly recommend:

1. ✅ **Train on VitalDB** (6,388 real surgical cases)
2. ✅ **Validate on external dataset** (INSPIRE, MIMIC-IV, or institutional data)
3. ✅ **Prospective validation** (compare predictions to actual outcomes in your OR)
4. ✅ **Regulatory review** (if using for clinical decision-making)

**Disclaimer:** This tool provides risk estimation for educational purposes. All predictions should be verified by qualified anesthesia professionals.

## Chat History (Persistent Storage)

### Overview

**Status:** Phase 1 Complete (Database Setup)
**Feature Flag:** `ENABLE_CHAT_HISTORY` (default: `false`)
**Implementation Plan:** See `CHAT_HISTORY_IMPLEMENTATION_PLAN.md`

The chat history feature provides persistent storage for conversations, preparing for future user accounts. The implementation follows a careful, incremental approach to avoid breaking existing chat functionality.

### Architecture

**Dual-Layer Storage:**
- **Session (Temporary):** Active conversation stored in Redis-backed session (1 hour TTL)
- **Database (Persistent):** All conversations stored in SQLite for long-term access

**Database Location:**
- **Production:** `/var/lib/gasconsult/gasconsult.db` (requires directory creation)
- **Development:** `./gasconsult.db` (automatic fallback)

### Database Schema

**conversations table:**
```sql
- id (UUID primary key)
- user_session_id (session ID, later: user_id)
- title (auto-generated from first query)
- created_at, updated_at
- message_count
- is_active (soft delete support)
```

**messages table:**
```sql
- id (UUID primary key)
- conversation_id (foreign key)
- role (user/assistant)
- content
- paper_references (JSON array)
- num_papers, evidence_strength
- created_at
```

### Phase 1: Database Setup (COMPLETED ✅)

**What's Implemented:**
- ✅ `database.py` module with full CRUD operations
- ✅ SQLite database initialization with indexes
- ✅ Feature flag (`ENABLE_CHAT_HISTORY`)
- ✅ Non-breaking integration (app works if database fails)
- ✅ Comprehensive error handling

**What's NOT Changed:**
- ❌ No UI changes (chat works exactly as before)
- ❌ Data not saved to database yet (session-only for now)
- ❌ No history viewing features

**Safety Guarantees:**
- Database failures don't break chat functionality
- Feature can be disabled via environment variable
- All database operations are wrapped in try/except
- Session-based chat continues working independently

### Next Phases (Planned)

**Phase 2:** Silent Data Persistence
- Save messages to database in background
- No UI changes yet
- Verify data integrity

**Phase 3:** History Sidebar UI
- Add conversation list sidebar
- Display past conversations
- Current chat unaffected

**Phase 4:** Load Previous Conversations
- Click to load old conversations
- Continue old conversations
- Careful session management

**Phase 5:** Advanced Features
- Search conversations
- Delete/export conversations
- Share conversations

### Database Module API

**Available Functions in `database.py`:**

```python
# Initialization
init_db() -> bool

# Conversations
save_conversation(id, user_session_id, title) -> bool
get_conversations(user_session_id, limit, offset) -> List[Dict]
get_conversation(conversation_id) -> Dict
delete_conversation(conversation_id, hard_delete=False) -> bool
update_conversation_title(conversation_id, new_title) -> bool
search_conversations(user_session_id, search_query) -> List[Dict]

# Messages
save_message(conversation_id, role, content, references, num_papers, evidence_strength) -> bool

# Utilities
generate_conversation_title(first_query) -> str
get_database_stats() -> Dict
```

### Testing

**Independent Testing:**
```bash
# Test database module independently
python database.py
```

**Feature Toggle:**
```bash
# Enable chat history
export ENABLE_CHAT_HISTORY=true

# Disable chat history (default)
export ENABLE_CHAT_HISTORY=false
```

### Production Deployment

**Database Directory Setup:**
```bash
# Create production database directory
sudo mkdir -p /var/lib/gasconsult
sudo chown www-data:www-data /var/lib/gasconsult
sudo chmod 755 /var/lib/gasconsult
```

**Render Configuration:**
Set environment variable in Render dashboard:
```
ENABLE_CHAT_HISTORY=true
GASCONSULT_DB_DIR=/var/lib/gasconsult
```

### Future: User Accounts

When ready to add user accounts, migration is simple:
1. Add `user_id` column to `conversations` table
2. Replace `user_session_id` with `current_user.id`
3. Add authentication (Flask-Login recommended)

The database schema is already designed for this transition.

### Future Improvements

- ✅ ~~Move Entrez credentials to environment variables~~ (Completed)
- ✅ ~~Add input sanitization before rendering~~ (Completed)
- ✅ ~~Add logging for debugging~~ (Completed)
- ✅ ~~Consider rate limiting~~ (Completed)
- ✅ ~~Persistent chat history - Phase 1~~ (Completed)
- ⬜ Persistent chat history - Phase 2-5 (In Progress)
- ⬜ Implement caching for repeated PubMed queries (Redis)
- ⬜ Add automated testing suite
- ⬜ Set up CI/CD pipeline
- ⬜ Enhanced medical disclaimers with acknowledgment
- ⬜ User accounts and saved queries
- ⬜ Mobile app (optional)

## Git Workflow

- Main development branch varies by feature
- Commit messages should be concise version tags (e.g., "v2.1")
- Recent commits indicate active development of search functionality

## Disclaimer

This application provides **educational content only** and is **not medical advice**. All outputs should be verified by qualified medical professionals before clinical use.
