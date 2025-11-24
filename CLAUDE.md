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
| Language | Python 3.x | - |

## Directory Structure

```
gasconsult/
├── app.py              # Main Flask application (single-file architecture)
├── requirements.txt    # Python dependencies
├── README.md           # Project description
├── CLAUDE.md           # This file - AI assistant guide
└── static/
    ├── favicon.ico     # Browser favicon
    └── logo.png        # Site logo
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

### Synonym Expansion (lines 80-85)

Automatically expands common clinical abbreviations:
- `txa` → `"tranexamic acid" OR TXA`
- `blood loss` → `"blood loss" OR hemorrhage OR transfusion`
- `spine surgery` → `"spine surgery" OR "spinal fusion" OR scoliosis`
- `pediatric` → `pediatric OR children OR peds`
- `ponv` → `PONV OR "postoperative nausea"`

## Environment Variables

| Variable | Required | Description |
|----------|----------|-------------|
| `OPENAI_API_KEY` | Yes | OpenAI API key for GPT-4o access |

**Note:** Entrez credentials are currently hardcoded in `app.py` (lines 8-9). Consider moving to environment variables for production.

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

### Security Notes

- API keys should be moved to environment variables
- Entrez email/API key are currently hardcoded (lines 8-9)
- OpenAI key properly loaded from environment
- User input is processed but not sanitized for display (potential XSS via `|safe` filter)

### Testing Queries

Good test queries for the application:
- "TXA in spine surgery"
- "blood loss scoliosis"
- "propofol peds"
- "PONV prevention"

### Potential Improvements

- Move Entrez credentials to environment variables
- Add input sanitization before rendering
- Implement caching for repeated queries
- Add logging for debugging
- Consider rate limiting

## Git Workflow

- Main development branch varies by feature
- Commit messages should be concise version tags (e.g., "v2.1")
- Recent commits indicate active development of search functionality

## Disclaimer

This application provides **educational content only** and is **not medical advice**. All outputs should be verified by qualified medical professionals before clinical use.
