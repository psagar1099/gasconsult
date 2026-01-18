# GasConsult.ai - Comprehensive Refactor Plan

**Date:** 2026-01-18
**Branch:** `claude/fix-chat-issues-N810c`
**Goal:** Transform monolithic single-file architecture into maintainable, modular codebase

---

## Current State Analysis

### Problems Identified
- **app.py:** 35,399 lines (unmanageable monolith)
- **18 HTML template constants** as string variables in app.py
- **43 duplicate `.nav` CSS definitions** across templates
- **12 duplicate `.chat-hints` CSS definitions** across templates
- **9 duplicate `.hint-chip` CSS definitions** across templates
- **No template inheritance** - every page duplicates header/footer/nav/CSS
- **Inline CSS/JS** scattered across 18 templates (thousands of duplicate lines)
- **No code splitting** - entire app in one file makes changes risky

### Impact on Development
- ❌ **Any navbar change requires editing 14 locations**
- ❌ **CSS changes require editing multiple templates**
- ❌ **Chat function breakage affects entire app.py**
- ❌ **Difficult to add new pages** (copy/paste 1000+ lines of boilerplate)
- ❌ **Testing is hard** - can't isolate routes or logic
- ❌ **Git conflicts** - everyone edits same massive file

---

## New Architecture

### Directory Structure
```
gasconsult/
├── app.py                          # Lightweight Flask app initialization (200 lines)
├── requirements.txt
├── .env.example
├── database.py                     # Existing database module ✓
├── config.py                       # NEW: Configuration management
├── wsgi.py                         # NEW: Production WSGI entry point
│
├── templates/                      # NEW: Jinja2 templates
│   ├── base.html                   # Base template with navbar/footer/CSS
│   ├── components/                 # Reusable UI components
│   │   ├── navbar.html
│   │   ├── footer.html
│   │   ├── evidence_badge.html
│   │   └── citation_preview.html
│   ├── chat.html                   # Homepage/chat interface
│   ├── preop.html                  # Pre-operative assessment
│   ├── hypotension.html            # IOH predictor
│   ├── dose_calc.html              # Dose calculator
│   ├── calculators.html            # Clinical calculators hub
│   ├── crisis.html                 # Crisis protocols
│   ├── difficult_airway.html       # Difficult airway
│   ├── informed_consent.html       # Informed consent
│   ├── evidence.html               # Evidence methodology
│   ├── library.html                # User library
│   ├── shared_response.html        # Shared response view
│   ├── auth/                       # Authentication pages
│   │   ├── login.html
│   │   ├── register.html
│   │   └── resend_verification.html
│   ├── admin/                      # Admin pages
│   │   └── dashboard.html
│   ├── legal/                      # Legal pages
│   │   ├── terms.html
│   │   ├── privacy.html
│   │   └── pricing.html
│   └── errors/                     # Error pages
│       ├── 404.html
│       ├── 500.html
│       └── 503.html
│
├── static/                         # Static assets
│   ├── css/                        # NEW: Organized CSS
│   │   ├── main.css                # Core styles ✓
│   │   ├── chat.css                # Chat-specific styles
│   │   ├── forms.css               # Form styles
│   │   └── components.css          # Component styles
│   ├── js/                         # NEW: Organized JavaScript
│   │   ├── utils.js                # Shared utilities ✓
│   │   ├── chat.js                 # Chat/streaming logic
│   │   ├── forms.js                # Form handling
│   │   └── calculators.js          # Calculator logic
│   ├── favicon.svg                 # Existing ✓
│   ├── logo.png                    # Existing ✓
│   ├── manifest.json               # Existing ✓
│   └── sw.js                       # Existing ✓
│
├── routes/                         # NEW: Route handlers
│   ├── __init__.py                 # Blueprint registration
│   ├── chat.py                     # Chat/streaming routes
│   ├── preop.py                    # Pre-op assessment routes
│   ├── hypotension.py              # IOH predictor routes
│   ├── calculators.py              # Calculator routes
│   ├── crisis.py                   # Crisis protocol routes
│   ├── auth.py                     # Authentication routes
│   ├── admin.py                    # Admin routes
│   └── api.py                      # API endpoints
│
├── services/                       # NEW: Business logic
│   ├── __init__.py
│   ├── pubmed.py                   # PubMed search logic
│   ├── gpt.py                      # GPT-4o integration
│   ├── evidence.py                 # Evidence strength calculation
│   ├── calculators.py              # Clinical calculator logic
│   ├── preop_agent.py              # Pre-op agentic AI
│   └── ioh_predictor.py            # IOH ML model
│
├── utils/                          # NEW: Utilities
│   ├── __init__.py
│   ├── sanitization.py             # Input sanitization
│   ├── query_processing.py         # Query expansion, abbreviations
│   ├── session.py                  # Session management helpers
│   └── validators.py               # Input validation
│
└── models/                         # NEW: Data models (future)
    ├── __init__.py
    ├── user.py                     # User model
    └── conversation.py             # Conversation model
```

---

## Implementation Phases

### Phase 1: Template Extraction (HIGHEST PRIORITY)
**Goal:** Eliminate 90% of code duplication

**Steps:**
1. Create `templates/` directory
2. Extract `base.html` with:
   - `<!DOCTYPE html>` through `<head>` setup
   - Navbar (single definition)
   - Footer (single definition)
   - All shared CSS (consolidated from all templates)
   - All shared JavaScript
   - Jinja2 blocks: `{% block title %}`, `{% block content %}`, `{% block scripts %}`
3. Create `templates/components/navbar.html` (with auth button logic)
4. Create `templates/components/footer.html`
5. Convert all 18 template constants to separate `.html` files
6. Update routes to use `render_template()` instead of `render_template_string()`

**Impact:**
- **app.py shrinks from 35,399 lines to ~8,000 lines** (77% reduction)
- **Navbar changes require editing 1 file** (not 14)
- **CSS changes in 1 location** (not 43)
- **New pages need ~50 lines** (not 1000+)

---

### Phase 2: CSS Consolidation
**Goal:** Single source of truth for styles

**Steps:**
1. Move all inline `<style>` blocks to external CSS files
2. Create `static/css/chat.css` (chat-specific styles)
3. Create `static/css/forms.css` (form styles)
4. Create `static/css/components.css` (badges, cards, etc.)
5. Update `base.html` to load CSS in correct order
6. Remove all inline CSS from templates

**Impact:**
- **Browser caching** - CSS loaded once, cached across pages
- **No duplicate CSS** - 100% consistency
- **Easier debugging** - all styles in organized files
- **Smaller HTML payloads** - faster page loads

---

### Phase 3: JavaScript Organization
**Goal:** Modular, maintainable JavaScript

**Steps:**
1. Extract chat/streaming JavaScript to `static/js/chat.js`
2. Extract form handling to `static/js/forms.js`
3. Extract calculator logic to `static/js/calculators.js`
4. Keep `static/js/utils.js` for shared utilities
5. Update `base.html` to load JS in correct order
6. Remove all inline JavaScript from templates

**Impact:**
- **No duplicate JavaScript** across templates
- **Easier testing** - isolated modules
- **Better browser caching**

---

### Phase 4: Route Splitting
**Goal:** app.py becomes lightweight orchestrator

**Steps:**
1. Create `routes/` directory with Flask Blueprints
2. Extract chat routes to `routes/chat.py`:
   - `@app.route("/", methods=["GET", "POST"])` → `@bp.route("/")`
   - `/stream` endpoint
   - `/clear`, `/clear-chat` endpoints
3. Extract preop routes to `routes/preop.py`
4. Extract hypotension routes to `routes/hypotension.py`
5. Extract calculator routes to `routes/calculators.py`
6. Extract auth routes to `routes/auth.py`
7. Update `app.py` to register blueprints

**Impact:**
- **app.py shrinks to ~200 lines** (initialization only)
- **Routes isolated** - chat changes don't affect preop
- **Easier testing** - test individual blueprints
- **Parallel development** - team members work on different route files

---

### Phase 5: Business Logic Extraction
**Goal:** Separate concerns - routes vs logic

**Steps:**
1. Create `services/pubmed.py`:
   - `search_pubmed()`, `fetch_papers()`, `extract_metadata()`
2. Create `services/gpt.py`:
   - `generate_response()`, `stream_response()`
3. Create `services/evidence.py`:
   - `get_evidence_strength()`, `calculate_confidence_percentage()`
4. Create `services/calculators.py`:
   - All calculator functions (IBW, MABL, BSA, QTc, etc.)
5. Create `services/preop_agent.py`:
   - `PreopAgent` class and helpers
6. Create `services/ioh_predictor.py`:
   - IOH ML model loading and prediction

**Impact:**
- **Reusable logic** - services can be called from multiple routes
- **Easier testing** - unit test individual services
- **Better organization** - clear separation of concerns

---

### Phase 6: Utility Organization
**Goal:** Shared utilities in logical modules

**Steps:**
1. Create `utils/query_processing.py`:
   - `expand_medical_abbreviations()`, `resolve_references()`, `detect_question_type()`
2. Create `utils/sanitization.py`:
   - `sanitize_user_query()`, `sanitize_input()`
3. Create `utils/session.py`:
   - `get_or_create_conversation_id()`, session helpers
4. Create `utils/validators.py`:
   - Input validation functions

**Impact:**
- **No more 500+ line helper functions in app.py**
- **Easier to find utilities** - organized by purpose
- **Reusable across services**

---

### Phase 7: Configuration Management
**Goal:** Environment-aware configuration

**Steps:**
1. Create `config.py`:
   - `Config`, `DevelopmentConfig`, `ProductionConfig` classes
   - Load from environment variables
2. Update `app.py` to use `app.config.from_object(config)`
3. Remove hardcoded configuration scattered in app.py

**Impact:**
- **Environment-specific settings** (dev vs prod)
- **Easier deployment** - no code changes needed
- **Security** - sensitive config not in code

---

## Migration Strategy

### Incremental Migration (SAFE)
We'll migrate **one page at a time** to avoid breaking the entire app:

1. ✅ Create new directory structure
2. ✅ Extract `base.html` with navbar/footer
3. ✅ Convert `/` (chat) route first → test thoroughly
4. ✅ Convert `/preop` route → test
5. ✅ Convert remaining routes one-by-one
6. ✅ Remove old template constants after all routes migrated
7. ✅ Extract services once routes stabilized

### Testing at Each Step
- Manual testing of converted routes
- Verify no regressions on unconverted routes
- Git commit after each successful page conversion
- Rollback individual commits if issues arise

---

## File Size Projections

### Before Refactor
| File | Lines | Size |
|------|-------|------|
| `app.py` | 35,399 | 1.2 MB |
| **Total** | **35,399** | **1.2 MB** |

### After Refactor
| File | Lines | Size | Notes |
|------|-------|------|-------|
| `app.py` | 200 | 8 KB | Initialization only |
| `routes/*.py` (8 files) | 3,000 | 120 KB | Route handlers |
| `services/*.py` (6 files) | 2,500 | 100 KB | Business logic |
| `templates/*.html` (18 files) | 4,000 | 150 KB | No duplicate CSS/JS |
| `static/css/*.css` (4 files) | 1,500 | 50 KB | Consolidated styles |
| `static/js/*.js` (4 files) | 1,000 | 40 KB | Organized JS |
| `utils/*.py` (4 files) | 800 | 32 KB | Utilities |
| **Total** | **13,000** | **500 KB** | 63% reduction |

**Key Wins:**
- **63% smaller codebase** (35K → 13K lines)
- **58% smaller file size** (1.2 MB → 500 KB)
- **Eliminates 22K lines of duplicate CSS/HTML/JS**

---

## Benefits Summary

### Developer Experience
✅ **Easier changes** - edit 1 file instead of 14
✅ **Faster development** - less boilerplate for new pages
✅ **Fewer bugs** - no more "forgot to update 1 of 14 navbars"
✅ **Better IDE support** - separate files = better autocomplete
✅ **Parallel work** - team members don't conflict on app.py
✅ **Easier testing** - isolate and test individual modules

### Performance
✅ **Faster page loads** - browser caches CSS/JS across pages
✅ **Smaller HTML payloads** - no inline CSS/JS in every page
✅ **Better minification** - external CSS/JS can be minified/compressed

### Maintenance
✅ **Easier debugging** - find bugs in isolated modules
✅ **Code reviews** - reviewers see only changed module, not entire app.py
✅ **Onboarding** - new devs understand structure faster
✅ **Refactoring** - safe to change isolated modules

---

## Risk Mitigation

### Potential Risks
1. **Breaking existing functionality** during migration
2. **Template rendering issues** (missing variables, incorrect paths)
3. **Import circular dependencies** when splitting modules
4. **Session handling changes** breaking chat
5. **CSS specificity conflicts** when consolidating styles

### Mitigation Strategies
1. ✅ **Incremental migration** - one page at a time
2. ✅ **Git commits after each page** - easy rollback
3. ✅ **Keep old code until verified** - dual implementation during migration
4. ✅ **Thorough testing** - manual QA on each converted page
5. ✅ **Preserve all existing functionality** - no feature changes during refactor
6. ✅ **Use same Flask patterns** - render_template instead of render_template_string

---

## Timeline Estimate

| Phase | Effort | Duration |
|-------|--------|----------|
| Phase 1: Template Extraction | High | 3-4 hours |
| Phase 2: CSS Consolidation | Medium | 2 hours |
| Phase 3: JavaScript Organization | Medium | 2 hours |
| Phase 4: Route Splitting | High | 3 hours |
| Phase 5: Business Logic Extraction | Medium | 2 hours |
| Phase 6: Utility Organization | Low | 1 hour |
| Phase 7: Configuration Management | Low | 1 hour |
| **Total** | - | **14-15 hours** |

**Note:** This is a one-time investment that will save hundreds of hours over the project's lifetime.

---

## Success Metrics

### Before Refactor
- ❌ Average time to add new page: **60 minutes** (copy/paste 1000+ lines, fix duplicates)
- ❌ Average time to update navbar: **30 minutes** (edit 14 locations)
- ❌ Lines of duplicate code: **22,000+**
- ❌ Git conflicts per week: **High** (everyone edits app.py)

### After Refactor (Target)
- ✅ Average time to add new page: **10 minutes** (extend base.html, write content)
- ✅ Average time to update navbar: **2 minutes** (edit 1 file)
- ✅ Lines of duplicate code: **0**
- ✅ Git conflicts per week: **Minimal** (separate files per feature)

---

## Next Steps

1. ✅ Get approval for refactor plan
2. ✅ Create new directory structure
3. ✅ Start Phase 1: Extract base template
4. ✅ Convert chat page (/)
5. ✅ Test thoroughly
6. ✅ Continue with remaining pages
7. ✅ Extract services and utilities
8. ✅ Final testing and deployment

---

## Questions & Decisions

**Q: Should we maintain backward compatibility with old template constants?**
A: Yes, during migration. Remove after all routes converted.

**Q: Should we migrate database.py to a models/ directory?**
A: Not in this refactor. Keep existing database.py working. Future enhancement.

**Q: Should we use Flask Blueprints or keep all routes in app.py?**
A: Use Blueprints for better organization and testability.

**Q: Should we keep existing functionality or improve during refactor?**
A: Keep existing functionality. Refactor first, improve later (separate PRs).

---

**STATUS:** Ready for implementation
**APPROVED BY:** Awaiting user confirmation
**IMPLEMENTATION START:** After approval
