# Testing Infrastructure - Prevent Template Refactoring Bugs

## Overview

This testing suite prevents bugs like the one where `chat_standalone.html` contained the wrong template (pre-op instead of home page). It provides multiple layers of automated validation.

## Quick Start

### Install Test Dependencies

```bash
pip install -r requirements.txt
```

This installs:
- `pytest==8.0.0` - Test framework
- `pytest-flask==1.3.0` - Flask testing utilities

### Run Tests

```bash
# Run all tests
./run_tests.sh

# Run only critical regression tests (fastest)
./run_tests.sh critical

# Run quick smoke test (home page only)
./run_tests.sh quick

# Run specific test suites
./run_tests.sh templates   # Template validation
./run_tests.sh snapshots   # HTML snapshots
./run_tests.sh routes      # Route smoke tests
```

## Test Suites

### 1. Route Smoke Tests (`test_routes.py`)

**Purpose:** Verify all routes return 200 and contain expected content

**Key Tests:**
- ✅ Home page returns 200 OK
- ✅ Home page is chat interface (NOT pre-op)
- ✅ All tool routes load correctly
- ✅ Static files are accessible
- ✅ Navigation JavaScript functions exist

**Example:**
```python
def test_home_is_chat_interface_not_preop(self, client):
    """
    CRITICAL: Home page should be chat interface, NOT pre-op assessment
    This catches the exact bug we just fixed
    """
    response = client.get('/')
    html = response.data.decode('utf-8')

    # Should contain chat interface elements
    assert 'chat-input' in html, \
        "Home page missing chat input (is this the pre-op form instead?)"
```

**Run:** `./run_tests.sh routes`

---

### 2. Template Validation Tests (`test_templates.py`)

**Purpose:** Validate template files contain correct content (catches template swap bugs)

**Key Tests:**
- ✅ `chat_standalone.html` is home page template
- ✅ `chat_standalone.html` NOT pre-op template
- ✅ Template loads correct CSS/JS
- ✅ Base template has Jinja2 blocks
- ✅ Navigation components call correct functions

**Example:**
```python
def test_chat_standalone_is_home_page(self):
    """chat_standalone.html should contain HOME PAGE content, not pre-op"""
    template_path = TEMPLATES_DIR / 'chat_standalone.html'
    content = template_path.read_text()

    # MUST contain home page title
    assert 'GasConsult.ai - AI-Powered Anesthesiology Assistant' in content, \
        "CRITICAL: chat_standalone.html missing home page title!"
```

**Run:** `./run_tests.sh templates`

---

### 3. HTML Snapshot Tests (`test_snapshots.py`)

**Purpose:** Detect unexpected HTML changes (regression testing)

**Key Tests:**
- ✅ Home page title unchanged
- ✅ Pre-op page title unchanged
- ✅ Home page has chat input
- ✅ Pre-op page has procedure field
- ✅ Pages don't have each other's content

**How it works:**
Stores "known-good" HTML snippets and compares current output:

```python
HOME_PAGE_SNAPSHOTS = {
    'title': 'GasConsult.ai - AI-Powered Anesthesiology Assistant',
    'hero_section_hint': 'chat-input',
}
```

**Run:** `./run_tests.sh snapshots`

---

## Pre-Commit Hook

**Automatic Protection:** Tests run before every commit that touches templates or `app.py`

### How It Works

1. You try to commit changes to templates or app.py
2. Git hook runs CRITICAL tests automatically
3. If tests fail → commit is BLOCKED
4. If tests pass → commit proceeds

### Example Output

```bash
$ git commit -m "Update template"

Running pre-commit tests...

Template files changed:
templates/chat_standalone.html

Running CRITICAL regression tests...

test_home_is_chat_interface_not_preop PASSED
test_preop_is_not_on_homepage PASSED
test_chat_standalone_is_home_page PASSED
test_home_is_not_preop PASSED

✓ Critical tests PASSED - commit allowed
```

### Bypass (NOT RECOMMENDED)

```bash
git commit --no-verify -m "Skip tests (dangerous)"
```

---

## What Each Test Catches

### The Bug We Fixed

**Problem:** `chat_standalone.html` contained pre-op HTML instead of home page HTML

**Tests that catch this:**
1. `test_home_is_chat_interface_not_preop` - Checks home page has chat input
2. `test_chat_standalone_is_home_page` - Validates template file content
3. `test_home_is_not_preop` - Ensures home ≠ pre-op
4. `test_preop_is_not_on_homepage` - Cross-validation

**All 4 tests would FAIL** if this bug happens again.

---

## Test Files Structure

```
gasconsult/
├── tests/
│   ├── __init__.py
│   ├── test_routes.py        # Route smoke tests
│   ├── test_templates.py     # Template validation
│   └── test_snapshots.py     # HTML snapshots
├── run_tests.sh               # Test runner script
├── pytest.ini                 # Pytest configuration
├── .git/hooks/pre-commit      # Pre-commit hook
└── TESTING.md                 # This file
```

---

## Adding New Tests

### When to Add Tests

**After every template refactor:**
1. Add snapshot for new route
2. Add validation for new template file
3. Update cross-page validation

### Example: Add Test for New Route

```python
# In test_routes.py
def test_new_tool_route_loads(self, client):
    """New tool should return 200 with expected content"""
    response = client.get('/new-tool')
    assert response.status_code == 200
    assert 'New Tool Title' in response.data.decode('utf-8')
```

### Example: Add Template Validation

```python
# In test_templates.py
def test_new_template_has_correct_content(self):
    """new_template.html should contain expected content"""
    template_path = TEMPLATES_DIR / 'new_template.html'
    content = template_path.read_text()
    assert 'Expected Title' in content
```

---

## CI/CD Integration (Future)

### GitHub Actions (Recommended)

```yaml
# .github/workflows/test.yml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: actions/setup-python@v2
        with:
          python-version: '3.10'
      - run: pip install -r requirements.txt
      - run: ./run_tests.sh all
```

### Render Deploy Hook

Add to `render.yaml`:
```yaml
services:
  - type: web
    buildCommand: "pip install -r requirements.txt && ./run_tests.sh critical"
```

---

## Troubleshooting

### Tests Fail After Template Change

**Good!** That means tests are working. Check:
1. Did you put the right content in the right template file?
2. Did you update a route to use a different template?
3. Do you need to update snapshots for intentional changes?

### Pre-Commit Hook Blocks Commit

**Options:**
1. Fix the failing tests (recommended)
2. Run `./run_tests.sh critical` to see detailed output
3. Bypass with `git commit --no-verify` (dangerous)

### Pytest Not Found

```bash
pip install pytest pytest-flask
```

### Import Errors

Make sure you're in the project root:
```bash
cd /home/user/gasconsult
pytest tests/
```

---

## Performance

### Test Speed

| Test Suite | Tests | Time | Use Case |
|------------|-------|------|----------|
| Quick | 4 | ~2s | Pre-commit (home page only) |
| Critical | 4 | ~3s | Pre-commit hook |
| Routes | 15+ | ~5s | Smoke testing |
| Templates | 10+ | ~1s | Template validation |
| Snapshots | 10+ | ~3s | Regression testing |
| **All** | **40+** | **~10s** | **Full validation** |

---

## Best Practices

### ✅ DO

- Run `./run_tests.sh critical` before committing template changes
- Add tests when creating new routes/templates
- Update snapshots when making intentional changes
- Review test failures carefully (they catch bugs!)

### ❌ DON'T

- Skip tests with `--no-verify` unless absolutely necessary
- Ignore failing tests
- Copy HTML to wrong template files
- Change routes without updating tests

---

## Summary

**3-Layer Protection:**

1. **Pre-commit Hook** - Blocks bad commits automatically
2. **Test Suites** - Comprehensive validation (routes, templates, snapshots)
3. **Manual Testing** - Run before refactoring (`./run_tests.sh all`)

**Result:** Template refactoring bugs caught BEFORE they reach production 🎯

---

## Questions?

**Run tests:** `./run_tests.sh`

**Check test output:** `./run_tests.sh critical`

**See what hook does:** `cat .git/hooks/pre-commit`

**Need help?** Check test output or read test file comments
