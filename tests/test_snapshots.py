"""
HTML Snapshot Tests - Detect Template Changes

Stores known-good HTML snippets and compares against current templates.
Alerts when templates change unexpectedly.
"""

import pytest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Set test environment
os.environ['FLASK_ENV'] = 'testing'
os.environ['OPENAI_API_KEY'] = 'test-key'
os.environ['ENTREZ_EMAIL'] = 'test@example.com'
os.environ['FLASK_SECRET_KEY'] = 'test-secret'
os.environ['TESTING'] = 'true'

from app import app


@pytest.fixture
def client():
    """Create test client"""
    app.config['TESTING'] = True
    app.config['WTF_CSRF_ENABLED'] = False
    with app.test_client() as client:
        yield client


# Known-good HTML snippets (extracted from working templates)
# These act as "golden masters" to detect unexpected changes

HOME_PAGE_SNAPSHOTS = {
    'title': 'GasConsult.ai - AI-Powered Anesthesiology Assistant',
    'meta_description': 'Evidence-based anesthesiology AI consultant',
    'hero_section_hint': 'chat-input',  # Should have chat input
    'static_css': '/static/main.css',
    'static_js': '/static/utils.js',
}

PREOP_PAGE_SNAPSHOTS = {
    'title': 'Pre-Operative Assessment',
    'form_element': 'procedure',  # Should have procedure field
    'submit_button': 'submit',
}

DOSE_CALC_SNAPSHOTS = {
    'title': 'Dose Calculator',
    'weight_input': 'weight',
}


class TestHomePageSnapshot:
    """Snapshot tests for home page"""

    def test_home_title_unchanged(self, client):
        """Home page title should match known-good snapshot"""
        response = client.get('/')
        html = response.data.decode('utf-8')
        assert HOME_PAGE_SNAPSHOTS['title'] in html, \
            "Home page title changed unexpectedly - template may have been swapped"

    def test_home_meta_description_unchanged(self, client):
        """Meta description should match known-good snapshot"""
        response = client.get('/')
        html = response.data.decode('utf-8')
        assert HOME_PAGE_SNAPSHOTS['meta_description'] in html, \
            "Home page meta description changed - wrong template?"

    def test_home_has_chat_input(self, client):
        """Home page must have chat input field"""
        response = client.get('/')
        html = response.data.decode('utf-8')
        assert HOME_PAGE_SNAPSHOTS['hero_section_hint'] in html, \
            "Home page missing chat-input - this is NOT the home page!"

    def test_home_loads_expected_css(self, client):
        """CSS reference should match snapshot"""
        response = client.get('/')
        html = response.data.decode('utf-8')
        assert HOME_PAGE_SNAPSHOTS['static_css'] in html, \
            "Home page CSS reference changed"

    def test_home_loads_expected_js(self, client):
        """JS reference should match snapshot"""
        response = client.get('/')
        html = response.data.decode('utf-8')
        assert HOME_PAGE_SNAPSHOTS['static_js'] in html, \
            "Home page JS reference changed"


class TestPreopPageSnapshot:
    """Snapshot tests for pre-op page"""

    def test_preop_title_unchanged(self, client):
        """Pre-op page title should match snapshot"""
        response = client.get('/preop')
        html = response.data.decode('utf-8')
        assert PREOP_PAGE_SNAPSHOTS['title'] in html, \
            "Pre-op page title changed unexpectedly"

    def test_preop_has_procedure_field(self, client):
        """Pre-op page must have procedure input"""
        response = client.get('/preop')
        html = response.data.decode('utf-8')
        assert PREOP_PAGE_SNAPSHOTS['form_element'] in html, \
            "Pre-op page missing procedure field - wrong template?"


class TestDoseCalcPageSnapshot:
    """Snapshot tests for dose calculator"""

    def test_dose_calc_title_unchanged(self, client):
        """Dose calc title should match snapshot"""
        response = client.get('/dose-calc')
        html = response.data.decode('utf-8')
        assert DOSE_CALC_SNAPSHOTS['title'] in html, \
            "Dose calc page title changed"

    def test_dose_calc_has_weight_input(self, client):
        """Dose calc must have weight input"""
        response = client.get('/dose-calc')
        html = response.data.decode('utf-8')
        assert DOSE_CALC_SNAPSHOTS['weight_input'] in html, \
            "Dose calc missing weight input"


class TestCrossPageValidation:
    """Ensure pages don't have each other's content"""

    def test_home_is_not_preop(self, client):
        """
        CRITICAL REGRESSION TEST:
        Home page should NOT contain pre-op specific content
        """
        home = client.get('/').data.decode('utf-8')
        preop = client.get('/preop').data.decode('utf-8')

        # Check that home is different from preop
        home_title = HOME_PAGE_SNAPSHOTS['title']
        preop_title = PREOP_PAGE_SNAPSHOTS['title']

        assert home_title in home, "Home missing its own title"
        assert preop_title in preop, "Preop missing its own title"

        # Home should NOT have preop title
        assert 'Pre-Operative Assessment - GasConsult.ai' not in home, \
            "CRITICAL: Home page contains pre-op title - WRONG TEMPLATE!"

    def test_preop_is_not_home(self, client):
        """Pre-op page should NOT be home page"""
        preop = client.get('/preop').data.decode('utf-8')

        # Should have pre-op content
        assert 'procedure' in preop.lower(), "Preop missing procedure field"

        # Should NOT be generic home page
        assert 'Pre-Operative Assessment' in preop, \
            "Preop doesn't have pre-op title - wrong template?"


# Quick test runner
if __name__ == '__main__':
    print("Running HTML snapshot tests...")
    pytest.main([__file__, '-v'])
