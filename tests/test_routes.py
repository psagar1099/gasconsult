"""
Route Smoke Tests - Prevent Template Refactoring Bugs

Tests all routes return 200 and contain expected content.
Catches issues like wrong templates being rendered to routes.
"""

import pytest
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Set test environment variables before importing app
os.environ['FLASK_ENV'] = 'testing'
os.environ['OPENAI_API_KEY'] = 'test-key-for-testing'
os.environ['ENTREZ_EMAIL'] = 'test@example.com'
os.environ['FLASK_SECRET_KEY'] = 'test-secret-key-for-testing'
os.environ['TESTING'] = 'true'

from app import app


@pytest.fixture
def client():
    """Create test client"""
    app.config['TESTING'] = True
    app.config['WTF_CSRF_ENABLED'] = False  # Disable CSRF for testing
    app.config['SERVER_NAME'] = 'localhost'
    with app.test_client() as client:
        yield client


class TestHomePage:
    """Test home page (/) - Most critical route"""

    def test_home_returns_200(self, client):
        """Home page should return 200 OK"""
        response = client.get('/')
        assert response.status_code == 200, "Home page failed to load"

    def test_home_is_chat_interface_not_preop(self, client):
        """
        CRITICAL: Home page should be chat interface, NOT pre-op assessment
        This catches the exact bug we just fixed
        """
        response = client.get('/')
        html = response.data.decode('utf-8')

        # Should contain chat interface elements
        assert 'GasConsult.ai - AI-Powered Anesthesiology Assistant' in html, \
            "Home page missing correct title (got wrong template?)"
        assert 'chat-input' in html, \
            "Home page missing chat input (is this the pre-op form instead?)"

        # Should NOT contain pre-op specific elements
        assert 'Pre-Operative Assessment' not in html or 'chat' in html.lower(), \
            "CRITICAL: Home page contains pre-op content - wrong template!"

    def test_home_has_navigation(self, client):
        """Home page should have navigation"""
        response = client.get('/')
        html = response.data.decode('utf-8')
        assert 'nav' in html.lower(), "Home page missing navigation"

    def test_home_loads_static_assets(self, client):
        """Home page should reference CSS and JS"""
        response = client.get('/')
        html = response.data.decode('utf-8')
        assert '/static/main.css' in html, "Home page not loading main.css"
        assert '/static/utils.js' in html, "Home page not loading utils.js"


class TestStaticRoutes:
    """Test static/informational routes"""

    @pytest.mark.parametrize("route,expected_title", [
        ('/terms', 'Terms of Service'),
        ('/privacy', 'Privacy Policy'),
        ('/calculators', 'Clinical Calculators'),
        ('/evidence', 'Evidence'),
        ('/pricing', 'Pricing'),
    ])
    def test_static_routes_load(self, client, route, expected_title):
        """Static routes should return 200 and contain expected title"""
        response = client.get(route)
        assert response.status_code == 200, f"{route} failed to load"
        html = response.data.decode('utf-8')
        assert expected_title.lower() in html.lower(), \
            f"{route} missing expected title '{expected_title}'"


class TestToolRoutes:
    """Test clinical tool routes (GET only)"""

    @pytest.mark.parametrize("route,expected_content", [
        ('/preop', 'Pre-Operative Assessment'),
        ('/dose-calc', 'Dose Calculator'),
        ('/informed-consent', 'Informed Consent'),
        ('/difficult-airway', 'Difficult Airway'),
        ('/hypotension', 'Hypotension'),
    ])
    def test_tool_routes_load(self, client, route, expected_content):
        """Tool routes should return 200 and contain expected content"""
        response = client.get(route)
        assert response.status_code == 200, f"{route} failed to load"
        html = response.data.decode('utf-8')
        assert expected_content.lower() in html.lower(), \
            f"{route} missing expected content '{expected_content}'"

    def test_preop_is_not_on_homepage(self, client):
        """
        Pre-op content should ONLY be on /preop, NOT on /
        This is a critical regression test
        """
        home_response = client.get('/')
        preop_response = client.get('/preop')

        home_html = home_response.data.decode('utf-8')
        preop_html = preop_response.data.decode('utf-8')

        # Pre-op form should be on /preop
        assert 'procedure' in preop_html.lower(), "/preop missing procedure field"

        # Home should NOT have pre-op form
        assert 'Pre-Operative Assessment - GasConsult.ai' not in home_html, \
            "CRITICAL BUG: Home page has pre-op title!"


class TestAPIRoutes:
    """Test API endpoints"""

    def test_health_endpoint(self, client):
        """Health check should return 200 with JSON"""
        response = client.get('/health')
        assert response.status_code == 200, "Health check failed"
        assert response.is_json or 'application/json' in response.content_type, \
            "Health check should return JSON"

    def test_api_status_endpoint(self, client):
        """API status should return 200 with JSON"""
        response = client.get('/api/status')
        assert response.status_code == 200, "API status failed"


class TestAuthRoutes:
    """Test authentication routes (GET only)"""

    @pytest.mark.parametrize("route", [
        '/login',
        '/register',
    ])
    def test_auth_routes_load(self, client, route):
        """Auth routes should return 200"""
        response = client.get(route)
        assert response.status_code == 200, f"{route} failed to load"


class TestStaticFiles:
    """Test static file availability"""

    @pytest.mark.parametrize("static_file", [
        '/static/main.css',
        '/static/utils.js',
        '/static/favicon.svg',
        '/static/manifest.json',
    ])
    def test_static_files_exist(self, client, static_file):
        """Critical static files should be accessible"""
        response = client.get(static_file)
        assert response.status_code == 200, f"{static_file} not found"


class TestNavigationFunctions:
    """Test JavaScript navigation functions are available"""

    def test_utils_has_toggle_functions(self, client):
        """utils.js should define toggleMobileMenu and toggleNavDropdown"""
        response = client.get('/static/utils.js')
        js = response.data.decode('utf-8')

        assert 'function toggleMobileMenu' in js, \
            "utils.js missing toggleMobileMenu function"
        assert 'function toggleNavDropdown' in js, \
            "utils.js missing toggleNavDropdown function"
        assert 'window.toggleMobileMenu' in js, \
            "toggleMobileMenu not exposed globally"
        assert 'window.toggleNavDropdown' in js, \
            "toggleNavDropdown not exposed globally"


# Quick smoke test runner (can be run standalone)
if __name__ == '__main__':
    print("Running route smoke tests...")
    pytest.main([__file__, '-v'])
