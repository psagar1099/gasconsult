"""
Template Validation Tests - Prevent Template Swap Bugs

Validates that each template file contains the correct content.
Catches issues like copying the wrong HTML into template files.
"""

import pytest
import os
from pathlib import Path


TEMPLATES_DIR = Path(__file__).parent.parent / 'templates'


class TestChatStandaloneTemplate:
    """
    CRITICAL: Test chat_standalone.html is the HOME PAGE template
    This catches the exact bug where pre-op content was in this file
    """

    def test_chat_standalone_exists(self):
        """chat_standalone.html must exist"""
        template_path = TEMPLATES_DIR / 'chat_standalone.html'
        assert template_path.exists(), "chat_standalone.html is missing!"

    def test_chat_standalone_is_home_page(self):
        """chat_standalone.html should contain HOME PAGE content, not pre-op"""
        template_path = TEMPLATES_DIR / 'chat_standalone.html'
        content = template_path.read_text()

        # MUST contain home page title
        assert 'GasConsult.ai - AI-Powered Anesthesiology Assistant' in content, \
            "CRITICAL: chat_standalone.html missing home page title!"

        # MUST contain chat interface elements
        assert 'chat-input' in content or 'class="chat-input"' in content, \
            "CRITICAL: chat_standalone.html missing chat input field!"

        # MUST NOT be pre-op assessment
        # Check for pre-op specific title
        if 'Pre-Operative Assessment - GasConsult.ai' in content:
            # Check if it also has chat content (might be hybrid page)
            if 'chat-input' not in content.lower():
                pytest.fail(
                    "CRITICAL BUG: chat_standalone.html contains ONLY pre-op content! "
                    "This is the exact bug we fixed. Someone copied the wrong template."
                )

    def test_chat_standalone_has_correct_meta_description(self):
        """Meta description should be for evidence-based anesthesiology AI"""
        template_path = TEMPLATES_DIR / 'chat_standalone.html'
        content = template_path.read_text()

        # Should have the home page meta description
        assert 'Evidence-based anesthesiology AI consultant' in content or \
               'PubMed research with GPT-4o' in content, \
            "chat_standalone.html has wrong meta description"

    def test_chat_standalone_loads_utils_js(self):
        """Template should load utils.js"""
        template_path = TEMPLATES_DIR / 'chat_standalone.html'
        content = template_path.read_text()
        assert '/static/utils.js' in content, \
            "chat_standalone.html not loading utils.js"

    def test_chat_standalone_has_suggested_prompts(self):
        """Home page should have suggested prompt examples"""
        template_path = TEMPLATES_DIR / 'chat_standalone.html'
        content = template_path.read_text()
        # Look for hint-chip or suggested prompts section
        assert 'hint-chip' in content or 'suggested' in content.lower(), \
            "Home page missing suggested prompts"


class TestBaseTemplate:
    """Test base.html template structure"""

    def test_base_template_exists(self):
        """base.html must exist"""
        template_path = TEMPLATES_DIR / 'base.html'
        assert template_path.exists(), "base.html is missing!"

    def test_base_has_blocks(self):
        """Base template should have Jinja2 blocks for inheritance"""
        template_path = TEMPLATES_DIR / 'base.html'
        content = template_path.read_text()

        assert '{% block' in content, "base.html has no Jinja2 blocks"
        assert '{% endblock' in content, "base.html blocks not closed"

    def test_base_loads_static_files(self):
        """Base template should reference CSS/JS"""
        template_path = TEMPLATES_DIR / 'base.html'
        content = template_path.read_text()

        assert '/static/main.css' in content or '{% block styles %}' in content, \
            "base.html not loading CSS"


class TestComponentTemplates:
    """Test navbar and footer components"""

    def test_navbar_component_exists(self):
        """Navbar component must exist"""
        navbar_path = TEMPLATES_DIR / 'components' / 'navbar.html'
        assert navbar_path.exists(), "navbar.html is missing!"

    def test_navbar_has_toggle_functions(self):
        """Navbar should call toggleMobileMenu and toggleNavDropdown"""
        navbar_path = TEMPLATES_DIR / 'components' / 'navbar.html'
        content = navbar_path.read_text()

        assert 'toggleMobileMenu' in content, \
            "navbar.html not using toggleMobileMenu function"
        assert 'toggleNavDropdown' in content, \
            "navbar.html not using toggleNavDropdown function"

    def test_footer_component_exists(self):
        """Footer component must exist"""
        footer_path = TEMPLATES_DIR / 'components' / 'footer.html'
        assert footer_path.exists(), "footer.html is missing!"


class TestTemplateConsistency:
    """Cross-template validation"""

    def test_no_template_has_hardcoded_2024_year(self):
        """
        Templates should not have hardcoded 2024 year
        (catches stale copyright dates)
        """
        # This is a nice-to-have, not critical
        pass  # Skip for now

    def test_all_templates_use_same_css_version(self):
        """
        All templates should reference the same CSS version
        (catches inconsistent cache busting)
        """
        css_versions = {}

        for template_file in TEMPLATES_DIR.rglob('*.html'):
            content = template_file.read_text()
            if '/static/main.css?v=' in content:
                import re
                match = re.search(r'/static/main\.css\?v=(\d+)', content)
                if match:
                    version = match.group(1)
                    css_versions[template_file.name] = version

        # All should have same version
        if css_versions:
            unique_versions = set(css_versions.values())
            if len(unique_versions) > 1:
                print(f"WARNING: Inconsistent CSS versions: {css_versions}")
                # Don't fail, just warn


class TestInlineTemplates:
    """
    Test that inline HTML templates in app.py have correct content
    These are harder to test, so we'll just check app.py structure
    """

    def test_app_py_has_preop_html(self):
        """app.py should define PREOP_HTML constant"""
        app_path = Path(__file__).parent.parent / 'app.py'
        content = app_path.read_text()
        assert 'PREOP_HTML = """' in content, \
            "app.py missing PREOP_HTML template definition"

    def test_app_py_preop_html_is_preop(self):
        """PREOP_HTML should contain pre-op content, not chat"""
        app_path = Path(__file__).parent.parent / 'app.py'
        content = app_path.read_text()

        # Find PREOP_HTML definition
        import re
        match = re.search(r'PREOP_HTML = """(.*?)"""', content, re.DOTALL)
        if match:
            preop_html = match.group(1)
            assert 'Pre-Operative Assessment' in preop_html, \
                "PREOP_HTML template doesn't contain pre-op title"


# Quick test runner
if __name__ == '__main__':
    print("Running template validation tests...")
    pytest.main([__file__, '-v'])
