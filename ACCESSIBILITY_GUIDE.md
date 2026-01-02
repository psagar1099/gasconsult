# Accessibility Implementation Guide - GasConsult.ai

## Current State
- **ARIA Coverage:** 91 attributes across 33,259 lines (0.27%)
- **WCAG Compliance:** Partial Level A
- **Major Gaps:** Missing labels, poor focus management, no skip links (except homepage)

---

## Priority 1: Critical ARIA Labels (Immediate)

### Mobile Menu Button
**Current:** Missing `aria-label`, `aria-expanded`, `aria-controls`

**Fix:**
```html
<button 
    class="mobile-menu-btn" 
    onclick="toggleMobileMenu()"
    aria-label="Toggle navigation menu"
    aria-expanded="false"
    aria-controls="mobileMenu">
    ☰
</button>
```

### Dropdown Menus
**Current:** No `aria-haspopup`, `aria-expanded`

**Fix:**
```html
<button 
    class="nav-dropdown-toggle"
    aria-haspopup="true"
    aria-expanded="false"
    aria-label="More options menu">
    More ▼
</button>
<div class="nav-dropdown-menu" role="menu">
    <a href="..." role="menuitem">Option</a>
</div>
```

### Form Inputs
**Current:** Some inputs lack proper label associations

**Fix:**
```html
<label for="email">Email Address</label>
<input 
    type="email" 
    id="email" 
    name="email"
    aria-required="true"
    aria-describedby="email-hint">
<span id="email-hint">We'll never share your email</span>
```

---

## Priority 2: Skip Links (All Pages)

**Current:** Only homepage has skip link  
**Needed:** All 18 templates

**Implementation:**
```html
<!-- At top of <body>, before navigation -->
<a href="#main-content" class="skip-to-content">Skip to main content</a>

<style>
.skip-to-content {
    position: absolute;
    left: -9999px;
    z-index: 999;
    padding: 1em;
    background: var(--blue-600);
    color: white;
    text-decoration: none;
}

.skip-to-content:focus {
    left: 50%;
    transform: translateX(-50%);
    top: 1em;
}
</style>

<!-- Add id to main content -->
<main id="main-content">...</main>
```

---

## Priority 3: Focus Management

### Visible Focus Indicators
```css
*:focus-visible {
    outline: 2px solid var(--blue-600);
    outline-offset: 2px;
    border-radius: 4px;
}

button:focus-visible,
a:focus-visible,
input:focus-visible {
    outline: 2px solid var(--blue-600);
    outline-offset: 2px;
}
```

### Focus Trap for Modals
```javascript
function trapFocus(element) {
    const focusableElements = element.querySelectorAll(
        'button, [href], input, select, textarea, [tabindex]:not([tabindex="-1"])'
    );
    const firstFocusable = focusableElements[0];
    const lastFocusable = focusableElements[focusableElements.length - 1];

    element.addEventListener('keydown', function(e) {
        if (e.key === 'Tab') {
            if (e.shiftKey) {
                if (document.activeElement === firstFocusable) {
                    lastFocusable.focus();
                    e.preventDefault();
                }
            } else {
                if (document.activeElement === lastFocusable) {
                    firstFocusable.focus();
                    e.preventDefault();
                }
            }
        }
    });
}
```

---

## Priority 4: Live Regions (Chat Interface)

**Current:** No `aria-live` for dynamic content

**Fix:**
```html
<div class="messages-container" 
     aria-live="polite" 
     aria-atomic="false"
     role="log">
    <!-- Chat messages appear here -->
</div>

<div class="streaming-indicator" 
     aria-live="assertive" 
     role="status">
    <!-- Loading states -->
</div>
```

---

## Priority 5: Form Validation

**Current:** No `aria-invalid`, `aria-describedby` for errors

**Fix:**
```html
<input 
    type="email" 
    id="email"
    aria-required="true"
    aria-invalid="false"
    aria-describedby="email-error">

<span id="email-error" role="alert" class="error-message">
    <!-- Error message appears here -->
</span>

<script>
function validateEmail(input) {
    const isValid = input.value.includes('@');
    input.setAttribute('aria-invalid', !isValid);
    
    const errorSpan = document.getElementById(input.id + '-error');
    if (!isValid) {
        errorSpan.textContent = 'Please enter a valid email address';
    } else {
        errorSpan.textContent = '';
    }
}
</script>
```

---

## Implementation Checklist

### Templates Requiring Updates (18 total)
- [ ] HTML (Homepage)
- [ ] PREOP_HTML
- [ ] LIBRARY_HTML
- [ ] SHARED_RESPONSE_HTML
- [ ] TERMS_HTML
- [ ] PRIVACY_POLICY_HTML
- [ ] EVIDENCE_HTML
- [ ] CRISIS_HTML
- [ ] QUICK_DOSE_HTML
- [ ] CALCULATORS_HTML
- [ ] HYPOTENSION_HTML
- [ ] INFORMED_CONSENT_HTML
- [ ] DIFFICULT_AIRWAY_HTML
- [ ] LOGIN_HTML
- [ ] REGISTER_HTML
- [ ] RESEND_VERIFICATION_HTML
- [ ] ADMIN_DASHBOARD_HTML
- [ ] PRICING_HTML

### Per-Template Tasks
1. Add skip link (`<a href="#main-content">`)
2. Add `id="main-content"` to main section
3. Add `aria-label` to mobile menu button
4. Add `aria-expanded` to all dropdowns
5. Add `aria-live` to dynamic content areas
6. Add `aria-invalid` to form inputs with validation
7. Add `role` attributes where semantic HTML isn't used

---

## Testing

### Manual Testing
1. **Keyboard Navigation:** Tab through entire page
2. **Screen Reader:** Test with NVDA/JAWS/VoiceOver
3. **Focus Visible:** Ensure all interactive elements show focus
4. **Skip Links:** Test skip-to-content functionality

### Automated Testing Tools
- **axe DevTools:** Browser extension for WCAG auditing
- **Lighthouse:** Built into Chrome DevTools
- **WAVE:** Web accessibility evaluation tool

**Target Metrics:**
- Lighthouse Accessibility Score: >90
- WCAG 2.1 Level AA compliance
- Zero critical issues in axe DevTools

---

## Resources

- [ARIA Authoring Practices Guide](https://www.w3.org/WAI/ARIA/apg/)
- [WebAIM WCAG Checklist](https://webaim.org/standards/wcag/checklist)
- [MDN ARIA Documentation](https://developer.mozilla.org/en-US/docs/Web/Accessibility/ARIA)

---

**Note:** This guide provides patterns and examples. Actual implementation across all 18 templates requires systematic updates to HTML structure and JavaScript event handlers.
