/**
 * GasConsult.ai - Shared JavaScript Utilities
 * Week 2 Performance Optimization: Consolidated duplicate form handlers
 * Reduces code duplication and improves maintainability
 */

/**
 * Initialize loading overlay for form submissions
 * Automatically shows loading overlay when forms are submitted
 * @param {string} formSelector - CSS selector for the form (default: 'form[method="POST"]')
 * @param {string} overlayId - ID of the loading overlay element (default: 'loadingOverlay')
 * @param {function} validateFn - Optional validation function, return true to show loading
 */
function initFormLoadingOverlay(formSelector = 'form[method="POST"]', overlayId = 'loadingOverlay', validateFn = null) {
    document.addEventListener('DOMContentLoaded', function() {
        const form = document.querySelector(formSelector);
        const loadingOverlay = document.getElementById(overlayId);

        if (form && loadingOverlay) {
            form.addEventListener('submit', function(e) {
                // If validation function provided, check it first
                if (validateFn && typeof validateFn === 'function') {
                    if (!validateFn(form)) {
                        return; // Don't show loading if validation fails
                    }
                }

                // Show loading overlay
                loadingOverlay.classList.add('active');
            });
        }
    });
}

/**
 * Close dropdowns when clicking outside
 * Reusable dropdown management for navigation and other dropdowns
 */
function initDropdownClickOutside() {
    document.addEventListener('click', function(e) {
        const dropdowns = document.querySelectorAll('.nav-dropdown-menu');
        dropdowns.forEach(dropdown => {
            if (!dropdown.parentElement.contains(e.target)) {
                dropdown.classList.remove('show');
            }
        });
    });
}

/**
 * Fill query from hint chips (homepage suggestions)
 * Called by onclick handlers on hint-chip elements
 * @param {Event} event - Click event from hint chip
 */
function fillQuery(event) {
    console.log('[fillQuery] Called with event:', event);
    const chip = event.currentTarget;
    const text = chip.textContent.trim();
    console.log('[fillQuery] Text:', text);
    const textarea = document.querySelector('.chat-input');
    console.log('[fillQuery] Textarea found:', !!textarea);
    if (textarea) {
        textarea.value = text;
        textarea.focus();
        console.log('[fillQuery] Value set successfully');
    } else {
        console.error('[fillQuery] No textarea with class .chat-input found');
    }
}

// Make fillQuery available globally for inline onclick handlers
window.fillQuery = fillQuery;

/**
 * Initialize all common UI utilities
 * Call this once per page for default behavior
 */
function initCommonUtilities() {
    initFormLoadingOverlay();
    initDropdownClickOutside();
}

// Auto-initialize if not explicitly disabled
if (typeof window.DISABLE_AUTO_INIT === 'undefined') {
    initCommonUtilities();
}
