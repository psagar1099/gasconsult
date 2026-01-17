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
