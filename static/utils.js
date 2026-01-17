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
 * Citation Preview Feature (Week 3 Innovation)
 * Shows paper abstracts in a popup when clicking citation numbers
 */
function initCitationPreviews() {
    let activePopup = null;
    let activeOverlay = null;

    function closePopup() {
        if (activePopup) {
            activePopup.remove();
            activePopup = null;
        }
        if (activeOverlay) {
            activeOverlay.remove();
            activeOverlay = null;
        }
    }

    function showCitationPreview(trigger, event) {
        event.preventDefault();
        event.stopPropagation();

        // Close any existing popup
        closePopup();

        const abstract = trigger.getAttribute('data-abstract');
        const pmid = trigger.getAttribute('data-pmid');
        const title = trigger.getAttribute('data-title');

        if (!abstract || !pmid) {
            console.error('Missing abstract or PMID data');
            return;
        }

        // Create overlay
        activeOverlay = document.createElement('div');
        activeOverlay.className = 'citation-preview-overlay';
        activeOverlay.addEventListener('click', closePopup);

        // Create popup
        activePopup = document.createElement('div');
        activePopup.className = 'citation-preview-popup';

        // Build popup HTML
        activePopup.innerHTML = `
            <div class="citation-preview-header">
                <div class="citation-preview-title">${title || 'Citation'}</div>
                <div class="citation-preview-close" title="Close">
                    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                        <line x1="18" y1="6" x2="6" y2="18"></line>
                        <line x1="6" y1="6" x2="18" y2="18"></line>
                    </svg>
                </div>
            </div>
            <a href="https://pubmed.ncbi.nlm.nih.gov/${pmid}/" target="_blank" rel="noopener noreferrer" class="citation-preview-pmid">
                <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <path d="M18 13v6a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V8a2 2 0 0 1 2-2h6"></path>
                    <polyline points="15 3 21 3 21 9"></polyline>
                    <line x1="10" y1="14" x2="21" y2="3"></line>
                </svg>
                View on PubMed (${pmid})
            </a>
            <div class="citation-preview-abstract-label">Abstract</div>
            <div class="citation-preview-abstract">${abstract}</div>
        `;

        // Add close button handler
        const closeBtn = activePopup.querySelector('.citation-preview-close');
        closeBtn.addEventListener('click', closePopup);

        // Position popup
        document.body.appendChild(activeOverlay);
        document.body.appendChild(activePopup);

        // Center popup on desktop, or use smart positioning
        if (window.innerWidth >= 769) {
            // Desktop: position near the trigger
            const rect = trigger.getBoundingClientRect();
            const popupWidth = 500;
            const popupHeight = Math.min(400, activePopup.offsetHeight);

            let left = rect.right + 12;
            let top = rect.top;

            // Keep popup in viewport
            if (left + popupWidth > window.innerWidth - 20) {
                left = rect.left - popupWidth - 12;
            }
            if (left < 20) {
                left = 20;
            }
            if (top + popupHeight > window.innerHeight - 20) {
                top = window.innerHeight - popupHeight - 20;
            }
            if (top < 20) {
                top = 20;
            }

            activePopup.style.left = `${left}px`;
            activePopup.style.top = `${top}px`;
        }
        // Mobile: CSS handles centering via transform

        // Prevent body scroll when popup is open
        document.body.style.overflow = 'hidden';
    }

    // Delegate click events to citation triggers
    document.addEventListener('click', function(e) {
        const trigger = e.target.closest('.citation-preview-trigger');
        if (trigger) {
            showCitationPreview(trigger, e);
        }
    });

    // Close popup on Escape key
    document.addEventListener('keydown', function(e) {
        if (e.key === 'Escape' && activePopup) {
            closePopup();
        }
    });

    // Restore body scroll when popup closes
    const originalClosePopup = closePopup;
    closePopup = function() {
        originalClosePopup();
        document.body.style.overflow = '';
    };
}

/**
 * PDF Export Feature (Week 3 Innovation)
 * Exports the current conversation as a PDF with citations
 */
function exportConversationPDF() {
    console.log('[exportConversationPDF] Starting export...');

    const exportBtn = document.querySelector('.export-pdf-btn');
    if (!exportBtn) {
        console.error('[exportConversationPDF] Export button not found');
        return;
    }

    // Show loading state
    exportBtn.classList.add('exporting');
    const originalHTML = exportBtn.innerHTML;
    exportBtn.innerHTML = `
        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
            <circle cx="12" cy="12" r="10"></circle>
        </svg>
        <span>Exporting...</span>
    `;

    // Collect conversation messages from the page
    const conversationHistory = [];
    const messageElements = document.querySelectorAll('.message');

    messageElements.forEach(function(msgEl) {
        const role = msgEl.classList.contains('user-message') ? 'user' : 'assistant';
        const contentEl = msgEl.querySelector('.message-content, .ai-response-content');
        const content = contentEl ? contentEl.innerHTML : '';

        const message = {
            role: role,
            content: content,
            references: [],
            num_papers: 0,
            evidence_strength: {}
        };

        // Extract references if this is an assistant message
        if (role === 'assistant') {
            const referenceItems = msgEl.querySelectorAll('.reference-item');
            referenceItems.forEach(function(refEl) {
                const titleEl = refEl.querySelector('.reference-link');
                const metaEl = refEl.querySelector('.reference-meta');
                const pmidMatch = titleEl ? titleEl.href.match(/\/(\d+)\/$/) : null;

                if (titleEl && metaEl && pmidMatch) {
                    const metaParts = metaEl.textContent.split(' - ');
                    const journalYear = metaParts[1] ? metaParts[1].split(', ') : ['Unknown', 'N/A'];

                    message.references.push({
                        title: titleEl.textContent.trim(),
                        authors: metaParts[0] || 'Unknown',
                        journal: journalYear[0] || 'Unknown',
                        year: journalYear[1] || 'N/A',
                        pmid: pmidMatch[1]
                    });
                }
            });

            message.num_papers = message.references.length;

            // Extract evidence strength if available
            const evidenceBadge = msgEl.querySelector('.evidence-badge');
            if (evidenceBadge) {
                const levelMatch = evidenceBadge.textContent.match(/✓\s*(\w+)\s*Confidence/);
                const confidenceMatch = evidenceBadge.textContent.match(/(\d+)%/);
                const studiesMatch = evidenceBadge.textContent.match(/(\d+)\s*studies/);

                message.evidence_strength = {
                    level: levelMatch ? levelMatch[1] : 'Unknown',
                    confidence_percentage: confidenceMatch ? parseInt(confidenceMatch[1]) : 0
                };
                message.num_papers = studiesMatch ? parseInt(studiesMatch[1]) : message.references.length;
            }
        }

        conversationHistory.push(message);
    });

    console.log('[exportConversationPDF] Collected', conversationHistory.length, 'messages');

    if (conversationHistory.length === 0) {
        alert('No conversation to export');
        exportBtn.classList.remove('exporting');
        exportBtn.innerHTML = originalHTML;
        return;
    }

    // Send to server for PDF generation
    fetch('/export-pdf', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({
            conversation_history: conversationHistory
        })
    })
    .then(function(response) {
        if (!response.ok) {
            throw new Error('Failed to generate PDF: ' + response.statusText);
        }
        return response.blob();
    })
    .then(function(blob) {
        // Create download link
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'gasconsult-conversation-' + new Date().toISOString().slice(0, 19).replace(/[:-]/g, '') + '.pdf';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        window.URL.revokeObjectURL(url);

        console.log('[exportConversationPDF] PDF downloaded successfully');

        // Reset button
        exportBtn.classList.remove('exporting');
        exportBtn.innerHTML = originalHTML;
    })
    .catch(function(error) {
        console.error('[exportConversationPDF] Error:', error);
        alert('Failed to export PDF: ' + error.message);

        // Reset button
        exportBtn.classList.remove('exporting');
        exportBtn.innerHTML = originalHTML;
    });
}

// Make exportConversationPDF available globally
window.exportConversationPDF = exportConversationPDF;

/**
 * Evidence Quality Chart Feature (Week 3 Innovation)
 * Shows interactive bar chart of study type breakdown
 */
function showEvidenceChart(badge, event) {
    event.preventDefault();
    event.stopPropagation();

    // Get data from badge
    const guidelines = parseInt(badge.getAttribute('data-guidelines')) || 0;
    const metaAnalyses = parseInt(badge.getAttribute('data-meta-analyses')) || 0;
    const systematicReviews = parseInt(badge.getAttribute('data-systematic-reviews')) || 0;
    const rcts = parseInt(badge.getAttribute('data-rcts')) || 0;
    const observational = parseInt(badge.getAttribute('data-observational')) || 0;
    const total = parseInt(badge.getAttribute('data-total')) || 0;
    const level = badge.getAttribute('data-level') || 'Unknown';
    const confidence = badge.getAttribute('data-confidence') || '0';

    if (total === 0) {
        return; // No data to show
    }

    // Create overlay
    const overlay = document.createElement('div');
    overlay.className = 'citation-preview-overlay';

    // Create chart popup
    const popup = document.createElement('div');
    popup.className = 'evidence-chart-popup';

    // Calculate percentages
    const studyTypes = [
        { label: 'Guidelines', count: guidelines, cssClass: 'guideline' },
        { label: 'Meta-Analyses', count: metaAnalyses, cssClass: 'meta-analysis' },
        { label: 'Systematic Reviews', count: systematicReviews, cssClass: 'systematic-review' },
        { label: 'RCTs', count: rcts, cssClass: 'rct' },
        { label: 'Observational', count: observational, cssClass: 'observational' }
    ];

    // Filter out zero counts and sort by count descending
    const nonZeroStudies = studyTypes.filter(s => s.count > 0).sort((a, b) => b.count - a.count);

    // Build chart HTML
    let chartHTML = `
        <div class="evidence-chart-header">
            <div class="evidence-chart-title">Evidence Quality Breakdown</div>
            <div class="evidence-chart-close" title="Close">
                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <line x1="18" y1="6" x2="6" y2="18"></line>
                    <line x1="6" y1="6" x2="18" y2="18"></line>
                </svg>
            </div>
        </div>
        <div class="evidence-chart-bars">
    `;

    nonZeroStudies.forEach(function(study) {
        const percentage = Math.round((study.count / total) * 100);
        chartHTML += `
            <div class="chart-bar-row">
                <div class="chart-bar-label">
                    <span>${study.label}</span>
                    <span class="chart-bar-count">${study.count}</span>
                </div>
                <div class="chart-bar-track">
                    <div class="chart-bar-fill ${study.cssClass}" style="width: ${percentage}%;">
                        <span class="chart-bar-percentage">${percentage}%</span>
                    </div>
                </div>
            </div>
        `;
    });

    // Add summary
    const topStudyType = nonZeroStudies[0];
    chartHTML += `
        </div>
        <div class="evidence-chart-summary">
            <strong>${level} Evidence Quality (${confidence}% confidence)</strong><br/>
            Based on ${total} ${total === 1 ? 'study' : 'studies'}, primarily ${topStudyType.label.toLowerCase()} (${topStudyType.count}/${total}).
            ${level === 'High' ? 'Strong evidence from high-quality sources.' : (level === 'Moderate' ? 'Moderate evidence - consider individual patient factors.' : 'Limited evidence - use caution and clinical judgment.')}
        </div>
    `;

    popup.innerHTML = chartHTML;

    // Add close handler
    const closeBtn = popup.querySelector('.evidence-chart-close');
    closeBtn.addEventListener('click', function() {
        popup.remove();
        overlay.remove();
        document.body.style.overflow = '';
    });

    // Close on overlay click
    overlay.addEventListener('click', function() {
        popup.remove();
        overlay.remove();
        document.body.style.overflow = '';
    });

    // Close on Escape key
    const escapeHandler = function(e) {
        if (e.key === 'Escape') {
            popup.remove();
            overlay.remove();
            document.body.style.overflow = '';
            document.removeEventListener('keydown', escapeHandler);
        }
    };
    document.addEventListener('keydown', escapeHandler);

    // Add to DOM
    document.body.appendChild(overlay);
    document.body.appendChild(popup);

    // Position popup on desktop
    if (window.innerWidth >= 769) {
        const rect = badge.getBoundingClientRect();
        const popupWidth = 450;
        const popupHeight = Math.min(400, popup.offsetHeight);

        let left = rect.right + 12;
        let top = rect.top;

        // Keep popup in viewport
        if (left + popupWidth > window.innerWidth - 20) {
            left = rect.left - popupWidth - 12;
        }
        if (left < 20) {
            left = 20;
        }
        if (top + popupHeight > window.innerHeight - 20) {
            top = window.innerHeight - popupHeight - 20;
        }
        if (top < 20) {
            top = 20;
        }

        popup.style.left = `${left}px`;
        popup.style.top = `${top}px`;
    }

    // Prevent body scroll
    document.body.style.overflow = 'hidden';

    console.log('[showEvidenceChart] Displayed chart with', total, 'studies');
}

// Make showEvidenceChart available globally
window.showEvidenceChart = showEvidenceChart;

/**
 * Initialize all common UI utilities
 * Call this once per page for default behavior
 */
function initCommonUtilities() {
    initFormLoadingOverlay();
    initDropdownClickOutside();
    initCitationPreviews();
}

// Auto-initialize if not explicitly disabled
if (typeof window.DISABLE_AUTO_INIT === 'undefined') {
    initCommonUtilities();
}
