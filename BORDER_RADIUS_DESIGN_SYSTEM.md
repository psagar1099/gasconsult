# Border-Radius Design System Implementation Guide

## Overview

This document outlines the standardization of border-radius values across the GasConsult.ai application using CSS custom properties (variables). The goal is to replace 403 inconsistent border-radius declarations across 18 HTML templates with a cohesive, maintainable design system.

## Design Tokens Created

**Location:** `/home/user/gasconsult/app.py` - Lines 7601-7612

```python
SHARED_DESIGN_TOKENS = """
        :root {
            /* Border Radius Scale */
            --radius-sm: 8px;     /* Small: badges, chips, small buttons */
            --radius-md: 12px;    /* Medium: buttons, inputs, small cards */
            --radius-lg: 16px;    /* Large: cards, panels, navigation */
            --radius-xl: 20px;    /* Extra large: feature cards, main containers */
            --radius-2xl: 24px;   /* Huge: hero sections, large panels */
            --radius-full: 100px; /* Pills, fully rounded elements */
        }
"""
```

## Current State Analysis

### Total Scope
- **Total border-radius declarations:** 403
- **Unique values:** 18 different pixel values
- **Templates affected:** 18 HTML templates
- **Shared CSS components:** 10 border-radius declarations in `SHARED_NAV_CSS`

### Value Distribution

| Current Value | Occurrences | Proposed Variable | Use Case |
|---------------|-------------|-------------------|----------|
| `12px` | 93 | `--radius-md` | Primary buttons, inputs, small cards |
| `16px` | 43 | `--radius-lg` | Large cards, panels, navigation |
| `20px` | 35 | `--radius-xl` | Feature cards, main containers |
| `8px` | 33 | `--radius-sm` | Badges, chips, small buttons |
| `14px` | 26 | `--radius-md` | **Standardize to 12px** |
| `10px` | 24 | `--radius-sm` | **Standardize to 8px** |
| `24px` | 22 | `--radius-2xl` | Hero sections, large panels |
| `100px` | 13 | `--radius-full` | Pills, fully rounded elements |
| `4px` | 11 | `--radius-sm` | **Reduce to 8px** (or create `--radius-xs: 4px`) |
| `6px` | 8 | `--radius-sm` | **Standardize to 8px** |
| `28px` | 7 | `--radius-2xl` | **Standardize to 24px** |
| `1px` | 7 | `--radius-sm` | **Very small, standardize to 8px** |
| `18px` | 7 | `--radius-lg` | **Standardize to 16px** |
| `2px` | 4 | `--radius-sm` | **Standardize to 8px** |
| `32px` | 3 | `--radius-2xl` | **Standardize to 24px** |
| `999px` | 2 | `--radius-full` | **Use 100px for consistency** |
| `3px` | 1 | `--radius-sm` | **Standardize to 8px** |
| `13px` | 1 | `--radius-md` | **Standardize to 12px** |

### Standardization Summary

**Primary Values (Keep As-Is):**
- 8px → `--radius-sm` (33 + 81 standardized = 114 total)
- 12px → `--radius-md` (93 + 27 standardized = 120 total)
- 16px → `--radius-lg` (43 + 7 standardized = 50 total)
- 20px → `--radius-xl` (35, keep as-is)
- 24px → `--radius-2xl` (22 + 13 standardized = 35 total)
- 100px → `--radius-full` (13 + 2 standardized = 15 total)

**Values to Standardize:**
- 1px, 2px, 3px, 4px, 6px, 10px → `--radius-sm` (8px)
- 13px, 14px → `--radius-md` (12px)
- 18px → `--radius-lg` (16px)
- 28px, 32px → `--radius-2xl` (24px)
- 999px → `--radius-full` (100px)

**Total Changes:**
- 322 values remain unchanged (just convert to CSS variable)
- 81 values will be slightly adjusted during standardization

## Templates by Priority (Border-Radius Usage)

### High Priority (50+ occurrences)
1. **PREOP_HTML** - 57 occurrences (Line 1597)
2. **HTML** (Homepage) - 53 occurrences (Line 7614)

### Medium Priority (20-50 occurrences)
3. **HYPOTENSION_HTML** - 36 occurrences (Line 20700)
4. **LIBRARY_HTML** - 31 occurrences (Line 10686)
5. **INFORMED_CONSENT_HTML** - 29 occurrences (Line 23374)
6. **SHARED_RESPONSE_HTML** - 25 occurrences (Line 11653)
7. **DIFFICULT_AIRWAY_HTML** - 25 occurrences (Line 24891)
8. **CRISIS_HTML** - 24 occurrences (Line 15061)
9. **CALCULATORS_HTML** - 23 occurrences (Line 18639)

### Low Priority (<20 occurrences)
10. **QUICK_DOSE_HTML** - 19 occurrences (Line 16692)
11. **ADMIN_DASHBOARD_HTML** - 15 occurrences (Line 29758)
12. **EVIDENCE_HTML** - 10 occurrences (Line 14072)
13. **REGISTER_HTML** - 9 occurrences (Line 28896)
14. **LOGIN_HTML** - 8 occurrences (Line 28387)
15. **RESEND_VERIFICATION_HTML** - 7 occurrences (Line 29402)
16. **TERMS_HTML** - 6 occurrences (Line 12469)
17. **PRIVACY_POLICY_HTML** - 5 occurrences (Line 13317)
18. **PRICING_HTML** - 0 occurrences (Line 30587)

### Shared Components
**SHARED_NAV_CSS** - 10 occurrences (Lines 7286-7600)

## Implementation Strategy

### Phase 1: Add Design Tokens to Templates (COMPLETED ✅)

**Status:** Design tokens constant created at line 7601-7612

### Phase 2: Inject Tokens into Each Template

Add `{SHARED_DESIGN_TOKENS}` to the `<style>` section of each template, immediately after the opening `<style>` tag:

```html
<style>
    {SHARED_DESIGN_TOKENS}

    /* Existing color variables */
    :root {
        --white: #FFFFFF;
        --gray-50: #F9FAFB;
        /* ... */
    }

    /* Rest of styles */
</style>
```

**Templates to Update (in order):**
1. HTML (Homepage)
2. PREOP_HTML
3. HYPOTENSION_HTML
4. LIBRARY_HTML
5. INFORMED_CONSENT_HTML
6. SHARED_RESPONSE_HTML
7. DIFFICULT_AIRWAY_HTML
8. CRISIS_HTML
9. CALCULATORS_HTML
10. QUICK_DOSE_HTML
11. ADMIN_DASHBOARD_HTML
12. EVIDENCE_HTML
13. REGISTER_HTML
14. LOGIN_HTML
15. RESEND_VERIFICATION_HTML
16. TERMS_HTML
17. PRIVACY_POLICY_HTML

**Note:** PRICING_HTML has no border-radius declarations, so it can be skipped or included for consistency.

### Phase 3: Replace Values (Gradual, Template-by-Template)

**Recommended Approach:** Start with high-impact, low-risk templates first.

#### Batch 1: Shared Components (Test First)
- **Target:** `SHARED_NAV_CSS` (10 occurrences)
- **Why First:** Changes affect all pages, catch any issues early
- **Risk:** Medium (global impact)
- **Test:** Homepage, all tool pages

**Example Replacement:**
```css
/* Before */
.nav-inner {
    border-radius: 16px;
}

/* After */
.nav-inner {
    border-radius: var(--radius-lg);
}
```

#### Batch 2: High-Volume Templates (Visual Consistency)
- **Targets:** PREOP_HTML (57), HTML/Homepage (53)
- **Why Next:** Highest usage, most visual impact
- **Risk:** Low (isolated to specific pages)
- **Test:** Homepage, preop tool

#### Batch 3: Medium-Volume Templates
- **Targets:** HYPOTENSION_HTML (36), LIBRARY_HTML (31), INFORMED_CONSENT_HTML (29), etc.
- **Why:** Moderate impact, manageable scope
- **Risk:** Low
- **Test:** Each tool individually

#### Batch 4: Low-Volume Templates
- **Targets:** AUTH pages (LOGIN_HTML, REGISTER_HTML), TERMS_HTML, etc.
- **Why Last:** Lowest usage, minimal visual changes
- **Risk:** Very low
- **Test:** Auth flows, static pages

### Phase 4: Validation & Cleanup

**After Each Batch:**
1. ✅ Visual regression testing (screenshot comparison)
2. ✅ Browser compatibility check (Chrome, Firefox, Safari)
3. ✅ Mobile responsiveness verification
4. ✅ Search for any missed `border-radius:` declarations

**Final Cleanup:**
```bash
# Verify no hardcoded values remain
grep -n "border-radius: [0-9]\+px" app.py | wc -l
# Should be 0 after completion
```

## Replacement Patterns (Search & Replace)

### Safe Automated Replacements (Exact Matches)

```python
# Use Edit tool with replace_all=True for these exact patterns:

# Primary values (no standardization needed)
"border-radius: 12px;" → "border-radius: var(--radius-md);"
"border-radius: 16px;" → "border-radius: var(--radius-lg);"
"border-radius: 20px;" → "border-radius: var(--radius-xl);"
"border-radius: 8px;" → "border-radius: var(--radius-sm);"
"border-radius: 24px;" → "border-radius: var(--radius-2xl);"
"border-radius: 100px;" → "border-radius: var(--radius-full);"
```

### Careful Manual Replacements (Standardization Required)

These require visual review as the actual pixel value will change:

```python
# Standardize to --radius-sm (8px)
"border-radius: 10px;" → "border-radius: var(--radius-sm);"  # -2px change
"border-radius: 6px;" → "border-radius: var(--radius-sm);"   # +2px change
"border-radius: 4px;" → "border-radius: var(--radius-sm);"   # +4px change
"border-radius: 1px;" → "border-radius: var(--radius-sm);"   # +7px change
"border-radius: 2px;" → "border-radius: var(--radius-sm);"   # +6px change
"border-radius: 3px;" → "border-radius: var(--radius-sm);"   # +5px change

# Standardize to --radius-md (12px)
"border-radius: 14px;" → "border-radius: var(--radius-md);"  # -2px change
"border-radius: 13px;" → "border-radius: var(--radius-md);"  # -1px change

# Standardize to --radius-lg (16px)
"border-radius: 18px;" → "border-radius: var(--radius-lg);"  # -2px change

# Standardize to --radius-2xl (24px)
"border-radius: 28px;" → "border-radius: var(--radius-2xl);" # -4px change
"border-radius: 32px;" → "border-radius: var(--radius-2xl);" # -8px change

# Standardize to --radius-full (100px)
"border-radius: 999px;" → "border-radius: var(--radius-full);" # -899px (visual impact minimal for pill shapes)
```

## Alternative: Create --radius-xs for Small Values

If the 1px-4px values serve a specific design purpose (subtle borders, etc.), consider adding:

```css
:root {
    --radius-xs: 4px;  /* Extra small: minimal rounding */
    --radius-sm: 8px;
    /* ... rest unchanged */
}
```

Then:
```python
"border-radius: 4px;" → "border-radius: var(--radius-xs);"
"border-radius: 1px;" → "border-radius: var(--radius-xs);"  # +3px change
"border-radius: 2px;" → "border-radius: var(--radius-xs);"  # +2px change
"border-radius: 3px;" → "border-radius: var(--radius-xs);"  # +1px change
```

This reduces visual impact but adds one more token to maintain.

## Testing Checklist

### Visual Testing
- [ ] Homepage hero section (--radius-xl, --radius-2xl)
- [ ] Navigation bar (--radius-lg)
- [ ] Primary buttons (--radius-md)
- [ ] Input fields (--radius-md)
- [ ] Cards/panels (--radius-lg, --radius-xl)
- [ ] Badges/chips (--radius-sm)
- [ ] Pills/tags (--radius-full)

### Functional Testing
- [ ] Click targets not affected by border-radius changes
- [ ] No visual overflow issues
- [ ] Hover states work correctly
- [ ] Focus indicators visible

### Browser/Device Testing
- [ ] Chrome (desktop + mobile)
- [ ] Firefox (desktop + mobile)
- [ ] Safari (desktop + iOS)
- [ ] Edge (desktop)

## Benefits of This Design System

1. **Consistency:** Unified visual language across all 18 templates
2. **Maintainability:** Single source of truth for border-radius values
3. **Scalability:** Easy to adjust global rounding with one change
4. **Accessibility:** Predictable visual patterns improve UX
5. **Developer Experience:** Clear semantic naming reduces guesswork
6. **Performance:** No impact (CSS variables are fast)
7. **Future-Proof:** Easy to add dark mode or theme variations

## Rollback Plan

If issues arise during implementation:

1. **Immediate:** Remove `{SHARED_DESIGN_TOKENS}` from affected template
2. **Quick Fix:** Revert specific template to previous commit
3. **Full Rollback:** Git revert to pre-design-system state

```bash
# Find which templates have been updated
git diff HEAD --name-only | grep app.py

# Revert specific changes
git checkout HEAD~1 app.py

# Or use git revert for specific commits
git revert <commit-hash>
```

## Next Steps

1. ✅ **COMPLETED:** Create `SHARED_DESIGN_TOKENS` constant
2. **TODO:** Inject `{SHARED_DESIGN_TOKENS}` into each template's `<style>` section
3. **TODO:** Start with Batch 1 (SHARED_NAV_CSS - 10 occurrences)
4. **TODO:** Proceed through Batches 2-4 based on priority
5. **TODO:** Final validation and cleanup
6. **TODO:** Update CHANGELOG.md with design system implementation

## Timeline Estimate

- **Phase 2 (Inject tokens):** 1-2 hours (18 templates)
- **Phase 3 Batch 1:** 30 minutes (shared components)
- **Phase 3 Batch 2:** 2 hours (high-volume templates)
- **Phase 3 Batch 3:** 3 hours (medium-volume templates)
- **Phase 3 Batch 4:** 1 hour (low-volume templates)
- **Phase 4 (Validation):** 2 hours (testing)

**Total:** ~9-10 hours of work

## Questions/Decisions Needed

1. **Keep --radius-xs (4px)?** Or standardize 1-4px to --radius-sm (8px)?
2. **Visual acceptance:** Is +7px change (1px → 8px) acceptable for very small elements?
3. **Testing scope:** Manual visual testing or automated screenshot diffing?
4. **Deployment strategy:** All at once or incremental template-by-template deployment?

---

**Document Version:** 1.0
**Created:** 2025-12-31
**Last Updated:** 2025-12-31
**Author:** Claude Code Agent
**Status:** Design Tokens Created, Ready for Phase 2
