# RAG Pipeline Critical Audit
**Date:** 2025-11-27
**Pipeline:** PubMed Search → Paper Fetch → GPT-4o Synthesis

## 🔴 CRITICAL ISSUES FOUND

### 1. **Environment Variables Not Set (BLOCKING)**
**Location:** Line 41-44
**Issue:** If environment variables are not set in Render, the app will fail silently
```python
Entrez.email = os.getenv("ENTREZ_EMAIL", "your-email@example.com")  # Invalid default!
Entrez.api_key = os.getenv("ENTREZ_API_KEY", "")  # Empty API key!
openai_client = openai.OpenAI(api_key=os.getenv("OPENAI_API_KEY"))  # Will fail if not set
```

**Impact:** PubMed returns HTTP 400 errors when email is invalid
**Status:** ✅ User just added env vars to Render - should be fixed once deployed

---

### 2. **PubMed Date Range Format (POTENTIAL ISSUE)**
**Location:** Lines 7055, 7061
**Current Format:** `("2015/01/01"[PDAT] : "3000"[PDAT])`

**Concern:** This format with quotes and slashes may be deprecated or invalid in newer PubMed API
**Alternative Format:** `2015:3000[pdat]` (simpler, more modern)

**Recommendation:** Monitor error logs. If 400 errors persist after env vars are set, try the simpler format.

---

### 3. **Sort Parameter May Be Invalid**
**Location:** Lines 7070, 7082, 7401
**Current:** `sort="relevance"`

**Issue:** PubMed E-utilities documentation suggests valid values might be:
- `"pub_date"` (publication date)
- `"author"` (first author)
- `"relevance"` (may work, but verify)

**Recommendation:** If searches continue to fail, try removing `sort` parameter entirely or use `sort="pub_date"`

---

### 4. **Error Handling Suppresses Important Info**
**Location:** Lines 7074-7076, 7086-7088
```python
except Exception as e:
    print(f"[ERROR] PubMed search failed (anesthesiology): {e}")
    ids = []  # Silently continues with empty results
```

**Issue:** Catches all exceptions without proper logging or user feedback
**Better Approach:**
```python
except Exception as e:
    logger.error(f"PubMed search failed: {e}", exc_info=True)
    # Re-raise critical errors, only catch known recoverable ones
```

---

### 5. **Bare Except Clause Hides Paper Processing Errors**
**Location:** Line 7175
```python
except:
    continue  # Silently skips papers with errors
```

**Issue:** Completely suppresses errors when processing papers - could hide critical bugs
**Impact:** Papers might be silently dropped without any indication

---

### 6. **OpenAI API Key Validation Missing**
**Location:** Line 44
```python
openai_client = openai.OpenAI(api_key=os.getenv("OPENAI_API_KEY"))
```

**Issue:** If `OPENAI_API_KEY` is not set, returns None and OpenAI client will fail at runtime (not startup)
**Better:**
```python
api_key = os.getenv("OPENAI_API_KEY")
if not api_key:
    raise ValueError("OPENAI_API_KEY environment variable is required")
openai_client = openai.OpenAI(api_key=api_key)
```

---

### 7. **No Rate Limiting for PubMed API**
**Location:** Entrez calls throughout

**Issue:** NCBI requires:
- Max 3 requests/second without API key
- Max 10 requests/second with API key

**Current:** No rate limiting implemented - could get IP banned
**Mitigation:** Biopython's Entrez module has built-in rate limiting if `Entrez.email` is set properly

---

### 8. **Session Data Could Grow Unbounded**
**Location:** Lines 7232-7238
```python
session[f'stream_data_{request_id}'] = {
    'prompt': prompt,  # Could be very large
    'refs': refs,
    'num_papers': num_papers,
    'raw_query': raw_query
}
```

**Issue:** Session stores full prompts and references - could exceed session storage limits
**Impact:** Session corruption on large conversations

---

## 🟡 MEDIUM PRIORITY ISSUES

### 9. **Synonym Expansion Could Break MeSH Terms**
**Location:** Lines 7038-7049

**Example Issue:**
```python
q = q.replace("etomidate", '"etomidate"[MeSH Terms] OR etomidate')
```

If user searches "etomidate induction", this becomes:
```
"etomidate"[MeSH Terms] OR etomidate induction
```
This is malformed! Should be:
```
("etomidate"[MeSH Terms] OR etomidate) AND induction
```

**Recommendation:** Use word boundary matching or more sophisticated parsing

---

### 10. **Magic Numbers Throughout Code**
- Line 7161: `papers[:8]` - Why 8?
- Line 7166: `abstract[:600]` - Why 600?
- Line 7168: `[:3]` authors - Why 3?

**Issue:** No explanation for these limits - makes tuning difficult

---

## 🟢 LOW PRIORITY / STYLE ISSUES

### 11. **Inconsistent Error Messages**
Some errors return HTML, some plain text, some redirect

### 12. **Debug Print Statements in Production**
All the `print(f"[DEBUG] ...")` should use proper logger

---

## ✅ THINGS THAT ARE CORRECT

1. ✅ Using `Entrez.read()` to parse XML responses
2. ✅ Streaming GPT responses for better UX
3. ✅ Conversation context management
4. ✅ CSRF protection enabled
5. ✅ Input sanitization working
6. ✅ Temperature settings (0.1 for evidence, 0.2 for follow-ups)

---

## 🔧 IMMEDIATE ACTION ITEMS

1. **VERIFY ENV VARS ARE SET IN RENDER:**
   - ENTREZ_EMAIL
   - ENTREZ_API_KEY
   - OPENAI_API_KEY
   - FLASK_SECRET_KEY

2. **TEST AFTER DEPLOYMENT:**
   - Try "etomidate" search
   - Check logs for any 400 errors
   - Verify papers are being returned

3. **IF STILL FAILING:**
   - Try removing `sort="relevance"` parameter
   - Try simpler date format: `2015:3000[pdat]`
   - Add detailed error logging to identify exact issue

4. **MONITORING:**
   - Watch for OpenAI API rate limits
   - Monitor session storage size
   - Check PubMed API usage

---

## 📊 CODE QUALITY SCORE

| Category | Score | Notes |
|----------|-------|-------|
| Functionality | 7/10 | Works when env vars set |
| Error Handling | 4/10 | Too many bare excepts |
| Logging | 6/10 | Has logging but uses print() |
| Security | 9/10 | Good - env vars, sanitization, CSRF |
| Maintainability | 6/10 | Single 7400 line file |
| **Overall** | **6.4/10** | **Functional but needs hardening** |

---

## 🎯 ROOT CAUSE OF CURRENT ISSUE

**The HTTP 400 errors are caused by:**
1. Missing/invalid `ENTREZ_EMAIL` environment variable in Render
2. Missing/invalid `ENTREZ_API_KEY` environment variable in Render

Once these are set correctly in Render and the app redeploys, searches should work.

**How to verify it's fixed:**
```bash
# In Render logs, you should see:
[DEBUG] Searching PubMed (anesthesiology)...
[DEBUG] Found 5 papers (anesthesiology)  # <-- Success!
[DEBUG] Fetching 5 papers from PubMed...
```

Instead of:
```bash
[ERROR] PubMed search failed (anesthesiology): HTTP Error 400: Bad Request
```
