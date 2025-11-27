# GasConsult.ai - Public Release Deployment Guide

## 🚀 Pre-Release Checklist

This guide covers everything you need to do before releasing gasconsult.ai to the public.

---

## ✅ Completed Security & Infrastructure Improvements

### 1. Security Hardening ✓
- [x] Moved hardcoded API credentials to environment variables
- [x] Added comprehensive input sanitization (XSS protection)
- [x] Implemented rate limiting (60 requests/minute default)
- [x] Added CSRF protection for all forms
- [x] Server-side session management

### 2. Monitoring & Logging ✓
- [x] Comprehensive logging system with configurable levels
- [x] Request/response logging middleware
- [x] Error tracking and exception handling
- [x] Health check endpoint (`/health`)
- [x] API status endpoint (`/api/status`)

### 3. Legal & Compliance ✓
- [x] Privacy Policy page (`/privacy`)
- [x] Existing Terms of Service page (`/terms`)
- [x] Medical disclaimer visible on site

### 4. Configuration Management ✓
- [x] `.env.example` file with all required variables
- [x] `python-dotenv` for environment variable loading
- [x] Proper `.gitignore` for sensitive files

---

## 📋 Remaining Tasks Before Public Release

### High Priority (Must Complete)

#### 1. Environment Configuration
**Action Required:** Create `.env` file with production values

```bash
# Copy the example file
cp .env.example .env

# Edit with your production values
nano .env
```

**Required Variables:**
- `OPENAI_API_KEY` - Get from https://platform.openai.com/api-keys
- `ENTREZ_EMAIL` - Your real email for NCBI
- `ENTREZ_API_KEY` - Get from https://www.ncbi.nlm.nih.gov/account/settings/
- `FLASK_SECRET_KEY` - Generate with: `python -c "import secrets; print(secrets.token_hex(32))"`

#### 2. Install Updated Dependencies

```bash
pip install -r requirements.txt
```

**New dependencies added:**
- `bleach==6.1.0` - HTML sanitization
- `Flask-Limiter==3.5.0` - Rate limiting
- `Flask-WTF==1.2.1` - CSRF protection
- `python-dotenv==1.0.0` - Environment variable loading

#### 3. Testing

**Manual Testing Checklist:**
- [ ] Test `/health` endpoint returns 200
- [ ] Test rate limiting (make 61 requests in 1 minute)
- [ ] Test XSS protection (try entering `<script>alert('xss')</script>`)
- [ ] Test all forms (chat, preop, hypotension)
- [ ] Test Privacy Policy page loads (`/privacy`)
- [ ] Test Terms of Service page loads (`/terms`)
- [ ] Verify no errors in logs

**Run the application:**
```bash
# Development
python app.py

# Production (recommended)
gunicorn app:app -b 0.0.0.0:8000 --workers 4 --timeout 120
```

#### 4. Deployment Configuration

**For Render/Heroku/Cloud Platforms:**

1. **Set Environment Variables** in your hosting dashboard:
   - `OPENAI_API_KEY`
   - `ENTREZ_EMAIL`
   - `ENTREZ_API_KEY`
   - `FLASK_SECRET_KEY`
   - `FLASK_ENV=production`
   - `LOG_LEVEL=INFO`
   - `RATE_LIMIT=60 per minute`

2. **Configure Build Settings:**
   - Build Command: `pip install -r requirements.txt`
   - Start Command: `gunicorn app:app -b 0.0.0.0:$PORT --workers 4`

3. **Enable HTTPS** (required for production)

4. **Set up monitoring** using the `/health` endpoint

#### 5. Domain & DNS
- [ ] Configure custom domain
- [ ] Set up SSL/TLS certificate
- [ ] Configure DNS records
- [ ] Test domain resolution

#### 6. Error Monitoring (Recommended)

**Option A: Sentry (Free tier available)**
```bash
pip install sentry-sdk[flask]
```

Add to `app.py`:
```python
import sentry_sdk
from sentry_sdk.integrations.flask import FlaskIntegration

sentry_sdk.init(
    dsn=os.getenv("SENTRY_DSN"),
    integrations=[FlaskIntegration()],
    environment=os.getenv("FLASK_ENV", "production")
)
```

Update `.env`:
```
SENTRY_DSN=your-sentry-dsn-here
```

### Medium Priority (Recommended)

#### 7. Performance Optimization
- [ ] Enable gzip compression
- [ ] Add CDN for static assets
- [ ] Implement PubMed query caching (Redis)
- [ ] Add database for analytics (optional)

#### 8. User Analytics (Privacy-Respecting)
- [ ] Add privacy-respecting analytics (Plausible, Fathom, or Simple Analytics)
- [ ] Track page views and feature usage
- [ ] Monitor error rates

#### 9. Backup & Disaster Recovery
- [ ] Set up automated backups
- [ ] Document recovery procedures
- [ ] Test restore process

### Low Priority (Nice to Have)

#### 10. Documentation
- [ ] API documentation
- [ ] User guide
- [ ] FAQ page
- [ ] Video tutorials

#### 11. Additional Features
- [ ] User accounts (optional)
- [ ] Saved queries
- [ ] Export chat history
- [ ] Mobile app

---

## 🔒 Security Best Practices

### Before Going Live:

1. **Never commit `.env` to git** (already in `.gitignore`)
2. **Use strong secret keys** (minimum 32 characters)
3. **Enable HTTPS only** (redirect HTTP to HTTPS)
4. **Keep dependencies updated** (`pip list --outdated`)
5. **Monitor security advisories** for Python packages
6. **Review logs regularly** for suspicious activity
7. **Set up alerts** for error spikes

### Ongoing Maintenance:

```bash
# Check for security vulnerabilities
pip install safety
safety check

# Update dependencies
pip install --upgrade pip
pip install -r requirements.txt --upgrade
```

---

## 📊 Monitoring Endpoints

### Health Check
```bash
curl https://yourdomain.com/health
```

**Expected Response (200 OK):**
```json
{
  "status": "healthy",
  "service": "gasconsult.ai",
  "version": "1.0.0",
  "python_version": "3.x.x",
  "checks": {
    "openai": true,
    "entrez_email": true,
    "entrez_api_key": true,
    "secret_key": true
  }
}
```

**Degraded Response (503):**
```json
{
  "status": "degraded",
  ...
}
```

### API Status
```bash
curl https://yourdomain.com/api/status
```

---

## 🚨 Troubleshooting

### Application Won't Start

1. **Check environment variables:**
   ```bash
   python -c "from dotenv import load_dotenv; import os; load_dotenv(); print('OPENAI_API_KEY:', bool(os.getenv('OPENAI_API_KEY')))"
   ```

2. **Check logs:**
   ```bash
   tail -f app.log  # if LOG_FILE is set
   ```

3. **Test imports:**
   ```bash
   python -c "import flask, bleach, flask_limiter, flask_wtf; print('All imports OK')"
   ```

### Rate Limiting Too Strict

Update `.env`:
```
RATE_LIMIT=120 per minute
```

Or in code, modify `app.py` line ~34.

### CSRF Errors

CSRF protection is now enabled. To exempt specific routes:
```python
@app.route("/your-route")
@csrf.exempt
def your_function():
    ...
```

---

## 📝 Legal Considerations

### Before Public Release:

1. **Review Privacy Policy** (`/privacy`)
   - Update contact email if needed
   - Ensure compliance with local laws

2. **Review Terms of Service** (`/terms`)
   - Consider adding limitation of liability
   - Clarify educational use only

3. **Medical Disclaimer**
   - Prominently display on all clinical pages
   - Consider requiring acknowledgment

4. **HIPAA Compliance**
   - You are NOT a covered entity
   - Do NOT accept PHI
   - Clearly state in disclaimers

5. **Professional Liability**
   - Consider consulting with legal counsel
   - Review malpractice insurance needs
   - Document educational purpose

---

## 🎯 Launch Checklist

### Day Before Launch:
- [ ] All environment variables set in production
- [ ] Dependencies installed and tested
- [ ] Health check endpoint returning 200
- [ ] HTTPS enabled and working
- [ ] Domain configured
- [ ] Error monitoring active (Sentry or similar)
- [ ] Backups configured
- [ ] Team has access to admin dashboard
- [ ] Emergency contact list prepared

### Launch Day:
- [ ] Final smoke test on production
- [ ] Monitor logs for first hour
- [ ] Test all critical paths
- [ ] Announce to target audience
- [ ] Monitor /health endpoint
- [ ] Watch for error spikes

### Post-Launch (First Week):
- [ ] Daily log review
- [ ] Monitor performance metrics
- [ ] Gather user feedback
- [ ] Address any bugs/issues
- [ ] Update documentation based on feedback

---

## 📞 Support & Maintenance

### Regular Maintenance Schedule:

**Daily:**
- Check `/health` endpoint
- Review error logs

**Weekly:**
- Review usage analytics
- Update dependencies if needed
- Check for security advisories

**Monthly:**
- Security audit
- Performance review
- User feedback analysis
- Backup testing

---

## 🔗 Useful Resources

- [Flask Security Best Practices](https://flask.palletsprojects.com/en/latest/security/)
- [OWASP Top 10](https://owasp.org/www-project-top-ten/)
- [HIPAA Compliance Guide](https://www.hhs.gov/hipaa/index.html)
- [GDPR Compliance](https://gdpr.eu/)
- [Python Security Best Practices](https://python.readthedocs.io/en/stable/library/security_warnings.html)

---

## 📧 Contact

For questions about deployment or security:
- Email: admin@gasconsult.ai
- GitHub: [Create an issue](https://github.com/yourusername/gasconsult/issues)

---

**Last Updated:** November 27, 2025
**Version:** 1.0.0
