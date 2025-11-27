# Changelog

All notable changes to gasconsult.ai will be documented in this file.

## [1.0.0] - 2025-11-27 - Public Release Preparation

### 🔒 Security Enhancements

#### Credentials & Secrets Management
- **Moved hardcoded API credentials to environment variables**
  - Entrez email and API key now loaded from `.env`
  - OpenAI API key properly configured
  - Flask secret key configurable via environment
  - Created comprehensive `.env.example` template

#### Input Validation & Sanitization
- **Added HTML sanitization using bleach library**
  - Created `sanitize_input()` function for safe HTML rendering
  - Created `sanitize_user_query()` for plain text input
  - Applied sanitization to all user input endpoints:
    - Chat queries (`/chat`)
    - Pre-op assessment form (`/preop`)
    - Hypotension predictor form (`/hypotension`)
  - Prevents XSS attacks while allowing safe GPT-generated HTML

#### Rate Limiting
- **Implemented Flask-Limiter for API protection**
  - Default: 60 requests per minute per IP
  - Configurable via `RATE_LIMIT` environment variable
  - Prevents abuse and DoS attacks
  - Uses in-memory storage (consider Redis for production scale)

#### CSRF Protection
- **Added Flask-WTF CSRF protection**
  - Protects all POST requests by default
  - Health check and status endpoints exempted
  - Automatic token validation

### 📊 Monitoring & Logging

#### Comprehensive Logging System
- **Implemented structured logging**
  - Configurable log levels (DEBUG, INFO, WARNING, ERROR, CRITICAL)
  - Optional file-based logging with rotation (10MB max, 5 backups)
  - Console output for development
  - Request/response logging middleware
  - Automatic exception logging
  - Startup configuration logging

#### Health & Status Endpoints
- **Created `/health` endpoint**
  - Returns JSON status with version info
  - Checks all critical services (OpenAI, Entrez, secrets)
  - Returns 200 (healthy) or 503 (degraded)
  - Perfect for deployment monitoring

- **Created `/api/status` endpoint**
  - Lists all available endpoints
  - Shows rate limit configuration
  - Useful for API discovery

### 📄 Legal & Compliance

#### Privacy Policy
- **Created comprehensive Privacy Policy (`/privacy`)**
  - GDPR/CCPA compliant
  - Clear data collection disclosure
  - User rights explanation
  - PHI/HIPAA guidance
  - Third-party service disclosure (OpenAI, PubMed)
  - Data retention policies
  - International user notice

#### Existing Legal Pages
- Terms of Service (`/terms`) - already present
- Medical disclaimers on clinical pages

### 🛠️ Infrastructure & Configuration

#### Environment Management
- **Added python-dotenv for configuration**
  - Automatic `.env` file loading
  - Clear separation of config from code
  - Documented all variables in `.env.example`

#### Dependencies Added
```
bleach==6.1.0           # HTML sanitization
Flask-Limiter==3.5.0    # Rate limiting
Flask-WTF==1.2.1        # CSRF protection
python-dotenv==1.0.0    # Environment variables
```

### 📚 Documentation

- **Created DEPLOYMENT_GUIDE.md**
  - Complete pre-release checklist
  - Environment setup instructions
  - Testing procedures
  - Deployment configuration
  - Security best practices
  - Monitoring guidelines
  - Troubleshooting guide

- **Created CHANGELOG.md** (this file)
  - Documents all changes for public release

- **Updated .env.example**
  - All required environment variables
  - Detailed comments and instructions
  - Default values where appropriate

### 🔧 Code Quality Improvements

#### Error Handling
- Global exception handler with logging
- Graceful error responses
- No stack traces exposed to users

#### Code Organization
- Added security functions section
- Added logging configuration section
- Improved code documentation
- Clear separation of concerns

### ⚠️ Breaking Changes

#### Environment Variables Required
Applications **must** now set these environment variables:
- `OPENAI_API_KEY` (required)
- `ENTREZ_EMAIL` (required)
- `ENTREZ_API_KEY` (recommended)
- `FLASK_SECRET_KEY` (recommended for production)

**Migration Steps:**
1. Copy `.env.example` to `.env`
2. Fill in your actual API keys
3. Restart the application

#### CSRF Protection
All POST requests now require CSRF tokens. If you're using the application via API or scripts:
- Either disable CSRF for specific routes using `@csrf.exempt`
- Or include CSRF tokens in requests

### 📝 Notes for Future Releases

#### High Priority (Not Yet Implemented)
- [ ] Automated testing suite (unit + integration tests)
- [ ] CI/CD pipeline
- [ ] Redis-backed rate limiting for multi-server deployments
- [ ] PubMed query caching layer
- [ ] Enhanced medical disclaimers on all clinical pages
- [ ] Cookie consent banner (if using analytics)

#### Medium Priority
- [ ] User accounts and saved queries
- [ ] Export functionality (PDF/JSON)
- [ ] Advanced analytics dashboard
- [ ] Mobile responsiveness improvements
- [ ] Progressive Web App (PWA) enhancements

#### Low Priority
- [ ] Internationalization (i18n)
- [ ] Dark mode toggle
- [ ] Accessibility audit (WCAG 2.1)
- [ ] Performance optimization (CDN, caching)

---

## Version History

### [0.x] - Pre-Release Development
- Initial implementation
- Core features: PubMed search, GPT-4 integration
- Pre-op assessment tool
- Hypotension predictor
- Quick dose calculator
- Basic UI/UX

---

**Versioning Scheme:** We use Semantic Versioning (SemVer)
- MAJOR version for incompatible API changes
- MINOR version for backwards-compatible functionality additions
- PATCH version for backwards-compatible bug fixes

**Release Dates:** All dates in YYYY-MM-DD format
