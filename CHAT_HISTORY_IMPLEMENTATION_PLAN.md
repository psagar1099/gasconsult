# Chat History Implementation Plan

## Executive Summary

This document outlines a safe, incremental approach to implementing persistent chat history for gasconsult.ai without breaking the existing chat functionality. The plan prioritizes backwards compatibility and includes comprehensive testing at each step.

---

## Current Architecture Analysis

### How Chat Currently Works

1. **Session-Based Storage (Temporary)**
   - Chat messages stored in `session['messages']` (Redis-backed)
   - Session expires after 1 hour
   - Format: `[{"role": "user", "content": "..."}, {"role": "assistant", "content": "...", "references": [...]}]`

2. **Single-Page Flow**
   - Homepage (`/`) handles both welcome screen AND chat interface
   - When `session['messages']` is empty → shows hero/welcome screen
   - When `session['messages']` has content → shows chat interface
   - Uses `session['chat_active']` flag to track state

3. **Message Processing Flow**
   ```
   User submits query (POST /)
   → Adds to session['messages']
   → Generates request_id
   → Stores stream_data in session
   → Redirects to /stream with request_id
   → Streams GPT response
   → Saves complete response to session['messages']
   → Displays in chat UI
   ```

### Why Previous Attempts Likely Failed

Based on the code analysis, the previous attempt to "transfer chat from home page to dedicated chat page" likely failed because:

1. **Session Data Not Properly Passed**
   - Moving to a separate route without carrying over `session['messages']`
   - The streaming endpoint (`/stream`) depends on data stored in the session during the POST request
   - If session state is lost during redirect/transfer, streaming breaks

2. **State Management Confusion**
   - The app uses `session['chat_active']` to determine UI state
   - Redirecting to a new route might have cleared this flag
   - The homepage route has complex logic to decide between hero view vs chat view

3. **Stream Data Synchronization**
   - The `/stream` endpoint looks for `stream_data_{request_id}` in the session
   - If the session is cleared or not properly shared between requests, streaming fails
   - Error: "Session expired or not found"

4. **Form Submission Conflicts**
   - Both homepage and chat page trying to POST to the same `/` route
   - Confusion about which route should handle the form submission

---

## Recommended Solution: SQLite-Based Persistent Storage

### Why SQLite?

✅ **No Additional Infrastructure** - Single file database, no server needed
✅ **Perfect for Future User Accounts** - Just add `user_id` column later
✅ **Persistent Across Sessions** - Survives server restarts
✅ **Lightweight & Fast** - Suitable for current scale
✅ **Easy to Migrate Later** - Can move to PostgreSQL if needed

### Architecture Overview

```
Current: Session (Temporary) → Lost after 1 hour
New:     Session (Current Chat) + SQLite (Persistent History)

         - Session continues to work exactly as before
         - After each message, also save to SQLite
         - User can browse/load old conversations from database
         - Current conversation still in session for performance
```

### Database Schema

```sql
-- conversations table
CREATE TABLE conversations (
    id TEXT PRIMARY KEY,                    -- UUID
    user_session_id TEXT,                   -- For now: session ID (later: user_id)
    title TEXT,                             -- Auto-generated from first query
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    message_count INTEGER DEFAULT 0,
    is_active BOOLEAN DEFAULT 1             -- Soft delete support
);

-- messages table
CREATE TABLE messages (
    id TEXT PRIMARY KEY,                    -- UUID
    conversation_id TEXT NOT NULL,          -- Foreign key to conversations
    role TEXT NOT NULL,                     -- 'user' or 'assistant'
    content TEXT NOT NULL,                  -- Message content
    references TEXT,                        -- JSON array of paper references
    num_papers INTEGER DEFAULT 0,           -- Number of papers cited
    evidence_strength TEXT,                 -- 'high', 'moderate', 'low'
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (conversation_id) REFERENCES conversations(id)
);

-- CREATE INDEX for performance
CREATE INDEX idx_conversations_session ON conversations(user_session_id, created_at DESC);
CREATE INDEX idx_messages_conversation ON messages(conversation_id, created_at ASC);
```

---

## Implementation Plan (Step-by-Step)

### Phase 1: Database Setup (Non-Breaking)

**Objective:** Add database layer WITHOUT changing any existing functionality

**Steps:**
1. Add SQLite dependency (already in Python stdlib, no new requirement)
2. Create `database.py` module with:
   - `init_db()` - Creates tables if they don't exist
   - `save_conversation()` - Saves conversation to database
   - `save_message()` - Saves individual message
   - `get_conversations()` - Retrieves conversation list
   - `get_conversation()` - Retrieves specific conversation
   - `delete_conversation()` - Soft delete

3. Initialize database on app startup
4. **NO CHANGES TO EXISTING ROUTES YET**

**Testing:**
- [ ] Database file created successfully
- [ ] Tables created with correct schema
- [ ] Can save/retrieve test data
- [ ] Existing chat functionality still works 100%

---

### Phase 2: Silent Data Persistence (Non-Breaking)

**Objective:** Start saving to database in background, no UI changes

**Steps:**
1. Modify `/` route to save messages to database AFTER they're added to session
2. On first message: Create new conversation record
3. On each message: Save to messages table
4. Use `session.get('conversation_id')` to track current conversation
5. **NO UI CHANGES - User doesn't see anything different yet**

**Code Changes:**
```python
# In index() route, after session['messages'].append(...)
if 'conversation_id' not in session:
    # First message - create new conversation
    conversation_id = str(uuid.uuid4())
    session['conversation_id'] = conversation_id
    save_conversation(
        id=conversation_id,
        user_session_id=session.get('session_id', 'anonymous'),
        title=generate_title(raw_query)  # Auto-generate from query
    )

# Save message to database (non-blocking)
save_message(
    conversation_id=session['conversation_id'],
    role='user',
    content=raw_query
)
```

**Testing:**
- [ ] Messages saved to database correctly
- [ ] Conversation records created properly
- [ ] Session still works exactly as before
- [ ] No performance degradation
- [ ] Existing chat UI unchanged

---

### Phase 3: Add History Sidebar (Additive Change)

**Objective:** Add UI to view history WITHOUT changing existing chat flow

**Steps:**
1. Add sidebar toggle button (hamburger menu icon)
2. Create sidebar component in HTML template
3. Add `/api/conversations` endpoint to fetch conversation list
4. Display conversation list with timestamps
5. **Current chat continues to work exactly as before**

**UI Changes:**
```html
<!-- Add to existing template -->
<button id="history-toggle" class="history-btn">
    <svg><!-- hamburger icon --></svg>
</button>

<div id="history-sidebar" class="history-sidebar hidden">
    <h3>Chat History</h3>
    <div id="conversation-list">
        <!-- Populated via AJAX -->
    </div>
</div>
```

**New Routes:**
```python
@app.route("/api/conversations")
def get_user_conversations():
    """Get list of user's past conversations"""
    session_id = session.get('session_id', 'anonymous')
    conversations = get_conversations(user_session_id=session_id, limit=50)
    return jsonify(conversations)
```

**Testing:**
- [ ] Sidebar appears/hides correctly
- [ ] Conversation list loads
- [ ] Current chat still works perfectly
- [ ] No interference between sidebar and main chat

---

### Phase 4: Load Previous Conversations (Carefully!)

**Objective:** Allow loading old conversations WITHOUT breaking current flow

**Steps:**
1. Add click handler to conversation list items
2. Create `/load-conversation/<id>` endpoint
3. **CRITICAL:** Clear session safely before loading
4. Load messages from database into `session['messages']`
5. Set `session['conversation_id']` to loaded conversation
6. Redirect to homepage (will show chat interface)

**Code Changes:**
```python
@app.route("/load-conversation/<conversation_id>")
def load_conversation(conversation_id):
    """Load a previous conversation from database"""
    # Get conversation from database
    conversation = get_conversation(conversation_id)

    if not conversation:
        return jsonify({'error': 'Conversation not found'}), 404

    # Clear current session safely
    session.pop('messages', None)
    session.pop('conversation_topic', None)
    session.pop('stream_data_*', None)  # Clear any stream data

    # Load messages from database
    session['messages'] = conversation['messages']
    session['conversation_id'] = conversation_id
    session['chat_active'] = True
    session.modified = True

    return redirect(url_for('index'))
```

**Testing:**
- [ ] Can load old conversations
- [ ] Loaded conversation displays correctly
- [ ] Can continue old conversation (new messages saved correctly)
- [ ] Can start new conversation after loading old one
- [ ] No corruption of data

---

### Phase 5: Polish & Features

**Objective:** Add nice-to-have features

**Features:**
1. Search conversations (by title/content)
2. Delete conversations
3. Share conversations (extend existing share feature)
4. Export conversations (PDF/text)
5. Conversation titles (auto-generated from first query)

---

## Risk Mitigation Strategies

### 1. Session Integrity Protection

**Risk:** Loading conversations might corrupt session state

**Mitigation:**
- Always use `session.pop()` instead of `del session[...]` (safer)
- Set `session.modified = True` after any changes
- Clear stream data keys explicitly
- Test with Redis session backend (production-like)

### 2. Streaming Synchronization

**Risk:** Stream endpoint can't find data after loading conversation

**Mitigation:**
- Don't allow new queries while old conversation is loading
- Clear all `stream_data_*` keys before loading
- Use request_id validation to ensure fresh stream data

### 3. Database Write Failures

**Risk:** Database errors break chat functionality

**Mitigation:**
- Wrap all database calls in try/except
- Log errors but don't fail the request
- Session-based chat continues working even if database fails

```python
try:
    save_message(conversation_id, role, content)
except Exception as e:
    logger.error(f"Failed to save message to database: {e}")
    # Continue - session-based chat still works
```

### 4. Performance Degradation

**Risk:** Database operations slow down chat

**Mitigation:**
- Use indexes on common queries
- Limit conversation list (e.g., last 50)
- Consider async database writes (future optimization)
- SQLite is fast for read/write operations at this scale

### 5. Session ID Persistence

**Risk:** Users lose history when session expires

**Mitigation:**
- Generate a persistent `user_tracking_id` (stored in cookie)
- Use this for conversation association instead of session ID
- Cookie lasts longer than session (e.g., 30 days)
- Later: Replace with user_id when accounts are added

---

## Testing Checklist

### Before Implementation
- [ ] Current chat works perfectly
- [ ] Streaming works correctly
- [ ] Session persistence works
- [ ] Clear/reset works

### After Phase 1 (Database Setup)
- [ ] Database initialized successfully
- [ ] Can write/read test data
- [ ] **All existing functionality still works**

### After Phase 2 (Silent Persistence)
- [ ] Messages saved to database
- [ ] Conversations created correctly
- [ ] **Chat UI unchanged and working**
- [ ] No performance impact

### After Phase 3 (History Sidebar)
- [ ] Sidebar appears/hides
- [ ] Conversation list loads
- [ ] **Current chat unaffected**
- [ ] Can continue chatting with sidebar open

### After Phase 4 (Load Conversations)
- [ ] Can load old conversations
- [ ] Loaded conversation displays correctly
- [ ] Can continue old conversations
- [ ] Can start new conversations
- [ ] Switching between conversations works
- [ ] **No data corruption**

### Edge Cases
- [ ] Load conversation while current chat active
- [ ] Clear chat then load old conversation
- [ ] Start new chat after loading old conversation
- [ ] Session expires while viewing history
- [ ] Database connection fails (chat still works)
- [ ] Multiple tabs/windows open

---

## User Account Preparation

When you're ready to add user accounts, the migration is simple:

### Database Changes:
```sql
-- Add user_id column
ALTER TABLE conversations ADD COLUMN user_id TEXT;

-- Migrate old data
UPDATE conversations
SET user_id = user_session_id
WHERE user_id IS NULL;

-- Update index
CREATE INDEX idx_conversations_user ON conversations(user_id, created_at DESC);
```

### Code Changes:
- Replace `session.get('session_id')` with `current_user.id`
- Add user authentication (Flask-Login)
- Update `/api/conversations` to filter by `user_id`

---

## Rollback Plan

If anything goes wrong at any phase:

1. **Immediate Rollback:**
   ```bash
   git revert <commit-hash>
   git push origin claude/add-chat-history-*
   ```

2. **Database Isolation:**
   - Database file: `gasconsult.db` (separate from app.py)
   - Can delete database file without affecting app
   - App works fine without database (session-based only)

3. **Feature Flags:**
   - Add `ENABLE_CHAT_HISTORY` environment variable
   - Can disable feature without code changes

---

## Timeline Estimate

| Phase | Description | Development Time | Testing Time |
|-------|-------------|------------------|--------------|
| 1 | Database setup | 2 hours | 1 hour |
| 2 | Silent persistence | 2 hours | 2 hours |
| 3 | History sidebar | 3 hours | 2 hours |
| 4 | Load conversations | 2 hours | 3 hours |
| 5 | Polish & features | 4 hours | 2 hours |
| **Total** | | **13 hours** | **10 hours** |

**Recommended Approach:** Implement one phase at a time, deploy to staging, test thoroughly, then proceed to next phase.

---

## Questions to Answer Before Starting

1. **Where should the database file live?**
   - Option A: `/home/user/gasconsult/gasconsult.db` (same directory as app.py)
   - Option B: `/var/lib/gasconsult/gasconsult.db` (dedicated data directory)
   - **Recommendation:** Option A for simplicity, Option B for production

2. **How long should conversations be kept?**
   - Option A: Forever (until user deletes)
   - Option B: Auto-delete after 90 days
   - **Recommendation:** Keep forever, add manual delete option

3. **Should we implement feature flag?**
   - Pros: Can disable quickly if issues arise
   - Cons: Extra complexity
   - **Recommendation:** Yes, add `ENABLE_CHAT_HISTORY` env var

4. **Do you want to implement all phases, or start with Phase 1-2?**
   - **Recommendation:** Start with Phase 1-2 (silent persistence), verify it works, then add UI

---

## Next Steps

Once you approve this plan, I will:

1. Create the database module (`database.py`)
2. Implement Phase 1 (Database Setup)
3. Test thoroughly
4. Commit and push
5. Wait for your approval before proceeding to Phase 2

**Please review this plan and let me know:**
- Do you approve this approach?
- Any concerns or questions?
- Should I proceed with Phase 1?
