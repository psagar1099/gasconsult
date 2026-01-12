# Phase 1: UI/UX Enhancements - Implementation Plan

## Overview
Implementing 4 quick-win features to outperform UpToDate AI and NYSORA AI:
1. Reasoning trace visibility
2. Predictive question suggestions
3. Enhanced evidence badges
4. Conversational context visualization

## Feature 1: Reasoning Trace Visibility

### Current State
- PreopAgent already has `self.reasoning_trace` list tracking 7 steps
- Not currently displayed in UI
- Data exists but hidden from users

### Implementation

#### Backend Changes (app.py)
1. **Extend reasoning trace to regular chat:**
   - Add `reasoning_steps` list to session messages
   - Track steps like: "Synonym expansion", "PubMed search (N papers)", "GPT synthesis"

2. **Pass reasoning trace in responses:**
   ```python
   # In index() route after processing
   session['messages'].append({
       "role": "assistant",
       "content": cleaned_response,
       "references": refs,
       "num_papers": num_papers,
       "reasoning_trace": reasoning_steps  # NEW
   })
   ```

3. **PreopAgent already has it:**
   ```python
   # Line 6829: self.reasoning_trace = []
   # Line 6836: self.reasoning_trace.append("✓ Step 1: Risk Stratification")
   # Just need to pass it to template
   ```

#### Frontend Changes (HTML template)
1. **Add collapsible reasoning section:**
   ```html
   {% if message.get('reasoning_trace') %}
   <div class="reasoning-section">
       <button class="reasoning-toggle" onclick="toggleReasoning(this)">
           <svg>...</svg> Show AI reasoning process
       </button>
       <div class="reasoning-content" style="display: none;">
           <div class="reasoning-timeline">
               {% for step in message.reasoning_trace %}
               <div class="reasoning-step">
                   <div class="step-indicator">{{ loop.index }}</div>
                   <div class="step-text">{{ step }}</div>
               </div>
               {% endfor %}
           </div>
       </div>
   </div>
   {% endif %}
   ```

2. **CSS styling:**
   - Timeline with connected dots
   - Smooth expand/collapse animation
   - Color-coded steps (search blue, analysis green, synthesis purple)

### Success Metrics
- Users can see exactly how the AI reached its conclusion
- Transparency builds trust (UpToDate doesn't show this level of detail)

---

## Feature 2: Predictive Question Suggestions

### Current State
- `.followup-container` already exists in HTML (line 10495)
- Empty placeholder waiting for content
- `data-message-index` attribute for tracking

### Implementation

#### Backend Changes
1. **Generate suggestions after each response:**
   ```python
   def generate_followup_questions(query, response_content, conversation_history):
       """Use GPT-4o to generate 3-5 smart follow-up questions"""
       prompt = f"""Based on this conversation:

       User asked: {query}
       Assistant answered: {response_content[:500]}...

       Generate 3-5 smart follow-up questions that:
       1. Address natural next questions (dosing → side effects → contraindications)
       2. Explore related topics not covered
       3. Consider clinical decision pathways

       Return as JSON array: ["Question 1?", "Question 2?", ...]
       """

       response = openai_client.chat.completions.create(
           model="gpt-4o",
           messages=[{"role": "user", "content": prompt}],
           temperature=0.3,
           response_format={"type": "json_object"}
       )

       return json.loads(response.choices[0].message.content)
   ```

2. **Store in message:**
   ```python
   suggested_questions = generate_followup_questions(raw_query, cleaned_response, session['messages'])

   session['messages'].append({
       ...
       "suggested_questions": suggested_questions  # NEW
   })
   ```

#### Frontend Changes
1. **Populate followup-container:**
   ```html
   {% if message.get('suggested_questions') %}
   <div class="followup-suggestions">
       <div class="followup-title">You might also want to know:</div>
       <div class="followup-grid">
           {% for question in message.suggested_questions %}
           <button class="followup-chip" onclick="askFollowup(this)">
               {{ question }}
           </button>
           {% endfor %}
       </div>
   </div>
   {% endif %}
   ```

2. **JavaScript for click-to-ask:**
   ```javascript
   function askFollowup(button) {
       const question = button.textContent.trim();
       document.getElementById('chatInput').value = question;
       document.getElementById('chatInput').focus();
       // Auto-scroll to input
       button.classList.add('selected');
   }
   ```

3. **CSS styling:**
   - Chip/pill design (rounded, clickable)
   - Hover effects
   - Selected state
   - Grid layout for 3-5 items

### Success Metrics
- Guides users through complex decision-making
- Reduces cognitive load
- No competitor has this feature

---

## Feature 3: Enhanced Evidence Badges

### Current State
- Basic evidence badge with High/Moderate/Low (lines 10437-10458)
- Shows confidence level + number of studies
- Static, non-interactive

### Implementation

#### Backend Changes
1. **Calculate detailed breakdown:**
   ```python
   def calculate_evidence_breakdown(references):
       """Generate detailed evidence analysis"""
       breakdown = {
           'paper_types': {'guideline': 0, 'meta-analysis': 0, 'systematic_review': 0, 'rct': 0, 'cohort': 0},
           'years': {'2024-2025': 0, '2022-2023': 0, '2020-2021': 0, 'older': 0},
           'consensus': 0.0,  # 0-1 score
           'recent_count': 0,  # Papers from last 2 years
           'high_quality_count': 0  # Guidelines + meta-analyses + systematic reviews
       }

       for ref in references:
           # Extract publication type and year
           pub_type = extract_publication_type(ref)
           year = extract_year(ref)

           breakdown['paper_types'][pub_type] += 1

           if year >= 2024:
               breakdown['years']['2024-2025'] += 1
               breakdown['recent_count'] += 1
           # ... etc

       # Calculate consensus (if most papers agree)
       breakdown['consensus'] = calculate_consensus_score(references)

       return breakdown
   ```

2. **Pass in evidence_strength:**
   ```python
   evidence_strength = {
       'level': 'High',  # or Moderate/Low
       'description': '...',
       'breakdown': breakdown,  # NEW
       'total_papers': num_papers
   }
   ```

#### Frontend Changes
1. **Make badge clickable:**
   ```html
   <div class="evidence-badge {{ level }}" onclick="toggleEvidenceDetails(this)">
       <div class="badge-summary">
           ✓ High Confidence • {{ message.num_papers }} studies
       </div>
       <div class="evidence-details-popup" style="display: none;">
           <!-- Detailed breakdown -->
           <div class="breakdown-section">
               <h4>Study Types</h4>
               <div class="breakdown-bar">
                   {% for type, count in breakdown.paper_types.items() if count > 0 %}
                   <div class="bar-segment" style="width: {{ (count / total) * 100 }}%">
                       {{ count }} {{ type }}
                   </div>
                   {% endfor %}
               </div>
           </div>

           <div class="breakdown-section">
               <h4>Publication Years</h4>
               <div class="year-timeline">
                   <!-- Visual timeline -->
               </div>
           </div>

           <div class="breakdown-section">
               <h4>Consensus</h4>
               <div class="consensus-meter">
                   <div class="consensus-fill" style="width: {{ consensus * 100 }}%"></div>
                   <span>{{ (consensus * 100)|int }}% agreement</span>
               </div>
           </div>
       </div>
   </div>
   ```

2. **CSS for popup:**
   - Modal or dropdown design
   - Animated expand
   - Visual charts (bar charts, timeline, meter)

### Success Metrics
- Users understand evidence quality at a glance
- Interactive exploration of evidence
- More transparent than UpToDate's static badges

---

## Feature 4: Conversational Context Visualization

### Current State
- No conversation overview
- Users lose track in long conversations
- No way to jump to previous Q&A pairs

### Implementation

#### Frontend Only (No Backend Changes)
1. **Add sidebar:**
   ```html
   <aside class="conversation-sidebar" id="conversationSidebar">
       <div class="sidebar-header">
           <h3>Conversation Map</h3>
           <button onclick="toggleSidebar()">×</button>
       </div>
       <div class="conversation-tree" id="conversationTree">
           <!-- Populated by JavaScript -->
       </div>
   </aside>

   <button class="sidebar-toggle" onclick="toggleSidebar()">
       <svg>...</svg> Show conversation map
   </button>
   ```

2. **JavaScript to build tree:**
   ```javascript
   function buildConversationTree() {
       const messages = document.querySelectorAll('.message');
       const tree = document.getElementById('conversationTree');
       tree.innerHTML = '';

       let currentTopic = null;

       messages.forEach((msg, index) => {
           if (msg.classList.contains('user-message')) {
               const text = msg.querySelector('.message-bubble').textContent;
               const keywords = extractKeywords(text);

               // Create tree node
               const node = document.createElement('div');
               node.className = 'tree-node';
               node.innerHTML = `
                   <div class="node-icon">🎯</div>
                   <div class="node-text">${text.substring(0, 40)}...</div>
               `;
               node.onclick = () => scrollToMessage(index);
               tree.appendChild(node);
           }
       });
   }

   function scrollToMessage(index) {
       const msg = document.querySelectorAll('.message')[index];
       msg.scrollIntoView({ behavior: 'smooth', block: 'center' });
       msg.classList.add('highlight');
       setTimeout(() => msg.classList.remove('highlight'), 2000);
   }
   ```

3. **CSS for sidebar:**
   - Slide-in animation from right
   - Tree structure with indentation
   - Highlight current message
   - Responsive (hidden on mobile by default)

### Success Metrics
- Easy navigation in long conversations
- Visual overview of conversation flow
- No competitor has this feature

---

## Implementation Order

1. **Reasoning Trace** (2-3 hours)
   - Simple backend change
   - Frontend template + CSS
   - Test with preop tool first

2. **Enhanced Evidence Badges** (3-4 hours)
   - Backend breakdown calculation
   - Frontend interactive popup
   - Visual charts

3. **Predictive Questions** (4-5 hours)
   - Backend GPT-4o integration
   - Frontend chip design
   - JavaScript interaction

4. **Conversation Sidebar** (3-4 hours)
   - Frontend only
   - JavaScript tree building
   - Responsive design

**Total estimated time: 12-16 hours**

---

## Testing Plan

1. **Reasoning Trace:**
   - Test with preop assessment (already has traces)
   - Test with regular chat (need to add traces)
   - Verify all steps display correctly

2. **Evidence Badges:**
   - Test with various paper counts (0, 1-2, 3-5, 6+)
   - Test with different paper types
   - Verify popup works on mobile

3. **Predictive Questions:**
   - Test question generation quality
   - Test click-to-ask functionality
   - Test with conversation context

4. **Sidebar:**
   - Test with short conversations (2-3 Q&A)
   - Test with long conversations (10+ Q&A)
   - Test responsive behavior
   - Test scroll-to-message

---

## Rollout Strategy

1. **Deploy to staging**
2. **A/B test with 20% of users**
3. **Collect feedback**
4. **Iterate based on feedback**
5. **Full rollout**

---

## Success Criteria

- [ ] Reasoning trace visible for all responses (preop + chat)
- [ ] Follow-up questions generated for 100% of responses
- [ ] Evidence badges interactive with detailed breakdown
- [ ] Conversation sidebar functional and responsive
- [ ] All features work on mobile
- [ ] Page load time < 2 seconds
- [ ] No JavaScript errors in console
- [ ] User feedback positive (>80% approval)
