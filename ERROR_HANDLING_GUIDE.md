# Error Handling Standardization Guide - GasConsult.ai

## Current State
- **Inconsistent Patterns:** Mix of flash messages, inline errors, alerts, and no handling
- **User Experience:** Varies across different routes and features
- **No Centralized System:** Error handling scattered throughout codebase

---

## Recommended Standard: Centralized Error Handler

### JavaScript Error Handler Class

```javascript
class ErrorHandler {
    static show(message, type = 'error', duration = 5000) {
        // Create toast notification
        const toast = document.createElement('div');
        toast.className = `toast toast-${type}`;
        toast.innerHTML = `
            <div class="toast-icon">${this.getIcon(type)}</div>
            <div class="toast-message">${message}</div>
            <button class="toast-close" onclick="this.parentElement.remove()">×</button>
        `;
        
        document.body.appendChild(toast);
        
        // Auto-dismiss after duration
        if (duration > 0) {
            setTimeout(() => toast.remove(), duration);
        }
    }
    
    static getIcon(type) {
        const icons = {
            error: '❌',
            warning: '⚠️',
            success: '✅',
            info: 'ℹ️'
        };
        return icons[type] || icons.info;
    }
    
    static async handleFetchError(response) {
        if (!response.ok) {
            const error = await response.json();
            this.show(error.message || 'Request failed', 'error');
            throw new Error(error.message);
        }
        return response;
    }
}
```

### CSS for Toast Notifications

```css
.toast {
    position: fixed;
    top: 20px;
    right: 20px;
    min-width: 300px;
    max-width: 500px;
    padding: 16px 20px;
    background: white;
    border-radius: var(--radius-md);
    box-shadow: 0 4px 12px rgba(0,0,0,0.15);
    display: flex;
    align-items: center;
    gap: 12px;
    z-index: 9999;
    animation: slideIn 0.3s ease;
}

@keyframes slideIn {
    from {
        transform: translateX(400px);
        opacity: 0;
    }
    to {
        transform: translateX(0);
        opacity: 1;
    }
}

.toast-error {
    border-left: 4px solid var(--red-500);
}

.toast-success {
    border-left: 4px solid var(--green-500);
}

.toast-warning {
    border-left: 4px solid var(--yellow-500);
}

.toast-info {
    border-left: 4px solid var(--blue-500);
}

.toast-message {
    flex: 1;
    font-size: 14px;
    color: var(--gray-800);
}

.toast-close {
    background: none;
    border: none;
    font-size: 20px;
    color: var(--gray-400);
    cursor: pointer;
    padding: 0;
    width: 24px;
    height: 24px;
}

.toast-close:hover {
    color: var(--gray-600);
}
```

---

## Usage Patterns

### 1. Fetch Requests (Recommended)

```javascript
async function submitForm(data) {
    try {
        const response = await fetch('/api/endpoint', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(data)
        });
        
        await ErrorHandler.handleFetchError(response);
        const result = await response.json();
        
        ErrorHandler.show('Success!', 'success');
        return result;
        
    } catch (error) {
        ErrorHandler.show(error.message || 'Request failed', 'error');
        throw error;
    }
}
```

### 2. Form Validation

```javascript
function validateForm(formData) {
    const errors = [];
    
    if (!formData.email) {
        errors.push('Email is required');
    }
    
    if (!formData.password || formData.password.length < 8) {
        errors.push('Password must be at least 8 characters');
    }
    
    if (errors.length > 0) {
        ErrorHandler.show(errors.join('<br>'), 'warning');
        return false;
    }
    
    return true;
}
```

### 3. Network Errors

```javascript
window.addEventListener('online', () => {
    ErrorHandler.show('Connection restored', 'success', 3000);
});

window.addEventListener('offline', () => {
    ErrorHandler.show('No internet connection', 'error', 0); // Don't auto-dismiss
});
```

---

## Backend Error Responses (Python/Flask)

### Standard Error Response Format

```python
def standardize_error(error_message, status_code=400):
    """Return standardized error response"""
    return jsonify({
        'success': False,
        'error': {
            'message': error_message,
            'code': status_code
        }
    }), status_code

# Usage in routes
@app.route('/api/example', methods=['POST'])
def example_route():
    try:
        # Process request
        data = request.json
        result = process_data(data)
        
        return jsonify({
            'success': True,
            'data': result
        }), 200
        
    except ValueError as e:
        return standardize_error(str(e), 400)
    except Exception as e:
        logging.error(f"Unexpected error: {str(e)}")
        return standardize_error("An unexpected error occurred", 500)
```

---

## Migration Plan

### Phase 1: Add Error Handler Class
- Add `ErrorHandler` class to shared JavaScript
- Add toast notification CSS to shared styles
- Test with one route

### Phase 2: Update Critical Routes
1. **Chat Interface** (`/stream`, `/chat`)
2. **Authentication** (`/login`, `/register`)
3. **Forms** (`/preop`, `/hypotension`)

### Phase 3: Replace Legacy Patterns
- Remove `alert()` calls → `ErrorHandler.show()`
- Replace inline error divs → Toast notifications
- Standardize Flask error responses

---

## Routes Requiring Updates

### High Priority (User-Facing)
- `/chat` - Chat submission errors
- `/stream` - SSE connection errors
- `/login`, `/register` - Authentication errors
- `/preop`, `/hypotension` - Form validation

### Medium Priority (Tools)
- `/quick-dose` - Calculation errors
- `/calculators` - Input validation
- `/library` - Save/delete errors

### Low Priority (Static)
- `/terms`, `/privacy` - Usually no errors
- `/` (homepage) - Minimal interaction

---

## Testing Checklist

- [ ] Test error toast appears and dismisses
- [ ] Test multiple errors (stacking)
- [ ] Test success notifications
- [ ] Test warning notifications
- [ ] Test network offline/online
- [ ] Test form validation errors
- [ ] Test API error responses
- [ ] Test mobile responsiveness

---

## Benefits

1. **Consistency:** Same error UX across all pages
2. **User-Friendly:** No jarring `alert()` popups
3. **Accessible:** Toast notifications with ARIA labels
4. **Maintainable:** Single error handling system
5. **Debuggable:** Centralized error logging

---

**Status:** Documentation complete. Implementation requires updates to shared JavaScript and systematic route-by-route migration.
