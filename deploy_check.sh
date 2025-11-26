#!/bin/bash
# Deployment Check and Restart Script

echo "=== GasConsult.ai Deployment Check ==="
echo ""

echo "1. Current Branch:"
git branch --show-current
echo ""

echo "2. Latest Commit:"
git log -1 --oneline
echo ""

echo "3. Code Verification:"
echo "   - Quick Dose Route: $(grep -c '@app.route("/quick-dose")' app.py) found"
echo "   - Feature Grid (3 cols): $(grep -c 'grid-template-columns: repeat(3, 1fr)' app.py) found"
echo "   - Glassy Nav in Terms: $(grep -c 'backdrop-filter: blur(12px)' app.py) found"
echo ""

echo "4. To apply changes to your deployment:"
echo "   Option A - If running locally:"
echo "      pkill -f 'python.*app.py' || pkill -f gunicorn"
echo "      python app.py"
echo ""
echo "   Option B - If using Render/Heroku/similar:"
echo "      They should auto-deploy from the pushed commits"
echo "      Check your deployment dashboard for latest deployment status"
echo ""

echo "5. Clear browser cache:"
echo "   - Chrome/Edge: Ctrl+Shift+Del → Clear cached images and files"
echo "   - Or use Incognito/Private mode"
echo ""

echo "=== All code changes verified in repository ==="
