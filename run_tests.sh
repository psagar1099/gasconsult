#!/bin/bash
#
# GasConsult.ai Test Suite Runner
# Runs automated tests to prevent template refactoring bugs
#

set -e  # Exit on first error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  GasConsult.ai Automated Test Suite${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Check if pytest is installed
if ! command -v pytest &> /dev/null; then
    echo -e "${RED}ERROR: pytest not installed${NC}"
    echo "Install with: pip install pytest"
    exit 1
fi

# Parse command line arguments
TEST_MODE="${1:-all}"

case "$TEST_MODE" in
    quick)
        echo -e "${YELLOW}Running QUICK smoke tests only...${NC}"
        pytest tests/test_routes.py::TestHomePage -v
        ;;

    critical)
        echo -e "${YELLOW}Running CRITICAL regression tests...${NC}"
        pytest tests/test_routes.py::TestHomePage::test_home_is_chat_interface_not_preop \
               tests/test_routes.py::TestToolRoutes::test_preop_is_not_on_homepage \
               tests/test_templates.py::TestChatStandaloneTemplate::test_chat_standalone_is_home_page \
               tests/test_snapshots.py::TestCrossPageValidation::test_home_is_not_preop \
               -v
        ;;

    templates)
        echo -e "${YELLOW}Running TEMPLATE validation tests...${NC}"
        pytest tests/test_templates.py -v
        ;;

    snapshots)
        echo -e "${YELLOW}Running HTML SNAPSHOT tests...${NC}"
        pytest tests/test_snapshots.py -v
        ;;

    routes)
        echo -e "${YELLOW}Running ROUTE smoke tests...${NC}"
        pytest tests/test_routes.py -v
        ;;

    all)
        echo -e "${YELLOW}Running FULL test suite...${NC}"
        pytest tests/ -v
        ;;

    *)
        echo -e "${RED}Unknown test mode: $TEST_MODE${NC}"
        echo ""
        echo "Usage: ./run_tests.sh [MODE]"
        echo ""
        echo "Available modes:"
        echo "  quick      - Quick smoke tests (home page only)"
        echo "  critical   - Critical regression tests (prevents template swap bugs)"
        echo "  templates  - Template validation tests"
        echo "  snapshots  - HTML snapshot tests"
        echo "  routes     - All route smoke tests"
        echo "  all        - Full test suite (default)"
        echo ""
        exit 1
        ;;
esac

# Check exit code
if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}  ✓ ALL TESTS PASSED${NC}"
    echo -e "${GREEN}========================================${NC}"
    exit 0
else
    echo ""
    echo -e "${RED}========================================${NC}"
    echo -e "${RED}  ✗ TESTS FAILED${NC}"
    echo -e "${RED}========================================${NC}"
    exit 1
fi
