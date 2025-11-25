#!/bin/bash

echo "========================================================================"
echo "ORCA Parser - Complete Test Suite Runner"
echo "========================================================================"
echo ""

# Run comprehensive tests
echo "📋 Running comprehensive test suite (36 tests)..."
echo "------------------------------------------------------------------------"
python tests/test_comprehensive.py

# Check exit code
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================================================"
    echo "✓ ALL TESTS PASSED!"
    echo "========================================================================"
    echo ""
    echo "Next steps:"
    echo "1. Start web UI: python app.py"
    echo "2. Open browser: http://localhost:5000"
    echo "3. Click 'Load Data' to see all parsed sections"
    echo ""
else
    echo ""
    echo "========================================================================"
    echo "✗ SOME TESTS FAILED"
    echo "========================================================================"
    echo ""
    echo "Check the output above for details."
    echo ""
    exit 1
fi
