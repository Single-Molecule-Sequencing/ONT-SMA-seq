#!/bin/bash
# Run the ONT-SMA-seq pipeline end-to-end test

set -e  # Exit on error

echo "Running ONT-SMA-seq Pipeline End-to-End Test..."
echo ""

# Check if Python 3 is available
if ! command -v python3 &> /dev/null; then
    echo "Error: python3 is required but not found"
    exit 1
fi

# Check if dependencies are installed
if ! python3 -c "import pysam, pod5, edlib" 2>/dev/null; then
    echo "Installing dependencies..."
    pip install -q -r requirements.txt
fi

# Run the test
python3 test_e2e.py
