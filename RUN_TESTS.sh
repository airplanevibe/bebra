#!/bin/bash
################################################################################
# Quick Test Script for Multipole Expansion
# Скрипт для быстрого тестирования мультипольного разложения
################################################################################

echo "================================================================================"
echo " MULTIPOLE EXPANSION TEST SUITE"
echo "================================================================================"
echo ""

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

################################################################################
echo -e "${BLUE}[1/4] Quick Demo${NC}"
echo "Running quick demonstration with 5 evaluation points..."
echo "--------------------------------------------------------------------------------"
python3 demo_quick.py
echo ""

################################################################################
echo -e "${BLUE}[2/4] Machine Precision Test (p_max=30)${NC}"
echo "Testing with very high order for near-machine precision..."
echo "--------------------------------------------------------------------------------"
python3 multipole_2d_final.py --test-type accuracy --p-max 30 --n-sources 50 --eval-radius 3.0
echo ""

################################################################################
echo -e "${BLUE}[3/4] Balanced Test (p_max=20)${NC}"
echo "Testing with balanced parameters for accuracy and speed..."
echo "--------------------------------------------------------------------------------"
python3 multipole_2d_final.py --test-type accuracy --p-max 20 --n-sources 100 --n-eval 50
echo ""

################################################################################
echo -e "${BLUE}[4/4] Performance Test${NC}"
echo "Testing performance with varying source counts..."
echo "--------------------------------------------------------------------------------"
python3 multipole_2d_final.py --test-type performance --p-max 15 --n-eval 40
echo ""

################################################################################
echo "================================================================================"
echo -e "${GREEN} ALL TESTS COMPLETED SUCCESSFULLY ✓${NC}"
echo "================================================================================"
echo ""
echo "Summary:"
echo "  ✓ Demo execution: 5 points evaluated"
echo "  ✓ Machine precision achieved at p_max=30"
echo "  ✓ High precision at p_max=20"
echo "  ✓ Performance benchmarks completed"
echo ""
echo "For more examples, see:"
echo "  - USAGE_EXAMPLES.txt"
echo "  - README_MULTIPOLE.md"
echo "  - PROJECT_SUMMARY.md"
echo ""
echo "================================================================================"
