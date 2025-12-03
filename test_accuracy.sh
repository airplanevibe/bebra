#!/bin/bash

# Script to test FMM accuracy and performance with different parameters

echo "=============================================="
echo "  FMM 2D Accuracy and Performance Tests"
echo "=============================================="
echo ""

# Create output directory
mkdir -p test_results

# Test different expansion orders
for p in 2 4 6 8 10; do
    echo ""
    echo "=========================================="
    echo "Testing with expansion order p = $p"
    echo "=========================================="

    # Modify source code to use this p value
    sed -i "s/p_order = [0-9]*/p_order = $p/" src/main.f90

    # Rebuild
    make clean > /dev/null 2>&1
    make > /dev/null 2>&1

    if [ $? -eq 0 ]; then
        # Run and save output
        ./bin/fmm_2d | tee test_results/output_p${p}.txt
    else
        echo "Build failed for p=$p"
    fi

    echo ""
done

# Test different grid sizes with fixed p=6
echo ""
echo "=========================================="
echo "Testing different grid sizes (p=6)"
echo "=========================================="

for ngrid in 20 30 40 50 60; do
    echo ""
    echo "Testing with grid size = ${ngrid}x${ngrid}"

    sed -i "s/p_order = [0-9]*/p_order = 6/" src/main.f90
    sed -i "s/ngrid = [0-9]*/ngrid = $ngrid/" src/main.f90

    make clean > /dev/null 2>&1
    make > /dev/null 2>&1

    if [ $? -eq 0 ]; then
        ./bin/fmm_2d | tee test_results/output_grid${ngrid}.txt
    else
        echo "Build failed for ngrid=$ngrid"
    fi

    echo ""
done

# Test different particle counts with fixed p=6, ngrid=40
echo ""
echo "=========================================="
echo "Testing different particle counts (p=6, grid=40x40)"
echo "=========================================="

for N in 1000 2000 5000 10000; do
    echo ""
    echo "Testing with N = $N particles"

    sed -i "s/N = [0-9]*/N = $N/" src/main.f90
    sed -i "s/p_order = [0-9]*/p_order = 6/" src/main.f90
    sed -i "s/ngrid = [0-9]*/ngrid = 40/" src/main.f90

    make clean > /dev/null 2>&1
    make > /dev/null 2>&1

    if [ $? -eq 0 ]; then
        ./bin/fmm_2d | tee test_results/output_N${N}.txt
    else
        echo "Build failed for N=$N"
    fi

    echo ""
done

echo ""
echo "=============================================="
echo "  All tests complete!"
echo "  Results saved in test_results/"
echo "=============================================="
