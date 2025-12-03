#!/bin/bash

echo "Testing different configurations..."
echo ""

# Test 1: Larger N
echo "=== Test 1: N=10000, p=6, grid=50x50 ==="
sed -i 's/N = [0-9]*/N = 10000/' src/main.f90
sed -i 's/p_order = [0-9]*/p_order = 6/' src/main.f90
sed -i 's/ngrid = [0-9]*/ngrid = 50/' src/main.f90
sed -i 's/nsteps = [0-9]*/nsteps = 1/' src/main.f90
make clean > /dev/null 2>&1 && make > /dev/null 2>&1
if [ $? -eq 0 ]; then
  ./bin/fmm_2d | grep -A 20 "PERFORMANCE SUMMARY"
fi
echo ""

# Test 2: Higher order
echo "=== Test 2: N=5000, p=8, grid=40x40 ==="
sed -i 's/N = [0-9]*/N = 5000/' src/main.f90
sed -i 's/p_order = [0-9]*/p_order = 8/' src/main.f90
sed -i 's/ngrid = [0-9]*/ngrid = 40/' src/main.f90
make clean > /dev/null 2>&1 && make > /dev/null 2>&1
if [ $? -eq 0 ]; then
  ./bin/fmm_2d | grep -A 20 "PERFORMANCE SUMMARY"
fi
echo ""

# Test 3: Lower order for speed
echo "=== Test 3: N=10000, p=4, grid=50x50 ==="
sed -i 's/N = [0-9]*/N = 10000/' src/main.f90
sed -i 's/p_order = [0-9]*/p_order = 4/' src/main.f90
sed -i 's/ngrid = [0-9]*/ngrid = 50/' src/main.f90
make clean > /dev/null 2>&1 && make > /dev/null 2>&1
if [ $? -eq 0 ]; then
  ./bin/fmm_2d | grep -A 20 "PERFORMANCE SUMMARY"
fi
echo ""

# Restore defaults
sed -i 's/N = [0-9]*/N = 5000/' src/main.f90
sed -i 's/p_order = [0-9]*/p_order = 6/' src/main.f90
sed -i 's/ngrid = [0-9]*/ngrid = 40/' src/main.f90
sed -i 's/nsteps = [0-9]*/nsteps = 3/' src/main.f90

echo "Tests complete!"
