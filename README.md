# Multipole Expansion Tests

This repository contains implementations and tests for multipole expansion methods for the 2D Laplace equation.

## Quick Start

```bash
# Run comprehensive tests with default parameters
python3 multipole_2d_final.py

# Test for machine precision accuracy (p_max=30)
python3 multipole_2d_final.py --test-type accuracy --p-max 30 --n-sources 50

# Custom parameters
python3 multipole_2d_final.py --p-max 25 --n-sources 100 --n-eval 50
```

## Files

- `multipole_2d_final.py` - Main implementation with command-line interface
- `test_multipole_expansion.py` - Initial version using spherical harmonics
- `test_multipole_expansion_v2.py` - Simplified version with Legendre polynomials
- `README_MULTIPOLE.md` - Detailed documentation (Russian)

## Key Results

✓ Achieves near-machine precision (relative error ~10⁻¹⁵) at p_max=30
✓ Direct series expansion using Addition Theorem for Legendre polynomials
✓ Full control over all parameters via main() or command-line arguments
✓ Comprehensive accuracy and performance tests

## Requirements

- Python 3
- NumPy
- SciPy

See `README_MULTIPOLE.md` for detailed documentation.

---
This repository was initialized by Terragon.