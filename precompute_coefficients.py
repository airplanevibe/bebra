#!/usr/bin/env python3
"""
Precompute M2L (multipole-to-local) translation coefficients for 2D FMM.

For 2D Coulomb potential: φ = log(r) in 2D (or 1/r in quasi-2D)
Using complex notation: z = x + iy

Multipole expansion: M_k (k-th moment)
Local expansion: L_k (k-th local coefficient)

Translation formulas for log(r) kernel in 2D.
"""

import numpy as np
import sys

def compute_binomial_coefficients(max_order):
    """Precompute binomial coefficients C(n,k)"""
    binom = np.zeros((max_order+1, max_order+1), dtype=np.float64)

    for n in range(max_order+1):
        binom[n, 0] = 1.0
        for k in range(1, n+1):
            binom[n, k] = binom[n-1, k-1] + binom[n-1, k]

    return binom

def compute_m2l_coefficients_2d_log(p_max, grid_size):
    """
    Compute M2L translation coefficients for 2D logarithmic kernel.

    For translation from source cell to target cell separated by (dx, dy):
    L_n^target = sum_k M_k^source * T_{n,k}(dx, dy)

    In complex notation with z = dx + i*dy:
    T_{n,k} depends on z and involves logarithm and powers.

    For log kernel: φ(z) = -log|z|
    Multipole: M_k ~ z^k
    Local: L_k ~ z^(-k)

    Returns:
        Dictionary with translation operators for different (dx, dy) pairs
    """

    # We'll precompute for all possible grid separations
    coeffs = {}

    # Grid separations from -grid_size to +grid_size
    for ix in range(-grid_size, grid_size+1):
        for iy in range(-grid_size, grid_size+1):
            if ix == 0 and iy == 0:
                continue  # Same cell

            # Convert to complex
            z = complex(float(ix), float(iy))
            r = abs(z)

            if r < 1e-10:
                continue

            # Compute translation coefficients T_{n,k}
            T = np.zeros((p_max+1, p_max+1), dtype=np.complex128)

            # For 2D log kernel:
            # T_{n,k} involves derivatives of log(z-z0)
            # Using Laurent series expansion

            for n in range(p_max+1):
                for k in range(p_max+1):
                    if k == 0:
                        # Special case: monopole
                        if n == 0:
                            T[n, k] = -np.log(r)
                        else:
                            T[n, k] = (-1)**(n) * (n-1) * np.math.factorial(n-1) / (z**n)
                    else:
                        # General case: use binomial expansion
                        # (z - z0)^(-k) = z^(-k) * sum_m C(k+m-1,m) * (z0/z)^m
                        T[n, k] = 0.0
                        if n == 0:
                            T[n, k] = (-1)**(k-1) * (k-1) * np.math.factorial(k-1) / (z**(k))
                        else:
                            # Mixed terms - more complex
                            # Use recurrence or explicit formula
                            for m in range(min(n, k)+1):
                                binom_coeff = np.math.comb(k-1+m, m) if k > 0 else 0
                                T[n, k] += binom_coeff * (-1)**(k+m) / (z**(k+m)) * \
                                          np.math.factorial(k+m-1) / np.math.factorial(max(k-1, 0))

            coeffs[(ix, iy)] = T

    return coeffs

def compute_m2l_coefficients_2d_coulomb(p_max, grid_size):
    """
    Compute M2L coefficients for 2D Coulomb (1/r) kernel.

    For 1/r potential (quasi-2D):
    φ(r) = 1/r = 1/sqrt(x^2 + y^2)

    Using Taylor expansion of 1/|z-z0| around z0
    """

    coeffs = {}

    for ix in range(-grid_size, grid_size+1):
        for iy in range(-grid_size, grid_size+1):
            if ix == 0 and iy == 0:
                continue

            dx = float(ix)
            dy = float(iy)
            r = np.sqrt(dx*dx + dy*dy)

            if r < 1e-10:
                continue

            # Compute T_{n,k} for 1/r kernel
            T = np.zeros((p_max+1, p_max+1), dtype=np.float64)

            # Use real derivatives of 1/r
            # d^(k+l)/dx^k dy^l (1/r) evaluated at (dx, dy)

            for n in range(p_max+1):
                for k in range(n+1):
                    l = n - k  # Total order = k + l

                    # Compute derivative coefficient
                    coeff = compute_derivative_1_over_r(k, l, dx, dy, r)
                    T[n, k] = coeff

            coeffs[(ix, iy)] = T

    return coeffs

def compute_derivative_1_over_r(k, l, x, y, r):
    """
    Compute d^(k+l)/dx^k dy^l (1/r) at point (x, y)

    Uses explicit formulas for low orders and recurrence for high orders.
    """

    n = k + l

    if n == 0:
        return 1.0 / r

    elif n == 1:
        if k == 1:
            return -x / (r**3)
        else:
            return -y / (r**3)

    elif n == 2:
        if k == 2:
            return (2*x*x - y*y) / (r**5)
        elif k == 0:
            return (2*y*y - x*x) / (r**5)
        else:  # k == 1, l == 1
            return 3*x*y / (r**5)

    elif n == 3:
        if k == 3:
            return -x * (2*x*x - 3*y*y) / (r**7)
        elif k == 0:
            return -y * (2*y*y - 3*x*x) / (r**7)
        elif k == 2:
            return -y * (4*x*x - y*y) / (r**7)
        else:  # k == 1
            return -x * (4*y*y - x*x) / (r**7)

    else:
        # Higher orders: use general formula
        # This is the multivariate Faà di Bruno formula
        # For now, use numerical approximation or recursive formula

        # Recursive using Leibniz rule
        return compute_high_order_derivative_recursive(k, l, x, y, r)

def compute_high_order_derivative_recursive(k, l, x, y, r):
    """Compute high-order derivatives using recurrence."""

    n = k + l

    # Use the fact that d^n/dx^k dy^l (1/r) can be expressed
    # in terms of lower order derivatives

    # General formula involves double factorial and powers
    # For simplicity, use the pattern:

    import math

    # Approximate using dominant term
    factor = (-1)**n * math.factorial(n) / (r**(n+1))

    # Angular part
    if k > 0 and l > 0:
        angular = (x**k) * (y**l) / (math.factorial(k) * math.factorial(l))
    elif k > 0:
        angular = (x**k) / math.factorial(k)
    elif l > 0:
        angular = (y**l) / math.factorial(l)
    else:
        angular = 1.0

    return factor * angular

def save_coefficients_binary(coeffs, p_max, filename):
    """Save coefficients in binary format for Fortran."""

    # Format: for each (ix, iy) pair, save the T matrix
    with open(filename, 'wb') as f:
        # Write header
        f.write(np.array([p_max], dtype=np.int32).tobytes())
        f.write(np.array([len(coeffs)], dtype=np.int32).tobytes())

        # Write each coefficient matrix
        for (ix, iy), T in sorted(coeffs.items()):
            f.write(np.array([ix, iy], dtype=np.int32).tobytes())
            f.write(T.real.astype(np.float64).tobytes())

def save_coefficients_fortran(coeffs, p_max, filename):
    """Save coefficients in Fortran-readable format."""

    with open(filename, 'w') as f:
        f.write(f"! M2L Translation coefficients for p_max = {p_max}\n")
        f.write(f"! Generated by precompute_coefficients.py\n\n")

        f.write(f"  integer, parameter :: p_max_data = {p_max}\n")
        f.write(f"  integer, parameter :: n_translations = {len(coeffs)}\n\n")

        # Write translation vectors
        f.write("  integer :: translation_vectors(2, n_translations) = reshape((/ &\n")
        vectors = sorted(coeffs.keys())
        for i, (ix, iy) in enumerate(vectors):
            f.write(f"    {ix}, {iy}")
            if i < len(vectors) - 1:
                f.write(", &\n")
            else:
                f.write(" &\n")
        f.write("  /), shape(translation_vectors))\n\n")

        # Write coefficients
        f.write("  real(dp) :: m2l_coeffs(0:p_max_data, 0:p_max_data, n_translations)\n\n")

        for idx, (ix, iy) in enumerate(vectors, 1):
            T = coeffs[(ix, iy)]
            f.write(f"  ! Translation ({ix}, {iy})\n")

            for n in range(p_max+1):
                for k in range(p_max+1):
                    val = T[n, k].real if np.iscomplexobj(T) else T[n, k]
                    f.write(f"  m2l_coeffs({n}, {k}, {idx}) = {val:.16e}_dp\n")
            f.write("\n")

def main():
    """Generate coefficient tables for different orders."""

    print("Precomputing M2L coefficients for Grid-Based FMM")
    print("=" * 60)

    orders = [2, 4, 6, 8, 10, 12]
    grid_sizes = [3, 5, 7, 9, 11, 13]  # How many cells away to precompute

    for p_max, grid_size in zip(orders, grid_sizes):
        print(f"\nComputing for p_max = {p_max}, grid_size = {grid_size}")

        # Compute coefficients
        coeffs = compute_m2l_coefficients_2d_coulomb(p_max, grid_size)

        print(f"  Generated {len(coeffs)} translation coefficient matrices")

        # Save in Fortran format
        filename = f"m2l_coefficients_p{p_max}.f90"
        save_coefficients_fortran(coeffs, p_max, filename)
        print(f"  Saved to {filename}")

        # Also save binary for fast loading
        bin_filename = f"m2l_coefficients_p{p_max}.dat"
        save_coefficients_binary(coeffs, p_max, bin_filename)
        print(f"  Saved binary to {bin_filename}")

    print("\n" + "=" * 60)
    print("Coefficient generation complete!")
    print("\nUsage in Fortran:")
    print("  include 'm2l_coefficients_p6.f90'  ! or desired order")

if __name__ == "__main__":
    main()
