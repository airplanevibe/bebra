#!/usr/bin/env python3
"""
Generate precomputed coefficients for Fortran multipole code.
Generates Legendre polynomial coefficients and other constants.
"""

import numpy as np
from scipy.special import legendre, lpmv, factorial
import sys


def generate_legendre_coefficients(p_max, num_angles=360):
    """
    Generate Legendre polynomial values for a set of angles.

    Args:
        p_max: Maximum order of Legendre polynomials
        num_angles: Number of angles to evaluate (0 to π)

    Returns:
        Array of shape (p_max+1, num_angles) with P_n(cos(theta))
    """
    # Angles from 0 to π
    angles = np.linspace(0, np.pi, num_angles)
    cos_angles = np.cos(angles)

    # Generate Legendre polynomial values
    P_values = np.zeros((p_max + 1, num_angles))

    for n in range(p_max + 1):
        P_n = legendre(n)
        P_values[n, :] = P_n(cos_angles)

    return P_values, angles


def generate_factorial_table(n_max):
    """Generate factorial values up to n_max."""
    factorials = np.zeros(n_max + 1)
    for n in range(n_max + 1):
        factorials[n] = float(factorial(n))
    return factorials


def write_fortran_module(p_max, num_angles=360, output_file='legendre_data.f90'):
    """
    Generate Fortran module with precomputed coefficients.
    """
    P_values, angles = generate_legendre_coefficients(p_max, num_angles)
    factorials = generate_factorial_table(2 * p_max)

    with open(output_file, 'w') as f:
        f.write("! Auto-generated module with Legendre polynomial coefficients\n")
        f.write("! Generated for p_max = {}\n".format(p_max))
        f.write("!\n")
        f.write("module legendre_data\n")
        f.write("    implicit none\n")
        f.write("    integer, parameter :: dp = selected_real_kind(15, 307)\n")
        f.write("    integer, parameter :: p_max = {}\n".format(p_max))
        f.write("    integer, parameter :: num_angles = {}\n".format(num_angles))
        f.write("    \n")

        # Write angle table
        f.write("    ! Angle table (0 to π)\n")
        f.write("    real(dp), parameter, dimension(num_angles) :: angle_table = (/ &\n")
        for i in range(num_angles):
            if i > 0 and i % 2 == 0:
                f.write(" &\n        ")
            f.write("{:.13e}_dp".format(angles[i]))
            if i < num_angles - 1:
                f.write(", ")
        f.write(" /)\n\n")

        # Write Legendre polynomial table
        f.write("    ! Legendre polynomial values P_n(cos(theta))\n")
        f.write("    ! Dimensions: (0:p_max, num_angles)\n")
        f.write("    real(dp), parameter, dimension(0:p_max, num_angles) :: legendre_table = reshape((/ &\n")

        count = 0
        for n in range(p_max + 1):
            for i in range(num_angles):
                if count > 0 and count % 2 == 0:
                    f.write(" &\n        ")
                f.write("{:.13e}_dp".format(P_values[n, i]))
                if n < p_max or i < num_angles - 1:
                    f.write(", ")
                count += 1
        f.write(" /), shape(legendre_table))\n\n")

        # Write factorial table
        f.write("    ! Factorial table 0! to {}!\n".format(2 * p_max))
        f.write("    real(dp), parameter, dimension(0:{}) :: factorial_table = (/ &\n".format(2 * p_max))
        for i in range(len(factorials)):
            if i > 0 and i % 2 == 0:
                f.write(" &\n        ")
            f.write("{:.13e}_dp".format(factorials[i]))
            if i < len(factorials) - 1:
                f.write(", ")
        f.write(" /)\n\n")

        f.write("end module legendre_data\n")

    print(f"Generated {output_file}")
    print(f"  p_max = {p_max}")
    print(f"  num_angles = {num_angles}")
    print(f"  Legendre table size: {(p_max+1) * num_angles} values")
    print(f"  Factorial table size: {len(factorials)} values")


def main():
    """Main function."""
    # Default parameters
    p_max = 20
    num_angles = 360

    if len(sys.argv) > 1:
        p_max = int(sys.argv[1])
    if len(sys.argv) > 2:
        num_angles = int(sys.argv[2])

    print(f"Generating Legendre coefficients...")
    print(f"  Maximum order: p_max = {p_max}")
    print(f"  Number of angles: {num_angles}")
    print()

    # Generate Fortran module
    write_fortran_module(p_max, num_angles, 'legendre_data.f90')

    print()
    print("Done! Use this module in your Fortran code:")
    print("  use legendre_data")


if __name__ == "__main__":
    main()
