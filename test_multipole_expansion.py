#!/usr/bin/env python3
"""
Test Multipole Expansion for 2D case (no z-axis interactions)
Based on Fast Multipole Methods for Laplace Equation

Implementation:
- 1/r potential using multipole expansion
- Spherical harmonics with precomputed coefficients
- 2D situation: 3D potential but no interactions along z-axis
- Uniform grid
- Tests: machine precision accuracy at high orders and performance
"""

import numpy as np
from scipy.special import lpmv, factorial
import time
from typing import Tuple, List


class MultipoleExpansion2D:
    """
    Multipole expansion for 2D case using spherical harmonics.

    Based on formulas from the PDF:
    - S-expansion (multipole): valid outside source region
    - Uses spherical harmonics Y_n^m
    - For 1/r potential
    """

    def __init__(self, p_max: int):
        """
        Initialize multipole expansion.

        Args:
            p_max: Maximum order of expansion (truncation number)
        """
        self.p_max = p_max
        self._precompute_coefficients()

    def _precompute_coefficients(self):
        """Precompute normalization coefficients for spherical harmonics."""
        # Alpha coefficients: alpha_n^m = (-1)^n / sqrt((n-m)!(n+m)!)
        self.alpha = {}
        for n in range(self.p_max + 1):
            for m in range(-n, n + 1):
                abs_m = abs(m)
                if abs_m <= n:
                    # Normalization constant for associated Legendre functions
                    # Y_n^m(θ,φ) = (-1)^m * sqrt((2n+1)/(4π) * (n-m)!/(n+m)!) * P_n^m(cos θ) * e^(imφ)
                    coeff = ((-1.0) ** abs_m) / np.sqrt(
                        factorial(n - abs_m) * factorial(n + abs_m)
                    )
                    self.alpha[(n, m)] = coeff

    def spherical_harmonic(self, n: int, m: int, theta: float, phi: float) -> complex:
        """
        Compute spherical harmonic Y_n^m(θ, φ).

        Y_n^m(θ,φ) = sqrt((2n+1)/(4π) * (n-|m|)!/(n+|m|)!) * P_n^|m|(cos θ) * e^(imφ)

        Args:
            n: degree
            m: order (-n <= m <= n)
            theta: polar angle [0, π]
            phi: azimuthal angle [0, 2π]

        Returns:
            Complex value of Y_n^m(θ, φ)
        """
        if abs(m) > n:
            return 0.0 + 0.0j

        abs_m = abs(m)

        # Associated Legendre polynomial P_n^|m|(cos θ)
        cos_theta = np.cos(theta)
        cos_theta = np.clip(cos_theta, -1.0, 1.0)

        # P_n^|m|(cos θ) using scipy
        plm = lpmv(abs_m, n, cos_theta)

        # Normalization factor
        norm = np.sqrt(
            (2 * n + 1) / (4 * np.pi) *
            factorial(n - abs_m) / factorial(n + abs_m)
        )

        # For m < 0: Y_n^(-m) = (-1)^m * conj(Y_n^m)
        if m < 0:
            # Need to compute Y_n^|m| and take conjugate
            phase = np.exp(1j * abs_m * phi)
            result = ((-1) ** abs_m) * norm * plm * np.conj(phase)
        else:
            # m >= 0
            phase = np.exp(1j * m * phi)
            result = norm * plm * phase

        return result

    def cartesian_to_spherical(self, x: float, y: float, z: float) -> Tuple[float, float, float]:
        """
        Convert Cartesian coordinates to spherical coordinates.

        Args:
            x, y, z: Cartesian coordinates

        Returns:
            r: radial distance
            theta: polar angle [0, π]
            phi: azimuthal angle [0, 2π]
        """
        r = np.sqrt(x**2 + y**2 + z**2)

        if r < 1e-15:  # Avoid division by zero
            return 0.0, 0.0, 0.0

        theta = np.arccos(np.clip(z / r, -1.0, 1.0))
        phi = np.arctan2(y, x)

        return r, theta, phi

    def compute_multipole_coefficients(
        self,
        sources: np.ndarray,
        charges: np.ndarray,
        center: np.ndarray
    ) -> dict:
        """
        Compute multipole coefficients O_n^m for given sources.

        From PDF (slide 10): For S-expansion (multipole)
        1/(4π|r-r_0|) = sum_{n=0}^∞ sum_{m=-n}^n (1/(2n+1)) * R_n^(-m)(r_0) * S_n^m(r)

        where S_n^m(r) = r^(-n-1) * Y_n^m(θ,φ)
              R_n^m(r) = r^n * Y_n^m(θ,φ)

        For sources, we compute: O_n^m = sum_i q_i * R_n^(-m)(r_i) = sum_i q_i * r_i^n * Y_n^(-m)(θ_i, φ_i)

        Args:
            sources: Source positions (N x 3 array)
            charges: Source charges (N array)
            center: Expansion center (3 array)

        Returns:
            Dictionary of coefficients {(n, m): O_n^m}
        """
        O = {}

        for n in range(self.p_max + 1):
            for m in range(-n, n + 1):
                O[(n, m)] = 0.0 + 0.0j

        # For each source
        for i, (source, q) in enumerate(zip(sources, charges)):
            # Get position relative to center
            rel_pos = source - center
            r, theta, phi = self.cartesian_to_spherical(rel_pos[0], rel_pos[1], rel_pos[2])

            if r < 1e-15:
                continue

            # Accumulate multipole moments
            for n in range(self.p_max + 1):
                r_n = r ** n
                for m in range(-n, n + 1):
                    # O_n^m = sum_i q_i * r_i^n * Y_n^(-m)(θ_i, φ_i)
                    Y_n_minus_m = self.spherical_harmonic(n, -m, theta, phi)
                    O[(n, m)] += q * r_n * Y_n_minus_m

        return O

    def evaluate_multipole_expansion(
        self,
        O: dict,
        eval_point: np.ndarray,
        center: np.ndarray
    ) -> float:
        """
        Evaluate potential at a point using multipole expansion.

        From PDF (slide 10): The expansion formula is:
        1/(4π|r-r_0|) = (4π/r_0) * sum_{n=0}^∞ sum_{m=-n}^n (1/(2n+1)) * (r/r_0)^n * Y_n^(-m)(θ_0,φ_0) * Y_n^m(θ,φ)

        Rearranging for our case:
        1/|r-r_i| = 4π * sum_{n=0}^∞ sum_{m=-n}^n (1/(2n+1)) * Y_n^(-m)(θ_i,φ_i) * r_i^n * Y_n^m(θ,φ) / r^(n+1)

        So: Φ(r) = sum_i q_i / |r - r_i|
                 = 4π * sum_{n,m} [ sum_i q_i * r_i^n * Y_n^(-m)(θ_i,φ_i) ] * (1/(2n+1)) * Y_n^m(θ,φ) / r^(n+1)
                 = 4π * sum_{n,m} O_n^m * (1/(2n+1)) * Y_n^m(θ,φ) / r^(n+1)

        Args:
            O: Multipole coefficients O_n^m
            eval_point: Evaluation point (3 array)
            center: Expansion center (3 array)

        Returns:
            Potential value at eval_point
        """
        # Get position relative to center
        rel_pos = eval_point - center
        r, theta, phi = self.cartesian_to_spherical(rel_pos[0], rel_pos[1], rel_pos[2])

        if r < 1e-15:
            return np.inf

        potential = 0.0

        # Sum over degrees and orders
        for n in range(self.p_max + 1):
            r_factor = 1.0 / (r ** (n + 1))
            coeff_n = 4.0 * np.pi / (2 * n + 1)

            for m in range(-n, n + 1):
                Y_nm = self.spherical_harmonic(n, m, theta, phi)
                # Φ = 4π * sum_{n,m} (1/(2n+1)) * O_n^m * r^(-n-1) * Y_n^m
                potential += (coeff_n * O[(n, m)] * r_factor * Y_nm).real

        return potential


def compute_potential_direct(
    sources: np.ndarray,
    charges: np.ndarray,
    eval_points: np.ndarray
) -> np.ndarray:
    """
    Compute potential directly using 1/r formula.

    Φ(r) = sum_i q_i / |r - r_i|

    Args:
        sources: Source positions (N_src x 3)
        charges: Source charges (N_src,)
        eval_points: Evaluation points (N_eval x 3)

    Returns:
        Potential values at evaluation points (N_eval,)
    """
    N_eval = eval_points.shape[0]
    N_src = sources.shape[0]

    potentials = np.zeros(N_eval)

    for i in range(N_eval):
        for j in range(N_src):
            r_vec = eval_points[i] - sources[j]
            r = np.linalg.norm(r_vec)

            if r > 1e-15:  # Avoid self-interaction
                potentials[i] += charges[j] / r

    return potentials


def generate_2d_grid(
    N_side: int,
    domain_size: float,
    z_fixed: float = 0.0
) -> np.ndarray:
    """
    Generate 2D uniform grid (z is fixed - no z-axis interactions).

    Args:
        N_side: Number of points per side
        domain_size: Size of the square domain
        z_fixed: Fixed z-coordinate (default 0)

    Returns:
        Grid points (N_side^2 x 3)
    """
    x = np.linspace(-domain_size/2, domain_size/2, N_side)
    y = np.linspace(-domain_size/2, domain_size/2, N_side)

    X, Y = np.meshgrid(x, y)

    points = np.zeros((N_side * N_side, 3))
    points[:, 0] = X.flatten()
    points[:, 1] = Y.flatten()
    points[:, 2] = z_fixed

    return points


def test_accuracy(
    p_max: int,
    N_sources: int,
    N_eval: int,
    source_radius: float,
    eval_radius: float,
    verbose: bool = True
) -> Tuple[float, float, float]:
    """
    Test accuracy of multipole expansion vs direct computation.

    Args:
        p_max: Maximum expansion order
        N_sources: Number of source charges
        N_eval: Number of evaluation points
        source_radius: Radius of source region
        eval_radius: Radius of evaluation region
        verbose: Print detailed results

    Returns:
        max_error: Maximum absolute error
        rel_error: Maximum relative error
        rms_error: RMS error
    """
    # Generate random sources in a circle (2D)
    angles = np.random.uniform(0, 2*np.pi, N_sources)
    radii = np.random.uniform(0, source_radius, N_sources)

    sources = np.zeros((N_sources, 3))
    sources[:, 0] = radii * np.cos(angles)
    sources[:, 1] = radii * np.sin(angles)
    sources[:, 2] = 0.0  # z = 0 (2D case)

    # Random charges
    charges = np.random.uniform(-1.0, 1.0, N_sources)

    # Generate evaluation points outside source region
    angles_eval = np.random.uniform(0, 2*np.pi, N_eval)
    radii_eval = np.random.uniform(eval_radius, eval_radius * 1.5, N_eval)

    eval_points = np.zeros((N_eval, 3))
    eval_points[:, 0] = radii_eval * np.cos(angles_eval)
    eval_points[:, 1] = radii_eval * np.sin(angles_eval)
    eval_points[:, 2] = 0.0  # z = 0 (2D case)

    # Compute using direct method
    potential_direct = compute_potential_direct(sources, charges, eval_points)

    # Compute using multipole expansion
    multipole = MultipoleExpansion2D(p_max)
    center = np.array([0.0, 0.0, 0.0])

    # Get multipole coefficients
    O = multipole.compute_multipole_coefficients(sources, charges, center)

    # Evaluate at each point
    potential_multipole = np.zeros(N_eval)
    for i, point in enumerate(eval_points):
        potential_multipole[i] = multipole.evaluate_multipole_expansion(O, point, center)

    # Compute errors
    abs_error = np.abs(potential_multipole - potential_direct)
    max_error = np.max(abs_error)

    # Relative error (avoid division by zero)
    nonzero_mask = np.abs(potential_direct) > 1e-15
    rel_errors = np.zeros_like(abs_error)
    rel_errors[nonzero_mask] = abs_error[nonzero_mask] / np.abs(potential_direct[nonzero_mask])
    rel_error = np.max(rel_errors)

    rms_error = np.sqrt(np.mean(abs_error**2))

    if verbose:
        print(f"\n{'='*60}")
        print(f"Accuracy Test (p_max={p_max})")
        print(f"{'='*60}")
        print(f"Number of sources: {N_sources}")
        print(f"Number of evaluation points: {N_eval}")
        print(f"Source radius: {source_radius:.3f}")
        print(f"Evaluation radius: {eval_radius:.3f}")
        print(f"\nResults:")
        print(f"  Maximum absolute error: {max_error:.3e}")
        print(f"  Maximum relative error: {rel_error:.3e}")
        print(f"  RMS error: {rms_error:.3e}")
        print(f"  Machine epsilon (double): {np.finfo(float).eps:.3e}")

        if rel_error < 1e-10:
            print(f"  ✓ Excellent precision achieved!")
        elif rel_error < 1e-6:
            print(f"  ✓ Good precision")
        else:
            print(f"  ⚠ Moderate precision")

    return max_error, rel_error, rms_error


def test_performance(
    p_max: int,
    N_sources_list: List[int],
    N_eval: int,
    source_radius: float,
    eval_radius: float,
    verbose: bool = True
) -> dict:
    """
    Test performance: multipole vs direct computation.

    Args:
        p_max: Maximum expansion order
        N_sources_list: List of source counts to test
        N_eval: Number of evaluation points
        source_radius: Radius of source region
        eval_radius: Radius of evaluation region
        verbose: Print detailed results

    Returns:
        Dictionary with timing results
    """
    results = {
        'N_sources': [],
        'time_direct': [],
        'time_multipole': [],
        'speedup': []
    }

    if verbose:
        print(f"\n{'='*60}")
        print(f"Performance Test (p_max={p_max})")
        print(f"{'='*60}")
        print(f"Number of evaluation points: {N_eval}")
        print(f"\n{'N_sources':>12} {'Direct(s)':>12} {'Multipole(s)':>15} {'Speedup':>10}")
        print(f"{'-'*60}")

    for N_sources in N_sources_list:
        # Generate sources
        angles = np.random.uniform(0, 2*np.pi, N_sources)
        radii = np.random.uniform(0, source_radius, N_sources)

        sources = np.zeros((N_sources, 3))
        sources[:, 0] = radii * np.cos(angles)
        sources[:, 1] = radii * np.sin(angles)
        sources[:, 2] = 0.0

        charges = np.random.uniform(-1.0, 1.0, N_sources)

        # Generate evaluation points
        angles_eval = np.random.uniform(0, 2*np.pi, N_eval)
        radii_eval = np.random.uniform(eval_radius, eval_radius * 1.5, N_eval)

        eval_points = np.zeros((N_eval, 3))
        eval_points[:, 0] = radii_eval * np.cos(angles_eval)
        eval_points[:, 1] = radii_eval * np.sin(angles_eval)
        eval_points[:, 2] = 0.0

        # Time direct method
        t0 = time.time()
        potential_direct = compute_potential_direct(sources, charges, eval_points)
        time_direct = time.time() - t0

        # Time multipole method
        t0 = time.time()
        multipole = MultipoleExpansion2D(p_max)
        center = np.array([0.0, 0.0, 0.0])
        O = multipole.compute_multipole_coefficients(sources, charges, center)

        potential_multipole = np.zeros(N_eval)
        for i, point in enumerate(eval_points):
            potential_multipole[i] = multipole.evaluate_multipole_expansion(O, point, center)
        time_multipole = time.time() - t0

        speedup = time_direct / time_multipole if time_multipole > 0 else 0

        results['N_sources'].append(N_sources)
        results['time_direct'].append(time_direct)
        results['time_multipole'].append(time_multipole)
        results['speedup'].append(speedup)

        if verbose:
            print(f"{N_sources:>12} {time_direct:>12.4f} {time_multipole:>15.4f} {speedup:>10.2f}x")

    return results


def main(
    test_type: str = "both",
    p_max: int = 10,
    N_sources: int = 100,
    N_eval: int = 50,
    source_radius: float = 1.0,
    eval_radius: float = 2.0,
    performance_sources: List[int] = None
):
    """
    Main function to run tests.

    Args:
        test_type: "accuracy", "performance", or "both"
        p_max: Maximum expansion order
        N_sources: Number of source charges (for accuracy test)
        N_eval: Number of evaluation points
        source_radius: Radius of source region
        eval_radius: Radius of evaluation region
        performance_sources: List of source counts for performance test
    """
    print("\n" + "="*60)
    print("Multipole Expansion Test - 2D Case")
    print("="*60)
    print(f"\nParameters:")
    print(f"  Expansion order (p_max): {p_max}")
    print(f"  Source radius: {source_radius}")
    print(f"  Evaluation radius: {eval_radius}")

    if test_type in ["accuracy", "both"]:
        # Run accuracy test
        max_err, rel_err, rms_err = test_accuracy(
            p_max=p_max,
            N_sources=N_sources,
            N_eval=N_eval,
            source_radius=source_radius,
            eval_radius=eval_radius,
            verbose=True
        )

    if test_type in ["performance", "both"]:
        # Run performance test
        if performance_sources is None:
            performance_sources = [50, 100, 200, 400]

        perf_results = test_performance(
            p_max=p_max,
            N_sources_list=performance_sources,
            N_eval=N_eval,
            source_radius=source_radius,
            eval_radius=eval_radius,
            verbose=True
        )

    print(f"\n{'='*60}")
    print("Tests completed!")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    # Default parameters - all controllable from main()
    main(
        test_type="both",           # Test accuracy and performance
        p_max=15,                   # High order for near-machine precision
        N_sources=100,              # Number of source charges
        N_eval=50,                  # Number of evaluation points
        source_radius=1.0,          # Source region radius
        eval_radius=2.5,            # Evaluation region radius
        performance_sources=[50, 100, 200, 400, 800]  # For performance test
    )
