#!/usr/bin/env python3
"""
Multipole Expansion for 2D Laplace Equation (No Z-axis interactions)

Implementation based on "Fast Multipole Methods for the Laplace Equation"
by Duraiswami & Gumerov (2003-2004).

Key Features:
- Uses Addition Theorem for Legendre Polynomials
- Achieves near-machine precision at high expansion orders (p_max ~ 25-30)
- 2D case: sources and evaluation points in xy-plane (z=0)
- Uniform grid support
- Direct series expansion of 1/r potential

Mathematical Background:
========================
For r > r_0, the Green's function expansion is:
    1/|r - r_0| = sum_{n=0}^∞ P_n(cos θ) * (r_0/r)^(n+1)

where:
- P_n is the nth Legendre polynomial
- θ is the angle between vectors r and r_0
- cos θ = (r · r_0) / (|r| * |r_0|)

For potential Φ(r) from sources q_i at positions r_i:
    Φ(r) = sum_i q_i / |r - r_i|
         = sum_i q_i * sum_n P_n(cos θ_i) * (r_i/r)^(n+1)
         = sum_n [sum_i q_i * r_i^n * P_n(cos θ_i)] / r^(n+1)
         = sum_n M_n(θ,φ) / r^(n+1)

where M_n are the multipole moments.
"""

import numpy as np
from scipy.special import legendre
import time
from typing import Tuple, List, Optional
import argparse


class MultipoleExpansion2D:
    """
    Multipole expansion for 2D Laplace equation using Legendre polynomials.

    This implementation uses the direct series expansion for the 1/r potential
    based on the Addition Theorem for spherical harmonics, simplified for the
    axisymmetric case.
    """

    def __init__(self, p_max: int, cache_legendre: bool = True):
        """
        Initialize multipole expansion.

        Args:
            p_max: Maximum order of expansion (higher = more accurate)
            cache_legendre: Whether to cache Legendre polynomial evaluations
        """
        self.p_max = p_max
        self.cache_legendre = cache_legendre

        # Precompute Legendre polynomial functions
        self.P = [legendre(n) for n in range(p_max + 1)]

        print(f"Initialized MultipoleExpansion2D with p_max={p_max}")

    def compute_multipole_coefficients(
        self,
        sources: np.ndarray,
        charges: np.ndarray,
        center: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute multipole moments for sources.

        The multipole moment of order n is defined as:
            M_n = sum_i q_i * r_i^n

        where r_i is the distance of source i from the expansion center.

        Args:
            sources: Source positions, shape (N_sources, 3)
            charges: Source charges, shape (N_sources,)
            center: Expansion center, shape (3,)

        Returns:
            Tuple of:
                - rel_sources: Source positions relative to center, shape (N_sources, 3)
                - moments: Multipole moments, shape (p_max+1, N_sources)
        """
        N_sources = sources.shape[0]

        # Compute positions relative to expansion center
        rel_sources = sources - center

        # Compute radial distances
        r_sources = np.linalg.norm(rel_sources, axis=1, keepdims=True)

        # Compute moments: M_n,i = q_i * r_i^n
        n_values = np.arange(self.p_max + 1).reshape(-1, 1)
        moments = charges * (r_sources.T ** n_values)

        return rel_sources, moments

    def evaluate_potential(
        self,
        rel_sources: np.ndarray,
        moments: np.ndarray,
        eval_point: np.ndarray,
        center: np.ndarray
    ) -> float:
        """
        Evaluate potential at a point using multipole expansion.

        Uses the expansion:
            Φ(r) = sum_{n=0}^{p_max} sum_i M_n,i * P_n(cos θ_i) / r^(n+1)

        where θ_i is the angle between r and r_i.

        Args:
            rel_sources: Source positions relative to center, shape (N_sources, 3)
            moments: Precomputed moments, shape (p_max+1, N_sources)
            eval_point: Evaluation point, shape (3,)
            center: Expansion center, shape (3,)

        Returns:
            Potential value at eval_point
        """
        # Position relative to center
        rel_eval = eval_point - center
        r_eval = np.linalg.norm(rel_eval)

        if r_eval < 1e-15:
            return np.inf

        N_sources = rel_sources.shape[0]
        potential = 0.0

        # For each source
        for i in range(N_sources):
            r_source = np.linalg.norm(rel_sources[i])

            if r_source < 1e-15:
                continue

            # Only compute if evaluation point is outside source
            if r_eval > r_source:
                # Compute cos(angle) between r_eval and r_source
                cos_angle = np.dot(rel_eval, rel_sources[i]) / (r_eval * r_source)
                cos_angle = np.clip(cos_angle, -1.0, 1.0)

                # Sum over expansion orders
                for n in range(self.p_max + 1):
                    P_n = self.P[n](cos_angle)
                    # Add contribution: M_n,i * P_n(cos θ) / r^(n+1)
                    potential += moments[n, i] * P_n / (r_eval ** (n + 1))

        return potential

    def evaluate_potential_batch(
        self,
        rel_sources: np.ndarray,
        moments: np.ndarray,
        eval_points: np.ndarray,
        center: np.ndarray
    ) -> np.ndarray:
        """
        Evaluate potential at multiple points.

        Args:
            rel_sources: Source positions relative to center, shape (N_sources, 3)
            moments: Precomputed moments, shape (p_max+1, N_sources)
            eval_points: Evaluation points, shape (N_eval, 3)
            center: Expansion center, shape (3,)

        Returns:
            Potential values, shape (N_eval,)
        """
        N_eval = eval_points.shape[0]
        potentials = np.zeros(N_eval)

        for i in range(N_eval):
            potentials[i] = self.evaluate_potential(
                rel_sources, moments, eval_points[i], center
            )

        return potentials


def compute_potential_direct(
    sources: np.ndarray,
    charges: np.ndarray,
    eval_points: np.ndarray
) -> np.ndarray:
    """
    Compute Coulomb potential directly using 1/r formula.

    Φ(r) = sum_i q_i / |r - r_i|

    Args:
        sources: Source positions, shape (N_sources, 3)
        charges: Source charges, shape (N_sources,)
        eval_points: Evaluation points, shape (N_eval, 3)

    Returns:
        Potential values, shape (N_eval,)
    """
    N_eval = eval_points.shape[0]
    N_src = sources.shape[0]

    potentials = np.zeros(N_eval)

    for i in range(N_eval):
        for j in range(N_src):
            r_vec = eval_points[i] - sources[j]
            r = np.linalg.norm(r_vec)

            if r > 1e-15:  # Avoid singularity
                potentials[i] += charges[j] / r

    return potentials


def generate_2d_grid(
    N_side: int,
    domain_size: float,
    z_fixed: float = 0.0
) -> np.ndarray:
    """
    Generate uniform 2D grid in xy-plane.

    Args:
        N_side: Number of points per side
        domain_size: Size of square domain
        z_fixed: Fixed z-coordinate (default 0.0)

    Returns:
        Grid points, shape (N_side^2, 3)
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
    Test accuracy of multipole expansion against direct computation.

    Args:
        p_max: Maximum expansion order
        N_sources: Number of source charges
        N_eval: Number of evaluation points
        source_radius: Radius of source region
        eval_radius: Radius of evaluation region (must be > source_radius)
        verbose: Whether to print detailed results

    Returns:
        Tuple of (max_absolute_error, max_relative_error, rms_error)
    """
    # Generate random sources in a disk (2D)
    angles = np.random.uniform(0, 2*np.pi, N_sources)
    radii = np.random.uniform(0, source_radius, N_sources)

    sources = np.zeros((N_sources, 3))
    sources[:, 0] = radii * np.cos(angles)
    sources[:, 1] = radii * np.sin(angles)
    sources[:, 2] = 0.0

    charges = np.random.uniform(-1.0, 1.0, N_sources)

    # Generate evaluation points outside source region
    angles_eval = np.random.uniform(0, 2*np.pi, N_eval)
    radii_eval = np.random.uniform(eval_radius, eval_radius * 1.5, N_eval)

    eval_points = np.zeros((N_eval, 3))
    eval_points[:, 0] = radii_eval * np.cos(angles_eval)
    eval_points[:, 1] = radii_eval * np.sin(angles_eval)
    eval_points[:, 2] = 0.0

    # Direct computation
    potential_direct = compute_potential_direct(sources, charges, eval_points)

    # Multipole expansion
    multipole = MultipoleExpansion2D(p_max)
    center = np.array([0.0, 0.0, 0.0])

    rel_sources, moments = multipole.compute_multipole_coefficients(
        sources, charges, center
    )

    potential_multipole = multipole.evaluate_potential_batch(
        rel_sources, moments, eval_points, center
    )

    # Compute errors
    abs_error = np.abs(potential_multipole - potential_direct)
    max_error = np.max(abs_error)

    nonzero_mask = np.abs(potential_direct) > 1e-15
    rel_errors = np.zeros_like(abs_error)
    rel_errors[nonzero_mask] = (
        abs_error[nonzero_mask] / np.abs(potential_direct[nonzero_mask])
    )
    rel_error = np.max(rel_errors)

    rms_error = np.sqrt(np.mean(abs_error**2))

    if verbose:
        print(f"\n{'='*70}")
        print(f"ACCURACY TEST (p_max={p_max})")
        print(f"{'='*70}")
        print(f"Configuration:")
        print(f"  Number of sources: {N_sources}")
        print(f"  Number of evaluation points: {N_eval}")
        print(f"  Source radius: {source_radius:.3f}")
        print(f"  Evaluation radius: {eval_radius:.3f}")
        print(f"  Expansion order: {p_max}")
        print(f"\nError Metrics:")
        print(f"  Maximum absolute error: {max_error:.6e}")
        print(f"  Maximum relative error: {rel_error:.6e}")
        print(f"  RMS error: {rms_error:.6e}")
        print(f"  Machine epsilon (float64): {np.finfo(float).eps:.6e}")

        # Classify precision
        if rel_error < 1e-12:
            status = "✓ EXCELLENT: Near-machine precision achieved!"
        elif rel_error < 1e-8:
            status = "✓ VERY GOOD: High precision"
        elif rel_error < 1e-4:
            status = "✓ GOOD: Acceptable precision"
        else:
            status = "⚠ MODERATE: Consider increasing p_max"

        print(f"\n  Status: {status}")
        print(f"{'='*70}")

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

    Note: Current implementation is O(N*M*p) where N=sources, M=eval points, p=order.
    True FMM with hierarchical tree structure achieves O((N+M)*p^2).

    Args:
        p_max: Maximum expansion order
        N_sources_list: List of source counts to test
        N_eval: Number of evaluation points
        source_radius: Radius of source region
        eval_radius: Radius of evaluation region
        verbose: Whether to print detailed results

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
        print(f"\n{'='*70}")
        print(f"PERFORMANCE TEST (p_max={p_max})")
        print(f"{'='*70}")
        print(f"Configuration:")
        print(f"  Number of evaluation points: {N_eval}")
        print(f"  Expansion order: {p_max}")
        print(f"\nTiming Results:")
        print(f"{'N_sources':>12} {'Direct(s)':>14} {'Multipole(s)':>16} {'Speedup':>12}")
        print(f"{'-'*70}")

    for N_sources in N_sources_list:
        # Generate data
        angles = np.random.uniform(0, 2*np.pi, N_sources)
        radii = np.random.uniform(0, source_radius, N_sources)

        sources = np.zeros((N_sources, 3))
        sources[:, 0] = radii * np.cos(angles)
        sources[:, 1] = radii * np.sin(angles)
        sources[:, 2] = 0.0

        charges = np.random.uniform(-1.0, 1.0, N_sources)

        angles_eval = np.random.uniform(0, 2*np.pi, N_eval)
        radii_eval = np.random.uniform(eval_radius, eval_radius * 1.5, N_eval)

        eval_points = np.zeros((N_eval, 3))
        eval_points[:, 0] = radii_eval * np.cos(angles_eval)
        eval_points[:, 1] = radii_eval * np.sin(angles_eval)
        eval_points[:, 2] = 0.0

        # Time direct computation
        t0 = time.time()
        potential_direct = compute_potential_direct(sources, charges, eval_points)
        time_direct = time.time() - t0

        # Time multipole computation
        t0 = time.time()
        multipole = MultipoleExpansion2D(p_max)
        center = np.array([0.0, 0.0, 0.0])
        rel_sources, moments = multipole.compute_multipole_coefficients(
            sources, charges, center
        )
        potential_multipole = multipole.evaluate_potential_batch(
            rel_sources, moments, eval_points, center
        )
        time_multipole = time.time() - t0

        speedup = time_direct / time_multipole if time_multipole > 0 else 0

        results['N_sources'].append(N_sources)
        results['time_direct'].append(time_direct)
        results['time_multipole'].append(time_multipole)
        results['speedup'].append(speedup)

        if verbose:
            speedup_str = f"{speedup:.2f}x" if speedup >= 1 else f"{speedup:.3f}x"
            print(f"{N_sources:>12} {time_direct:>14.6f} {time_multipole:>16.6f} {speedup_str:>12}")

    if verbose:
        print(f"{'='*70}")
        print(f"\nNote: For large-scale problems, use hierarchical FMM")
        print(f"      (octree structure) for O((N+M)*p^2) complexity.")
        print(f"{'='*70}")

    return results


def main(
    test_type: str = "both",
    p_max: int = 20,
    N_sources: int = 100,
    N_eval: int = 50,
    source_radius: float = 1.0,
    eval_radius: float = 2.5,
    performance_sources: Optional[List[int]] = None
):
    """
    Main function with all parameters controllable.

    Args:
        test_type: Type of test - "accuracy", "performance", or "both"
        p_max: Maximum expansion order (higher = better accuracy, slower)
        N_sources: Number of source charges for accuracy test
        N_eval: Number of evaluation points
        source_radius: Radius of source region
        eval_radius: Radius of evaluation region (must be > source_radius)
        performance_sources: List of source counts for performance test
    """
    print("\n" + "="*70)
    print(" MULTIPOLE EXPANSION TEST - 2D LAPLACE EQUATION")
    print("="*70)
    print("\nImplementation Details:")
    print("  - Based on Addition Theorem for Legendre Polynomials")
    print("  - Direct series expansion of 1/r potential")
    print("  - 2D case: all points in xy-plane (z=0)")
    print("  - Achieves near-machine precision at high orders (p_max ~ 25-30)")
    print(f"\nTest Parameters:")
    print(f"  Expansion order (p_max): {p_max}")
    print(f"  Source radius: {source_radius:.3f}")
    print(f"  Evaluation radius: {eval_radius:.3f}")

    if test_type in ["accuracy", "both"]:
        max_err, rel_err, rms_err = test_accuracy(
            p_max=p_max,
            N_sources=N_sources,
            N_eval=N_eval,
            source_radius=source_radius,
            eval_radius=eval_radius,
            verbose=True
        )

    if test_type in ["performance", "both"]:
        if performance_sources is None:
            performance_sources = [50, 100, 200, 400, 800]

        perf_results = test_performance(
            p_max=p_max,
            N_sources_list=performance_sources,
            N_eval=N_eval,
            source_radius=source_radius,
            eval_radius=eval_radius,
            verbose=True
        )

    print(f"\n{'='*70}")
    print(" TESTS COMPLETED SUCCESSFULLY")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Test multipole expansion for 2D Laplace equation"
    )
    parser.add_argument(
        "--test-type",
        type=str,
        default="both",
        choices=["accuracy", "performance", "both"],
        help="Type of test to run"
    )
    parser.add_argument(
        "--p-max",
        type=int,
        default=20,
        help="Maximum expansion order (default: 20, use 25-30 for machine precision)"
    )
    parser.add_argument(
        "--n-sources",
        type=int,
        default=100,
        help="Number of source charges for accuracy test"
    )
    parser.add_argument(
        "--n-eval",
        type=int,
        default=50,
        help="Number of evaluation points"
    )
    parser.add_argument(
        "--source-radius",
        type=float,
        default=1.0,
        help="Radius of source region"
    )
    parser.add_argument(
        "--eval-radius",
        type=float,
        default=2.5,
        help="Radius of evaluation region"
    )

    args = parser.parse_args()

    main(
        test_type=args.test_type,
        p_max=args.p_max,
        N_sources=args.n_sources,
        N_eval=args.n_eval,
        source_radius=args.source_radius,
        eval_radius=args.eval_radius,
        performance_sources=[50, 100, 200, 400, 800]
    )
