#!/usr/bin/env python3
"""
Test Multipole Expansion for 2D case - Version 2
Using Addition Theorem for Legendre Polynomials directly

Based on slide 7 from PDF:
G(r - r_0) = 1/|r - r_0| = sum_{n=0}^∞ P_n(cos θ) * (r/r_0)^n / r_0   for r < r_0
           = sum_{n=0}^∞ P_n(cos θ) * (r_0/r)^(n+1)                    for r > r_0

where θ is the angle between vectors r and r_0
"""

import numpy as np
from scipy.special import legendre
import time
from typing import Tuple, List


class MultipoleExpansion2DSimple:
    """
    Simplified multipole expansion for 2D using Legendre polynomials.
    """

    def __init__(self, p_max: int):
        """
        Initialize multipole expansion.

        Args:
            p_max: Maximum order of expansion
        """
        self.p_max = p_max
        # Precompute Legendre polynomials
        self.P = [legendre(n) for n in range(p_max + 1)]

    def compute_multipole_coefficients(
        self,
        sources: np.ndarray,
        charges: np.ndarray,
        center: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute multipole moments M_n for sources.

        For each source at r_i with charge q_i:
        M_n = sum_i q_i * r_i^n * (x_i/r_i, y_i/r_i, z_i/r_i)^n

        We'll store the contribution more carefully using spherical coordinates.

        Args:
            sources: Source positions (N x 3)
            charges: Source charges (N,)
            center: Expansion center (3,)

        Returns:
            Tuple of (radial_moments, angular_info)
        """
        N_sources = sources.shape[0]

        # Store source information relative to center
        rel_sources = sources - center
        r_sources = np.linalg.norm(rel_sources, axis=1)

        # Compute moments: for multipole expansion
        # We need to store q_i * r_i^n and direction vectors
        moments = np.zeros((self.p_max + 1, N_sources), dtype=np.float64)

        for n in range(self.p_max + 1):
            for i in range(N_sources):
                if r_sources[i] > 1e-15:
                    moments[n, i] = charges[i] * (r_sources[i] ** n)

        return rel_sources, moments

    def evaluate_potential(
        self,
        rel_sources: np.ndarray,
        moments: np.ndarray,
        eval_point: np.ndarray,
        center: np.ndarray
    ) -> float:
        """
        Evaluate potential using multipole expansion with Addition Theorem.

        From PDF slide 7:
        1/|r - r_0| = sum_{n=0}^∞ P_n(μ) * (r_0/r)^(n+1)  for r > r_0
        where μ = cos(θ) = (r · r_0) / (|r| * |r_0|)

        Args:
            rel_sources: Source positions relative to center (N x 3)
            moments: Precomputed moments (p_max+1 x N)
            eval_point: Evaluation point (3,)
            center: Expansion center (3,)

        Returns:
            Potential at eval_point
        """
        # Position relative to center
        rel_eval = eval_point - center
        r_eval = np.linalg.norm(rel_eval)

        if r_eval < 1e-15:
            return np.inf

        potential = 0.0
        N_sources = rel_sources.shape[0]

        # For each source
        for i in range(N_sources):
            r_source = np.linalg.norm(rel_sources[i])

            if r_source < 1e-15:
                continue

            # Check if eval point is outside source region
            if r_eval > r_source:
                # cos(angle) between r_eval and r_source
                cos_angle = np.dot(rel_eval, rel_sources[i]) / (r_eval * r_source)
                cos_angle = np.clip(cos_angle, -1.0, 1.0)

                # Sum over Legendre polynomials
                contrib = 0.0
                r_ratio = r_source / r_eval

                for n in range(self.p_max + 1):
                    # P_n(cos θ)
                    P_n = self.P[n](cos_angle)
                    # moments[n, i] = q_i * r_source^n
                    # Contribution: q_i * r_source^n * P_n(cos θ) / r_eval^(n+1)
                    contrib += moments[n, i] * P_n / (r_eval ** (n + 1))

                potential += contrib

        return potential


def compute_potential_direct(
    sources: np.ndarray,
    charges: np.ndarray,
    eval_points: np.ndarray
) -> np.ndarray:
    """Direct computation of Coulomb potential."""
    N_eval = eval_points.shape[0]
    N_src = sources.shape[0]

    potentials = np.zeros(N_eval)

    for i in range(N_eval):
        for j in range(N_src):
            r_vec = eval_points[i] - sources[j]
            r = np.linalg.norm(r_vec)

            if r > 1e-15:
                potentials[i] += charges[j] / r

    return potentials


def test_accuracy(
    p_max: int,
    N_sources: int,
    N_eval: int,
    source_radius: float,
    eval_radius: float,
    verbose: bool = True
) -> Tuple[float, float, float]:
    """Test accuracy of multipole expansion."""
    # Generate sources in a circle (2D)
    angles = np.random.uniform(0, 2*np.pi, N_sources)
    radii = np.random.uniform(0, source_radius, N_sources)

    sources = np.zeros((N_sources, 3))
    sources[:, 0] = radii * np.cos(angles)
    sources[:, 1] = radii * np.sin(angles)
    sources[:, 2] = 0.0

    charges = np.random.uniform(-1.0, 1.0, N_sources)

    # Evaluation points outside
    angles_eval = np.random.uniform(0, 2*np.pi, N_eval)
    radii_eval = np.random.uniform(eval_radius, eval_radius * 1.5, N_eval)

    eval_points = np.zeros((N_eval, 3))
    eval_points[:, 0] = radii_eval * np.cos(angles_eval)
    eval_points[:, 1] = radii_eval * np.sin(angles_eval)
    eval_points[:, 2] = 0.0

    # Direct computation
    potential_direct = compute_potential_direct(sources, charges, eval_points)

    # Multipole expansion
    multipole = MultipoleExpansion2DSimple(p_max)
    center = np.array([0.0, 0.0, 0.0])

    rel_sources, moments = multipole.compute_multipole_coefficients(sources, charges, center)

    potential_multipole = np.zeros(N_eval)
    for i, point in enumerate(eval_points):
        potential_multipole[i] = multipole.evaluate_potential(rel_sources, moments, point, center)

    # Errors
    abs_error = np.abs(potential_multipole - potential_direct)
    max_error = np.max(abs_error)

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
    """Test performance comparison."""
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

        # Time direct
        t0 = time.time()
        potential_direct = compute_potential_direct(sources, charges, eval_points)
        time_direct = time.time() - t0

        # Time multipole
        t0 = time.time()
        multipole = MultipoleExpansion2DSimple(p_max)
        center = np.array([0.0, 0.0, 0.0])
        rel_sources, moments = multipole.compute_multipole_coefficients(sources, charges, center)

        potential_multipole = np.zeros(N_eval)
        for i, point in enumerate(eval_points):
            potential_multipole[i] = multipole.evaluate_potential(rel_sources, moments, point, center)
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
    Main function - all parameters controllable.

    Args:
        test_type: "accuracy", "performance", or "both"
        p_max: Maximum expansion order
        N_sources: Number of source charges
        N_eval: Number of evaluation points
        source_radius: Radius of source region
        eval_radius: Radius of evaluation region
        performance_sources: List of source counts for performance test
    """
    print("\n" + "="*60)
    print("Multipole Expansion Test - 2D Case (Simplified)")
    print("="*60)
    print(f"\nParameters:")
    print(f"  Expansion order (p_max): {p_max}")
    print(f"  Source radius: {source_radius}")
    print(f"  Evaluation radius: {eval_radius}")

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
    main(
        test_type="both",
        p_max=20,                   # Higher order for better precision
        N_sources=100,
        N_eval=50,
        source_radius=1.0,
        eval_radius=2.5,
        performance_sources=[50, 100, 200, 400, 800]
    )
