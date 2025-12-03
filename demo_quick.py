#!/usr/bin/env python3
"""
Быстрая демонстрация мультипольного разложения
Quick demonstration of multipole expansion
"""

import numpy as np
from multipole_2d_final import MultipoleExpansion2D, compute_potential_direct

def demo():
    """Быстрая демонстрация с примером"""

    print("\n" + "="*70)
    print(" QUICK DEMO: Multipole Expansion for 2D Laplace Equation")
    print("="*70)

    # Параметры / Parameters
    p_max = 25
    N_sources = 50
    source_radius = 1.0

    print(f"\nConfiguration:")
    print(f"  Expansion order: p_max = {p_max}")
    print(f"  Number of sources: {N_sources}")
    print(f"  Source radius: {source_radius}")

    # Создать источники в круге / Create sources in a disk
    print(f"\nGenerating {N_sources} random sources in a disk...")
    angles = np.random.uniform(0, 2*np.pi, N_sources)
    radii = np.random.uniform(0, source_radius, N_sources)

    sources = np.zeros((N_sources, 3))
    sources[:, 0] = radii * np.cos(angles)
    sources[:, 1] = radii * np.sin(angles)
    sources[:, 2] = 0.0  # 2D: z = 0

    charges = np.random.uniform(-1.0, 1.0, N_sources)
    print(f"  Charges range: [{charges.min():.2f}, {charges.max():.2f}]")

    # Точки оценки / Evaluation points
    eval_points = np.array([
        [3.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [-3.0, 0.0, 0.0],
        [0.0, -3.0, 0.0],
        [2.5, 2.5, 0.0]
    ])
    print(f"\nEvaluation at {len(eval_points)} points:")
    for i, pt in enumerate(eval_points):
        print(f"  Point {i+1}: ({pt[0]:.1f}, {pt[1]:.1f}, {pt[2]:.1f})")

    # Прямое вычисление / Direct computation
    print(f"\nComputing potential directly...")
    potential_direct = compute_potential_direct(sources, charges, eval_points)

    # Мультипольное разложение / Multipole expansion
    print(f"Computing potential using multipole expansion (p_max={p_max})...")
    multipole = MultipoleExpansion2D(p_max)
    center = np.array([0.0, 0.0, 0.0])

    rel_sources, moments = multipole.compute_multipole_coefficients(
        sources, charges, center
    )

    potential_multipole = multipole.evaluate_potential_batch(
        rel_sources, moments, eval_points, center
    )

    # Сравнение результатов / Compare results
    print(f"\n" + "="*70)
    print(" RESULTS")
    print("="*70)
    print(f"\n{'Point':>7} {'Direct':>15} {'Multipole':>15} {'Error':>15}")
    print("-"*70)

    for i in range(len(eval_points)):
        error = abs(potential_multipole[i] - potential_direct[i])
        rel_error = error / abs(potential_direct[i]) if abs(potential_direct[i]) > 1e-15 else 0
        print(f"{i+1:>7} {potential_direct[i]:>15.10f} {potential_multipole[i]:>15.10f} {error:>15.2e}")

    # Статистика ошибок / Error statistics
    abs_errors = np.abs(potential_multipole - potential_direct)
    max_error = np.max(abs_errors)

    nonzero_mask = np.abs(potential_direct) > 1e-15
    rel_errors = np.zeros_like(abs_errors)
    rel_errors[nonzero_mask] = abs_errors[nonzero_mask] / np.abs(potential_direct[nonzero_mask])
    max_rel_error = np.max(rel_errors)

    print("-"*70)
    print(f"\nError Statistics:")
    print(f"  Maximum absolute error: {max_error:.6e}")
    print(f"  Maximum relative error: {max_rel_error:.6e}")
    print(f"  Machine epsilon (float64): {np.finfo(float).eps:.6e}")

    if max_rel_error < 1e-10:
        print(f"\n  ✓ Excellent precision! Very close to machine accuracy.")
    elif max_rel_error < 1e-6:
        print(f"\n  ✓ Good precision!")
    else:
        print(f"\n  ⚠ Moderate precision. Consider increasing p_max.")

    print("\n" + "="*70)
    print(" DEMO COMPLETED")
    print("="*70)
    print(f"\nTry running with different parameters:")
    print(f"  - Higher p_max (25-30) for better accuracy")
    print(f"  - More sources for testing larger systems")
    print(f"  - Different evaluation points")
    print(f"\nSee USAGE_EXAMPLES.txt for more examples.")
    print("="*70 + "\n")


if __name__ == "__main__":
    # Установить seed для воспроизводимости / Set seed for reproducibility
    np.random.seed(42)
    demo()
