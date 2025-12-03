#!/usr/bin/env python3
"""
Visualization script for Fortran multipole simulation output.
Reads CSV files and creates animation.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import sys
import os


def load_frame(filename):
    """Load a single frame from CSV."""
    try:
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return data
    except:
        return None


def visualize_frames(frame_pattern='frame_*.csv', save_animation=False, output_file='animation.mp4'):
    """
    Visualize all frames matching the pattern.

    Args:
        frame_pattern: Glob pattern for frame files
        save_animation: Whether to save animation to file
        output_file: Output filename for animation
    """
    # Find all frame files
    frames = sorted(glob.glob(frame_pattern))

    if not frames:
        print(f"No frames found matching pattern: {frame_pattern}")
        return

    print(f"Found {len(frames)} frames")

    # Load first frame to get domain info
    first_data = load_frame(frames[0])
    if first_data is None:
        print("Could not load first frame")
        return

    # Determine domain bounds
    all_x = []
    all_y = []
    for f in frames:
        data = load_frame(f)
        if data is not None:
            all_x.extend(data[:, 0])
            all_y.extend(data[:, 1])

    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = min(all_y), max(all_y)

    # Add margin
    margin = 0.1
    xrange = xmax - xmin
    yrange = ymax - ymin
    xmin -= margin * xrange
    xmax += margin * xrange
    ymin -= margin * yrange
    ymax += margin * yrange

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 8))

    def update_plot(frame_num):
        """Update plot for given frame number."""
        if frame_num >= len(frames):
            return

        filename = frames[frame_num]
        data = load_frame(filename)

        if data is None:
            return

        ax.clear()

        x = data[:, 0]
        y = data[:, 1]

        if data.shape[1] >= 3:
            charges = data[:, 2]
            # Color by charge
            colors = ['red' if q > 0 else 'blue' for q in charges]
            sizes = [abs(q) * 50 + 10 for q in charges]
        else:
            colors = 'blue'
            sizes = 20

        ax.scatter(x, y, c=colors, s=sizes, alpha=0.7, edgecolors='black', linewidth=0.5)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title(f'Frame {frame_num} / {len(frames)-1}\n{os.path.basename(filename)}')

        # Add statistics
        if data.shape[1] >= 3:
            total_charge = np.sum(charges)
            ax.text(0.02, 0.98, f'N={len(x)}, Q_tot={total_charge:.2f}',
                   transform=ax.transAxes, va='top', fontsize=10,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Interactive mode
    if not save_animation:
        print("Interactive visualization (close window or Ctrl+C to exit)")
        print("Press any key to advance frame...")

        for i in range(len(frames)):
            update_plot(i)
            plt.pause(0.1)

        plt.show()

    # Save animation
    else:
        print(f"Creating animation: {output_file}")
        anim = animation.FuncAnimation(
            fig, update_plot, frames=len(frames),
            interval=50, repeat=True
        )

        # Save
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=20, metadata=dict(artist='Fortran Multipole'), bitrate=1800)
        anim.save(output_file, writer=writer)
        print(f"Animation saved to {output_file}")


def plot_energy_log(logfile='energy.log'):
    """Plot energy evolution from log file."""
    if not os.path.exists(logfile):
        print(f"Log file not found: {logfile}")
        return

    # Read log file
    data = []
    with open(logfile, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                try:
                    parts = line.split()
                    step = int(parts[0])
                    time = float(parts[1])
                    energy = float(parts[2])
                    data.append([step, time, energy])
                except:
                    pass

    if not data:
        print("No valid data in log file")
        return

    data = np.array(data)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Energy vs time
    ax1.plot(data[:, 1], data[:, 2], 'b-', linewidth=2)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Total Energy')
    ax1.set_title('Energy Evolution')
    ax1.grid(True, alpha=0.3)

    # Energy vs step
    ax2.plot(data[:, 0], data[:, 2], 'r-', linewidth=2)
    ax2.set_xlabel('Step')
    ax2.set_ylabel('Total Energy')
    ax2.set_title('Energy vs Step')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('energy_evolution.png', dpi=150)
    print("Saved energy_evolution.png")
    plt.show()


def main():
    """Main function."""
    import argparse

    parser = argparse.ArgumentParser(description='Visualize Fortran multipole simulation')
    parser.add_argument('--pattern', type=str, default='frame_*.csv',
                       help='Glob pattern for frame files')
    parser.add_argument('--save', action='store_true',
                       help='Save animation to file')
    parser.add_argument('--output', type=str, default='animation.mp4',
                       help='Output filename for animation')
    parser.add_argument('--energy', action='store_true',
                       help='Plot energy evolution')
    parser.add_argument('--energy-file', type=str, default='energy.log',
                       help='Energy log file')

    args = parser.parse_args()

    if args.energy:
        plot_energy_log(args.energy_file)
    else:
        visualize_frames(args.pattern, args.save, args.output)


if __name__ == "__main__":
    main()
