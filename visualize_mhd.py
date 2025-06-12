#!/usr/bin/env python3
"""
Plot MHD solver results exactly like Figure 4 from Dedner et al. 2002
Shows: ρ, u_y, B_x, B_y
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import os

def load_final_solution(output_dir="Result"):
    """Load the final solution file"""
    files = sorted(glob.glob(os.path.join(output_dir, "out_state_*.csv")))
    if not files:
        print(f"No output files found in {output_dir}")
        return None
    
    # Load the last file (final time)
    final_file = files[-1]
    print(f"Loading final solution from: {final_file}")
    return pd.read_csv(final_file)

def plot_dedner_figure4_exact(output_dir="Result", save_name="figure4_reproduction.png"):
    """Create a plot exactly like Figure 4 in Dedner et al. 2002"""
    
    # Load data
    data = load_final_solution(output_dir)
    if data is None:
        return
    
    # Create figure with 2x2 subplots matching paper layout
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    
    # Common plot settings
    line_color = 'blue'
    line_width = 2
    
    # Top left: Density (ρ)
    ax = axes[0, 0]
    ax.plot(data['x'], data['rho'], color=line_color, linewidth=line_width)
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(0.5, 4.0)  # Expected range from paper
    ax.set_xlabel('x')
    ax.set_ylabel(r'$\rho$', fontsize=12)
    ax.grid(True, alpha=0.3)
    
    # Top right: y-velocity (u_y)
    ax = axes[0, 1]
    ax.plot(data['x'], data['v'], color=line_color, linewidth=line_width)  # 'v' is u_y
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(-1.0, 1.0)  # Expected range from paper
    ax.set_xlabel('x')
    ax.set_ylabel(r'$u_y$', fontsize=12)
    ax.grid(True, alpha=0.3)
    
    # Bottom left: x-magnetic field (B_x)
    ax = axes[1, 0]
    ax.plot(data['x'], data['bx'], color=line_color, linewidth=line_width)
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(1.0, 2.0)  # Adjust based on your Bx values
    ax.set_xlabel('x')
    ax.set_ylabel(r'$B_x$', fontsize=12)
    ax.grid(True, alpha=0.3)
    # Add reference line for initial Bx
    Bx_init = 5.0 / np.sqrt(4.0 * np.pi)
    ax.axhline(y=Bx_init, color='red', linestyle='--', alpha=0.5, linewidth=1)
    
    # Bottom right: y-magnetic field (B_y)
    ax = axes[1, 1]
    ax.plot(data['x'], data['by'], color=line_color, linewidth=line_width)
    ax.set_xlim(-0.5, 0.5)
    # Set ylim around the expected constant value
    By_init = 5.0 / np.sqrt(4.0 * np.pi)
    ax.set_ylim(0.5, 6.0)  # Adjust to see the actual range
    ax.set_xlabel('x')
    ax.set_ylabel(r'$B_y$', fontsize=12)
    ax.grid(True, alpha=0.3)
    # Add reference line for constant By
    ax.axhline(y=By_init, color='red', linestyle='--', alpha=0.5, linewidth=1,
               label=f'Initial: {By_init:.3f}')
    ax.legend(loc='best', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(save_name, dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print diagnostics
    print(f"\nSolution at t = 0.2:")
    print(f"ρ: [{data['rho'].min():.3f}, {data['rho'].max():.3f}]")
    print(f"u_y: [{data['v'].min():.3f}, {data['v'].max():.3f}]")
    print(f"B_x: [{data['bx'].min():.3f}, {data['bx'].max():.3f}] (initial: {Bx_init:.3f})")
    print(f"B_y: [{data['by'].min():.3f}, {data['by'].max():.3f}] (initial: {By_init:.3f})")
    
    # Check conservation
    print(f"\nConservation check:")
    print(f"B_x should be constant: variation = {data['bx'].max() - data['bx'].min():.6f}")
    print(f"B_y should be constant: variation = {data['by'].max() - data['by'].min():.6f}")

def plot_reference_comparison(output_dir="Result"):
    """Plot with reference lines showing expected behavior"""
    
    data = load_final_solution(output_dir)
    if data is None:
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Initial values
    Bx_init = 5.0 / np.sqrt(4.0 * np.pi)
    By_init = 5.0 / np.sqrt(4.0 * np.pi)
    
    # Density
    ax = axes[0, 0]
    ax.plot(data['x'], data['rho'], 'b-', linewidth=2, label='Numerical')
    ax.axhline(y=1.0, color='gray', linestyle=':', alpha=0.5, label='Initial')
    ax.set_xlim(-0.5, 0.5)
    ax.set_xlabel('x')
    ax.set_ylabel(r'$\rho$')
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_title('Density')
    
    # u_y (should be zero initially)
    ax = axes[0, 1]
    ax.plot(data['x'], data['v'], 'b-', linewidth=2, label='Numerical')
    ax.axhline(y=0.0, color='gray', linestyle=':', alpha=0.5, label='Initial')
    ax.set_xlim(-0.5, 0.5)
    ax.set_xlabel('x')
    ax.set_ylabel(r'$u_y$')
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_title('y-Velocity')
    
    # B_x (should remain constant)
    ax = axes[1, 0]
    ax.plot(data['x'], data['bx'], 'b-', linewidth=2, label='Numerical')
    ax.axhline(y=Bx_init, color='red', linestyle='--', label=f'Expected: {Bx_init:.3f}')
    ax.set_xlim(-0.5, 0.5)
    ax.set_xlabel('x')
    ax.set_ylabel(r'$B_x$')
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_title('x-Magnetic Field')
    
    # B_y (should remain constant)
    ax = axes[1, 1]
    ax.plot(data['x'], data['by'], 'b-', linewidth=2, label='Numerical')
    ax.axhline(y=By_init, color='red', linestyle='--', label=f'Expected: {By_init:.3f}')
    ax.set_xlim(-0.5, 0.5)
    ax.set_xlabel('x')
    ax.set_ylabel(r'$B_y$')
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_title('y-Magnetic Field')
    
    plt.suptitle('1D MHD Riemann Problem - Comparison with Expected Values', fontsize=14)
    plt.tight_layout()
    plt.savefig('reference_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()

def check_dedner_values():
    """Print expected values from Dedner Figure 4"""
    print("\nExpected values from Dedner et al. Figure 4:")
    print("="*50)
    print("ρ: Should show shock-rarefaction structure")
    print("   - Left state: 1.0")
    print("   - Peak value: ~3.8")
    print("   - Complex wave structure")
    print("\nu_y: Should develop from initially zero")
    print("   - Range approximately [-0.5, 0.3]")
    print("\nB_x: Should remain constant at", 5.0/np.sqrt(4*np.pi))
    print("\nB_y: Should remain constant at", 5.0/np.sqrt(4*np.pi))
    print("\nNote: If B_x and B_y are not constant, there's a")
    print("fundamental issue with the flux calculation.")

if __name__ == "__main__":
    # Create plots
    plot_dedner_figure4_exact()
    plot_reference_comparison()
    check_dedner_values()
    
    print("\nPlots saved as:")
    print("  - figure4_reproduction.png (matches paper layout)")
    print("  - reference_comparison.png (with expected values)")
