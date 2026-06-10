"""
analyze_curvature.py — Analyze local curvature radius of complex (_3) PWB centerlines.

Computes the local curvature radius ρ(s) = 1/κ(s) along the Bezier bending path
and plots centerline shape + curvature profile. No Lumerical needed.

Uses bend-only path (no taper1 straight section) for cleaner curvature analysis.

Usage:
    python PD-PWB-SMF/analyze_curvature.py
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))       # project root → config/
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # PD-PWB-SMF/ → pwb_core

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from pwb_core import PWBParameters, generate_pwb_bend_path, path_arc_length


# ---------------------------------------------------------------------------
# Curvature computation with smoothing
# ---------------------------------------------------------------------------

def gaussian_smooth(data, sigma=2.0):
    """1D Gaussian smoothing via convolution (no scipy dependency)."""
    data = np.asarray(data, dtype=float)
    # kernel size: ±3σ
    radius = int(np.ceil(3 * sigma))
    x = np.arange(-radius, radius + 1)
    kernel = np.exp(-0.5 * (x / sigma) ** 2)
    kernel /= kernel.sum()
    # reflect edges to avoid boundary artifacts
    padded = np.pad(data, radius, mode='reflect')
    smoothed = np.convolve(padded, kernel, mode='valid')
    # trim to original length (may be off by 1 due to even/odd kernel)
    return smoothed[:len(data)]


def compute_curvature_radius(path, smooth_sigma=2.0):
    """Compute local curvature radius along a 2D x-z path.

    Applies Gaussian smoothing before differentiation to suppress
    numerical noise from discrete sampling and C1-only Bezier junctions.

    Returns (s, rho) where s is cumulative arc length and rho is
    curvature radius, both in meters. rho = 1/kappa.
    """
    path = np.asarray(path, dtype=float)
    x_raw, z_raw = path[:, 0], path[:, 2]

    # Smooth coordinates to reduce derivative noise
    x = gaussian_smooth(x_raw, sigma=smooth_sigma)
    z = gaussian_smooth(z_raw, sigma=smooth_sigma)

    # discrete derivatives
    dx = np.gradient(x)
    dz = np.gradient(z)
    ddx = np.gradient(dx)
    ddz = np.gradient(dz)

    # curvature kappa = |x'z'' - z'x''| / (x'^2 + z'^2)^(3/2)
    numerator = np.abs(dx * ddz - dz * ddx)
    denominator = np.clip((dx**2 + dz**2) ** 1.5, 1e-30, None)
    kappa = numerator / denominator

    rho = 1.0 / np.clip(kappa, 1e-30, None)
    s = path_arc_length(path)

    return s, rho


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

RHO_DISPLAY_MAX = 200.0   # cap curvature radius display at 200 μm


def plot_single_case(params, label=None, ax_path=None, ax_curv=None, color=None):
    """Compute curvature and add plots. Returns min(rho) in μm."""
    path = generate_pwb_bend_path(params)
    s, rho = compute_curvature_radius(path)

    s_um = s * 1e6
    rho_um = np.clip(rho * 1e6, 0, RHO_DISPLAY_MAX)  # cap display
    x_um = path[:, 0] * 1e6
    z_um = path[:, 2] * 1e6

    if ax_path is not None:
        ax_path.plot(x_um, z_um, color=color, linewidth=1.5, label=label or "")
    if ax_curv is not None:
        ax_curv.plot(s_um, rho_um, color=color, linewidth=1.5)

    return np.nanmin(rho_um)


def plot_curvature_comparison(cases, title="Curvature Radius Comparison"):
    """Compare multiple parameter sets side by side.

    cases: list of (params, label) tuples
    """
    colors = cm.tab10(np.linspace(0, 1, len(cases)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    min_rhos = []
    for (params, label), color in zip(cases, colors):
        min_rho = plot_single_case(params, label, ax_path=ax1, ax_curv=ax2, color=color)
        min_rhos.append(min_rho)
        print(f"  {label:35s}  min ρ = {min_rho:.1f} μm")

    # centerline plot
    ax1.set_xlabel("X (μm)")
    ax1.set_ylabel("Z (μm)")
    ax1.set_title("Centerline Shape (X-Z plane)")
    ax1.set_aspect("equal")
    ax1.legend(fontsize=8, loc="lower left")
    ax1.invert_yaxis()
    ax1.grid(True, alpha=0.3)

    # curvature plot — focus on the low-ρ range
    ax2.set_xlabel("Arc length s (μm)")
    ax2.set_ylabel("Curvature radius ρ (μm)")
    ax2.set_title("Local Curvature Radius along Path")
    ax2.set_ylim(0, RHO_DISPLAY_MAX)
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    fig.suptitle(title, fontsize=13, fontweight='bold')
    plt.tight_layout()

    safe_title = title.replace(" ", "_").replace(":", "").replace("/", "_")
    out_dir = Path(__file__).resolve().parent.parent / "results" / "curvature_analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"curvature_{safe_title}.png"
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    print(f"\nSaved: {out_path}")
    plt.close(fig)

    return min_rhos


# ---------------------------------------------------------------------------
# Analysis routines
# ---------------------------------------------------------------------------

BASE = PWBParameters()
BASE.L = 250e-6; BASE.h = 100e-6; BASE.l1 = 100e-6
BASE.bend_shape = 0.5; BASE.drop_shape = 0.7


def _make_params(**overrides):
    p = PWBParameters()
    for attr in ['L', 'h', 'l1', 'bend_shape', 'drop_shape', 'complex_segments']:
        setattr(p, attr, getattr(BASE, attr))
    for k, v in overrides.items():
        setattr(p, k, v)
    return p


def analyze_bend_lift():
    print("=" * 60)
    print("Varying bend_lift (arch height)")
    print("=" * 60)
    cases = [(_make_params(bend_lift=v), f"bend_lift = {v:.2f}")
             for v in [0.0, 0.05, 0.12, 0.20, 0.35]]
    plot_curvature_comparison(cases, "Effect of bend_lift on Curvature Radius")


def analyze_arch_position():
    print("=" * 60)
    print("Varying arch_position (peak position)")
    print("=" * 60)
    cases = [(_make_params(arch_position=v), f"arch_position = {v:.2f}")
             for v in [0.25, 0.35, 0.45, 0.55, 0.65, 0.75]]
    plot_curvature_comparison(cases, "Effect of arch_position on Curvature Radius")


def analyze_bend_shape():
    print("=" * 60)
    print("Varying bend_shape (arch transition smoothness)")
    print("=" * 60)
    cases = [(_make_params(bend_shape=v), f"bend_shape = {v:.2f}")
             for v in [0.1, 0.3, 0.5, 0.7, 0.9]]
    plot_curvature_comparison(cases, "Effect of bend_shape on Curvature Radius")


def analyze_drop_shape():
    print("=" * 60)
    print("Varying drop_shape (drop transition smoothness)")
    print("=" * 60)
    cases = [(_make_params(drop_shape=v), f"drop_shape = {v:.2f}")
             for v in [0.1, 0.3, 0.5, 0.7, 0.9]]
    plot_curvature_comparison(cases, "Effect of drop_shape on Curvature Radius")


def analyze_default():
    """Single detailed analysis with default parameters."""
    print("=" * 60)
    print("Default _3 bend-only centerline detailed analysis")
    print("=" * 60)
    params = _make_params()
    path = generate_pwb_bend_path(params)
    s, rho = compute_curvature_radius(path)

    x_um = path[:, 0] * 1e6
    z_um = path[:, 2] * 1e6
    s_um = s * 1e6
    rho_um = rho * 1e6
    rho_clip = np.clip(rho_um, 0, RHO_DISPLAY_MAX)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # centerline
    ax = axes[0, 0]
    ax.plot(x_um, z_um, 'b-', linewidth=1.5)
    ax.set_xlabel("X (μm)"); ax.set_ylabel("Z (μm)")
    ax.set_title("Centerline (X-Z)")
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.grid(True, alpha=0.3)

    # curvature radius vs arc length
    ax = axes[0, 1]
    ax.plot(s_um, rho_clip, 'r-', linewidth=1.5)
    ax.set_xlabel("Arc length s (μm)"); ax.set_ylabel("ρ (μm)")
    ax.set_title("Curvature Radius vs Arc Length")
    ax.set_ylim(0, RHO_DISPLAY_MAX)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=np.nanmin(rho_clip), color='gray', linestyle='--', alpha=0.5,
               label=f"min ρ = {np.nanmin(rho_clip):.1f} μm")
    ax.legend()

    # curvature radius color-coded on centerline
    ax = axes[1, 0]
    vmax = min(np.nanmax(rho_um), RHO_DISPLAY_MAX)
    norm = plt.Normalize(np.nanmin(rho_um), vmax)
    points = ax.scatter(x_um, z_um, c=rho_um, cmap='plasma',
                        s=3, norm=norm, vmin=0, vmax=vmax)
    plt.colorbar(points, ax=ax, label="ρ (μm)")
    ax.set_xlabel("X (μm)"); ax.set_ylabel("Z (μm)")
    ax.set_title("Curvature Radius along Path")
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.grid(True, alpha=0.3)

    # histogram
    ax = axes[1, 1]
    finite = rho_um[np.isfinite(rho_um)]
    clipped = finite[finite < RHO_DISPLAY_MAX]
    ax.hist(clipped, bins=50, color='steelblue', edgecolor='white', alpha=0.8)
    ax.axvline(x=np.nanmin(rho_um), color='r', linestyle='--',
               label=f"min = {np.nanmin(rho_um):.1f} μm")
    ax.axvline(x=np.nanmedian(rho_um), color='orange', linestyle='--',
               label=f"median = {np.nanmedian(rho_um):.1f} μm")
    ax.set_xlabel("ρ (μm)"); ax.set_ylabel("Count")
    ax.set_xlim(0, RHO_DISPLAY_MAX)
    ax.set_title("Distribution of Curvature Radius")
    ax.legend()

    fig.suptitle("Default _3 Bend-only Centerline Analysis", fontweight='bold')
    plt.tight_layout()

    out_dir = Path(__file__).resolve().parent.parent / "results" / "curvature_analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "curvature_default.png"
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    print(f"\nSaved: {out_path}")
    print(f"  Min curvature radius: {np.nanmin(rho_um):.1f} μm")
    print(f"  Median curvature radius: {np.nanmedian(rho_um):.1f} μm")
    print(f"  Max curvature radius (clipped): {np.nanmax(clipped):.1f} μm")
    plt.close(fig)


if __name__ == "__main__":
    print("PD-PWB-SMF Complex Path Curvature Analysis\n")

    analyze_default()
    analyze_bend_lift()
    analyze_arch_position()
    analyze_bend_shape()
    analyze_drop_shape()
