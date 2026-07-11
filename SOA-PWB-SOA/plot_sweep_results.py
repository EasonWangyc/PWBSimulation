"""
Plot sweep results: 6 panels showing taper optimisation.

Plots 1-4 (width sweeps):  Overlap (%)  vs  r_in, r_in_2, r_out, r_out_2
  - r_in / r_out      -> blue dots + blue cubic-spline fit
  - r_in_2 / r_out_2  -> red  dots + red  cubic-spline fit
  For each slice, the *other* parameter is held at the value that
  maximises T_forward.

Plots 5-6 (length sweeps):  Loss (dB) vs  taper1_length, taper2_length.

Output: results/Pictures/sweep_summary.png
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

RESULTS_DIR = Path(__file__).resolve().parent / "results"
OUTPUT      = RESULTS_DIR / "Pictures" / "sweep_summary.png"
OUTPUT.parent.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
w = pd.read_csv(RESULTS_DIR / "sweep_width_results.csv")
L = pd.read_csv(RESULTS_DIR / "sweep_length_results.csv")

t1 = w[w.direction == "taper1"].copy()
t2 = w[w.direction == "taper2"].copy()
L1 = L[L.direction == "taper1"].dropna(subset=["T_forward"])
L2 = L[L.direction == "taper2"].dropna(subset=["T_forward"])

# ---------------------------------------------------------------------------
# Helper: slice a 2D grid along one axis, fixing the other at best T
# ---------------------------------------------------------------------------
def _slice_1d(df, x_col, fix_col):
    """Return (x_vals_um, y_overlap_pct, best_fix_um)."""
    grouped = df.groupby(fix_col)["T_forward"].mean()
    best_fix = grouped.idxmax()
    sub = df[df[fix_col] == best_fix].sort_values(x_col)
    return (sub[x_col].values,
            sub["T_forward"].values * 100,
            best_fix)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
plt.rcParams["font.family"]  = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"
FONT = 14

fig, axes = plt.subplots(2, 3, figsize=(21, 12))
ax = axes.flatten()

# ---- 1. r_in (blue) — fix at best r_in_2 ----
x, y, best = _slice_1d(t1, "r1_um", "r2_um")
xs = np.linspace(x.min(), x.max(), 200)
cs = CubicSpline(x, y)
ax[0].plot(xs, cs(xs), color="blue", linewidth=1.5)
ax[0].scatter(x, y, color="blue", s=40, zorder=5)
ax[0].set_xlabel(r"$r_{\rm in}$ (μm)", fontsize=FONT)
ax[0].set_ylabel("Overlap (%)", fontsize=FONT)
ax[0].set_title(f"Taper-1:  $r_{{\\rm in}}$  "
                f"($r_{{\\rm in2}}={best:.1f}$ μm)", fontsize=FONT)
ax[0].tick_params(labelsize=FONT - 2)
ax[0].grid(True, alpha=0.3)

# ---- 2. r_in_2 (red) — fix at best r_in ----
x, y, best = _slice_1d(t1, "r2_um", "r1_um")
xs = np.linspace(x.min(), x.max(), 200)
cs = CubicSpline(x, y)
ax[1].plot(xs, cs(xs), color="red", linewidth=1.5)
ax[1].scatter(x, y, color="red", s=40, zorder=5)
ax[1].set_xlabel(r"$r_{\rm in2}$ (μm)", fontsize=FONT)
ax[1].set_ylabel("Overlap (%)", fontsize=FONT)
ax[1].set_title(f"Taper-1:  $r_{{\\rm in2}}$  "
                f"($r_{{\\rm in}}={best:.1f}$ μm)", fontsize=FONT)
ax[1].tick_params(labelsize=FONT - 2)
ax[1].grid(True, alpha=0.3)

# ---- 3. r_out (blue) — fix at best r_out_2 ----
x, y, best = _slice_1d(t2, "r1_um", "r2_um")
xs = np.linspace(x.min(), x.max(), 200)
cs = CubicSpline(x, y)
ax[2].plot(xs, cs(xs), color="blue", linewidth=1.5)
ax[2].scatter(x, y, color="blue", s=40, zorder=5)
ax[2].set_xlabel(r"$r_{\rm out}$ (μm)", fontsize=FONT)
ax[2].set_ylabel("Overlap (%)", fontsize=FONT)
ax[2].set_title(f"Taper-2:  $r_{{\\rm out}}$  "
                f"($r_{{\\rm out2}}={best:.1f}$ μm)", fontsize=FONT)
ax[2].tick_params(labelsize=FONT - 2)
ax[2].grid(True, alpha=0.3)

# ---- 4. r_out_2 (red) — fix at best r_out ----
x, y, best = _slice_1d(t2, "r2_um", "r1_um")
xs = np.linspace(x.min(), x.max(), 200)
cs = CubicSpline(x, y)
ax[3].plot(xs, cs(xs), color="red", linewidth=1.5)
ax[3].scatter(x, y, color="red", s=40, zorder=5)
ax[3].set_xlabel(r"$r_{\rm out2}$ (μm)", fontsize=FONT)
ax[3].set_ylabel("Overlap (%)", fontsize=FONT)
ax[3].set_title(f"Taper-2:  $r_{{\\rm out2}}$  "
                f"($r_{{\\rm out}}={best:.1f}$ μm)", fontsize=FONT)
ax[3].tick_params(labelsize=FONT - 2)
ax[3].grid(True, alpha=0.3)

# ---- 5. taper1 length -> Loss (dB) ----
x1 = L1["length_um"].values
y1 = -10 * np.log10(L1["T_forward"].values.astype(float))
xs1 = np.linspace(x1.min(), x1.max(), 200)
cs1 = CubicSpline(x1, y1)
ax[4].plot(xs1, cs1(xs1), color="green", linewidth=1.5)
ax[4].scatter(x1, y1, color="green", s=40, zorder=5)
ax[4].set_xlabel("Taper-1 Length (μm)", fontsize=FONT)
ax[4].set_ylabel("Loss (dB)", fontsize=FONT)
ax[4].set_title("Taper-1 Length", fontsize=FONT)
ax[4].tick_params(labelsize=FONT - 2)
ax[4].grid(True, alpha=0.3)

# ---- 6. taper2 length -> Loss (dB) ----
x2 = L2["length_um"].values
y2 = -10 * np.log10(L2["T_forward"].values.astype(float))
xs2 = np.linspace(x2.min(), x2.max(), 200)
cs2 = CubicSpline(x2, y2)
ax[5].plot(xs2, cs2(xs2), color="green", linewidth=1.5)
ax[5].scatter(x2, y2, color="green", s=40, zorder=5)
ax[5].set_xlabel("Taper-2 Length (μm)", fontsize=FONT)
ax[5].set_ylabel("Loss (dB)", fontsize=FONT)
ax[5].set_title("Taper-2 Length", fontsize=FONT)
ax[5].tick_params(labelsize=FONT - 2)
ax[5].grid(True, alpha=0.3)

# ---------------------------------------------------------------------------
fig.tight_layout()
fig.savefig(OUTPUT, dpi=200, bbox_inches="tight")
print(f"Saved -> {OUTPUT}")
plt.show()
