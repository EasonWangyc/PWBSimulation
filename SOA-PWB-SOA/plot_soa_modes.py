r"""
Plot SOA waveguide fundamental TE mode profiles for publication.

Loads the completed AR-coated simulation, extracts source mode profile
(SOA output WG mode) and plots |E|² with logarithmic colour scale using
the inferno colormap.

Usage
-----
.. code-block:: powershell

   & "D:\Program Files\Lumerical\python\python.exe" SOA-PWB-SOA\plot_soa_modes.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.colors import LogNorm

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import SOA_RESULTS_DIR, add_lumerical_api_path

add_lumerical_api_path()
import lumapi  # noqa: E402

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
FSP_PATH = Path(__file__).resolve().parent / "entire with ar.fsp"
PICTURES_DIR = SOA_RESULTS_DIR / "Pictures"
PICTURES_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Plot style
# ---------------------------------------------------------------------------
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["xtick.labelsize"] = 11
plt.rcParams["ytick.labelsize"] = 11
plt.rcParams["axes.labelsize"] = 13
plt.rcParams["axes.titlesize"] = 14

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print(f"Loading {FSP_PATH} ...")
fdtd = lumapi.FDTD()
fdtd.load(str(FSP_PATH))

source_data = fdtd.getresult("source", "mode profile")
fdtd.close()

# ---------------------------------------------------------------------------
# Compute |E|^2 intensity
# ---------------------------------------------------------------------------
E = source_data["E"]                         # (1, Ny, Nz, 1, 3)
Ex = E[0, :, :, 0, 0]                        # Re(Ex)
Ey = E[0, :, :, 0, 1]                        # Re(Ey)
Ez = E[0, :, :, 0, 2]                        # Re(Ez)
intensity = np.abs(Ex) ** 2 + np.abs(Ey) ** 2 + np.abs(Ez) ** 2

y = np.squeeze(source_data["y"]) * 1e6        # μm
z = np.squeeze(source_data["z"]) * 1e6        # μm

# Normalise to unity peak for visualisation
intensity_norm = intensity / np.nanmax(intensity)

# ---------------------------------------------------------------------------
# Plot: both SOA modes side by side
# ---------------------------------------------------------------------------
fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(12, 5))

# --- shared plotting helper ---
def _plot_mode(ax, yy, zz, data, title):
    vmin = max(np.nanmin(data[data > 0]), 1e-5)  # floor for log scale
    vmax = np.nanmax(data)
    norm = LogNorm(vmin=vmin, vmax=vmax)

    Y, Z = np.meshgrid(yy, zz, indexing="ij")
    pcm = ax.pcolormesh(Y, Z, data, cmap="inferno", shading="auto", norm=norm)
    cbar = fig.colorbar(pcm, ax=ax, pad=0.02,
                        format=ticker.LogFormatterSciNotation())
    cbar.set_label("|E|$^{2}$ (norm.)", fontsize=11)
    ax.set_xlabel("Y (μm)")
    ax.set_ylabel("Z (μm)")
    ax.set_title(title)
    ax.set_aspect("equal")

# Left: SOA output waveguide mode  (source position)
_plot_mode(ax_left, y, z, intensity_norm,
           "SOA Output WG\nFundamental TE Mode")

# Right: SOA input waveguide mode  (identical cross-section → same data)
_plot_mode(ax_right, y, z, intensity_norm,
           "SOA Input WG\nFundamental TE Mode")

fig.suptitle("SOA Ridge Waveguide — Fundamental TE Mode Profile  (@ 1550 nm)",
             fontsize=15, fontweight="bold")
plt.tight_layout()

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out_path = PICTURES_DIR / "SOA_fundamental_modes.png"
fig.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Saved → {out_path}")
plt.show()
