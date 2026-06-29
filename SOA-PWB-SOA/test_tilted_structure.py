"""
Test building tilted SOA waveguide geometry in FDTD.

Verifies whether a tilted ridge waveguide (SOA) can be constructed and
connected to a straight PWB.  This is a pure geometry test — no FDTD
simulation is run; inspect the saved .fsp in the Lumerical CAD.

Tilt model
----------
Two independent angles (both rotation around the global Y axis):

  theta_chip  — SOA chip body tilt relative to the PWB axis (+X)
  theta_ridge — ridge waveguide tilt *relative to the chip body*

The effective ridge angle w.r.t. the PWB axis is::

    theta_eff = theta_chip + theta_ridge

Typical values (from literature / qa.md discussion):
  - theta_chip  ≈ 7°   (facet anti-reflection cleave)
  - theta_ridge ≈ 7°   (additional ridge tilt within the SOA chip)

The tilted SOA sections are built using Lumerical structure groups:
each SOA's substrate, core and ridge are placed in a group at local
coordinates, then the group is rotated and positioned so the output
facet lands at the correct x coordinate.

Output file
-----------
SOA-PWB-SOA/test_tilted_structure.fsp
"""

import sys
from pathlib import Path
import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import SOA_DIR, SOA_BASE_FSP, MATERIAL_DB, add_lumerical_api_path

add_lumerical_api_path()
import lumapi  # noqa: E402

SAVE_PATH = SOA_DIR / "test_tilted_structure.fsp"

# =========================================================================
# Configurable angles (degrees) — edit these to explore different tilts
# =========================================================================
THETA_CHIP = 7.0        # SOA chip body tilt (XZ plane, around Y)
THETA_RIDGE = 7.0       # ridge waveguide extra tilt relative to chip body

# =========================================================================
# SOA waveguide dimensions (metres)
# =========================================================================
SOA_SIM_LEN = 15e-6     # length of SOA section inside the simulation domain
SOA_CORE_W  = 5e-6      # waveguide core width  (Y-span, before rotation)
SOA_CORE_H  = 2e-6      # waveguide core height (Z-span, before rotation)
SOA_RIDGE_W = 3e-6      # ridge width
SOA_RIDGE_H = 0.2e-6    # ridge (active) height
SOA_SUB_PAD = 2e-6      # extra substrate thickness below core

# =========================================================================
# Refractive indices @ 1550 nm (custom dielectrics — not in material DB)
# =========================================================================
N_InP      = 3.17       # InP substrate / cladding
N_InGaAsP  = 3.40       # InGaAsP active region (approximate)

# =========================================================================
# PWB (straight section between the two SOA facets)
# =========================================================================
PWB_LENGTH = 250e-6
PWB_RADIUS = 1.1e-6
PWB_N_SEG   = 200

# =========================================================================
# Helpers
# =========================================================================

def _rot_y(theta_deg):
    """Return (cos, sin) for a rotation around +Y by *theta_deg* degrees."""
    th = np.deg2rad(theta_deg)
    return np.cos(th), np.sin(th)


def _build_soa_group(fdtd, theta_deg, x_facet, sign, prefix):
    """Build ONE tilted SOA section (substrate + core + ridge) inside a group.

    Parameters
    ----------
    fdtd : lumapi.FDTD handle
    theta_deg : float
        Total effective tilt angle of this SOA section (degrees, around Y).
    x_facet : float
        Global X coordinate where the SOA facet meets the PWB.
    sign : int
        -1 for the output SOA (waveguide extends in the **-X** direction
        from the facet); +1 for the input SOA (extends in **+X**).
    prefix : str
        Name prefix for the group and its child structures.
    """
    cos_th, sin_th = _rot_y(theta_deg)
    half = SOA_SIM_LEN / 2.0

    # Group centre in global coordinates — offset by half the SOA length
    # *along the tilted axis* away from the facet.
    gx = x_facet + sign * half * cos_th
    gy = 0.0
    gz = sign * half * sin_th

    # ---- create group ----
    fdtd.addgroup()
    fdtd.set("name", prefix)

    # ---- substrate (InP, wider & taller than core) ----
    fdtd.addrect()
    fdtd.set("name", prefix + "_sub")
    fdtd.set("material", "<Object defined dielectric>")
    fdtd.set("index", N_InP)
    fdtd.set("x", 0.0)
    fdtd.set("y", 0.0)
    fdtd.set("z", 0.0)
    fdtd.set("x span", SOA_SIM_LEN)
    fdtd.set("y span", SOA_CORE_W + 4e-6)
    fdtd.set("z span", SOA_CORE_H + SOA_SUB_PAD)
    fdtd.addtogroup(prefix)

    # ---- core (InGaAsP active region) ----
    fdtd.addrect()
    fdtd.set("name", prefix + "_core")
    fdtd.set("material", "<Object defined dielectric>")
    fdtd.set("index", N_InGaAsP)
    fdtd.set("x", 0.0)
    fdtd.set("y", 0.0)
    fdtd.set("z", SOA_SUB_PAD / 2.0)
    fdtd.set("x span", SOA_SIM_LEN)
    fdtd.set("y span", SOA_CORE_W)
    fdtd.set("z span", SOA_CORE_H)
    fdtd.addtogroup(prefix)

    # ---- ridge (InP cladding on top of active region) ----
    fdtd.addrect()
    fdtd.set("name", prefix + "_ridge")
    fdtd.set("material", "<Object defined dielectric>")
    fdtd.set("index", N_InP)
    fdtd.set("x", 0.0)
    fdtd.set("y", 0.0)
    fdtd.set("z", SOA_SUB_PAD / 2.0 + SOA_CORE_H / 2.0 + SOA_RIDGE_H / 2.0)
    fdtd.set("x span", SOA_SIM_LEN)
    fdtd.set("y span", SOA_RIDGE_W)
    fdtd.set("z span", SOA_RIDGE_H)
    fdtd.addtogroup(prefix)

    # ---- position & rotate the group ----
    fdtd.set("x", gx)
    fdtd.set("y", gy)
    fdtd.set("z", gz)
    fdtd.set("first axis", "y")
    fdtd.set("rotation 1", theta_deg)


def _build_pwb(fdtd, x_start, x_end, radius, n_seg, material):
    """Place straight circular segments along +X."""
    span = x_end - x_start
    for i in range(n_seg):
        t0 = i / n_seg
        t1 = (i + 1) / n_seg
        x0 = x_start + t0 * span
        x1 = x_start + t1 * span
        cx = (x0 + x1) / 2.0
        seg_len = x1 - x0

        fdtd.addcircle()
        fdtd.set("name", f"PWB_{i}")
        fdtd.set("material", material)
        fdtd.set("make ellipsoid", 0)
        fdtd.set("x", cx)
        fdtd.set("y", 0.0)
        fdtd.set("z", 0.0)
        fdtd.set("radius", radius)
        fdtd.set("z span", seg_len)
        fdtd.set("first axis", "y")
        fdtd.set("rotation 1", 90)


# =========================================================================
# Main
# =========================================================================

def main():
    fdtd = lumapi.FDTD(hidden=False)   # visible so user can inspect
    try:
        fdtd.deleteall()

        # Import materials so "InP" / "Vancore B" etc. are resolvable
        if MATERIAL_DB.exists():
            fdtd.importmaterialdb(str(MATERIAL_DB))

        # ---- optionally load the existing SOA .fsp for reference ----
        if SOA_BASE_FSP.exists():
            print(f"Loading reference SOA from: {SOA_BASE_FSP}")
            # We don't load it here because we're building SOA from scratch;
            # uncomment the line below to overlay the reference .fsp:
            # fdtd.load(str(SOA_BASE_FSP))

        # Compute effective ridge angles
        theta_out = THETA_CHIP + THETA_RIDGE   # output SOA
        theta_in  = THETA_CHIP + THETA_RIDGE   # input SOA (same sign → symmetric)

        # ---- SOA output section (facet at x = 0, extends -X) ----
        _build_soa_group(fdtd, theta_out, x_facet=0.0, sign=-1, prefix="SOA_out")

        # ---- SOA input section (facet at x = PWB_LENGTH, extends +X) ----
        _build_soa_group(fdtd, theta_in, x_facet=PWB_LENGTH, sign=+1, prefix="SOA_in")

        # ---- PWB (straight, along +X) ----
        _build_pwb(fdtd, 0.0, PWB_LENGTH, PWB_RADIUS, PWB_N_SEG, "Vancore B")

        # ---- save ----
        fdtd.save(str(SAVE_PATH))
        print(f"\nStructure saved → {SAVE_PATH}")
        print(f"  theta_chip  = {THETA_CHIP}°")
        print(f"  theta_ridge = {THETA_RIDGE}°")
        print(f"  theta_eff   = {theta_out}°  (output & input)")
        print(f"  SOA sim length = {SOA_SIM_LEN * 1e6:.1f} μm per side")
        print(f"  PWB length     = {PWB_LENGTH * 1e6:.1f} μm")
        print(f"\nOpen the .fsp in Lumerical FDTD to inspect the geometry.")

    finally:
        # Keep open for visual inspection; close with fdtd.close() if desired
        print("\nFDTD session still open — close manually or call fdtd.close().")


if __name__ == "__main__":
    main()
