"""
Test building tilted SOA waveguide geometry in FDTD.

Verifies that a tilted ridge waveguide (SOA) can be built in the XY plane
and connected to a straight PWB.  Pure geometry test — no FDTD simulation
is run; inspect the saved .fsp in the Lumerical CAD.

Tilt model (top-down / XY-plane view)
--------------------------------------
All tilts are rotations around the global Z axis (looking from above).

  theta_chip  — SOA chip body tilt relative to +X (PWB axis)
  theta_ridge — ridge waveguide extra tilt *relative to the chip body*

Effective ridge angle:  theta_eff = theta_chip + theta_ridge

Typical values: theta_chip ≈ 7°, theta_ridge ≈ 7° → theta_eff ≈ 14°.

Geometry (looking from +Z down onto the XY plane)::

    SOA out (tilted)          PWB (straight)        SOA in (tilted)
         /                                           /
        / θ                                         / θ
       /_______ … _______________________________… /
       facet @ 0        250 um             facet @ L

Each SOA is built as 3 independent rotated rectangles — no structure
groups — to avoid Lumerical API selection / group-coordinate pitfalls.

Output file: SOA-PWB-SOA/test_tilted_structure.fsp
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
THETA_CHIP  = 7.0        # SOA chip body tilt (XY plane, around Z)
THETA_RIDGE = 7.0        # ridge extra tilt relative to chip body

# =========================================================================
# SOA waveguide dimensions (metres)
# =========================================================================
SOA_SIM_LEN = 15e-6      # length of SOA section inside the simulation domain
SOA_CORE_W  = 5e-6       # waveguide core width (Y-span, before rotation)
SOA_CORE_H  = 2e-6       # waveguide core height (Z-span)
SOA_RIDGE_W = 3e-6       # ridge width
SOA_RIDGE_H = 0.2e-6     # ridge (active) height
SOA_SUB_PAD = 2e-6       # extra substrate thickness below core (Z)

# =========================================================================
# Refractive indices @ 1550 nm (custom dielectrics — not in material DB)
# =========================================================================
N_InP     = 3.17         # InP substrate / cladding
N_InGaAsP = 3.40         # InGaAsP active region (approximate)

# =========================================================================
# PWB (straight section between the two SOA facets)
# =========================================================================
PWB_LENGTH = 250e-6
PWB_RADIUS = 1.1e-6
PWB_N_SEG  = 200

# =========================================================================
# Helpers
# =========================================================================

def _build_tilted_soa(fdtd, theta_body_deg, theta_ridge_deg,
                      x_facet, sign, prefix):
    """Build ONE tilted SOA waveguide section (substrate + core + ridge).

    All rotations are around the **global Z axis** (tilt in XY plane).

    Parameters
    ----------
    fdtd : lumapi.FDTD handle
    theta_body_deg : float
        Tilt of the chip body (substrate + core) in XY plane [degrees].
    theta_ridge_deg : float
        Tilt of the ridge in XY plane [degrees] (typically
        theta_body_deg + extra internal ridge tilt).
    x_facet : float
        Global X coordinate of the facet (where SOA meets PWB).
    sign : int
        -1 → output side (waveguide extends in **-X** from facet).
        +1 → input  side (waveguide extends in **+X** from facet).
    prefix : str
        Name prefix, e.g. ``"SOA_out"`` or ``"SOA_in"``.
    """
    theta_body  = np.deg2rad(theta_body_deg)
    theta_ridge = np.deg2rad(theta_ridge_deg)
    half = SOA_SIM_LEN / 2.0

    # ---- direction vectors in XY plane ----
    # For output (sign = -1): waveguide goes -X  with a +Y tilt
    #   direction = (-cos θ, +sin θ)
    #   rotation  = 180° - θ   (local X axis points along direction)
    #
    # For input  (sign = +1): waveguide goes +X  with a +Y tilt
    #   direction = (+cos θ, +sin θ)
    #   rotation  = +θ

    # body (substrate + core)
    if sign < 0:
        body_dir_x = -np.cos(theta_body)
        body_dir_y = +np.sin(theta_body)
        body_rot   = 180.0 - theta_body_deg
    else:
        body_dir_x = +np.cos(theta_body)
        body_dir_y = +np.sin(theta_body)
        body_rot   = theta_body_deg

    body_cx = x_facet + half * body_dir_x
    body_cy = half * body_dir_y

    # ridge (may have additional internal tilt)
    if sign < 0:
        ridge_dir_x = -np.cos(theta_ridge)
        ridge_dir_y = +np.sin(theta_ridge)
        ridge_rot   = 180.0 - theta_ridge_deg
    else:
        ridge_dir_x = +np.cos(theta_ridge)
        ridge_dir_y = +np.sin(theta_ridge)
        ridge_rot   = theta_ridge_deg

    ridge_cx = x_facet + half * ridge_dir_x
    ridge_cy = half * ridge_dir_y

    # ---- Z offsets (same for body & ridge — Z is rotation axis) ----
    z_core  = SOA_SUB_PAD / 2.0
    z_ridge = SOA_SUB_PAD / 2.0 + SOA_CORE_H / 2.0 + SOA_RIDGE_H / 2.0

    # =================================================================
    # 1. Substrate  (at body angle)
    # =================================================================
    fdtd.addrect()
    fdtd.set("name", prefix + "_sub")
    fdtd.set("material", "<Object defined dielectric>")
    fdtd.set("index", N_InP)
    fdtd.set("x", body_cx)
    fdtd.set("y", body_cy)
    fdtd.set("z", 0.0)
    fdtd.set("x span", SOA_SIM_LEN)
    fdtd.set("y span", SOA_CORE_W + 4e-6)
    fdtd.set("z span", SOA_CORE_H + SOA_SUB_PAD)
    fdtd.set("first axis", "z")
    fdtd.set("rotation 1", body_rot)

    # =================================================================
    # 2. Core / active region  (at body angle)
    # =================================================================
    fdtd.addrect()
    fdtd.set("name", prefix + "_core")
    fdtd.set("material", "<Object defined dielectric>")
    fdtd.set("index", N_InGaAsP)
    fdtd.set("x", body_cx)
    fdtd.set("y", body_cy)
    fdtd.set("z", z_core)
    fdtd.set("x span", SOA_SIM_LEN)
    fdtd.set("y span", SOA_CORE_W)
    fdtd.set("z span", SOA_CORE_H)
    fdtd.set("first axis", "z")
    fdtd.set("rotation 1", body_rot)

    # =================================================================
    # 3. Ridge  (at ridge angle, may differ from body)
    # =================================================================
    fdtd.addrect()
    fdtd.set("name", prefix + "_ridge")
    fdtd.set("material", "ridge")                # from material DB
    fdtd.set("x", ridge_cx)
    fdtd.set("y", ridge_cy)
    fdtd.set("z", z_ridge)
    fdtd.set("x span", SOA_SIM_LEN)
    fdtd.set("y span", SOA_RIDGE_W)
    fdtd.set("z span", SOA_RIDGE_H)
    fdtd.set("first axis", "z")
    fdtd.set("rotation 1", ridge_rot)


def _build_pwb(fdtd, x_start, x_end, radius, n_seg, material):
    """Place straight circular segments along +X (PWB bridge)."""
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
    fdtd = lumapi.FDTD(hidden=False)
    try:
        fdtd.deleteall()

        if MATERIAL_DB.exists():
            fdtd.importmaterialdb(str(MATERIAL_DB))

        # Effective angles
        theta_eff = THETA_CHIP + THETA_RIDGE

        print(f"THETA_CHIP  = {THETA_CHIP}°   (body: sub + core)")
        print(f"THETA_RIDGE = {THETA_RIDGE}°   (extra on top of chip)")
        print(f"theta_eff   = {theta_eff}°  (ridge total)")

        # ---- SOA output (facet at x=0, extends -X) ----
        _build_tilted_soa(
            fdtd,
            theta_body_deg=THETA_CHIP,
            theta_ridge_deg=theta_eff,
            x_facet=0.0,
            sign=-1,
            prefix="SOA_out",
        )

        # ---- SOA input (facet at x=PWB_LENGTH, extends +X) ----
        _build_tilted_soa(
            fdtd,
            theta_body_deg=THETA_CHIP,
            theta_ridge_deg=theta_eff,
            x_facet=PWB_LENGTH,
            sign=+1,
            prefix="SOA_in",
        )

        # ---- PWB ----
        _build_pwb(fdtd, 0.0, PWB_LENGTH, PWB_RADIUS, PWB_N_SEG, "Vancore B")

        fdtd.save(str(SAVE_PATH))
        print(f"\nSaved → {SAVE_PATH}")
        print(f"SOA sim length = {SOA_SIM_LEN * 1e6:.1f} μm each side")
        print(f"PWB length     = {PWB_LENGTH * 1e6:.1f} μm")
        print(f"Inspect the .fsp in Lumerical FDTD (top-down view recommended).")

    finally:
        print("\nFDTD session open — close manually or fdtd.close().")


if __name__ == "__main__":
    main()
