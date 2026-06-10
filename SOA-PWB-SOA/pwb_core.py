"""
Core helpers for SOA-PWB-SOA simulations.

Models the on-chip SOA output → PWB → external SOA input coupling path.
Structure is a straight in-line configuration:
  SOA output waveguide → Taper-1 (mode expansion) → PWB straight
  → Taper-2 (mode compression) → external SOA input waveguide

This module does not import lumapi directly; callers create and pass FDTD/MODE handles.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import MATERIAL_DB, SOA_RESULTS_DIR

PICTURES_DIR = SOA_RESULTS_DIR / "Pictures"


class SOAPWBParams:
    """Parameters for SOA-PWB-SOA simulation.

    All geometric dimensions in SI units (meters). Use X * 1e-6 for microns.
    """

    def __init__(self):
        # ---- Wavelength ----
        self.wavelength = 1.55e-6          # center wavelength [m]

        # ---- On-chip SOA output waveguide (InP ridge) ----
        self.wg1_width = 2.0e-6            # waveguide width [m]
        self.wg1_height = 1.5e-6           # waveguide height [m]
        self.wg1_active_thickness = 0.3e-6 # active region thickness [m]
        self.wg1_length = 10e-6            # SOA output segment length [m]

        # ---- Taper-1: SOA output → PWB (mode expansion) ----
        self.taper1_w_in = 2.0e-6          # input width (match SOA output) [m]
        self.taper1_w_out = 4.0e-6         # output width (match PWB) [m]
        self.taper1_h_in = 1.5e-6          # input height [m]
        self.taper1_h_out = 4.0e-6         # output height [m]
        self.taper1_length = 30e-6         # taper length [m]

        # ---- PWB polymer waveguide (middle straight section) ----
        self.pwb_width = 4.0e-6            # PWB cross-section width [m]
        self.pwb_height = 4.0e-6           # PWB cross-section height [m]
        self.pwb_length = 50e-6            # PWB straight section length [m]

        # ---- Taper-2: PWB → external SOA input (mode compression) ----
        self.taper2_w_in = 4.0e-6          # input width (match PWB) [m]
        self.taper2_w_out = 2.0e-6         # output width (match ext SOA) [m]
        self.taper2_h_in = 4.0e-6          # input height [m]
        self.taper2_h_out = 1.5e-6         # output height [m]
        self.taper2_length = 30e-6         # taper length [m]

        # ---- External SOA input waveguide ----
        self.wg2_width = 2.0e-6            # waveguide width [m]
        self.wg2_height = 1.5e-6           # waveguide height [m]
        self.wg2_length = 10e-6            # external SOA input segment length [m]

        # ---- Materials ----
        self.material_soa = "InP"          # SOA waveguide material
        self.material_pwb = "Vancore B"    # PWB polymer material
        self.material_substrate = "SiO2 (Glass) - Palik"

        # ---- Simulation settings ----
        self.mesh_accuracy = 2
        self.simulation_time = 2000e-15    # [s]


def compute_positions(params):
    """Compute cumulative x-positions for each structural section.

    Returns a dict with keys: soa_end, taper1_end, pwb_end, taper2_end, total
    """
    x_soa_end = params.wg1_length
    x_taper1_end = x_soa_end + params.taper1_length
    x_pwb_end = x_taper1_end + params.pwb_length
    x_taper2_end = x_pwb_end + params.taper2_length
    x_total = x_taper2_end + params.wg2_length

    return {
        "soa_start": 0.0,
        "soa_end": x_soa_end,
        "taper1_start": x_soa_end,
        "taper1_end": x_taper1_end,
        "pwb_start": x_taper1_end,
        "pwb_end": x_pwb_end,
        "taper2_start": x_pwb_end,
        "taper2_end": x_taper2_end,
        "ext_soa_start": x_taper2_end,
        "ext_soa_end": x_total,
        "total": x_total,
    }


def create_pwb_structure_in_fdtd(fdtd, params):
    """Build the full SOA-PWB-SOA structure in the FDTD workspace.

    Structure (all along +x axis, centered at y=0, z=0):
      1. On-chip SOA output waveguide (InP ridge, rectangular)
      2. Taper-1: linear transition SOA → PWB (mode expansion)
      3. PWB polymer waveguide (straight)
      4. Taper-2: linear transition PWB → external SOA (mode compression)
      5. External SOA input waveguide (rectangular)

    Args:
        fdtd: lumapi.FDTD handle (already opened by caller).
        params: SOAPWBParams instance.

    Returns:
        dict: positions dict from compute_positions().
    """
    pos = compute_positions(params)

    fdtd.deleteall()
    fdtd.importmaterialdb(str(MATERIAL_DB))

    # --- 1. On-chip SOA output waveguide (InP ridge) ---
    fdtd.addrect()
    fdtd.set("name", "SOA_output_ridge")
    fdtd.set("material", params.material_soa)
    fdtd.set("x", pos["soa_start"] + params.wg1_length / 2)
    fdtd.set("y", 0)
    fdtd.set("z", 0)
    fdtd.set("x span", params.wg1_length)
    fdtd.set("y span", params.wg1_width)
    fdtd.set("z span", params.wg1_height)

    # --- 2. Taper-1: SOA → PWB (linear expansion via pyramid) ---
    fdtd.addpyramid()
    fdtd.set("name", "PWB_taper1")
    fdtd.set("material", params.material_pwb)
    fdtd.set("override mesh order from material database", 1)
    fdtd.set("mesh order", 3)
    fdtd.set("override color opacity from material database", 1)
    fdtd.set("alpha", 0.5)
    fdtd.set("x", pos["taper1_start"] + params.taper1_length / 2)
    fdtd.set("y", 0)
    fdtd.set("z", 0)
    fdtd.set("x span", params.taper1_length)
    fdtd.set("y span bottom", params.taper1_w_in)
    fdtd.set("y span top", params.taper1_w_out)
    fdtd.set("z span bottom", params.taper1_h_in)
    fdtd.set("z span top", params.taper1_h_out)
    fdtd.set("first axis", "y")
    fdtd.set("rotation 1", 90)

    # --- 3. PWB polymer waveguide (straight section) ---
    fdtd.addrect()
    fdtd.set("name", "PWB_straight")
    fdtd.set("material", params.material_pwb)
    fdtd.set("override mesh order from material database", 1)
    fdtd.set("mesh order", 3)
    fdtd.set("override color opacity from material database", 1)
    fdtd.set("alpha", 0.5)
    fdtd.set("x", pos["pwb_start"] + params.pwb_length / 2)
    fdtd.set("y", 0)
    fdtd.set("z", 0)
    fdtd.set("x span", params.pwb_length)
    fdtd.set("y span", params.pwb_width)
    fdtd.set("z span", params.pwb_height)

    # --- 4. Taper-2: PWB → external SOA (linear compression via pyramid) ---
    fdtd.addpyramid()
    fdtd.set("name", "PWB_taper2")
    fdtd.set("material", params.material_pwb)
    fdtd.set("override mesh order from material database", 1)
    fdtd.set("mesh order", 3)
    fdtd.set("override color opacity from material database", 1)
    fdtd.set("alpha", 0.5)
    fdtd.set("x", pos["taper2_start"] + params.taper2_length / 2)
    fdtd.set("y", 0)
    fdtd.set("z", 0)
    fdtd.set("x span", params.taper2_length)
    fdtd.set("y span bottom", params.taper2_w_in)
    fdtd.set("y span top", params.taper2_w_out)
    fdtd.set("z span bottom", params.taper2_h_in)
    fdtd.set("z span top", params.taper2_h_out)
    fdtd.set("first axis", "y")
    fdtd.set("rotation 1", 90)

    # --- 5. External SOA input waveguide ---
    fdtd.addrect()
    fdtd.set("name", "Ext_SOA_input")
    fdtd.set("material", params.material_soa)
    fdtd.set("x", pos["ext_soa_start"] + params.wg2_length / 2)
    fdtd.set("y", 0)
    fdtd.set("z", 0)
    fdtd.set("x span", params.wg2_length)
    fdtd.set("y span", params.wg2_width)
    fdtd.set("z span", params.wg2_height)

    return pos


def setup_fdtd_simulation(fdtd, params, positions):
    """Configure the FDTD simulation region, source, and monitors.

    Args:
        fdtd: lumapi.FDTD handle.
        params: SOAPWBParams instance.
        positions: dict from compute_positions().
    """
    # --- Lateral margin for source/monitor spans ---
    y_margin = max(params.pwb_width, params.wg1_width,
                   params.taper1_w_out, params.taper2_w_in) * 1.5
    z_margin = max(params.pwb_height, params.wg1_height,
                   params.taper1_h_out, params.taper2_h_in) * 1.5

    x_min = -params.wg1_length * 0.5
    x_max = positions["total"] + params.wg2_length * 0.5

    fdtd.setresource("FDTD", "GPU", True)
    fdtd.addfdtd()
    fdtd.set("dimension", "3D")
    fdtd.set("x min", x_min)
    fdtd.set("x max", x_max)
    fdtd.set("y min", -y_margin)
    fdtd.set("y max", y_margin)
    fdtd.set("z min", -z_margin)
    fdtd.set("z max", z_margin)
    fdtd.set("x min bc", "PML")
    fdtd.set("x max bc", "PML")
    fdtd.set("y min bc", "PML")
    fdtd.set("y max bc", "PML")
    fdtd.set("z min bc", "PML")
    fdtd.set("z max bc", "PML")
    fdtd.set("mesh accuracy", params.mesh_accuracy)
    fdtd.set("mesh type", "auto non-uniform")
    fdtd.set("simulation time", params.simulation_time)

    # --- Mode source: inject fundamental TE mode from SOA output ---
    source_x = positions["soa_start"] + params.wg1_length * 0.3
    fdtd.addmode()
    fdtd.set("name", "source")
    fdtd.set("injection axis", "x-axis")
    fdtd.set("direction", "forward")
    fdtd.set("x", source_x)
    fdtd.set("y", 0)
    fdtd.set("z", 0)
    fdtd.set("y span", 2 * y_margin)
    fdtd.set("z span", 2 * z_margin)
    fdtd.set("mode selection", "fundamental TE mode")
    fdtd.set("wavelength start", params.wavelength)
    fdtd.set("wavelength stop", params.wavelength)

    # --- X-normal power monitors at key positions ---
    monitor_positions = [
        ("input_monitor", positions["soa_end"]),
        ("taper1_monitor", positions["taper1_end"]),
        ("pwb_monitor", positions["pwb_end"]),
        ("output_monitor", positions["ext_soa_start"]),
    ]
    for name, mx in monitor_positions:
        fdtd.addpower()
        fdtd.set("name", name)
        fdtd.set("monitor type", "2D X-normal")
        fdtd.set("x", mx)
        fdtd.set("y", 0)
        fdtd.set("z", 0)
        fdtd.set("y span", 2 * y_margin)
        fdtd.set("z span", 2 * z_margin)

    # --- Mode expansion at output ---
    fdtd.addmodeexpansion()
    fdtd.set("name", "mode_expansion")
    fdtd.setexpansion("input", "output_monitor")
    fdtd.set("mode selection", "fundamental TE mode")
    fdtd.set("x", positions["ext_soa_start"])
    fdtd.set("y", 0)
    fdtd.set("z", 0)
    fdtd.set("y span", 2 * y_margin)
    fdtd.set("z span", 2 * z_margin)

    # --- Y-normal transmission monitor (side view for field propagation) ---
    fdtd.addpower()
    fdtd.set("name", "transmission_monitor")
    fdtd.set("monitor type", "2D Y-normal")
    fdtd.set("y", 0)
    fdtd.set("x min", x_min)
    fdtd.set("x max", x_max)
    fdtd.set("z", 0)
    fdtd.set("z span", 2 * z_margin)


def get_data(fdtd, params):
    """Read simulation results from all monitors.

    Returns:
        dict with keys: source_E, transmission_E, monitors dict, T_forward.
    """
    source_E = fdtd.getresult("source", "mode profile")
    transmission_E = fdtd.getresult("transmission_monitor", "E")
    input_E = fdtd.getresult("input_monitor", "E")
    taper1_E = fdtd.getresult("taper1_monitor", "E")
    pwb_E = fdtd.getresult("pwb_monitor", "E")
    output_E = fdtd.getresult("output_monitor", "E")

    mode_data = fdtd.getresult("mode_expansion", "expansion for input")
    T_forward = mode_data.get("T_forward", None)

    return {
        "source_E": source_E,
        "transmission_E": transmission_E,
        "input_E": input_E,
        "taper1_E": taper1_E,
        "pwb_E": pwb_E,
        "output_E": output_E,
        "T_forward": T_forward,
    }


def _plot_x_normal_monitor(ax, monitor_data, title):
    """Plot intensity on an X-normal (Y-Z) monitor plane."""
    e_field = monitor_data["E"]
    ex = e_field[0, :, :, 0, 0]
    ey = e_field[0, :, :, 0, 1]
    ez = e_field[0, :, :, 0, 2]
    intensity = np.abs(ex) ** 2 + np.abs(ey) ** 2 + np.abs(ez) ** 2

    y = np.squeeze(monitor_data["y"]) * 1e6
    z = np.squeeze(monitor_data["z"]) * 1e6
    y_grid, z_grid = np.meshgrid(y, z, indexing="ij")
    ax.pcolormesh(y_grid, z_grid, intensity, cmap="jet", shading="auto")
    ax.set_xlabel("Y (um)")
    ax.set_ylabel("Z (um)")
    ax.set_title(title)


def visualize_and_save_results(fdtd, params, output_path=None):
    """Generate field visualization figure and return T_forward.

    Creates a 2x3 subplot: source mode + 4 cross-section monitors + transmission side view.
    """
    results = get_data(fdtd, params)

    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["xtick.labelsize"] = 12
    plt.rcParams["ytick.labelsize"] = 12

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    (ax_src, ax_in, ax_t1), (ax_pwb, ax_out, ax_trans) = axes

    _plot_x_normal_monitor(ax_src, results["source_E"], "Source Mode")
    _plot_x_normal_monitor(ax_in, results["input_E"], "Input (SOA End)")
    _plot_x_normal_monitor(ax_t1, results["taper1_E"], "Taper-1 End")
    _plot_x_normal_monitor(ax_pwb, results["pwb_E"], "PWB End")
    _plot_x_normal_monitor(ax_out, results["output_E"], "Output (Ext SOA)")

    # Transmission monitor (side view: X-Z plane, Y=0)
    trans_e = results["transmission_E"]
    e_field = trans_e["E"]
    ex = e_field[:, 0, :, 0, 0]
    ey = e_field[:, 0, :, 0, 1]
    ez = e_field[:, 0, :, 0, 2]
    intensity = np.abs(ex) ** 2 + np.abs(ey) ** 2 + np.abs(ez) ** 2

    x = np.squeeze(trans_e["x"]) * 1e6
    z = np.squeeze(trans_e["z"]) * 1e6
    x_grid, z_grid = np.meshgrid(x, z, indexing="ij")
    ax_trans.pcolormesh(x_grid, z_grid, intensity, cmap="jet", shading="auto")
    ax_trans.set_xlabel("X (um)")
    ax_trans.set_ylabel("Z (um)")
    ax_trans.set_title("Field Propagation (X-Z)")

    plt.suptitle("SOA-PWB-SOA Field Distribution", fontsize=14)
    plt.tight_layout()

    if output_path is None:
        PICTURES_DIR.mkdir(parents=True, exist_ok=True)
        output_path = PICTURES_DIR / "single_run.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()

    return results["T_forward"]


# ============================================================
# FDE Mode Analysis (requires Lumerical MODE Solutions)
# ============================================================


def run_fde_mode_analysis(params, mode_obj=None):
    """Compute fundamental mode field distributions at three cross-sections.

    Uses Lumerical MODE Solutions FDE (Finite-Difference Eigenmode) solver.

    Sections analyzed:
      1. SOA output waveguide cross-section
      2. PWB polymer waveguide cross-section
      3. External SOA input waveguide cross-section

    Args:
        params: SOAPWBParams instance.
        mode_obj: lumapi.MODE handle. If None, creates a new one.

    Returns:
        dict with keys: soa_mode, pwb_mode, ext_soa_mode.
        Each contains: E (3D array), x, y, neff.
    """
    close_on_exit = mode_obj is None
    if mode_obj is None:
        try:
            import lumapi
            mode_obj = lumapi.MODE()
        except AttributeError:
            print("WARNING: MODE Solutions not available. FDE analysis skipped.")
            return None

    try:
        results = {}

        # --- 1. SOA output waveguide mode ---
        results["soa_mode"] = _compute_single_mode(
            mode_obj, params,
            width=params.wg1_width,
            height=params.wg1_height,
            material=params.material_soa,
            label="SOA_output",
        )

        # --- 2. PWB waveguide mode ---
        results["pwb_mode"] = _compute_single_mode(
            mode_obj, params,
            width=params.pwb_width,
            height=params.pwb_height,
            material=params.material_pwb,
            label="PWB",
        )

        # --- 3. External SOA input waveguide mode ---
        results["ext_soa_mode"] = _compute_single_mode(
            mode_obj, params,
            width=params.wg2_width,
            height=params.wg2_height,
            material=params.material_soa,
            label="Ext_SOA_input",
        )

        return results

    finally:
        if close_on_exit and mode_obj is not None:
            mode_obj.close()


def _compute_single_mode(mode, params, width, height, material, label):
    """Compute fundamental TE mode for a single rectangular waveguide cross-section.

    Args:
        mode: lumapi.MODE handle.
        params: SOAPWBParams instance.
        width, height: waveguide cross-section dimensions [m].
        material: material name string.
        label: name for the solver region.

    Returns:
        dict with E, x, y, neff, or None on failure.
    """
    margin = max(width, height) * 2

    mode.addfde()
    mode.set("name", label)
    mode.set("solver type", "2D X normal")
    mode.set("x", 0)
    mode.set("y", 0)
    mode.set("z", 0)
    mode.set("y span", width + margin)
    mode.set("z span", height + margin)
    mode.set("wavelength", params.wavelength)
    mode.set("number of trial modes", 5)
    mode.set("mesh refinement", "conformal variant 1")

    # Add waveguide structure in MODE
    mode.addrect()
    mode.set("name", f"{label}_core")
    mode.set("material", material)
    mode.set("x", 0)
    mode.set("y", 0)
    mode.set("z", 0)
    mode.set("x span", 1e-9)  # thin in propagation direction
    mode.set("y span", width)
    mode.set("z span", height)

    # Add substrate
    mode.addrect()
    mode.set("name", f"{label}_substrate")
    mode.set("material", params.material_substrate)
    mode.set("x", 0)
    mode.set("y", 0)
    mode.set("z min", -margin * 2)
    mode.set("z max", 0)
    mode.set("x span", 1e-9)
    mode.set("y span", width + margin * 2)

    mode.run()
    mode.findmodes()

    try:
        neff = mode.getdata("FDE::data::mode1", "neff")
        E = mode.getresult("FDE::data::mode1", "mode profile")
        x = mode.getdata("FDE::data::mode1", "x")
        y = mode.getdata("FDE::data::mode1", "y")
        mode.deleteall()
        return {"E": E, "x": x, "y": y, "neff": neff}
    except Exception:
        mode.deleteall()
        return None


def compute_mode_overlap(mode1, mode2):
    """Compute the mode overlap factor eta between two modes.

    eta = |integral E1(x,y) . E2*(x,y) dxdy|^2 / [integral|E1|^2 dxdy * integral|E2|^2 dxdy]

    Args:
        mode1, mode2: dicts from _compute_single_mode() with keys E, x, y.

    Returns:
        float: overlap factor eta (0 to 1), or None on failure.
    """
    if mode1 is None or mode2 is None:
        return None

    E1 = mode1["E"]  # shape: (3, ny, nx, ...)
    E2 = mode2["E"]

    # Extract transverse components: Ex, Ey, Ez at the central slice
    e1_x = E1[0, :, :, 0, 0].flatten()
    e1_y = E1[1, :, :, 0, 0].flatten()
    e1_z = E1[2, :, :, 0, 0].flatten()

    e2_x = E2[0, :, :, 0, 0].flatten()
    e2_y = E2[1, :, :, 0, 0].flatten()
    e2_z = E2[2, :, :, 0, 0].flatten()

    # Full vector overlap
    numerator = np.abs(
        np.dot(e1_x, np.conj(e2_x))
        + np.dot(e1_y, np.conj(e2_y))
        + np.dot(e1_z, np.conj(e2_z))
    ) ** 2

    denom1 = (
        np.dot(e1_x, np.conj(e1_x))
        + np.dot(e1_y, np.conj(e1_y))
        + np.dot(e1_z, np.conj(e1_z))
    )
    denom2 = (
        np.dot(e2_x, np.conj(e2_x))
        + np.dot(e2_y, np.conj(e2_y))
        + np.dot(e2_z, np.conj(e2_z))
    )

    if denom1 == 0 or denom2 == 0:
        return 0.0

    return float(numerator / (denom1 * denom2))


def plot_mode_profiles(fde_results, output_path=None):
    """Plot normalized mode field distributions from FDE analysis.

    Args:
        fde_results: dict from run_fde_mode_analysis().
        output_path: optional path for saving the figure.
    """
    if fde_results is None:
        print("No FDE results to plot.")
        return

    plt.rcParams["font.family"] = "Times New Roman"
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    labels = ["SOA Output Mode", "PWB Mode", "Ext SOA Input Mode"]
    keys = ["soa_mode", "pwb_mode", "ext_soa_mode"]

    for ax, label, key in zip(axes, labels, keys):
        mode = fde_results.get(key)
        if mode is None:
            ax.set_title(f"{label}\n(not available)")
            continue
        E = mode["E"]
        intensity = (
            np.abs(E[0, :, :, 0, 0]) ** 2
            + np.abs(E[1, :, :, 0, 0]) ** 2
            + np.abs(E[2, :, :, 0, 0]) ** 2
        )
        y = np.squeeze(mode["y"]) * 1e6
        z = np.squeeze(mode["x"]) * 1e6  # MODE may swap x/y axes
        y_grid, z_grid = np.meshgrid(y, z, indexing="ij")
        ax.pcolormesh(y_grid, z_grid, intensity / np.max(intensity),
                       cmap="hot", shading="auto")
        ax.set_xlabel("Y (um)")
        ax.set_ylabel("Z (um)")
        neff = mode.get("neff", "?")
        ax.set_title(f"{label}\nn_eff = {neff}")

    plt.suptitle("Fundamental Mode Profiles (|E|^2 normalized)", fontsize=14)
    plt.tight_layout()

    if output_path is not None:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()
