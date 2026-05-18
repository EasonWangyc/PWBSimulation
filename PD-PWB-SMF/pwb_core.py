"""
pwb_core.py — PD-PWB-SMF 仿真核心库
包含参数类、几何生成、仿真设置、数据读取及可视化函数。
"""

import os
import numpy as np
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd


# ---------------------------------------------------------------------------
# 参数类
# ---------------------------------------------------------------------------

class PWBParameters:
    def __init__(self):
        self.L = 250e-6
        self.h = 150e-6
        self.r1 = 6.2e-6
        self.r = 0.8e-6
        self.R = 60e-6
        self.r2 = 7e-6
        self.l1 = 100e-6
        # self.l2 = 60e-6
        self.curve_points = 100
        self.wavelength = 1.55e-6
        self.mesh_accuracy = 1
        self.simulation_time = 2000e-15


# ---------------------------------------------------------------------------
# 路径生成
# ---------------------------------------------------------------------------

def generate_pwb_path_1(params):
    curve_points = params.curve_points
    L = params.L
    h = params.h
    r1 = params.r1
    r = params.r
    R = params.R
    r2 = params.r2
    l1 = params.l1
    l2 = h - R
    l = L - l1 - R

    # taper1
    x1 = np.linspace(0, l1, curve_points)
    y1 = np.zeros(curve_points)
    z1 = np.zeros(curve_points)

    # bending
    angle_start = 0
    angle_end = np.pi / 2
    t = np.linspace(angle_start, angle_end, 2 * curve_points)
    center_x = l1 + l
    center_z = -R
    x2 = center_x + R * np.sin(t)
    y2 = np.zeros_like(x2)
    z2 = center_z + R * np.cos(t)

    # taper2
    z3 = np.linspace(z2[-1], z2[-1] - l2, curve_points)
    y3 = np.zeros(curve_points)
    x3 = np.zeros(curve_points) + x2[-1]

    x = np.concatenate((x1, x2, x3))
    y = np.zeros_like(x)
    z = np.concatenate((z1, z2, z3))
    path = np.column_stack((x, y, z))
    return path


def generate_pwb_path_2(params):
    curve_points = params.curve_points
    L = params.L
    h = params.h
    r1 = params.r1
    r = params.r
    R = params.R
    r2 = params.r2
    l1 = params.l1
    l2 = h - R
    l = L - l1 - R

    # taper1
    x1 = np.linspace(0, l1, curve_points)
    y1 = np.zeros(curve_points)
    z1 = np.zeros(curve_points)

    # bending
    angle_start = 0
    angle_end = np.pi / 2
    t = np.linspace(angle_start, angle_end, 2 * curve_points)
    center_x = l1 + l
    center_z = -R
    x2 = center_x + R * np.sin(t)
    y2 = np.zeros_like(x2)
    z2 = center_z + R * np.cos(t)

    # taper2
    z3 = np.linspace(z2[-1], z2[-1] - l2, curve_points)
    y3 = np.zeros(curve_points)
    x3 = np.zeros(curve_points) + x2[-1]

    x = np.concatenate((x1, x2, x3))
    y = np.zeros_like(x)
    z = np.concatenate((z1, z2, z3))
    path = np.column_stack((x, y, z))
    return path


# ---------------------------------------------------------------------------
# 结构生成
# ---------------------------------------------------------------------------

def generate_pwb_structure_1(fdtd, params):
    """按照参数生成对应结构（Section 1：bending 段）"""
    L = params.L
    h = params.h
    r1 = params.r1
    r = params.r
    r2 = params.r2
    l1 = params.l1
    R = params.R
    l2 = h - R
    l = L - l1 - R
    path = generate_pwb_path_1(params)

    fdtd.deleteall()
    fdtd.load("D:/simulation/Simulation Project/simulation/PD-PWB-SMF/SMF.fsp")
    fdtd.importmaterialdb("D:/simulation/Simulation Project/simulation/database.mdf")

    # 创建 bending
    num_segments_2 = 200
    path_length_2 = 200
    segment_length_2 = path_length_2 // num_segments_2
    for i in range(num_segments_2):
        start_idx = i * segment_length_2 + 100
        end_idx = min((i + 1) * segment_length_2, path_length_2) + 100

        if start_idx >= end_idx:
            continue

        start_pos = path[start_idx]
        end_pos = path[end_idx]

        center_x = (start_pos[0] + end_pos[0]) / 2
        center_y = (start_pos[1] + end_pos[1]) / 2
        center_z = (start_pos[2] + end_pos[2]) / 2

        segment_len = np.linalg.norm(end_pos - start_pos)

        direction = end_pos - start_pos
        norm = np.linalg.norm(direction)
        if norm > 0:
            direction = direction / norm
        else:
            direction = np.array([1, 0, 0])
        dz, dy, dx = direction[2], direction[1], direction[0]
        theta = np.arctan2(dx, dz) * 180 / np.pi  # 绕 y 轴
        phi = np.arcsin(dy) * 180 / np.pi          # 绕 x 轴

        fdtd.addcircle()
        fdtd.set("name", f"Bending_{i}")
        fdtd.set("material", "Vancore B")
        fdtd.set("make ellipsoid", 0)
        fdtd.set("x", center_x)
        fdtd.set("y", center_y)
        fdtd.set("z", center_z)
        fdtd.set("radius", r)
        fdtd.set("z span", segment_len)
        fdtd.set("first axis", "y")
        fdtd.set("rotation 1", theta)
        fdtd.set("second axis", "x")
        fdtd.set("rotation 2", phi)


def generate_pwb_structure_2(fdtd, params):
    """按照参数生成对应结构（Section 2：bending + taper2 段）"""
    L = params.L
    h = params.h
    r1 = params.r1
    r = params.r
    r2 = params.r2
    l1 = params.l1
    R = params.R
    l2 = h - R
    l = L - l1 - R
    path = generate_pwb_path_2(params)

    fdtd.deleteall()
    fdtd.load("D:/simulation/Simulation Project/simulation/PD-PWB-SMF/SMF.fsp")
    fdtd.importmaterialdb("D:/simulation/Simulation Project/simulation/database.mdf")

    # 创建 bending
    num_segments_2 = 200
    path_length_2 = 200
    segment_length_2 = path_length_2 // num_segments_2
    for i in range(num_segments_2):
        start_idx = i * segment_length_2 + 100
        end_idx = min((i + 1) * segment_length_2, path_length_2) + 100

        if start_idx >= end_idx:
            continue

        start_pos = path[start_idx]
        end_pos = path[end_idx]

        center_x = (start_pos[0] + end_pos[0]) / 2
        center_y = (start_pos[1] + end_pos[1]) / 2
        center_z = (start_pos[2] + end_pos[2]) / 2

        segment_len = np.linalg.norm(end_pos - start_pos)

        direction = end_pos - start_pos
        norm = np.linalg.norm(direction)
        if norm > 0:
            direction = direction / norm
        else:
            direction = np.array([1, 0, 0])
        dz, dy, dx = direction[2], direction[1], direction[0]
        theta = np.arctan2(dx, dz) * 180 / np.pi
        phi = np.arcsin(dy) * 180 / np.pi

        fdtd.addcircle()
        fdtd.set("name", f"Bending_{i}")
        fdtd.set("material", "Vancore B")
        fdtd.set("make ellipsoid", 0)
        fdtd.set("x", center_x)
        fdtd.set("y", center_y)
        fdtd.set("z", center_z)
        fdtd.set("radius", r)
        fdtd.set("z span", segment_len)
        fdtd.set("first axis", "y")
        fdtd.set("rotation 1", theta)
        fdtd.set("second axis", "x")
        fdtd.set("rotation 2", phi)

    # 创建 taper2
    num_segments_3 = 100
    for i in range(num_segments_3):
        radius_2 = r + (r2 - r) * (i + 1) / 100
        fdtd.addcircle()
        fdtd.set("name", f"taper2_{i}")
        fdtd.set("material", "Vancore B")
        fdtd.set("make ellipsoid", 0)
        fdtd.set("x", L)
        fdtd.set("y", 0)
        fdtd.set("z", -R - (i + 1 / 2) * l2 / 100)
        fdtd.set("radius", radius_2)
        fdtd.set("z span", l2 / 100)


# ---------------------------------------------------------------------------
# 仿真设置
# ---------------------------------------------------------------------------

def setup_fdtd_simulation_1(fdtd, params):
    """设置 Section 1 仿真区域、光源和监视器"""
    L = params.L
    h = params.h
    r1 = params.r1
    r = params.r
    R = params.R
    r2 = params.r2
    l1 = params.l1
    l2 = h - R
    l = L - l1 - R

    x_min, x_max = L - R, L
    y_min, y_max = 0, 0
    z_min, z_max = -R, 0
    margin = 6e-6

    fdtd.setresource("FDTD", "GPU", True)
    fdtd.addfdtd()
    fdtd.set("dimension", "3D")
    fdtd.set("x min", x_min)
    fdtd.set("x max", x_max + margin)
    fdtd.set("y min", y_min - margin)
    fdtd.set("y max", y_max + margin)
    fdtd.set("z min", z_min)
    fdtd.set("z max", z_max + margin)
    fdtd.set("express mode", True)
    fdtd.set("x min bc", "PML")
    fdtd.set("x max bc", "PML")
    fdtd.set("y min bc", "PML")
    fdtd.set("y max bc", "PML")
    fdtd.set("z min bc", "PML")
    fdtd.set("z max bc", "PML")
    fdtd.set("mesh accuracy", params.mesh_accuracy)
    fdtd.set("mesh type", "auto non-uniform")
    fdtd.set("simulation time", params.simulation_time)

    fdtd.addmode()
    fdtd.set("name", "source")
    fdtd.set("injection axis", "x-axis")
    fdtd.set("direction", "forward")
    fdtd.set("x", x_min + 1e-6)
    fdtd.set("y", 0)
    fdtd.set("z", 0)
    fdtd.set("y span", 2 * margin)
    fdtd.set("z span", 2 * margin)
    fdtd.set("wavelength start", params.wavelength)
    fdtd.set("wavelength stop", params.wavelength)

    fdtd.addpower()
    fdtd.set("name", "monitor_1")
    fdtd.set("monitor type", "2D Z-normal")
    fdtd.set("x", L)
    fdtd.set("y", 0)
    fdtd.set("z", -R)
    fdtd.set("y span", 2 * margin)
    fdtd.set("x span", 2 * margin)

    fdtd.addpower()
    fdtd.set("name", "transmission_monitor")
    fdtd.set("monitor type", "2D Y-normal")
    fdtd.set("y", 0)
    fdtd.set("x", (2 * L - R) / 2)
    fdtd.set("x span", R + 6e-6)
    fdtd.set("z", -R / 2)
    fdtd.set("z span", R + 6e-6)


def setup_fdtd_simulation_2(fdtd, params):
    """设置 Section 2 仿真区域、光源和监视器"""
    L = params.L
    h = params.h
    r1 = params.r1
    r = params.r
    R = params.R
    r2 = params.r2
    l1 = params.l1
    l2 = h - R
    l = L - l1 - R

    x_min, x_max = L - R, L
    y_min, y_max = 0, 0
    z_min, z_max = -h, 0
    margin = 10e-6

    fdtd.setresource("FDTD", "GPU", True)
    fdtd.addfdtd()
    fdtd.set("dimension", "3D")
    fdtd.set("x min", x_min - margin)
    fdtd.set("x max", x_max + margin)
    fdtd.set("y min", y_min - margin)
    fdtd.set("y max", y_max + margin)
    fdtd.set("z min", z_min - margin)
    fdtd.set("z max", z_max + margin)
    fdtd.set("express mode", True)
    fdtd.set("x min bc", "PML")
    fdtd.set("x max bc", "PML")
    fdtd.set("y min bc", "PML")
    fdtd.set("y max bc", "PML")
    fdtd.set("z min bc", "PML")
    fdtd.set("z max bc", "PML")
    fdtd.set("mesh accuracy", params.mesh_accuracy)
    fdtd.set("mesh type", "auto non-uniform")
    fdtd.set("simulation time", params.simulation_time)

    fdtd.addmode()
    fdtd.set("name", "source")
    fdtd.set("injection axis", "x-axis")
    fdtd.set("direction", "forward")
    fdtd.set("x", L - R + 2e-6)
    fdtd.set("y", 0)
    fdtd.set("z", z_max)
    fdtd.set("y span", 3e-6)
    fdtd.set("z span", 3e-6)
    fdtd.set("wavelength start", params.wavelength)
    fdtd.set("wavelength stop", params.wavelength)

    fdtd.addpower()
    fdtd.set("name", "monitor_1")
    fdtd.set("monitor type", "2D Z-normal")
    fdtd.set("x", L)
    fdtd.set("y", 0)
    fdtd.set("z", -h)
    fdtd.set("y span", 2 * margin)
    fdtd.set("x span", 2 * margin)

    fdtd.addpower()
    fdtd.set("name", "transmission_monitor")
    fdtd.set("monitor type", "2D Y-normal")
    fdtd.set("y", 0)
    fdtd.set("x", (2 * L - R) / 2)
    fdtd.set("x span", R + 2 * margin)
    fdtd.set("z", -h / 2)
    fdtd.set("z span", h + 2 * margin)


# ---------------------------------------------------------------------------
# 数据读取
# ---------------------------------------------------------------------------

def get_data_1(fdtd, params):
    """读取 Section 1 仿真结果"""
    source_E = fdtd.getresult("source", "mode profile")
    transmission_E = fdtd.getresult("transmission_monitor", "E")
    transmission_P = fdtd.getresult("transmission_monitor", "P")
    monitor_1_E = fdtd.getresult("monitor_1", "E")
    T_total = fdtd.getresult("monitor_1", "T")

    results = {
        'source_E': source_E,
        'monitor_1_E': monitor_1_E,
        'transmission_E': transmission_E,
        'transmission_P': transmission_P,
        'T_total': T_total["T"][0],
    }
    return results


def get_data_2(fdtd, params):
    """读取 Section 2 仿真结果"""
    source_E = fdtd.getresult("source", "mode profile")
    transmission_E = fdtd.getresult("transmission_monitor", "E")
    transmission_P = fdtd.getresult("transmission_monitor", "P")
    monitor_1_E = fdtd.getresult("monitor_1", "E")
    T_total = fdtd.getresult("monitor_1", "T")

    results = {
        'source_E': source_E,
        'monitor_1_E': monitor_1_E,
        'transmission_E': transmission_E,
        'transmission_P': transmission_P,
        'T_total': T_total["T"][0],
    }
    return results


# ---------------------------------------------------------------------------
# 可视化
# ---------------------------------------------------------------------------

def visualize_and_save_results_1(fdtd, params):
    """可视化 Section 1 电场分布并保存图片，返回 T_total"""
    results = get_data_1(fdtd, params)

    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    fig, ax = plt.subplots(figsize=(6, 5))
    trans_E = results['transmission_E']
    E_field = trans_E['E']
    Ex = E_field[:, 0, :, 0, 0]
    Ey = E_field[:, 0, :, 0, 1]
    Ez = E_field[:, 0, :, 0, 2]
    E_intensity_trans = np.abs(Ex) ** 2 + np.abs(Ey) ** 2 + np.abs(Ez) ** 2
    x_trans = np.squeeze(trans_E['x']) * 1e6
    z_trans = np.squeeze(trans_E['z']) * 1e6
    X_trans, Z_trans = np.meshgrid(x_trans, z_trans, indexing='ij')

    vmin = 1e-5
    vmax = np.nanmax(E_intensity_trans)
    norm = LogNorm(vmin=vmin, vmax=vmax)
    pcm = ax.pcolormesh(X_trans, Z_trans, E_intensity_trans, cmap="inferno", shading='auto', norm=norm)
    fig.colorbar(pcm, ax=ax, pad=0.02, format=ticker.LogFormatterSciNotation())

    ax.set_xlabel('X (μm)')
    ax.set_ylabel('Z (μm)')
    ax.set_title('Transmission Monitor')
    plt.tight_layout()

    out_path = "D:/simulation/Simulation Project/simulation/PD-PWB-SMF/results/Pictures/Section1.png"
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.show()

    return results['T_total']


def plot_Ttotal_loss_heatmap_rR(data_file, output_dir=None, cmap='viridis_r', eps=1e-12, show=True):
    """
    读取 r, R, T_total 数据文件，计算 loss = -10*log10(|T_total|) 并绘制热力图。
    返回 (loss_matrix, r_vals, R_vals)
    """
    for enc in ('utf-8', 'gbk', 'latin-1'):
        try:
            df = pd.read_csv(data_file, sep=r'\s*\t\s*', engine='python', encoding=enc)
            break
        except Exception:
            continue

    cols = [c.strip() for c in df.columns]
    df.columns = cols

    if len(cols) == 3:
        r_col, R_col, T_col = cols[0], cols[1], cols[2]
    else:
        lower = [c.lower() for c in cols]
        r_col = cols[lower.index(next(x for x in lower if x.startswith('r')))] if any(x.startswith('r') for x in lower) else cols[0]
        R_col = cols[lower.index(next(x for x in lower if x.startswith('r') and x != r_col.lower()))] if len(cols) > 1 else cols[1]
        T_col = next((c for c in cols if 't' in c.lower()), cols[-1])

    r_vals = np.sort(np.unique(df[r_col].values.astype(float)))
    R_vals = np.sort(np.unique(df[R_col].values.astype(float)))

    loss_matrix = np.full((len(r_vals), len(R_vals)), np.nan)

    for _, row in df.iterrows():
        try:
            i = int(np.argwhere(r_vals == float(row[r_col]))[0])
            j = int(np.argwhere(R_vals == float(row[R_col]))[0])
        except Exception:
            i = int(np.argmin(np.abs(r_vals - float(row[r_col]))))
            j = int(np.argmin(np.abs(R_vals - float(row[R_col]))))
        T = float(row[T_col])
        T_pos = np.clip(np.abs(-1.0 * T), eps, None)
        loss_matrix[i, j] = -10.0 * np.log10(T_pos)

    plt.rcParams['font.family'] = 'Times New Roman'
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(loss_matrix, origin='lower', aspect='auto', cmap=cmap)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('loss (dB)', fontsize=12)

    ax.set_xticks(np.arange(len(R_vals)))
    ax.set_xticklabels([f"{v:.0f}" if v == int(v) else f"{v:.1f}" for v in R_vals], rotation=45)
    ax.set_yticks(np.arange(len(r_vals)))
    ax.set_yticklabels([f"{v:.2f}" if v < 10 else f"{v:.1f}" for v in r_vals])

    ax.set_xlabel('R (μm)', fontsize=12)
    ax.set_ylabel('r (μm)', fontsize=12)
    ax.set_title('Loss heatmap', fontsize=14)
    plt.tight_layout()

    if output_dir is None:
        output_dir = os.path.dirname(data_file) or '.'
    os.makedirs(output_dir, exist_ok=True)
    base = os.path.splitext(os.path.basename(data_file))[0]
    out_path = os.path.join(output_dir, f"{base}_loss_heatmap.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)

    print("Saved heatmap to:", out_path)
    return loss_matrix, r_vals, R_vals
