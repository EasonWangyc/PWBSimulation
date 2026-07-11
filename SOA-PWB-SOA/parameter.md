# SOA-PWB-SOA 参数说明

本目录用于 片上SOA输出端 → PWB → 外置SOA输入端 的耦合仿真。

两端 SOA 波导结构由 `SOA_base_with_ar&cladding.fsp` 定义（脊形波导，含200nm有源区+增透膜+包层）。
PWB 以 `addcircle()` + ellipsoid 构建，截面为椭圆。

提供：
- `test_setup.py`：只生成并保存结构，不运行 FDTD。
- `run_single.py`：生成结构、运行一次仿真并绘图。
- `pwb_core.py`：核心模块，包含参数、几何、FDTD设置、FDE分析和可视化。
- `sweep_taper_width.py`：2D 网格扫描 taper 截面尺寸，隔离式小窗口 FDTD。
- `sweep_taper_length.py`：1D 扫描 taper 长度，mode expansion 直读 T_forward/T_backward。
- `plot_sweep_results.py`：纯数据分析脚本（无需 Lumerical），读取 CSV 生成 6 面板汇总图。

## 0. 路径配置

脚本中的路径统一从根目录 `config/sim_config.py` 获取：
- `SOA_DIR`：当前场景目录。
- `SOA_RESULTS_DIR`：结果输出目录。
- `SOA_BASE_FSP`：基础工程文件（含 SOA 波导 + 增透膜 + 包层）。
- `SOA_TILTED_FSP`：倾斜 SOA 工程文件（用于 `run_tilted.py`）。
- `MATERIAL_DB`：材料库文件。
- `add_lumerical_api_path()`：将 Lumerical Python API 路径加入 `sys.path`。

## 1. 几何结构

结构沿 x 轴直线排列，PWB 由 addcircle() 段组成，总长 250μm：

```
x=0 ──────────────────────────── x=250μm
  [SOA_out, .fsp]  Taper-1    PWB直段    Taper-2  [SOA_in, .fsp]
                   (r_in→r_pwb) (r_pwb恒定) (r_pwb→r_out)
```

PWB 截面为椭圆形（`use_ellipsoid=True`），两个轴半径独立可调。

核心几何参数（单位：SI，以 `* 1e-6` 表示微米）：

| 参数 | 说明 | 默认值 | 最优值 |
|------|------|--------|--------|
| `total_length` | PWB 总长度（含 taper） | 250e-6 | — |
| `taper1_length` | Taper-1 长度 | 100e-6 | 待定 |
| `taper2_length` | Taper-2 长度 | 100e-6 | 待定 |

PWB 直段长度 = `total_length - taper1_length - taper2_length`（自动计算）。

### 半径剖面（椭圆截面，两个独立轴）

| 参数 | 说明 | 默认值 | 最优值（宽度扫描） |
|------|------|--------|-------------------|
| `r_in` | Taper-1 SOA侧半径 (axis-1) | 1.3e-6 | **0.5e-6** |
| `r_in_2` | Taper-1 SOA侧半径 (axis-2) | 2.2e-6 | **2.8e-6** |
| `r_pwb` | PWB 中间段半径 (axis-1) | 0.9e-6 | **0.8e-6**（固定） |
| `r_pwb_2` | PWB 中间段半径 (axis-2) | 0.9e-6 | **0.8e-6**（固定） |
| `r_out` | Taper-2 SOA侧半径 (axis-1) | 0.8e-6 | **0.5e-6** |
| `r_out_2` | Taper-2 SOA侧半径 (axis-2) | 1.6e-6 | **1.5e-6** |
| `use_ellipsoid` | 椭圆截面 | True | — |

半径沿 x 轴的演化（线性插值）：
- Taper-1（0 → taper1_length）：r_in / r_in_2 → r_pwb / r_pwb_2
- PWB 直段：r_pwb / r_pwb_2 恒定
- Taper-2（total_length−taper2_length → total_length）：r_pwb / r_pwb_2 → r_out / r_out_2

### 两端 SOA 波导参数（由 .fsp 定义）

| 参数 | PE-SOA（左侧输入） | PG-SOA（右侧输出） |
|------|-------------------|-------------------|
| 波导类型 | 脊形（600μm渐变 2→5μm） | 脊形（固定 1μm×3μm） |
| 脊高度 | 2 μm | 3 μm |
| 有源层 | 200 nm | — |

### 截面扫描范围（sweep_taper_width.py）

| 扫描组 | 参数 | 范围 | 步长 | 点数 |
|--------|------|------|------|------|
| Taper-1 | r_in | 0.7–1.3 μm | 0.1 | 7 |
| Taper-1 | r_in_2 | 2.4–3.6 μm | 0.2 | 7 |
| Taper-2 | r_out | 0.2–0.8 μm | 0.1 | 7 |
| Taper-2 | r_out_2 | 0.9–2.1 μm | 0.2 | 7 |

采用 **2D 网格扫描**（r_in × r_in_2 全部 7×7=49 组合 + r_out × r_out_2 全部 7×7=49 组合）。

### 长度扫描范围（sweep_taper_length.py）

| 扫描组 | 范围 | 步长 | 点数 |
|--------|------|------|------|
| taper1_length | 10–100 μm | 10 | 11 |
| taper2_length | 10–100 μm | 10 | 11 |

## 2. 材料参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `material_pwb` | PWB聚合物材料 | "Vancore B" |

SOA 材料在 .fsp 中定义（InP 脊形 + SiO2 包层 + 增透膜）。

## 3. 仿真参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `wavelength` | 中心波长 | 1.55e-6 |
| `mesh_accuracy` | 网格精度 (1-8) | 1 |
| `simulation_time` | 仿真时间 [s] | 2000e-15（宽度扫描）/ 4000e-15（长度扫描） |
| `curve_points` | 路径采样点数 | 200 |

## 4. 监视器与仿真策略

### 宽度扫描（隔离模式，sweep_taper_width.py）

每个 taper 独立测试，FDTD 窗口极小（~20 μm）：

**Taper-1（forward）：**
- FDTD x ∈ [-10, 10] μm
- 源：SOA 输出波导内 x=-5 μm，forward
- Monitor + Mode Expansion：x=5 μm（进入 taper 5 μm 处）
- 读取：`T_forward`

**Taper-2（backward）：**
- FDTD x ∈ [240, 260] μm
- 源：SOA 输入波导内 x=255 μm，**backward**
- Monitor + Mode Expansion：x=245 μm（进入 taper 5 μm 处）
- 读取：`T_backward`

### 长度扫描（隔离模式，sweep_taper_length.py）

FDTD 窗口随 taper 长度 L 动态缩放：

**Taper-1（forward）：**
- FDTD x ∈ [-10, L+8] μm
- 源：x=-5 μm，forward
- Monitor + Mode Expansion：x=L（taper 结束处）
- 读取：`T_forward`

**Taper-2（backward）：**
- FDTD x ∈ [242−L, 260] μm
- 源：x=255 μm，backward
- Monitor + Mode Expansion：x=250−L（taper 结束处）
- 读取：`T_backward`

## 5. 运行方式

只生成结构：
```powershell
& "D:\Program Files\Lumerical\python\python.exe" SOA-PWB-SOA\test_setup.py
```

单次仿真：
```powershell
& "D:\Program Files\Lumerical\python\python.exe" SOA-PWB-SOA\run_single.py
```

宽度扫描（2D grid，~98 次仿真，耗时较长）：
```powershell
& "D:\Program Files\Lumerical\python\python.exe" SOA-PWB-SOA\sweep_taper_width.py
```

长度扫描（~22 次仿真）：
```powershell
& "D:\Program Files\Lumerical\python\python.exe" SOA-PWB-SOA\sweep_taper_length.py
```

结果绘图（纯数据分析，无需 Lumerical）：
```powershell
& "D:\Program Files\Lumerical\python\python.exe" SOA-PWB-SOA\plot_sweep_results.py
```

## 6. 输出文件

| 脚本 | 输出 |
|------|------|
| `sweep_taper_width.py` | `results/sweep_width_results.csv` |
| `sweep_taper_length.py` | `results/sweep_length_results.csv` |
| `plot_sweep_results.py` | `results/Pictures/sweep_summary.png`（6 面板汇总图） |
