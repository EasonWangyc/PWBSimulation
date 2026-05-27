# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 项目概要

PWB（光子引线键合）仿真项目，通过 Python 调用 Lumerical/Ansys FDTD API（`lumapi`）。按耦合场景分目录：`PD-PWB-SMF/`、`LD-PWB-SMF/`、`LNOI-PWB-SMF/`，每个目录有各自的 `pwb_core.py`、运行脚本和结果。详见 `AGENTS.md`。

## 常用命令

```bash
# 语法检查（不需要 Lumerical）
python -m py_compile sim_config.py
python -m py_compile PD-PWB-SMF/pwb_core.py

# 纯几何单元测试（不需要 Lumerical）
python PD-PWB-SMF/test_complex_geometry.py

# 后处理分析（不需要 Lumerical）
python PD-PWB-SMF/analyze_results.py

# 仅构建结构，不运行 FDTD（需要 Lumerical）
python PD-PWB-SMF/test_setup.py
python LD-PWB-SMF/test_setup.py
python LNOI-PWB-SMF/test_setup.py

# 单次 FDTD 仿真（需要 Lumerical，耗时长）
python PD-PWB-SMF/run_single.py

# 参数扫描（需要 Lumerical，耗时极长——必须先征得同意）
python PD-PWB-SMF/sweep_r_R.py
python PD-PWB-SMF/sweep_h_R.py
python LNOI-PWB-SMF/sweep_h1_h2_w1_w2.py
```

## 架构要点

### 路径体系

所有路径通过 `sim_config.py` 统一管理，脚本从中导入常量，不硬编码路径。三个环境变量可覆盖默认值：`SIM_PROJECT_ROOT`、`LUMERICAL_API_PATH`、`SIM_MATERIAL_DB`。

### 三个场景的统一结构

每个场景目录（`PD-PWB-SMF/`、`LD-PWB-SMF/`、`LNOI-PWB-SMF/`）结构一致：

| 文件 | 作用 |
|------|------|
| `pwb_core.py` | 参数类 + 几何生成 + FDTD 设置 + 结果读取 + 可视化 |
| `test_setup.py` | 仅构建并保存 `.fsp`，不运行仿真 |
| `run_single.py` | 一次完整仿真 |
| `parameter.md` | 该场景的物理参数说明 |
| `*.ipynb` | 历史探索记录，不要大改，可复用逻辑应抽取到 `.py` |

### PD-PWB-SMF 的三种中心线模式

`pwb_core.py` 中用 `_1`/`_2`/`_3` 后缀区分：

- **`_1`**：简单 90° 四分之一圆弧弯曲
- **`_2`**：弯曲段 + 垂直 taper2
- **`_3`**：分段三次贝塞尔复杂中心线（arch + drop 两段），曲率由 `bend_lift`、`arch_position`、`bend_shape`、`drop_shape` 控制。**注意：`_3` 中的 `R` 参数不直接控制曲率**，仅用于估算输出 taper 长度。`_1`/`_2`/`_3` 中的弯曲半径参数物理含义不同，不要混淆。

纯几何辅助函数（不依赖 `lumapi`）：`cubic_bezier()`、`path_arc_length()`、`normalized_path_position()`、`centerline_to_segments()`、`generate_pwb_path_3()`、`generate_radius_profile_3()`。

### LD-PWB-SMF：导入光源

使用导入光源数据集 `Section2 output.mat`。文件缺失时自动回退到普通 mode source，由 `params.use_imported_source` 控制。

### LNOI-PWB-SMF：pyramid + 直波导

用 `addpyramid()` 构建 taper，`addrect()` 构建直波导。参数 `h1/h2` 是半厚度（旋转 90° 后的 x span 一半），不是直接的高度值。

## 开发约定

- **单位**：几何参数统一用 SI（米），通常以 `数值 * 1e-6` 表示微米
- **T_total 与 loss**：`loss = -10 * log10(abs(-1.0 * T_total))`，修改前需确认物理约定
- **结果输出**：放到对应场景的 `results/` 子目录，图片文件名需编码扫描参数值
- **Notebook**：保留不动，可复用逻辑抽取到 `.py` 文件
- **修改范围**：限制在用户提到的具体场景目录内，不要擅自重组目录结构
- **FDTD 代价高**：优先做代码检查和后处理，未经明确要求不要运行完整仿真或参数扫描
