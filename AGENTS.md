# AGENTS.md

## 项目背景

这个工作区是一个 PWB（Photonic Wire Bonding，光子引线键合）仿真项目。核心流程是通过 Python 调用 Lumerical/Ansys FDTD Python API（`lumapi`）来构建仿真结构、运行 FDTD 仿真，并对导出的数据、图片等结果进行分析。

项目按耦合场景组织：

- `PD-PWB-SMF/`：当前脚本化程度最高的目录。包含可复用 Python 代码 `pwb_core.py`、单次运行脚本、参数扫描脚本、`.fsp` 工程文件、notebook 和已生成结果。
- `LD-PWB-SMF/`：LD 到 PWB 到 SMF 的仿真资源，包含 notebook、`.fsp` 工程、LSF 脚本片段和扫描结果。
- `LNOI-PWB-SMF/`：LNOI 到 PWB 到 SMF 的仿真资源，以及参数扫描图片和 CSV 结果。
- `SMF-PWB-SMF/`：SMF 到 PWB 到 SMF 的仿真工程和导出结果。
- `Taper/`：用于构建 taper 的 Lumerical 脚本片段。
- `U-Turn/`：用途需先检查后再判断，不要直接假设其角色。
- 根目录下的 `data analysis.ipynb` 和 `results_analysis.ipynb` 主要用于数据分析和绘图。

大型二进制文件和生成结果是研究流程的一部分：

- `.fsp` 文件是 Lumerical FDTD 工程文件。
- `results/` 下的 `.h5`、`.png`、`.jpg`、`.txt`、`.csv` 文件通常是仿真输出或分析产物。
- 除非用户明确要求，不要删除、重新生成或批量改动这些结果文件。

## 运行前提

运行仿真脚本需要本机安装 Lumerical/Ansys，并能访问 FDTD Python API。现有脚本使用如下路径引入 `lumapi`：

```python
sys.path.append("D:\\Program Files\\Lumerical\\v241\\api\\python\\")
import lumapi
```

脚本还依赖常见科学计算包，例如 `numpy`、`matplotlib` 和 `pandas`。

FDTD 仿真通常耗时较长，并依赖商业软件许可、GPU 和本机资源。不要随意运行完整参数扫描；除非用户明确要求，优先做代码检查、轻量分析或后处理。

## 路径注意事项

部分脚本当前使用绝对 Windows 路径，例如：

```text
D:/simulation/Simulation Project/simulation/PD-PWB-SMF/...
```

当前工作区根目录是：

```text
D:/simulation/Simulation Project
```

在运行或修改脚本前，需要确认路径中额外的 `simulation/` 目录层级是否是用户机器上的真实结构，还是旧目录布局遗留。若需要改路径，应保持改动范围很小，并说明所依据的假设。

## 关键 Python 入口

`PD-PWB-SMF/` 中的主要脚本：

- `pwb_core.py`：共享参数类、PWB 路径生成、结构构建、FDTD 设置、结果读取和绘图辅助函数。
- `run_single.py`：构建 Section 1，保存 `.fsp`，运行 FDTD，绘制电场图片，并打印 `T_total`。
- `sweep_r_R.py`：对波导半径 `r` 和弯曲半径 `R` 做二维扫描；输出 `results/r_R_scan/T_total_sweep_results.txt` 和对应电场图。
- `sweep_h_R.py`：对总高度 `h` 和弯曲半径 `R` 做二维扫描；输出 `results/h_R_scan/T_total_sweep_results.txt`。
- `analyze_results.py`：读取已有扫描结果并绘制 loss 热力图。它不启动 FDTD，是最适合轻量运行的脚本。
- `test_setup.py`：只构建并保存 FDTD 结构，不执行完整仿真；但仍然依赖 `lumapi`。

## 开发约定

- 修改应尽量限制在用户提到的具体场景目录和脚本内。
- 若要增强可复现性，优先参数化路径；不要在未经明确同意的情况下重组目录结构。
- 保持物理单位一致。代码中的几何参数大多使用 SI 单位，通常以微米数值乘以 `1e-6` 表示。
- 谨慎处理 `T_total` 的符号和 loss 换算。现有分析使用 `loss = -10 * log10(abs(-1.0 * T_total))`；修改公式前需确认物理约定。
- 新增分析代码时，输出文件应放入对应场景目录的 `results/` 子目录。
- 对 notebook 不做大规模重写。若逻辑需要复用，优先抽取到 `.py` 文件。
- 生成图片的文件名应清楚包含扫描参数值。

## 验证方式

轻量验证：

- 仅运行后处理：`python PD-PWB-SMF/analyze_results.py`
- 对 Python 文件做语法检查：`python -m py_compile PD-PWB-SMF/pwb_core.py`

仿真相关验证仅在用户要求且本机 Lumerical 可用时执行：

- 只构建结构、不完整运行：`python PD-PWB-SMF/test_setup.py`
- 单次 FDTD 仿真：`python PD-PWB-SMF/run_single.py`
- 参数扫描：`python PD-PWB-SMF/sweep_r_R.py` 或 `python PD-PWB-SMF/sweep_h_R.py`
