# SOA-PWB-SOA 参数说明

本目录用于 片上SOA输出端 → PWB → 外置SOA输入端 的耦合仿真。

提供：
- `test_setup.py`：只生成并保存结构，不运行 FDTD。
- `run_single.py`：生成结构、运行一次仿真并绘图。
- `pwb_core.py`：核心模块，包含参数、几何、FDTD设置、FDE分析和可视化。

## 0. 路径配置

脚本中的路径统一从根目录 `config/sim_config.py` 获取：
- `SOA_DIR`：当前场景目录。
- `SOA_RESULTS_DIR`：结果输出目录。
- `MATERIAL_DB`：材料库文件。
- `add_lumerical_api_path()`：将 Lumerical Python API 路径加入 `sys.path`。

## 1. 几何结构

结构沿 x 轴直线排列（无弯曲），包含五段：

```
x=0 ─────── x=L1 ─────── x=L2 ─────── x=L3 ─────── x=L4
  SOA输出     Taper-1      PWB直段     Taper-2     外置SOA输入
  (InP脊形)   (模式扩展)   (聚合物)    (模式压缩)   (待定结构)
```

核心几何参数（单位：SI，以 `* 1e-6` 表示微米）：

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `wg1_width` | 片上SOA输出波导宽度 | 待提供 |
| `wg1_height` | 片上SOA输出波导高度 | 待提供 |
| `wg1_active_thickness` | 片上SOA有源区厚度 | 待提供 |
| `wg1_length` | 片上SOA输出段长度 | 10e-6 |
| `taper1_w_in` | Taper-1输入端宽度 | = wg1_width |
| `taper1_w_out` | Taper-1输出端宽度 | 待提供 |
| `taper1_h_in` | Taper-1输入端高度 | = wg1_height |
| `taper1_h_out` | Taper-1输出端高度 | 待提供 |
| `taper1_length` | Taper-1长度 | 待提供 |
| `pwb_length` | PWB中间直段长度 | 待提供 |
| `taper2_w_in` | Taper-2输入端宽度 | = taper1_w_out |
| `taper2_w_out` | Taper-2输出端宽度 | 待提供 |
| `taper2_h_in` | Taper-2输入端高度 | = taper1_h_out |
| `taper2_h_out` | Taper-2输出端高度 | 待提供 |
| `taper2_length` | Taper-2长度 | 待提供 |
| `wg2_width` | 外置SOA输入波导宽度 | 待提供 |
| `wg2_height` | 外置SOA输入波导高度 | 待提供 |
| `wg2_length` | 外置SOA输入段长度 | 10e-6 |

## 2. 材料参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `material_soa` | SOA波导材料 | "InP" (需在材料库中定义) |
| `material_pwb` | PWB聚合物材料 | "Vancore B" |
| `material_clad` | 包层材料 | 空气 (index=1.0) |
| `pwb_index` | PWB折射率 | 待提供 |

## 3. 仿真参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `wavelength` | 中心波长 | 1.55e-6 |
| `mesh_accuracy` | 网格精度 (1-8) | 2 |
| `simulation_time` | 仿真时间 [s] | 2000e-15 |

## 4. FDE 模式分析

使用 Lumerical MODE Solutions 的 FDE (Finite-Difference Eigenmode) 求解器，
计算三个截面的基模场分布和模式重叠因子：

- SOA 输出端截面 → 基模 E₁(x,y)
- PWB 中间截面 → 基模 E₂(x,y)
- 外置 SOA 输入端截面 → 基模 E₃(x,y)

模式重叠因子 η：
```
η = |∫∫ E₁(x,y)·E₂*(x,y) dxdy|² / [∫∫|E₁|²dxdy · ∫∫|E₂|²dxdy]
```

## 5. 运行方式

只生成结构：
```powershell
& "D:\Program Files\Lumerical\python\python.exe" SOA-PWB-SOA\test_setup.py
```

单次仿真：
```powershell
& "D:\Program Files\Lumerical\python\python.exe" SOA-PWB-SOA\run_single.py
```
