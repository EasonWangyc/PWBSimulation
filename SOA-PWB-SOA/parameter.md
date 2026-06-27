# SOA-PWB-SOA 参数说明

本目录用于 片上SOA输出端 → PWB → 外置SOA输入端 的耦合仿真。

两端 SOA 波导结构由 `SOA_base.fsp` 定义（脊形波导，含200nm有源区）。
PWB 以 `addcircle()` + ellipsoid 构建，截面为圆或椭圆。

提供：
- `test_setup.py`：只生成并保存结构，不运行 FDTD。
- `run_single.py`：生成结构、运行一次仿真并绘图。
- `pwb_core.py`：核心模块，包含参数、几何、FDTD设置、FDE分析和可视化。

## 0. 路径配置

脚本中的路径统一从根目录 `config/sim_config.py` 获取：
- `SOA_DIR`：当前场景目录。
- `SOA_RESULTS_DIR`：结果输出目录。
- `SOA_BASE_FSP`：包含两端SOA波导的基座文件。
- `MATERIAL_DB`：材料库文件。
- `add_lumerical_api_path()`：将 Lumerical Python API 路径加入 `sys.path`。

## 1. 几何结构

结构沿 x 轴直线排列（无弯曲），PWB 由 addcircle() 段组成，总长 250μm：

```
x=0 ──────────────────────────── x=250μm
  [SOA_out, .fsp]  Taper-1    PWB直段    Taper-2  [SOA_in, .fsp]
                   (扩半径)   (等半径)   (缩半径)
```

PWB 截面为圆形（`use_ellipsoid=False`）或椭圆形（`use_ellipsoid=True`）。

核心几何参数（单位：SI，以 `* 1e-6` 表示微米）：

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `total_length` | PWB 总长度（含 taper） | 250e-6 |
| `taper1_length` | Taper-1 长度 | 30e-6 |
| `taper2_length` | Taper-2 长度 | 30e-6 |

PWB 直段长度 = `total_length - taper1_length - taper2_length`（自动计算）。

### 半径剖面

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `r_in` | Taper-1 起始半径 | 1.0e-6 |
| `r_pwb` | PWB 中间段半径 | 1.5e-6 |
| `r_out` | Taper-2 终止半径 | 1.0e-6 |
| `use_ellipsoid` | 是否使用椭圆截面 | True |
| `r_in_2` | Taper-1 起始第二轴半径 | 2.0e-6 |
| `r_pwb_2` | PWB 中间段第二轴半径 | 1.5e-6 |
| `r_out_2` | Taper-2 终止第二轴半径 | 2.0e-6 |

半径沿 x 轴的演化：
- Taper-1：r_in → r_pwb（线性，第一个1/3段）
- PWB 直段：r_pwb 恒定（中间1/3段）
- Taper-2：r_pwb → r_out（线性，最后1/3段）

### 两端 SOA 波导参数（由 .fsp 定义）

| 参数 | 值 |
|------|-----|
| 片上 SOA 输出端宽度 | 5 μm（脊形，由2μm渐变而来） |
| 片上 SOA 输出端高度 | 2 μm |
| 有源区厚度 | 200 nm |
| 外置 SOA 输入宽度 | 5 μm（两端）/ 2 μm（中间） |
| 实际 PWB 打线尺寸 | 5.5 μm × 1.8 μm（供参考） |

## 2. 材料参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `material_pwb` | PWB聚合物材料 | "Vancore B" |

SOA 材料在 .fsp 中定义。

## 3. 仿真参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `wavelength` | 中心波长 | 1.55e-6 |
| `mesh_accuracy` | 网格精度 (1-8) | 2 |
| `simulation_time` | 仿真时间 [s] | 2000e-15 |
| `curve_points` | 路径采样点数 | 200 |

## 4. 监视器与模式扩展

FDTD 中共设 4 个 X-normal 功率监视器 + 2 个模式扩展：

| 监视器 | 位置 | 用途 |
|--------|------|------|
| `input_monitor` | Taper-1 中段 | 确认源注入 |
| `pwb_in_monitor` | Taper-1 末端 | PWB 入口截面 |
| `pwb_out_monitor` | Taper-2 始端 | PWB 出口截面 |
| `output_monitor` | Taper-2 末端 | 输出截面 |
| `mode_exp_pwb_in` | Taper-1 末端 | 基模传输效率（过 Taper-1 后） |
| `mode_exp_output` | Taper-2 末端 | 总基模传输效率 |
| `transmission_monitor` | Y-normal 全视场 | 侧面传播分布 |

## 5. FDE 模式分析（辅助）

使用 Lumerical MODE Solutions 的 FDE 求解器计算三个截面的基模：
- SOA 输出端截面（5μm × 2μm 脊形）
- PWB 中间截面（r_pwb × r_pwb_2 近似）
- 外置 SOA 输入端截面

**主要优化手段是 FDTD + mode expansion monitor 直读 T_forward，FDE 仅作 sanity check。**

## 6. 运行方式

只生成结构：
```powershell
& "D:\Program Files\Lumerical\python\python.exe" SOA-PWB-SOA\test_setup.py
```

单次仿真：
```powershell
& "D:\Program Files\Lumerical\python\python.exe" SOA-PWB-SOA\run_single.py
```
