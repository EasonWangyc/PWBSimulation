"""
test_setup.py — 结构搭建测试
仅生成结构并保存 .fsp 文件，不执行仿真，用于快速验证几何和仿真设置是否正确。
"""

import sys
sys.path.append("D:\\Program Files\\Lumerical\\v241\\api\\python\\")
import lumapi

from pwb_core import (
    PWBParameters,
    generate_pwb_structure_1,
    generate_pwb_structure_2,
    setup_fdtd_simulation_1,
    setup_fdtd_simulation_2,
)

SAVE_PATH = "D:/simulation/Simulation Project/simulation/PD-PWB-SMF/temp.fsp"

# --- 测试 Section 1（bending 段 + 仿真设置）---
params = PWBParameters()
fdtd = lumapi.FDTD()

generate_pwb_structure_1(fdtd, params)
setup_fdtd_simulation_1(fdtd, params)
fdtd.save(SAVE_PATH)
print("Section 1 结构已保存:", SAVE_PATH)
fdtd.close()

# --- 取消注释以测试 Section 2（bending + taper2 段）---
# params2 = PWBParameters()
# fdtd2 = lumapi.FDTD()
# generate_pwb_structure_2(fdtd2, params2)
# setup_fdtd_simulation_2(fdtd2, params2)
# fdtd2.save(SAVE_PATH)
# print("Section 2 结构已保存:", SAVE_PATH)
# fdtd2.close()
