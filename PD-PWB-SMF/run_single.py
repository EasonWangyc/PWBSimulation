"""
run_single.py — 单次完整仿真
搭建 Section 1 结构，运行仿真，可视化电场分布并打印 T_total。
"""

import sys
sys.path.append("D:\\Program Files\\Lumerical\\v241\\api\\python\\")
import lumapi

from pwb_core import (
    PWBParameters,
    generate_pwb_structure_1,
    setup_fdtd_simulation_1,
    visualize_and_save_results_1,
)

SAVE_PATH = "D:/simulation/Simulation Project/simulation/PD-PWB-SMF/test/test.fsp"

params = PWBParameters()
fdtd = lumapi.FDTD()

generate_pwb_structure_1(fdtd, params)
setup_fdtd_simulation_1(fdtd, params)
fdtd.save(SAVE_PATH)

fdtd.run()

T_total = visualize_and_save_results_1(fdtd, params)
print(f"T_total = {T_total}")

fdtd.close()
