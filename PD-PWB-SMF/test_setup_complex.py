"""
Build and save the complex planar PWB structure without running FDTD.

Use this script to inspect the generated geometry in Lumerical before running
full simulations or sweeps.
"""

import sys

sys.path.append("D:\\Program Files\\Lumerical\\v241\\api\\python\\")
import lumapi

from pwb_core import (
    PWBParameters,
    generate_pwb_structure_3,
    setup_fdtd_simulation_3,
)


SAVE_PATH = "D:/simulation/Simulation Project/PD-PWB-SMF/temp_complex.fsp"


params = PWBParameters()
fdtd = lumapi.FDTD()

generate_pwb_structure_3(fdtd, params)
setup_fdtd_simulation_3(fdtd, params)
fdtd.save(SAVE_PATH)
print("Complex planar PWB structure saved:", SAVE_PATH)
fdtd.close()
