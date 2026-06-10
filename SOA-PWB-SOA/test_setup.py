"""Build and save the SOA-PWB-SOA structure without running FDTD."""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import SOA_DIR, add_lumerical_api_path

add_lumerical_api_path()
import lumapi

from pwb_core import SOAPWBParams, create_pwb_structure_in_fdtd, setup_fdtd_simulation

SAVE_PATH = SOA_DIR / "SOA_PWB_SOA_temp.fsp"

params = SOAPWBParams()
fdtd = lumapi.FDTD()

try:
    positions = create_pwb_structure_in_fdtd(fdtd, params)
    setup_fdtd_simulation(fdtd, params, positions)
    fdtd.save(str(SAVE_PATH))
    print("SOA-PWB-SOA structure saved:", SAVE_PATH)
    print(f"Total structure length: {positions['total'] * 1e6:.1f} um")
finally:
    fdtd.close()
