"""Shared path configuration for the PWB simulation workspace."""

import os
import sys
from pathlib import Path


PROJECT_ROOT = Path(os.environ.get("SIM_PROJECT_ROOT", Path(__file__).resolve().parent.parent)).resolve()

LUMERICAL_API_PATH = Path(
    os.environ.get("LUMERICAL_API_PATH", r"D:\Program Files\Lumerical\api\python")
)
MATERIAL_DB = Path(os.environ.get("SIM_MATERIAL_DB", PROJECT_ROOT / "config" / "database.mdf"))

PD_DIR = PROJECT_ROOT / "PD-PWB-SMF"
if str(PD_DIR) not in sys.path:
    sys.path.insert(0, str(PD_DIR))
LD_DIR = PROJECT_ROOT / "LD-PWB-SMF"
LNOI_DIR = PROJECT_ROOT / "LNOI-PWB-SMF"
SOA_DIR = PROJECT_ROOT / "SOA-PWB-SOA"
SOA_RESULTS_DIR = SOA_DIR / "results"
SOA_BASE_FSP = SOA_DIR / "SOA_base.fsp"

PD_RESULTS_DIR = PD_DIR / "results"
LD_RESULTS_DIR = LD_DIR / "results"
LNOI_RESULTS_DIR = LNOI_DIR / "results"

PD_SMF_FSP = PD_DIR / "SMF.fsp"
PD_SECTION1_FSP = PD_DIR / "Section1.fsp"

LD_BASE_FSP = LD_DIR / "LD.fsp"
LD_SOURCE_DATASET = LD_DIR / "Section2 output.mat"

LNOI_BASE_FSP = LNOI_DIR / "LNOI.fsp"


def add_lumerical_api_path():
    """Add the configured Lumerical Python API path to sys.path."""
    api_path = str(LUMERICAL_API_PATH)
    if api_path not in sys.path:
        sys.path.append(api_path)
    return LUMERICAL_API_PATH


def project_path(*parts):
    """Return a path under the project root."""
    return PROJECT_ROOT.joinpath(*parts)
