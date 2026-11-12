# Copyright (c) 2026 Jaden John Allison (Fringeneering). All Rights Reserved.
# Contact: fringeneering@gmail.com | GitHub: github.com/fringeneering-jpg
# Dual License: GNU AGPLv3 / Fringeneering Commercial Telemetry License.

"""wurral — Unified Quantum-to-Macro Zero-Parameter Astrodynamics & Telemetry Engine."""

__version__ = "1.0.0"
__author__ = "Jaden John Allison (Fringeneering)"
__license__ = "AGPLv3 / Commercial Telemetry License"

from wurral.constants import (
    S_DEBT, PWC_COEFF, PWC_EXP, CHI_CHIRAL, H0_SI, DRHO_GEOM, A_ENTROPY,
    DGM_1AU_KM3S2, M_PENNING_FESI, R_JOIN_AU, A_JOIN, DIPOLE, SHEAR,
    R_ENTROPY_CROSSOVER_AU, R_JOIN_CROSSOVER_AU, C_SI, C_KMS, G_SI, AU_SI, AU_KM,
)
from wurral.subatomic import run_sector1_audit
from wurral.spin_structure import run_sector2_audit
from wurral.integrator import propagate_solar_system, load_solar_system_state, run_sector3_audit
from wurral.telemetry import (
    dynamic_one_way_two_way_link_budget,
    viking_lander_turnaround_offset,
    run_link_budget_audit,
    run_telemetry_audit,
)
from wurral.etno_wake import run_sector4_audit

__all__ = [
    "__version__", "S_DEBT", "PWC_COEFF", "PWC_EXP", "CHI_CHIRAL", "H0_SI",
    "DRHO_GEOM", "A_ENTROPY", "DGM_1AU_KM3S2", "M_PENNING_FESI", "R_JOIN_AU",
    "A_JOIN", "DIPOLE", "SHEAR", "R_ENTROPY_CROSSOVER_AU", "R_JOIN_CROSSOVER_AU",
    "C_SI", "C_KMS", "G_SI", "AU_SI", "AU_KM",
    "run_sector1_audit", "run_sector2_audit", "propagate_solar_system",
    "load_solar_system_state", "run_sector3_audit",
    "dynamic_one_way_two_way_link_budget", "viking_lander_turnaround_offset",
    "run_link_budget_audit", "run_telemetry_audit", "run_sector4_audit",
]
