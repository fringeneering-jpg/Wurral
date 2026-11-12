# Copyright (c) 2026 Jaden John Allison (Fringeneering). All Rights Reserved.
# Contact: fringeneering@gmail.com | GitHub: github.com/fringeneering-jpg
# Dual License: GNU AGPLv3 / Fringeneering Commercial Telemetry License.

"""Strictly locked zero-parameter physical, constitutive, and derived constants."""

import math

S_DEBT: float = 1.0 / 19.0
PWC_COEFF: float = 18.0 / 19.0
PWC_EXP: float = 2.0 / 33.0
CHI_CHIRAL: float = 12.0 / 11.0

C_SI: float = 299792458.0
C_KMS: float = C_SI / 1000.0
G_SI: float = 6.67430e-11
AU_SI: float = 149597870700.0
AU_KM: float = AU_SI / 1000.0
DAY_SEC: float = 86400.0
YEAR_DAYS: float = 365.25
YEAR_SEC: float = YEAR_DAYS * DAY_SEC
JULIAN_CENTURY_DAYS: float = 36525.0
J2000_JD: float = 2451545.0

H0_SI: float = 73.935 / 3.08567758149e19
DRHO_GEOM: float = (1.5 * math.gamma(1.0 / 3.0) - 3.0) / (4.0 * math.pi / 3.0)
A_ENTROPY: float = C_SI * H0_SI * (1.0 + DRHO_GEOM)
DGM_1AU_KM3S2: float = A_ENTROPY * (AU_SI ** 2) / 1.0e9
DGM_1AU_SI: float = DGM_1AU_KM3S2 * 1.0e9

M_PENNING_FESI: float = math.sqrt((1.0 / (1.0 - 0.00733)) / (1.0 / (1.0 - 0.00088)))

M_E_KEV: float = 510.998950
ALPHA_EM: float = 1.0 / 137.035999084
M_P_ISO_MEV: float = 938.27208816
M_P_SEAL_MEV: float = 937.190397372
Q_WELD_KEV: float = 91.846216503
RHO_MAX_GW: float = 1.3040143e15
R_P0_FM: float = 0.841400

GM_SUN_KM3S2: float = 132712440041.93938
GM_SUN_SI: float = GM_SUN_KM3S2 * 1.0e9
SUN_REQ_KM: float = 695700.0
OMEGA_SPIN_SUN: float = 2.8653e-6
K_SUN: float = (2.0 * OMEGA_SPIN_SUN * (SUN_REQ_KM * 1000.0)) / C_SI
SUN_J2: float = 2.20e-7
SUN_POLE_RA_DEG: float = 286.13
SUN_POLE_DEC_DEG: float = 63.87

EARTH_J2: float = 1.08262668e-3
EARTH_REQ_KM: float = 6378.137

JIBA_DR1_M: float = 100.55
JIBA_R1_KM: float = 7190.137
JIBA_ALPHA: float = 2.635
R_MOON_ORBIT_KM: float = 384399.0

R_JOIN_AU: float = 122.0
A_JOIN: float = (3.0 / 19.0) * A_ENTROPY
DIPOLE: float = 7.0 / 19.0
SHEAR: float = (3.0 / 19.0) * math.cos(math.radians(30.1))

R_ENTROPY_CROSSOVER_AU: float = math.sqrt(GM_SUN_SI / A_ENTROPY) / AU_SI
R_JOIN_CROSSOVER_AU: float = math.sqrt(GM_SUN_SI / A_JOIN) / AU_SI
