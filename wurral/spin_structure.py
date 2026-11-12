# Copyright (c) 2026 Jaden John Allison (Fringeneering). All Rights Reserved.
# Contact: fringeneering@gmail.com | GitHub: github.com/fringeneering-jpg
# Dual License: GNU AGPLv3 / Fringeneering Commercial Telemetry License.

"""Sector 2 Audit: Planetary Rotation Lock, Lunar Viscosity & Jiba LEO Boundary Decay."""

from typing import Dict, Any, List
import math
from wurral.constants import (
    PWC_EXP, M_PENNING_FESI, JIBA_DR1_M, JIBA_R1_KM, JIBA_ALPHA, R_MOON_ORBIT_KM,
)


def core_steepening_factor(rho_c: float, rho_mean: float) -> float:
    return float(1.0 + (PWC_EXP / 3.0) * math.log(rho_c / rho_mean))


def flattening_relief_factor(r_pol_m: float, r_eq_m: float, j2: float) -> float:
    return float(((r_pol_m / r_eq_m) ** (1.0 / 6.0)) * math.sqrt(1.0 + 1.5 * j2))


def planetary_rotation_lock() -> Dict[str, Any]:
    planets_data = [
        ("Earth",   1.00000011, 6.3781e6, 6.3568e6, 1.08263e-3, 3.986004e14, 23.9345, 12.5, 5.514, (12.0/209.0)*(1.00000011**(1.0/3.0)), True),
        ("Mars",    1.52366231, 3.3962e6, 3.3762e6, 1.9555e-3,  4.282838e13, 24.6229,  6.2, 3.933, (12.0/209.0)*(1.52366231**(1.0/3.0)), True),
        ("Jupiter", 5.20260,    7.1492e7, 6.6854e7, 1.46966e-2, 1.266865e17,  9.9250, 15.0, 1.326, 11.0/38.0, False),
        ("Saturn",  9.55491,    6.0268e7, 5.4364e7, 1.62907e-2, 3.793119e16, 10.5730, 12.0, 0.687,  7.0/19.0, False),
        ("Uranus",  19.21845,   2.5559e7, 2.4973e7, 3.5107e-3,  5.793940e15, 17.2400,  8.0, 1.270, (3.0/19.0)*(19.0/18.0), False),
        ("Neptune", 30.11039,   2.4764e7, 2.4341e7, 3.4084e-3,  6.836529e15, 16.1100,  9.0, 1.638,  3.0/19.0, False),
    ]
    results: Dict[str, Any] = {}
    errs: List[float] = []
    for name, a_au, req, rpol, j2, gm_p, p_obs, rhoc, rhom, q_gate, is_rocky in planets_data:
        cs = core_steepening_factor(rhoc, rhom)
        fr = flattening_relief_factor(rpol, req, j2)
        pm = M_PENNING_FESI if is_rocky else 1.0
        w_eq = math.sqrt(gm_p / (req**3)) * pm * cs * fr
        t_sh = (2.0 * math.pi / w_eq) / 3600.0
        p_pred = t_sh / q_gate
        signed_pct = (p_pred / p_obs - 1.0) * 100.0
        errs.append(abs(signed_pct))
        results[name] = {
            "t_sh_hr": t_sh, "obs_hr": p_obs, "pred_hr": p_pred,
            "signed_error_pct": signed_pct, "abs_error_pct": abs(signed_pct),
        }
    mean_res = float(sum(errs) / len(errs))
    return {
        "planets": results,
        "mean_residual_pct": mean_res,
        "mercury_3_2_ratchet_hr": (2.0 / 3.0) * 87.969 * 24.0,
        "verified": abs(mean_res - 1.20) < 0.05,
    }


def lunar_libration_viscosity() -> Dict[str, Any]:
    tau_rest = 1.5 * (3.986004418e14 / (384399000.0**3)) * 1.98e31 * math.sin(2.0 * math.radians(100.0 / 3600.0))
    eta_pa_s = tau_rest / ((1022.0 / 384399000.0) * (7.0e5**3))
    return {"restoring_torque_nm": tau_rest, "viscosity_pa_s": eta_pa_s, "verified": abs(eta_pa_s - 2.216e5) < 500.0}


def jiba_leo_dipole_decay() -> Dict[str, Any]:
    decayed_mm = JIBA_DR1_M * ((JIBA_R1_KM / R_MOON_ORBIT_KM) ** JIBA_ALPHA) * 1000.0
    return {"decayed_amplitude_mm": decayed_mm, "verified": abs(decayed_mm - 2.81) < 0.02}


def run_sector2_audit() -> Dict[str, Any]:
    sp = planetary_rotation_lock()
    lv = lunar_libration_viscosity()
    jd = jiba_leo_dipole_decay()
    return {
        "sector": "Sector 2: Planetary Spin Flattening Lock & Boundary Decay",
        "passed": sp["verified"] and lv["verified"] and jd["verified"],
        "spin_lock": sp, "lunar_viscosity": lv, "leo_dipole_decay": jd,
    }
