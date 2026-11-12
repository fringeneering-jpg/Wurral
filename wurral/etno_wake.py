# Copyright (c) 2026 Jaden John Allison (Fringeneering). All Rights Reserved.
# Contact: fringeneering@gmail.com | GitHub: github.com/fringeneering-jpg
# Dual License: GNU AGPLv3 / Fringeneering Commercial Telemetry License.

"""Sector 4 Engine: Vectorized 4th-order Yoshida Symplectic Zero-Knob Anderson eTNO Wake."""

from typing import Dict, Any
import math
import numpy as np
from wurral.constants import PWC_COEFF, PWC_EXP, R_JOIN_AU, A_JOIN, DIPOLE, SHEAR, AU_SI, YEAR_SEC, GM_SUN_SI

ETNO_CATALOG = [
    {"name": "90377 Sedna",         "a": 506.8,  "e": 0.8495, "i": 11.93, "Omega": 144.40, "omega": 311.29},
    {"name": "2012 VP113 (Biden)",  "a": 261.5,  "e": 0.6922, "i": 24.05, "Omega": 90.80,  "omega": 293.72},
    {"name": "2015 TG387 (Goblin)", "a": 1170.0, "e": 0.9444, "i": 11.67, "Omega": 300.99, "omega": 118.01},
    {"name": "2004 VN112",          "a": 327.4,  "e": 0.8555, "i": 25.55, "Omega": 66.02,  "omega": 327.08},
    {"name": "2007 TG422",          "a": 502.6,  "e": 0.9292, "i": 18.59, "Omega": 112.91, "omega": 285.66},
    {"name": "2013 RF98",           "a": 363.6,  "e": 0.8965, "i": 29.54, "Omega": 67.64,  "omega": 311.76},
]

GM_SUN_AUYR2 = 4.0 * (math.pi ** 2)
GM_NEP_AUYR2 = GM_SUN_AUYR2 * (6.836529e15 / GM_SUN_SI)
_lon_r, _lat_r, _tilt_r = map(math.radians, (71.4, -5.2, 60.2))
WAKE_AXIS = np.array([math.cos(_lat_r)*math.cos(_lon_r), math.cos(_lat_r)*math.sin(_lon_r), math.sin(_lat_r)])
GAL_NORM  = np.array([math.sin(_tilt_r)*math.sin(_lon_r), -math.sin(_tilt_r)*math.cos(_lon_r), math.cos(_tilt_r)])

_W1 = 1.0 / (2.0 - 2.0 ** (1.0 / 3.0))
_W0 = -(2.0 ** (1.0 / 3.0)) * _W1
_YC = [0.5 * _W1, 0.5 * (_W0 + _W1), 0.5 * (_W0 + _W1), 0.5 * _W1]
_YD = [_W1, _W0, _W1, 0.0]


def _elem_to_rv(a: float, e: float, i_d: float, Om_d: float, om_d: float, M_d: float = 180.0):
    i, Om, om, M = map(math.radians, (i_d, Om_d, om_d, M_d))
    E = M
    for _ in range(15):
        E -= (E - e * math.sin(E) - M) / (1.0 - e * math.cos(E))
    xo, yo = a * (math.cos(E) - e), a * math.sqrt(1.0 - e**2) * math.sin(E)
    n = math.sqrt(GM_SUN_AUYR2 / (a**3)); r = a * (1.0 - e * math.cos(E))
    vxo, vyo = -(a**2 * n / r) * math.sin(E), (a**2 * n / r) * math.sqrt(1.0 - e**2) * math.cos(E)
    cO, sO, cw, sw, ci, si = math.cos(Om), math.sin(Om), math.cos(om), math.sin(om), math.cos(i), math.sin(i)
    P = np.array([cO*cw - sO*sw*ci, sO*cw + cO*sw*ci, sw*si])
    Q = np.array([-cO*sw - sO*cw*ci, -sO*sw + cO*cw*ci, cw*si])
    return P * xo + Q * yo, P * vxo + Q * vyo


def _rv_to_q_a(pos: np.ndarray, vel: np.ndarray):
    r = float(np.linalg.norm(pos))
    en = 0.5 * float(np.dot(vel, vel)) - GM_SUN_AUYR2 / r
    if en >= 0:
        return float("inf"), float("nan")
    a = -GM_SUN_AUYR2 / (2.0 * en)
    h = np.cross(pos, vel)
    e = float(np.linalg.norm(np.cross(vel, h) / GM_SUN_AUYR2 - pos / r))
    return a, a * (1.0 - e)


def _etno_accel(pos: np.ndarray) -> np.ndarray:
    r2 = np.einsum("ij,ij->i", pos, pos)
    r = np.sqrt(r2)
    r_hat = pos / r[:, None]
    acc = -GM_SUN_AUYR2 * pos / (r2 * r)[:, None]
    a_nep = 30.07
    mask = r > a_nep
    if np.any(mask):
        z2 = pos[mask, 2] ** 2
        qf = 0.75 * GM_NEP_AUYR2 * (a_nep**2) / (r[mask]**5)
        acc[mask, 0] -= qf * pos[mask, 0] * (1.0 - 5.0 * z2 / r2[mask])
        acc[mask, 1] -= qf * pos[mask, 1] * (1.0 - 5.0 * z2 / r2[mask])
        acc[mask, 2] -= qf * pos[mask, 2] * (3.0 - 5.0 * z2 / r2[mask])
    rat = r / R_JOIN_AU
    f_join = PWC_COEFF * ((rat**2) / (1.0 + rat**2)) / (1.0 + rat**(-PWC_EXP))
    a_c = A_JOIN * (YEAR_SEC**2) / AU_SI
    cos_t = np.dot(pos, WAKE_AXIS) / r
    acc_rad = -(a_c * f_join)[:, None] * r_hat
    acc_dip = -(DIPOLE * a_c * f_join)[:, None] * (WAKE_AXIS[None, :] - 0.5 * cos_t[:, None] * r_hat)
    tang = np.cross(GAL_NORM[None, :], r_hat)
    tn = np.linalg.norm(tang, axis=1, keepdims=True)
    tang = np.where(tn > 1e-12, tang / np.maximum(tn, 1e-12), 0.0)
    acc_shr = (SHEAR * a_c * f_join * (1.0 - cos_t))[:, None] * tang
    return acc + acc_rad + acc_dip + acc_shr


def run_sector4_audit(sim_years: float = 400000.0, dt_yr: float = 25.0) -> Dict[str, Any]:
    p_l, v_l = [], []
    for idx, obj in enumerate(ETNO_CATALOG):
        off = 65.0 if idx % 2 == 0 else -60.0
        p0, v0 = _elem_to_rv(obj["a"], obj["e"], obj["i"], obj["Omega"], (obj["omega"] + off) % 360.0, 180.0)
        p_l.append(p0); v_l.append(v0)
    p_s, v_s = _elem_to_rv(506.8, 1.0 - 30.5 / 506.8, 11.93, 144.40, 311.29, 0.0)
    p_l.append(p_s); v_l.append(v_s)

    pw, vw = np.array(p_l), np.array(v_l)
    steps = int(sim_years / dt_yr)
    for _ in range(steps):
        for k in range(4):
            pw += _YC[k] * vw * dt_yr
            if _YD[k] != 0.0:
                vw += _YD[k] * _etno_accel(pw) * dt_yr

    bound_cnt = sum(1 for i in range(6) if math.isfinite(_rv_to_q_a(pw[i], vw[i])[0]))
    _, q_sedna_final = _rv_to_q_a(pw[6], vw[6])
    return {
        "sector": "Sector 4: Zero-Knob Symplectic Anderson eTNO Wake Engine",
        "passed": bound_cnt == 6 and q_sedna_final >= 65.0,
        "sim_years": sim_years,
        "orbits_bound": f"{bound_cnt}/6",
        "proto_sedna_q_init_au": 30.5,
        "proto_sedna_q_final_au": q_sedna_final,
    }
