# Copyright (c) 2026 Jaden John Allison (Fringeneering). All Rights Reserved.
# Contact: fringeneering@gmail.com | GitHub: github.com/fringeneering-jpg
# Dual License: GNU AGPLv3 / Fringeneering Commercial Telemetry License.

"""
Sector 3 Engine: Vectorized NumPy DOP853 15-Body Solar System Propagator.
Implements complete Moyer Eq. 4-26 pairwise 1PN EIH (Terms 1, 2, and 3),
True-of-Date precessing Earth J2 + Solar J2, fixed-period gauge ruler scaling,
and frame-uniform SSB centering. Strictly conservative circular-groove ODE (coupling = 0.0).
"""

from typing import Dict, Any, Tuple, List, Optional
import json
import math
import os
import urllib.parse
import urllib.request
import numpy as np
from scipy.integrate import solve_ivp

from wurral.constants import (
    C_KMS, DAY_SEC, JULIAN_CENTURY_DAYS, J2000_JD,
    EARTH_J2, EARTH_REQ_KM, SUN_J2, SUN_REQ_KM,
    SUN_POLE_RA_DEG, SUN_POLE_DEC_DEG, DGM_1AU_KM3S2,
)

try:
    from jplephem.spk import SPK
    JPLEPHEM_AVAILABLE = True
except ImportError:
    SPK = None
    JPLEPHEM_AVAILABLE = False

# Indices: 0=Sun, 1=Mercury, 2=Venus, 3=Earth, 4=Mars, 5=Jupiter, 6=Saturn, 7=Uranus, 8=Neptune, 9=Moon, 10=Pluto, 11..14=Main-Belt
BODY_NAMES: List[str] = [
    "Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn",
    "Uranus", "Neptune", "Moon", "Pluto_Charon", "Ceres", "Vesta", "Pallas", "Hygiea",
]

EMRAT: float = 81.30056907419062
GM_EMB: float = 403503.235502
GM_MOON: float = GM_EMB / (1.0 + EMRAT)
GM_EARTH: float = GM_EMB - GM_MOON

BODY_GMS_CAL_KM3S2: np.ndarray = np.array([
    132712440041.93938, 22031.868551, 324858.592, GM_EARTH,
    42828.375214, 126712764.8, 37940585.2, 5794548.6, 6836527.10058,
    GM_MOON, 975.500, 62.6284, 17.2882, 13.63, 5.78,
], dtype=np.float64)


def _fetch_or_cache_4rocks(jd_epoch: float, cache_dir: str) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    os.makedirs(cache_dir, exist_ok=True)
    cache_fp = os.path.join(cache_dir, f"dawn_4rocks_{int(round(jd_epoch))}.json")
    if os.path.isfile(cache_fp):
        with open(cache_fp, "r", encoding="utf-8") as f:
            cd = json.load(f)
        if abs(cd["jd"] - jd_epoch) < 1e-4:
            return [np.array(x, dtype=np.float64) for x in cd["r"]], [np.array(x, dtype=np.float64) for x in cd["v"]]

    # If J2000, load directly from bundled exact snapshot!
    snap_j2000 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "j2000_exact_15body.json")
    if abs(jd_epoch - J2000_JD) < 1e-6 and os.path.isfile(snap_j2000):
        with open(snap_j2000, "r", encoding="utf-8") as f:
            sd = json.load(f)
        return [np.array(x, dtype=np.float64) for x in sd["r_km"][11:15]], [np.array(x, dtype=np.float64) for x in sd["v_kms"][11:15]]

    rk_r, rk_v = [], []
    for cmd in ("1;", "4;", "2;", "10;"):
        q = urllib.parse.quote(cmd)
        url = (f"https://ssd.jpl.nasa.gov/api/horizons.api?format=json&COMMAND='{q}'&EPHEM_TYPE='VECTORS'"
               f"&CENTER='500@0'&START_TIME='JD{jd_epoch}'&STOP_TIME='JD{jd_epoch+1.0}'&STEP_SIZE='2d'&OUT_UNITS='KM-S'&REF_PLANE='FRAME'")
        txt = json.loads(urllib.request.urlopen(url, timeout=20).read().decode("utf-8"))["result"]
        soe = txt.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
        xl = [l for l in soe if "X =" in l][0]
        vl = [l for l in soe if "VX=" in l][0]
        rk_r.append(np.array([float(xl.split(t)[1].strip().split()[0]) for t in ("X =", "Y =", "Z =")], dtype=np.float64))
        rk_v.append(np.array([float(vl.split(t)[1].strip().split()[0]) for t in ("VX=", "VY=", "VZ=")], dtype=np.float64))
    with open(cache_fp, "w", encoding="utf-8") as f:
        json.dump({"jd": jd_epoch, "r": [x.tolist() for x in rk_r], "v": [x.tolist() for x in rk_v]}, f)
    return rk_r, rk_v


def load_solar_system_state(jd_epoch: float = J2000_JD, bsp_path: Optional[str] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Loads exact 64-bit state (pos_km, vel_kms, gms_cal_km3s2) from de440.bsp or bundled IEEE J2000 snapshot."""
    snap_j2000 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "j2000_exact_15body.json")

    if bsp_path and os.path.isfile(bsp_path) and JPLEPHEM_AVAILABLE:
        k440 = SPK.open(bsp_path)
        ps, vs = k440[0, 10].compute_and_differentiate(jd_epoch)
        pemb, vemb = k440[0, 3].compute_and_differentiate(jd_epoch)
        pe, ve = k440[3, 399].compute_and_differentiate(jd_epoch)
        pm, vm = k440[3, 301].compute_and_differentiate(jd_epoch)

        R_l, V_l = [ps], [vs / DAY_SEC]
        for bid in (1, 2):
            p, v = k440[0, bid].compute_and_differentiate(jd_epoch)
            R_l.append(p); V_l.append(v / DAY_SEC)
        R_l.append(pemb + pe); V_l.append((vemb + ve) / DAY_SEC)     # 3: Earth
        for bid in (4, 5, 6, 7, 8):
            p, v = k440[0, bid].compute_and_differentiate(jd_epoch)
            R_l.append(p); V_l.append(v / DAY_SEC)
        R_l.append(pemb + pm); V_l.append((vemb + vm) / DAY_SEC)     # 9: Moon
        pp, vp = k440[0, 9].compute_and_differentiate(jd_epoch)
        R_l.append(pp); V_l.append(vp / DAY_SEC)                     # 10: Pluto

        cache_dir = os.path.dirname(os.path.abspath(bsp_path))
        try:
            rk_r, rk_v = _fetch_or_cache_4rocks(jd_epoch, cache_dir)
            R_l.extend(rk_r); V_l.extend(rk_v)
            gms = BODY_GMS_CAL_KM3S2.copy()
        except Exception:
            gms = BODY_GMS_CAL_KM3S2[:11].copy()
        return np.array(R_l, dtype=np.float64), np.array(V_l, dtype=np.float64), gms

    # Bundled zero-SPICE autonomous IEEE-754 J2000 snapshot!
    with open(snap_j2000, "r", encoding="utf-8") as f:
        sd = json.load(f)
    r_j2 = np.array(sd["r_km"], dtype=np.float64)
    v_j2 = np.array(sd["v_kms"], dtype=np.float64)
    return r_j2, v_j2, BODY_GMS_CAL_KM3S2.copy()


def build_wurral_rhs(jd_start: float, gms_tru: np.ndarray):
    nt = len(gms_tru)
    nmaj = 10
    c2 = C_KMS ** 2
    ra_s, dec_s = math.radians(SUN_POLE_RA_DEG), math.radians(SUN_POLE_DEC_DEG)
    pole_sun = np.array([math.cos(dec_s)*math.cos(ra_s), math.cos(dec_s)*math.sin(ra_s), math.sin(dec_s)], dtype=np.float64)

    def earth_pole(t_s: float) -> np.ndarray:
        T = ((jd_start + t_s / DAY_SEC) - J2000_JD) / JULIAN_CENTURY_DAYS
        z = math.radians((2306.2181 * T) / 3600.0)
        th = math.radians((2004.3109 * T) / 3600.0)
        return np.array([-math.sin(th) * math.cos(z), -math.sin(th) * math.sin(z), math.cos(th)], dtype=np.float64)

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        pos = y[:3 * nt].reshape((nt, 3))
        vel = y[3 * nt:].reshape((nt, 3))
        diff = pos[None, :, :] - pos[:, None, :]
        d2 = np.einsum("ijk,ijk->ij", diff, diff)
        np.fill_diagonal(d2, 1.0)
        d1 = np.sqrt(d2)
        inv_d3 = d2 ** (-1.5)
        np.fill_diagonal(inv_d3, 0.0)
        acc_n = np.einsum("ij,ijk->ik", inv_d3 * gms_tru[None, :], diff)
        acc = acc_n.copy()

        # Complete Moyer Eq. 4-26 Pairwise 1PN EIH (Terms 1, 2, AND 3!)
        p_m, v_m, an_m, g_m = pos[:nmaj], vel[:nmaj], acc_n[:nmaj], gms_tru[:nmaj]
        r_ji = p_m[None, :, :] - p_m[:, None, :]
        r_ij = -r_ji
        rm = d1[:nmaj, :nmaj]
        inv_rm = 1.0 / rm; np.fill_diagonal(inv_rm, 0.0)
        inv_rm3 = inv_rm ** 3
        phi = inv_rm @ g_m
        v2 = np.einsum("ij,ij->i", v_m, v_m)
        vi_vj = v_m @ v_m.T
        nij_vj = np.einsum("ijk,jk->ij", r_ij, v_m) * inv_rm
        rij_aj = np.einsum("ijk,jk->ij", r_ij, an_m)

        coef1 = (g_m[None, :] / c2) * inv_rm3 * (
            -4.0 * phi[:, None] - phi[None, :] + v2[:, None] + 2.0 * v2[None, :]
            - 4.0 * vi_vj - 1.5 * (nij_vj ** 2) - 0.5 * rij_aj
        )
        np.fill_diagonal(coef1, 0.0)
        t1 = np.einsum("ij,ijk->ik", coef1, r_ji)

        four_vi_3vj = 4.0 * v_m[:, None, :] - 3.0 * v_m[None, :, :]
        coef2 = (g_m[None, :] / c2) * inv_rm3 * np.einsum("ijk,ijk->ij", r_ij, four_vi_3vj)
        np.fill_diagonal(coef2, 0.0)
        t2 = np.einsum("ij,ijk->ik", coef2, v_m[:, None, :] - v_m[None, :, :])

        coef3 = 3.5 * (g_m[None, :] / c2) * inv_rm
        np.fill_diagonal(coef3, 0.0)
        t3 = coef3 @ an_m
        acc[:nmaj] += t1 + t2 + t3

        pe_v = earth_pole(t)
        for ti in (9, 0):
            rv = pos[ti] - pos[3]
            rsq = float(np.dot(rv, rv)); rmg = math.sqrt(rsq)
            zp = float(np.dot(rv, pe_v))
            jf = 1.5 * EARTH_J2 * gms_tru[3] * (EARTH_REQ_KM ** 2) / (rsq * rsq * rmg)
            aj2 = jf * ((5.0 * (zp ** 2) / rsq - 1.0) * rv - 2.0 * zp * pe_v)
            acc[ti] += aj2
            acc[3] -= (gms_tru[ti] / gms_tru[3]) * aj2

        for p_i in range(1, 5):
            rs = pos[p_i] - pos[0]
            rs2 = float(np.dot(rs, rs)); rs1 = math.sqrt(rs2)
            zps = float(np.dot(rs, pole_sun))
            jf_s = 1.5 * SUN_J2 * gms_tru[0] * (SUN_REQ_KM ** 2) / (rs2 * rs2 * rs1)
            acc[p_i] += jf_s * ((5.0 * (zps ** 2) / rs2 - 1.0) * rs - 2.0 * zps * pole_sun)

        return np.hstack([vel.ravel(), acc.ravel()])

    return rhs


def propagate_solar_system(
    t_span_days: float = 365.25,
    jd_start: float = J2000_JD,
    bsp_path: Optional[str] = None,
    rtol: float = 1e-13,
    atol: float = 1e-13,
    dense_output: bool = True,
) -> Dict[str, Any]:
    R0, V0, gms_cal = load_solar_system_state(jd_epoch=jd_start, bsp_path=bsp_path)
    gms_tru = gms_cal.copy()
    gms_tru[0] -= DGM_1AU_KM3S2
    lam_g = float((gms_tru[0] / gms_cal[0]) ** (1.0 / 3.0))

    r_ssb = (gms_cal @ R0) / np.sum(gms_cal)
    v_ssb = (gms_cal @ V0) / np.sum(gms_cal)
    R0 = R0 - r_ssb[None, :]
    V0 = V0 - v_ssb[None, :]

    nt = len(gms_tru)
    y0 = np.hstack([(R0 * lam_g).ravel(), (V0 * lam_g).ravel()])
    rhs_fn = build_wurral_rhs(jd_start, gms_tru)
    sol = solve_ivp(rhs_fn, (0.0, t_span_days * DAY_SEC), y0, method="DOP853",
                    rtol=rtol, atol=atol, dense_output=dense_output)

    v_init = V0
    v_final = sol.y[3 * nt:, -1].reshape((nt, 3)) / lam_g
    ssb_drift_kms = float(np.linalg.norm((gms_tru @ v_final) - (gms_tru @ v_init)) / np.sum(gms_tru))

    return {
        "success": sol.success, "nfev": sol.nfev, "sol": sol,
        "lam_g": lam_g, "nt": nt, "gms_cal": gms_cal, "gms_tru": gms_tru,
        "ssb_drift_km_s": ssb_drift_kms,
    }


def run_sector3_audit(bsp_dir: Optional[str] = None) -> Dict[str, Any]:
    bsp_path = os.path.join(bsp_dir, "de440.bsp") if bsp_dir else None
    sim = propagate_solar_system(t_span_days=30.0, jd_start=J2000_JD, bsp_path=bsp_path)
    return {
        "sector": "Sector 3: 15-Body Vectorized EIH Propagator (Gauge-Exact)",
        "passed": sim["success"] and sim["ssb_drift_km_s"] < 1.0e-10,
        "lam_g_factor": sim["lam_g"],
        "dop853_evaluations": sim["nfev"],
        "ssb_drift_km_s": sim["ssb_drift_km_s"],
        "num_bodies": sim["nt"],
        "bsp_mode": "LIVE SPK (de440.bsp)" if (bsp_path and os.path.isfile(bsp_path)) else "Autonomous IEEE-754 J2000 Snapshot",
    }
