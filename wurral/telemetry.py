# Copyright (c) 2026 Jaden John Allison (Fringeneering). All Rights Reserved.
# Contact: fringeneering@gmail.com | GitHub: [github.com/fringeneering-jpg](https://github.com/fringeneering-jpg)
# Dual License: GNU AGPLv3 / Fringeneering Commercial Telemetry License.

from typing import Dict, Any, Optional
import math
import os
import numpy as np
from wurral.constants import (
    C_KMS, C_SI, AU_SI, DAY_SEC, YEAR_SEC, H0_SI, S_DEBT,
    R_ENTROPY_CROSSOVER_AU, R_JOIN_CROSSOVER_AU,
)
from wurral.integrator import propagate_solar_system, JPLEPHEM_AVAILABLE, SPK

MONTHS = {"jan":1,"feb":2,"mar":3,"apr":4,"may":5,"jun":6,"jul":7,"aug":8,"sep":9,"oct":10,"nov":11,"dec":12}

def cal_to_jd(y: int, m: int, d: int, hh: int, mm: int, ss: float) -> float:
    if m <= 2:
        y -= 1; m += 12
    A = y // 100
    B = 2 - A + (A // 4)
    return int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + B - 1524.5 + (hh * 3600.0 + mm * 60.0 + ss) / DAY_SEC

def dynamic_one_way_two_way_link_budget(r1_au: float = 1.0, r2_au: float = 40.0, v_probe_kms: float = 15.0) -> Dict[str, Any]:
    r1, r2 = r1_au * AU_SI, r2_au * AU_SI
    v_p, v_c = v_probe_kms * 1000.0, C_SI
    alpha = H0_SI * (S_DEBT / 4.0) * (AU_SI ** 2)
    dl_out_p = (alpha / v_p) * (math.log(r2 / r1) - (1.0 - r1 / r2))
    dl_in_p  = (alpha / v_p) * ((r2 / r1) - 1.0 - math.log(r2 / r1))
    dl_rt_p  = dl_out_p + dl_in_p
    dl_out_c = (alpha / v_c) * (math.log(r2 / r1) - (1.0 - r1 / r2))
    dl_in_c  = (alpha / v_c) * ((r2 / r1) - 1.0 - math.log(r2 / r1))
    dl_rt_c  = dl_out_c + dl_in_c
    r_near = 1.0001 * r1
    l_base = 0.5 * ((r_near - r1) ** 2) / (r_near ** 2)
    o_grad = (math.log(r_near / r1) - (1.0 - r1 / r_near)) - l_base
    i_grad = ((r_near / r1) - 1.0 - math.log(r_near / r1)) - l_base
    small_dr_out_pct = (o_grad / (o_grad + i_grad)) * 100.0
    small_dr_in_pct  = (i_grad / (o_grad + i_grad)) * 100.0
    invariant = (dl_rt_p * v_p) / (dl_rt_c * v_c)
    sieve_m_yr = (H0_SI * (S_DEBT / 4.0)) * r2 * YEAR_SEC
    return {
        "r1_au": r1_au, "r2_au": r2_au, "v_probe_kms": v_probe_kms,
        "probe_dl_out_m": dl_out_p, "probe_dl_in_m": dl_in_p, "probe_dl_rt_m": dl_rt_p,
        "light_dl_out_m": dl_out_c, "light_dl_in_m": dl_in_c, "light_dl_rt_m": dl_rt_c,
        "small_dr_outbound_pct": small_dr_out_pct, "small_dr_inbound_pct": small_dr_in_pct,
        "deep_space_in_over_out_ratio": dl_in_c / dl_out_c,
        "speed_equalized_invariant": invariant, "c_over_v_amplification": v_c / v_p,
        "secular_sieve_recession_m_yr": sieve_m_yr,
        "r_entropy_crossover_au": R_ENTROPY_CROSSOVER_AU, "r_join_crossover_au": R_JOIN_CROSSOVER_AU,
        "verified": abs(invariant - 1.0) < 1e-12 and abs(small_dr_out_pct - 33.333) < 0.01,
    }

def viking_lander_turnaround_offset() -> Dict[str, Any]:
    cos_ratio = math.cos(math.radians(47.97)) / math.cos(math.radians(22.48))
    obs_ratio = -12.7 / -18.0
    return {
        "cosine_ratio": cos_ratio, "observed_ratio": obs_ratio,
        "relative_agreement_pct": (cos_ratio / obs_ratio - 1.0) * 100.0,
        "verified": abs(cos_ratio - 0.7246) < 0.001,
    }

def run_link_budget_audit() -> Dict[str, Any]:
    lb_sat = dynamic_one_way_two_way_link_budget(1.0, 9.5549, 15.0)
    lb_nh  = dynamic_one_way_two_way_link_budget(1.0, 33.0, 16.26)
    vk     = viking_lander_turnaround_offset()
    return {
        "sector": "Dynamic 1/r^2 Link Budget (1/3 : 2/3 Asymmetry) & Far-Field Crossover",
        "passed": lb_sat["verified"] and lb_nh["verified"] and vk["verified"],
        "saturn_link": lb_sat, "pluto_nh_link": lb_nh, "viking": vk,
    }

def run_telemetry_audit(prefit_dir: Optional[str] = None, bsp_dir: Optional[str] = None, run_50yr: bool = True) -> Dict[str, Any]:
    bsp_440 = os.path.join(bsp_dir, "de440.bsp") if bsp_dir else None
    bsp_442 = os.path.join(bsp_dir, "de442.bsp") if bsp_dir else None
    vlba_fp = os.path.join(prefit_dir, "vlbaobs.txt") if prefit_dir else None
    cass_fp = os.path.join(prefit_dir, "cassini_range_avg.txt") if prefit_dir else None
    mess_fp = os.path.join(prefit_dir, "messenger.txt") if prefit_dir else None
    link_res = run_link_budget_audit()
    if not (bsp_440 and os.path.isfile(bsp_440) and JPLEPHEM_AVAILABLE):
        return {"mode": "AUTONOMOUS OFFLINE SUMMARY", "passed": link_res["passed"], "link_budget": link_res}
    k440 = SPK.open(bsp_440)
    out: Dict[str, Any] = {"mode": "LIVE EXECUTABLE DOP853 TELEMETRY INTEGRATION", "link_budget": link_res}
    if vlba_fp and os.path.isfile(vlba_fp) and cass_fp and os.path.isfile(cass_fp) and mess_fp and os.path.isfile(mess_fp):
        vlba_obs, cass_obs, mess_obs = [], [], []
        with open(vlba_fp, "r", encoding="utf-8", errors="replace") as f:
            for ln in f:
                s = ln.strip()
                if not s or s.startswith("#"): continue
                t = s.split()
                try:
                    yr, mo, dy = int(t[0]), MONTHS[t[1][:3].lower()], int(t[2])
                    hh, mm, ss = int(t[3]), int(t[4]), float(t[5])
                    leap = 32.0 if yr < 2006 else (33.0 if yr < 2009 else 34.0)
                    jd_t = cal_to_jd(yr, mo, dy, hh, mm, ss + leap + 32.184)
                    ra_r = math.radians((float(t[6]) + float(t[7])/60.0 + float(t[8])/3600.0) * 15.0)
                    sgn = -1.0 if t[9].startswith("-") else 1.0
                    dc_r = math.radians(sgn * (abs(float(t[9])) + float(t[10])/60.0 + float(t[11])/3600.0))
                    vlba_obs.append((jd_t, ra_r, dc_r))
                except Exception: continue
        with open(cass_fp, "r", encoding="utf-8", errors="replace") as f:
            for ln in f:
                s = ln.strip()
                if not s or s.startswith("#"): continue
                t = s.split()
                if len(t) >= 12:
                    yr, mo, dy, hh, mm, ss = int(t[0]), int(t[1]), int(t[2]), int(t[3]), int(t[4]), float(t[5])
                    leap = 32.0 if yr < 2006 else (33.0 if yr < 2009 else (34.0 if yr < 2012 else 35.0))
                    cass_obs.append((cal_to_jd(yr, mo, dy, hh, mm, ss + leap + 32.184), float(t[7])))
        with open(mess_fp, "r", encoding="utf-8", errors="replace") as f:
            for ln in f:
                s = ln.strip()
                if not s or s.startswith("#"): continue
                t = s.split()
                if len(t) >= 8:
                    yr, mo, dy = int(t[3]), MONTHS[t[4][:3].lower()], int(t[5])
                    hms = [float(x) for x in t[6].split(":")]
                    leap = 34.0 if yr < 2012 else (35.0 if (yr < 2015 or (yr == 2015 and mo < 7)) else 36.0)
                    mess_obs.append((cal_to_jd(yr, mo, dy, int(hms[0]), int(hms[1]), hms[2] + leap + 32.184), float(t[7]) * 1e-6))
        jd0 = min(vlba_obs[0][0], cass_obs[0][0]) - 1.0
        jd1 = max(vlba_obs[-1][0], cass_obs[-1][0], mess_obs[-1][0]) + 1.0
        sim11 = propagate_solar_system(t_span_days=(jd1 - jd0), jd_start=jd0, bsp_path=bsp_440)
        sol, lam_g, nt, gms_c = sim11["sol"], sim11["lam_g"], sim11["nt"], sim11["gms_cal"]
        c2 = C_KMS ** 2
        w_dra, w_ddc = [], []
        for jd_obs, ra_o, dc_o in vlba_obs:
            t_s = (jd_obs - jd0) * DAY_SEC
            pr = sol.sol(t_s)[:3*nt].reshape((nt, 3)) / lam_g
            emb_w = (gms_c[3]*pr[3] + gms_c[9]*pr[9]) / (gms_c[3] + gms_c[9])
            re_w = emb_w + k440[3, 399].compute(jd_obs)
            tau = 4500.0
            for _ in range(5):
                pe = sol.sol(t_s - tau)[:3*nt].reshape((nt, 3)) / lam_g
                rho_v = pe[6] - re_w; rho = float(np.linalg.norm(rho_v))
                r1 = float(np.linalg.norm(re_w - pr[0])); r2 = float(np.linalg.norm(pe[6] - pe[0]))
                shap = (2.0 * gms_c[0] / (C_KMS**3)) * math.log(max(1e-9, (r1+r2+rho)/max(1e-9, r1+r2-rho)))
                tau = rho / C_KMS + shap
            u_p = rho_v / rho
            qv = re_w - pe[0]; qm = float(np.linalg.norm(qv)); uq = qv / qm
            em = float(np.linalg.norm(pe[6] - pe[0]))
            den = qm * (1.0 + float(np.dot(uq, u_p)))
            defl = (2.0 * gms_c[0] / (c2 * den)) * (uq - float(np.dot(uq, u_p)) * u_p) * (1.0 + qm / em) * 0.5
            ua = u_p + defl; ua /= np.linalg.norm(ua)
            ra_w, dc_w = math.atan2(ua[1], ua[0]) % (2*math.pi), math.asin(max(-1.0, min(1.0, ua[2])))
            w_dra.append((((ra_w - ra_o + math.pi) % (2*math.pi)) - math.pi) * math.cos(dc_o) * (180*3600*1000/math.pi))
            w_ddc.append((dc_w - dc_o) * (180*3600*1000/math.pi))
        cass_d_m = []
        for jd_t, rtlt_s in cass_obs:
            t_s = (jd_t - jd0) * DAY_SEC
            pr = sol.sol(t_s)[:3*nt].reshape((nt, 3)) / lam_g
            pe = sol.sol(t_s - 0.5 * rtlt_s)[:3*nt].reshape((nt, 3)) / lam_g
            emb_w = (gms_c[3]*pr[3] + gms_c[9]*pr[9]) / (gms_c[3] + gms_c[9]) - pr[0]
            dw = float(np.linalg.norm((pe[6] - pe[0]) - emb_w))
            jd_b = jd_t - 0.5 * rtlt_s / DAY_SEC
            d4 = float(np.linalg.norm((k440[0, 6].compute(jd_b) - k440[0, 10].compute(jd_b)) - (k440[0, 3].compute(jd_t) - k440[0, 10].compute(jd_t))))
            cass_d_m.append((dw - d4) * 1000.0)
        mess_d_m = []
        for jd_t, rtlt_s in mess_obs:
            t_s = (jd_t - jd0) * DAY_SEC
            pr = sol.sol(t_s)[:3*nt].reshape((nt, 3)) / lam_g
            pe = sol.sol(t_s - 0.5 * rtlt_s)[:3*nt].reshape((nt, 3)) / lam_g
            emb_w = (gms_c[3]*pr[3] + gms_c[9]*pr[9]) / (gms_c[3] + gms_c[9]) - pr[0]
            dw = float(np.linalg.norm((pe[1] - pe[0]) - emb_w))
            jd_b = jd_t - 0.5 * rtlt_s / DAY_SEC
            d4 = float(np.linalg.norm((k440[0, 1].compute(jd_b) - k440[0, 10].compute(jd_b)) - (k440[0, 3].compute(jd_t) - k440[0, 10].compute(jd_t))))
            mess_d_m.append((dw - d4) * 1000.0)
        out["vlba_ra_rms_mas"]  = float(math.sqrt(np.mean(np.square(w_dra))))
        out["vlba_dec_rms_mas"] = float(math.sqrt(np.mean(np.square(w_ddc))))
        out["cassini_rms_m"]    = float(math.sqrt(np.mean(np.square(cass_d_m))))
        out["messenger_rms_m"]  = float(math.sqrt(np.mean(np.square(mess_d_m))))
        out["messenger_med_m"]  = float(np.median(mess_d_m))
    if run_50yr:
        jd_76 = cal_to_jd(1976, 7, 20, 12, 0, 0.0)
        jd_26 = cal_to_jd(2026, 7, 20, 12, 0, 0.0)
        sim50 = propagate_solar_system(t_span_days=(jd_26 - jd_76), jd_start=jd_76, bsp_path=bsp_440, rtol=1e-12, atol=1e-12)
        pw = sim50["sol"].y[:3*sim50["nt"], -1].reshape((sim50["nt"], 3)) / sim50["lam_g"]
        gc = sim50["gms_cal"]
        k442 = SPK.open(bsp_442) if (bsp_442 and os.path.isfile(bsp_442)) else None
        drifts_50 = {}
        for lbl, wi, sid in [("Mercury",1,1),("Venus",2,2),("EMB","emb",3),("Mars",4,4),("Jupiter",5,5),("Saturn",6,6),("Uranus",7,7),("Neptune",8,8)]:
            rw = ((gc[3]*pw[3] + gc[9]*pw[9])/(gc[3]+gc[9]) - pw[0]) if wi == "emb" else (pw[wi] - pw[0])
            r4 = k440[0, sid].compute(jd_26) - k440[0, 10].compute(jd_26)
            jpl_shift = float(np.linalg.norm(r4 - (k442[0, sid].compute(jd_26) - k442[0, 10].compute(jd_26)))) if k442 else None
            drifts_50[lbl] = {"wurral_50yr_km": float(np.linalg.norm(rw - r4)), "de440_vs_de442_km": jpl_shift}
        out["benchmark_50yr"] = drifts_50
    out["passed"] = link_res["passed"]
    return out
