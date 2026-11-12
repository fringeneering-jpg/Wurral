#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2026 Jaden John Allison (Fringeneering). All Rights Reserved.
# Contact: fringeneering@gmail.com | GitHub: github.com/fringeneering-jpg
# Dual License: GNU AGPLv3 / Fringeneering Commercial Telemetry License.

"""wurral_cli.py — Frictionless Command-Line Audit Entry Point for Wurral v1.0.0."""

import argparse
import sys
from wurral import (
    __version__,
    run_sector1_audit,
    run_sector2_audit,
    run_sector3_audit,
    run_telemetry_audit,
    run_sector4_audit,
)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Wurral Zero-Parameter Quantum-to-Macro Astrodynamics CLI")
    ap.add_argument("--all", action="store_true", help="Run all 4 sectors + Link Budget & Telemetry audit")
    ap.add_argument("--sector", type=int, choices=[1, 2, 3, 4], help="Run a single sector (1..4)")
    ap.add_argument("--bsp-dir", default=None, help="Directory containing de440.bsp / de442.bsp")
    ap.add_argument("--prefit-dir", default=None, help="Directory containing vlbaobs.txt, cassini_range_avg.txt, messenger.txt")
    ap.add_argument("--skip-50yr", action="store_true", help="Skip the 50-year 1976-2026 SPK integration")
    ap.add_argument("--out", default=None, help="Write audit ledger to file")
    args = ap.parse_args(argv)

    if not args.all and args.sector is None:
        args.all = True

    lines = [
        "=" * 118,
        f"  WURRAL ASTRODYNAMICS & TELEMETRY AUDIT ENGINE v{__version__} | Dual Licensed: AGPLv3 / Commercial",
        "=" * 118,
    ]
    overall_ok = True

    if args.all or args.sector == 1:
        s1 = run_sector1_audit()
        overall_ok &= s1["passed"]
        lines.extend([
            f"\n[{'PASS' if s1['passed'] else 'FAIL'}] {s1['sector']}",
            f"  * Harmonic V(q) Equilibrium Pressure : {s1['jacobian']['analytic_pressure_mev_fm3']:.6f} MeV/fm^3 (spread: {s1['jacobian']['max_spread_mev_fm3']:.2e})",
            f"  * Master EOS Ceiling P(rho_max)      : {s1['eos']['p_max_pa']:.4e} Pa (c_s^2/c_0^2 = 3/19)",
            f"  * Neutron P*dV Closure               : {s1['neutron']['q_pred_kev']:.3f} keV vs 782.333 keV ({s1['neutron']['residual_ev']:+.2f} eV | Algebraic n-1H: {s1['neutron']['residual_algebraic_ev']:+.1f} eV)",
            f"  * Deuteron State-Transition Closure  : {s1['deuteron']['b_d_pred_mev']:.6f} MeV vs 2.224573 MeV ({s1['deuteron']['delta_ev']:+.2f} eV)",
            f"  * Finite Knot Compressibility Bounds : Volumetric +/-{s1['compressibility']['volumetric_bound_pct']:.1f}% | Mass-Impedance Span = {s1['compressibility']['mass_impedance_span_pct']:.3f}%",
        ])

    if args.all or args.sector == 2:
        s2 = run_sector2_audit()
        overall_ok &= s2["passed"]
        lines.extend([
            f"\n[{'PASS' if s2['passed'] else 'FAIL'}] {s2['sector']}",
            f"  * 6 Unlocked Planets Mean Spin Residual : {s2['spin_lock']['mean_residual_pct']:.2f}% (Zero Free Parameters!)",
        ])
        for nm, pd in s2["spin_lock"]["planets"].items():
            lines.append(f"    - {nm:<8}: Obs = {pd['obs_hr']:7.3f} h | Pred = {pd['pred_hr']:7.3f} h ({pd['signed_error_pct']:+6.2f}%)")
        lines.extend([
            f"  * Exact Lunar Libration Viscosity       : {s2['lunar_viscosity']['viscosity_pa_s']:.3e} Pa*s",
            f"  * Jiba LEO Dipole Decayed at Moon       : {s2['leo_dipole_decay']['decayed_amplitude_mm']:.2f} mm (Local Sun-Earth throat verified)",
        ])

    if args.all or args.sector == 3:
        s3 = run_sector3_audit(args.bsp_dir)
        overall_ok &= s3["passed"]
        lines.extend([
            f"\n[{'PASS' if s3['passed'] else 'FAIL'}] {s3['sector']}",
            f"  * Ephemeris State Source  : {s3['bsp_mode']} ({s3['num_bodies']} bodies, separate Earth & Moon)",
            f"  * Gauge Ruler Factor      : lam_g = {s3['lam_g_factor']:.12f} | SSB Momentum Drift = {s3['ssb_drift_km_s']:.3e} km/s",
        ])
        tel = run_telemetry_audit(args.prefit_dir, args.bsp_dir, run_50yr=not args.skip_50yr)
        overall_ok &= tel['passed']
        lb = tel["link_budget"]["saturn_link"]
        nh = tel["link_budget"]["pluto_nh_link"]
        vk = tel["link_budget"]["viking"]
        lines.extend([
            f"\n[{'PASS' if tel['passed'] else 'FAIL'}] Dynamic 1/r^2 Link Budget & Telemetry Pre-Fit Audit ({tel['mode']})",
            f"  * Small-dr Gradient Split : Outbound = {lb['small_dr_outbound_pct']:.2f}% (1/3) | Inbound = {lb['small_dr_inbound_pct']:.2f}% (2/3)",
            f"  * Speed-Equalized Ratio   : (dL_probe * v_p)/(dL_light * c_0) = {lb['speed_equalized_invariant']:.6f} (Amp = {nh['c_over_v_amplification']:.1f}x at Pluto)",
            f"  * Far-Field Crossovers    : r_entropy = {lb['r_entropy_crossover_au']:.1f} AU | r_join = {lb['r_join_crossover_au']:.1f} AU | Viking Cos Ratio = {vk['cosine_ratio']:.4f}",
        ])
        if "vlba_ra_rms_mas" in tel:
            lines.extend([
                f"  * Live 6.45-Yr VLBA Saturn RMS   : RA*cosD = {tel['vlba_ra_rms_mas']:.3f} mas | Dec = {tel['vlba_dec_rms_mas']:.3f} mas",
                f"  * Live 10.0-Yr Cassini Range RMS : {tel['cassini_rms_m']:.2f} m vs DE440 (Zero Kuiper Ring)",
                f"  * Live 10.8-Yr MESSENGER Blind   : {tel['messenger_rms_m']:.2f} m RMS (Median: {tel['messenger_med_m']:+.2f} m across 45 Mercury orbits)",
            ])
        if "benchmark_50yr" in tel:
            lines.append("  * Live 50.0-Yr Continuous Unpinned Drift (1976 Viking -> 2026 vs DE440/DE442):")
            for b_nm, d_inf in tel["benchmark_50yr"].items():
                sh = f" (DE440->DE442 shift: {d_inf['de440_vs_de442_km']:.3f} km)" if d_inf["de440_vs_de442_km"] is not None else ""
                lines.append(f"    - {b_nm:<8}: {d_inf['wurral_50yr_km']:8.3f} km{sh}")

    if args.all or args.sector == 4:
        s4 = run_sector4_audit()
        overall_ok &= s4["passed"]
        lines.extend([
            f"\n[{'PASS' if s4['passed'] else 'FAIL'}] {s4['sector']}",
            f"  * 0.40 Myr Symplectic Run : {s4['orbits_bound']} eTNOs Bound | Proto-Sedna q lifted {s4['proto_sedna_q_init_au']:.1f} -> {s4['proto_sedna_q_final_au']:.2f} AU",
        ])

    lines.extend(["", "=" * 118, f"  OVERALL VERDICT: {'PASS (EXIT 0)' if overall_ok else 'FAIL (EXIT 1)'}", "=" * 118])
    report_str = "\n".join(lines)
    print(report_str)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(report_str + "\n")
    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
