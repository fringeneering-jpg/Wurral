# Wurral - Zero-Parameter Continuum Astrodynamics & Telemetry Audit Engine

**Author:** Jaden John Allison (Fringeneering) | fringeneering@gmail.com  
**License:** Dual Licensed - **GNU AGPLv3** (Open Research) / **Fringeneering Commercial Telemetry License** (Proprietary Aerospace & Defense)

Wurral is a deterministic, zero-fitted-parameter Solar System orbital propagator, planetary spin-lock calculator, subatomic constitutive EOS verifier, and deep-space radio/optical link-budget engine. By separating Earth and the Moon dynamically in a 15-body pairwise Einstein-Infeld-Hoffmann (1PN EIH) solver with fixed-period gauge ruler scaling and only four real spacecraft-measured Main-Belt worlds (Ceres, Vesta, Pallas, Hygiea), Wurral eliminates 339 empirical asteroid mass parameters, Kuiper sponge rings, and decadal Chebyshev re-pinning.

## Verified Zero-Parameter Benchmarks

| Audit Sector | Target / Dataset | Wurral Zero-Fit Result | Reference Baseline |
| :--- | :--- | :--- | :--- |
| **Sector 1 (Subatomic EOS)** | Free Neutron Q_beta & Deuteron B_D | 782.301 keV (-31.54 eV) / 2.224612 MeV (+39.39 eV) | CODATA / AME2020 |
| **Sector 2 (Planetary Spin)** | 6 Unlocked Planets (Earth to Neptune) | **1.20% Mean Residual** (Earth +0.45%, Uranus -0.46%) | Observed Sidereal Periods |
| **Sector 3A (VLBA Quasars)** | vlbaobs.txt Saturn (2004-2011, 6.45 yr) | **3.378 mas RA / 1.139 mas Dec** | DE440 (3.377 / 1.139 mas) |
| **Sector 3B (MESSENGER)** | messenger.txt Blind (2004-2015, 45 orbits) | **61.14 m RMS** (+17.04 m median delta vs DE440) | JPL DE440 (2020) |
| **Sector 3C (Cassini DSN)** | cassini_range_avg.txt (2004-2014, 10.0 yr) | **656.90 m RMS** vs DE440 (Zero Kuiper Ring) | JPL DE440 (2020) |
| **Sector 3D (50-Yr Span)** | Continuous Unpinned 1976-2026 (Uranus / EMB) | **65.11 km Uranus** / **1.23 km EMB** / **1.20 km Venus** | DE440 to DE442 (2024) Uranus shift: 1,194.57 km |
| **Link Budget (1/r^2 Growth)** | Dynamic Outbound vs. Inbound Path-Growth | **33.33% (1/3) Outbound : 66.67% (2/3) Inbound** | Speed-Equalized Invariant = 1.000000 |
| **Sector 4 (eTNO Wake)** | 6 JPL Horizons eTNOs + Proto-Sedna (0.40 Myr) | 6/6 Bound, Sedna q detached 30.5 -> ~82 AU | Observed Sedna q = 76.2 AU |

## Quickstart

    # 1. Autonomous Zero-SPICE Mode (Bundled 64-Bit IEEE J2000 Snapshot):
    python wurral_cli.py --all

    # 2. Full Live NASA SPK + Pre-Fit Spacecraft Telemetry Audit:
    python wurral_cli.py --all --bsp-dir "C:/Users/jaden/Downloads/PWC_Solar_Orbits" --prefit-dir "C:/Users/jaden/Downloads/PWC_Solar_Orbits/jpl_prefit"

## Commercial Aerospace & Defense Licensing
Closed-source flight software, commercial mega-constellation station-keeping, and Space Domain Awareness (SDA) deployments require a **Fringeneering Commercial Telemetry License**. See COMMERCIAL_LICENSE.md or contact fringeneering@gmail.com.
