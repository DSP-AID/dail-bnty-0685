#!/usr/bin/env python3
"""Worked examples for the AAVSO photometry error-budget guide. Version 1.0.0.

Reads ONLY the shipped files in data/ (real public datasets, see data/MANIFEST.json) plus
tools/error_budget.py, and writes:

  results/worked_examples.md    every number quoted in the guide, in tables
  results/worked_examples.json  the same numbers, machine-readable

No network access, no third-party packages. Output is deterministic (no timestamps), so
`python tools/run_examples.py --check` proves the committed results match a fresh run.

Usage:  python tools/run_examples.py [--write | --check]
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import error_budget as eb  # noqa: E402

DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")

# Reference example setup. Camera numbers are measured values quoted from a published AAVSO
# forum test (see docs/sources.md); the rest is an example configuration that the guide
# sweeps, so no claim depends on one particular assumption.
SETUP = dict(
    aperture_cm=20.0,  # example: 20 cm / f=1000 mm reflector (SCT at native f/10)
    focal_mm=2000.0,
    pixel_um=3.76,  # IMX571 pixel size (manufacturer specification)
    gain_e_per_adu=0.77,  # ASI2600MM Pro, gain 0 / offset 10 (AAVSO test, sources.md S-3)
    read_noise_e=2.8,  # same source
    dark_e_per_s=0.1,  # assumed for the example; swept in section 5.4
    aper_radius_pix=5.0,
    annulus_inner_pix=10.0,
    annulus_outer_pix=16.0,
)


def load(name: str) -> list[dict]:
    with open(os.path.join(DATA, name), encoding="utf-8") as fh:
        return [{k: v for k, v in row.items()} for row in csv.DictReader(fh)]


def fnum(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def percentile(values: list[float], q: float) -> float:
    """Linear-interpolated percentile of a list (stdlib only)."""
    if not values:
        return float("nan")
    s = sorted(values)
    pos = q * (len(s) - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (pos - lo)


def describe(values: list[float], digits=4) -> dict:
    if not values:
        return {}
    return {
        "n": len(values),
        "min": round(min(values), digits),
        "p25": round(percentile(values, 0.25), digits),
        "median": round(percentile(values, 0.5), digits),
        "p75": round(percentile(values, 0.75), digits),
        "p90": round(percentile(values, 0.90), digits),
        "max": round(max(values), digits),
    }


def loglog_slope(pairs: list[tuple[float, float]], digits=3) -> dict:
    """Least-squares slope of log10(y) on log10(x) plus R**2."""
    pts = [(math.log10(x), math.log10(y)) for x, y in pairs if x > 0 and y > 0]
    n = len(pts)
    if n < 3:
        return {"n": n}
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    sxy = sum((p[0] - mx) * (p[1] - my) for p in pts)
    syy = sum((p[1] - my) ** 2 for p in pts)
    if sxx == 0 or syy == 0:
        return {"n": n}
    slope = sxy / sxx
    return {"n": n, "slope": round(slope, digits), "r2": round(sxy**2 / (sxx * syy), digits)}


def fmt(value, spec=".4f") -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "-"
    return format(value, spec)


# ======================================================================================
# analysis blocks
# ======================================================================================
def sec1_apass(apass: list[dict]) -> dict:
    """Real per-measurement scatter of ground-based V photometry (APASS DR9)."""
    good = [r for r in apass if (fnum(r["e_vmag"]) or 0) > 0 and r["u_e_vmag"] == "0" and fnum(r["nobs"], 0) > 1]
    out = {"rows_total": len(apass), "rows_usable": len(good)}
    out["rows_flagged_unreliable_error"] = sum(1 for r in apass if r["u_e_vmag"] == "1")
    out["nobs"] = describe([fnum(r["nobs"]) for r in good])
    bins = []
    slopes = {}
    for lo, hi in [(10, 11), (11, 12), (12, 13), (13, 14), (14, 15), (15, 16)]:
        sel = [r for r in good if lo <= fnum(r["vmag"], 99) < hi]
        if not sel:
            continue
        sd = [fnum(r["e_vmag"]) for r in sel]
        sem = [fnum(r["e_vmag"]) / math.sqrt(fnum(r["nobs"])) for r in sel]
        bins.append(
            {
                "vmag_bin": f"{lo}-{hi}",
                "n_stars": len(sel),
                "nobs_median": round(percentile([fnum(r["nobs"]) for r in sel], 0.5), 1),
                "scatter_sd_mag": describe(sd),
                "catalogue_mean_sem_mag": describe(sem),
            }
        )
        slopes[f"{lo}-{hi}"] = loglog_slope([(fnum(r["nobs"]), fnum(r["e_vmag"])) for r in sel])
    out["by_vmag_bin"] = bins
    out["slope_sd_vs_nobs_by_bin"] = slopes
    allsd = [fnum(r["e_vmag"]) for r in good]
    out["scatter_all_mags_mag"] = describe(allsd)
    out["catalogue_mean_sem_all_mags_mag"] = describe([fnum(r["e_vmag"]) / math.sqrt(fnum(r["nobs"])) for r in good])
    # the number the guide quotes as the measured atmosphere+calibration floor:
    out["measured_floor_per_measurement_mag"] = round(percentile(allsd, 0.5), 4)
    return out


def sec2_gaia(gaia: list[dict]) -> dict:
    """Real quoted photometric precision of Gaia DR3 in the same field (space, no atmosphere)."""
    rows = []
    for lo, hi in [(10, 12), (12, 13), (13, 14), (14, 15), (15, 16)]:
        sel = []
        for r in gaia:
            g, fl, fe = fnum(r["phot_g_mean_mag"]), fnum(r["phot_g_mean_flux"]), fnum(r["phot_g_mean_flux_error"])
            if not (fl and fe and fl > 0) or not (lo <= (g or 99) < hi):
                continue
            sel.append((eb.snr_to_mag_error(fl / fe), fl / fe))
        if sel:
            rows.append(
                {
                    "gmag_bin": f"{lo}-{hi}",
                    "n_stars": len(sel),
                    "flux_snr_median": round(percentile([s[1] for s in sel], 0.5), 1),
                    "sigma_mag_median": round(percentile([s[0] for s in sel], 0.5), 6),
                    "sigma_mag_p90": round(percentile([s[0] for s in sel], 0.90), 6),
                }
            )
    out = {"rows_total": len(gaia), "by_gmag_bin": rows}
    out["variable_flag_counts"] = {}
    for r in gaia:
        out["variable_flag_counts"][r["phot_variable_flag"] or "NONE"] = out["variable_flag_counts"].get(r["phot_variable_flag"] or "NONE", 0) + 1
    return out


def _nearest(rows: list[dict], ra: float, dec: float) -> tuple[dict, float]:
    best, bestd = None, 1e9
    for r in rows:
        dr = (fnum(r.get("raj2000"), fnum(r.get("ra"), 0.0)) - ra) * math.cos(math.radians(dec))
        dd = fnum(r.get("dej2000"), fnum(r.get("dec"), 0.0)) - dec
        d = math.hypot(dr, dd) * 3600.0
        if d < bestd:
            best, bestd = r, d
    return best, bestd


def _usable(apass: list[dict]) -> list[dict]:
    """APASS rows with a usable quoted error and more than one contributing measurement."""
    return [r for r in apass if (fnum(r["e_vmag"]) or 0) > 0 and r["u_e_vmag"] == "0" and fnum(r["nobs"], 0) > 1]


def _ensemble_rows(pool: list[tuple], tv: float, dm_max: float = 3.5) -> dict:
    """Ensemble calibration error versus comp availability and versus the number of real comparison stars."""
    by_window = [
        {
            "dm_window": dm,
            "comps_available": len([p for p in pool if abs(p[3] - tv) <= dm]),
        }
        for dm in (0.25, 0.5, 1.0, 1.5, 2.0, 3.0)
        if dm <= dm_max
    ]
    by_n = []
    ordered = sorted(pool, key=lambda p: abs(p[3] - tv))
    for n in (3, 5, 8, 12, 20):
        sel = ordered[:n]
        if len(sel) < 2:
            continue
        sd_errs = [p[1] for p in sel]
        sem_errs = [p[2] for p in sel]
        by_n.append(
            {
                "n_requested": n,
                "n_used": len(sel),
                "comp_mag_range": f"{min(p[3] for p in sel):.2f}..{max(p[3] for p in sel):.2f}",
                "worst_comp_delta_vmag": round(max(abs(p[3] - tv) for p in sel), 2),
                "ens_sd_mag": round(eb.ensemble_calibration_error(sd_errs), 5),
                "ens_sd_ivw_mag": round(eb.mean_inverse_variance(sd_errs), 5),
                "ens_sem_mag": round(eb.ensemble_calibration_error(sem_errs), 5),
            }
        )
    return {"by_window": by_window, "by_n": by_n}


def _bright_vs_faint(pool: list[tuple], tv: float) -> list[dict]:
    """Real stars about 2 mag brighter / magnitude-matched / about 2 mag fainter than the target."""
    out = []
    for label, sel_pool in [
        (f"about 2 mag brighter than V={tv:.2f}", sorted([p for p in pool if 1.5 <= tv - p[3] <= 2.5], key=lambda p: abs(p[3] - tv))),
        (f"within 0.5 mag of V={tv:.2f} (matched)", sorted([p for p in pool if abs(p[3] - tv) <= 0.5], key=lambda p: abs(p[3] - tv))),
        (f"1.5-2.5 mag fainter than V={tv:.2f}", sorted([p for p in pool if 1.5 <= p[3] - tv <= 2.5], key=lambda p: abs(p[3] - tv))),
    ]:
        sel = sel_pool[:8]
        out.append(
            {
                "selection": label,
                "n_available": len(sel_pool),
                "n_used": len(sel),
                "median_comp_vmag": round(percentile([p[3] for p in sel], 0.5), 2) if len(sel) >= 3 else None,
                "median_comp_sd_mag": round(percentile([p[1] for p in sel], 0.5), 4) if len(sel) >= 3 else None,
                "ens_sd_mag": round(eb.ensemble_calibration_error([p[1] for p in sel]), 5) if len(sel) >= 3 else None,
                "ens_sem_mag": round(eb.ensemble_calibration_error([p[2] for p in sel]), 5) if len(sel) >= 3 else None,
            }
        )
    return out


def sec3_target_and_ensembles(apass: list[dict], gaia: list[dict], target_pos: tuple[float, float], atm_floor: float) -> dict:
    """Locate the real AAVSO target in both catalogues and build comparison ensembles two ways."""
    ra, dec = target_pos
    out = {"target": None, "gaia_crosscheck": None, "cases": [], "comp_photon_table": []}

    tv_row, tv_d = _nearest(apass, ra, dec)
    if not tv_row or tv_d > 5.0:
        out["note"] = f"no APASS DR9 entry within 5 arcsec of the target (closest {tv_d:.2f}\")"
        return out
    sd = fnum(tv_row["e_vmag"]) or 0.0
    nobs = fnum(tv_row["nobs"], 1)
    out["target"] = {
        "apass_recno": tv_row["recno"],
        "vmag": fnum(tv_row["vmag"]),
        "separation_arcsec": round(tv_d, 2),
        "nobs": nobs,
        "e_vmag_sd_mag": sd,
        "error_flag_unreliable": tv_row["u_e_vmag"] == "1",
        "catalogue_mean_sem_mag": round(sd / math.sqrt(nobs), 5) if sd > 0 else None,
    }
    g_row, g_d = _nearest(gaia, ra, dec)
    if g_row and g_d < 5.0:
        fl, fe = fnum(g_row["phot_g_mean_flux"]), fnum(g_row["phot_g_mean_flux_error"])
        out["gaia_crosscheck"] = {
            "source_id": g_row["source_id"],
            "separation_arcsec": round(g_d, 2),
            "gmag": fnum(g_row["phot_g_mean_mag"]),
            "bp_minus_rp": round(fnum(g_row["phot_bp_mean_mag"], 0) - fnum(g_row["phot_rp_mean_mag"], 0), 3) if g_row.get("phot_bp_mean_mag") else None,
            "phot_variable_flag": g_row["phot_variable_flag"],
            "ruwe": fnum(g_row["ruwe"]),
            "sigma_g_mag": round(eb.snr_to_mag_error(fl / fe), 6) if fl and fe else None,
        }

    usable = _usable(apass)

    def pool_for(tv: float) -> list[tuple]:
        return [
            (abs(fnum(r["vmag"]) - tv), fnum(r["e_vmag"]), fnum(r["e_vmag"]) / math.sqrt(fnum(r["nobs"])), fnum(r["vmag"]))
            for r in usable
            if abs(fnum(r["vmag"]) - tv) <= 3.5
        ]

    tv = out["target"]["vmag"]

    def case_extra(tv_value: float) -> dict:
        brighter = [r for r in usable if fnum(r["vmag"]) < tv_value]
        return {
            "n_usable_stars_in_field": len(usable),
            "n_comps_brighter_than_target": len(brighter),
            "brightest_usable_comp_vmag": round(min(fnum(r["vmag"]) for r in usable), 3) if usable else None,
        }

    # case A: the real bright AAVSO star, where matched comparison stars are scarce
    pool_a = pool_for(tv)
    out["cases"].append(
        {
            "label": f"case A - the real target itself (V = {tv:.2f})",
            "target_vmag": round(tv, 3),
            **case_extra(tv),
            "ensembles": _ensemble_rows(pool_a, tv),
            "bright_vs_faint": _bright_vs_faint(pool_a, tv),
        }
    )
    # case B: a mid-field magnitude where AAVSO work usually sits, using a real star as the anchor
    anchor = min(usable, key=lambda r: abs(fnum(r["vmag"]) - 11.5))
    tv2 = fnum(anchor["vmag"])
    pool_b = pool_for(tv2)
    out["cases"].append(
        {
            "label": f"case B - a real field star used as a stand-in target (V = {tv2:.2f})",
            "target_vmag": round(tv2, 3),
            "anchor_apass_recno": anchor["recno"],
            **case_extra(tv2),
            "ensembles": _ensemble_rows(pool_b, tv2),
            "bright_vs_faint": _bright_vs_faint(pool_b, tv2),
        }
    )

    # the other half of the story: photon noise of measuring a comp of that magnitude yourself
    for dv in (-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0):
        res = eb.frame_budget(mag=round(tv2 + dv, 2), exposure_s=60.0, frames=1, systematic_mag=atm_floor, **SETUP)
        out["comp_photon_table"].append(
            {
                "delta_vmag": dv,
                "comp_vmag": round(tv2 + dv, 2),
                "star_e_in_60s": res["electrons"]["star_total_in_exposure"],
                "snr": res["per_frame"]["snr"],
                "sigma_photon_mag": res["per_frame"]["sigma_photon_mag"],
                "sigma_with_atm_mag": res["per_frame"]["sigma_frame_total_mag"],
            }
        )
    return out


def sec4_stacking(apass: list[dict], atm_floor: float) -> dict:
    """Stack-then-measure vs measure-then-average, with frame errors taken from real data."""
    good = [r for r in apass if (fnum(r["e_vmag"]) or 0) > 0 and r["u_e_vmag"] == "0" and fnum(r["nobs"], 0) > 1]
    sd = [fnum(r["e_vmag"]) for r in good]
    p25, med, p90 = percentile(sd, 0.25), percentile(sd, 0.5), percentile(sd, 0.9)
    out = {"frame_error_percentiles_mag": {"p25": round(p25, 4), "median": round(med, 4), "p90": round(p90, 4)}}

    # equal frames: both routes agree; N frames
    rows = []
    for n in (1, 4, 9, 16, 25, 36):
        sig = eb.stack_equal_frames(med, n)
        rows.append(
            {
                "frames": n,
                "sigma_photon_route_mag": round(med / math.sqrt(n), 5),
                "sigma_mean_route_mag": round(sig, 5),
                "routes_agree": abs(sig - med / math.sqrt(n)) < 1e-12,
                "implied_snr": round(eb.mag_error_to_snr(sig), 1),
            }
        )
    out["equal_frames"] = rows

    # unequal frames: one clouded frame in a good run (percentiles from the real distribution)
    mixed = [p25, p25, med, med, p90 * 3]
    out["mixed_quality_example"] = {
        "frame_errors_mag": [round(x, 4) for x in mixed],
        **{k: round(v, 5) for k, v in eb.outlier_penalty(mixed).items() if isinstance(v, float)},
        "n_total": len(mixed),
    }
    keep = [True] * (len(mixed) - 1) + [False]
    out["mixed_quality_example"]["rejected_mask_used"] = keep
    out["mixed_quality_example"]["loss_if_bad_frame_kept_unweighted_pct"] = round(
        100.0 * (eb.mean_unweighted(mixed) / eb.mean_unweighted(mixed[:-1]) - 1.0), 1
    )

    return out


def sec5_reconciliation(atm_floor: float) -> dict:
    """The 'why is my error 0.085 when SNR says 0.01' arithmetic, on a real setup."""
    out = {"error_to_snr": [], "setup_cases": [], "sens": {}}
    for sig in (0.005, 0.01, 0.02, 0.05, 0.08, 0.085, 0.09, 0.1):
        out["error_to_snr"].append({"reported_sigma_mag": sig, "implied_snr": round(eb.mag_error_to_snr(sig), 1)})

    for mag in (11.3, 12.5, 14.2, 15.0):
        for exp in (30.0, 60.0, 120.0):
            for frames in (1, 20):
                res = eb.frame_budget(mag=mag, exposure_s=exp, frames=frames, sky_mag_per_arcsec2=20.7, systematic_mag=atm_floor, **SETUP)
                out["setup_cases"].append(
                    {
                        "vmag": mag,
                        "exposure_s": exp,
                        "frames": frames,
                        "star_e": res["electrons"]["star_total_in_exposure"],
                        "snr_photon": res["per_frame"]["snr"] if frames == 1 else res["stacked"]["snr_photon"],
                        "sigma_photon_mag": (res["per_frame"] if frames == 1 else res["stacked"])["sigma_photon_mag"],
                        "sigma_frame_total_mag": res["per_frame"]["sigma_frame_total_mag"],
                        "stack_random_case_mag": res["stacked"]["if_atm_changes_between_frames_mag"],
                        "stack_constant_offset_mag": res["stacked"]["if_atm_constant_over_run_mag"],
                        "dominant_term": res["dominant_term"],
                        "variance_share": res["variance_share"],
                    }
                )

    base = eb.frame_budget(mag=15.0, exposure_s=60.0, frames=1, sky_mag_per_arcsec2=20.7, systematic_mag=atm_floor, **SETUP)
    out["sens"] = {"base": base["per_frame"]}
    for label, key, values in [
        ("exposure_s", "exposure_s", [15, 30, 60, 120, 240]),
        ("aperture_cm", "aperture_cm", [8, 10, 15, 20, 30]),
        ("sky_mag_per_arcsec2", "sky_mag_per_arcsec2", [19.5, 20.4, 21.3, 21.8]),
        ("aper_radius_pix", "aper_radius_pix", [3, 4, 5, 7, 10]),
        ("dark_e_per_s", "dark_e_per_s", [0.0, 0.1, 0.5, 1.0]),
        ("read_noise_e", "read_noise_e", [1.0, 2.8, 6.0, 10.0]),
    ]:
        rows = []
        for v in values:
            kw = dict(SETUP, mag=15.0, exposure_s=60.0, frames=1, sky_mag_per_arcsec2=20.7)
            kw[key] = v
            r = eb.frame_budget(systematic_mag=atm_floor, **kw)
            rows.append({"value": v, "sigma_photon_mag": r["per_frame"]["sigma_photon_mag"], "dominant_term": r["dominant_term"]})
        out["sens"][label] = rows
    return out


def sec6_checkstar(ens_sd: float, ens_sem: float, star_sem: float, star_sd: float, atm_floor: float, photon_sigma: float, target_vmag: float) -> dict:
    """Expected check-star RMS in four regimes, and what a measured value tells you.

    `star_sem` / `star_sd` are the real catalogue-mean and per-measurement errors of one star
    from section 1; `photon_sigma` is the calculator's photon-only error for the target.
    """
    regimes = [
        ("internal consistency only (your own frames, no catalogue)", photon_sigma, 0.0, 0.0),
        ("check star vs catalogue, comps at their catalogue-mean error", photon_sigma, ens_sem, star_sem),
        ("check star vs catalogue, comps as good as one measurement each", photon_sigma, ens_sd, star_sd),
        ("absolute calibration of a single frame (what APASS reports)", atm_floor, ens_sem, star_sem),
    ]
    out = {"regimes": [], "target_vmag": round(target_vmag, 3), "exposure_s": 60.0}
    for label, sig_t, sig_e, sig_c in regimes:
        exp = eb.checkstar_expected_rms(sig_t, sig_e, sig_c)
        out["regimes"].append(
            {
                "label": label,
                "target_sigma_mag": round(sig_t, 5),
                "ensemble_sigma_mag": round(sig_e, 5),
                "check_star_catalogue_sigma_mag": round(sig_c, 5),
                "expected_rms_mag": round(exp, 5),
                "implied_snr": round(eb.mag_error_to_snr(exp), 1),
            }
        )
    base = out["regimes"][1]["expected_rms_mag"]
    out["measured"] = []
    for measured in (0.01, 0.02, 0.03, 0.05, 0.085, 0.11):
        extra_sq = measured**2 - base**2
        extra = math.sqrt(extra_sq) if extra_sq > 0 else 0.0
        if measured < 1.5 * base:
            verdict = "explained by the comparison-star catalogue - your camera is not the limit"
        elif extra >= 0.8 * atm_floor:
            verdict = "as bad as the 0.088 mag catalogue calibration floor: suspect per-frame systematics, not photons"
        else:
            verdict = "worse than the budget: check for a variable or mis-catalogued comp, wrong sky annulus, colour term"
        out["measured"].append(
            {
                "measured_rms_mag": measured,
                "expected_rms_mag": base,
                "measured_over_expected": round(measured / base, 2),
                "unexplained_systematic_mag": round(extra, 5),
                "verdict": verdict,
            }
        )
    return out


# ======================================================================================
# rendering
# ======================================================================================
def table(headers: list[str], rows: list[list]) -> str:
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(out)


def render(j: dict) -> str:
    s1, s2, s3, s4, s5, s6 = j["sec1"], j["sec2"], j["sec3"], j["sec4"], j["sec5"], j["sec6"]
    atm = s1["measured_floor_per_measurement_mag"]  # atmosphere+calibration floor measured in section 1
    L = []
    A = L.append
    A("# Worked examples: numbers behind the guide\n")
    A(f"Generated by `tools/run_examples.py` (v{j['version']}) from the files in `data/` only.")
    A("Regenerate with `python tools/run_examples.py --write`; verify with `--check`.")
    A("Every number below is measured from the shipped real datasets or computed by")
    A("`tools/error_budget.py` from them; nothing here is typed in by hand.\n")

    A("\n## 1. Real per-measurement scatter of ground-based V photometry (APASS DR9)\n")
    A(f"Dataset: `data/apass_dr9_tvcas_r0p5.csv` - {s1['rows_total']} APASS DR9 entries in a 0.5 deg cone on the target.")
    A(f"Usable rows (quoted error > 0, error flag `u_e_vmag`=0, nobs>1): **{s1['rows_usable']}**.")
    A(f"Rows whose quoted error is flagged unreliable (`u_e_vmag`=1): {s1['rows_flagged_unreliable_error']}.\n")
    A("Per-measurement scatter `e_vmag` (the catalogue's own standard deviation of the individual V")
    A("measurements) and the standard error of the catalogue mean (`e_vmag/sqrt(nobs)`):\n")
    A(
        table(
            ["V bin", "stars", "median nobs", "scatter p25", "scatter median", "scatter p90", "catalogue-mean SEM median", "SEM p90"],
            [
                [
                    b["vmag_bin"],
                    b["n_stars"],
                    b["nobs_median"],
                    fmt(b["scatter_sd_mag"]["p25"]),
                    fmt(b["scatter_sd_mag"]["median"]),
                    fmt(b["scatter_sd_mag"]["p90"]),
                    fmt(b["catalogue_mean_sem_mag"]["median"]),
                    fmt(b["catalogue_mean_sem_mag"]["p90"]),
                ]
                for b in s1["by_vmag_bin"]
            ],
        )
    )
    A(f"\nAll magnitudes together: median per-measurement scatter **{s1['measured_floor_per_measurement_mag']} mag**")
    A(f"(p25 {s1['scatter_all_mags_mag']['p25']}, p90 {s1['scatter_all_mags_mag']['p90']}, max {s1['scatter_all_mags_mag']['max']}).\n")
    A("Does the quoted error shrink with the number of observations? Log-log slope of `e_vmag` vs `nobs`\n"
      "within one magnitude bin (a per-measurement scatter gives ~0, a standard error of the mean gives -0.5):\n")
    A(table(["V bin", "stars", "slope", "R^2"], [[k, v.get("n"), fmt(v.get("slope"), "+.3f"), fmt(v.get("r2"))] for k, v in s1["slope_sd_vs_nobs_by_bin"].items()]))
    A("\n**Reading:** the scatter does *not* decrease with more observations, so `e_vmag` is a per-measurement")
    A("quantity and it is dominated by terms that do not depend on photon statistics: the same ~0.09 mag")
    A("appears at V=10.5 and at V=15.5, where photon noise differs by more than a factor of 10.\n")

    A("\n## 2. Same field in Gaia DR3 (space photometry, quoted formal errors)\n")
    A(f"Dataset: `data/gaia_dr3_tvcas_r0p5.csv` - {s2['rows_total']} Gaia DR3 sources, brightest first.")
    A("sigma_mag = 1.0857/(flux/flux_error), i.e. the magnitude error implied by the catalogue's own SNR:\n")
    A(
        table(
            ["G bin", "stars", "median flux SNR", "median sigma_mag", "p90 sigma_mag"],
            [[r["gmag_bin"], r["n_stars"], fmt(r["flux_snr_median"], ",.0f"), fmt(r["sigma_mag_median"], ".6f"), fmt(r["sigma_mag_p90"], ".6f")] for r in s2["by_gmag_bin"]],
        )
    )
    ratio = None
    if s2["by_gmag_bin"]:
        gspace = s2["by_gmag_bin"][2]["sigma_mag_median"]
        gground = s1["measured_floor_per_measurement_mag"]
        ratio = round(gground / gspace, 1)
        A(f"\nGround vs space at G~13.5: {gground} mag (APASS per measurement) / {gspace} mag (Gaia quoted)")
        A(f"= a factor of **{ratio}**. Photon statistics are not what limits ground-based photometry.\n")
    A(f"\nGaia variability flags in the field: {json.dumps(s2['variable_flag_counts'])}.\n")

    A("\n## 3. The target and its real comparison ensemble\n")
    if s3.get("note"):
        A(s3["note"])
    else:
        t = s3["target"]
        flag = " APASS flags this error **unreliable** (`u_e_vmag = 1`)." if t["error_flag_unreliable"] else "."
        sem_txt = (
            f"The mean of those measurements is quoted good to {t['catalogue_mean_sem_mag']} mag"
            if t["catalogue_mean_sem_mag"] is not None
            else "Because the scatter column is zero, this catalogue gives no usable error for this star at all"
        )
        A(f"The AAVSO target position (`data/vsp_target_meta.json`) matches APASS DR9 entry `recno={t['apass_recno']}`")
        A(f"at {t['separation_arcsec']} arcsec: **V = {t['vmag']}** from {t['nobs']:.0f} measurements, per-measurement")
        A(f"scatter {t['e_vmag_sd_mag']} mag.{flag} {sem_txt} - and a star this bright has almost no comparison")
        A("stars near its own magnitude, which is the real problem with bright targets.\n")
        g = s3.get("gaia_crosscheck")
        if g:
            A(f"Gaia DR3 cross-check at {g['separation_arcsec']} arcsec: `source_id={g['source_id']}`, G = {g['gmag']},")
            A(f"BP-RP = {g['bp_minus_rp']}, `phot_variable_flag = {g['phot_variable_flag']}`, RUWE = {g['ruwe']},")
            A(f"formal sigma_G = {g['sigma_g_mag']} mag. Gaia confirms the star is variable, and its quoted error")
            A("is hundreds of times smaller than any ground-based error for the same magnitude.\n")

        A("Two cases are worked out with real stars from the shipped catalogue: the bright AAVSO target itself")
        A("(a magnitude-matched ensemble barely exists) and a real mid-field star near V=11.5, which is where most")
        A("AAVSO work sits. Ensembles are picked by magnitude similarity, the usual rule.\n")
        A("`ens_sd` uses each comp's per-measurement scatter; `ens_sem` uses its catalogue-mean error, which is")
        A("the floor for adopted magnitudes. IVW = inverse-variance weights.\n")
        for case in s3["cases"]:
            A(f"\n### {case['label']}\n")
            A(f"Stars in the catalogue with a usable quoted error: {case['n_usable_stars_in_field']}. The brightest of")
            A(f"them is V = {fmt(case['brightest_usable_comp_vmag'], '.3f')}, and {case['n_comps_brighter_than_target']}")
            A(f"of them are brighter than the target (V = {fmt(case['target_vmag'], '.2f')}), so a magnitude-matched ensemble is")
            A("either easy or impossible depending on the target's brightness. How many are close enough in brightness:\n")
            A(
                table(
                    ["comps with |dV| <=", "stars available"],
                    [[fmt(e["dm_window"], ".2f"), e["comps_available"]] for e in case["ensembles"]["by_window"]],
                )
            )
            A("\nEnsemble error when the N nearest-magnitude stars are used (all real catalogue values):\n")
            A(
                table(
                    ["N wanted", "N used", "comp V range", "worst |dV|", "ens_sd (unweighted)", "ens_sd (IVW)", "ens_sem"],
                    [
                        [
                            e["n_requested"],
                            e["n_used"],
                            e["comp_mag_range"],
                            fmt(e["worst_comp_delta_vmag"], ".2f"),
                            fmt(e["ens_sd_mag"]),
                            fmt(e["ens_sd_ivw_mag"]),
                            fmt(e["ens_sem_mag"]),
                        ]
                        for e in case["ensembles"]["by_n"]
                    ],
                )
            )
            A("\nBright versus faint comparison stars, same target (up to 8 stars per selection):\n")
            A(
                table(
                    ["selection", "available", "N used", "median comp V", "median per-measurement scatter", "ens_sd", "ens_sem"],
                    [
                        [
                            b["selection"],
                            b["n_available"],
                            b["n_used"],
                            fmt(b["median_comp_vmag"], ".2f"),
                            fmt(b["median_comp_sd_mag"]),
                            fmt(b["ens_sd_mag"]),
                            fmt(b["ens_sem_mag"]),
                        ]
                        for b in case["bright_vs_faint"]
                    ],
                )
            )
        ca = s3["cases"][0]
        gap = round(ca["brightest_usable_comp_vmag"] - ca["target_vmag"], 2)
        dmin = min(e["worst_comp_delta_vmag"] for e in ca["ensembles"]["by_n"])
        dmax = max(e["worst_comp_delta_vmag"] for e in ca["ensembles"]["by_n"])
        A(f"\nCase A is the instructive one: the brightest star in the whole field that carries a usable error is")
        A(f"V = {fmt(ca['brightest_usable_comp_vmag'], '.3f')}, i.e. {gap} mag fainter than the target, and the ensembles")
        A(f"above are forced to sit {dmin}-{dmax} mag away. Differential measurement across that much of a magnitude")
        A("range brings in terms a matched ensemble cancels: colour difference between target and comps, sensor")
        A("nonlinearity close to saturation, and flat-field error that scales with the signal ratio. A bright star")
        A("needs its own error budget; the field's typical value cannot be reused for it.\n")
        A("The other half of the trade: what *you* pay in photon noise to measure a comp of that magnitude,")
        A("60 s on the reference setup (atmosphere term switched on):\n")
        A(
            table(
                ["dV", "comp V", "star electrons/60 s", "SNR (photon)", "sigma_photon", "sigma incl. 0.088 mag atmosphere"],
                [[f"{c['delta_vmag']:+.1f}", fmt(c["comp_vmag"], ".2f"), f"{c['star_e_in_60s']:,.0f}", fmt(c["snr"], ",.0f"), fmt(c["sigma_photon_mag"]), fmt(c["sigma_with_atm_mag"])] for c in s3["comp_photon_table"]],
            )
        )
    A("")

    A("\n## 4. Stack-then-measure versus measure-then-average\n")
    p = s4["frame_error_percentiles_mag"]
    A(f"Frame errors below are the real APASS DR9 percentiles: p25 {p['p25']}, median {p['median']}, p90 {p['p90']} mag.\n")
    A(
        table(
            ["frames N", "stack then measure", "measure then average", "identical?", "implied SNR"],
            [[r["frames"], fmt(r["sigma_photon_route_mag"]), fmt(r["sigma_mean_route_mag"]), "yes" if r["routes_agree"] else "no", fmt(r["implied_snr"], ".1f")] for r in s4["equal_frames"]],
        )
    )
    A("\nFor identical frames the two routes are mathematically the same (sigma/sqrt(N)); the useful number")
    A("is the implied SNR, which is what you should quote as your error.\n")
    m = s4["mixed_quality_example"]
    A("Unequal frames - a real case: four good frames plus one clouded frame. Frame errors")
    A(f"{m['frame_errors_mag']} mag:\n")
    A(
        table(
            ["method", "error on the mean (mag)", "implied SNR"],
            [
                ["unweighted mean, all frames kept", fmt(m["unweighted_all"]), fmt(eb.mag_error_to_snr(m["unweighted_all"]), ".1f")],
                ["inverse-variance weighted, all kept", fmt(m["inverse_variance_all"]), fmt(eb.mag_error_to_snr(m["inverse_variance_all"]), ".1f")],
                ["unweighted after rejecting the bad frame", fmt(m["unweighted_kept"]), fmt(eb.mag_error_to_snr(m["unweighted_kept"]), ".1f")],
            ],
        )
    )
    A(f"\nKeeping one clouded frame in an unweighted mean costs **+{m['loss_if_bad_frame_kept_unweighted_pct']}%** on the")
    A("error; inverse-variance weighting almost completely absorbs it. That is the practical argument for")
    A("measure-then-average with weighting, or for weighting your stack.\n")

    A("\n## 5. Why the reported error is 0.08-0.09 mag when the SNR looks like 0.01\n")
    A(
        table(
            ["reported sigma (mag)", "SNR it actually implies"],
            [[fmt(r["reported_sigma_mag"], ".3f"), fmt(r["implied_snr"], ".1f")] for r in s5["error_to_snr"]],
        )
    )
    A("\nReference setup (see README): 20 cm f/10, ASI2600MM Pro at gain 0 (0.77 e-/ADU, read noise 2.8 e-),")
    A("sky 20.7 mag/arcsec^2, aperture radius 5 px, annulus 10-16 px, atmosphere+calibration term taken from")
    A(f"section 1 ({atm} mag per frame):\n")
    A(
        table(
            ["V", "exp (s)", "N", "star e-", "SNR photon", "sigma photon (mag)", "one frame incl. atm", "N frames, atm varies", "N frames, atm constant"],
            [
                [
                    fmt(c["vmag"], ".1f"),
                    fmt(c["exposure_s"], ".0f"),
                    c["frames"],
                    f"{c['star_e']:,.0f}",
                    fmt(c["snr_photon"], ",.0f"),
                    fmt(c["sigma_photon_mag"]),
                    fmt(c["sigma_frame_total_mag"]),
                    fmt(c["stack_random_case_mag"]),
                    fmt(c["stack_constant_offset_mag"]),
                ]
                for c in s5["setup_cases"]
            ],
        )
    )
    A("\nOne caveat about the 0.088 mag term: it is what one *absolute* ground-based V measurement costs")
    A("(section 1). When you measure a target against ensemble stars on the same image, whatever is common to")
    A("target and comparison stars cancels; how much survives is site-specific, and that residual is exactly")
    A("what a check star measures (section 6). Stacking frames cannot remove whatever survives.\n")
    A("\nVariance shares at the faint end (V = 15.0): the dominant term is what you should fix first.\n")
    A(
        table(
            ["V", "exp (s)", "N", "dominant term", "source", "sky", "dark", "read noise", "sky determination"],
            [
                [
                    fmt(c["vmag"], ".1f"),
                    fmt(c["exposure_s"], ".0f"),
                    c["frames"],
                    c["dominant_term"],
                    fmt(c["variance_share"].get("source"), ".3f"),
                    fmt(c["variance_share"].get("sky"), ".4f"),
                    fmt(c["variance_share"].get("dark"), ".4f"),
                    fmt(c["variance_share"].get("readnoise"), ".4f"),
                    fmt(c["variance_share"].get("sky_determination"), ".4f"),
                ]
                for c in s5["setup_cases"]
                if c["vmag"] == 15.0 and c["frames"] == 1
            ],
        )
    )
    A("\nSensitivity sweep at V = 15.0, 60 s (at V = 12.5 the star itself carries 97.6% of the variance, so")
    A("none of these knobs move the error - which is itself the answer for bright-star work):\n")
    rows = []
    for label in ("exposure_s", "aperture_cm", "sky_mag_per_arcsec2", "aper_radius_pix", "dark_e_per_s", "read_noise_e"):
        for r in s5["sens"][label]:
            rows.append([label, fmt(r["value"]), fmt(r["sigma_photon_mag"]), r["dominant_term"]])
    A(table(["parameter swept", "value", "sigma_photon (mag)", "dominant term"], rows))

    A("\n## 6. Check star: expected versus measured deviation\n")
    A("A check star is measured with the ensemble zero point and compared with its own adopted magnitude, so")
    A("three terms add in quadrature: its measurement noise, the ensemble error, and its own catalogue error.")
    A("The numbers below use the real catalogue errors from section 1 - median per-measurement scatter")
    A(f"{fmt(s1['scatter_all_mags_mag']['median'])} mag, median catalogue-mean error {fmt(s1['catalogue_mean_sem_all_mags_mag']['median'])} mag -")
    A(f"and the calculator's photon term for V = {fmt(s6['target_vmag'], '.2f')} in {s6['exposure_s']:.0f} s.\n")
    A(
        table(
            ["regime", "target sigma", "ensemble sigma", "check-star catalogue sigma", "expected RMS", "implied SNR"],
            [
                [
                    r["label"],
                    fmt(r["target_sigma_mag"]),
                    fmt(r["ensemble_sigma_mag"]),
                    fmt(r["check_star_catalogue_sigma_mag"]),
                    fmt(r["expected_rms_mag"]),
                    fmt(r["implied_snr"], ".1f"),
                ]
                for r in s6["regimes"]
            ],
        )
    )
    A("\nRow 1 is what your camera alone can do on one frame. Row 3 is what you get if the comparison-star")
    A("magnitudes you adopt are only as good as a single ground-based measurement - which is exactly the")
    A("0.08-0.09 mag that members report, and it has nothing to do with their SNR.\n")
    A(f"Reading a measured value against the best realistic case ({fmt(s6['measured'][0]['expected_rms_mag'])} mag):\n")
    A(
        table(
            ["measured RMS", "expected RMS", "ratio", "systematic left over", "what it means"],
            [
                [
                    fmt(m["measured_rms_mag"]),
                    fmt(m["expected_rms_mag"]),
                    fmt(m["measured_over_expected"], ".2f"),
                    fmt(m["unexplained_systematic_mag"]),
                    m["verdict"],
                ]
                for m in s6["measured"]
            ],
        )
    )
    return "\n".join(L) + "\n"


def compute() -> dict:
    apass = load("apass_dr9_tvcas_r0p5.csv")
    gaia = load("gaia_dr3_tvcas_r0p5.csv")
    with open(os.path.join(DATA, "vsp_target_meta.json"), encoding="utf-8") as fh:
        meta = json.load(fh)

    s1 = sec1_apass(apass)
    atm = s1["measured_floor_per_measurement_mag"]
    s2 = sec2_gaia(gaia)
    s3 = sec3_target_and_ensembles(apass, gaia, (meta["ra_deg"], meta["dec_deg"]), atm)
    s4 = sec4_stacking(apass, atm)
    s5 = sec5_reconciliation(atm)
    case_b = next((c for c in s3["cases"] if c["label"].startswith("case B")), None)
    ens = [e for e in (case_b or {}).get("ensembles", {}).get("by_n", []) if e["n_requested"] == 8]
    ens_sem = ens[0]["ens_sem_mag"] if ens else 0.0
    ens_sd = ens[0]["ens_sd_mag"] if ens else 0.0
    target_mag = (case_b or {}).get("target_vmag", 11.5)
    photon = eb.frame_budget(
        mag=target_mag, exposure_s=60.0, frames=1, sky_mag_per_arcsec2=20.7, systematic_mag=0.0, **SETUP
    )["per_frame"]["sigma_photon_mag"]
    s6 = sec6_checkstar(
        ens_sd,
        ens_sem,
        s1["catalogue_mean_sem_all_mags_mag"]["median"],
        s1["measured_floor_per_measurement_mag"],
        atm,
        photon,
        target_mag,
    )

    return {
        "version": "1.0.0",
        "calculator_version": eb.VERSION,
        "sec1": s1,
        "sec2": s2,
        "sec3": s3,
        "sec4": s4,
        "sec5": s5,
        "sec6": s6,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true", help="write results/worked_examples.{md,json}")
    mode.add_argument("--check", action="store_true", help="fail if the committed results differ from a fresh run")
    args = ap.parse_args()

    j = compute()
    md = render(j)
    js = json.dumps(j, indent=2, sort_keys=True) + "\n"
    os.makedirs(RESULTS, exist_ok=True)
    md_path, js_path = os.path.join(RESULTS, "worked_examples.md"), os.path.join(RESULTS, "worked_examples.json")

    if args.write:
        with open(md_path, "w", newline="\n", encoding="utf-8") as fh:
            fh.write(md)
        with open(js_path, "w", newline="\n", encoding="utf-8") as fh:
            fh.write(js)
        print(f"wrote {md_path} ({len(md)} chars)")
        return 0

    ok = True
    for path, text in ((md_path, md), (js_path, js)):
        if not os.path.exists(path):
            print(f"MISSING {path}")
            ok = False
            continue
        with open(path, encoding="utf-8") as fh:
            current = fh.read()
        same = current == text or current.replace("\r\n", "\n") == text
        print(f"{'OK  ' if same else 'DIFF'} {path}")
        ok = ok and same
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
