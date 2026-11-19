#!/usr/bin/env python3
"""AAVSO photometry error-budget calculator. Version 1.0.0 (2026-11-19).

Pure Python standard library, no third-party dependencies, no network access.
All arithmetic is exact or documented; nothing in this module is fitted or tuned.

What it computes
----------------
1. Aperture signal-to-noise of one exposure from first principles (electrons):
   source shot noise + sky shot noise + dark current + read noise
   (the standard aperture SNR equation, Howell 1989; see docs/sources.md).
2. Conversion of SNR to magnitude error: sigma_mag = (2.5/ln10)/SNR = 1.085736/SNR.
3. Stacking: N exposures, stack-then-measure vs measure-then-average, unweighted vs
   inverse-variance weighting, and the effect of one bad frame.
4. Comparison-star ensemble calibration error, and the predicted check-star deviation.
5. A variance breakdown that says which term dominates -> the diagnostic checklist.

Conventions
-----------
* Magnitudes are Vega-system V magnitudes unless stated otherwise.
* `zeropoint_e_per_s` = electrons recorded per second by the whole system
  (telescope + optics + filter + sensor + QE) for a V = 0 star. It is computed
  from published constants by `zeropoint_from_flux()`; it can also be measured on
  your own data (that is the recommended route, see the guide).
* Sky brightness is given as V magnitudes per square arcsecond (`sky_mag_per_arcsec2`).
* `sigma_systematic_mag` is an *optional extra variance term* for atmosphere
  (scintillation + transparency) plus calibration: it is NOT invented here - the
  guide measures it from real APASS DR9 data (0.086-0.095 mag per measurement in
  the TV Cas field) and lets the user substitute their own value. Set it to 0 to
  get the pure photon-noise answer.

Usage examples
--------------
  python tools/error_budget.py snr --mag 11.3 --exposure 60 --frames 20 \
      --aperture-cm 20 --focal-mm 2000 --pixel-um 3.76 --sky 20.7 \
      --gain 0.77 --rdnoise 2.8 --dark 0.1 --aper-radius 5 --atm 0.09
  python tools/error_budget.py ensemble --comps data/apass_dr9_tvcas_r0p5.csv \
      --target-mag 11.46 --n 8 --dm 1.5
  python tools/error_budget.py stack --sigmas 0.03,0.03,0.03,0.03,0.25
  python tools/error_budget.py budget --mag 11.3 --exposure 60 --sky 20.7 --atm 0.09
  python tools/error_budget.py selftest
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys

VERSION = "1.0.0"
MAG_PER_NEPER = 2.5 / math.log(10.0)  # 1.0857362...: d(mag) = 2.5/ln10 * d(F)/F


# ----------------------------------------------------------------------------- helpers
def mag_to_flux_ratio(mag: float) -> float:
    """Flux ratio relative to a V = 0 star."""
    return 10.0 ** (-0.4 * mag)


def snr_to_mag_error(snr: float) -> float:
    """1-sigma magnitude error of a measurement with signal-to-noise `snr`."""
    if snr <= 0:
        raise ValueError("snr must be > 0")
    return MAG_PER_NEPER / snr


def mag_error_to_snr(sigma_mag: float) -> float:
    """The SNR that a reported magnitude error actually implies (the reality check)."""
    if sigma_mag <= 0:
        raise ValueError("sigma_mag must be > 0")
    return MAG_PER_NEPER / sigma_mag


def zeropoint_from_flux(
    photon_flux_per_angstrom: float = 970.0,
    bandwidth_angstrom: float = 880.0,
    aperture_cm: float = 20.0,
    obstruction: float = 0.15,
    throughput: float = 0.45,
) -> float:
    """Electrons recorded per second for a V = 0 star (documented in docs/sources.md).

    photon_flux_per_angstrom: Vega photon flux density at the V effective wavelength
    bandwidth_angstrom:        Johnson V effective width
    obstruction:               linear fraction of the *diameter* blocked by secondary + spider,
                               so the clear area is pi*(D/2)^2 * (1 - obstruction^2)
    throughput:                total (mirror x filter x window x QE), electrons per incident photon
    """
    area_cm2 = math.pi * (aperture_cm / 2.0) ** 2
    clear = area_cm2 * (1.0 - obstruction**2)
    return photon_flux_per_angstrom * bandwidth_angstrom * clear * throughput


def sky_e_per_pixel_per_s(sky_mag_per_arcsec2: float, zeropoint_e_per_s: float, pixel_area_arcsec2: float) -> float:
    """Sky background in e-/pixel/s for a sky of `sky_mag_per_arcsec2` mag/arcsec^2."""
    return zeropoint_e_per_s * mag_to_flux_ratio(sky_mag_per_arcsec2) * pixel_area_arcsec2


def pixel_area_arcsec2(pixel_um: float, focal_mm: float) -> float:
    """Plate scale area per pixel in arcsec^2 (206.265 arcsec per radian)."""
    scale = 206.2648 * pixel_um / focal_mm
    return scale * scale


def exposure_electrons(
    mag: float,
    exposure_s: float,
    zeropoint_e_per_s: float,
    extinction_mag: float = 0.0,
) -> float:
    """Star electrons above background recorded in the *whole star* (not just aperture)."""
    return zeropoint_e_per_s * mag_to_flux_ratio(mag + extinction_mag) * exposure_s


def aperture_snr(
    star_e: float,
    sky_e_per_pix_s: float,
    exposure_s: float,
    dark_e_per_pix_s: float,
    read_noise_e: float,
    aper_pixels: float,
    annulus_pixels: float = 0.0,
    n_dark_frames: float = 0.0,
    n_flat_frames: float = 0.0,
) -> float:
    """Standard aperture SNR (Howell 1989 / Merline & Howell 1995 form), all terms in electrons.

    `aper_pixels` is the area of the photometry aperture in pixels; `annulus_pixels` the
    area of the sky annulus (used to propagate the uncertainty of the sky estimate).
    Star flux in the aperture must already be aperture-corrected (pass it as `star_e`).
    """
    src = star_e
    var_pix = sky_e_per_pix_s * exposure_s + dark_e_per_pix_s * exposure_s + read_noise_e**2
    sky = aper_pixels * sky_e_per_pix_s * exposure_s
    dark = aper_pixels * dark_e_per_pix_s * exposure_s
    read = aper_pixels * read_noise_e**2
    var = src + sky + dark + read
    if annulus_pixels > 0:
        # the mean sky per aperture pixel comes from `annulus_pixels` pixels, and is multiplied
        # by aper_pixels in the aperture sum -> its variance enters as (n_ap/n_ann) * n_ap * var_pix
        var += (aper_pixels**2 / annulus_pixels) * var_pix
    if n_dark_frames > 0:
        var += (aper_pixels**2 / n_dark_frames) * (dark + read)
    if n_flat_frames > 0 and src > 0:
        var += src**2 / n_flat_frames  # photon-transfer noise of a finite flat set
    return src / math.sqrt(var) if var > 0 else float("inf")


def variance_terms(
    star_e: float,
    sky_e_per_pix_s: float,
    exposure_s: float,
    dark_e_per_pix_s: float,
    read_noise_e: float,
    aper_pixels: float,
    annulus_pixels: float = 0.0,
) -> dict:
    """Variance contributions to the aperture sum (electrons^2), for diagnosis."""
    out = {
        "source": star_e,
        "sky": aper_pixels * sky_e_per_pix_s * exposure_s,
        "dark": aper_pixels * dark_e_per_pix_s * exposure_s,
        "readnoise": aper_pixels * read_noise_e**2,
    }
    if annulus_pixels > 0:
        var_pix = sky_e_per_pix_s * exposure_s + dark_e_per_pix_s * exposure_s + read_noise_e**2
        out["sky_determination"] = (aper_pixels**2 / annulus_pixels) * var_pix
    return out


# ----------------------------------------------------------------------------- stacking
def stack_equal_frames(sigma_frame_mag: float, n_frames: int) -> float:
    """Error of a stack/mean of N identical frames: sigma/sqrt(N) (both routes agree)."""
    if n_frames <= 0:
        raise ValueError("n_frames must be > 0")
    return sigma_frame_mag / math.sqrt(n_frames)


def mean_unweighted(sigmas: list[float]) -> float:
    """Error of the plain (unweighted) mean of N measurements with individual errors."""
    n = len(sigmas)
    if n == 0:
        raise ValueError("no frames")
    return math.sqrt(sum(s * s for s in sigmas)) / n


def mean_inverse_variance(sigmas: list[float]) -> float:
    """Error of the inverse-variance weighted mean = the error of an optimally built stack."""
    w = [1.0 / (s * s) for s in sigmas]
    return 1.0 / math.sqrt(sum(w))


def stack_snr_from_frame_snr(frame_snr: list[float]) -> float:
    """SNR of a summed stack of frames with different SNRs: sqrt(sum(SNR_i^2))."""
    return math.sqrt(sum(s * s for s in frame_snr))


def outlier_penalty(sigmas: list[float], keep: list[bool] | None = None) -> dict:
    """Compare keeping every frame with rejecting outliers, on the same data."""
    if keep is None:
        med = sorted(sigmas)[len(sigmas) // 2]
        keep = [s <= 2.0 * med for s in sigmas]
    good = [s for s, k in zip(sigmas, keep) if k]
    return {
        "n_total": len(sigmas),
        "n_kept": len(good),
        "unweighted_all": mean_unweighted(sigmas),
        "unweighted_kept": mean_unweighted(good),
        "inverse_variance_all": mean_inverse_variance(sigmas),
        "inverse_variance_kept": mean_inverse_variance(good),
    }


# ----------------------------------------------------------------------------- ensembles
def ensemble_calibration_error(comp_errors: list[float], weights: list[float] | None = None) -> float:
    """Error of the ensemble zero point built from comparison stars.

    The calibration constant is a weighted mean of the per-comparison offsets, so for
    unweighted use it is sqrt(sum(sigma_i^2))/N; with inverse-variance weights it reduces
    to 1/sqrt(sum(1/sigma_i^2)) (the same arithmetic as `mean_inverse_variance`).
    """
    if not comp_errors:
        raise ValueError("no comparison stars")
    if weights is None:
        n = len(comp_errors)
        return math.sqrt(sum(e * e for e in comp_errors)) / n
    wsum = sum(weights)
    return math.sqrt(sum((w / wsum) ** 2 * e * e for w, e in zip(weights, comp_errors)))


def measurement_total_error(sigma_frame_mag: float, sigma_ensemble_mag: float) -> float:
    """Error on the calibrated magnitude of one target measurement."""
    return math.sqrt(sigma_frame_mag**2 + sigma_ensemble_mag**2)


def checkstar_expected_rms(
    sigma_target_mag: float, sigma_ensemble_mag: float, sigma_check_catalogue_mag: float | None = None
) -> float:
    """Expected RMS of a check star's measured minus adopted magnitude.

    Three independent terms add in quadrature: the check star's own measurement noise, the
    ensemble zero point it is calibrated with, and its own adopted catalogue magnitude error.
    The ensemble term does NOT cancel in the difference - it enters with the opposite sign.
    If `sigma_check_catalogue_mag` is omitted it is taken equal to the ensemble error, i.e. a
    check star excluded from the ensemble but of comparable quality.
    """
    sc = sigma_ensemble_mag if sigma_check_catalogue_mag is None else sigma_check_catalogue_mag
    return math.sqrt(sigma_target_mag**2 + sigma_ensemble_mag**2 + sc**2)


def read_comp_csv(path: str, mag_col: str, err_col: str, nobs_col: str | None = None) -> list[dict]:
    """Read comparison stars (magnitude + error) from a catalogue CSV.

    `err_col` is the *standard deviation of the individual measurements* in catalogues such as
    APASS DR9; pass `nobs_col` to convert it to the standard error of the catalogue mean
    (err / sqrt(nobs)), which is the quantity that limits an adopted comparison magnitude.
    """
    rows = []
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            try:
                mag = float(row[mag_col])
                err = float(row[err_col])
                if nobs_col:
                    n = float(row[nobs_col])
                    if n > 1:
                        err = err / math.sqrt(n)
            except (KeyError, TypeError, ValueError):
                continue
            if err > 0:
                rows.append({"mag": mag, "err": err})
    return rows


# ----------------------------------------------------------------------------- top level
def frame_budget(
    mag: float,
    exposure_s: float,
    frames: int = 1,
    aperture_cm: float = 20.0,
    focal_mm: float = 2000.0,
    pixel_um: float = 3.76,
    gain_e_per_adu: float = 0.77,
    read_noise_e: float = 2.8,
    dark_e_per_s: float = 0.1,
    sky_mag_per_arcsec2: float = 20.7,
    aper_radius_pix: float = 5.0,
    annulus_inner_pix: float = 10.0,
    annulus_outer_pix: float = 16.0,
    aperture_fraction: float = 1.0,
    extinction_mag: float = 0.0,
    zeropoint_e_per_s: float | None = None,
    systematic_mag: float = 0.0,
) -> dict:
    """Full per-frame and stacked error budget for one target + setup."""
    zp = zeropoint_e_per_s or zeropoint_from_flux(aperture_cm=aperture_cm)
    parea = pixel_area_arcsec2(pixel_um, focal_mm)
    skypix = sky_e_per_pixel_per_s(sky_mag_per_arcsec2, zp, parea)
    star_total = exposure_electrons(mag, exposure_s, zp, extinction_mag)
    star = star_total * aperture_fraction
    n_aper = math.pi * aper_radius_pix**2
    n_ann = math.pi * (annulus_outer_pix**2 - annulus_inner_pix**2)

    snr1 = aperture_snr(star, skypix, exposure_s, dark_e_per_s, read_noise_e, n_aper, n_ann)
    sig_photon = snr_to_mag_error(snr1)
    sig_frame = math.sqrt(sig_photon**2 + systematic_mag**2)

    terms = variance_terms(star, skypix, exposure_s, dark_e_per_s, read_noise_e, n_aper, n_ann)
    total = sum(terms.values())
    share = {k: (v / total if total else 0.0) for k, v in sorted(terms.items(), key=lambda kv: -kv[1])}

    # Stacking N frames. Two cases for the atmosphere+calibration term, because that is the
    # question members ask: it averages down only if it changes between frames.
    n = max(1, int(frames))
    snr_stack = snr1 * math.sqrt(n)
    sig_photon_stack = sig_photon / math.sqrt(n)
    sig_frame_total = sig_frame  # sqrt(photon^2 + systematic^2) for one frame
    case_random = sig_frame_total / math.sqrt(n)  # every frame independent (best case)
    case_systematic = math.sqrt(sig_photon_stack**2 + systematic_mag**2)  # offset common to the run
    return {
        "inputs": {
            "version": VERSION,
            "mag": mag,
            "exposure_s": exposure_s,
            "frames": frames,
            "aperture_cm": aperture_cm,
            "focal_mm": focal_mm,
            "pixel_um": pixel_um,
            "arcsec_per_pixel": round(math.sqrt(parea), 4),
            "gain_e_per_adu": gain_e_per_adu,
            "read_noise_e": read_noise_e,
            "dark_e_per_s": dark_e_per_s,
            "sky_mag_per_arcsec2": sky_mag_per_arcsec2,
            "sky_e_per_pix_s": round(skypix, 4),
            "aperture_radius_pix": aper_radius_pix,
            "aperture_pixels": round(n_aper, 2),
            "annulus_pixels": round(n_ann, 2),
            "aperture_fraction": aperture_fraction,
            "extinction_mag": extinction_mag,
            "zeropoint_e_per_s": round(zp, 3),
            "systematic_mag": systematic_mag,
        },
        "electrons": {
            "star_total_in_exposure": round(star_total, 1),
            "star_in_aperture": round(star, 1),
            "sky_in_aperture": round(skypix * n_aper * exposure_s, 1),
            "dark_in_aperture": round(dark_e_per_s * n_aper * exposure_s, 1),
            "read_noise_in_aperture": round(math.sqrt(n_aper) * read_noise_e, 1),
            "star_adu_in_aperture": round(star / gain_e_per_adu, 1),
        },
        "per_frame": {
            "snr": round(snr1, 2),
            "sigma_photon_mag": round(sig_photon, 5),
            "sigma_frame_total_mag": round(sig_frame, 5),
            "implied_snr_from_reported_error": round(mag_error_to_snr(sig_frame), 2),
        },
        "stacked": {
            "frames": n,
            "snr_photon": round(snr_stack, 2),
            "sigma_photon_mag": round(sig_photon_stack, 5),
            "if_atm_changes_between_frames_mag": round(case_random, 5),
            "if_atm_constant_over_run_mag": round(case_systematic, 5),
        },
        "variance_share": {k: round(v, 4) for k, v in share.items()},
        "dominant_term": next(iter(share)),
    }


# ----------------------------------------------------------------------------- self test
def _selftest() -> int:
    """Assertive checks of every formula above; run with `error_budget.py selftest`."""
    ok = 0

    def check(name, cond):
        nonlocal ok
        assert cond, f"FAILED: {name}"
        ok += 1

    check("mag-per-neper constant", abs(MAG_PER_NEPER - 1.085736) < 1e-6)
    check("snr->mag roundtrip", abs(snr_to_mag_error(mag_error_to_snr(0.05)) - 0.05) < 1e-12)
    check("SNR 100 -> 0.0109 mag", abs(snr_to_mag_error(100.0) - 0.0108574) < 1e-6)
    check("error 0.085 implies SNR 12.8", abs(mag_error_to_snr(0.085) - 12.7734) < 1e-3)

    # background-limited: SNR grows as sqrt(exposure time)
    a = aperture_snr(1000.0, 50.0, 60.0, 0.0, 0.0, 78.5)
    b = aperture_snr(1000.0 * 4, 50.0, 240.0, 0.0, 0.0, 78.5)
    check("background-limited sqrt(t)", abs(b / a - 2.0) < 1e-9)

    # read-noise-limited: SNR grows almost linearly with time (source shot noise negligible)
    c = aperture_snr(1.0, 0.0, 60.0, 0.0, 10.0, 78.5)
    d = aperture_snr(4.0, 0.0, 240.0, 0.0, 10.0, 78.5)
    check("read-noise-limited ~linear in t", abs(d / c - 4.0) < 0.01)

    # source shot noise only: SNR = sqrt(N)
    check("pure shot noise", abs(aperture_snr(10000.0, 0.0, 1.0, 0.0, 0.0, 1.0) - 100.0) < 1e-9)

    # stacking: equal frames
    check("stack sqrt(N)", abs(stack_equal_frames(0.09, 20) - 0.09 / math.sqrt(20)) < 1e-15)

    # inverse-variance weighting is never worse than the unweighted mean
    s = [0.02, 0.03, 0.05, 0.2, 0.04]
    check("IVW <= unweighted", mean_inverse_variance(s) < mean_unweighted(s))
    check("IVW == stack SNR route", abs(mean_inverse_variance([0.1, 0.1]) - 0.1 / math.sqrt(2)) < 1e-15)
    pen = outlier_penalty(s)
    check("rejecting the bad frame helps", pen["unweighted_kept"] < pen["unweighted_all"])

    # stack SNR from per-frame SNRs equals the IVW result for equal-quality frames
    check("stack SNR route", abs(stack_snr_from_frame_snr([50, 50, 50]) - 50 * math.sqrt(3)) < 1e-9)

    # ensemble: identical comps -> error falls as 1/sqrt(N); doubling N cuts by sqrt2
    e4 = ensemble_calibration_error([0.09] * 4)
    e16 = ensemble_calibration_error([0.09] * 16)
    check("ensemble 1/sqrt(N)", abs(e4 / e16 - math.sqrt(4)) < 1e-12)
    check("ensemble N=1 identity", abs(ensemble_calibration_error([0.09]) - 0.09) < 1e-15)

    # one bad comparison star inflates the ensemble far above the good-star limit
    good = ensemble_calibration_error([0.09] * 8)
    mixed = ensemble_calibration_error([0.09] * 7 + [0.6])
    check("one bad comp dominates", mixed > 2.0 * good)

    # check star: ensemble term enters with 2x variance, so RMS > target error alone
    check("checkstar rms composition", abs(checkstar_expected_rms(0.03, 0.02) - math.sqrt(0.03**2 + 2 * 0.02**2)) < 1e-15)
    check(
        "checkstar rms uses the check star's own catalogue error",
        abs(checkstar_expected_rms(0.03, 0.02, 0.04) - math.sqrt(0.03**2 + 0.02**2 + 0.04**2)) < 1e-15,
    )

    # plate scale: 2.0 um on a 1000 mm camera -> 0.4125296 arcsec/pixel (206.2648"/rad)
    check("plate scale", abs(math.sqrt(pixel_area_arcsec2(2.0, 1000.0)) - 0.4125296) < 1e-6)

    # zeropoint: documented default value is reproducible
    expected_zp = 970.0 * 880.0 * (math.pi * 10.0**2) * (1.0 - 0.15**2) * 0.45
    check("zeropoint default", abs(zeropoint_from_flux() - expected_zp) < 1e-6)

    # a fainter target must have a larger error, and more frames must reduce it
    f1 = frame_budget(mag=13.0, exposure_s=60.0, frames=1)
    f2 = frame_budget(mag=11.0, exposure_s=60.0, frames=1)
    check("brighter is better", f2["per_frame"]["sigma_photon_mag"] < f1["per_frame"]["sigma_photon_mag"])
    f3 = frame_budget(mag=13.0, exposure_s=60.0, frames=16)
    check("stacking helps", f3["stacked"]["sigma_photon_mag"] < f1["per_frame"]["sigma_photon_mag"])
    check(
        "stack photon scaling exact",
        abs(f3["stacked"]["sigma_photon_mag"] - f1["per_frame"]["sigma_photon_mag"] / 4.0) < 1e-6,
    )

    # sky brightness must matter: brighter sky -> larger error for the same star
    g1 = frame_budget(mag=13.0, exposure_s=60.0, sky_mag_per_arcsec2=21.7)
    g2 = frame_budget(mag=13.0, exposure_s=60.0, sky_mag_per_arcsec2=19.5)
    check("brighter sky worse", g2["per_frame"]["sigma_photon_mag"] > g1["per_frame"]["sigma_photon_mag"])

    print(f"selftest PASS: {ok} checks")
    return 0


# ----------------------------------------------------------------------------- CLI
def _load_comps(args):
    if args.comps:
        nobs_col = args.nobs_col if args.err_mode == "sem" else None
        rows = read_comp_csv(args.comps, args.mag_col, args.err_col, nobs_col)
    else:
        rows = [{"mag": float(p.split(":")[0]), "err": float(p.split(":")[1])} for p in args.comp_list]
    if args.target_mag is not None and args.dm:
        rows = [r for r in rows if abs(r["mag"] - args.target_mag) <= args.dm]
    rows.sort(key=lambda r: (abs(r["mag"] - (args.target_mag or 0.0))))
    return rows[: args.n]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0] + f" (v{VERSION})")
    sub = ap.add_subparsers(dest="cmd", required=True)

    common = argparse.ArgumentParser(add_help=False)
    for flag, kw in [
        ("--mag", dict(type=float, default=11.3, help="target V magnitude")),
        ("--exposure", dict(type=float, default=60.0, help="single exposure length, seconds")),
        ("--frames", dict(type=int, default=1, help="number of exposures to stack/average")),
        ("--aperture-cm", dict(type=float, default=20.0, help="telescope aperture, cm (diameter)")),
        ("--focal-mm", dict(type=float, default=2000.0, help="image scale focal length, mm")),
        ("--pixel-um", dict(type=float, default=3.76, help="pixel size, microns")),
        ("--gain", dict(type=float, default=0.77, help="system gain, e-/ADU")),
        ("--rdnoise", dict(type=float, default=2.8, help="read noise, electrons/pixel")),
        ("--dark", dict(type=float, default=0.1, help="dark current, e-/pixel/s")),
        ("--sky", dict(type=float, default=20.7, help="sky brightness, V mag/arcsec^2")),
        ("--aper-radius", dict(type=float, default=5.0, help="aperture radius, pixels")),
        ("--annulus", dict(type=float, nargs=2, default=[10.0, 16.0], metavar=("IN", "OUT"))),
        ("--aper-frac", dict(type=float, default=1.0, help="fraction of the star flux inside the aperture")),
        ("--extinction", dict(type=float, default=0.0, help="extinction at your airmass, mag")),
        ("--zeropoint", dict(type=float, default=None, help="override e-/s for V=0 (measured on your data)")),
        ("--atm", dict(type=float, default=0.0, help="per-frame atmosphere+calibration floor, mag (see guide)")),
    ]:
        common.add_argument(flag, **kw)

    p_snr = sub.add_parser("snr", parents=[common], help="per-frame and stacked SNR / magnitude error")
    p_budget = sub.add_parser("budget", parents=[common], help="which noise term dominates")
    p_sens = sub.add_parser("sens", parents=[common], help="sweep the parameters that could fix it")

    p_stack = sub.add_parser("stack", help="measure-then-average vs stack-then-measure")
    p_stack.add_argument("--sigmas", default=None, help="comma-separated per-frame magnitude errors")
    p_stack.add_argument("--frame-sigma", type=float, default=None, help="identical frames: error of one frame")
    p_stack.add_argument("--n", type=int, default=None, help="identical frames: how many")
    p_stack.add_argument("--keep", default=None, help="optional 0/1 mask, e.g. 1,1,1,1,0")

    p_ens = sub.add_parser("ensemble", help="comparison-star ensemble calibration error")
    p_ens.add_argument("--comps", help="CSV with a magnitude column and an error column")
    p_ens.add_argument("--mag-col", default="vmag")
    p_ens.add_argument("--err-col", default="e_vmag")
    p_ens.add_argument(
        "--err-mode",
        choices=["sd", "sem"],
        default="sd",
        help="'sd' = the catalogue error is the scatter of one measurement (APASS DR9 default); "
        "'sem' = convert it to the error of the catalogue mean using --nobs-col",
    )
    p_ens.add_argument("--nobs-col", default="nobs", help="number of measurements behind each catalogue magnitude")
    p_ens.add_argument("--comp-list", nargs="*", default=[], help="fallback: mag:err pairs")
    p_ens.add_argument("--target-mag", type=float, default=None)
    p_ens.add_argument("--dm", type=float, default=0.0, help="keep comps within this magnitude difference")
    p_ens.add_argument("--n", type=int, default=8)

    p_chk = sub.add_parser("checkstar", help="expected check-star deviation RMS")
    p_chk.add_argument("--target-sigma", type=float, required=True)
    p_chk.add_argument("--ensemble-sigma", type=float, required=True)
    p_chk.add_argument("--measured-rms", type=float, default=None, help="your measured value -> verdict")
    p_chk.add_argument(
        "--check-cat-sigma",
        type=float,
        default=None,
        help="the check star's own adopted catalogue magnitude error; defaults to the ensemble sigma",
    )

    sub.add_parser("selftest", help="run the built-in formula checks")
    for p in (p_snr, p_budget, p_sens):
        p.add_argument("--json", action="store_true")

    args = ap.parse_args(argv)

    if args.cmd == "selftest":
        return _selftest()

    def kwargs(extra=None):
        kw = dict(
            mag=args.mag,
            exposure_s=args.exposure,
            frames=args.frames,
            aperture_cm=args.aperture_cm,
            focal_mm=args.focal_mm,
            pixel_um=args.pixel_um,
            gain_e_per_adu=args.gain,
            read_noise_e=args.rdnoise,
            dark_e_per_s=args.dark,
            sky_mag_per_arcsec2=args.sky,
            aper_radius_pix=args.aper_radius,
            annulus_inner_pix=args.annulus[0],
            annulus_outer_pix=args.annulus[1],
            aperture_fraction=args.aper_frac,
            extinction_mag=args.extinction,
            zeropoint_e_per_s=args.zeropoint,
            systematic_mag=getattr(args, "atm", 0.0),
        )
        if extra:
            kw.update(extra)
        return kw

    if args.cmd in ("snr", "budget"):
        res = frame_budget(**kwargs())
        print(json.dumps(res, indent=2))
        return 0

    if args.cmd == "sens":
        base = frame_budget(**kwargs())
        out = {"base": base["per_frame"], "sweeps": {}}
        for label, key, values in [
            ("exposure_s", "exposure_s", [10, 30, 60, 120, 300]),
            ("frames", "frames", [1, 4, 9, 16, 25]),
            ("aperture_cm", "aperture_cm", [8, 10, 15, 20, 30]),
            ("sky_mag_per_arcsec2", "sky_mag_per_arcsec2", [19.5, 20.4, 21.3, 21.8]),
            ("aper_radius_pix", "aper_radius_pix", [3, 4, 5, 7, 10]),
        ]:
            rows = []
            for v in values:
                r = frame_budget(**kwargs({key: v}))
                rows.append(
                    {
                        "value": v,
                        "snr": r["per_frame"]["snr"] if key != "frames" else r["stacked"]["snr_photon"],
                        "sigma_photon_mag": (r["stacked"] if key == "frames" else r["per_frame"])["sigma_photon_mag"],
                        "dominant_term": r["dominant_term"],
                    }
                )
            out["sweeps"][label] = rows
        print(json.dumps(out, indent=2))
        return 0

    if args.cmd == "stack":
        if args.sigmas:
            sigmas = [float(x) for x in args.sigmas.split(",")]
        elif args.frame_sigma is not None and args.n:
            sigmas = [args.frame_sigma] * args.n
        else:
            print("stack needs --sigmas, or --frame-sigma together with --n", file=sys.stderr)
            return 2
        keep = None
        if args.keep:
            keep = [bool(int(x)) for x in args.keep.split(",")]
        res = outlier_penalty(sigmas, keep)
        res["frame_snr_equivalent"] = {k: round(mag_error_to_snr(v), 2) for k, v in res.items() if isinstance(v, float)}
        print(json.dumps(res, indent=2))
        return 0

    if args.cmd == "ensemble":
        comps = _load_comps(args)
        if not comps:
            print("no comparison stars selected", file=sys.stderr)
            return 2
        errors = [c["err"] for c in comps]
        ens = ensemble_calibration_error(errors)
        res = {
            "error_mode": getattr(args, "err_mode", "sd"),
            "n_comps": len(comps),
            "comp_mags": [round(c["mag"], 3) for c in comps],
            "comp_errors_mag": [round(c["err"], 4) for c in comps],
            "ensemble_error_unweighted_mag": round(ens, 5),
            "ensemble_error_inverse_variance_mag": round(mean_inverse_variance(errors), 5),
            "floor_if_every_comp_as_good_as_the_best_one_mag": round(min(errors) / math.sqrt(len(comps)), 5),
        }
        print(json.dumps(res, indent=2))
        return 0

    if args.cmd == "checkstar":
        exp = checkstar_expected_rms(args.target_sigma, args.ensemble_sigma, args.check_cat_sigma)
        res = {"expected_rms_mag": round(exp, 5)}
        if args.measured_rms is not None:
            ratio = args.measured_rms / exp
            res["measured_rms_mag"] = args.measured_rms
            res["measured_over_expected"] = round(ratio, 2)
            res["verdict"] = (
                "consistent with the budget (look at photon noise first)"
                if ratio < 1.5
                else "systematic problem: a bad/variable comparison star or a wrong sky/aperture setup"
            )
        print(json.dumps(res, indent=2))
        return 0

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
