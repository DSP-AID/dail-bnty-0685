#!/usr/bin/env python3
"""Formula tests for tools/error_budget.py plus consistency checks on the shipped data.

    python -m unittest discover -s tests -v      (or: python tests/test_error_budget.py)

Standard library only, no network. The data file used by the last test class ships in
`data/`, so these tests run on a clean machine with no API keys.
"""
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))

import error_budget as eb  # noqa: E402

DATA_CSV = os.path.join(ROOT, "data", "apass_dr9_tvcas_r0p5.csv")


def median(values):
    return sorted(values)[len(values) // 2]


class TestMagSnr(unittest.TestCase):
    def test_constant(self):
        self.assertAlmostEqual(eb.MAG_PER_NEPER, 1.085736, places=6)

    def test_round_trip(self):
        for snr in (5, 50, 460, 6170):
            self.assertAlmostEqual(eb.mag_error_to_snr(eb.snr_to_mag_error(snr)), snr, places=9)

    def test_exact_value(self):
        self.assertAlmostEqual(eb.snr_to_mag_error(100.0), 0.01085736, places=8)

    def test_rejects_nonpositive(self):
        for bad in (0.0, -3.0):
            with self.assertRaises(ValueError):
                eb.snr_to_mag_error(bad)
            with self.assertRaises(ValueError):
                eb.mag_error_to_snr(bad)


class TestFluxAndScale(unittest.TestCase):
    def test_flux_ratio(self):
        self.assertAlmostEqual(eb.mag_to_flux_ratio(5.0), 0.01, places=12)
        self.assertAlmostEqual(eb.mag_to_flux_ratio(-2.5), 10.0, places=12)

    def test_zeropoint_scales_with_area(self):
        # doubling the diameter must quadruple the light grasp
        ratio = eb.zeropoint_from_flux(aperture_cm=40.0) / eb.zeropoint_from_flux(aperture_cm=20.0)
        self.assertAlmostEqual(ratio, 4.0, places=12)

    def test_zeropoint_exact(self):
        zp = eb.zeropoint_from_flux(
            photon_flux_per_angstrom=1000.0, bandwidth_angstrom=900.0, aperture_cm=10.0, obstruction=0.2, throughput=0.5
        )
        clear = math.pi * 5.0**2 * (1.0 - 0.04)
        self.assertAlmostEqual(zp, 1000.0 * 900.0 * clear * 0.5, delta=1.0)

    def test_pixel_area(self):
        # 3.76 um on a 2000 mm scope is 0.3878 arcsec/pixel
        self.assertAlmostEqual(math.sqrt(eb.pixel_area_arcsec2(3.76, 2000.0)), 0.387778, places=5)

    def test_sky_e_per_pixel(self):
        zp = 1.0e8
        area = 0.25
        got = eb.sky_e_per_pixel_per_s(21.0, zp, area)
        self.assertAlmostEqual(got, zp * 10.0 ** (-0.4 * 21.0) * area, places=9)

    def test_extinction(self):
        zp = 1.0e8
        a = eb.exposure_electrons(12.0, 60.0, zp)
        b = eb.exposure_electrons(12.0, 60.0, zp, extinction_mag=5.0)
        self.assertAlmostEqual(a / b, 100.0, places=9)


class TestSnrModel(unittest.TestCase):
    def test_pure_poisson(self):
        # a bare Poisson source must give sqrt(N) noise exactly
        self.assertAlmostEqual(eb.aperture_snr(10000.0, 0.0, 10.0, 0.0, 0.0, 10.0), 100.0, places=9)

    def test_hand_case_with_sky(self):
        # var = star + n_ap * sky_per_pix (no dark, no read noise, no annulus term)
        snr = eb.aperture_snr(1000.0, 1.0, 10.0, 0.0, 0.0, 10.0)
        self.assertAlmostEqual(snr, 1000.0 / math.sqrt(1100.0), places=9)

    def test_snr_matches_variance_terms(self):
        args = (4466.0, 0.93, 60.0, 0.1, 2.8, 78.5, 659.7)
        snr = eb.aperture_snr(*args)
        terms = eb.variance_terms(*args)
        self.assertAlmostEqual(snr, args[0] / math.sqrt(sum(terms.values())), places=9)

    def test_more_flats_reduce_noise(self):
        base = eb.aperture_snr(5000.0, 1.0, 60.0, 0.2, 3.0, 80.0, 600.0)
        withflats = eb.aperture_snr(5000.0, 1.0, 60.0, 0.2, 3.0, 80.0, 600.0, n_dark_frames=20.0, n_flat_frames=100.0)
        self.assertLess(withflats, base)

    def test_annulus_term_shrinks_with_more_annulus_pixels(self):
        small = eb.aperture_snr(5000.0, 1.0, 60.0, 0.2, 3.0, 80.0, 200.0)
        big = eb.aperture_snr(5000.0, 1.0, 60.0, 0.2, 3.0, 80.0, 2000.0)
        self.assertGreater(big, small)


class TestStacking(unittest.TestCase):
    def test_equal_frames(self):
        self.assertAlmostEqual(eb.stack_equal_frames(0.09, 9), 0.03, places=12)

    def test_two_routes_agree_for_equal_frames(self):
        # stack-then-measure == measure-then-average when the frames are identical
        sig = 0.09
        self.assertAlmostEqual(eb.mean_unweighted([sig] * 16), eb.stack_equal_frames(sig, 16), places=12)

    def test_inverse_variance_beats_or_equals_unweighted(self):
        s = [0.02, 0.03, 0.09, 0.25]
        self.assertLessEqual(eb.mean_inverse_variance(s), eb.mean_unweighted(s))

    def test_inverse_variance_equals_unweighted_when_equal(self):
        s = [0.04] * 6
        self.assertAlmostEqual(eb.mean_inverse_variance(s), eb.mean_unweighted(s) , places=12)

    def test_rejecting_a_bad_frame_helps(self):
        s = [0.02, 0.02, 0.03, 0.03, 0.30]
        p = eb.outlier_penalty(s)
        self.assertAlmostEqual(p["unweighted_all"], eb.mean_unweighted(s), places=12)
        self.assertLess(p["unweighted_kept"], p["unweighted_all"])
        self.assertLessEqual(p["inverse_variance_all"], p["unweighted_all"])

    def test_stack_snr_from_frame_snr(self):
        # combining 4 frames with SNR 100, 120, 80, 90 is quadrature in SNR space
        want = math.sqrt(100**2 + 120**2 + 80**2 + 90**2)
        self.assertAlmostEqual(eb.stack_snr_from_frame_snr([100.0, 120.0, 80.0, 90.0]), want, places=9)


class TestEnsemble(unittest.TestCase):
    def test_equal_errors_give_sem(self):
        e = [0.05] * 4
        self.assertAlmostEqual(eb.ensemble_calibration_error(e), 0.05 / 2.0, places=12)

    def test_more_comps_help(self):
        e = [0.088] * 3
        f = [0.088] * 12
        self.assertLess(eb.ensemble_calibration_error(f), eb.ensemble_calibration_error(e))

    def test_inverse_variance_prefer_lighter_errors(self):
        errs = [0.02, 0.02, 0.09, 0.09]
        self.assertLessEqual(eb.mean_inverse_variance(errs), eb.ensemble_calibration_error(errs))

    def test_total_and_checkstar_composition(self):
        self.assertAlmostEqual(eb.measurement_total_error(0.03, 0.04), 0.05, places=12)
        self.assertAlmostEqual(eb.checkstar_expected_rms(0.03, 0.04), math.sqrt(0.03**2 + 2 * 0.04**2), places=12)
        self.assertAlmostEqual(
            eb.checkstar_expected_rms(0.03, 0.04, 0.05), math.sqrt(0.03**2 + 0.04**2 + 0.05**2), places=12
        )


class TestCompReader(unittest.TestCase):
    def test_reads_and_converts(self):
        with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False, encoding="utf-8") as fh:
            fh.write("Vmag,e_Vmag,nobs\n11.0,0.09,9\n12.0,0.12,4\n13.0,0.00,20\nbad,x,x\n14.0,0.08,16\n")
            path = fh.name
        try:
            sd = eb.read_comp_csv(path, "Vmag", "e_Vmag")
            self.assertEqual([round(r["err"], 4) for r in sd], [0.09, 0.12, 0.08])
            sem = eb.read_comp_csv(path, "Vmag", "e_Vmag", nobs_col="nobs")
            self.assertEqual([round(r["err"], 5) for r in sem], [0.03, 0.06, 0.02])
        finally:
            os.unlink(path)


class TestFrameBudget(unittest.TestCase):
    def test_frame_total_includes_systematic(self):
        r = eb.frame_budget(mag=12.5, exposure_s=60.0, frames=1, systematic_mag=0.088)
        self.assertGreater(r["per_frame"]["sigma_frame_total_mag"], 0.088)
        self.assertAlmostEqual(
            r["per_frame"]["sigma_frame_total_mag"],
            math.sqrt(r["per_frame"]["sigma_photon_mag"] ** 2 + 0.088**2),
            places=4,
        )

    def test_stacking_two_limiting_cases(self):
        # an offset that is constant over the run does not average down; one that changes between
        # frames does. The stacked photon term must always be photon/sqrt(N).
        r = eb.frame_budget(mag=12.5, exposure_s=60.0, frames=16, systematic_mag=0.088)
        self.assertAlmostEqual(r["stacked"]["sigma_photon_mag"], r["per_frame"]["sigma_photon_mag"] / 4.0, places=4)
        self.assertLess(
            r["stacked"]["if_atm_changes_between_frames_mag"], r["stacked"]["if_atm_constant_over_run_mag"]
        )
        self.assertAlmostEqual(r["stacked"]["if_atm_constant_over_run_mag"], 0.088, places=3)

    def test_dominant_term_flips_faint_to_bright(self):
        bright = eb.frame_budget(mag=9.0, exposure_s=60.0)
        faint = eb.frame_budget(mag=17.0, exposure_s=20.0, sky_mag_per_arcsec2=19.0)
        self.assertEqual(bright["dominant_term"], "source")
        self.assertIn(faint["dominant_term"], ("sky", "readnoise"))


class TestShippedData(unittest.TestCase):
    """The claims the guide makes must be true of the CSV that ships in this repo."""

    @unittest.skipUnless(os.path.exists(DATA_CSV), "shipped APASS DR9 CSV missing")
    def test_apass_errors_are_per_measurement_and_flat_in_magnitude(self):
        sd = eb.read_comp_csv(DATA_CSV, "vmag", "e_vmag")
        sem = eb.read_comp_csv(DATA_CSV, "vmag", "e_vmag", nobs_col="nobs")
        self.assertGreater(len(sd), 2000)
        self.assertTrue(all(r["err"] > 0 for r in sd))
        # per-measurement scatter: ~0.088 mag, and essentially magnitude independent
        med = median([r["err"] for r in sd])
        self.assertLess(abs(med - 0.088), 0.004)
        bright = [r["err"] for r in sd if r["mag"] < 12.0]
        faint = [r["err"] for r in sd if r["mag"] > 14.5]
        self.assertLess(abs(median(bright) - median(faint)), 0.02)
        # catalogue means are ~5x better, which is the ratio sqrt(<nobs>) with nobs ~ 20
        ratio = med / median([r["err"] for r in sem])
        self.assertTrue(4.0 < ratio < 6.5, f"unexpected sd/sem ratio {ratio}")

    def test_reported_error_of_0085_implies_snr_13(self):
        self.assertAlmostEqual(eb.mag_error_to_snr(0.085), 12.8, places=1)

    def test_cli_selftest_runs_offline(self):
        out = subprocess.run(
            [sys.executable, os.path.join(ROOT, "tools", "error_budget.py"), "budget", "--mag", "12.5", "--exposure", "60"],
            capture_output=True,
            text=True,
            check=True,
        )
        doc = json.loads(out.stdout)
        self.assertIn("per_frame", doc)
        self.assertGreater(doc["per_frame"]["snr"], 1.0)


class TestCLI(unittest.TestCase):
    """The documented command lines must work offline and print the numbers in the guide."""

    def run_cli(self, *argv):
        out = subprocess.run(
            [sys.executable, os.path.join(ROOT, "tools", "error_budget.py"), *argv],
            capture_output=True,
            text=True,
            check=True,
        )
        return json.loads(out.stdout)

    def test_stack_identical_frames(self):
        doc = self.run_cli("stack", "--frame-sigma", "0.088", "--n", "16")
        self.assertAlmostEqual(doc["unweighted_all"], 0.022, places=5)
        self.assertAlmostEqual(doc["frame_snr_equivalent"]["unweighted_all"], 49.35, places=2)

    def test_stack_needs_some_input(self):
        proc = subprocess.run(
            [sys.executable, os.path.join(ROOT, "tools", "error_budget.py"), "stack"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(proc.returncode, 2)

    def test_checkstar_uses_the_check_star_catalogue_error(self):
        doc = self.run_cli(
            "checkstar", "--target-sigma", "0.0026", "--ensemble-sigma", "0.0067",
            "--check-cat-sigma", "0.0216", "--measured-rms", "0.085",
        )
        self.assertAlmostEqual(doc["expected_rms_mag"], math.sqrt(0.0026**2 + 0.0067**2 + 0.0216**2), places=4)
        self.assertAlmostEqual(doc["measured_over_expected"], 3.73, places=2)

    def test_ensemble_from_shipped_catalogue(self):
        doc = self.run_cli(
            "ensemble", "--comps", DATA_CSV, "--target-mag", "11.5", "--dm", "0.5", "--n", "8", "--err-mode", "sem",
        )
        self.assertEqual(doc["n_comps"], 8)
        self.assertAlmostEqual(doc["ensemble_error_unweighted_mag"], 0.0067, places=4)


if __name__ == "__main__":
    unittest.main(verbosity=2)
