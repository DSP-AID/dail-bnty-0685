# Photometry error budget: stacking versus measuring, and why your error is 0.085 mag

Companion document for AAVSO HELP thread 3/3 ("Photometry error budget - stack vs measure, and
diagnosing high errors"). Written for CCD observers of any experience level: every formula is
written out, every number is either measured from public data shipped in this repository or a
published instrument test, and `bash verify.sh` reproduces the whole thing offline on a machine
with no API keys.

Subject star: **TV Cas** (AAVSO AUID `000-BBB-443`), V = 7.24 - a bright, variable target in a
sparse field, which turns out to be the worst possible case for comparison-star work.

---

## The short answers

1. **A single ground-based V measurement costs about 0.088 mag, and that cost is almost
   independent of your telescope, your exposure time and the star's brightness.** Measured on
   2 772 real APASS DR9 stars in the TV Cas field: median per-measurement scatter 0.088 mag
   (p25 0.061, p90 0.139), and the same value, 0.086-0.095 mag, in every magnitude bin from
   V = 10 to V = 15.5 - a range over which photon noise alone should differ by more than a
   factor of ten. The atmosphere and the calibration are what you are paying for, not photons.
2. **A quoted error of 0.085 mag means SNR 12.8, nothing more and nothing less**
   (SNR = 1.0857 / sigma). If your aperture SNR is 460 (V = 11.3, 60 s on a 20 cm scope), then
   the 0.085 mag you reported is *not* photon noise: it is the calibration and comparison-star
   term, and no amount of extra exposure or extra frames will remove the part of it that does
   not change between frames.
3. **Stack-then-measure and measure-then-average are the same arithmetic** for identical frames:
   sigma / sqrt(N). 16 frames at 0.088 mag each give 0.022 mag, i.e. SNR 49. What you actually
   report is the implied SNR of the combination, never the per-frame error.
4. **Frames of unequal quality are where "measure each frame" wins**: on a realistic run of four
   good frames plus one clouded frame (frame errors taken from the real APASS percentiles), an
   unweighted mean over all five costs 0.0887 mag against 0.0379 mag if the bad frame is dropped
   - a **+134 %** penalty - while inverse-variance weighting on all five gives 0.0353 mag, better
   than dropping it by hand.
5. **The comparison-star ensemble sets a floor you cannot buy your way out of.** With 8
   magnitude-matched stars whose catalogue magnitudes are good to their catalogue-mean error, the
   zero point is good to 0.0067 mag; if those adopted magnitudes are only as good as a single
   ground-based measurement (APASS-grade), the same 8 stars give 0.0322 mag, and your check star
   then wanders by about 0.094 mag - exactly the "0.08-0.09" that observers report as their
   personal error.
6. **Bright targets are a special failure mode.** In a 0.5 degree APASS DR9 cone on TV Cas there
   is not one star with a usable quoted error within 2 mag of V = 7.24; the brightest usable
   comparison is V = 10.011, and any ensemble you build sits 2.9-3.5 mag below the target. A
   bright star needs its own error budget; the field's typical value does not apply to it.

---

## 1. What limits one ground-based measurement (measured, not modelled)

The catalogue column `e_vmag` in APASS DR9 is the scatter of the individual V measurements that
were averaged into the catalogue magnitude. Section 1 of
[`results/worked_examples.md`](../results/worked_examples.md) computes it per magnitude bin from
the shipped CSV:

| V bin | stars | median `nobs` | median per-measurement scatter | error of the catalogue mean |
|---|---|---|---|---|
| 10-11 | 24 | 23.5 | 0.0870 mag | 0.0180 mag |
| 11-12 | 57 | 22.0 | 0.0860 mag | 0.0181 mag |
| 12-13 | 141 | 22.0 | 0.0890 mag | 0.0185 mag |
| 13-14 | 330 | 22.0 | 0.0880 mag | 0.0186 mag |
| 14-15 | 552 | 20.0 | 0.0950 mag | 0.0216 mag |
| 15-16 | 1027 | 11.0 | 0.0780 mag | 0.0242 mag |

Two things follow from this table, and both are tested rather than asserted:

* **The scatter is flat in magnitude.** V = 10.5 and V = 15.5 give the same ~0.09 mag even though
  the photon rate differs by a factor of about 100. A photon-noise-limited measurement cannot do
  that. Whatever dominates does not care about signal level: it is transparency change,
  scintillation, colour terms and plate/zero-point calibration.
* **The scatter does not shrink with the number of observations.** Within one magnitude bin the
  log-log slope of `e_vmag` against `nobs` is +0.07, -0.26, +0.04, +0.09, +0.15, +0.02 (R^2 at
  most 0.30, mostly below 0.03). A standard error of the mean would give a slope of -0.5. This
  matters because APASS does not document which one it is; the test settles it. `e_vmag` is a
  per-measurement scatter, and the error on an *adopted* catalogue magnitude is
  `e_vmag / sqrt(nobs)`, the last column above.

The same field measured from space is the control experiment. Gaia DR3's own quoted flux errors
(`sigma_mag = 1.0857 / (flux / flux_error)`) work out to 0.00018-0.00028 mag for G = 10-16, so at
G ~ 13.5 the ground-based per-measurement scatter is larger by a factor of **476**. Photon
statistics are not what limits ground-based photometry; the atmosphere plus the calibration is.

For the rest of this document that measured 0.088 mag is used as one lumped "atmosphere +
calibration" systematic, which is deliberate: it is the quantity an observer can check against
their own data, and it is what the catalogue says one measurement costs. Krisciunas & Schaefer
(1991) would split it into scintillation, transparency and extinction-coefficient terms, and that
split is where this analysis stops short - see section 7.

---

## 2. The arithmetic you need, in four lines

| what you have | what to report |
|---|---|
| aperture SNR of one measurement | sigma_mag = **1.0857 / SNR** mag |
| a reported sigma, wanting the honest SNR | SNR = **1.0857 / sigma** |
| N independent frames of equal quality | sigma_mean = **sigma_frame / sqrt(N)** |
| ensemble zero point from N comparison stars with errors e_i | **sqrt(sum(e_i^2)) / N** unweighted, **1 / sqrt(sum(1/e_i^2))** inverse-variance weighted |

The second line is the reality check that answers most HELP posts:

| reported sigma (mag) | SNR it actually implies |
|---|---|
| 0.005 | 217 |
| 0.010 | 109 |
| 0.020 | 54 |
| 0.050 | 22 |
| 0.085 | **12.8** |
| 0.100 | 11 |

Variance adds in quadrature, never linearly: a measurement against an ensemble is
`sqrt(sigma_frame^2 + sigma_ensemble^2)`, and a check star deviates by
`sqrt(sigma_own^2 + sigma_ensemble^2 + sigma_check_star_catalogue^2)` (the ensemble term does not
cancel in that difference - it enters with the opposite sign).

The calculator, `tools/error_budget.py`, implements all of this plus the full aperture model
(Howell 1989; sky-annulus term after Merline & Howell 1995). These are the commands used to
produce every number in this document:

```bash
python3 tools/error_budget.py budget --mag 12.5 --exposure 60 --frames 20 --atm 0.088
    # full budget, variance shares, and both stacking limits (0.0197 mag vs 0.0880 mag)

python3 tools/error_budget.py sens --mag 15 --exposure 60
    # which parameter would actually help

python3 tools/error_budget.py stack --frame-sigma 0.088 --n 16
    # both routes: 0.0220 mag, implied SNR 49

python3 tools/error_budget.py ensemble --comps data/apass_dr9_tvcas_r0p5.csv \
    --target-mag 11.5 --dm 0.5 --n 8 --err-mode sem
    # ensemble floor from a real catalogue: 0.0067 mag (drop --err-mode for the 0.0322 mag case)

python3 tools/error_budget.py checkstar --target-sigma 0.0026 --ensemble-sigma 0.0067 \
    --check-cat-sigma 0.0216 --measured-rms 0.085
    # expected 0.0228 mag, measured/expected = 3.73
```

---

## 3. Stack then measure, or measure then average?

For identical frames the two routes are algebraically identical: summing N frames multiplies the
signal and adds the variances, so SNR grows as sqrt(N), which is the same thing as averaging N
magnitude measurements each good to sigma_frame. The measured APASS percentiles make it concrete
(frame error 0.088 mag, the real median):

| frames N | stack then measure | measure then average | identical? | implied SNR |
|---|---|---|---|---|
| 1 | 0.0880 | 0.0880 | yes | 12.3 |
| 4 | 0.0440 | 0.0440 | yes | 24.7 |
| 9 | 0.0293 | 0.0293 | yes | 37.0 |
| 16 | 0.0220 | 0.0220 | yes | 49.4 |
| 25 | 0.0176 | 0.0176 | yes | 61.7 |
| 36 | 0.0147 | 0.0147 | yes | 74.0 |

So the choice is not "which is more accurate" - they are the same - but *which one exposes your
bad frames*. Frames of unequal quality, using real percentiles from the same catalogue (p25 0.061,
median 0.088, and one clouded frame at 3 x p90 = 0.417 mag):

| method | error on the mean |
|---|---|
| unweighted, all five kept | 0.0887 mag |
| unweighted, clouded frame dropped | 0.0379 mag |
| inverse-variance weighted, all five kept | **0.0353 mag** |

One bad frame in an unweighted mean costs **+134 %**; weighting absorbs it automatically and even
beats manual rejection. That is the practical case for measuring every frame and combining with
inverse-variance weights (or, equivalently, weighting the stack) instead of taking one unweighted
mean. Keep the summed stack when you need maximum detection depth on a faint object; report the
weighted mean of per-frame magnitudes when you need an *error bar*, because then each frame's
quality is visible and rejectable.

One warning that follows from section 1: if your per-frame error is dominated by the 0.088 mag
atmosphere+calibration term, stacking N frames only helps if that term changes between frames.
The calculator prints both limits; for V = 12.5, 60 s, 20 frames it gives **0.0197 mag** if the
systematic varies frame to frame (it averages down) and **0.0880 mag** if it is a constant offset
over the run (it does not). Re-deriving the zero point on every frame - which is what
measure-then-average with a fresh ensemble does - is the only way to attack the varying case.

---

## 4. "My SNR says 0.01, why is my error 0.085?"

Reference setup used throughout (camera numbers from the published AAVSO ASI2600MM Pro test; all
other inputs are labelled assumptions swept in section 5 of the worked examples): 20 cm f/10
reflector, 3.76 um pixels, gain 0 (0.77 e-/ADU), read noise 2.80 e-, dark 0.1 e-/s, sky
20.7 mag/arcsec^2, aperture radius 5 px with a 10-16 px annulus.

| V | exposure | star electrons | photon SNR | sigma from photons alone | one frame incl. the measured floor |
|---|---|---|---|---|---|
| 11.3 | 60 s | 213 740 | 460 | 0.0024 mag | 0.0880 mag |
| 12.5 | 60 s | 70 776 | 263 | 0.0041 mag | 0.0881 mag |
| 15.0 | 60 s | 7 078 | 75 | 0.0144 mag | 0.0892 mag |

The first row is the whole mystery: photon statistics say 0.0024 mag, the measurement costs
0.088 mag, a factor of 37. The reason is section 1 - that 0.088 mag is what a real ground-based V
measurement costs at *any* magnitude, measured on 2 772 stars. Your SNR really is 460; your error
is not photon noise.

At the faint end the detector terms become visible again, and that is where effort pays:

| V = 15.0 | dominant term | source | sky | dark | read noise | sky determination |
|---|---|---|---|---|---|---|
| 30 s | source | 0.740 | 0.046 | 0.049 | **0.129** | 0.036 |
| 60 s | source | 0.800 | 0.050 | 0.053 | 0.070 | 0.028 |
| 120 s | source | 0.834 | 0.052 | 0.056 | 0.036 | 0.023 |

Sensitivity of the photon error at V = 15.0, 60 s (full table in section 5 of the worked
examples): exposure 15 -> 240 s moves sigma from 0.0322 to 0.0070 mag; aperture 8 -> 30 cm moves
it from 0.0477 to 0.0092 mag; sky brightness from 19.5 to 21.8 mag/arcsec^2 only moves it from
0.0152 to 0.0142 mag; read noise from 1 to 10 e- moves it from 0.0139 to 0.0202 mag and makes
read noise the dominant term. So at V = 15 with this camera, exposure time and aperture are worth
money and moonlight is nearly free - while at V = 12.5 the star itself carries 97.6 % of the
variance and none of those knobs do anything at all.

**One subtlety that trips people up.** The 0.088 mag floor is the cost of an *absolute*
measurement - one frame calibrated onto catalogue magnitudes. When you measure a target against
ensemble stars on the same image, whatever is common to target and comparison stars cancels: that
is why differential photometry can beat the absolute floor. How much survives is site-specific
(scintillation decorrelates with separation and with aperture, transparency does not, flat-field
and colour terms depend on your setup), and the quantity that tells you is the check star -
section 6. Stacking frames cannot remove whatever survives; it only averages down the part that is
independent per frame.

---

## 5. Comparison-star ensembles: the floor you cannot buy out of

Two cases from the shipped catalogue, both real stars. Case A is TV Cas itself (V = 7.239, 16
measurements; APASS quotes `e_vmag = 0` with the error flagged unreliable, so the catalogue gives
no usable error for this star at all). Gaia confirms it as variable (`phot_variable_flag =
VARIABLE`, G = 7.203, RUWE 0.826) and quotes a formal sigma_G of 0.0017 mag - still hundreds of
times smaller than the ground-based floor. Case B is a real field star at V = 11.50, which is
where most AAVSO work sits.

How many stars in the same 0.5 degree APASS cone can serve as comparison stars:

| comps with |dV| <= | case A, target V = 7.24 | case B, target V = 11.50 |
|---|---|---|
| 0.25 mag | 0 | 27 |
| 0.50 mag | 0 | 57 |
| 1.00 mag | 0 | 131 |
| 2.00 mag | 0 | 348 |
| 3.00 mag | 4 | 797 |

For case A the *brightest star in the whole field with a usable error* is V = 10.011, i.e. 2.77
mag fainter than the target, so every ensemble is forced to sit 2.9-3.5 mag away from it. Across
that much of a magnitude range you carry colour difference between target and comps, sensor
nonlinearity close to saturation on the target while the comps are photon-starved, and flat-field
error that scales with the signal ratio. A mag 7 target needs its own error budget; the field's
typical value does not describe it.

Ensemble zero-point error versus number of stars used (all values from the catalogue):

| N used | case A comp V range | case A ensemble (per-measurement basis) | case B comp V range | case B ensemble (per-measurement) | case B ensemble (catalogue-mean basis) |
|---|---|---|---|---|---|
| 3 | 10.01..10.17 | 0.0761 mag | 11.49..11.50 | 0.0483 mag | 0.0103 mag |
| 8 | 10.01..10.55 | 0.0393 mag | 11.42..11.58 | 0.0322 mag | 0.0067 mag |
| 20 | 10.01..10.69 (only 15 exist) | 0.0261 mag | 11.38..11.65 | 0.0197 mag | 0.0041 mag |

Read the two right-hand columns of case B as your realistic range: if you adopt comparison
magnitudes at their catalogue-mean quality, eight matched stars give a zero point good to
**0.0067 mag**; if the magnitudes you use are only as good as one ground-based measurement each,
the same eight stars give **0.0322 mag**. Five times better for doing the bookkeeping properly -
that is what "use a properly calibrated ensemble" is worth in practice.

Magnitude-matched versus faint comps (8 stars each, case B): matched comps (median V = 11.49)
give 0.0322 / 0.0067 mag; comps 1.5-2.5 mag fainter (median V = 13.03) give 0.0648 / 0.0415 mag,
twice as bad, because their own measurements were noisier.

And the other half of the trade - what *you* pay in photon noise to measure a comparison star of a
given magnitude yourself, 60 s on the reference setup:

| comp V | star electrons | SNR | sigma_photon | incl. 0.088 mag floor |
|---|---|---|---|---|
| 8.50 | 2 817 636 | 1 678 | 0.0006 | 0.0880 |
| 11.50 | 177 781 | 420 | 0.0026 | 0.0880 |
| 13.50 | 28 176 | 163 | 0.0067 | 0.0882 |
| 15.50 | 4 466 | 57 | 0.0192 | 0.0901 |

Going faint for comparison stars is cheap in photon terms down to about V = 14 and only then
starts to bite; what it really costs you is the colour match, the saturation behaviour on the
bright target, and a worse catalogue error per comp. That is the whole trade-off, and both sides
of it are now numbers.

---

## 6. Diagnosing a high error: what your check star should do

A check star is measured with the ensemble zero point and compared against its own adopted
magnitude, so three terms add in quadrature: its measurement noise, the ensemble error, and its
own catalogue error. With the real catalogue numbers from section 1 (median per-measurement
scatter 0.088 mag, median catalogue-mean error 0.0216 mag) and the calculator's photon term for a
V = 11.50 star in 60 s:

| regime | target sigma | ensemble sigma | check-star catalogue sigma | expected RMS | implied SNR |
|---|---|---|---|---|---|
| your own frames only, no catalogue involved | 0.0026 | - | - | **0.0026** | 419 |
| check star vs catalogue, comps at their catalogue-mean error | 0.0026 | 0.0067 | 0.0216 | **0.0228** | 48 |
| check star vs catalogue, comps as good as one measurement each | 0.0026 | 0.0322 | 0.0880 | **0.0937** | 12 |
| absolute calibration of a single frame (what APASS itself reports) | 0.0880 | 0.0067 | 0.0216 | **0.0909** | 12 |

Row 1 is what your camera can do on one frame; row 3 is what you get when the comparison-star
magnitudes you adopted are only APASS-per-measurement quality. If your check star RMS is in the
0.08-0.09 range and your ensemble is built from single-measurement-grade catalogue magnitudes,
**your error budget is already correct** - the fault is in the comparison-star calibration, not in
your aperture, your SNR or your focusing.

Reading a measured value against the best realistic case (0.0228 mag):

| measured check-star RMS | ratio to expectation | systematic left over | what it means |
|---|---|---|---|
| 0.010 | 0.44 | - | better than expected: your comps are correlated with each other, or your quoted errors are too big |
| 0.020 | 0.88 | - | explained by the comparison-star catalogue; your camera is not the limit |
| 0.030 | 1.32 | 0.0195 mag | still within the catalogue floor; do not chase it with exposure time |
| 0.050 | 2.20 | 0.0445 mag | worse than the budget: suspect a variable or mis-catalogued comp, a wrong sky annulus, or a colour term |
| 0.085 | 3.73 | 0.0819 mag | as large as the measured 0.088 mag calibration floor: per-frame systematics, not photons |
| 0.110 | 4.83 | 0.1076 mag | something specific is wrong (variable comp, blend, saturation); start with per-comp residual plots |

### The checklist, in order

1. Compute `sigma_photon` for your setup and exposure (`error_budget.py budget`). If it is much
   smaller than your reported error, stop adding exposure time: you are not photon limited.
2. Compute the ensemble floor from your comparison stars' *adopted* magnitude errors
   (`error_budget.py ensemble --err-mode sem` for catalogue-mean quality). That is your hard floor.
3. Expected check-star RMS = `sqrt(sigma_own^2 + sigma_ensemble^2 + sigma_check_catalogue^2)`.
4. Measure the actual check-star RMS over your run and take the ratio. Below ~1.5: the budget
   explains you; above it, go hunting.
5. When the ratio is high, in this order: plot each comparison star's residual against time (a
   variable or mis-catalogued comp shows up as one line drifting or jumping); re-check the sky
   annulus against a neighbouring star and against the seeing; check the target's peak ADU against
   the linear range (bright targets saturate while comps do not); check the colour difference
   between target and mean comp.
6. Bright-target special case (TV Cas): with no comp within 2 mag, quote the larger error rather
   than the field average, use short exposures or defocus/ND to keep the target unsaturated, and
   prefer ensembles calibrated on Gaia/Sloan-based photometry so that the comps' catalogue errors
   are not the dominant term.

---

## 7. What this analysis does not cover

* The measured 0.088 mag floor is used as one term. Scintillation, transparency change and
  extinction-coefficient error are not separated; Krisciunas & Schaefer (1991) give the formulas
  if you want the split, and their scintillation term scales as D^(-2/3) t^(-1/2) X^(3/2), which
  is why small telescopes hit this floor early.
* The differential regime is not modelled term by term: how much of the floor cancels between
  target and comps depends on separation, aperture and weather, so section 4 tells you to measure
  that residual with a check star rather than guess it.
* Colour terms, crowding/blending and nonlinearity are discussed qualitatively only; the numbers
  in section 5 for bright targets are about comparison-star availability and catalogue errors, not
  a colour-term model.
* Sky brightness (20.7 mag/arcsec^2), dark current, throughput, central obstruction and the Vega
  photon flux used to set the electron scale are labelled assumptions, swept over in section 5 of
  the worked examples so that no conclusion depends on them. Camera gain and read noise come from
  one published AAVSO instrument test, not from your camera - re-measure them at home.
* No AAVSO light curves were used: the AID API needs a token, so the "measured check-star RMS"
  values in section 6 are labelled probe values spanning the range observers report, not shipped
  data. Everything else is shipped data.

---

## 8. Repository layout and how to reproduce it

```
data/     apass_dr9_tvcas_r0p5.csv, gaia_dr3_tvcas_r0p5.csv, vsp_target_meta.json, MANIFEST.json
tools/    fetch_data.py (v1.0.0, the only script that needs a network),
          error_budget.py (v1.0.0, calculator + `selftest`), run_examples.py (v1.0.0, generator)
results/  worked_examples.md (every table quoted above), worked_examples.json (same, machine readable)
docs/     sources.md (provenance, exact queries, column semantics, measured-versus-assumed table)
tests/    test_error_budget.py (36 unit tests, including checks that the shipped CSV really
          supports the 0.088 mag flat-floor claim, and that the documented commands work)
guide/    this document
verify.sh one command: selftest + unit tests + regenerate-and-diff + sha256 against MANIFEST.json
```

```bash
bash verify.sh                            # everything, offline, stdlib Python 3 only
python3 tools/error_budget.py selftest    # 23 formula checks
python3 -m unittest discover -s tests     # 36 tests
python3 tools/run_examples.py --check     # regenerates and byte-compares results/
```

Versions: `error_budget.py` 1.0.0, `run_examples.py` 1.0.0, `fetch_data.py` 1.0.0; data retrieved
2026-11-18T23:33:55Z (hashes in `data/MANIFEST.json`). See `docs/sources.md` for the exact ADQL
queries, the catalogue column definitions relied on, and the papers behind each formula.
