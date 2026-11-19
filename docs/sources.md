# Sources, provenance and what is measured versus assumed

Every number in `guide/aavso-error-budget-guide.md` comes from one of three places: a public
catalogue shipped in `data/`, a published instrument test or paper cited here, or an explicitly
labelled assumption that section 5 of the worked examples sweeps. Nothing else was used.

## 1. Datasets shipped in this repository

They are frozen copies, so the whole analysis runs offline with no API keys or tokens. The fetcher
(`tools/fetch_data.py`, v1.0.0) is only needed to re-download them; the analysis scripts never
touch the network. All three were retrieved on **2026-11-18T23:33:55Z**.

| file | what it is | source | rows | sha256 (full value in `data/MANIFEST.json`) |
|---|---|---|---|---|
| `data/vsp_target_meta.json` | target identity for the star that the AAVSO HELP thread is about, TV Cas (AUID `000-BBB-443`, chartid `X43374GGF`, J2000 00:19:18.74 +59:08:20.6) | AAVSO VSP chart API, `https://app.aavso.org/vsp/api/chart?fov=30&star=TV+Cas&maglimit=14.5&format=json` | - | `5dc2ebe9...e7d8d8f` |
| `data/apass_dr9_tvcas_r0p5.csv` | AAVSO APASS DR9 all-sky photometry, 0.5 deg cone on the target | ESA Gaia archive TAP sync endpoint `https://gea.esac.esa.int/tap-server/tap/sync`, schema `external.apassdr9` (the archive's mirror of VizieR II/336/apass9) | 3262 | `19e7c139...cf1b271b` |
| `data/gaia_dr3_tvcas_r0p5.csv` | Gaia DR3 sources in the same cone, brightest first (`TOP 2000`, G < 19) | same TAP endpoint, schema `gaiadr3.gaia_source` | 2000 | `f32ff04d...845a8708` |

`data/MANIFEST.json` holds the exact ADQL text of both queries plus the hashes, so a reader can
re-run the identical query and diff the result. `verify.sh` re-checks these hashes on every run.

## 2. Column semantics that the analysis depends on

### APASS DR9 (`external.apassdr9`, VizieR II/336/apass9)

| column | published definition used here |
|---|---|
| `vmag` | Johnson V-band magnitude of the catalogue mean |
| `e_vmag` | "Vmag uncertainty" / "Uncertainty in Vmag" (VizieR II/336 and its mirrors) |
| `u_e_vmag` | "[0/1] Uncertainty flag on e_Vmag" - 1 means the catalogue itself does not trust that error |
| `nobs` | number of contributing observations behind the mean magnitude |

The published definition does **not** say whether `e_vmag` is the standard deviation of the
individual measurements or the standard error of their mean, and the difference matters a great
deal. It is settled empirically in section 1 of `results/worked_examples.md`: within a single
magnitude bin the log-log slope of `e_vmag` against `nobs` is +0.02 to +0.15 with R^2 at most
0.30, whereas a standard error of the mean would give -0.5. So `e_vmag` is treated as a
**per-measurement** scatter, and divided by `sqrt(nobs)` whenever the error of an *adopted*
catalogue magnitude is needed.

Rows with `u_e_vmag = 1` (480 of 3262 here, 456 of which carry `e_vmag = 0`) are dropped from every
statistic; that leaves 2772 usable stars.

### Gaia DR3 (`gaiadr3.gaia_source`)

`phot_g_mean_flux` and `phot_g_mean_flux_error` are converted with
sigma_mag = 2.5 / ln(10) * (sigma_F / F) = 1.0857 / (F / sigma_F), the standard large-SNR
conversion used by the Gaia photometry pipeline documentation
(`gea.esac.esa.int/archive/documentation/GDR3/`, release notes, heliocentric photometry section).
`phot_variable_flag` and `ruwe` are used as published.

## 3. Method references

| reference | what it provides | where it is used |
|---|---|---|
| Howell, S. B. 1989, PASP 101, 616 (aperture photometry SNR) | the aperture SNR equation with source, sky, dark and read-noise terms in electrons | `tools/error_budget.py: aperture_snr()` |
| Merline, J. J. & Howell, S. B. 1995, PASP 107, 544 | the extra variance from estimating the sky from a finite annulus, (n_ap^2 / n_ann) * var_pixel | same function, `annulus_pixels` branch; reported as the `sky_determination` variance share |
| Krisciunas, K. & Schaefer, B. E. 1991, PASP 103, 1033 | scintillation noise scaling as D^(-2/3) * t^(-1/2) * X^(3/2), and the fact that it is only partly cancelled by differential photometry | used qualitatively in the guide; **not** implemented as a separate term (see the table below) |
| Bessell, M. S. 1990, PASP 102, 1181 (Vega spectrophotometric fluxes) | Vega photon flux density near the V effective wavelength, about 970 photons cm^-2 s^-1 A^-1, and the Johnson V width of 880 A | `tools/error_budget.py: zeropoint_from_flux()` defaults |
| AAVSO instrument test report, ZWO ASI2600MM Pro, `https://archive.aavso.org/zwo-asi2600mm-pro-review` | gain 0 / offset 10 gives 0.77 e-/ADU with 2.80 e- read noise; 2x2 binning gives 3.08 e-/ADU with 6.63 e- | the reference setup in `tools/run_examples.py` (`SETUP`) and section 4 of the guide |
| AAVSO VSP chart API | target position, AUID and chart id used to identify the star in the catalogues | `data/vsp_target_meta.json`, `tools/fetch_data.py` |

## 4. Measured, published, or assumed

| quantity | value used | status |
|---|---|---|
| per-measurement V scatter floor of ground-based photometry | **0.088 mag** (p25 0.061, p90 0.139) | **measured** from the shipped APASS CSV, median of 2772 rows |
| same, as a function of magnitude | flat from V = 10 to V = 15.5 (0.086-0.095 mag per bin) | **measured**, section 1 table |
| error of an adopted catalogue magnitude | 0.018-0.024 mag | **measured** (`e_vmag / sqrt(nobs)`) |
| Gaia DR3 formal sigma_G at G = 10-16 | 0.00018-0.00028 mag | **measured** from the shipped Gaia CSV |
| TV Cas brightness and comparison-star availability | V = 7.239 (APASS), G = 7.203 (Gaia, `phot_variable_flag = VARIABLE`); brightest usable comp in the field V = 10.011 | **measured** |
| ensemble calibration errors in guide section 5 | computed from the rows above | **derived**, no new input |
| camera gain and read noise | 0.77 e-/ADU, 2.80 e- | **published third-party test** (AAVSO report), not reproduced here |
| sky brightness | 20.7 mag/arcsec^2 in V | **assumption**; swept over 19.5-21.8 in section 5 of the worked examples |
| dark current | 0.1 e-/s | **assumption**; swept over 0-1.0 |
| total throughput, central obstruction | 0.45, and 0.15 linear | **assumption**; only sets the absolute electron scale, and section 5 sweeps aperture diameter 8-30 cm instead |
| Vega photon flux, V bandwidth | 970 ph cm^-2 s^-1 A^-1, 880 A | **published** (Bessell 1990); affects only the absolute electron scale - sigma_photon goes as sqrt(zeropoint) |
| atmosphere + calibration term used as one lumped 0.088 mag systematic | measured size, single term | **modelled coarsely on purpose**: Krisciunas & Schaefer would split it into scintillation, transparency and extinction-coefficient terms. The guide says so: the size of the term is measured, its internal split is not |
| "measured" check-star RMS values in guide section 6 (0.01-0.11 mag) | probe values | **illustrative inputs**, clearly labelled: the range an observer might measure, not shipped data |

## 5. What was tried and rejected (so that nobody repeats it)

| source | outcome |
|---|---|
| AAVSO AID photometry API | HTTP 401 without a token, so it cannot be part of an offline verification path |
| `www.aavso.org` help/forum HTML pages | HTTP 403 from the WAF for non-browser clients |
| VSX object detail pages | JavaScript shell, no data in the response body |
| Gaia `epoch_photometry` | not populated for this field in the current archive release (empty result) |
| ASAS-SN light-curve endpoints | connection failed from this network at retrieval time |
| MAST public APSS download | HTTP 403 without credentials |

The two keyless services that did work - ESA Gaia TAP sync and the AAVSO VSP chart API - are the
ones shipped in `data/`.

## 6. Reproducing everything from a clean checkout

```bash
bash verify.sh                                  # all of the below + sha256 check against MANIFEST
python3 tools/error_budget.py selftest          # 23 formula checks
python3 -m unittest discover -s tests           # 36 unit tests, incl. shipped-data and CLI checks
python3 tools/run_examples.py --write           # regenerates results/worked_examples.{md,json}
python3 tools/run_examples.py --check           # fails if the committed results differ
```

`tools/fetch_data.py` is the only script that needs a network, and it writes a fresh
`MANIFEST.json`; nothing in `verify.sh` calls it.
