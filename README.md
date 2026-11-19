# AAVSO HELP 3/3 - photometry error budget: stack vs measure, and diagnosing high errors

An error-budget guide for AAVSO CCD observers, built on real public data instead of rule-of-thumb
numbers, with a calculator you can run against your own comparison-star list.

Start with **[`guide/aavso-error-budget-guide.md`](guide/aavso-error-budget-guide.md)**. It answers
the three questions the HELP thread asks:

| question | where the answer is |
|---|---|
| where do photometric errors actually come from | guide sections 1-2, measured on 2 772 APASS DR9 stars and cross-checked against Gaia DR3 in the same field |
| stack-then-measure versus measure-then-average | guide section 3 (they are algebraically identical for equal frames; weighting is what matters for unequal ones) |
| why a reported error is 0.085 mag when the SNR looks like 0.01, and how to diagnose it | guide sections 4 and 6, with a check-star decision table and a calculator command for each step |

## Verify it in one command

```bash
bash verify.sh
```

Runs offline on a clean machine with no API keys: 23 formula checks (`error_budget.py selftest`),
36 unit tests, sha256 of every shipped data file against `data/MANIFEST.json`, and a byte-for-byte
regeneration of `results/worked_examples.md`. Standard library Python 3 only - no numpy, no
network.

## Headline results (all measured, all reproducible)

* One ground-based V measurement costs **0.088 mag**, flat from V = 10 to V = 15.5: median
  per-measurement scatter of APASS DR9 in a 0.5 deg cone on TV Cas. The log-log slope of that
  scatter against `nobs` is ~0, which proves the column is a per-measurement standard deviation
  and not a standard error - something the catalogue itself does not state.
* The same field from space: Gaia DR3 formal sigma_G = 0.00018-0.00028 mag, i.e. ground-based
  photometry is worse by a factor of **476** at G ~ 13.5. Photon statistics are not the limit.
* A reported 0.085 mag error is SNR 12.8; on the reference setup V = 11.3 gives photon SNR 460 and
  sigma_photon = 0.0024 mag, so the extra factor of 37 is atmosphere + calibration + comparison
  stars, not light grasps or exposure time.
* Stacking N equal frames and averaging N measurements both give sigma/sqrt(N) (16 x 0.088 ->
  0.022 mag, SNR 49). For unequal frames, inverse-variance weighting on all five frames of a
  clouded run gives 0.0353 mag where an unweighted mean gives 0.0887 mag (**+134 %** penalty).
* Ensemble floors from the real catalogue: 8 magnitude-matched comps at their catalogue-mean error
  give **0.0067 mag**; the same 8 comps at per-measurement quality give **0.0322 mag**, and a check
  star then wanders by **0.094 mag** - which is exactly the "0.08-0.09" observers report.
* The target itself, TV Cas (V = 7.239), has **zero** comparison stars with a usable quoted error
  within 2 mag in that whole field; the brightest usable one is V = 10.011.

## Layout

```
guide/aavso-error-budget-guide.md   the deliverable: explanation, tables, decision checklist
results/worked_examples.md          every table quoted in the guide, generated from the shipped data
results/worked_examples.json        the same numbers, machine readable
docs/sources.md                     provenance, exact ADQL queries, column semantics,
                                    measured-versus-published-versus-assumed table, dead ends tried
data/                               APASS DR9 + Gaia DR3 cones on TV Cas + VSP target metadata,
                                    with MANIFEST.json holding the queries and sha256 hashes
tools/error_budget.py               calculator: SNR <-> mag error, aperture model, stacking both
                                    routes, ensemble floors, check-star diagnosis (`selftest` built in)
tools/run_examples.py               regenerates results/ from data/ and fails if they drift (--check)
tools/fetch_data.py                 the only script that needs a network (re-download the datasets)
tests/test_error_budget.py          36 tests: formulas, the documented CLI commands, and the claim
                                    that the shipped CSV really shows a flat 0.088 mag floor
verify.sh                           all of the above in one command
```

## Try the calculator on your own data

```bash
python3 tools/error_budget.py budget --mag 12.5 --exposure 60 --frames 20 --atm 0.088
python3 tools/error_budget.py ensemble --comps data/apass_dr9_tvcas_r0p5.csv \
    --target-mag 11.5 --dm 0.5 --n 8 --err-mode sem
python3 tools/error_budget.py checkstar --target-sigma 0.0026 --ensemble-sigma 0.0067 \
    --check-cat-sigma 0.0216 --measured-rms 0.085
```

`--comps` takes any CSV with a magnitude column and an error column, so pointing it at your own
comparison-star list is the intended use (`--help` on any subcommand lists the options).

## Provenance

Datasets were retrieved on 2026-11-18T23:33:55Z from the ESA Gaia archive TAP service (schemas
`external.apassdr9` and `gaiadr3.gaia_source`) and from the AAVSO VSP chart API; the exact queries,
row counts and sha256 hashes are in `data/MANIFEST.json`, and `verify.sh` re-checks them. Camera
numbers come from the published AAVSO ZWO ASI2600MM Pro test report; formulas from Howell (1989),
Merline & Howell (1995) and Krisciunas & Schaefer (1991). `docs/sources.md` records which numbers
are measured, which are published third-party values and which are labelled assumptions that the
sweeps show no conclusion depends on.

Versions: `error_budget.py` 1.0.0, `run_examples.py` 1.0.0, `fetch_data.py` 1.0.0.
