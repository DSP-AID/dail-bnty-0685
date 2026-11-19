#!/usr/bin/env python3
"""Fetch the real public datasets used by the error-budget guide. Version: 1.0.0.

Datasets (all keyless, all cited in docs/sources.md):
  A. AAVSO VSP chart metadata for the target star -> coordinates + chart id + AUID.
     https://app.aavso.org/vsp/api/chart
  B. APASS DR9 (AAVSO Photometric All-Sky Survey) rows in a cone around the target,
     via the Gaia archive TAP mirror table `external.apassdr9`.
     Columns of interest: vmag, e_vmag (= standard deviation of the individual V
     measurements that went into the mean), nobs (number of observations that went
     into vmag), u_e_vmag (error-estimate reliability flag).
  C. Gaia DR3 `gaiadr3.gaia_source` rows in the same cone: G magnitudes with
     phot_g_mean_flux / phot_g_mean_flux_error (fluxes in electron/s) and
     phot_variable_flag.

Re-running is optional: the CSVs are shipped in data/ so every example can be
reproduced offline. This script needs network access; the verification scripts do not.

Usage:  python tools/fetch_data.py [--outdir data] [--radius 0.5]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone

TAP_URL = "https://gea.esac.esa.int/tap-server/tap/sync"
VSP_URL = "https://app.aavso.org/vsp/api/chart"
USER_AGENT = "dail-bnty-0685-error-budget/1.0 (+https://github.com/DSP-AID)"

TARGET_NAME = "TV Cas"  # real AAVSO target, coordinates taken from VSP (not hardcoded)
CONE_RADIUS_DEG = 0.5
GAIA_MAX_ROWS = 2000


def http_get(url: str, timeout: int = 180) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "*/*"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def tap_query(query: str, timeout: int = 300) -> str:
    """Run an ADQL query against the Gaia archive TAP sync endpoint, return CSV text."""
    body = urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": " ".join(query.split())}
    ).encode()
    req = urllib.request.Request(
        TAP_URL,
        data=body,
        headers={"User-Agent": USER_AGENT, "Content-Type": "application/x-www-form-urlencoded"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", "replace")


def hms_to_deg(text: str) -> float:
    parts = [float(x) for x in text.split(":")]
    return 15.0 * (parts[0] + parts[1] / 60.0 + parts[2] / 3600.0)


def dms_to_deg(text: str) -> float:
    sign = -1.0 if text.strip().startswith("-") else 1.0
    parts = [float(x) for x in text.replace("-", "").split(":")]
    return sign * (parts[0] + parts[1] / 60.0 + parts[2] / 3600.0)


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: str, text: str) -> int:
    rows = [ln for ln in text.splitlines() if ln.strip()]
    # validate that it parses as CSV before shipping it
    list(csv.DictReader(rows))
    with open(path, "w", newline="", encoding="utf-8") as fh:
        fh.write("\n".join(rows) + "\n")
    return len(rows) - 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__), "..", "data"))
    ap.add_argument("--radius", type=float, default=CONE_RADIUS_DEG)
    args = ap.parse_args()
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    retrieved = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    chart = json.loads(http_get(f"{VSP_URL}?fov=30&star={urllib.parse.quote(TARGET_NAME)}&maglimit=14.5&format=json").decode())
    ra, dec = hms_to_deg(chart["ra"]), dms_to_deg(chart["dec"])
    chart_meta = {
        "source_url": f"{VSP_URL}?fov=30&star={TARGET_NAME}&maglimit=14.5&format=json",
        "retrieved_utc": retrieved,
        "star": chart["star"],
        "auid": chart["auid"],
        "chartid": chart["chartid"],
        "ra_hms": chart["ra"],
        "dec_dms": chart["dec"],
        "ra_deg": round(ra, 7),
        "dec_deg": round(dec, 7),
    }
    with open(os.path.join(outdir, "vsp_target_meta.json"), "w", encoding="utf-8") as fh:
        json.dump(chart_meta, fh, indent=2)

    apass_q = (
        "SELECT recno, raj2000, dej2000, nobs, mobs, b_v, vmag, e_vmag, u_e_vmag, "
        "bmag, e_bmag, r_mag, e_r_mag, i_mag, e_i_mag FROM external.apassdr9 WHERE "
        f"CONTAINS(POINT('ICRS',raj2000,dej2000),CIRCLE('ICRS',{ra},{dec},{args.radius}))=1 "
        "AND vmag IS NOT NULL ORDER BY vmag"
    )
    gaia_q = (
        "SELECT TOP %d source_id, ra, dec, phot_g_mean_mag, phot_g_mean_flux, phot_g_mean_flux_error, "
        "phot_bp_mean_mag, phot_rp_mean_mag, ruwe, phot_variable_flag FROM gaiadr3.gaia_source WHERE "
        f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{args.radius}))=1 "
        "AND phot_g_mean_mag IS NOT NULL AND phot_g_mean_mag < 19 ORDER BY phot_g_mean_mag" % GAIA_MAX_ROWS
    )

    datasets = [
        ("apass_dr9_tvcas_r0p5.csv", apass_q, "AAVSO APASS DR9 (Gaia archive TAP mirror: external.apassdr9)"),
        ("gaia_dr3_tvcas_r0p5.csv", gaia_q, "Gaia DR3 gaiadr3.gaia_source (ESA Gaia archive TAP)"),
    ]
    manifest = {"retrieved_utc": retrieved, "target": chart_meta, "datasets": []}
    for fname, query, label in datasets:
        text = tap_query(query)
        nrow = write_csv(os.path.join(outdir, fname), text)
        manifest["datasets"].append(
            {
                "file": fname,
                "label": label,
                "service": TAP_URL,
                "query": " ".join(query.split()),
                "rows": nrow,
                "sha256": sha256_file(os.path.join(outdir, fname)),
            }
        )
        print(f"{fname}: {nrow} rows")
    manifest["vsp_target_meta_sha256"] = sha256_file(os.path.join(outdir, "vsp_target_meta.json"))

    with open(os.path.join(outdir, "MANIFEST.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
    print("wrote", os.path.join(outdir, "MANIFEST.json"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
