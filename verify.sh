#!/usr/bin/env bash
# One command, fully offline (stdlib Python 3 only, no API keys): proves the formulas,
# the shipped data, and that the published results are what the code produces today.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"

echo "== formula selftest =="
"$PY" tools/error_budget.py selftest

echo
echo "== unit tests =="
"$PY" -m unittest discover -s tests

echo
echo "== data integrity (sha256 against data/MANIFEST.json) =="
"$PY" - <<'EOF'
import hashlib, json, os, sys

man = json.load(open(os.path.join("data", "MANIFEST.json"), encoding="utf-8"))
targets = [(d["file"], d["sha256"]) for d in man["datasets"]]
targets.append(("vsp_target_meta.json", man["vsp_target_meta_sha256"]))
bad = []
for name, want in targets:
    got = hashlib.sha256(open(os.path.join("data", name), "rb").read()).hexdigest()
    if got != want:
        bad.append("%s: %s != %s" % (name, got, want))
if bad:
    print("\n".join(bad))
    sys.exit(1)
print("OK   %d data files match data/MANIFEST.json" % len(targets))
EOF

echo
echo "== worked examples still reproduce byte-for-byte =="
"$PY" tools/run_examples.py --check

echo
echo "ALL CHECKS PASSED"
