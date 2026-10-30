#!/usr/bin/env python3
"""Recompute statistics from the original Borlido/Botti 2019 workbook.

Download DOI:10.1021/acs.jctc.9b00322.s002 separately. Requires openpyxl.
No new electronic-structure calculations are performed by this script.
"""
import argparse
import hashlib
import json
from pathlib import Path
import openpyxl

SHA256 = "21fc5c597411093c98e6209f3a94e907f6e81b6e27ee87611bcd73bd59ac2db0"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("workbook", type=Path)
    args = ap.parse_args()
    if hashlib.sha256(args.workbook.read_bytes()).hexdigest() != SHA256:
        ap.error("Workbook differs from the pinned original 2019 supplement")
    workbook = openpyxl.load_workbook(args.workbook, data_only=True, read_only=True)
    rows = list(workbook.active.values)[2:]
    result = {"source": "https://doi.org/10.1021/acs.jctc.9b00322.s002", "sha256": SHA256}
    for name, col in [("PBE", 5), ("BJ", 10), ("mBJ", 11), ("HSE06", 13)]:
        pairs = [(r[col], r[18]) for r in rows
                 if isinstance(r[col], (int, float)) and isinstance(r[18], (int, float))]
        errors = [a-b for a, b in pairs]
        result[name] = {"n": len(pairs), "MAE_eV": sum(abs(e) for e in errors)/len(pairs),
                        "ME_eV": sum(errors)/len(pairs),
                        "MAPE_percent": 100*sum(abs(a-b)/b for a, b in pairs)/len(pairs),
                        "false_metals_below_0.01_eV": sum(a < .01 for a, b in pairs)}
    result["Si"] = [{"PBE_eV": r[5], "mBJ_eV": r[11], "experiment_eV": r[18]}
                    for r in rows if r[0] == "Si"]
    result["missing_mBJ"] = [r[0] for r in rows if not isinstance(r[11], (int, float))]
    workbook.close()
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
