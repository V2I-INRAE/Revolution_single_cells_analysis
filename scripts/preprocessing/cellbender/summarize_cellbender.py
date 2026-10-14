#!/usr/bin/env python3
"""Summarize every run attempt, retaining pending/failed samples rather than selecting latest runs."""

import argparse
import csv
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    rows = []
    for source in sorted((args.root / "data/raw_data").glob("*/*_RSEC_MolsPerCell_Unfiltered_MEX.zip")):
        sample = source.parent.name
        prep_path = args.root / "results/cellbender-preparation" / sample / "preparation.json"
        prep = json.loads(prep_path.read_text()) if prep_path.exists() else {}
        runs = sorted((args.root / "results/cellbender" / sample).glob("*/run.json"))
        for run_path in runs or [None]:
            config = json.loads(run_path.read_text()) if run_path else {}
            validation_path = run_path.parent / "validation.json" if run_path else None
            validation = json.loads(validation_path.read_text()) if validation_path and validation_path.exists() else {}
            review_path = run_path.parent / "qc_review.json" if run_path else None
            review = json.loads(review_path.read_text()) if review_path and review_path.exists() else {}
            warnings = validation.get("log_warning_lines", []) + [
                f"FPR {fpr}: {line}" for fpr, lines in validation.get("html_warning_lines", {}).items() for line in lines
            ]
            if config.get("validation_error"):
                warnings.append(config["validation_error"])
            rows.append({
                "sample": sample, "run_id": config.get("run_id"),
                "expected_cells": config.get("expected_cells", prep.get("expected_cells")),
                "total_droplets_included": config.get("total_droplets_included", prep.get("selected_total_droplets_included")),
                "candidate_total_droplets_included": prep.get("candidate_total_droplets_included"),
                "cells_called": validation.get("cells_called"),
                "pct_rhapsody_calls_retained": validation.get("pct_rhapsody_calls_retained"),
                "pct_cellbender_calls_in_rhapsody": validation.get("pct_cellbender_calls_in_rhapsody"),
                "jaccard_percent": validation.get("jaccard_percent"),
                "elbo_converged": review.get("elbo_converged"),
                "status": config.get("status", "prepared_not_run" if prep else "preparation_pending"),
                "qc_status": review.get("qc_status", validation.get("qc_status", "not_assessed")),
                "warnings": " | ".join(warnings), "output_dir": config.get("output_dir"),
            })
    output = args.root / "results/cellbender"
    output.mkdir(parents=True, exist_ok=True)
    (output / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    with (output / "summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {len(rows)} sample/run rows to {output / 'summary.tsv'}")


if __name__ == "__main__":
    main()
