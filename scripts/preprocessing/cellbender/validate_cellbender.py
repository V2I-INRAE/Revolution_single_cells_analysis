#!/usr/bin/env python3
"""Validate one completed run; numerical success is not automatic QC acceptance."""

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import zipfile

from bs4 import BeautifulSoup
from cellbender.remove_background.data.io import load_data
import h5py
import numpy as np
import pandas as pd
import scanpy as sc


AMBIENT_GENES = "COX1 COX2 COX3 ATP6 ND1 ND2 ND4 ND5 CYTB JCHAIN IGHM SFTPC SFTPA1 CXCL8 S100A8 S100A12 CSF3R GPR84 RETN G0S2 SDS STEAP4 TREM1".split()


def regenerate_reports(output):
    from cellbender.remove_background import report

    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Regenerate reports inside a Slurm allocation")
    config_path = output / "run.json"
    config = json.loads(config_path.read_text())
    if config.get("returncode") != 0:
        raise ValueError("Cannot repair reports for an unsuccessful CellBender run")
    repair_dir = output / f"report-regeneration-{os.environ['SLURM_JOB_ID']}"
    repair_dir.mkdir(exist_ok=False)
    hashes = {}
    for path in output.glob(f"{config['sample']}_cellbender_FPR_*.h5"):
        with path.open("rb") as handle:
            hashes[path.name] = hashlib.file_digest(handle, "sha256").hexdigest()
    os.environ["INPUT_FILE"] = config["input_mex"]
    os.environ.pop("TRUTH_FILE", None)
    original_cwd = Path.cwd()
    try:
        os.chdir(repair_dir)
        for fpr in config["fpr"]:
            stem = f"{config['sample']}_cellbender_FPR_{fpr}"
            shutil.copy2(output / f"{stem}_report.html", repair_dir / f"{stem}_report.original.html")
            os.environ["OUTPUT_FILE"] = str(output / f"{stem}.h5")
            regenerated = repair_dir / f"{stem}_report.html"
            report.run_notebook_make_html(str(Path(report.__file__).with_suffix(".ipynb")), str(regenerated))
            html = BeautifulSoup(regenerated.read_text(), "html.parser")
            text = html.get_text(" ", strip=True)
            if "Traceback (most recent call last)" in text or "Summary of warnings:" not in text or not html.select('img[src^="data:image/"]'):
                raise ValueError(f"FPR {fpr}: regenerated report still failed")
    finally:
        os.chdir(original_cwd)
    for name, expected_hash in hashes.items():
        with (output / name).open("rb") as handle:
            if hashlib.file_digest(handle, "sha256").hexdigest() != expected_hash:
                raise ValueError(f"Report regeneration changed count file {name}")
    for fpr in config["fpr"]:
        name = f"{config['sample']}_cellbender_FPR_{fpr}_report.html"
        shutil.copy2(repair_dir / name, output / name)
    if "report_repair" in config:
        config.setdefault("report_repair_history", []).append(config["report_repair"])
    config["report_repair"] = {
        "previous_validation_error": config.pop("validation_error", None),
        "report_source_sha256": hashlib.sha256(Path(report.__file__).read_bytes()).hexdigest(),
        "report_notebook_sha256": hashlib.sha256(Path(report.__file__).with_suffix(".ipynb").read_bytes()).hexdigest(),
        "counts_sha256_unchanged": hashes, "archive_directory": str(repair_dir),
    }
    config["status"] = "completed_pending_qc"
    config_path.write_text(json.dumps(config, indent=2) + "\n")


def validate(output):
    config = json.loads((output / "run.json").read_text())
    if config["status"] not in ["completed_pending_qc", "validated_pending_qc", "accepted"]:
        raise ValueError("CellBender did not complete successfully")
    sample = config["sample"]
    prefix = f"{sample}_cellbender"
    preparation = json.loads(Path(config["preparation_json"]).read_text())
    called = (output / f"{prefix}_cell_barcodes.csv").read_text().splitlines()
    if not called or len(set(called)) != len(called):
        raise ValueError("Empty or duplicated CellBender calls")
    with zipfile.ZipFile(preparation["inputs"]["filtered_mex"]["path"]) as archive:
        with gzip.open(archive.open("barcodes.tsv.gz"), "rt") as handle:
            bd_called = {line.strip() for line in handle}
    overlap = len(set(called) & bd_called)
    log = (output / f"{prefix}.log").read_text()
    priors = re.search(r"Using (\d+) probable cell barcodes, plus an additional (\d+) barcodes, and (\d+) empty droplets", log)
    if not priors or any(int(value) == 0 for value in priors.groups()):
        raise ValueError("Missing or zero probable cells/additional barcodes/empty droplets")
    if "Completed remove-background." not in log or not (output / f"{prefix}.pdf").is_file():
        raise ValueError("Missing completion message or PDF diagnostic")

    # Official CellBender input reader: rows are barcodes, columns are features.
    raw = load_data(config["input_mex"])
    raw_barcodes = np.asarray(raw["barcodes"]).astype(str)
    raw_ids = np.asarray(raw["gene_ids"]).astype(str)
    raw_names = np.asarray(raw["gene_names"]).astype(str)
    if raw["matrix"].shape != (preparation["barcodes"], preparation["features"]) or int(raw["matrix"].sum()) != preparation["total_molecules"]:
        raise ValueError("Independent input-reader check disagrees with preparation")
    barcode_index = {barcode: i for i, barcode in enumerate(raw_barcodes)}
    raw_cells = raw["matrix"].tocsr()[[barcode_index[b] for b in called]]
    del raw
    raw_totals = np.asarray(raw_cells.sum(axis=0)).ravel()
    bd_mask = np.asarray([b in bd_called for b in called])
    group_masks = {"rhapsody_called": bd_mask, "additional_calls": ~bd_mask}
    raw_by_group = {group: np.asarray(raw_cells[mask].sum(axis=0)).ravel()
                    for group, mask in group_masks.items()}
    per_fpr = []
    per_group = []
    ambient_rows = []
    ambient_group_rows = []
    report_warnings = {}
    for fpr in config["fpr"]:
        stem = f"{prefix}_FPR_{fpr}"
        filtered_path = output / f"{stem}_filtered.h5"
        corrected = sc.read_10x_h5(filtered_path)
        if list(corrected.obs_names) != called or not np.array_equal(corrected.var["gene_ids"].astype(str), raw_ids):
            raise ValueError(f"FPR {fpr}: filtered barcodes/features not aligned with input")
        if not np.array_equal(np.asarray(corrected.var_names).astype(str), raw_names):
            raise ValueError(f"FPR {fpr}: gene names changed")
        if not np.all(np.isfinite(corrected.X.data)) or np.any(corrected.X.data < 0) or np.any(corrected.X.data != np.floor(corrected.X.data)):
            raise ValueError(f"FPR {fpr}: counts are not finite nonnegative integers")
        matrix = corrected.X.tocsr().astype(np.int64)
        if np.any((matrix - raw_cells).data > 0):
            raise ValueError(f"FPR {fpr}: correction added counts")
        full = load_data(str(output / f"{stem}.h5"))
        if not np.array_equal(np.asarray(full["barcodes"]).astype(str), raw_barcodes) or not np.array_equal(np.asarray(full["gene_ids"]).astype(str), raw_ids):
            raise ValueError(f"FPR {fpr}: full H5 changed original barcode/feature order")
        if (full["matrix"].tocsr()[[barcode_index[b] for b in called]] != matrix).nnz:
            raise ValueError(f"FPR {fpr}: full and filtered H5 disagree")
        del full
        metrics = dict(csv.reader((output / f"{stem}_metrics.csv").open(newline="")))
        if int(float(metrics["found_cells"])) != len(called):
            raise ValueError(f"FPR {fpr}: metrics and barcode calls disagree")
        report = BeautifulSoup((output / f"{stem}_report.html").read_text(), "html.parser")
        report_text = report.get_text("\n", strip=True)
        if "Traceback (most recent call last)" in report_text or report.select(".output_error"):
            raise ValueError(f"FPR {fpr}: notebook errors in HTML report")
        if "Summary of warnings:" not in report_text:
            raise ValueError(f"FPR {fpr}: HTML report did not reach its warning summary")
        if not report.select('img[src^="data:image/"]'):
            raise ValueError(f"FPR {fpr}: HTML report contains no rendered diagnostic plots")
        summary_text = report_text.split("Summary of warnings:", 1)[1]
        report_warnings[str(fpr)] = [line for line in summary_text.splitlines() if line not in ["¶", "None."]]
        totals = np.asarray(matrix.sum(axis=0)).ravel()
        per_fpr.append({"fpr": fpr, "raw_molecules_in_called_cells": int(raw_totals.sum()),
                        "corrected_molecules_in_called_cells": int(totals.sum()),
                        "fraction_removed": float(1 - totals.sum() / raw_totals.sum())})
        for group, mask in group_masks.items():
            before_group = raw_by_group[group]
            after_group = np.asarray(matrix[mask].sum(axis=0)).ravel()
            per_group.append({"fpr": fpr, "barcode_group": group, "barcodes": int(mask.sum()),
                              "raw_molecules": int(before_group.sum()),
                              "corrected_molecules": int(after_group.sum()),
                              "fraction_removed": float(1 - after_group.sum() / before_group.sum()) if before_group.sum() else None})
            for gene in AMBIENT_GENES:
                matches = np.flatnonzero(raw_names == gene)
                before = int(before_group[matches].sum()) if len(matches) else None
                after = int(after_group[matches].sum()) if len(matches) else None
                ambient_group_rows.append({"fpr": fpr, "barcode_group": group, "gene": gene,
                                           "matched_features": len(matches), "feature_ids": ";".join(raw_ids[matches]),
                                           "raw_molecules": before, "corrected_molecules": after,
                                           "fraction_removed": 1 - after / before if before else None})
        for gene in AMBIENT_GENES:
            matches = np.flatnonzero(raw_names == gene)
            # Do not infer pig gene aliases or silently substitute other identifiers.
            before = int(raw_totals[matches].sum()) if len(matches) else None
            after = int(totals[matches].sum()) if len(matches) else None
            ambient_rows.append({"fpr": fpr, "gene": gene, "matched_features": len(matches),
                                 "feature_ids": ";".join(raw_ids[matches]), "raw_molecules": before,
                                 "corrected_molecules": after,
                                 "fraction_removed": 1 - after / before if before else None})
        pd.DataFrame({"feature_id": raw_ids, "gene": raw_names, "raw_molecules": raw_totals,
                      "corrected_molecules": totals}).to_csv(output / f"gene_counts_FPR_{fpr}.tsv.gz", sep="\t", index=False)
        if fpr == config["main_fpr"]:
            pd.DataFrame({"barcode": called, "raw_molecules": np.asarray(raw_cells.sum(axis=1)).ravel(),
                          "corrected_molecules": np.asarray(matrix.sum(axis=1)).ravel(),
                          "raw_genes_detected": np.asarray((raw_cells > 0).sum(axis=1)).ravel(),
                          "corrected_genes_detected": np.asarray((matrix > 0).sum(axis=1)).ravel(),
                          "rhapsody_called": bd_mask}).to_csv(output / "called_cell_counts.tsv.gz", sep="\t", index=False)

    with h5py.File(output / f"{prefix}_FPR_0.01.h5") as handle:
        probability = handle["droplet_latents/cell_probability"][:]
        analyzed_indices = handle["droplet_latents/barcode_indices_for_latents"][:]
        if len(probability) != len(analyzed_indices) or set(raw_barcodes[analyzed_indices[probability > 0.5]]) != set(called):
            raise ValueError("Cell probabilities disagree with called barcode CSV")
        if np.any(~np.isfinite(probability)) or np.any((probability < 0) | (probability > 1)):
            raise ValueError("Invalid cell probabilities")
        # The writer flattens learning_curve as learning_curve_train_elbo etc.
        train_elbo = handle["metadata/learning_curve_train_elbo"][:]
        train_epoch = handle["metadata/learning_curve_train_epoch"][:]
        test_elbo = handle["metadata/learning_curve_test_elbo"][:]
        test_epoch = handle["metadata/learning_curve_test_epoch"][:]
    if len(train_elbo) != config["epochs"] or np.any(~np.isfinite(train_elbo)) or np.any(~np.isfinite(test_elbo)):
        raise ValueError("Incomplete or nonfinite ELBO history")
    pd.DataFrame({"epoch": train_epoch, "train_elbo": train_elbo}).to_csv(output / "train_elbo.tsv", sep="\t", index=False)
    pd.DataFrame({"epoch": test_epoch, "test_elbo": test_elbo}).to_csv(output / "test_elbo.tsv", sep="\t", index=False)
    pd.DataFrame(ambient_rows).to_csv(output / "ambient_gene_validation.tsv", sep="\t", index=False)
    pd.DataFrame(ambient_group_rows).to_csv(output / "ambient_gene_by_barcode_group.tsv", sep="\t", index=False)
    pd.DataFrame(per_group).to_csv(output / "barcode_group_validation.tsv", sep="\t", index=False)
    result = {
        "sample": sample, "run_id": config["run_id"], "numerical_validation": "passed",
        "expected_cells": config["expected_cells"], "total_droplets_included": config["total_droplets_included"],
        "cells_called": len(called), "rhapsody_cells": len(bd_called), "overlap_cells": overlap,
        "pct_rhapsody_calls_retained": 100 * overlap / len(bd_called),
        "pct_cellbender_calls_in_rhapsody": 100 * overlap / len(called),
        "jaccard_percent": 100 * overlap / len(set(called) | bd_called),
        "priors": dict(zip(["probable_cells", "additional_barcodes", "empty_droplets"], map(int, priors.groups()))),
        "probabilities_between_0.1_and_0.9": int(((probability > 0.1) & (probability < 0.9)).sum()),
        "per_fpr": per_fpr, "per_fpr_by_barcode_group": per_group, "html_warning_lines": report_warnings,
        "log_warning_lines": [line for line in log.splitlines() if "warning" in line.lower()],
        "elbo_converged": None, "qc_status": "pending_visual_and_warning_review",
        "biological_validation": "aggregate gene comparisons only; cell-type-specific preservation requires validated labels",
        "checkpoint_policy": "retained pending QC and FPR review",
    }
    (output / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    # Explicitly select the agreed main FPR without copying large H5 files.
    for suffix in [".h5", "_filtered.h5", "_metrics.csv", "_report.html"]:
        alias = output / f"{prefix}{suffix}"
        target = f"{prefix}_FPR_0.01{suffix}"
        if not alias.exists():
            alias.symlink_to(target)
    config["status"] = "validated_pending_qc"
    (output / "run.json").write_text(json.dumps(config, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--regenerate-reports", action="store_true",
                        help="Repair HTML reports without retraining or changing counts, then validate")
    args = parser.parse_args()
    output = args.output.resolve()
    try:
        if args.regenerate_reports:
            regenerate_reports(output)
        validate(output)
    except Exception as error:
        config_path = output / "run.json"
        config = json.loads(config_path.read_text())
        config.update(status="validation_failed", validation_error=str(error))
        config_path.write_text(json.dumps(config, indent=2) + "\n")
        raise


if __name__ == "__main__":
    main()
