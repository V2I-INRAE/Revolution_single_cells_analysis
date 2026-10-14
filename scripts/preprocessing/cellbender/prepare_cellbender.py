#!/usr/bin/env python3
"""Read-only BD Rhapsody MEX checks and per-sample barcode-rank diagnostics."""

import argparse
import csv
import gzip
import itertools
import json
from pathlib import Path
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot_rank(sample, ranked, expected, droplets, label, output):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    rank = np.arange(1, len(ranked) + 1)
    for ax in axes:
        ax.plot(rank, ranked, color="#3A5BA0", linewidth=1)
        ax.set_yscale("log")
        ax.axvline(expected, color="#D4753C", linestyle="--", label=f"BD putative cells: {expected:,}")
        ax.axvline(droplets, color="#5A8F5A", linestyle=":",
                   label=f"{label} M: {droplets:,} ({ranked[droplets-1]:,} molecules)")
        ax.set_xlabel("Barcode rank (descending molecules)")
        ax.set_ylabel("RSEC molecules per barcode (log scale)")
        ax.grid(alpha=0.15)
    axes[0].set_xscale("log")
    axes[0].set_title("All unfiltered barcodes; log rank")
    axes[0].legend(fontsize=8)
    axes[1].set_xlim(1, min(len(ranked), max(50000, int(droplets * 1.5))))
    axes[1].set_title("Cell / background transition; linear rank")
    fig.suptitle(f"{sample}: unfiltered BD Rhapsody WTA barcode-rank diagnostic")
    fig.tight_layout()
    fig.savefig(output / "barcode_rank.png", dpi=300)
    plt.close(fig)


def prepare_sample(source, destination):
    sample = source.name
    paths = {
        "unfiltered_mex": source / f"{sample}_RSEC_MolsPerCell_Unfiltered_MEX.zip",
        "filtered_mex": source / f"{sample}_RSEC_MolsPerCell_MEX.zip",
        "metrics": source / f"{sample}_Metrics_Summary.csv",
    }
    with paths["metrics"].open(newline="") as handle:
        rows = list(csv.reader(handle))
    section = next(i for i, row in enumerate(rows) if row == ["#Cells#"])
    header = rows[section + 1]
    cell_rows = [row for row in rows[section + 2:] if row and not row[0].startswith("#")]
    metrics = [dict(zip(header, row)) for row in cell_rows if len(row) == len(header)]
    metrics = [row for row in metrics if row.get("Bioproduct_Type") == "mRNA"]
    if len(metrics) != 1:
        raise ValueError(f"{sample}: expected one mRNA row in #Cells#")
    expected = int(metrics[0]["Putative_Cell_Count"])

    with zipfile.ZipFile(paths["unfiltered_mex"]) as archive:
        if sorted(archive.namelist()) != ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]:
            raise ValueError(f"{sample}: unexpected unfiltered archive members")
        with gzip.open(archive.open("features.tsv.gz"), "rt") as handle:
            features = [line.rstrip("\r\n").split("\t") for line in handle]
        if not all(len(row) == 3 and row[2] == "Gene Expression" for row in features):
            raise ValueError(f"{sample}: expected three-column Gene Expression features")
        if len({row[0] for row in features}) != len(features):
            raise ValueError(f"{sample}: duplicate feature IDs")
        with gzip.open(archive.open("barcodes.tsv.gz"), "rt") as handle:
            barcodes = [line.strip() for line in handle]
        if len(set(barcodes)) != len(barcodes) or not all(b.isdigit() for b in barcodes):
            raise ValueError(f"{sample}: duplicate or nonnumeric bead barcodes")
        with gzip.open(archive.open("matrix.mtx.gz"), "rt") as handle:
            if handle.readline().strip() != "%%MatrixMarket matrix coordinate integer general":
                raise ValueError(f"{sample}: unexpected MatrixMarket header")
            line = handle.readline()
            while line.startswith("%"):
                line = handle.readline()
            nfeatures, nbarcodes, nnz = map(int, line.split())
            if (nfeatures, nbarcodes) != (len(features), len(barcodes)):
                raise ValueError(f"{sample}: matrix dimensions do not match feature/barcode files")
            counts = np.zeros(nbarcodes, dtype=np.int64)
            entries, total = 0, 0
            while lines := list(itertools.islice(handle, 200000)):
                values = np.loadtxt(lines, dtype=np.int64, ndmin=2)
                if values.shape[1] != 3 or np.any(values[:, 2] < 0):
                    raise ValueError(f"{sample}: invalid matrix entries")
                if np.any((values[:, 0] < 1) | (values[:, 0] > nfeatures)) or np.any(
                    (values[:, 1] < 1) | (values[:, 1] > nbarcodes)
                ):
                    raise ValueError(f"{sample}: out-of-bounds coordinates")
                np.add.at(counts, values[:, 1] - 1, values[:, 2])
                entries += len(values)
                total += int(values[:, 2].sum())
    if entries != nnz or int(counts.sum()) != total or np.any(counts < 1):
        raise ValueError(f"{sample}: nnz/column sums inconsistent with unfiltered input")
    with zipfile.ZipFile(paths["filtered_mex"]) as archive:
        with gzip.open(archive.open("barcodes.tsv.gz"), "rt") as handle:
            filtered = [line.strip() for line in handle]
    if len(set(filtered)) != len(filtered) or len(filtered) != expected:
        raise ValueError(f"{sample}: filtered barcode count differs from metrics or is duplicated")
    if not set(filtered).issubset(set(barcodes)):
        raise ValueError(f"{sample}: filtered calls missing from unfiltered input")

    # Starting point only: the curve must be inspected before selecting M.
    candidate = min(nbarcodes - 1, int(np.ceil(2 * expected / 1000)) * 1000)
    if not 0 < expected < candidate < nbarcodes:
        raise ValueError(f"{sample}: insufficient empty barcodes for this starting point")
    order = np.argsort(-counts, kind="stable")
    ranked = counts[order]
    ranks = sorted({expected, candidate, *range(5000, min(nbarcodes, 60000) + 1, 5000)})
    summary = {
        "sample": sample,
        "inputs": {key: {"path": str(path.resolve()), "bytes": path.stat().st_size,
                          "mtime_ns": path.stat().st_mtime_ns} for key, path in paths.items()},
        "validation": "passed",
        "features": nfeatures,
        "barcodes": nbarcodes,
        "nnz": nnz,
        "total_molecules": total,
        "filtered_barcodes": len(filtered),
        "duplicate_feature_names": nfeatures - len({row[1] for row in features}),
        "expected_cells": expected,
        "expected_cells_source": "Metrics Summary #Cells# mRNA Putative_Cell_Count",
        "candidate_total_droplets_included": candidate,
        "selected_total_droplets_included": None,
        "parameter_review": "pending visual review; candidate is approximately 2x expected cells",
        "molecules_at_rank": {str(rank): int(ranked[rank - 1]) for rank in ranks},
    }
    output = destination / sample
    output.mkdir(parents=True, exist_ok=False)
    with gzip.open(output / "barcode_counts.tsv.gz", "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["rank", "barcode", "molecules", "rhapsody_called"])
        called = set(filtered)
        writer.writerows((rank, barcodes[i], int(counts[i]), int(barcodes[i] in called))
                         for rank, i in enumerate(order, 1))
    (output / "preparation.json").write_text(json.dumps(summary, indent=2) + "\n")
    plot_rank(sample, ranked, expected, candidate, "Candidate", output)
    print(json.dumps(summary), flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-data", type=Path, default=Path("data/raw_data"))
    parser.add_argument("--outdir", type=Path, default=Path("results/cellbender-preparation"))
    parser.add_argument("--samples", nargs="+")
    parser.add_argument("--review-total-droplets", type=int)
    parser.add_argument("--review-note")
    args = parser.parse_args()
    if args.review_total_droplets is not None:
        if not args.samples or len(args.samples) != 1 or not args.review_note:
            parser.error("Review requires exactly one --samples value and --review-note")
        output = args.outdir / args.samples[0]
        path = output / "preparation.json"
        summary = json.loads(path.read_text())
        selected = args.review_total_droplets
        if not summary["expected_cells"] < selected < summary["barcodes"]:
            raise ValueError("Reviewed rank must exceed expected cells and leave surely-empty barcodes")
        ranked = np.loadtxt(output / "barcode_counts.tsv.gz", delimiter="\t", skiprows=1, usecols=2, dtype=np.int64)
        summary.update(selected_total_droplets_included=selected, selected_rank_molecules=int(ranked[selected - 1]),
                       parameter_review=args.review_note)
        plot_rank(summary["sample"], ranked, summary["expected_cells"], selected, "Selected", output)
        path.write_text(json.dumps(summary, indent=2) + "\n")
        print(json.dumps(summary, indent=2))
        return
    if args.review_note:
        parser.error("--review-note requires --review-total-droplets")
    sources = sorted(path.parent for path in args.raw_data.glob("*/*_RSEC_MolsPerCell_Unfiltered_MEX.zip"))
    if args.samples:
        if set(args.samples) - {source.name for source in sources}:
            raise ValueError("Unknown requested sample")
        sources = [source for source in sources if source.name in args.samples]
    summaries = [prepare_sample(source, args.outdir) for source in sources]
    manifest = args.outdir / ("manifest.tsv" if not args.samples else "manifest-selected.tsv")
    with manifest.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "expected_cells", "candidate_total_droplets_included", "preparation_json"])
        writer.writerows((row["sample"], row["expected_cells"], row["candidate_total_droplets_included"],
                          str((args.outdir / row["sample"] / "preparation.json").resolve())) for row in summaries)


if __name__ == "__main__":
    main()
