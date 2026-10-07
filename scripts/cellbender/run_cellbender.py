#!/usr/bin/env python3
"""Run one reviewed sample with CellBender inside its own Slurm allocation."""

import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import zipfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample")
    parser.add_argument("--epochs", type=int, default=150)
    parser.add_argument("--learning-rate", type=float, default=1e-4)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    preparation_path = root / "results/cellbender-preparation" / args.sample / "preparation.json"
    preparation = json.loads(preparation_path.read_text())
    expected = preparation["expected_cells"]
    total_droplets = preparation["selected_total_droplets_included"]
    if total_droplets is None:
        raise ValueError("Review the barcode-rank curve and select total droplets in preparation.json first")
    if preparation["validation"] != "passed" or not expected < total_droplets < preparation["barcodes"]:
        raise ValueError("Preparation failed or reviewed total-droplets-included is invalid")
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Submit this script through Slurm; do not train on a login node")
    threads = int(os.environ["SLURM_CPUS_PER_TASK"])
    task = os.environ.get("SLURM_ARRAY_TASK_ID")
    run_id = f"job-{os.environ['SLURM_JOB_ID']}" + (f"-task-{task}" if task is not None else "")
    output = root / "results/cellbender" / args.sample / run_id
    work = root / "data/cellbender" / args.sample / run_id / f"{args.sample}_unfiltered_MEX"
    output.mkdir(parents=True, exist_ok=False)
    work.mkdir(parents=True, exist_ok=False)
    source = Path(preparation["inputs"]["unfiltered_mex"]["path"])
    if source.stat().st_size != preparation["inputs"]["unfiltered_mex"]["bytes"] or source.stat().st_mtime_ns != preparation["inputs"]["unfiltered_mex"]["mtime_ns"]:
        raise ValueError("Input archive changed since preparation")
    with zipfile.ZipFile(source) as archive:
        if sorted(archive.namelist()) != ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]:
            raise ValueError("Unexpected archive members")
        archive.extractall(work)
    command = [str(Path(sys.executable).parent / "cellbender"), "remove-background",
               "--input", str(work), "--output", str(output / f"{args.sample}_cellbender.h5"),
               "--expected-cells", str(expected), "--total-droplets-included", str(total_droplets),
               "--fpr", "0.0", "0.01", "0.05", "--epochs", str(args.epochs),
               "--learning-rate", str(args.learning_rate), "--cpu-threads", str(threads)]
    config = {
        "sample": args.sample, "run_id": run_id, "preparation_json": str(preparation_path),
        "input_mex": str(work), "output_dir": str(output), "expected_cells": expected,
        "total_droplets_included": total_droplets, "parameter_review": preparation["parameter_review"],
        "epochs": args.epochs, "learning_rate": args.learning_rate, "cpu_threads": threads,
        "fpr": [0.0, 0.01, 0.05], "main_fpr": 0.01, "command": command,
        "versions": {p: importlib.metadata.version(p) for p in ["cellbender", "torch", "scvi-tools", "numpy", "pandas"]},
        "report_source_sha256": hashlib.sha256(
            Path(importlib.util.find_spec("cellbender.remove_background.report").origin).read_bytes()).hexdigest(),
        "report_notebook_sha256": hashlib.sha256(
            Path(importlib.util.find_spec("cellbender.remove_background.report").origin).with_suffix(".ipynb").read_bytes()).hexdigest(),
        "status": "running", "started_unix": time.time(),
    }
    config_path = output / "run.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n")
    environment = os.environ.copy()
    environment["PATH"] = str(Path(sys.executable).parent) + os.pathsep + environment["PATH"]
    print(json.dumps(config, indent=2), flush=True)
    result = subprocess.run(command, cwd=output, env=environment)
    config.update(status="completed_pending_qc" if result.returncode == 0 else "failed",
                  returncode=result.returncode, elapsed_seconds=time.time() - config["started_unix"])
    config_path.write_text(json.dumps(config, indent=2) + "\n")
    if result.returncode:
        raise SystemExit(result.returncode)
    subprocess.run([sys.executable, str(root / "scripts/cellbender/validate_cellbender.py"), str(output)],
                   env=environment, check=True)


if __name__ == "__main__":
    main()
