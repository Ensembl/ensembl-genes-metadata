#!/usr/bin/env python3
"""Merge BAM files by tissue group and output a list of dictionaries."""

import argparse
import json
import os
from pathlib import Path
import subprocess

import pandas as pd


def merge_bams(tissue_group, output_dir):
    """Merge BAM files for a given tissue group."""
    tissue = tissue_group["tissue"].iloc[0]
    bam_files = tissue_group["bamFile"].tolist()
    merged_bam = os.path.join(output_dir, f"{tissue}_merged.bam")

    # samtools merge
    subprocess.run(["samtools", "merge", "-f", merged_bam] + bam_files, check=True)

    return {
        "taxon_id": tissue_group["taxon_id"].iloc[0],
        "genomeDir": tissue_group["genomeDir"].iloc[0],
        "gca": tissue_group["gca"].iloc[0],
        "platform": tissue_group["platform"].iloc[0],
        "paired": tissue_group["paired"].iloc[0],
        "tissue": tissue,
        "bamFile": merged_bam,
    }


def process(input_json, output_dir):
    """Process the input JSON file, group by tissue, and merge BAM files."""
    with open(input_json, "r", encoding="utf-8") as f:
        input_data = json.load(f)

    df = pd.DataFrame(input_data)

    # Group by tissue and write one CSV per group
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    grouped = df.groupby("tissue")

    output_list = []
    for tissue, group_df in grouped:
        csv_path = output_dir / f"{tissue}.csv"
        group_df.to_csv(csv_path, index=False)
        merged_entry = merge_bams(group_df, output_dir)
        output_list.append(merged_entry)

    # Final output: emit list of dictionaries
    output_json = output_dir / "merged_output.json"
    with open(output_json, "w", encoding="utf-8") as f:
        json.dump(output_list, f, indent=2)

    print(f"Output written to {output_json}")
    return output_list


__version__ = "1.0.0"

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="list of dictionaries")
    parser.add_argument(
        "--output", required=True, help="Directory for CSVs and merged BAMs"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    args = parser.parse_args()

    process(args.input, args.output)
