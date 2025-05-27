import os
import subprocess
import pandas as pd
from pathlib import Path
import argparse
import json

def merge_bams(tissue_group, output_dir):
    tissue = tissue_group['tissue'].iloc[0]
    bam_files = tissue_group['bamFile'].tolist()
    merged_bam = os.path.join(output_dir, f"{tissue}_merged.bam")

    # samtools merge
    subprocess.run(['samtools', 'merge', '-f', merged_bam] + bam_files, check=True)

    return {
        "taxon_id": tissue_group['taxon_id'].iloc[0],
        "genomeDir": tissue_group['genomeDir'].iloc[0],
        "gca": tissue_group['gca'].iloc[0],
        "platform": tissue_group['platform'].iloc[0],
        "paired": tissue_group['paired'].iloc[0],
        "tissue": tissue,
        "bamFile": merged_bam
    }

def process(input_json, output_dir):
    with open(input_json, 'r') as f:
        input_data = json.load(f)

    df = pd.DataFrame(input_data)

    # Group by tissue and write one CSV per group
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    grouped = df.groupby('tissue')

    output_list = []
    for tissue, group_df in grouped:
        csv_path = output_dir / f"{tissue}.csv"
        group_df.to_csv(csv_path, index=False)
        merged_entry = merge_bams(group_df, output_dir)
        output_list.append(merged_entry)

    # Final output: emit list of dictionaries
    output_json = output_dir / "merged_output.json"
    with open(output_json, 'w') as f:
        json.dump(output_list, f, indent=2)

    print(f"Output written to {output_json}")
    return output_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="list of dictionaries")
    parser.add_argument("--output", required=True, help="Directory for CSVs and merged BAMs")
    args = parser.parse_args()

    process(args.input, args.output)
