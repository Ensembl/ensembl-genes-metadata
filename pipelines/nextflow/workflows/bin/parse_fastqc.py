#!/usr/bin/env python3
# pylint: disable=missing-module-docstring
#  See the NOTICE file distributed with this work for additional information
#  regarding copyright ownership.
#
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

from pathlib import Path
import re
import json
import logging
import argparse
from typing import Any, Dict, Optional


def parse_fastqc_summary(summary_path: Path) -> Dict[str, str]:
    """Parse the summary.txt file from FASTQC output.

    Args:
        summary_path (str): Path to the summary.txt file.

    Returns:
        Dict[str, str]: Dictionary containing module statuses.
    """
    summary_data: Dict[str, str] = {}
    if not summary_path.exists():
        logging.error("Summary file not found: %s", summary_path)
        return summary_data

    with summary_path.open("r") as summary_file:
        for line in summary_file:
            test_results = line.strip().split("\t")
            if len(test_results) == 3:
                status, module, _ = test_results
                module = module.replace(" ", "_").lower()
                summary_data[module] = status

    return summary_data


def parse_fastqc_data(fastqc_data_path: Path) -> Dict[str, int]:
    """Parse the fastqc_data.txt file from FASTQC output.

    Args:
        fastqc_data_path (str): Path to the fastqc_data.txt file.

    Returns:
        Dict[str, int]: Dictionary containing total sequences,
        sequence length, and GC content.
    """
    fastqc_data: Dict[str, int] = {}
    if not fastqc_data_path.exists():
        logging.error("FASTQC data file not found: %s", fastqc_data_path)
        return fastqc_data

    with fastqc_data_path.open("r") as f:
        data = f.read()

    total_sequences_pattern = r"Total Sequences\s+(\d+)"
    sequence_length_pattern = r"Sequence length\s+(\d+)"
    gc_content_pattern = r"%GC\s+(\d+)"

    total_sequences_match = re.search(total_sequences_pattern, data)
    if total_sequences_match:
        fastqc_data["total_sequences"] = int(total_sequences_match.group(1))

    gc_content_match = re.search(gc_content_pattern, data)
    if gc_content_match:
        fastqc_data["gc_content"] = int(gc_content_match.group(1))

    sequence_length_match = re.search(sequence_length_pattern, data)
    if sequence_length_match:
        fastqc_data["sequence_length"] = int(sequence_length_match.group(1))
    print(fastqc_data)
    return fastqc_data


def convert_to_json(run_accession_dir: Path, data_file_json: str, run_id: int) -> None:
    """Convert FASTQC output to JSON format.

    Args:
        fastqc_dir (str): Path to the folder containing FASTQC results.
    """
    table_data_files: Dict[str, list[Dict[str, str]]] = {"data_files": []}
    with open(data_file_json, "r") as input_file:
        data = json.load(input_file)
    fastq_files = [str(file.name) for file in run_accession_dir.iterdir() if file.suffix in ("_fastqc")]

    # fastq_files = [file for file in os.listdir(run_accession_dir) if file.endswith((".fastq", ".gz"))]
    for fastq_file in fastq_files:
        summary_path = Path(run_accession_dir) / fastq_file / "summary.txt"
        # f'{fastq_file.replace(".fastq.gz", "_fastqc")}'
        fastqc_data_path = Path(
            Path(run_accession_dir) / fastq_file / "fastqc_data.txt",
        )
        # f'{fastq_file.replace(".fastq.gz", "_fastqc")}'
        # Find the dictionary with the matching name
        run_accession_dict: Optional[Dict[str, str]] = next(
            (item for item in data["data_files"] if item["file_name"] == fastq_file.replace("_fastqc", "")),
            None,
        )
        if run_accession_dict is None:
            logging.warning("No matching entry found for %s. Skipping this file.", fastq_file)
            continue

        summary_data = parse_fastqc_summary(summary_path)
        fastqc_data = parse_fastqc_data(fastqc_data_path)

        data_file: Dict[str, Any] = {
            "run_id": run_id,
            **run_accession_dict,
            **summary_data,
            **fastqc_data,
        }
        print(data_file)
        table_data_files["data_files"].append(data_file)
    json_data_files = json.dumps(table_data_files)

    with open("complete_insert_into_data_file.json", "w") as file:
        file.write(json_data_files)

    logging.info("Summary data converted to JSON.")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description="Parse FASTQC output and create JSON with specified keys")
    parser.add_argument("--fastqc_results_path", required=True, type=str, help="Path to the FASTQC data file")
    parser.add_argument("--data_file_json", required=True, type=str, help="Output directory")
    parser.add_argument("--run_id", required=True, type=int, help="Table identifier for run_accession")
    return parser.parse_args()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    args = parse_args()
    convert_to_json(Path(args.fastqc_results_path), args.data_file_json, args.run_id)
