# See the NOTICE file distributed with this work for additional information #pylint: disable=missing-module-docstring
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Select the best subset of short-read transcriptomic data to align to the genome"""

import argparse
from collections import Counter
import csv
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import json
import random
import pymysql


def mysql_fetch_data(
    query: str, database: str, host: str, port: int, user: str, password: str
) -> Optional[List[Tuple]]:
    """
    Fetch data from MySQL database based on the provided query.

    Args:
        query: SQL query to be executed.
        database: Name of the database.
        host: Database host.
        port: Port number for the connection.
        user: Username for the database connection.
        password: Password for the database connection.

    Returns:
        A list of tuples representing the rows fetched from the database, or None if an error occurs.
    """
    try:
        conn = pymysql.connect(host=host, user=user, passwd=password, port=port, database=database.strip())

        with conn.cursor() as cursor:
            cursor.execute(query)
            info: List[Tuple] = cursor.fetchall()


    except pymysql.Error as err:
        print(err)

    finally:
        conn.close()
    return info

def select_data(
    taxon_id: str,
    reads_mapped_cutoff: float,
    prioritise_tissues: bool,
    max_num_runs: int,
    config: Dict[str, Any],
) -> Dict[str, Dict[str, Any]]:
    """
    Select the best data to align to the genome.

    Args:
        taxon_id: The taxon ID to filter by.
        reads_mapped_cutoff: Minimum percent of reads mapped to pass.
        prioritise_tissues: Whether to prioritise certain tissue types.
        max_num_runs: Maximum number of runs to select.
        config: Configuration details including database connection info and tissue prioritisation.

    Returns:
        A list of selected run accession IDs.
    """

    selected_runs: Dict[str, Dict[str, Any]] = {}
    # select the runs from the database that have passed QC and the percent_mapped reads
    # is greater than reads_mapped_cutoff
    data_query = (
        f"SELECT  run.run_accession, run.tissue, align.uniquely_mapped_reads_percentage FROM "
        f"run INNER JOIN align ON run.run_id=align.run_id WHERE run.qc_status='ALIGNED' AND "
        f"align.uniquely_mapped_reads_percentage>={reads_mapped_cutoff} AND run.taxon_id={taxon_id};"
    )
    data_fetch = mysql_fetch_data(
        data_query,
        config["server_details"]["db_name"],
        config["server_details"]["db_host"],
        config["server_details"]["db_port"],
        config["server_details"]["db_user"],
        config["server_details"]["db_pass"],
    )
    if not data_fetch:
        print("No data fetched or an error occurred.")
        return {}

    # Prepare a dictionary to store the fetched run data
    run_dict: Dict[str, Dict[str, Any]] = {
        row[0]: {"tissue": row[1], "percent_mapped": row[2]} for row in data_fetch
    }
    if prioritise_tissues:
        prioritised_tissues: List[str] = config["tissue_types"]["prioritise"]
        tissue_counter: Counter = Counter()
        prioritised_runs: Dict[str, Dict[str, Any]] = {}
        for key, value in run_dict.items():
            if value["tissue"] in prioritised_tissues:
                least_common_tissue, _ = tissue_counter.most_common()[-1] if tissue_counter else (None, None)
                if value["tissue"] == least_common_tissue or not tissue_counter[value["tissue"]]:
                    prioritised_runs[key] = value
                    tissue_counter[value["tissue"]] += 1
                    if len(prioritised_runs) == max_num_runs:
                        break
        # Fill up remaining slots
        if len(prioritised_runs) < max_num_runs:
            remaining_runs: Dict[str, Dict[str, Any]] = {
                k: v for k, v in run_dict.items() if k not in prioritised_runs
            }
            sample_size: int = min(max_num_runs - len(prioritised_runs), len(remaining_runs))
            additional_runs = random.sample(list(remaining_runs.items()), sample_size)
            prioritised_runs.update(dict(additional_runs))

        #selected_runs = prioritised_runs.keys()[:max_num_runs]
        selected_runs = {k: prioritised_runs[k] for k in list(prioritised_runs.keys())[:max_num_runs]}
    else:
        #sorted_runs: List[str] = sorted(run_dict, key=lambda x: run_dict[x]["percent_mapped"], reverse=True)
        #selected_runs = sorted_runs[:max_num_runs]
        sorted_runs = dict(sorted(run_dict.items(), key=lambda item: item[1]["percent_mapped"], reverse=True))
        selected_runs = {k: sorted_runs[k] for k in list(sorted_runs.keys())[:max_num_runs]}  # Slice to get top runs

    print(list(selected_runs))
    return selected_runs


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Parameters")
    parser.add_argument("-t", "--taxon_id", required=True, help="Taxon ID")
    parser.add_argument("--output_dir",  help="Output directory")
    parser.add_argument(
        "--reads_mapped_cutoff", type=int, default=50, help="The minimum allowed for percent_mapped reads."
    )

    parser.add_argument(
        "-p",
        "--prioritise_tissues",
        action="store_true",
        default=False,
        help="Prioritise runs from tissues that are specified in config.json.",
    )

    parser.add_argument(
        "--max_num_runs",
        type=int,
        default=100,
        help="The maximum number of runs to be included in the output.",
    )

    return parser.parse_args()


def main() -> None:
    """Entrypoint"""
    args = parse_args()
    # Get the directory where the current script is located
    script_dir = Path(__file__).parent.resolve()

    # Define the path to 'ensembl-genes-metadata/config.json' relative to your script's location
    config_path = script_dir.parents[3] / "ensembl-genes-metadata" / "conf" / "config.json"

    # Open the config file
    with open(config_path, "r") as f:
        config = json.load(f)

    runs_to_use = select_data(
        args.taxon_id, args.reads_mapped_cutoff, args.prioritise_tissues, args.max_num_runs, config
    )

    # Ensure output directory exists
    output_csv = Path(args.output_dir) / args.taxon_id / "run_accession.csv"
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    # Write the selected runs to a CSV file
    with open(output_csv, "w", newline="") as csvfile:
        fieldnames = ["Run Accession", "Tissue", "Percent Mapped"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for run_accession, run_data in runs_to_use.items():
            writer.writerow({
                "Run Accession": run_accession,
                "Tissue": run_data["tissue"],
                "Percent Mapped": run_data["percent_mapped"]
            })

if __name__ == "__main__":
    main()
