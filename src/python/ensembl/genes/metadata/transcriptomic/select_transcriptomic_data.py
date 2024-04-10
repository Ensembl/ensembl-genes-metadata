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
import pymysql
import json
import random
from collections import Counter

#need to add path to config
with open(os.environ["ENSCODE"] + "/ensembl-genes-metadata/config.json", "r") as f:
    config = json.load(f)

def mysql_fetch_data(query, database, host, port, user, password):
    try:
        conn = pymysql.connect(
            host=host, user=user, passwd=password, port=port, database=database.strip()
        )

        cursor = conn.cursor()
        cursor.execute(query)
        info = cursor.fetchall()

    except pymysql.Error as err:
        print(err)

    cursor.close()
    conn.close()
    return info

def select_data(taxon_id, reads_mapped_cutoff, prioritise_tissues, max_num_runs):
    """Select the best data to align to the genome."""
    
    selected_runs = []
    #select the runs from the database that have passed QC and the percent_mapped reads is greater than reads_mapped_cutoff
    data_query = (
        "SELECT  run.run_accession, run.tissue, align.percent_mapped FROM run INNER JOIN align on run.run_id=align.run_id WHERE run.qc_status='qc_pass' AND align.percent_mapped>=" + str(reads_mapped_cutoff) + " and run.taxon_id=" + taxon_id + ";"
    )
    data_fetch = mysql_fetch_data(
        data_query,
        config["server_details"]["db_name"],
        config["server_details"]["db_host"],
        config["server_details"]["db_port"],
        config["server_details"]["db_user"],
        config["server_details"]["db_pass"],
    )
    run_dict = {}
    for tuple in data_fetch:
        run_dict[tuple[0]] = {"tissue":tuple[1],
                               "percent_mapped":tuple[2]}
    if prioritise_tissues:
        prioritised_tissues = config["tissue_types"]["prioritise"]

        tissue_counter = Counter()
        prioritised_runs = {}

        for key, value in run_dict.items():
            if value["tissue"] in prioritised_tissues:
                least_common_tissue, _ = tissue_counter.most_common()[-1] if tissue_counter else (None, None)
                if value["tissue"] == least_common_tissue or not tissue_counter[value["tissue"]]:
                    prioritised_runs[key] = value
                    tissue_counter[value["tissue"]] += 1
                    if len(prioritised_runs) == max_num_runs:
                        break
        #fill up remaining slots
        if len(prioritised_runs) < max_num_runs:
            remaining_runs = {k: v for k, v in run_dict.items() if k not in prioritised_runs}
            sample_size = min(max_num_runs - len(prioritised_runs), len(remaining_runs))
            additional_runs = random.sample(remaining_runs.items(), sample_size)
            prioritised_runs.update(additional_runs)

        selected_runs = list(prioritised_runs)[:max_num_runs]
            
    else:
        sorted_runs = sorted(run_dict,key=lambda x:run_dict[x]['percent_mapped'])
        selected_runs = sorted_runs[:max_num_runs]

    return selected_runs    
        
        
def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Parameters")
    parser.add_argument(
        "-t",
        "--taxon_id",
        required=True,
        help="Taxon ID"
    )

    parser.add_argument(
        "--reads_mapped_cutoff",
        type=int, default=50, help="The minimum allowed for percent_mapped reads."
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
        type=int, default=100, help="The maximum number of runs to be included in the output."
    )

    return parser.parse_args()

def main() -> None:
    """Entrypoint"""
    args = parse_args()
    runs_to_use=select_data(args.taxon_id, args.reads_mapped_cutoff, args.prioritise_tissues, args.max_num_runs)
    print(runs_to_use)
    
if __name__ == "__main__":
    main()
