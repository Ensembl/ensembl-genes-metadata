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
"""Check the availability for short and long read data from ENA website given a taxon id"""
import os.path
from pathlib import Path
from typing import List
import argparse
import requests

def ena_rest_api(query:str, taxon_id: int, batching_option:bool, batch_size:int , output_dir_path:Path)->int:
    search_url = f"https://www.ebi.ac.uk/ena/portal/api/search?display=report&query={query}&domain=read&result=read_run&fields=run_accession"  # pylint: disable=line-too-long
    #def searchUrl = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=${taxonQuery}%20AND%20instrument_platform=ILLUMINA%20AND%20library_layout=PAIRED%20AND%20library_source=TRANSCRIPTOMIC%20AND%20
    #first_created%3E=${lastCheckedDate.trim()}&domain=read&fields=run_accession"

    search_result = requests.get(search_url)
    results = search_result.text.strip().split("\n")[1:]
     # If batching is enabled, split the results based on batch_size
    if batching_option and batch_size > 0:
        batches = [results[i:i + batch_size] for i in range(0, len(results), batch_size)]
        output_dir = output_dir_path / taxon_id / "batch"
        # Create the directory structure if it doesn't exist
        output_dir.mkdir(parents=True, exist_ok=True)
        # Save each batch into a separate file
        for idx, batch in enumerate(batches):
            batch_file = output_dir / f"batch_{idx + 1}.txt"
            with open(batch_file, 'w') as f:
                f.write("\n".join(batch))
    return len(results)

def check_data_from_ena(# pylint: disable=too-many-locals
        taxon_id: int,
        tree: bool,
        batching_option:bool, batch_size:int, output_dir:Path
) -> None:
    """Query ENA API to get short or long read data"""

    if tree:
        query = f"tax_tree({taxon_id})"
    else:
        query = f"tax_eq({taxon_id})"
    
    query_short_paired=query+f" AND instrument_platform=ILLUMINA AND library_layout=PAIRED AND library_source=TRANSCRIPTOMIC AND first_created>=2019-01-01"
    query_short_single=query+f" AND instrument_platform=ILLUMINA AND library_layout=SINGLE AND library_source=TRANSCRIPTOMIC AND first_created>=2019-01-01"
    query_pacbio=query+f" AND instrument_platform=PACBIO_SMRT AND library_source=TRANSCRIPTOMIC AND first_created>=2019-01-01"
    query_onp=query+f" AND instrument_platform=OXFORD_NANOPORE AND library_source=TRANSCRIPTOMIC AND first_created>=2019-01-01"
    
    short_paired_runs=ena_rest_api(query_short_paired, taxon_id, batching_option, batch_size, output_dir)
    short_single_runs=ena_rest_api(query_short_single, taxon_id, batching_option, batch_size, output_dir)
    pacbio_read_runs =ena_rest_api(query_pacbio, taxon_id, batching_option, batch_size, output_dir)
    onp_read_runs=ena_rest_api(query_onp, taxon_id, batching_option, batch_size, output_dir)

    print(f"{taxon_id};Short-read paired-end illumina;{short_paired_runs};Short-read single-end illumina;{short_single_runs};Long-read PacBio;{pacbio_read_runs};Long_read ONP;{onp_read_runs}")
    #print (text.BOLD+f"Short-read paired-end illumina data available! "+text.END+f"Found {short_paired_runs} runs.")
    #print (text.BOLD+f"Short-read single-end illumina data available! "+text.END+f"Found {short_single_runs} runs.")
    #print (text.BOLD+f"Long-read PacBio data available! "+text.END+f"Found {pacbio_read_runs} runs.")
    #print (text.BOLD+f"Long_read ONP data available! "+text.END+f"Found {onp_read_runs} runs.")
        
class text:
    """formatting set"""

    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
    END = "\033[0m"


class InputSchema(argparse.ArgumentParser):
    """Input arguments"""
    def __init__(self):
        super().__init__()

        self.add_argument(
            "-t", "--taxon_id", type=str, required=False, help="Taxon id"
        )

        self.add_argument(
            "--tree", action='store_true', required=False, help="Turn on the 'Include subordinate taxa' option in your query to ENA"
        )
        self.add_argument(
            "--output_dir",   required=False, help="Output directory path"
        )
        self.add_argument(
            "--batching_option", type=bool, required=False, help="Batch run accession"
        )
        self.add_argument(
            "--batch_size", type=int, required=False, help="Batch size"
        )
        self.add_argument(
            "-f", "--file", type=str, required=False, help="Path to the file containing a list of taxon ids"
        )

def main() -> None:
    """Entrypoint"""
    parser=InputSchema()
    args = parser.parse_args()
    if args.file:
        with open(args.file, 'r') as input_file:
            taxon_ids=input_file.read().splitlines()
            for taxon_id in taxon_ids:
                check_data_from_ena(taxon_id, args.tree, args.batching_option, args.batch_size, Path(args.output_dir))
    else:
        check_data_from_ena(args.taxon_id, args.tree, args.batching_option, args.batch_size, Path(args.output_dir))

if __name__ == "__main__":
    main()

