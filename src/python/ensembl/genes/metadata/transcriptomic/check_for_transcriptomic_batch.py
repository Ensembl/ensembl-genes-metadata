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
from pathlib import Path
from typing import List
import argparse
import requests


def ena_rest_api(query: str) -> List:
    """Retrieve list of run accession from ENA

    Args:
        query (str): query 

    Returns:
        List: list of run accessions per taxon_id
    """
    search_url = f"https://www.ebi.ac.uk/ena/portal/api/search?display=report&query={query}&domain=read&result=read_run&fields=run_accession"  # pylint: disable=line-too-long

    search_result = requests.get(search_url, timeout=20)
    results = search_result.text.strip().split("\n")[1:]

    return results


def check_data_from_ena(  # pylint: disable=too-many-locals
    taxon_id: str, batching_option: bool, batch_size: int, output_dir: Path
) -> None:
    """Get run list of run accession from Ena
    Only for short read if batching option is enabled, the list of run accession
    is splitted in multiple list according to the batching size and stored in txt files.

    Args:
        batching_option (bool): enable the batching for short reads only
        batch_size (int): according to the size the list of run accession is splitted and 
        saved in different files
        output_dir (Path): output dir
    """

    query = f"tax_eq({taxon_id})"

    query_short_paired = (
        query
        + " AND instrument_platform=ILLUMINA AND library_layout=PAIRED AND library_source=TRANSCRIPTOMIC \
            AND first_created>=2019-01-01"
    )
    query_short_single = (
        query
        + " AND instrument_platform=ILLUMINA AND library_layout=SINGLE AND library_source=TRANSCRIPTOMIC \
            AND first_created>=2019-01-01"
    )
    query_pacbio = (
        query
        + " AND instrument_platform=PACBIO_SMRT AND library_source=TRANSCRIPTOMIC \
            AND first_created>=2019-01-01"
    )
    query_onp = (
        query
        + " AND instrument_platform=OXFORD_NANOPORE AND library_source=TRANSCRIPTOMIC \
            AND first_created>=2019-01-01"
    )

    short_paired_runs = ena_rest_api(query_short_paired)
    short_single_runs = ena_rest_api(query_short_single)
    pacbio_read_runs = ena_rest_api(query_pacbio)
    onp_read_runs = ena_rest_api(query_onp)

    print(
        f"{taxon_id};Short-read paired-end illumina;{len(short_paired_runs)};Short-read single-end illumina;\
        {len(short_single_runs)};Long-read PacBio;{len(pacbio_read_runs)};Long_read ONP;{len(onp_read_runs)}"
    )

    # ONLY FORE PAIRED SHORT READS
    # If batching is enabled, split the results based on batch_size
    if batching_option and batch_size > 0:
        batches = [
            short_paired_runs[i : i + batch_size] for i in range(0, len(short_paired_runs), batch_size)
        ]
        output_dir = output_dir / taxon_id / "batch"
        # Create the directory structure if it doesn't exist
        output_dir.mkdir(parents=True, exist_ok=True)
        # Save each batch into a separate file
        for idx, batch in enumerate(batches):
            batch_file = output_dir / f"batch_{idx + 1}.txt"
            with open(batch_file, "w") as f:
                f.write("\n".join(batch))


class InputSchema(argparse.ArgumentParser):
    """Input arguments"""

    def __init__(self):
        super().__init__()

        self.add_argument("-t", "--taxon_id", type=str, required=False, help="Taxon id")
        self.add_argument("--output_dir", required=False, help="Output directory path")
        self.add_argument("--batching_option", type=bool, required=False, help="Batch run accession")
        self.add_argument("--batch_size", type=int, required=False, help="Batch size")
        self.add_argument(
            "-f", "--file", type=str, required=False, help="Path to the file containing a list of taxon ids"
        )


def main() -> None:
    """Entrypoint"""
    parser = InputSchema()
    args = parser.parse_args()
    if args.file:
        with open(args.file, "r") as input_file:
            taxon_ids = input_file.read().splitlines()
            for taxon_id in taxon_ids:
                check_data_from_ena(taxon_id, args.batching_option, args.batch_size, Path(args.output_dir))
    else:
        check_data_from_ena(args.taxon_id, args.batching_option, args.batch_size, Path(args.output_dir))


if __name__ == "__main__":
    main()
