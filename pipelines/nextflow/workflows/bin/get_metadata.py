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
"""Script to produce a json file per given run, containing all accessions, in the format:

{
    table : [
        metakey : metavalue,
        ...
        ]
},...

Returns:
    str : json element containing metadata for given run
"""

import argparse
import re
from typing import Dict
import json
import requests
from tenacity import retry, stop_after_attempt, wait_fixed


def request_data(run_accession: str, fields: list) -> str:
    """Make an HTTP request for the metadata of the given run_accession.

    Args:
        run_accession (str): Unique run accession.
        fields (list): List of fields of interest.

    Returns:
        str: TSV file containing metadata.
    """

    query = 'run_accession="' + run_accession + '"'

    data = {
        "result": "read_run",
        "query": query,
        "fields": ",".join(fields),
        "format": "tsv",
    }

    @retry(stop=stop_after_attempt(3), wait=wait_fixed(5))
    def make_request():
        try:
            response = requests.post("https://www.ebi.ac.uk/ena/portal/api/search", data=data, timeout=20)
            response.raise_for_status()
            return response.text
        except requests.RequestException as e:
            raise RuntimeError(f"Error fetching metadata: {e}") from e

    try:
        result = make_request()
        return result
        # Further processing of result can be done here
    except RuntimeError as e:
        print("Error retriving data")
        return str(e)


def json_parse(response: str, fields: list):
    """Parse response from HTTP request into a list of JSON strings.

    Args:
        response (str): TSV from the HTTP request.
        fields (list): List of fields of interest.

    Returns:
        list: List of JSON formatted dictionaries of metadata and values.
    """

    table = response.split("\n")
    meta_data = table[1].split("\t")
    output_data = {}
    for key, value in zip(fields, meta_data):
        output_data[key] = value

    table_run = {
        "run": {
            "taxon_id": output_data["tax_id"],
            "run_accession": output_data["run_accession"],
            "qc_status": "NOT_CHECKED" if len(output_data["fastq_ftp"].split(";")) >= 2 else "FILE_ISSUE",
            "sample_accession": output_data["sample_accession"],
            "study_accession": output_data["study_accession"],
            "read_type": output_data["library_strategy"],
            "platform": output_data["instrument_platform"],
            "paired": output_data["library_layout"] == "PAIRED",
            "experiment": "; ".join(
                value
                for value in [
                    re.sub(r"[!\"#$%&()*\+,\-\'.\/:;<=>?@\[\]^`{|}~]", "", output_data["experiment_alias"]),
                    re.sub(r"[!\"#$%&()*\+,\-\'.\/:;<=>?@\[\]^`{|}~]", "", output_data["experiment_title"]),
                ]
                if value is not None
            ).rstrip("; "),
            "run_description": re.sub(
                r"[!\"#$%&()*\+,\-\'.\/:;<=>?@\[\]^`{|}~]", "", output_data["description"]
            )[:250],
            "library_name": "; ".join(
                value
                for value in [output_data["library_source"], output_data["library_name"]]
                if value is not None
            ).rstrip("; "),
            "library_selection": output_data["library_selection"],
            "tissue": "; ".join(
                value
                for value in [
                    re.sub(r"[!\"#$%&()*\+,\-\'.\/:;<=>?@\[\]^`{|}~]", "", output_data["tissue_type"]),
                    re.sub(r"[!\"#$%&()*\+,\-\'.\/:;<=>?@\[\]^`{|}~]", "", output_data["tissue_lib"]),
                ]
                if value != ""
            ).rstrip("; "),
            "cell_line": re.sub(r"[!\"#$%&()*\+,\-\'.\/:;<=>?@\[\]^`{|}~]", "", output_data["cell_line"]),
            "cell_type": output_data["cell_type"],
            "strain": "; ".join(
                value
                for value in [
                    output_data["strain"],
                    output_data["cultivar"],
                    output_data["ecotype"],
                    output_data["isolate"],
                ]
                if value != ""
            ).rstrip("; "),
        }
    }

    table_study = {
        "study": {
            "study_accession": output_data["study_accession"],
            "center_name": output_data["center_name"],
        }
    }

    table_data_files: Dict[str, list[Dict[str, str]]] = {"data_files": []}

    file_url = output_data["fastq_ftp"].split(";")
    file_md5 = output_data["fastq_md5"].split(";")

    # Assuming the expected length is 3
    if len(file_url) == len(file_md5) == 3:
        file_entries = [f for f in file_url if "_1.fastq.gz" in f or "_2.fastq.gz" in f]
        md5_entries = [
            md5 for md5, file in zip(file_md5, file_url) if "_1.fastq.gz" in file or "_2.fastq.gz" in file
        ]
        # md5_entries=[]
        # for md5, file in zip(file_md5, file_url):
        #    if "_1.fastq.gz" in file or "_2.fastq.gz" in file:
        #       md5_entries.append(md5)
        file_url = file_entries
        file_md5 = md5_entries

    for url, md5 in zip(file_url, file_md5):
        extension_name = url.split("/")[-1]
        base_name = extension_name.split(".")[0]
        read = {"file_name": base_name, "file_url": url, "md5": md5}
        table_data_files["data_files"].append(read)

    json_run = json.dumps(table_run)
    json_study = json.dumps(table_study)
    json_data_files = json.dumps(table_data_files)

    with open("insert_into_run.json", "w") as file:
        file.write(json_run)
    with open("insert_into_study.json", "w") as file:
        file.write(json_study)
    with open("insert_into_data_file.json", "w") as file:
        file.write(json_data_files)


def main() -> None:
    """Module's entry-point."""

    fields = [
        "run_accession",
        "fastq_ftp",  # uses ftp, and contains both reads (if paired) separated by a semi colon (;)
        "fastq_md5",  # both separated by semicolon (;)
        "tax_id",
        "sample_accession",
        "study_accession",
        "library_strategy",  # read_type (it refers to wether a run is rna-seq, wgs, ...)
        "instrument_platform",  # platform
        "library_layout",  # paired
        "experiment_alias",
        "experiment_title",
        "description",
        "library_source",  # library_name (name is requested later, but it's often empty \
        # so we'll catch it still to fill tissue-related info if available but best to use \
        # source for this field)
        "library_selection",
        "tissue_type",
        "tissue_lib",
        "cell_line",
        "cell_type",
        "library_name",  # see library_source note
        "strain",
        "cultivar",
        "ecotype",
        "isolate",  #  different ways of expression "strains"
        # 'cage_protocol',       # this could be important to some farmed data, so added
        "center_name",
    ]

    parser = argparse.ArgumentParser(
        prog="get_metadata.py", description="Get metadata for a given run accession from the ENA site"
    )
    parser.add_argument(
        "--run",
        "-r",
        default="SRR10059726",
        help="Run accession to grab metadata from",
    )
    args = parser.parse_args()

    data = request_data(args.run, fields)
    json_parse(data, fields)


if __name__ == "__main__":
    main()
