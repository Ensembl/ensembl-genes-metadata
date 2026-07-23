#!/usr/bin/env python3

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


import logging
import argparse
import requests  # type: ignore
from datetime import datetime
import json
from tenacity import retry, stop_after_attempt, wait_random


@retry(stop=stop_after_attempt(5), wait=wait_random(min=1, max=10))
def connection_ncbi(uri: str) -> requests.Response:
    """Connects to NCBI API and retrieves data for a given URI
    Args:
        uri (str): URI to request data from NCBI API
    Returns:
        requests.Response: response object from NCBI API
    """
    response = requests.get(uri, timeout=60)
    response.raise_for_status()
    return response


def get_ncbi_json(accession, ncbi_url, attempt_update):
    uri = f"{ncbi_url}/genome/accession/{accession}/dataset_report?filters.exclude_atypical=false&filters.assembly_version=all_assemblies"
    response = connection_ncbi(uri)
    data = response.json()

    if data == {}:
        logging.info(f"Assembly {accession} is deleted from NCBI")
        query_update_status = f"UPDATE assembly a SET is_current = 'suppressed' WHERE CONCAT(a.gca_chain, '.', a.gca_version) = '{accession}';"
        print(query_update_status)
        attempt_update = False

    return data, attempt_update


def main():
    parser = argparse.ArgumentParser(description="Fetch assembly metadata from NCBI API")
    parser.add_argument(
        "--accession",
        type=str,
        required=True,
        help="GCA accession to retrieve metadata",
    )
    parser.add_argument(
        "--ncbi_url",
        type=str,
        required=True,
        help="NCBI API base URL",
    )

    args = parser.parse_args()
    logging.info(args)

    accession = args.accession.strip()
    ncbi_url = args.ncbi_url

    attempt_update = True
    data, attempt_update = get_ncbi_json(accession, ncbi_url, attempt_update)

    with open(f"{accession}_metadata.json", "w") as json_file:
        json.dump(data, json_file, indent=4)
    logging.info(f"Metadata for assembly {accession} saved to {accession}_metadata.json")

    if attempt_update == True:
        print("true")
    elif attempt_update == False:
        print("false")


if __name__ == "__main__":
    main()
