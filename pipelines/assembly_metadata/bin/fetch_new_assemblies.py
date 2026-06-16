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


"""The module will connect to NCBI API and the assembly metadata database to determine the list of
GCA accessions to register in the assembly metadata database.

Args:
    taxon (int): Taxon ID, by default it is 2759 (Eukaryote domain)
    date_update (str): date used to retrieve new assemblies from NCBI (optional)
    metadata (str): path to the metadata database parameters in json format
    ncbi (str): path to the NCBI API parameters in json format
    ncbi_url (str): NCBI API URL
Returns:
    stdout: prints to the standard output a list of GCA accessions
    path: assemblies_to_register.txt: file containing the list of GCA accessions to register in the assembly metadata database
"""

import requests # type: ignore
import json
from datetime import datetime
import pymysql # type: ignore
import argparse
import logging
from tenacity import retry, stop_after_attempt, wait_random # type: ignore
from typing import Dict, List, Any, Tuple, Optional

DEFAULT_ASSEMBLY_DATE = "01/01/2019"

def set_date(ncbi_params: Dict[str, Any], date_update: Optional[str] = None) -> Tuple[Dict[str, str], str]:
    """
    Set a date filter for retrieving assemblies from NCBI.

    Args:
        ncbi_params (dict): NCBI API parameters to update in place.
        date_update (str, optional): Release date filter (MM/DD/YYYY). 
                                     Defaults to DEFAULT_ASSEMBLY_DATE if not provided.

    Returns:
        tuple:
            - dict: Updated NCBI API parameters with the date filter applied.
            - str: The date used (either provided or default).
    """
    logging.info('Checking date to retrieve assemblies')

    effective_date = date_update or DEFAULT_ASSEMBLY_DATE
    ncbi_params['filters.first_release_date'] = effective_date

    if date_update:
        logging.info("Using provided date %r to retrieve assemblies", effective_date)
    else:
        logging.info("No date provided. Using default date %r", effective_date)

    return ncbi_params, effective_date

@retry(stop=stop_after_attempt(3), wait=wait_random(min=1, max=3))
def connection_ncbi(uri: str, params: Dict[str, str]) -> requests.Response:
    """Make a GET request to the NCBI API.

    Args:
        uri (str): The URI for the NCBI API endpoint.
        params (dict): Parameters for the GET request.

    Returns:
        requests.Response: The response object from the GET request.
    """
    logging.debug("Requesting NCBI URL: %s with params: %s", uri, params)
    response = requests.get(uri, params=params)
    response.raise_for_status()
    return response

def fetch_gca_list(taxon: int, ncbi_params: Dict[str, str], ncbi_url: str) -> set[str]:
    """
    Fetch a list of GCA accessions from NCBI API based on the taxon ID.

    Args:
        taxon (int): Taxon ID, by default it is 2759 (Eukaryote domain)
        ncbi_params (dict): NCBI API's parameters
        ncbi_url (str): NCBI API base URL

    Returns:
        set: set of GCA accessions
    """
    gca_list: list[str] = []
    page_token=None
    uri = f"{ncbi_url}/genome/taxon/{str(taxon)}/dataset_report"
    next_page = True

    while next_page:
        response = connection_ncbi(uri, ncbi_params)
        logging.info(f"URL: {response.url}")
        data = response.json()

        if 'reports' in data:
            gca_list.extend(report['accession'] for report in data['reports'])

            page_token = data.get('next_page_token')
            if page_token:
                ncbi_params['page_token'] = page_token
            elif page_token is None:
                next_page = False
            else:
                raise ValueError("Page token is not being set correctly")

        else:
            logging.info("No assemblies found. Setting next_page as false")
            next_page = False

    return set(gca_list)

def build_db_query(release_date: str) -> str:
    """Build the MySQL query to retrieve assemblies from the metadata database.

    Args:
        release_date (str): date to be used to retrieve assemblies (mm/dd/yyyy)

    Returns:
        str: SQL query string
    """
    release_date_sql = datetime.strptime(release_date, '%m/%d/%Y').strftime('%Y-%m-%d')
    return f"""
        SELECT CONCAT(gca_chain, '.', gca_version)
        FROM assembly
        WHERE release_date >= DATE('{release_date_sql}')
    """

def fetch_records_db(db_params: Dict[str, Any], query: str) -> List[str]:
    """Fetch assemblies that have been registered after the last update.

    Args:
        db_params (dict): database connection parameters
        query (str): mysql query to retrieve data

    Returns:
        list: list of GCA accessions recorded after the last update date
    """
    with pymysql.connect(**db_params) as conn:
        with conn.cursor() as cur:
            cur.execute(query)
            output = cur.fetchall()
            reg_gca = [row[0] for row in output]

    return reg_gca

def get_gca_to_register(gca_list: set, query: str, metadata_params: Dict[str, Any]) -> List[str]:
    """Return GCA accessions from NCBI not yet present in the metadata database.

    Args:
        gca_list (set): GCA accessions retrieved from NCBI API
        query (str): SQL query to fetch already-registered accessions
        metadata_params (dict): metadata database connection parameters

    Returns:
        list: accessions to register
    """
    logging.info('Getting assemblies from assembly metadata database')
    records_metadata = fetch_records_db(metadata_params, query=query)
    return list(gca_list - set(records_metadata))

def main():
    """module's entry-point
    """

    logging.basicConfig(filename="fetch_new_assemblies.log", level=logging.DEBUG, filemode='w',
                    format="%(asctime)s:%(levelname)s:%(message)s")

    parser = argparse.ArgumentParser(prog='fetch_new_assemblies.py',
                                    description='Identify new assemblies to register in the assembly metadata database. Assemblies are retrieved from NCBI API and compared with the assembly metadata database.')

    parser.add_argument('--taxon',
                        default=2759,
                        type=int,
                        help='Valid Taxon id: Eukaryota - 2759')
    parser.add_argument('--date_update',
                        type=str,
                        help="Last update date")
    parser.add_argument('--metadata',
                        type=str,
                        required=True,
                        help="Path to the metadata database params in json format")
    parser.add_argument('--ncbi',
                        type=str,
                        required=True,
                        help="Path to the NCBI API params in json format")
    parser.add_argument('--ncbi_url',
                        type=str,
                        required=True,
                        help="NCBI API URL")

    args = parser.parse_args()
    logging.info(args)

    with open(args.metadata, 'r') as file:
        metadata_params = json.load(file)

    with open(args.ncbi, 'r') as file:
        ncbi_params = json.load(file)

    if args.date_update:
        try:
            datetime.strptime(args.date_update, '%m/%d/%Y')
        except ValueError:
            raise ValueError("Please enter a valid date format (mm/dd/yyyy)")
    else:
        logging.info(f"Default date ({DEFAULT_ASSEMBLY_DATE}) will be used to retrieve assemblies")

    ncbi_params, release_date = set_date(ncbi_params, args.date_update)
    gca_list = fetch_gca_list(args.taxon, ncbi_params, args.ncbi_url)

    if len(gca_list) > 0:
        query = build_db_query(release_date)
        accessions_to_register = get_gca_to_register(gca_list, query, metadata_params)

        with open("assemblies_to_register.txt", 'w') as file:
            for accession in accessions_to_register:
                print(accession)
                file.write(accession + '\n')

        logging.info(f'Accessions to register: {len(accessions_to_register)}. Please note that some assemblies might belong to unspecified species')

    else:
        logging.info(f'No assemblies found since {release_date}')

if __name__ == '__main__':
    main()
