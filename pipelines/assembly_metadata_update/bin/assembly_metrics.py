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

import argparse
import json
import logging
from typing import Dict, Any, List, Tuple
import re
import pymysql # type: ignore

def execute_query(query, db_params):
    conn = pymysql.connect(**db_params)
    cursor = conn.cursor()
    cursor.execute(query)
    result = cursor.fetchall()
    cursor.close()
    conn.close()
    return result

def execute_write(query: str, db_params: Dict[str, Any]) -> int:
    """
    Execute INSERT/UPDATE/DELETE and commit.
    Returns number of affected rows.
    """
    conn = pymysql.connect(**db_params)
    try:
        with conn.cursor() as cursor:
            affected = cursor.execute(query)
        conn.commit()
        return affected
    finally:
        conn.close()

def getting_metrics_registry(assembly_id, metadata_params):

    query_asm_metrics = f"SELECT metrics_name, metrics_value FROM assembly_metrics WHERE assembly_id = '{assembly_id}';"
    logging.info(query_asm_metrics)

    asm_metrics_registry_tuple  = execute_query(query_asm_metrics, metadata_params)

    asm_metrics_registry = {}
    for metric in asm_metrics_registry_tuple:
        asm_metrics_registry[metric[0]] = metric[1]

    return asm_metrics_registry

def normalise_value(metric_name: str, value: Any) -> Any:
    """
    Normalise values so comparisons don't fail just because of type/format
    differences (e.g., '13' vs 13, '37.0x' vs '37').
    """
    if value is None:
        return None

    # If it's already numeric, keep it (but unify ints/floats where sensible)
    if isinstance(value, (int, float)):
        return float(value) if isinstance(value, float) else int(value)

    # Strings: strip whitespace
    if isinstance(value, str):
        s = value.strip()

        # Special case: genome_coverage like "37.0x" -> 37.0
        if metric_name == "genome_coverage":
            # keep only the first numeric token
            m = re.search(r"[-+]?\d*\.?\d+", s)
            return float(m.group(0)) if m else s

        # Generic numeric parsing:
        # - if it looks like an int -> int
        # - if it looks like a float -> float
        # Otherwise keep as string
        if re.fullmatch(r"[-+]?\d+", s):
            return int(s)
        if re.fullmatch(r"[-+]?\d*\.\d+", s):
            return float(s)

        return s

    # Fallback: compare as-is
    return value

def sql_escape(value: Any) -> str:
    """
    Minimal SQL string escaping for single quotes.
    (If you're using a DB driver, prefer parameterised queries instead.)
    """
    if value is None:
        return "NULL"
    return str(value).replace("'", "''")

def generate_metric_upserts(
    accession: str,
    ncbi: Dict[str, Any],
    registry: Dict[str, Any],
    assembly_id: int,
    metadata_params: Dict[str, Any],
) -> List[str]:
    output_line_list: List[str] = []

    ncbi_keys = set(ncbi.keys())
    reg_keys = set(registry.keys())

    for metric_name in sorted(ncbi_keys - reg_keys):
        new_raw_value = ncbi[metric_name]
        new_value = normalise_value(metric_name, new_raw_value)

        insert_query = (
            f"INSERT INTO assembly_metrics (assembly_id, metrics_name, metrics_value) "
            f"VALUES ({assembly_id}, '{sql_escape(metric_name)}', '{sql_escape(new_value)}');"
        )
        logging.info(insert_query)
        affected = execute_write(insert_query, metadata_params)
        output_line = f"{accession}, asm_metrics, false, NA, {metric_name}:{new_value}"
        output_line_list.append(output_line)

    for metric_name in sorted(ncbi_keys & reg_keys):
        new_raw = ncbi[metric_name]
        old_raw = registry[metric_name]

        new_norm = normalise_value(metric_name, new_raw)
        old_norm = normalise_value(metric_name, old_raw)

        if new_norm != old_norm:
            update_query = (
                f"UPDATE assembly_metrics "
                f"SET metrics_value = '{sql_escape(new_norm)}' "
                f"WHERE assembly_id = {assembly_id} AND metrics_name = '{sql_escape(metric_name)}';"
            )
            logging.info(update_query)
            affected = execute_write(update_query, metadata_params)
            output_line = f"{accession}, asm_metrics, false, {metric_name}:{old_raw}, {metric_name}:{new_raw}"
            output_line_list.append(output_line)

    return output_line_list

def main():
    """ Module's entry point
    """

    logging.basicConfig(filename="update_assembly_metrics.log", level=logging.DEBUG, filemode='w',
                    format="%(asctime)s:%(levelname)s:%(message)s")

    parser = argparse.ArgumentParser(prog='update_assembly_metrics.py',
                                    description="Compare and update assembly metrics between NCBI data and Registry database.")

    parser.add_argument('--accession_json',
                        type=str,
                        required=True,
                        help='Path to JSON file with assembly metadata retrieved from NCBI API')
    parser.add_argument('--accession',
                        type=str,
                        required=True,
                        help='GCA accession to retrieve metadata')
    parser.add_argument('--metadata_params',
                        type=str,
                        required=True,
                        help='Database connection parameters for metadata database in JSON format')
    

    args = parser.parse_args()
    logging.info(args)

    accession = args.accession.strip()

    with open(args.accession_json, 'r') as json_file:
        data = json.load(json_file)
    
    with open(args.metadata_params, 'r') as params_file:
        metadata_params = json.load(params_file)

    # Getting assembly_id from registry to use in metrics queries
    query_assembly_id = f"SELECT assembly_id FROM assembly WHERE CONCAT(gca_chain, '.', gca_version) = '{accession}'"
    assembly_id = execute_query(query_assembly_id, metadata_params)[0][0]
    logging.info(f"Assembly {accession} has assembly_id: {assembly_id}")

    # Getting existing metrics from Registry
    asm_metrics_registry = getting_metrics_registry(assembly_id, metadata_params)

    # Getting metrics from NCBI report
    assembly_metrics_ncbi = data['reports'][0]['assembly_stats']
    
    output_line_list = generate_metric_upserts(
    accession=accession,
    ncbi=assembly_metrics_ncbi,
    registry=asm_metrics_registry,
    assembly_id=assembly_id,
    metadata_params=metadata_params)

    for output_line in output_line_list:
        print(output_line)

if __name__ == '__main__':
    main()
