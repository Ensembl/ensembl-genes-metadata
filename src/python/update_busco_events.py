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

"""
Updates the assembly_events table with the BUSCO status/priority of each assembly, ahead of
deciding which ones need a BUSCO run.

Two independent passes:

1. mark_existing_busco_done: assemblies that already have a BUSCO metric (assembly_metrics,
   metrics_name='assembly.busco') but no assembly_events row yet are marked
   genome_busco.status=done, genome_busco.priority=high.
2. insert_busco_candidates: assemblies with no assembly_events row yet that are current,
   haploid, chromosome/complete-genome level, have a genus-level taxon and a
   total_sequence_length metric are marked genome_busco.status=pending, with
   genome_busco.priority set from a tier based on genome size and transcriptomic read
   coverage (from the ENA transcriptomics cache API): large_genome, high, medium, or low.

Both passes are read-only unless --execute is passed; without it, the insert queries that
would have run are logged but not sent to the database. Nothing is printed to stdout/stderr --
everything (including every query and the run summary) goes to update_busco_events.log.

Usage
-----

    python update_busco_events.py --metadata-params '{"host": "...", "port": 4527,
        "user": "...", "password": "...", "database": "gb_assembly_metadata"}' [--execute]
"""

import argparse
import json
import logging
from typing import Any, Dict

from gb_metadata.db_utils import execute_query, execute_write
from gb_metadata.utils import connection_api

TRANSCRIPTOMICS_URI = "http://genebuild-metadata.ebi.ac.uk:8000/api/transcriptomics/ena-cache"

STATUS_EVENT = "genome_busco.status"
PRIORITY_EVENT = "genome_busco.priority"

LARGE_GENOME_THRESHOLD_BP = 8_000_000_000
HIGH_PRIORITY_MIN_READS = 10

EXISTING_BUSCO_QUERY = """
select asm.assembly_id, CONCAT(asm.gca_chain, '.', asm.gca_version)
from assembly asm
JOIN assembly_metrics am ON am.assembly_id = asm.assembly_id
LEFT JOIN assembly_events ae ON asm.assembly_id = ae.assembly_id
where am.metrics_name = 'assembly.busco'
AND ae.assembly_id IS NULL;
"""

CANDIDATE_ASSEMBLIES_QUERY = """
select DISTINCT asm.assembly_id, CONCAT(asm.gca_chain, '.', asm.gca_version),
    asm.lowest_taxon_id, t.taxon_class_id, am.metrics_value
from assembly asm
LEFT JOIN assembly_events ae ON asm.assembly_id = ae.assembly_id
JOIN assembly_metrics am ON asm.assembly_id = am.assembly_id
JOIN bioproject b ON asm.assembly_id = b.assembly_id
INNER JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
INNER JOIN taxonomy t ON asm.lowest_taxon_id = t.lowest_taxon_id
WHERE ae.assembly_id is NULL
AND asm.is_current = 'current'
AND asm_type = 'haploid'
AND asm_level IN ('chromosome', 'Complete genome')
AND t.taxon_class = 'genus'
AND am.metrics_name = 'total_sequence_length'
;
"""


def insert_busco_event(
    assembly_id: int, status: str, priority: str, db_params: Dict[str, Any], execute: bool
) -> None:
    """Insert/no-op-update the genome_busco.status and genome_busco.priority events for one assembly."""
    query = f"""
    insert into assembly_events (assembly_id, event, status)
    values ({assembly_id}, '{STATUS_EVENT}', '{status}'),
        ({assembly_id}, '{PRIORITY_EVENT}', '{priority}')
    ON DUPLICATE KEY UPDATE assembly_id = {assembly_id};
    """

    if not execute:
        logging.info("Execution skipped for assembly_id %s:\n%s", assembly_id, query)
        return

    execute_write(query, db_params)


def mark_existing_busco_done(db_params: Dict[str, Any], execute: bool) -> int:
    """Mark assemblies that already have a BUSCO metric but no assembly_events row as done/high."""
    rows = execute_query(EXISTING_BUSCO_QUERY, db_params)
    logging.info("Total assemblies with BUSCO metrics fetched: %s", len(rows))

    for assembly_id, gca in rows:
        logging.info("Marking existing BUSCO result as done: %s (assembly_id=%s)", gca, assembly_id)
        insert_busco_event(assembly_id, "done", "high", db_params, execute)

    return len(rows)


def classify_priority(genome_size: int, num_reads: int, num_reads_genus: int) -> str:
    """Return the BUSCO priority tier for one candidate assembly."""
    if genome_size > LARGE_GENOME_THRESHOLD_BP:
        return "large_genome"
    if num_reads > HIGH_PRIORITY_MIN_READS:
        return "high"
    if num_reads <= HIGH_PRIORITY_MIN_READS and num_reads_genus > 0:
        return "medium"
    return "low"


def fetch_transcriptomics_cache() -> Dict[str, Any]:
    """Fetch the ENA transcriptomics read-count cache used to prioritize candidates."""
    response = connection_api(TRANSCRIPTOMICS_URI)
    return response.json()


def insert_busco_candidates(db_params: Dict[str, Any], execute: bool) -> Dict[str, int]:
    """Classify and mark newly-eligible assemblies as pending, with a priority tier."""
    rows = execute_query(CANDIDATE_ASSEMBLIES_QUERY, db_params)
    logging.info("Total candidate assemblies fetched: %s", len(rows))

    transcriptomics_cache = fetch_transcriptomics_cache()

    tier_counts = {"large_genome": 0, "high": 0, "medium": 0, "low": 0}
    for assembly_id, gca, taxon_id, genus_taxon_id, genome_size in rows:
        num_reads = (
            transcriptomics_cache.get(f"{taxon_id}:True", {})
            .get("data", {})
            .get("short_read_paired_end_illumina", 0)
        )
        num_reads_genus = (
            transcriptomics_cache.get(f"{genus_taxon_id}:True", {})
            .get("data", {})
            .get("short_read_paired_end_illumina", 0)
        )

        priority = classify_priority(int(genome_size), num_reads, num_reads_genus)
        tier_counts[priority] += 1

        logging.info(
            "Candidate %s (assembly_id=%s) classified as %s", gca, assembly_id, priority
        )
        insert_busco_event(assembly_id, "pending", priority, db_params, execute)

    return tier_counts


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        prog="update_busco_events.py",
        description="Update assembly_events with BUSCO status/priority ahead of deciding what needs running.",
    )
    parser.add_argument(
        "--metadata-params",
        type=json.loads,
        required=True,
        help="JSON string with the metadata DB connection parameters "
        "(host, port, user, password, database).",
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="If this option is added, the insert queries are actually executed. "
        "Without it, queries are logged but not run.",
    )
    return parser.parse_args()


def main() -> None:
    """Module entry point."""
    logging.basicConfig(
        filename="update_busco_events.log",
        level=logging.DEBUG,
        filemode="w",
        format="%(asctime)s:%(levelname)s:%(message)s",
    )

    args = parse_args()
    db_params = args.metadata_params

    logging.info("Execute mode: %s", args.execute)

    done_count = mark_existing_busco_done(db_params, args.execute)
    tier_counts = insert_busco_candidates(db_params, args.execute)

    logging.info("Marked as already done: %s", done_count)
    logging.info("Large genome: %s", tier_counts["large_genome"])
    logging.info("High priority: %s", tier_counts["high"])
    logging.info("Medium priority: %s", tier_counts["medium"])
    logging.info("Low priority: %s", tier_counts["low"])
    logging.info("Total candidates: %s", sum(tier_counts.values()))


if __name__ == "__main__":
    main()
