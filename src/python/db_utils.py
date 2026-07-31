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

"""Shared MySQL connection helpers for the assembly_metadata and assembly_metadata_update
pipeline scripts.
"""

import logging

from typing import Any, Dict, Tuple

import pymysql  # type: ignore


def execute_query(query: str, db_params: Dict[str, Any]) -> Tuple[Tuple[Any, ...], ...]:
    """Run a query and return every row."""
    logging.info("QUERY: %s", query)
    with pymysql.connect(**db_params) as conn:
        with conn.cursor() as cursor:
            cursor.execute(query)
            return cursor.fetchall()


def execute_write(query: str, db_params: Dict[str, Any]) -> int:
    """Execute an INSERT/UPDATE/DELETE and commit. Returns the number of affected rows."""
    logging.info("QUERY: %s", query)
    with pymysql.connect(**db_params) as conn:
        with conn.cursor() as cursor:
            affected = cursor.execute(query)
        conn.commit()
        return affected


def fetch_one_row(query: str, db_params: Dict[str, Any], context: str) -> Tuple[Any, ...]:
    """Run query and return its single result row, or raise if the row count isn't exactly 1."""
    results = execute_query(query, db_params)
    if len(results) != 1:
        raise ValueError(f"Expected one result for {context}, but got {len(results)} results")
    return results[0]
