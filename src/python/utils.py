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
""" Utility functions.
"""
import logging

import requests  # type: ignore
from tenacity import retry, stop_after_attempt, wait_random  # type: ignore


@retry(stop=stop_after_attempt(10), wait=wait_random(min=1, max=20))
def connection_api(uri: str) -> requests.Response:
    """Connect to the API and return the HTTP response, retrying transient failures."""
    logging.info("Connecting to API: %s", uri)
    response = requests.get(uri, timeout=30)
    response.raise_for_status()
    return response
