# See the NOTICE file distributed with this work for additional information
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
"""Utilities for parsing BUSCO summary strings."""

import re
from typing import Dict, Optional, Union

BuscoValue = Optional[Union[float, int]]

_BUSCO_COMPONENT_PATTERN = re.compile(r"(?P<component>[CSDFMn]):(?P<value>\d+(?:\.\d+)?)%?")

BUSCO_COMPONENTS = (
    "complete",
    "single_copy",
    "duplicated",
    "fragmented",
    "missing",
    "n_genes",
)

_COMPONENT_NAMES = {
    "C": "complete",
    "S": "single_copy",
    "D": "duplicated",
    "F": "fragmented",
    "M": "missing",
    "n": "n_genes",
}


def parse_busco_string(busco_string: Optional[str]) -> Dict[str, BuscoValue]:
    """Parse standard BUSCO components without discarding composition fields.

    Args:
        busco_string: BUSCO summary such as
            ``C:98.2%[S:97.5%,D:0.7%],F:0.8%,M:1.0%,n:255``.

    Returns:
        A dictionary containing ``complete``, ``single_copy``, ``duplicated``,
        ``fragmented``, ``missing``, and ``n_genes``. Missing or malformed
        components are returned as ``None``. Percentages are floats and ``n``
        is an integer.
    """
    parsed: Dict[str, BuscoValue] = {component: None for component in BUSCO_COMPONENTS}

    if not isinstance(busco_string, str):
        return parsed

    for match in _BUSCO_COMPONENT_PATTERN.finditer(busco_string):
        component = _COMPONENT_NAMES[match.group("component")]
        value = match.group("value")
        parsed[component] = int(value) if component == "n_genes" else float(value)

    return parsed
