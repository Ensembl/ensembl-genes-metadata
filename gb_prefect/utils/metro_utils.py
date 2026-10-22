"""Where the nf-metro live-progress server is, from Prefect (see gb_prefect/nextflow/metro.py).

Live maps are on when the Prefect Variable METRO_URL_VARIABLE holds the server's URL, e.g.
http://genebuild-metadata-prefect.ebi.ac.uk:8081; delete the Variable to switch them off. The
server's token, if it uses one, is the Secret block METRO_TOKEN_SECRET_BLOCK.
"""

from dataclasses import dataclass
from typing import Optional

from prefect.blocks.system import Secret  # type: ignore
from prefect.variables import Variable  # type: ignore

METRO_URL_VARIABLE = "nf-metro-url"
METRO_TOKEN_SECRET_BLOCK = "nf-metro-token"


@dataclass(frozen=True)
class MetroConfig:
    """The nf-metro server's URL and token (None if the server runs without one)."""

    url: str
    token: Optional[str] = None


def load_metro_config() -> Optional[MetroConfig]:
    """Return the nf-metro server configured in Prefect, or None if live maps are off.

    A missing token Secret block means the server runs without a token. Errors talking to
    Prefect are raised; the caller decides whether that should stop the run.
    """
    url = Variable.get(METRO_URL_VARIABLE, default=None)
    if not url:
        return None
    try:
        token = Secret.load(METRO_TOKEN_SECRET_BLOCK).get()
    except ValueError:  # no such block
        token = None
    return MetroConfig(url=str(url), token=str(token) if token else None)
