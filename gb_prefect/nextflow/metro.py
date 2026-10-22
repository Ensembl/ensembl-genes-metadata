"""Register a Nextflow run with an nf-metro live-progress server -- standard library only, no Prefect.

nf-metro (https://github.com/seqeralabs/nf-metro) `serve-multi` draws each registered run's
pipeline map and lights up its stations from the run's Nextflow weblog events:

1. register_run POSTs the pipeline's map (.mmd) to <server>/maps and gets the run's URLs back;
2. Nextflow is launched with -with-weblog <events URL>;
3. the live map is at the view URL.

The server keeps runs in memory only (its newest 100), so a view URL stops working after a
server restart or once 100 newer runs were registered.
"""

import json
import urllib.parse
import urllib.request
from dataclasses import dataclass
from typing import Optional

# Seconds to wait for the server; registering is skipped rather than holding up the run.
REGISTER_TIMEOUT_S = 10


@dataclass(frozen=True)
class MetroRun:
    """A run registered with the nf-metro server.

    Attributes:
        id: The server's id for the run.
        view_url: The run's live map, for people.
        events_url: Where Nextflow sends its events (-with-weblog); includes the token if one
            is set, so treat it as a secret.
    """

    id: str
    view_url: str
    events_url: str


def register_run(
    server_url: str,
    map_text: str,
    name: str,
    token: Optional[str] = None,
    timeout: float = REGISTER_TIMEOUT_S,
) -> MetroRun:
    """Register a run with the nf-metro server at server_url and return its URLs.

    The token goes in a header for registering, and in the events URL's query string, since
    Nextflow's weblog can't send headers. Raises urllib.error.URLError (server unreachable),
    HTTPError (e.g. 401 bad token, 400 invalid map) or ValueError (unexpected response).
    """
    base = server_url.rstrip("/") + "/"
    request = urllib.request.Request(
        urllib.parse.urljoin(base, "maps") + "?" + urllib.parse.urlencode({"name": name}),
        data=map_text.encode("utf-8"),
        headers={"Content-Type": "text/plain; charset=utf-8", **({"X-Metro-Token": token} if token else {})},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:  # nosec B310 -- configured URL
        body = json.loads(response.read().decode("utf-8"))
    try:
        run_id, view, events = body["id"], body["view"], body["events"]
    except (KeyError, TypeError) as err:
        raise ValueError(f"Unexpected response from {base}maps: {body!r}") from err

    events_url = urllib.parse.urljoin(base, events)
    if token:
        separator = "&" if urllib.parse.urlparse(events_url).query else "?"
        events_url += separator + urllib.parse.urlencode({"token": token})
    return MetroRun(id=str(run_id), view_url=urllib.parse.urljoin(base, view), events_url=events_url)
