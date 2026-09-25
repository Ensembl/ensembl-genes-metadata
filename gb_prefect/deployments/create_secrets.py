"""Create (or update) the Prefect Secret blocks holding the assembly registry credentials.

The registry update flows load their credentials from these blocks at run time, so the
database password and Slack token never appear in a deployment's parameters.

Each JSON file is read once, validated, and stored in Prefect; the file itself is not needed
afterwards (delete it, or keep it somewhere owner-readable only).

Example (run from the repository root, with PREFECT_API_URL pointing at the Prefect server):

    PYTHONPATH=. python gb_prefect/deployments/create_secrets.py \
        --metadata-json ~/secrets/metadata_params.json \
        --slack-json ~/secrets/slack_params.json

metadata_params.json holds the metadata DB connection, same keys the pipeline expects:
    {"host": "...", "port": 1234, "user": "...", "password": "...", "database": "..."}
slack_params.json holds the Slack bot connection:
    {"slack_bot_token": "xoxb-..."}
"""

import argparse
import json
from pathlib import Path

from prefect.blocks.system import Secret  # type: ignore

from gb_prefect.utils.credentials_utils import DEFAULT_METADATA_SECRET_BLOCK, DEFAULT_SLACK_SECRET_BLOCK


def save_secret_from_json(json_path: str, block_name: str, overwrite: bool) -> None:
    """Validate a JSON file and store its content as a Prefect Secret block."""
    value = json.loads(Path(json_path).expanduser().read_text(encoding="utf-8"))
    Secret(value=value).save(block_name, overwrite=overwrite)
    print(f"Saved Secret block '{block_name}' with keys: {sorted(value)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata-json", help="JSON file with metadata DB connection parameters.")
    parser.add_argument("--slack-json", help="JSON file with Slack bot connection parameters.")
    parser.add_argument("--metadata-block", default=DEFAULT_METADATA_SECRET_BLOCK)
    parser.add_argument("--slack-block", default=DEFAULT_SLACK_SECRET_BLOCK)
    parser.add_argument(
        "--overwrite", action="store_true", help="Replace the blocks if they already exist."
    )
    args = parser.parse_args()

    if not args.metadata_json and not args.slack_json:
        parser.error("Pass at least one of --metadata-json / --slack-json.")
    if args.metadata_json:
        save_secret_from_json(args.metadata_json, args.metadata_block, args.overwrite)
    if args.slack_json:
        save_secret_from_json(args.slack_json, args.slack_block, args.overwrite)
