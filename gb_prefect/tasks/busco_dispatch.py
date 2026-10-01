from typing import Any, Dict, List

from prefect import task  # type: ignore
from prefect.deployments import run_deployment  # type: ignore
from prefect.runtime import flow_run  # type: ignore

from gb_metadata.update_busco_events import get_busco_status


@task(log_prints=True)
def dispatch_rows(
    rows: List[Dict[str, str]],
    db_params: Dict[str, Any],
    force: bool,
    dry_run: bool,
    deployment_name: str,
) -> Dict[str, List[str]]:
    """Trigger deployment_name for each row {gca, taxon_id, busco_dataset} not already done.

    A GCA already marked genome_busco.status=done is skipped unless force=True. Shared by
    every BUSCO_dispatch input mode (CSV, bare GCA list, ...) -- only how `rows` gets built
    differs between them; this is the part that actually decides what runs.
    """
    triggered: List[str] = []
    skipped: List[str] = []

    for row in rows:
        gca = row["gca"]
        status = get_busco_status(gca, db_params)

        if status == "done" and not force:
            print(f"Skipping {gca}: already done (use --force to override)")
            skipped.append(gca)
            continue

        if dry_run:
            print(f"[dry run] Would trigger {deployment_name} for {gca}")
        else:
            # Scoped to this dispatch run, not the bare gca -- an unscoped key would dedupe
            # against ANY past run for this gca forever (including from a previous dispatch,
            # hours or days ago), silently returning the old run instead of creating a new
            # one, with no error. Scoping to flow_run.id still protects against the same row
            # being double-triggered within this one dispatch run, without blocking a later,
            # separate dispatch from ever re-triggering the same gca again.
            run_deployment(
                name=deployment_name,
                parameters={
                    "gca": gca,
                    "taxon_id": row["taxon_id"],
                    "busco_dataset": row.get("busco_dataset") or None,
                },
                timeout=0,
                idempotency_key=f"{flow_run.id}-{gca}",
            )
            print(f"Triggered {deployment_name} for {gca}")

        triggered.append(gca)

    print(f"Triggered: {len(triggered)}, skipped: {len(skipped)}")
    return {"triggered": triggered, "skipped": skipped}
