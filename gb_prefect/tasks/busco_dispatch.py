from typing import Any, Dict, List

from prefect import task  # type: ignore
from prefect.client.orchestration import get_client  # type: ignore
from prefect.deployments import run_deployment  # type: ignore
from prefect.runtime import flow_run  # type: ignore

from gb_metadata.update_busco_events import get_busco_status, mark_busco_in_progress


@task(log_prints=True)
def get_deployment_pool(deployment_name: str) -> str:
    """Return the work pool deployment_name ("flow_name/deployment_name") runs on.

    Looked up at run time rather than configured separately, so the pool the automatic
    dispatcher sizes its batch against can't drift from where BUSCO_gca actually runs.
    """
    with get_client(sync_client=True) as client:
        deployment = client.read_deployment_by_name(deployment_name)

    if not deployment.work_pool_name:
        raise ValueError(f"Deployment '{deployment_name}' has no work pool.")
    print(f"Deployment '{deployment_name}' runs on work pool '{deployment.work_pool_name}'")
    return deployment.work_pool_name


@task(log_prints=True)
def get_free_slots(pool_name: str) -> int:
    """Return free concurrency slots for pool_name (concurrency_limit - active_slots).

    Raises if the pool has no concurrency_limit set -- automatic dispatch refuses to run
    without one, rather than guessing a fallback cap, to avoid flooding the cluster.
    """
    with get_client(sync_client=True) as client:
        status = client.read_work_pool_concurrency_status(pool_name)

    if status.concurrency_limit is None:
        raise ValueError(
            f"Work pool '{pool_name}' has no concurrency_limit set -- automatic dispatch "
            "refuses to run without one."
        )

    free = max(status.concurrency_limit - status.active_slots, 0)
    print(
        f"Pool '{pool_name}': {status.active_slots}/{status.concurrency_limit} slots active, "
        f"{free} free"
    )
    return free


@task(log_prints=True)
def dispatch_rows(
    rows: List[Dict[str, str]],
    db_params: Dict[str, Any],
    force: bool,
    dry_run: bool,
    deployment_name: str,
    update_registry: bool = False,
    post_run_cleanup: bool = False,
) -> Dict[str, List[str]]:
    """Trigger deployment_name for each row {gca, taxon_id, busco_dataset} not already done.

    A GCA already marked genome_busco.status=done is skipped unless force=True. Shared by
    every BUSCO_dispatch input mode (CSV, bare GCA list, ...) -- only how `rows` gets built
    differs between them; this is the part that actually decides what runs.

    update_registry=False (the default) just triggers the run, with no DB write here. When
    True, genome_busco.status is updated to in_progress right after a successful trigger,
    and update_registry=True is also passed on to the triggered run, so the Nextflow
    pipeline's --update_registry loads the results and marks it done on completion.

    post_run_cleanup is passed on to each triggered run (see busco_gca_flow), to mark
    failures and clean run directories afterwards: always on in the automatic dispatcher,
    optional (default off) in the CSV / GCA-list ones.
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
            if update_registry:
                print(f"[dry run] Would mark {gca} as in_progress")
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
                    "update_registry": update_registry,
                    "post_run_cleanup": post_run_cleanup,
                },
                timeout=0,
                idempotency_key=f"{flow_run.id}-{gca}",
            )
            print(f"Triggered {deployment_name} for {gca}")

            if update_registry:
                # Only after a successful trigger -- if run_deployment() itself had failed,
                # marking in_progress first would leave a stuck ghost status with no job
                # actually submitted.
                mark_busco_in_progress(gca, db_params, execute=True)
                print(f"Marked {gca} as in_progress")

        triggered.append(gca)

    print(f"Triggered: {len(triggered)}, skipped: {len(skipped)}")
    return {"triggered": triggered, "skipped": skipped}
