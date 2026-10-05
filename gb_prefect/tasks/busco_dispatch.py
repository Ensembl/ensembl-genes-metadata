from typing import Any, Dict, List

from prefect import task  # type: ignore
from prefect.client.orchestration import get_client  # type: ignore
from prefect.client.schemas.filters import (  # type: ignore
    DeploymentFilter,
    DeploymentFilterId,
    FlowRunFilter,
    FlowRunFilterState,
    FlowRunFilterStateType,
)
from prefect.client.schemas.objects import StateType  # type: ignore
from prefect.deployments import run_deployment  # type: ignore
from prefect.runtime import flow_run  # type: ignore

from gb_metadata.update_busco_events import claim_busco, get_busco_status, release_busco_claim


@task(log_prints=True)
def get_free_slots(deployment_name: str) -> int:
    """Return how many more runs of deployment_name ("flow_name/deployment_name") can be
    dispatched now: its work pool's concurrency_limit, minus active slots (Pending/Running),
    minus deployment_name's runs still waiting to start (Scheduled, including Late).

    The pool is looked up from the deployment at run time, so it can't drift from where
    BUSCO_gca actually runs. Waiting runs are subtracted because the pool doesn't count them:
    without that, a full pool would get more runs queued on every dispatch.

    Raises if the pool has no concurrency_limit set -- automatic dispatch refuses to run
    without one, rather than guessing a fallback cap, to avoid flooding the cluster.
    """
    with get_client(sync_client=True) as client:
        deployment = client.read_deployment_by_name(deployment_name)
        pool_name = deployment.work_pool_name
        if not pool_name:
            raise ValueError(f"Deployment '{deployment_name}' has no work pool.")
        status = client.read_work_pool_concurrency_status(pool_name)
        waiting = client.count_flow_runs(
            deployment_filter=DeploymentFilter(id=DeploymentFilterId(any_=[deployment.id])),
            flow_run_filter=FlowRunFilter(
                state=FlowRunFilterState(type=FlowRunFilterStateType(any_=[StateType.SCHEDULED]))
            ),
        )

    if status.concurrency_limit is None:
        raise ValueError(
            f"Work pool '{pool_name}' has no concurrency_limit set -- automatic dispatch "
            "refuses to run without one."
        )

    free = max(status.concurrency_limit - status.active_slots - waiting, 0)
    print(
        f"Pool '{pool_name}' (for {deployment_name}): {status.active_slots}/"
        f"{status.concurrency_limit} slots active, {waiting} run(s) waiting to start, {free} free"
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

    A GCA already marked genome_busco.status=done or in_progress is skipped unless
    force=True. Shared by every BUSCO_dispatch input mode (CSV, bare GCA list, automatic)
    -- only how `rows` gets built differs between them; this is the part that actually
    decides what runs.

    update_registry=False (the default here) just triggers the run, with no DB write. When
    True, the GCA is first *claimed* -- genome_busco.status set to in_progress only if it
    still has the status read a moment before (claim_busco) -- and only then triggered; if
    the trigger fails, the claim is released back to that status. A GCA whose status changed
    in between (e.g. claimed by another dispatcher) is skipped, so it can't be dispatched
    twice, and a late write can't overwrite a done. update_registry=True is also passed on
    to the triggered run, so the Nextflow pipeline's --update_registry loads the results and
    marks it done on completion. A GCA with no status row yet is triggered without a claim.

    post_run_cleanup is passed on to each triggered run (see busco_gca_flow), to mark
    failures and clean run directories afterwards: always on in the automatic dispatcher,
    optional (default off) in the CSV / GCA-list ones.
    """
    triggered: List[str] = []
    skipped: List[str] = []

    for row in rows:
        gca = row["gca"]
        status = get_busco_status(gca, db_params)

        if status in ("done", "in_progress") and not force:
            print(f"Skipping {gca}: already {status} (use --force to override)")
            skipped.append(gca)
            continue

        # in_progress can't be claimed again (the UPDATE would change nothing); with --force
        # it's re-triggered as is.
        claim_from = status if update_registry and status not in (None, "in_progress") else None

        if dry_run:
            if claim_from:
                print(f"[dry run] Would claim {gca} ({claim_from} -> in_progress)")
            print(f"[dry run] Would trigger {deployment_name} for {gca}")
            triggered.append(gca)
            continue

        if claim_from:
            if not claim_busco(gca, (claim_from,), db_params):
                print(f"Skipping {gca}: status changed from {claim_from} before it could be claimed")
                skipped.append(gca)
                continue
            print(f"Claimed {gca} ({claim_from} -> in_progress)")

        try:
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
        except Exception:
            if claim_from:
                # No run was created, so don't leave a ghost in_progress behind.
                release_busco_claim(gca, claim_from, db_params)
                print(f"Trigger failed for {gca}; released claim back to {claim_from}")
            raise
        print(f"Triggered {deployment_name} for {gca}")
        triggered.append(gca)

    print(f"Triggered: {len(triggered)}, skipped: {len(skipped)}")
    return {"triggered": triggered, "skipped": skipped}
