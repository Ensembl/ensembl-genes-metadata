# Operations

How the Prefect setup on Codon is put together, and what to do after common changes. Commands run
from the repository root, with the venv that has this package installed (`pip install .`) active and
`PREFECT_API_URL` pointing at the server
(`http://genebuild-metadata-prefect.ebi.ac.uk:4200/api`). Every step can be re-run.

## Setting up from scratch

1. **Work pools.** The `slurm-cli` type comes from this package (registered as a Prefect plugin when
   it's installed):

    ```bash
    prefect work-pool create genebuild-metadata --type slurm-cli
    ```
    ```bash
    prefect work-pool create genebuild-busco --type slurm-cli
    ```
    ```bash
    prefect work-pool update genebuild-busco --concurrency-limit 10
    ```
    ```bash
    prefect work-pool create genebuild-dispatch --type process
    ```

    The BUSCO pool's concurrency limit is how many BUSCO runs go at once; the automatic dispatcher
    sizes each batch from it and refuses to run without one.

2. **Workers**, one per pool, kept running (e.g. in `tmux` or as a service), **as the user that
   owns the output directories**: everything a run does runs as that user.

    ```bash
    prefect worker start --pool genebuild-metadata --type slurm-cli
    ```
    ```bash
    prefect worker start --pool genebuild-dispatch --type process
    ```

    (and the same for `genebuild-busco`).

3. **Secret blocks**: see [Credentials](credentials.md).

4. **Deployments**: see [Deploying flows](deploying.md).

5. **Live maps (optional).** An [nf-metro](https://github.com/seqeralabs/nf-metro) server where both
   the Codon compute nodes and your browser can reach it; we use the Prefect VM, port 8081. It needs
   Python 3.11+ in its own venv (`pip install "nf-metro==2.1.0"`) and a random token
   (`python3 -c "import secrets; print(secrets.token_urlsafe(32))"`):

    ```bash
    nf-metro serve-multi --host 0.0.0.0 --port 8081 --token <token>
    ```

    Then point Prefect at it, and store the same token:

    ```bash
    prefect variable set nf-metro-url http://genebuild-metadata-prefect.ebi.ac.uk:8081
    ```
    ```bash
    python gb_prefect/setup/create_secrets.py --metro-token
    ```

    Check from a compute node that it is reachable:
    `srun -t 5 curl -s -o /dev/null -w "%{http_code}\n" http://genebuild-metadata-prefect.ebi.ac.uk:8081/`
    should print `200`.

    For a long-lived setup, run it as a service with the token in a root-readable environment file
    (`EnvironmentFile=` in a systemd unit, `--token ${NF_METRO_TOKEN}`), since nf-metro only takes the
    token on its command line.

## What uses what

| Thing | Comes from | Updated by |
|---|---|---|
| Flow code (`gb_prefect`) | the deployment's Git branch, cloned for each run | pushing the branch |
| Pipelines, scripts, nf-metro maps | the **ENSCODE** checkout (`--enscode`) | `git pull` in ENSCODE |
| `gb_metadata` (`src/python/`), imported by flows | the **installed package in the venv** (`--asm-venv`) | reinstalling it in the venv |
| Nextflow | the module loaded by the job (`--nextflow-module`) | redeploying with another module |
| Pools, schedules, Secret blocks, Variables | the Prefect server | the CLI or UI |

## Runbook

| After | Do |
|---|---|
| a change to flow code only | push the branch; nothing else |
| a change to a pipeline, script or map | merge it, then `git pull` the ENSCODE checkout |
| a change to `src/python/` | reinstall it in the venv (`pip install -e <checkout>`; `git pull` there) |
| a new, moved or renamed flow; new resources, pools or default parameters | redeploy (`deploy_flows.py`) |
| redeploying `busco-dispatch-automatic` | check its schedule is still set in the UI |
| a new pool | `prefect work-pool create ...` and start its worker |
| a credential change | `create_secrets.py ... --overwrite` |
| rotating the nf-metro token | new token on the server *and* `create_secrets.py --metro-token --overwrite`; restart the server (live maps of runs in progress stop updating) |
| switching live maps off | `prefect variable unset nf-metro-url` |
| a Slurm module removed from the cluster | update the pipeline config or `--nextflow-module`, then pull ENSCODE / redeploy |
