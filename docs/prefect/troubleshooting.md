# Troubleshooting

Where to look, from first to last:

1. **The flow run in the Prefect UI**: its state, logs (they include Nextflow's or the script's
   output, live) and artifacts (the run summary with the command and the end of the log).
2. **The run folder** (`<outdir>/<folder_name>/`, or `<outdir>/<gca>/` for BUSCO):
   `nextflow_<name>.log` or `script_<name>.log`, and the command record `*_command.sh`.
3. **The job folder** (`<outdir>/slurm_jobs/`): the worker's sbatch script and Slurm's
   `slurm_<job>.out` / `.err`. Look here when the run crashed before the flow started.
4. **Nextflow's own files** in the run folder: `.nextflow.log`, and a failed task's
   `work/<xx>/<hash>/.command.log`, `.command.err` and `.exitcode`.

## Common problems

| Symptom | Cause | Fix |
|---|---|---|
| Run **Crashed** at once; worker log: `PermissionError ... slurm_cli_<id>.sh` | The worker's Unix user can't write the deployment's `outdir`. | Use an `outdir` that user can write, or run your own worker as yourself for tests. |
| Run stays **Scheduled/Late** | No worker is polling the pool, or a concurrency limit is full. | Start the pool's worker; check limits in the UI. |
| Nextflow `failed with return code -9`, and the worker log says `ProcessWorker` | A Nextflow flow was deployed to a **process** pool: the head job ran on the worker's machine and was killed. | Deploy it to a `slurm-cli` pool. |
| `Could not start nextflow ... check the work pool's setup_commands` (return code 127) | `nextflow` isn't on `PATH` in the job: the module isn't loaded. | Check `--nextflow-module` and the module name; redeploy. |
| Run fails at once with a missing file for the entrypoint | Flows were moved or renamed after deploying. | Redeploy. |
| `ImportError` / old behaviour from `gb_metadata` | The venv has an old install of `src/python`. | Reinstall in the venv and `git pull` its checkout. |
| A task fails with `module(s) are unknown: "singularity-..."` | The pipeline loads a module the cluster no longer has. | Update the module in the pipeline's config. |
| Run **Failed**: `exited 0, but N task(s) failed and were ignored by the pipeline` | The pipeline's `errorStrategy 'ignore'` let tasks fail; Nextflow said "Successful", but those tasks produced nothing. | Look at the listed tasks' `work/` folders. Only set `fail_on_ignored_errors=False` if such failures are really harmless. |
| `sbatch: error: Batch job submission failed: Socket timed out` in the log | Slurm controller trouble. | Usually transient: Nextflow retries; rerun if it gave up. |
| `ENSCODE is required when dry_run=False` | No `enscode` parameter and no `$ENSCODE`. | Set the deployment's `enscode` (`--enscode` when deploying). |
| `Input file must be an absolute path` / `Input file not found` | The input file is relative or not on the cluster filesystem. | Give an absolute path visible from Codon. |
| `WARNING: Live map not available (...)` | The nf-metro server is down, the token is wrong, or the map isn't in ENSCODE. The run continues without a live map. | Check the server, the `nf-metro-token` block, and `git pull` ENSCODE. |
| A "Live pipeline map" link shows nothing | The nf-metro server restarted, or has registered 100 newer runs since (it keeps its newest 100, in memory). | Nothing to fix for old runs. |
| A GCA stays `in_progress` in BUSCO | Its Slurm job was killed outright (time limit, out of memory), so the flow's failure step never ran. | Check the job's `.err`; raise the deployment's `time_limit` or memory; reset the status. |

## Re-running by hand

The command record in the run folder is the exact command that ran, as a bash script:

```bash
bash /path/to/run_folder/nextflow_<name>_command.sh
```

For a Nextflow run, credentials are in the `-params-file` it points to (owner-only, removed by some
flows' cleanup); for a script run, set the environment variables the record lists as not recorded.
Run it from a Slurm job or an interactive session with `nextflow` and the venv loaded, as the
deployed job would.
