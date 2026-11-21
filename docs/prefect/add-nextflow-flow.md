# Add a Nextflow flow

This walks through adding a flow for a new Nextflow pipeline, using a made-up pipeline,
`assembly_qc`, that checks the assemblies in a GCA list. Each step shows the code to write; the
real flow it is modelled on is
[`gb_registry_by_gca.py`](https://github.com/Ensembl/ensembl-genes-metadata/blob/main/gb_prefect/flows/assembly_metadata/gb_registry_by_gca.py)
with its spec
[`assembly_metadata.py`](https://github.com/Ensembl/ensembl-genes-metadata/blob/main/gb_prefect/pipelines/assembly_metadata.py).
Read the [overview](index.md) first if Prefect is new to you.

A Nextflow flow is three small pieces:

1. a **spec**: plain Python that knows the pipeline's command line and returns a `NextflowRun`;
2. a **flow**: resolves the run's inputs (ENSCODE, credentials, folders) and calls the spec, then
   the task;
3. a **deployment**: one entry in `deploy_flows.py`.

You never build a command line, an sbatch script, a log or an artifact yourself: the generic task
`run_nextflow_pipeline` does all of that.

## 1. Write the spec

Create `gb_prefect/pipelines/assembly_qc.py`. One module per pipeline, or per pipeline *mode* when
one pipeline has several (e.g. `statistics_busco_genome.py`). Specs import nothing from Prefect, so
they are easy to test.

```python
"""Spec for the assembly_qc pipeline."""

from typing import Optional

from gb_prefect.models.nextflow import NextflowRun


def asm_qc_run(
    outdir: str,
    enscode: str,
    metadata_params_string: str,
    gca_list: str,
    min_contig_n50: Optional[int] = None,
) -> NextflowRun:
    """Return the NextflowRun of the assembly_qc pipeline for the accessions in gca_list."""
    return NextflowRun(
        name="assembly_qc",                                         # names log, record, artifact
        enscode=enscode,
        pipeline="ensembl-genes-metadata/pipelines/assembly_qc/main.nf",   # relative to ENSCODE
        config="ensembl-genes-metadata/pipelines/assembly_qc/nextflow.config",
        outdir=outdir,                                              # launch dir and outputs
        params={                                                    # --key value on the command line
            "output_dir": outdir,
            "gca_list": gca_list,
            "min_contig_n50": min_contig_n50,                       # None: left out
            "strict": True,                                         # True: bare flag --strict
        },
        secret_params={"metadata_params_string": metadata_params_string},  # -params-file, 0600
        with_report=True,
        dag_file="assembly_qc_dag_%Y%m%d-%H%M%S.mmd",
        metro_map="ensembl-genes-metadata/pipelines/assembly_qc/assets/metro_map.mmd",
    )
```

What each `NextflowRun` field does (full list in
[`models/nextflow.py`](https://github.com/Ensembl/ensembl-genes-metadata/blob/main/gb_prefect/models/nextflow.py)):

| Field | Notes |
|---|---|
| `name` | Short and unique within `outdir`: files are `nextflow_<name>.log`, `nextflow_<name>_command.sh`, `nextflow_<name>_params.json`. |
| `pipeline`, `config` | Paths **relative to ENSCODE**, so the version that runs is the ENSCODE checkout's. `config` is passed with `-C`, which makes Nextflow use *only* that file: pass the pipeline's own `nextflow.config`, or leave `config` out. |
| `params` | Rendered as `--key value`. `True` becomes a bare flag, `False` becomes `--key false`, `None` is left out. Strings, numbers and paths only. |
| `secret_params` | Anything secret (DB and Slack credentials, structured values): written to an owner-only `-params-file`, never on the command line, in logs or in artifacts. A key can't be in both `params` and `secret_params`. |
| `csv_input` | When the *flow* builds the input rows itself (e.g. one GCA per run): a `CsvInput(param="csvFile", columns=[...], rows=[...])` is written to `<outdir>/<name>.csv` and passed as `--csvFile`. A file the *user* gives (like `gca_list` here) is just a path in `params`. |
| `dag_file` | Written with `-with-dag`; `%Y%m%d-%H%M%S` is filled in at launch so reruns in the same folder don't clash. Use `.mmd`. |
| `metro_map` | The pipeline's nf-metro map, for live maps (step 4). Optional. |
| `fail_on_ignored_errors` | Default `True`: the run fails if Nextflow exits 0 but ignored failed tasks (`errorStrategy 'ignore'`). Keep it unless the pipeline's ignored failures are really harmless. |

Validate input combinations here, in plain Python, and raise `ValueError`: e.g.
`asm_update_run` refuses both or neither of `gca_list` and `screen_date`.

## 2. Write the flow

Create `gb_prefect/flows/assembly_metadata/gb_asm_qc.py` (or a new group folder: see the note at the
end of this step).

```python
from datetime import datetime
from typing import Optional

from prefect import flow  # type: ignore

from gb_prefect.models.pipeline_options import TaskRunOptions
from gb_prefect.pipelines.assembly_qc import asm_qc_run
from gb_prefect.tasks.nextflow import run_nextflow_pipeline
from gb_prefect.utils.credentials_utils import PipelineCredentials, resolve_credentials
from gb_prefect.utils.enscode_utils import resolve_enscode
from gb_prefect.utils.path_utils import check_input_file, run_folder


@flow(name="gb_asm_qc", log_prints=True)
def gb_asm_qc_flow(
    gca_list: str,
    outdir: str,
    min_contig_n50: Optional[int] = None,
    folder_name: Optional[str] = None,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: Optional[str] = None,
):
    """Run the assembly_qc Nextflow pipeline for the GCA accessions in gca_list.

    gca_list: absolute path to a file with one GCA accession per line.
    outdir, folder_name: the run goes to <outdir>/<folder_name>; folder_name defaults to
        asm_qc_<today>.
    min_contig_n50: passed to the pipeline; omitted, the pipeline's default.
    enscode: ENSCODE directory; defaults to $ENSCODE.
    dry_run: write and log the command without running Nextflow or loading credentials.
    credentials, metadata_secret_block: DB credentials, explicit or from the Secret block.
    """
    credentials = resolve_credentials(credentials, dry_run=dry_run, metadata_secret_block=metadata_secret_block)
    run = asm_qc_run(
        outdir=run_folder(outdir, folder_name, f"asm_qc_{datetime.now().strftime('%Y-%m-%d')}"),
        enscode=resolve_enscode(enscode, dry_run),
        metadata_params_string=credentials.metadata_params_string,
        gca_list=check_input_file(gca_list, dry_run),
        min_contig_n50=min_contig_n50,
    )
    return run_nextflow_pipeline(run, TaskRunOptions(dry_run=dry_run))
```

The conventions every flow follows, so deployments and users can rely on them:

- **Parameters are the flow's arguments.** Give every optional one a default; the docstring is the
  parameter documentation users see. Keep `outdir` and `enscode` as names: the deploy script sets
  them for every Nextflow flow.
- **`dry_run`** must work without credentials, ENSCODE or input files on the machine:
  `resolve_credentials(..., dry_run=dry_run)`, `resolve_enscode(enscode, dry_run)` and
  `check_input_file(path, dry_run)` all handle it.
- **`run_folder(outdir, folder_name, default)`** gives each run its own folder; the default should
  be dated.
- **`check_input_file`** rejects relative paths (Nextflow runs from the run folder, where a relative
  path means something else) and, unless it's a dry run, missing files.
- **Return the task's result.** `run_nextflow_pipeline` raises `NextflowRunError` on failure, which
  fails the flow run; don't catch it unless you need a post-run step (see `busco_gca_flow`, which
  marks the GCA failed and cleans up before re-raising).
- `@flow(name=...)` is the flow's name on the server; deployments are addressed as
  `<flow name>/<deployment name>`.

!!! note "A new folder under `flows/`"
    A new group folder needs an `__init__.py` and an entry in the `packages` list in
    `pyproject.toml` (e.g. `"gb_prefect.flows.my_group"`): packages are listed explicitly, not
    discovered.

## 3. Add the deployment

In
[`gb_prefect/setup/deploy_flows.py`](https://github.com/Ensembl/ensembl-genes-metadata/blob/main/gb_prefect/setup/deploy_flows.py),
add an entry to `FLOWS`:

```python
    "asm_qc": FlowDeployment(
        name="gb_asm_qc",                                   # deployment name
        entrypoint="gb_prefect/flows/assembly_metadata/gb_asm_qc.py:gb_asm_qc_flow",
        work_pool="metadata",                               # a key of POOLS
        description="Assembly QC",
        tags=["metadata_qc"],
        cpu=1,
        memory=4,                                           # GB, for the Nextflow head job
        time_limit=6,                                       # hours: must cover the whole pipeline
        concurrency_limit=1,                                # if two runs must never overlap
    ),
```

The job's `cpu`, `memory` and `time_limit` are the **Nextflow head job's**: the pipeline's own
processes get their resources from the pipeline's config. Set `time_limit` to the longest the whole
pipeline can take. See [Deploying flows](deploying.md) for pools, output folders and the other
options.

## 4. Optional: a live map

If the pipeline has an [nf-metro](https://github.com/seqeralabs/nf-metro) map, set `metro_map` in
the spec (step 1) and every real run gets a live map of its progress, linked from the flow run in
the Prefect UI. Maps live next to the pipeline:
`pipelines/<pipeline>/assets/metro_map.mmd` in this repo, or one map per mode in
`ensembl-genes-nf` (`assets/metro_map_<mode>.mmd`).

A map lists stations, the connections between them, and which Nextflow process each station stands
for; start from an existing one, e.g.
[`assembly_metadata/assets/metro_map.mmd`](https://github.com/Ensembl/ensembl-genes-metadata/blob/main/pipelines/assembly_metadata/assets/metro_map.mmd):

```text
%%metro title: Assembly QC
%%metro line: qc | QC | #2db572
%%metro process: fetch | (^|:)FETCH_ASSEMBLY$
%%metro process: qc_check | (^|:)QC_CHECK$

graph LR
    subgraph checks [Checks]
        fetch([Fetch assembly])
        qc_check([QC check])
        fetch -->|qc| qc_check
    end
```

Anchor the `%%metro process` patterns (`(^|:)NAME$`): `WRITE2DB` must not also match
`INTEGRITY_WRITE2DB`. Check the map with `pip install "nf-metro==2.1.0"`, then:

```bash
nf-metro validate pipelines/assembly_qc/assets/metro_map.mmd
```

```bash
nf-metro check-mapping pipelines/assembly_qc/assets/metro_map.mmd --processes processes.txt
```

where `processes.txt` lists the pipeline's fully qualified process names (e.g.
`ASSEMBLY_QC:FETCH_ASSEMBLY`), one per line. For this repo's pipelines, also render the README
diagram with `src/bash/render_metro_maps.sh`.

## 5. Test it

See [Testing](testing.md). At least:

- a test of the spec: the arguments it renders for each input mode;
- a dry run of the flow, which needs nothing but Python:

```bash
python -c "from gb_prefect.flows.assembly_metadata.gb_asm_qc import gb_asm_qc_flow as f; f(gca_list='/abs/gcas.txt', outdir='/tmp/qc', dry_run=True)"
```

The dry run writes `nextflow_assembly_qc_command.sh` in the run folder: check it is the command you
expect.

## 6. Deploy and run it

Commit and push your branch (the worker clones it), then deploy a test copy and run it. Use
`--name-suffix _test` so the production deployment isn't touched, and `--metadata-pool` to send it to
your own test pool if you have one:

```bash
python gb_prefect/setup/deploy_flows.py --branch <your-branch> --outdir <test_outdir> --enscode <enscode> --asm-venv <venv> --flows asm_qc --name-suffix _test
```

```bash
prefect deployment run 'gb_asm_qc/gb_asm_qc_test' --param gca_list=/abs/path/gcas.txt --param dry_run=true
```

Then without `dry_run`. In the Prefect UI, the flow run shows the live Nextflow output in its logs,
a run-summary artifact and, with a map, a "Live pipeline map" link. If something goes wrong, see
[Troubleshooting](troubleshooting.md). When it works, deploy without `--name-suffix`.

## Checklist

- [ ] Spec in `gb_prefect/pipelines/`, paths relative to ENSCODE, secrets in `secret_params`.
- [ ] Flow in `gb_prefect/flows/<group>/`: `outdir`, `folder_name`, `enscode`, `dry_run`,
      credentials; docstring documents every parameter.
- [ ] New group folder: `__init__.py` and `pyproject.toml` `packages` entry.
- [ ] `FlowDeployment` in `deploy_flows.py` with a realistic `time_limit`.
- [ ] Optional map in the pipeline's `assets/`, validated and checked against the process names.
- [ ] Tests for the spec and a dry run of the flow.
- [ ] A row for the flow in the `gb_prefect` README's Flows table.
- [ ] ENSCODE checkout updated (`git pull`) if the pipeline or map changed.
