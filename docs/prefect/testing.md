# Testing flows

Flows can be tested without a Prefect server, a cluster or Nextflow. Tests live in `tests/` and run
with pytest:

```bash
pip install ".[cicd]"
```

```bash
python -m pytest tests
```

Use a **fresh venv** with the package installed: stale installs give misleading import errors.

## Calling flows and tasks without Prefect

A `@flow` or `@task` keeps the plain function as `.fn`. Calling `.fn` runs your code without Prefect's
engine or a server:

```python
result = run_nextflow_pipeline.fn(run, TaskRunOptions(create_artifact=False))
output = gb_is_reference_flow.fn(file_path="/abs/gcas.txt", outdir=str(tmp_path), dry_run=True)
```

Pass `TaskRunOptions(create_artifact=False)`, or replace the artifact function, since creating an
artifact needs a server:

```python
monkeypatch.setattr(artifact_utils, "create_markdown_artifact", lambda **kw: published.append(kw))
```

## What to test

**The spec**, for each input mode: it is plain Python, so check the command it renders.

```python
from gb_prefect.nextflow.command import nextflow_argv

run = asm_qc_run("/out", "/enscode", "{}", gca_list="/in/gcas.txt")
argv = nextflow_argv(run, None, None)
assert argv[:3] == ["nextflow", "-C", "/enscode/ensembl-genes-metadata/pipelines/assembly_qc/nextflow.config"]
assert "--gca_list" in argv and "--min_contig_n50" not in argv   # None is left out
```

**The flow**, with the task replaced, to check what it builds from its parameters (defaults, folder
names, credentials, input checks):

```python
def test_flow_builds_the_run(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(flow_module, "run_nextflow_pipeline", lambda run, options: calls.append(run))
    monkeypatch.setenv("ENSCODE", "/enscode")
    flow_module.gb_asm_qc_flow.fn(gca_list=str(gcas), outdir=str(tmp_path),
                                  credentials=PipelineCredentials('{"host": "db"}'))
    [run] = calls
    assert run.outdir.startswith(str(tmp_path / "asm_qc_"))
```

**A dry run**, which needs no credentials or input files on the machine: pass `dry_run=True` and no
`credentials`; the run uses placeholders.

**End to end** with the real task and a fake `nextflow`, from the `fake_nextflow` fixture in
`tests/conftest.py`. It puts a `nextflow` on `PATH` that prints its arguments and exits; control it
with environment variables:

| Variable | Effect |
|---|---|
| `FAKE_NF_RC` | its exit code (default 0) |
| `FAKE_NF_EXTRA` | extra output to print (e.g. a Nextflow `-- Error is ignored` line) |
| `FAKE_NF_ARGS_FILE` | a file to write its exact arguments to, one per line (the log is redacted) |

```python
def test_end_to_end(tmp_path, monkeypatch, fake_nextflow):
    monkeypatch.setattr(flow_module, "run_nextflow_pipeline", run_nextflow_pipeline.fn)
    result = flow_module.gb_asm_qc_flow.fn(..., credentials=CREDENTIALS)
    assert result.returncode == 0
    assert "s3cr3t" not in Path(result.log_file).read_text()   # credentials never logged
```

For script flows, run real commands (`[sys.executable, "-c", "..."]`) through `run_script.fn`, or put
a stand-in script under a fake ENSCODE; see `tests/test_script.py`.

## Things the test setup does for you

- `gb_metadata` imports resolve to this checkout's `src/python`, not an installed copy.
- Live maps are switched off (an autouse fixture), so no test asks a server for `nf-metro-url`; to
  test them, set `load_metro_config` yourself, as `tests/test_metro.py` does with a fake server.

## Before deploying

1. The tests pass in a fresh venv.
2. A dry run of the flow shows the command you expect in the run folder's `*_command.sh`.
3. A `_test` deployment ([Deploying flows](deploying.md)) runs for real on Codon, first with
   `dry_run=true`.
