# Output

The pipeline's primary output is the populated Genebuild assembly metadata
database: new `assembly`, `species`, `organism`, `taxonomy` and ToL ID
records for every accession registered during the run. Alongside that, it
writes report and intermediate files to `<output_dir>`.

## Final outputs

Located in `<output_dir>/nextflow_output/`:

| File | Description |
|------|-------------|
| `report.txt` | Summary report of newly registered assemblies: counts by assembly type/level, assemblies from relevant BioProjects, RefSeq availability, and flagged issues (invalid taxon IDs, missing BioSample or submitter info) |
| `gca_to_run_ncbi.csv` | Input for the BUSCO Nextflow pipeline, listing shortlisted assemblies for genome quality assessment |
| `gca_list_to_report.txt` | List of all GCA accessions processed in the run |

## Intermediate outputs per accession

Located in `<output_dir>/nextflow_output/<GCA>/`:

| File | Description |
|------|-------------|
| `<GCA>_assembly.json` | Raw assembly metadata retrieved from NCBI |
| `<GCA>_metadata.json` | Processed metadata with updated database keys |
| `<GCA>_species.json` | Species information validated against NCBI Taxonomy |
| `<GCA>_tolid.json` | ToL ID information from the Darwin Tree of Life API |
| `<GCA>_assembly.last_id` | Database ID assigned to the assembly record |
| `<GCA>_metadata.last_id` | Database ID assigned to the metadata record |
| `<GCA>_species.last_id` | Database ID assigned to the species record |
| `<GCA>_tolid.last_id` | Database ID assigned to the ToL ID record |

`gca_to_run_ncbi.csv` is meant to be fed directly into the downstream BUSCO
statistics pipeline — see [Parameters](parameters.md) for how the shortlist
is chosen.
