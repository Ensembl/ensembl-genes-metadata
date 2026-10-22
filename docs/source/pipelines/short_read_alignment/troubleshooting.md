## Troubleshooting

### Invalid input samplesheet

If the pipeline fails during startup with a schema validation error, check that your input CSV:

- Contains all required columns.
- Uses the correct column names.
- Does not contain empty values in required fields.
- Uses `true` or `false` for the `paired` column.

You can validate the samplesheet against the schema in `assets/schema_input.json`.

---

### Genome download fails

If the reference genome cannot be downloaded:

- Verify that the `gca` accession is valid.
- Check your internet connection.
- Ensure the NCBI Datasets API is reachable.
- If the genome is already available locally, specify it using the `genome_file` column.

---

### FASTQ download fails

Common causes include:

- Incorrect FASTQ URLs or paths.
- FTP server temporarily unavailable.
- Network connectivity problems.

If downloads repeatedly fail, verify that the files are accessible manually:

```bash
wget <FASTQ_URL>
```

---

### Corrupted BAM files

The pipeline validates BAM files using `samtools quickcheck`.

If a BAM fails validation, it is excluded from downstream analyses.

Common causes are:

- Interrupted alignment jobs.
- Filesystem issues.
- Incomplete file transfers.
- Missing BAM EOF marker.

Inspect the corresponding process log in the Nextflow work directory for details.

---

### BAM indexing fails

Large genomes may contain chromosomes exceeding the size limit supported by BAI indexes.

The pipeline automatically falls back to creating a CSI index when BAI indexing is not possible.

---

### Tissue merge fails

The merge step requires that all BAM files within a tissue:

- are coordinate sorted,
- originate from the same reference genome,
- have valid BAM headers.

If `samtools merge` reports header incompatibilities, verify that all samples were aligned against the same genome assembly.

---

### BigWig generation fails

`bamCoverage` requires:

- an indexed BAM file,
- a coordinate-sorted BAM,
- sufficient disk space.

If stranded BigWigs are requested, ensure that the stranded BAM generation step completed successfully.

---

### CRAM conversion fails

CRAM conversion requires access to the reference genome used during alignment.

If CRAM creation fails, verify that:

- the reference FASTA exists,
- the FASTA index (`.fai`) is available,
- the BAM header matches the reference.

---

### Running out of disk space

Alignment files can be large.

Consider:

- enabling FASTQ deletion after alignment,
- storing the Nextflow work directory on a filesystem with sufficient capacity,
- periodically cleaning the work directory after successful runs.

---

### Filesystem latency

On networked filesystems (e.g. Lustre, GPFS, NFS), recently created files may not be immediately visible.

Increase the value of:

```bash
--files_latency
```

if you observe intermittent "file not found" errors.

---

### Resuming failed runs

After fixing the underlying issue, rerun the workflow using:

```bash
nextflow run main.nf -resume
```

Nextflow will reuse completed processes and execute only the missing or failed steps.

---

### Inspecting logs

Useful files include:

- `.nextflow.log` — overall workflow log.
- `work/<hash>/.command.log` — process stdout/stderr.
- `work/<hash>/.command.err` — process errors.
- `work/<hash>/.command.sh` — executed script.
- `work/<hash>/.exitcode` — process exit status.

These files are usually sufficient to diagnose most failures.