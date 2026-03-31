# Backed FastAPI services for the GB metadata app

Make sure you add password to conf/db_config.json

Launch backend with:
```bash
uvicorn metadata_app.backend.app.main:app --reload
```

## Genebuild Metadata API – Endpoint Overview

The backend exposes a FastAPI-based REST API for querying and managing **genebuild metadata** (assemblies, annotations, taxonomy, transcriptomics, projects, reports, and handovers), as well as some utility/maintenance endpoints.

All API routes are prefixed with `/api`.

---

### Health & Root

#### `GET /api/`

Simple root endpoint for the API.

- **Response 200**
  - `message`: Welcome string.
- **Example**
  ```json
  {
    "message": "Welcome to the Genebuild API"
  }
  ```

#### `GET /api/health`

Basic health check.

- **Response 200**
  - `status`: `"healthy"`.
- **Example**
  ```json
  {
    "status": "healthy"
  }
  ```

---

### Assemblies

Router: `metadata_app.backend.app.api.routes.assemblies` mounted at `/api/assemblies`.

This group serves endpoints for querying assembly-level metadata (e.g. genome builds).

Typical patterns (exact paths may differ; see `assemblies.py` for full details):

- `GET /api/assemblies`
  - List assemblies with basic metadata and optional filters (species, status, etc.).
- `GET /api/assemblies/{assembly_id}`
  - Detailed metadata for a single assembly.
- `GET /api/assemblies/{assembly_id}/stats`
  - Statistics or summary metrics (e.g. scaffold counts, N50) if implemented.

---

### Annotations

Router: `metadata_app.backend.app.api.routes.annotations` mounted at `/api/annotations`.

Endpoints related to **annotation runs** and gene builds.

Example patterns:

- `GET /api/annotations`
  - List annotation runs / builds with filters (assembly, bioproject, status, pipeline version).
- `GET /api/annotations/{annotation_id}`
  - Detailed metadata for a single annotation run.
- `GET /api/annotations/{annotation_id}/logs`
  - Optional logs / provenance information if provided.

---

### Taxonomy

Router: `metadata_app.backend.app.api.routes.taxonomy` mounted at `/api/taxonomy`.

Endpoints for taxonomy metadata associated with genomes and projects.

Example patterns:

- `GET /api/taxonomy/{taxon_id}`
  - Metadata for a specific taxon (scientific name, rank, lineage).
- `GET /api/taxonomy/{taxon_id}/assemblies`
  - Assemblies associated with a given taxon.

---

### Transcriptomics

Router: `metadata_app.backend.app.api.routes.transcriptomics` mounted at `/api/transcriptomics`.

Transcriptomics-related metadata, e.g. RNA-seq evidence used during annotation.

Example patterns:

- `GET /api/transcriptomics`
  - List transcriptomics datasets (bioprojects, experiments, libraries).
- `GET /api/transcriptomics/{dataset_id}`
  - Detailed metadata for a specific dataset.
- `GET /api/transcriptomics/{dataset_id}/usage`
  - How this dataset was used in genebuild (assemblies / annotations that reference it).

---

### Home Page Summary

Router: `metadata_app.backend.app.api.routes.home_page` mounted at `/api/home_page`.

Aggregated metadata for the **frontend home/dashboard**.

Common patterns:

- `GET /api/home_page`
  - Returns summary stats and “cards” for the UI, e.g.:
    - counts of assemblies, annotations, transcriptomics datasets,
    - lists of “recently updated” projects / assemblies,
    - high-level status indicators.

---

### Reports – Annotations

Router: `metadata_app.backend.app.api.routes.report_annotations` mounted at `/api/report/anno`.

Endpoints focused on reporting / summarising annotation data.

Example patterns:

- `GET /api/report/anno`
  - Returns a table-like report of annotations (one row per build).
- `GET /api/report/anno/{annotation_id}`
  - Detailed report for a specific annotation, including QC metrics and flags.

---

### Reports – Assemblies

Router: `metadata_app.backend.app.api.routes.report_assemblies` mounted at `/api/report/asm`.

Report-style endpoints for **assemblies**.

Example patterns:

- `GET /api/report/asm`
  - Summary table of assemblies, including statuses and key metrics.
- `GET /api/report/asm/{assembly_id}`
  - Detailed report for a specific assembly, potentially aggregating multiple sub‑records.

---

### Bioproject Search

Router: `metadata_app.backend.app.api.routes.bioproject_search` mounted at `/api/bioproject_search`.

Search endpoints keyed around **BioProject IDs**.

Likely patterns:

- `GET /api/bioproject_search`
  - Query assemblies/annotations by `bioproject`, free‑text search, or partial IDs.
- `GET /api/bioproject_search/{bioproject_id}`
  - All metadata linked to a given BioProject.

---

### Taxonomy Search

Router: `metadata_app.backend.app.api.routes.taxonomy_search` mounted at `/api/taxonomy_search`.

Search endpoints driven by taxonomy and free text.

Example patterns:

- `GET /api/taxonomy_search`
  - Search by scientific name, common name, or taxon ID.
- `GET /api/taxonomy_search/{query}`
  - Convenience endpoint to match a single query string to taxa and related assemblies.

---

### Projects

Router: `metadata_app.backend.app.api.routes.project` mounted at `/api/project`.

Higher-level **project** concept; a project may link assemblies, annotations, datasets, and handovers.

Likely patterns:

- `GET /api/project`
  - List projects, with filtering (status, owner, species, etc.).
- `GET /api/project/{project_id}`
  - Detailed project view aggregating all related metadata.
- `POST /api/project`
  - (If implemented) create or register a new project.
- `PATCH /api/project/{project_id}`
  - (If implemented) update project metadata.

---

### Handovers

Router: `metadata_app.backend.app.api.routes.handover` mounted at `/api/handover`.

Endpoints describing **handover** of data between stages (e.g. from genebuild to downstream consumers).

Example patterns:

- `GET /api/handover`
  - List handover records (who handed over what, when, and to where).
- `GET /api/handover/{handover_id}`
  - Detailed metadata for a single handover.
- `POST /api/handover`
  - (If implemented) register a new handover event.

---

### Database Cleaning / Maintenance

Router: `metadata_app.backend.app.api.routes.db_clean` mounted at `/api/clean`.

Utility endpoints to perform database clean‑up tasks (used with care, likely internal‑only).

Example patterns:

- `POST /api/clean/run`
  - Run a cleaning job (e.g. remove orphan / stale records).
- `GET /api/clean/status`
  - Inspect status/log of the most recent clean‑up operation.

These endpoints should be **restricted in production** (e.g. behind auth, internal network only).

---

### Error Handling

All endpoints share common error handling:

- **HTTPException**:
  - Returned as:
    ```json
    {
      "message": "<error detail>"
    }
    ```
  - With the appropriate HTTP status (e.g. 404, 400).

- **Unhandled exceptions**:
  - Logged at `ERROR` level.
  - Returned as:
    ```json
    {
      "message": "Internal server error"
    }
    ```
  - With status code `500`.

---

### Frontend Integration

The Next.js frontend is served separately from the API under `/` (root) and static asset routes. Non‑`/api/...` requests are routed to the frontend, while all JSON endpoints described here live under `/api`.

For concrete request/response models and query parameters, refer to the individual router modules in `metadata_app/backend/app/api/routes/`.