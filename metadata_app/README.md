# Metadata App

The **Metadata App** is a web application for exploring and managing metadata associated with Ensembl Genebuild pipelines. It provides a simple API and frontend to surface information about genomes, genebuilds, and associated analyses.

---

## Features

- Browse and search metadata for genomes, genebuilds, and assemblies.
- Simple REST API for programmatic access to metadata.
- Docker-based deployment for reproducible environments and easy local development.

---

## Repository Structure

The app is split into a backend API and a frontend UI:

- `backend/` – backend service (API, data access, validation logic).
- `frontend/` – frontend web application (UI for browsing and querying metadata).

Refer to the `README` files inside `backend/` and `frontend/` for component‑specific details.

---

## Getting Started

### Prerequisites

- Git
- Docker and Docker Compose
- (Optional, for direct development) Python / Node.js as required by the backend and frontend

### DB configs
You need to add a db configs to metadata_app/backend/conf. Add read user everywhere exept where it specifies write in the name
```bash
{
    "meta": {
        "host": "mysql-ens-genebuild-prod-1",
        "user": "",
        "password": "",
        "database": "gb_assembly_metadata",
        "port": },

    "meta_write": {
        "host": "mysql-ens-genebuild-prod-1",
        "user": "",
        "password": "",
        "database": "gb_assembly_metadata",
        "port": },

    "transcriptomic": {
        "host": "mysql-ens-genebuild-prod-1",
        "user": "",
        "password": "",
        "database": "gb_transcriptomic_registry",
        "port": 
    },

    "prod": {
    "host": "mysql-ens-production-1",
    "user": "",
    "password": "",
    "database": "ensembl_genome_metadata",
    "port": 
  },

  "main": {
    "host": "mysql-ens-meta-prod-1",
    "user": "",
    "password": "",
    "database": "ensembl_metadata",
    "port":  }
}

```




## Local Development (Without Docker)

You can run the backend and frontend directly on your machine.

### Backend

```bash
cd backend
# Create and activate a virtual environment, then install dependencies
# python -m venv .venv
# source .venv/bin/activate
# pip install -r requirements.txt

# Run the backend server
# uvicorn app.main:app --reload
```

### Frontend

```bash
cd frontend
# Install dependencies
# npm install

# Start the development server
# npm run dev 
```

### Docker

Build from the repository root with the app Dockerfile under `metadata_app/`:

```bash
docker build -f metadata_app/Dockerfile .
```
