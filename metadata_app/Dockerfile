FROM python:3.10-slim

WORKDIR /ensembl-genes-metadata

RUN apt-get update && \
    apt-get install -y curl && \
    curl -fsSL https://deb.nodesource.com/setup_22.x | bash - && \
    apt-get install -y nodejs && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY metadata_app/backend/requirements.txt ./metadata_app/backend/requirements.txt
RUN pip install --no-cache-dir -r metadata_app/backend/requirements.txt

COPY metadata_app/backend ./metadata_app/backend/
COPY metadata_app/frontend ./metadata_app/frontend/

WORKDIR /ensembl-genes-metadata/metadata_app/frontend
RUN npm install -g npm@latest
RUN npm install
RUN npm run build

WORKDIR /ensembl-genes-metadata

EXPOSE 8000

CMD ["uvicorn", "metadata_app.backend.app.main:app", "--host", "0.0.0.0", "--port", "8000"]
