ARG BUILDPLATFORM

FROM --platform=$BUILDPLATFORM node:20-bookworm-slim AS frontend-build

WORKDIR /app/metadata_app/frontend

COPY metadata_app/frontend/package.json metadata_app/frontend/package-lock.json ./
ENV NEXT_TELEMETRY_DISABLED=1
RUN npm ci

COPY metadata_app/frontend/ ./
RUN npm run build


FROM python:3.10-slim

WORKDIR /ensembl-genes-metadata

COPY metadata_app/backend/requirements.txt ./metadata_app/backend/requirements.txt
RUN pip install --no-cache-dir -r metadata_app/backend/requirements.txt

COPY metadata_app/backend ./metadata_app/backend/
RUN mkdir -p ./metadata_app/frontend
COPY --from=frontend-build /app/metadata_app/frontend/out ./metadata_app/frontend/out

EXPOSE 8000

CMD ["uvicorn", "metadata_app.backend.app.main:app", "--host", "0.0.0.0", "--port", "8000"]
