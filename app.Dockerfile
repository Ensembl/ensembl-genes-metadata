# Use Python base image
FROM python:3.10-slim

# Set working directory
WORKDIR /ensembl-genes-metadata


# Install Node.js for frontend build
RUN apt-get update && \
    apt-get install -y curl && \
    curl -fsSL https://deb.nodesource.com/setup_18.x | bash - && \
    apt-get install -y nodejs && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Copy and install backend dependencies first
COPY metadata_app/backend/requirements.txt ./
RUN pip install fastapi uvicorn
RUN pip install -r requirements.txt

# Copy backend code
COPY metadata_app/backend ./metadata_app/backend/

# Copy logging
COPY logs/ ./logs

# Copy frontend code and build it
COPY metadata_app/frontend/ ./metadata_app/frontend/
WORKDIR /ensembl-genes-metadata/metadata_app/frontend

# Make sure you have the correct next.config.js for static export
COPY metadata_app/frontend/next.config.ts ./

# Install frontend dependencies and build static export
RUN npm install
RUN npm run build

# Verify the build worked
RUN ls -la out/

# Move back to main workdir
WORKDIR /ensembl-genes-metadata

# Only expose port 8000 (FastAPI will serve everything)
EXPOSE 8000

# Start only FastAPI (it will serve both API and frontend)
CMD ["uvicorn", "metadata_app.backend.app.main:app", "--host", "0.0.0.0", "--port", "8000"]