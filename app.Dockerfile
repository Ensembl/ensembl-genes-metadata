# Use Python base image
FROM python:3.10-slim

# Set working directory
WORKDIR /ensembl-genes-metadata

# Copy backend
COPY metadata_app/backend ./metadata_app/backend/

# Copy logging
COPY logs/ ./logs

# Install backend dependencies
RUN pip install fastapi uvicorn
COPY metadata_app/backend/requirements.txt .
RUN pip install -r requirements.txt

# Install Node.js for frontend
RUN apt-get update && \
    apt-get install -y curl && \
    curl -fsSL https://deb.nodesource.com/setup_18.x | bash - && \
    apt-get install -y nodejs && \
    apt-get clean

# Copy frontend and install dependencies
COPY metadata_app/frontend/ ./metadata_app/frontend/
WORKDIR metadata_app/frontend
RUN npm install && npm run build

# Move back to main workdir
WORKDIR /ensembl-genes-metadata

# Expose ports
EXPOSE 8000 3000

# Start both frontend and backend
CMD ["sh", "-c", "uvicorn metadata_app.backend.app.main:app --host 0.0.0.0 --port 8000 & npm --prefix metadata_app/frontend start"]