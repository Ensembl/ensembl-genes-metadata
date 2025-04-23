# app/main.py
from fastapi import FastAPI, Depends, HTTPException, Query
from typing import List, Optional, Dict, Any
import pandas as pd
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
import logging

from metadata_app.backend.app.api.routes import assemblies, taxonomy, transcriptomics
from metadata_app.backend.app.core.database import setup_logging

# Initialize FastAPI app
app = FastAPI(
    title="Genebuild Metadata API",
    description="API for querying and filtering genebuild metadata",
    version="1.0.0"
)

# Set up CORS for frontend access
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Set to specific origins in production !!!!!!
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Set up logging
setup_logging()

# Include routers from different modules
app.include_router(assemblies.router, prefix="/api/assemblies", tags=["assemblies"])
app.include_router(taxonomy.taxonomy, prefix="/api/taxonomy", tags=["taxonomy"])
app.include_router(transcriptomics.transcriptomics, prefix="/api/transcriptomics", tags=["transcriptomics"])

@app.get("/")
async def root():
    return {"message": "Welcome to the Genebuild API"}

@app.get("/health")
async def health_check():
    return {"status": "healthy"}

# Error handlers
@app.exception_handler(HTTPException)
async def http_exception_handler(request, exc):
    return JSONResponse(
        status_code=exc.status_code,
        content={"message": exc.detail},
    )

@app.exception_handler(Exception)
async def general_exception_handler(request, exc):
    logging.error(f"Unhandled exception: {str(exc)}")
    return JSONResponse(
        status_code=500,
        content={"message": "Internal server error"},
    )