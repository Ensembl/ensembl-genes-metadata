# app/main.py
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, FileResponse
from fastapi.staticfiles import StaticFiles
import logging
from pathlib import Path

from metadata_app.backend.app.api.routes import assemblies, taxonomy, transcriptomics, home_page, annotations, \
	report_annotations, report_assemblies, bioproject_search, taxonomy_search, project
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
app.include_router(annotations.annotations, prefix="/api/annotations", tags=["annotations"])
app.include_router(taxonomy.taxonomy, prefix="/api/taxonomy", tags=["taxonomy"])
app.include_router(transcriptomics.transcriptomics, prefix="/api/transcriptomics", tags=["transcriptomics"])
app.include_router(home_page.home_page, prefix="/api/home_page", tags=["home_page"])
app.include_router(report_annotations.report, prefix="/api/report/anno", tags=["report_anno"])
app.include_router(report_assemblies.report, prefix="/api/report/asm", tags=["report_asm"])
app.include_router(bioproject_search.router, prefix="/api/bioproject_search", tags=["bioproject_search"])
app.include_router(taxonomy_search.router, prefix="/api/taxonomy_search", tags=["taxonomy_search"])
app.include_router(project.project_router, prefix="/api/project", tags=["project"])


# API health check endpoints
@app.get("/api/")
async def api_root():
	return {"message": "Welcome to the Genebuild API"}


@app.get("/api/health")
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


# Next.js frontend serving
# Next.js static export creates an 'out' directory
REACT_BUILD_PATH = Path("metadata_app/frontend/out")  # Next.js static export directory

# Check if React build directory exists
if REACT_BUILD_PATH.exists():
	# Mount Next.js static files (CSS, JS, images, etc.)
	# Next.js uses _next directory for static assets, not static
	nextjs_static_path = REACT_BUILD_PATH / "_next"
	if nextjs_static_path.exists():
		app.mount("/_next", StaticFiles(directory=str(nextjs_static_path)), name="nextjs_static")


	# Serve common static assets
	@app.get("/favicon.ico")
	async def favicon():
		favicon_path = REACT_BUILD_PATH / "favicon.ico"
		if favicon_path.exists():
			return FileResponse(str(favicon_path))
		raise HTTPException(status_code=404, detail="Favicon not found")


	@app.get("/manifest.json")
	async def manifest():
		manifest_path = REACT_BUILD_PATH / "manifest.json"
		if manifest_path.exists():
			return FileResponse(str(manifest_path))
		raise HTTPException(status_code=404, detail="Manifest not found")


	@app.get("/robots.txt")
	async def robots():
		robots_path = REACT_BUILD_PATH / "robots.txt"
		if robots_path.exists():
			return FileResponse(str(robots_path))
		raise HTTPException(status_code=404, detail="Robots.txt not found")


	# Root endpoint - serve Next.js app
	@app.get("/")
	async def serve_react_root():
		index_file = REACT_BUILD_PATH / "index.html"
		if index_file.exists():
			return FileResponse(str(index_file))
		raise HTTPException(status_code=404, detail="Next.js app not found")


	# Catch-all handler: serve Next.js app for all non-API routes
	@app.get("/{full_path:path}")
	async def serve_react_app(full_path: str):
		# Don't serve Next.js app for API routes or Next.js static assets
		if full_path.startswith("api/") or full_path.startswith("_next/"):
			raise HTTPException(status_code=404, detail="Resource not found")

		# Check if there's an HTML file for this route (Next.js static export)
		html_file = REACT_BUILD_PATH / f"{full_path}.html"
		if html_file.exists():
			return FileResponse(str(html_file))

		# Check if there's an index.html in a directory for this route
		dir_index = REACT_BUILD_PATH / full_path / "index.html"
		if dir_index.exists():
			return FileResponse(str(dir_index))

		# Fallback to main index.html for client-side routing
		index_file = REACT_BUILD_PATH / "index.html"
		if index_file.exists():
			return FileResponse(str(index_file))
		else:
			raise HTTPException(status_code=404, detail="Next.js app not found")
else:
	logging.warning("Next.js build directory not found. Make sure to run 'npm run build' first.")


	# Fallback root endpoint if Next.js build doesn't exist
	@app.get("/")
	async def root():
		return {"message": "Welcome to the Genebuild API - Next.js frontend not built yet"}