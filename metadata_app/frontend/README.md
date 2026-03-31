# Frontend for the genebuild metadat app

Run the development server:

```bash
npm run dev
```

To build the static files (accessible via the beckend)
```bash
npm run build
```

To run with backend 
```bash
uvicorn metadata_app.backend.app.main:app --reload
```