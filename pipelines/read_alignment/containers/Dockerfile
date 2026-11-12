FROM python:3.10-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        wget \
        gzip \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    PyMySQL==1.2.3 \
    SQLAlchemy==2.0.54 \
    pandas==2.3.3 \
    requests==2.34.2

WORKDIR /app

CMD ["python"]