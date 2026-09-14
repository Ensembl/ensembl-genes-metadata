# app/core/database.py
import json
import logging
import os
from contextlib import contextmanager
from pathlib import Path

import pymysql
from pymysql.cursors import DictCursor


BACKEND_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DB_CONFIG_PATH = BACKEND_ROOT / "conf" / "db_config.json"
LOG_DIR = BACKEND_ROOT / "logs"
LOG_FILE = LOG_DIR / "app.log"


def setup_logging():
    """Configure logging for the application."""
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(
        level=logging.INFO,
        format=log_format,
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(LOG_FILE),
        ],
    )
    return logging.getLogger(__name__)


def load_db_config():
    """Load database configuration from file or environment."""
    config_path = Path(
        os.environ.get("DB_CONFIG_PATH", str(DEFAULT_DB_CONFIG_PATH))
    ).expanduser()

    if not config_path.exists():
        logging.error("Config file not found: %s", config_path)
        raise FileNotFoundError(f"Database config file '{config_path}' not found.")

    with config_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


@contextmanager
def get_db_connection(config_key):
    """Context manager for database connections."""
    db_config = load_db_config()

    if config_key not in db_config:
        raise KeyError(
            f"Database config key '{config_key}' not found. Available keys: {list(db_config.keys())}"
        )

    connection = None
    try:
        connection = pymysql.connect(
            **db_config[config_key],
            cursorclass=DictCursor,
        )
        yield connection
    except pymysql.Error as exc:
        logging.error("Database connection error: %s", exc)
        raise
    finally:
        if connection:
            connection.close()
