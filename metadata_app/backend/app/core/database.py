# app/core/database.py
import os
import json
import logging
import pymysql
from pymysql.cursors import DictCursor
from contextlib import contextmanager


def setup_logging():
	"""Configure logging for the application"""
	log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
	logging.basicConfig(
		level=logging.INFO,
		format=log_format,
		handlers=[
			logging.StreamHandler(),
			logging.FileHandler("app.log")
		]
	)
	return logging.getLogger(__name__)


def load_db_config():
	"""Load database configuration from file or environment"""
	config_path = os.environ.get("DB_CONFIG_PATH", "conf/db_config.json")

	if not os.path.exists(config_path):
		logging.error(f"Config file not found: {config_path}")
		raise FileNotFoundError(f"Database config file '{config_path}' not found.")

	with open(config_path, "r") as f:
		return json.load(f)


@contextmanager
def get_db_connection(config_key):
	"""Context manager for database connections"""
	db_config = load_db_config()

	if config_key not in db_config:
		raise KeyError(f"Database config key '{config_key}' not found")

	connection = None
	try:
		connection = pymysql.connect(
			**db_config[config_key],
			cursorclass=DictCursor
		)
		yield connection
	except pymysql.Error as e:
		logging.error(f"Database connection error: {e}")
		raise
	finally:
		if connection:
			connection.close()