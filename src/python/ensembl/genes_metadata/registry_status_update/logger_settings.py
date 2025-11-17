import logging
from datetime import datetime


def get_logger(name: str = __name__) -> logging.Logger:
	"""Configure and return a logger with today's date in the filename."""

	# Generate filename with current date
	today_str = datetime.now().strftime("%Y-%m-%d")
	log_filename = f"genebuild_status_updater_{today_str}.log"

	# Create logger
	logger = logging.getLogger(name)
	logger.setLevel(logging.INFO)

	# Prevent adding multiple handlers if called multiple times
	if not logger.handlers:
		# File handler
		fh = logging.FileHandler(log_filename)
		fh.setLevel(logging.INFO)
		fh_formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
		fh.setFormatter(fh_formatter)
		logger.addHandler(fh)

		# Stream handler (console)
		sh = logging.StreamHandler()
		sh.setLevel(logging.INFO)
		sh.setFormatter(fh_formatter)
		logger.addHandler(sh)

	logger.info(f"Logging started. Writing to {log_filename}")
	return logger