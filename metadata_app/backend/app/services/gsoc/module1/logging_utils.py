"""
logging_utils.py
Single shared helper for obtaining a module-level logger.

Previously, individual Module 1 files configured logging inconsistently -
some called `logging.getLogger(__name__)` (correct, named logger), others
called the bare `logging.info(...)` / `logging.error(...)` module-level
functions directly, which silently uses the root logger and loses the
per-module name in log output. This module gives every file a single,
consistent way to obtain a properly named logger.
"""

import logging


def get_logger(name: str) -> logging.Logger:
    """
    Return a module-level logger for the given name.

    Args:
        name: Usually the caller's `__name__`, so log records are
            attributed to the module that emitted them.

    Returns:
        A standard library `logging.Logger` instance.
    """
    return logging.getLogger(name)
