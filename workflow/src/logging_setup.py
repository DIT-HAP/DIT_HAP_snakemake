"""Logging configuration for DIT-HAP scripts."""

# =============================================================================
# IMPORTS
# =============================================================================
import sys

from loguru import logger

# =============================================================================
# CORE LOGIC
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )
