"""Logging utilities."""

import logging
import sys


def setup_logging(level: str = "INFO") -> logging.Logger:
    """
    Setup logging configuration.
    
    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR)
        
    Returns:
        Configured logger
    """
    logger = logging.getLogger("crispr_offtarget")
    logger.setLevel(getattr(logging, level.upper()))
    
    # Clear existing handlers
    logger.handlers.clear()
    
    # Console handler
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(getattr(logging, level.upper()))
    
    # Format
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    handler.setFormatter(formatter)
    
    logger.addHandler(handler)
    
    return logger
