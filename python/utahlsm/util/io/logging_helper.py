# 
# UtahLSM
# 
# Copyright (c) 2017–2025 Jeremy A. Gibbs
# Copyright (c) 2017–2025 Rob Stoll
# Copyright (c) 2017–2025 Eric Pardyjak
# Copyright (c) 2017–2025 Pete Willemsen
# 
# This file is part of UtahLSM.
# 
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
# 

# logging_helper.py
import logging
import logging.handlers

_buffer_handler = None
_initialized = False

def _ensure_buffered():
    """Start buffered logging if not already active."""
    global _buffer_handler, _initialized
    if _initialized:
        return
    root_logger = logging.getLogger()
    if not root_logger.handlers:
        _buffer_handler = logging.handlers.MemoryHandler(capacity=100)
        root_logger.addHandler(_buffer_handler)
        root_logger.setLevel(logging.DEBUG)
    _initialized = True

def get_logger(name: str = None):
    """
    Ensures buffered logging is active if no config yet.
    """
    _ensure_buffered()
    return logging.getLogger(name)

def finalize_logging(level_str: str = "info"):
    """Replace buffer with real handlers and flush stored logs."""
    global _buffer_handler
    root_logger = logging.getLogger()

    LOG_LEVELS = {"info": logging.INFO, "debug": logging.DEBUG}
    log_level  = LOG_LEVELS.get(level_str.lower(), logging.INFO)
    filler = "_"
    log_format = logging.Formatter(
        "{asctime} [{levelname:^8s}] {name:.>10s}: {message}",
        datefmt="%Y-%m-%d %H:%M:%S", style="{"
    )

    # file handlers
    file_handler = logging.FileHandler("utahlsm.log", mode="w")
    file_handler.setLevel(log_level)
    file_handler.setFormatter(log_format)
    
    # console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(log_format)

    # reset root handlers
    root_logger.handlers = []
    root_logger.setLevel(log_level)
    root_logger.addHandler(file_handler)
    root_logger.addHandler(console_handler)

    # flush buffered messages
    if _buffer_handler is not None:
        _buffer_handler.setTarget(root_logger)
        _buffer_handler.flush()
        root_logger.removeHandler(_buffer_handler)
        _buffer_handler = None
