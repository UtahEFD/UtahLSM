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
"""A helper module for configuring logging.

This module provides a two-stage logging setup. It allows modules to
start logging immediately upon import, buffering the messages in memory.
Once the main configuration is loaded, `finalize_logging` is called to
set up the final file and console handlers and flush all buffered messages.
"""
from typing import Optional, Dict
from pathlib import Path
import logging
import logging.handlers

_buffer_handler: Optional[logging.handlers.MemoryHandler] = None
_initialized: bool = False

def _ensure_buffered() -> None:
    """Starts buffered logging if it is not already active.

    This internal function sets up a `MemoryHandler` on the root logger
    to capture all log messages generated before the main logging
    configuration is finalized. This ensures no messages are lost during
    the initial setup phase of the model.
    """
    global _buffer_handler, _initialized  # pylint: disable=global-statement
    if _initialized:
        return
    root_logger = logging.getLogger()
    if not root_logger.handlers:
        _buffer_handler = logging.handlers.MemoryHandler(capacity=100)
        root_logger.addHandler(_buffer_handler)
        root_logger.setLevel(logging.DEBUG)
    _initialized = True

def get_logger(name: Optional[str] = None) -> logging.Logger:
    """Gets a logger instance and ensures buffering is active.

    This is the main function that should be called by other modules to get
    a logger. It guarantees that the buffering handler is in place before
    returning the logger.

    Args:
        name: The name of the logger, typically `__name__`. If None, returns
            the root logger [Optional[str]].

    Returns:
        A `logging.Logger` instance with the specified name or the root logger
        if name is None.
    """
    _ensure_buffered()
    return logging.getLogger(name)

def finalize_logging(level_str: str = 'info') -> None:
    """Replaces the buffer with final handlers and flushes stored logs.

    This function should be called once after the main configuration has been
    read. It removes the temporary memory handler and replaces it with
    configured file and console handlers. It then flushes any messages that
    were buffered during startup to the new handlers.

    Args:
        level_str: The desired logging level as a string (e.g., 'info',
            'debug'). Defaults to 'info'.
    """
    global _buffer_handler  # pylint: disable=global-statement
    root_logger = logging.getLogger()

    log_levels: Dict[str, int] = {'info': logging.INFO, 'debug': logging.DEBUG}
    log_level = log_levels.get(level_str.lower(), logging.INFO)
    log_format = logging.Formatter(
        '{asctime} [{levelname:^8s}] {name:.>10s}: {message}',
        datefmt='%Y-%m-%d %H:%M:%S', style='{'
    )

    # Create logs directory in the python folder (go up 3 levels from this file)
    log_dir: Path = Path(__file__).resolve().parents[3] / 'logs'
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file: Path = log_dir / 'utahlsm.log'

    # file handler
    file_handler = logging.FileHandler(log_file, mode='w')
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
