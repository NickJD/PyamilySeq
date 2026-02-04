import logging
import os
import sys
import argparse
from datetime import datetime
from io import StringIO
import re

PyamilySeq_Version = 'v1.3.3'
WELCOME = f"Thank you for using PyamilySeq {PyamilySeq_Version} - A tool for gene clustering and pangenome analysis."
CITATION = "Please Cite PyamilySeq: https://doi.org/10.1093/nargab/lqaf198"
ISSUE = "Please report any issues to: https://github.com/NickJD/PyamilySeq/issues"

def configure_logger(logger_name, enable_file=False, log_dir=None, level=logging.INFO, verbose=False):
    """
    Create and return a configured logger.
    - logger_name: full logger name (e.g. "PyamilySeq.Group_Splitter")
    - enable_file: if True, create a timestamped logfile in log_dir
    - log_dir: directory for logfile (defaults to cwd)
    - level: console log level (default INFO)
    - verbose: if True, sets console level to DEBUG and file to DEBUG
    """
    logger = logging.getLogger(logger_name)
    # Clear previous handlers to avoid duplicate logs on repeated imports/runs
    if logger.hasHandlers():
        logger.handlers.clear()

    # Determine levels
    console_level = logging.DEBUG if verbose else level
    logger.setLevel(logging.DEBUG if verbose else level)

    # Formatter without logger name (keeps output clean)
    formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")

    # Console handler -> write to stdout by default
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(console_level)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Optional file handler
    file_handler = None
    if enable_file:
        if not log_dir:
            log_dir = os.getcwd()
        os.makedirs(log_dir, exist_ok=True)
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        # Use short tool name derived from logger_name for the filename
        safe_name = logger_name.split('.')[-1]
        file_name = f"{safe_name}-{ts}.log"
        fh = logging.FileHandler(os.path.join(log_dir, file_name))
        fh.setLevel(logging.DEBUG)  # file always capture debug for diagnostics
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        file_handler = fh
        logger.debug("File logging enabled: %s", os.path.join(log_dir, file_name))

    # Standard startup banner for all tools (printed once per logger instance)
    # If banner hasn't been printed at all, log it normally (will go to console and file if present).
    if not getattr(logger, "_welcome_printed", False):
        logger.info("%s", WELCOME)
        logger.info("%s", CITATION)
        logger.info("%s", ISSUE)
        setattr(logger, "_welcome_printed", True)
        # Mark that banner also written to file if file handler exists
        if file_handler:
            # Also write formatted lines directly into the file to guarantee presence
            try:
                for msg in (WELCOME, CITATION, ISSUE):
                    rec = logging.LogRecord(logger.name, logging.INFO, "", 0, msg, None, None)
                    formatted = formatter.format(rec)
                    # write to file handler's stream and flush
                    try:
                        file_handler.stream.write(formatted + "\n")
                        file_handler.stream.flush()
                    except Exception:
                        # Best-effort; ignore write errors
                        pass
                setattr(logger, "_welcome_file_written", True)
            except Exception:
                pass
    else:
        # Banner already printed (likely to console by an early logger). Ensure it is written to file
        # if file logging was just enabled and it hasn't yet been written to file.
        if file_handler and not getattr(logger, "_welcome_file_written", False):
            # Write banner lines directly into the file handler's stream (avoid duplicating console output).
            try:
                for msg in (WELCOME, CITATION, ISSUE):
                    rec = logging.LogRecord(logger.name, logging.INFO, "", 0, msg, None, None)
                    formatted = formatter.format(rec)
                    try:
                        file_handler.stream.write(formatted + "\n")
                        file_handler.stream.flush()
                    except Exception:
                        pass
                setattr(logger, "_welcome_file_written", True)
            except Exception:
                pass

    return logger

#  ArgumentParser subclass that logs usage/help/errors via the logger
class LoggingArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args, logger_name=None, **kwargs):
        super().__init__(*args, **kwargs)
        # If logger_name provided, use that logger; otherwise use root logger
        self._logger = logging.getLogger(logger_name) if logger_name else logging.getLogger()
        # Emit the parser description immediately on creation so it appears for normal runs
        # (tools create an early console-only logger before constructing the parser).
        if getattr(self, 'description', None):
            try:
                self._logger.info("%s", str(self.description))
            except Exception:
                # If logging fails for any reason, swallow to avoid breaking parser creation.
                pass

    def print_usage(self, file=None):
        # Preserve default usage printing to console; description already logged at init.
        super().print_usage(file)

    def print_help(self, file=None):
        # Capture help output, strip description (already logged), and print the rest to console.
        sio = StringIO()
        super().print_help(sio)
        help_text = sio.getvalue()
        if self.description:
            pattern = re.escape(str(self.description)) + r'(\r?\n){1,2}'
            help_text = re.sub(pattern, '', help_text, count=1)
        out_file = file if file is not None else sys.stdout
        out_file.write(help_text)

    def exit(self, status=0, message=None):
        # Preserve argparse behaviour by writing any exit message to stderr and exiting.
        if message:
            sys.stderr.write(message)
        raise SystemExit(status)

    def error(self, message):
        # Print usage to stderr (as argparse does) and log a concise error message via logger.
        super().print_usage(sys.stderr)
        prog = self.prog if hasattr(self, 'prog') else ''
        self._logger.error("%s: error: %s", prog, message)
        self.exit(2)
