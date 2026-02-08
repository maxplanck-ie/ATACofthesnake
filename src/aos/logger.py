import logging
from rich.logging import RichHandler
from pathlib import Path

def setup_logger(logfile: Path) -> logging.Logger:
    _logger = logging.getLogger()
    rich_handler = RichHandler(rich_tracebacks=True, show_time=False, show_level=True, show_path=False)
    _fmt = logging.Formatter('%(levelname)s - %(asctime)s - %(message)s', datefmt="%H:%M:%S")
    file_handler = logging.FileHandler(logfile, mode='w')
    
    file_handler.setFormatter(_fmt)
    _logger.setLevel(logging.DEBUG)
    file_handler.setLevel(logging.DEBUG)
    rich_handler.setLevel(logging.DEBUG)
        

    _logger.addHandler(rich_handler)
    _logger.addHandler(file_handler)
    return _logger