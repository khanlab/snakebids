import logging
from typing import Union


def setup_logging(level: Union[str, int] = logging.DEBUG):
    logger = logging.getLogger("snakebids")
    stream_handler = logging.StreamHandler()
    formatter = logging.Formatter("%(levelname)-8s %(message)s")
    stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.setLevel(level)
