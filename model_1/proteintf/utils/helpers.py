import os
import sys
import copy
import random
import logging
import numpy as np

import tensorflow as tf

from proteintf.common.tee_logger import TeeLogger


def with_fallback(preferred, fallback):
    """
    Code taken from https://github.com/allenai/allennlp/blob/master/allennlp/common/params.py#L121

    Deep merge two dicts, preferring values from `preferred`.
    """
    preferred_keys = set(preferred.keys())
    fallback_keys = set(fallback.keys())
    common_keys = preferred_keys & fallback_keys

    merged = {}

    for key in preferred_keys - fallback_keys:
        merged[key] = copy.deepcopy(preferred[key])
    for key in fallback_keys - preferred_keys:
        merged[key] = copy.deepcopy(fallback[key])

    for key in common_keys:
        preferred_value = preferred[key]
        fallback_value = fallback[key]

        if isinstance(preferred_value, dict) and isinstance(fallback_value, dict):
            merged[key] = with_fallback(preferred_value, fallback_value)
        else:
            merged[key] = copy.deepcopy(preferred_value)

    return merged


def prepare_environment(seed):
    random.seed(seed)
    np.random.seed(seed)
    tf.set_random_seed(seed)


def prepare_global_logging(serialization_dir: str, log_prefix='std') -> logging.FileHandler:
    """
    Code taken from https://github.com/allenai/allennlp/blob/master/allennlp/common/util.py#L209

    This function configures 3 global logging attributes - streaming stdout and stderr
    to a file as well as the terminal, setting the formatting for the python logging
    library and setting the interval frequency for the Tqdm progress bar.
    Note that this function does not set the logging level, which is set in ``allennlp/run.py``.
    Parameters
    ----------
    serialization_dir : ``str``, required.
        The directory to stream logs to.
    Returns
    -------
    ``logging.FileHandler``
        A logging file handler that can later be closed and removed from the global logger.
    """
    std_out_file = os.path.join(serialization_dir, log_prefix + "out.log")
    sys.stdout = TeeLogger(std_out_file,    # type: ignore
                           sys.stdout,
                           file_friendly_terminal_output=False)
    sys.stderr = TeeLogger(os.path.join(serialization_dir, log_prefix + "stderr.log"),  # type: ignore
                           sys.stderr,
                           file_friendly_terminal_output=False)

    stdout_handler = logging.FileHandler(std_out_file)
    stdout_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s'))
    logging.getLogger().addHandler(stdout_handler)

    return stdout_handler


def cleanup_global_logging(stdout_handler: logging.FileHandler) -> None:
    """
    Code taken from https://github.com/allenai/allennlp/blob/master/allennlp/common/util.py#L253

    This function closes any open file handles and logs set up by `prepare_global_logging`.
    Parameters
    ----------
    stdout_handler : ``logging.FileHandler``, required.
        The file handler returned from `prepare_global_logging`, attached to the global logger.
    """
    stdout_handler.close()
    logging.getLogger().removeHandler(stdout_handler)

    if isinstance(sys.stdout, TeeLogger):
        sys.stdout = sys.stdout.cleanup()
    if isinstance(sys.stderr, TeeLogger):
        sys.stderr = sys.stderr.cleanup()
