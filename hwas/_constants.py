"""Default configuration parameters

By: Robert Vogel
Affiliation: Palmer Lab at UCSD

"""

import os
from . import __version__


DEFAULT_META_PREFIX = "##"
DEFAULT_CONFIG_FILENAME = "config"
DEFAULT_DIRMODE = 0o770

DEFAULT_DB_HOST = "localhost"
DEFAULT_DB_PORT = "5432"
DEFAULT_DB_NAME = "hsrats"
DEFAULT_DB_PASSWORD_ENV_VAR = "PALMER_LAB_DB_PASSWORD"

DEFAULT_DB_USER = None
if "USER" in os.environ:
    DEFAULT_DB_USER = os.environ["USER"]

DEFAULT_LOG_DIR = "logs"


VERSION = __version__.version


