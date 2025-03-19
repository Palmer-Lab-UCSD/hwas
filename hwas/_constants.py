"""Default configuration parameters

By: Robert Vogel
Affiliation: Palmer Lab at UCSD

"""

import os
import re
from . import __version__


VERSION = __version__.version

DEFAULT_META_PREFIX = "##"
DEFAULT_HEADER_PREFIX = "#"
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



IS_COVARIATE = re.compile("covariate")
PHENOTYPE_TABLENAME = "gwas_phenotypes"
METADATA_TABLENAME = "descriptions"

SAMPLE_COLNAME = "rfid"

COVARIATE_DELIMITER = ','
COVARIATE_TYPE_TOKEN = '%covariate%'

MEASURE_TYPE_COLNAME = 'trait_covariate'
MEASURE_NAME_COLNAME = "measure"

OUTPUT_DELIMITER = ','
