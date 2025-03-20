"""Default configuration parameters

By: Robert Vogel
Affiliation: Palmer Lab at UCSD

"""

import os
import re
from . import __version__


VERSION = __version__.version



DEFAULT_CONFIG_PATH = "hwas.templates"
DEFAULT_META_PREFIX = "##"
DEFAULT_HEADER_PREFIX = "#"
DEFAULT_CONFIG_FILENAME = "config"
DEFAULT_DIRMODE = 0o770

OUTPUT_DELIMITER = ','

# DEFAULT_DB_HOST = "localhost"
# DEFAULT_DB_PORT = "5432"
# DEFAULT_DB_NAME = "hsrats"


DEFAULT_DB_USER = None
if "USER" in os.environ:
    DEFAULT_DB_USER = os.environ["USER"]

DEFAULT_LOG_DIR = "logs"


# environmental variables to check value
ENV_DB_PASSWORD = "PALMER_LAB_DB_PASSWORD"
ENV_BIN = "PALMER_BIN"
ENV_DB_USERNAME = "PALMER_DB_USERNAME"
ENV_DB_HOST = "PALMER_DB_HOST"
ENV_DB_PORT = "PALMER_DB_PORT"
ENV_DB_NAME = "PALMER_DB_NAME"


# Palmer lab database specific parameters
IS_COVARIATE = re.compile("covariate")
PHENOTYPE_TABLENAME = "gwas_phenotypes"
METADATA_TABLENAME = "descriptions"

SAMPLE_COLNAME = "rfid"

COVARIATE_DELIMITER = ','
COVARIATE_FILE_SUFFIX = "_covariates.csv"
COVARIATE_TYPE_TOKEN = '%covariate%'

PHENOTYPE_FILE_SUFFIX = ".csv"

MEASURE_TYPE_COLNAME = 'trait_covariate'
MEASURE_NAME_COLNAME = "measure"

