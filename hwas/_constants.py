"""Default configuration parameters and environment variables"""

import os
import re
from . import __version__


VERSION = __version__.version


GITHUB_URL = "https://github.com/Palmer-Lab-UCSD/hwas"

DEFAULT_META_PREFIX = "##"
DEFAULT_HEADER_PREFIX = "#"
DEFAULT_DIRMODE = 0o770

OUTPUT_DELIMITER = ','



DEFAULT_DB_USER = None
if "USER" in os.environ:
    DEFAULT_DB_USER = os.environ["USER"]

DEFAULT_LOG_DIR = "logs"


# environmental variables to check value
ENV_BIN = "PALMER_BIN"
ENV_SLURM_QOS = "PALMER_QOS"
ENV_SLURM_ACCOUNT = "PALMER_ACCOUNT"
ENV_SLURM_PARTITION = "PALMER_PARTITION"
ENV_DB_NAME = "PALMER_DB_NAME"
ENV_DB_USERNAME = "PALMER_DB_USERNAME"
ENV_DB_HOST = "PALMER_DB_HOST"
ENV_DB_PORT = "PALMER_DB_PORT"
ENV_DB_PW = "PALMER_DB_PW"


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


TEMPLATES_PATH = "hwas.templates"

TEMPLATES_CONFIG = "config.template"
FILENAME_CONFIG = "config"

TEMPLATES_HWAS_SBATCH="hwas.sh.template"
FILENAME_HWAS_SBATCH="hwas.sh"


FILENAME_PHENOTYPE="phenotype.csv"
FILENAME_COVARIATES="covariates.csv"
