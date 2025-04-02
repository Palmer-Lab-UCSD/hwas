"""Default configuration parameters and environment variables"""

import os
import re
from . import __version__


MAX_PHENOTYPE_VERSIONS = 1000

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



# COMMON 
VERSION = __version__.version
GITHUB_URL = "https://github.com/Palmer-Lab-UCSD/hwas"
FILENAME_PHENOTYPE="phenotype.csv"
FILENAME_COVARIATES="covariates.csv"
DEFAULT_LOG_DIR = "logs"


# Directory setup, file information 
DEFAULT_META_PREFIX = "##"
DEFAULT_HEADER_PREFIX = "#"
DEFAULT_DIRMODE = 0o770

OUTPUT_DELIMITER = ','



# VCF file pars
CONTIG_REGEX  = "chr[0-9]+"


# Palmer lab database specific parameters
DEFAULT_DB_USER = os.environ.get("USER")


PHENOTYPE_TABLENAME = "gwas_phenotypes"
METADATA_TABLENAME = "descriptions"
SAMPLE_COLNAME = "rfid"
COVARIATE_DELIMITER = ','
MEASURE_TYPE_COLNAME = 'trait_covariate'
MEASURE_NAME_COLNAME = "measure"

COVARIATE_TYPE_TOKEN = '%covariate%'
IS_COVARIATE = re.compile("covariate")


# Templates
TEMPLATES_PATH = "hwas.templates"

TEMPLATES_CONFIG = "config.template"
FILENAME_CONFIG = "config"

TEMPLATES_HWAS_SBATCH="hwas.sh.template"
FILENAME_HWAS_SBATCH="hwas.sh"

TEMPLATES_HGRM_SBATCH="hgrm.sh.template"
FILENAME_HWAS_SBATCH="hgrm.sh"



# hgrm computation

HGRM_CPUS_PER_TASK="1"
HGRM_ALLOC_TIME="4:00:00"
HGRM_CONTIG_FILE="contigs"
