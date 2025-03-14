"""Default configuration parameters

By: Robert Vogel
Affiliation: Palmer Lab at UCSD

"""

import os
import re

DEFAULT_DB = ("palmerlab-main-database-c2021-08-02.c6sgfwysomht"
                ".us-west-2.rds.amazonaws.com")

DEFAULT_DB_NAME = "PalmerLab_Datasets"
DEFAULT_DB_USER = "postgres"
DEFAULT_DB_PORT = "5432"


OUTPUT_META_PREFIX = "##"

IS_COVARIATE = re.compile("covariate")

with open(os.path.join(os.path.abspath(__file__),"VERSION"), "r") as fid:
    VERSION = fd.readline()
