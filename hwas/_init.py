"""Initialize parameter settings and directory structure.




Functions:
    interface: interface with command line to set up directories
        and configuration file in the required format
"""

import os
import logging


from . import _constants
from . import _config

logger = logging.getLogger(__name__)


def interface(schema_name: str,
              phenotype: str) -> None:

    
    # Use template to initialize configuration file.  A subset of
    # these values will be initialized to those sepecied by environmental
    # variables
    cfg = _config.init()


    # set up path and directories for phenotype and
    # project (schema_name)
    schema_path = os.path.join(os.path.expanduser(os.getcwd()),
                               schema_name)

    if not os.path.isdir(schema_path):
        os.mkdir(schema_name, mode = _constants.DEFAULT_DIRMODE)

    path = os.path.join(schema_path, phenotype)
    if not os.path.isdir(path):
        os.mkdir(path, mode = _constants.DEFAULT_DIRMODE)


    for i in range(1, _constants.MAX_PHENOTYPE_VERSIONS):
        tmp = os.path.join(path, f"{i:03d}")

        if not os.path.isdir(tmp):
            path = tmp
            os.mkdir(path, mode = _constants.DEFAULT_DIRMODE)
            break

    if not os.path.isdir(path):
        raise IsADirectoryError("Couldn't specify version of this directory")

    logger.info("Made direcotory %s", path)


    # update the configuration parameters to be consistent with the directory
    # structure created
    cfg.set("common", "path", path)
    cfg.set("common", "schema", schema_name)
    cfg.set("common", "phenotype", phenotype)
    cfg.set("common", "covariates_file", _constants.FILENAME_COVARIATES)
    cfg.set("common", "phenotype_file", _constants.FILENAME_PHENOTYPE)
    cfg.set("query", "db_pw_env", _constants.ENV_DB_PW)


    # # if another configuration file is provided, read and and sections,
    # # options and values to the initialized configuration.
    # if config is not None and os.path.isfile(config):

    #     user_cfg = _config.ConfigParser()
    #     user_cfg.read(config)
    #     _config.merge_cfg(cfg, user_cfg)

    # elif config is not None and not os.path.isfile(config):
    #     raise FileNotFoundError(f"The input file {config} is not found.")


    config_filename = os.path.join(path, _constants.FILENAME_CONFIG)
    with open(config_filename, "w") as fid:
        cfg.write(fid)

    logger.info("Wrote config file to %s", config_filename)
    
