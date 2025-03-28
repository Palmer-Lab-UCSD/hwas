"""Initialize parameter settings and directory structure.



Typical usage example:
    init.run(config, account_numer, quality_of_service_code,
                '/home/me/.local/bin',
                database_name, environment_variable_for_password,
                schema_name, phenotype_measurement_name)

Functions:
    run: set up directories and configuration file
"""

import os
import logging

from . import _constants
from . import _config
from . import _templates

logger = logging.getLogger(__name__)


def run(schema_name: str,
        phenotype: str,
        partition: str | None,
        account: str | None,
        qos: str | None,
        config: str | None = None, 
        bin: str | None = None) -> None:

    
    # Use template to initialize configuration file.  A subset of
    # these values will be initialized to those sepecied by environmental
    # variables
    cfg = _config.init()


    # set specified parameter values to corresponding configuration
    # option
    if account is not None:
        cfg.set("slurm", "account", account)

    if qos is not None:
        cfg.set("slurm", "qos", qos)

    if bin is not None:
        cfg.set("common", "bin", os.path.expanduser(bin))


    # set up path and directories for phenotype and
    # project (schema_name)

    schema_path = os.path.join(os.path.expanduser(os.getcwd()),
                               schema_name)

    if not os.path.isdir(schema_path):
        os.mkdir(schema_name, mode = _constants.DEFAULT_DIRMODE)

    path = os.path.join(schema_path, phenotype)

    if os.path.isdir(path):
        raise FileExistsError(f"Analysis for {schema_name}/{phenotype}"
                         " already exists in the current directory."
                         " Either delete these data or change directories.")

    os.mkdir(path, mode = _constants.DEFAULT_DIRMODE)
    logger.info("Made direcotory %s", path)


    # update the configuration parameters to be consistent with the directory
    # structure created
    cfg.set("common", "path", path)
    cfg.set("common", "schema", schema_name)
    cfg.set("common", "phenotype", phenotype)
    cfg.set("common", "phenotype_file",
            os.path.join(path, _constants.FILENAME_PHENOTYPE))
    cfg.set("common", "covariates_file",
            os.path.join(path, _constants.FILENAME_COVARIATES))

    
    os.mkdir(os.path.join(path, cfg.get("common","logs")),
             mode = _constants.DEFAULT_DIRMODE)
    logger.info("Made direcotory %s", cfg.get("common", "logs"))



    # if another configuration file is provided, read and and sections,
    # options and values to the initialized configuration.
    if config is not None and os.path.isfile(config):

        user_cfg = _config.ConfigParser()
        user_cfg.read(config)
        _config.merge_cfg(cfg, user_cfg)

    elif config is not None and not os.path.isfile(config):
        raise FileNotFoundError(f"The input file {config} is not found.")


    config_filename = os.path.join(path, _constants.FILENAME_CONFIG)
    with open(config_filename, "w") as fid:
        cfg.write(fid)

    logger.info("Wrote config file to %s", config_filename)


    try: 

        tname = _templates.get_template_filename(_constants.TEMPLATES_HWAS_SBATCH)

    except FileNotFoundError as err:

        raise FileNotFoundError(f"Template {_constants.TEMPLATES_HWAS_SBATCH}"
                                " not found.  Please submit a GitHub issue"
                                " at https://github.com/Palmer-Lab-UCSD/hwas") from None


    # write out script files
    hwas_tmp = _templates.render(tname, cfg.section_to_dict("query"))

    with open(os.path.join(cfg.get("common", "path"), "hwas.sh"), "w") as fid:
        fid.write(hwas_tmp)
        

    hgrm_tmp = _templates.render(tname, cfg.section_to_dict("hgrm"))

    with open(os.path.join(cfg.get("common", "path"), "hwas.sh"), "w") as fid:
        fid.write(tmp)
