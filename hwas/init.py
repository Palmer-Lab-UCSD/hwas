

import os
import configparser
import string
import importlib.resources

from . import __version__
from . import _constants



def run(config: str | None, 
        account: str | None,
        qos: str | None,
        bin: str | None,
        dbname: str | None,
        env_pw: str | None,
        schema_name: str,
        phenotype: str) -> None:

    
    cfg = configparser.ConfigParser()
    if (len(cfg.read(importlib.resources
               .files(_constants.DEFAULT_CONFIG_PATH)
               .joinpath('config'))) != 1):
        raise FileNotFoundError("Default config file for package not found."
                                " Please submit an issue on GitHub.")


    if config is not None and os.path.exists(config):
        cfg.read(config)
    elif config is not None and not os.path.exists(config):
        raise FileNotFoundError(f"The input file {config} is not found.")


    # Set up directories for results
    path = os.path.join(os.getcwd(), schema_name)

    if not os.path.exists(path):
        os.mkdir(path, mode = _constants.DEFAULT_DIRMODE)

    path = os.path.join(path, phenotype)

    if os.path.exists(path):
        raise FileExistsError(f"Analysis for {schema_name}/{phenotype}"
                         " already exists in the current directory."
                         " Either delete these data or change directories.")


    os.mkdir(path, mode = _constants.DEFAULT_DIRMODE)

    log_path = os.path.join(path, _constants.DEFAULT_LOG_DIR)
    os.mkdir(log_path, mode = _constants.DEFAULT_DIRMODE)


    cfg["common"]["version"] = __version__.version
    cfg["common"]["path"] = path
    cfg["common"]["schema"] = schema_name
    cfg["common"]["phenotype"] = phenotype
    cfg["common"]["user"] = os.environ["USER"]
    cfg["common"]["logs"] = log_path
    cfg["common"]["bin"] = bin


    if account is not None:
        cfg["slurm"]["account"] = account

    if qos is not None:
        cfg["slurm"]["qos"] = qos

    cfg["query"]["user"] = _constants.DEFAULT_DB_USER
    cfg["query"]["env_pw"] = _constants.ENV_DB_PASSWORD

    if _constants.ENV_DB_USERNAME in os.environ:
        cfg["query"]["user"] = os.environ[_constants.ENV_DB_USERNAME]

    if _constants.ENV_DB_HOST in os.environ:
        cfg["query"]["host"] = os.environ[_constants.ENV_DB_HOST]

    if _constants.ENV_DB_PORT in os.environ:
        cfg["query"]["port"] = os.environ[_constants.ENV_DB_PORT]

    if _constants.ENV_DB_NAME in os.environ:
        cfg["query"]["dbname"] = os.environ[_constants.ENV_DB_NAME]

    if bin is None and _constants.ENV_BIN in os.environ:
        bin = os.environ[_constants.ENV_BIN]

    if bin is None:
        raise FileNotFoundError("Software bin directory not found.")

    cfg["common"]["bin"] = bin

    cfg["output"] = dict(
            meta_prefix = _constants.DEFAULT_META_PREFIX
            )
        
    if dbname is not None:
        cfg["query"]["dbname"] = dbname 

    


    with (open(os.path.join(path, _constants.DEFAULT_CONFIG_FILENAME), "w")
          as fid):
        cfg.write(fid)


    # TODO
    # with (open(os.path.join(path, ), "r") as fin,
    #       open(os.path.join(path, _constants.HWAS_BASH), "w") as fout):

    #     tmp = string.Template(fin.read())
        
