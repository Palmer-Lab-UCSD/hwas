

import os
import configparser

from . import __version__
from . import _constants



def run(config: str | None, 
        account: str | None,
        qos: str | None,
        dbname: str | None,
        db_pw_env_var: str | None,
        schema_name: str,
        phenotype: str) -> None:

    cfg = configparser.ConfigParser()

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


    cfg["common"] = dict(
            version = __version__.version,
            path = path,
            schema = schema_name,
            phenotype = phenotype,
            user = os.environ["USER"],
            logs = log_path
            )

    cfg["slurm"] = dict()
    if account is not None:
        cfg["slurm"]["account"] = account

    if qos is not None:
        cfg["slurm"]["qos"] = qos

    cfg["query"] = dict(
            host = _constants.DEFAULT_DB_HOST,
            port = _constants.DEFAULT_DB_PORT,
            user = _constants.DEFAULT_DB_USER,
            dbname = _constants.DEFAULT_DB_NAME,
            db_password_env_var = _constants.DEFAULT_DB_PASSWORD_ENV_VAR,
            outdir = '${common:path}',
            schema = '${common:schema}',
            phenotype = '${common:phenotype}',
            )

    if "PALMER_DB_USERNAME" in os.environ:
        cfg["query"]["username"] = os.environ["PALMER_DB_USERNAME"]

    if "PALMER_DB_HOST" in os.environ:
        cfg["query"]["host"] = os.environ["PALMER_DB_HOST"]

    if "PALMER_DB_PORT" in os.environ:
        cfg["query"]["port"] = os.environ["PALMER_DB_PORT"]

    if "PALMER_DB_NAME" in os.environ:
        cfg["query"]["dbname"] = os.environ["PALMER_DB_NAME"]

    if "PALMER_BIN" in os.environ:
        cfg["common"]["bin"] = os.environ["PALMER_BIN"]


    cfg["output"] = dict(
            meta_prefix = _constants.DEFAULT_META_PREFIX
            )
        
    if dbname is not None:
        cfg["query"]["dbname"] = dbname 


    with (open(os.path.join(path, _constants.DEFAULT_CONFIG_FILENAME), "w")
          as fid):
        cfg.write(fid)

