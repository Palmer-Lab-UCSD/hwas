

import os
import configparser

from . import __version__



def run(config: str | None, 
        account: str | None,
        qos: str | None,
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
        os.mkdir(path, mode = 0o770)

    path = os.path.join(path, phenotype)

    if os.path.exists(path):
        raise FileExistsError(f"Analysis for {schema_name}/{phenotype}"
                         " already exists in the current directory."
                         " Either delete these data or change directories.")


    os.mkdir(path, mode = 0o750)

    log_path = os.path.join(path, "logs")
    os.mkdir(log_path, mode = 0o750)

    cfg["hwas_pkg"] = dict(
            version = __version__.version
            )

    cfg["common"] = dict(
            path = path,
            schema = schema_name,
            phenotype = phenotype,
            user = os.environ["USER"],
            logs = os.path.join(path, "logs")
            )

    cfg["slurm"] = dict()
    if account is not None:
        cfg["slurm"]["account"] = account

    if qos is not None:
        cfg["slurm"]["qos"] = qos

    cfg["db"] = dict(
            host = "localhost",
            port = 5432,
            username = os.environ["USER"],
            dbname = "data"
            )

    if "PALMER_DB_USERNAME" in os.environ:
        cfg["db"]["username"] = os.environ["PALMER_DB_USERNAME"]

    if "PALMER_DB_HOST" in os.environ:
        cfg["db"]["host"] = os.environ["PALMER_DB_HOST"]

    if "PALMER_DB_PORT" in os.environ:
        cfg["db"]["port"] = os.environ["PALMER_DB_PORT"]

    if "PALMER_DB_DBNAME" in os.environ:
        cfg["db"]["dbname"] = os.environ["PALMER_DB_DBNAME"]

    if "PALMER_BIN" in os.environ:
        cfg["common"]["bin"] = os.environ["PALMER_BIN"]


    cfg["output"] = dict(
            meta_prefix = "##"
            )
        
    with open(os.path.join(path, "config"), "w") as fid:
        cfg.write(fid)

