"""Query data from the Palmer Lab database

Query covariates and phenotypic measurements from the Palmer Lab 
database.

"""

# Query palmer lab database 
#
import sys
import re
import os
import logging
import psycopg as pg


from . import _constants
from . import _config
from . import _db
from . import _io


logger = logging.getLogger(__name__)


def run(dbname: str | None,
        host: str | None,
        port: str | int | None,
        db_user: str | None,
        cmd: str) -> None:

    if isinstance(port, int):
        port = str(port)

    pars = _config.get_config_section(_constants.FILENAME_CONFIG,
                                    "query")

    if dbname is not None:
        pars.dbname = dbname
    if host is not None:
        pars.host = host
    if port is not None:
        pars.port = port
    if db_user is not None:
        pars.db_user = db_user

    if pars.env_pw is None:
        pw = None
    else:
        pw = os.environ.get(pars.env_pw)

    with (_db.connect(pars.dbname,
                      pars.host,
                      pars.port,
                      pars.db_user,
                      pw) as conn,
          _db.data_cur(conn) as data_cur):



        # determine wither schema exists
        if not _db.is_schema_unique(data_cur, pars.schema):
            raise ValueError(f"The schema, {pars.schema}, is not"
                " uniquely defined in db")
    

        # determine whether schema has required tables
        if not _db.is_table_unique(data_cur, pars.schema, _constants.METADATA_TABLENAME):
            raise ValueError(f"The table, {_constants.METADATA_TABLENAME} is not"
                    f" uniquely defined in schema {pars.schema}.")


        if not _db.is_table_unique(data_cur, pars.schema, _constants.PHENOTYPE_TABLENAME):
            raise ValueError(f"The table, {_constants.PHENOTYPE_TABLENAME} is not"
                    f" uniquely defined in schema {pars.schema}.")
        

        # covariate data
        covariate_names = _db.get_covariate_names(data_cur, pars.schema, pars.phenotype)

        logger.info("Identified covariates: %s",covariate_names)
        colnames, cov_out = _db.get_records(data_cur,
                                            pars.schema,
                                            _constants.PHENOTYPE_TABLENAME,  
                                            covariate_names,
                                            _constants.SAMPLE_COLNAME)

        logger.info("Covariate table columns: %s", colnames)

        meta_data = _io.make_output_metadata(pars.dbname,
                                               pars.schema,
                                               pars.phenotype,
                                               cmd)

        _io.write_to_file(pars.covariates_file,
                          cov_out,
                          colnames,
                          meta_data = meta_data)


        # phenotype data
        colnames, pheno_out = _db.get_records(data_cur,
                                              pars.schema,
                                              _constants.PHENOTYPE_TABLENAME,
                                              [pars.phenotype],
                                              _constants.SAMPLE_COLNAME)

        logger.info("Phenotype table columns: %s", colnames)

        _io.write_to_file(pars.phenotype_file,
                          pheno_out,
                          colnames,
                          meta_data = meta_data)
