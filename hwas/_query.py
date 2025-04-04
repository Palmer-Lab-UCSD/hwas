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


def interface(**kwargs) -> None:

    pars = _config.get_config_section(_constants.FILENAME_CONFIG,
                                    "query")

    pars.update(**kwargs)

    if not pars.is_specification_complete():
        print(pars)
        raise ValueError("A required value wasn't specified, see above list.")

    _config.update_config_section(pars)

    pw = os.environ.get(pars.db_pw_env)         # type: ignore[attr-defined]

    with (_db.connect(pars.dbname,              # type: ignore[attr-defined]
                      pars.host,                # type: ignore[attr-defined]
                      pars.port,                # type: ignore[attr-defined]
                      pars.db_user,             # type: ignore[attr-defined]
                      pw) as conn,              
          _db.data_cur(conn) as data_cur):



        # determine wither schema exists
        if not _db.is_schema_unique(data_cur, pars.schema):         # type: ignore[attr-defined]
            raise ValueError(f"The schema, {pars.schema}, is not"   # type: ignore[attr-defined]
                " uniquely defined in db")
    

        # determine whether schema has required tables
        if not _db.is_table_unique(data_cur,
                                   pars.schema,                     # type: ignore[attr-defined]
                                   _constants.METADATA_TABLENAME):  
            raise ValueError(f"The table, {_constants.METADATA_TABLENAME} is not"
                    f" uniquely defined in schema {pars.schema}.")  # type: ignore[attr-defined]



        if not _db.is_table_unique(data_cur,
                                   pars.schema,                     # type: ignore[attr-defined]
                                   _constants.PHENOTYPE_TABLENAME):
            raise ValueError(f"The table, {_constants.PHENOTYPE_TABLENAME} is not"
                    f" uniquely defined in schema {pars.schema}.")  # type: ignore[attr-defined]
        

        # covariate data
        covariate_names = _db.get_covariate_names(data_cur,
                                                  pars.schema,      # type: ignore[attr-defined]
                                                  pars.phenotype)   # type: ignore[attr-defined]

        logger.info("Identified covariates: %s",covariate_names)
        colnames, cov_out = _db.get_records(data_cur,
                                            pars.schema,            # type: ignore[attr-defined]
                                            _constants.PHENOTYPE_TABLENAME,  
                                            covariate_names,
                                            _constants.SAMPLE_COLNAME)

        logger.info("Covariate table columns: %s", colnames)

        meta_data = _io.make_output_metadata(pars.dbname,           # type: ignore[attr-defined]
                                             pars.schema,           # type: ignore[attr-defined]
                                             pars.phenotype,        # type: ignore[attr-defined]
                                             kwargs["cmd"])

        _io.write_to_file(pars.covariates_file,                     # type: ignore[attr-defined]
                          cov_out,
                          colnames,
                          meta_data = meta_data)


        # phenotype data
        colnames, pheno_out = _db.get_records(data_cur,
                                              pars.schema,          # type: ignore[attr-defined]
                                              _constants.PHENOTYPE_TABLENAME,
                                              [pars.phenotype],     # type: ignore[attr-defined]
                                              _constants.SAMPLE_COLNAME)

        logger.info("Phenotype table columns: %s", colnames)

        _io.write_to_file(pars.phenotype_file,                      # type: ignore[attr-defined]
                          pheno_out,
                          colnames,
                          meta_data = meta_data)
