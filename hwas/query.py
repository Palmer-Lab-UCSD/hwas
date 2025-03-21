"""Query data from the Palmer Lab database

By: Robert Vogel
Affiliation: Palmer Lab at UCSD



DESCRIPTION

Query covariates and phenotypic measurements from the Palmer Lab 
database.



ACKNOWLEDGEMENT:

This code has been reviewed by Claude, the AI assistant from Anthropic.
The code was designed and implemented by Robert Vogel, code recommendations
that were provided by Claude were adapted and implemented by Robert Vogel.
"""

# Query palmer lab database 
#
import sys
import re
import os
import logging
import psycopg as pg


from . import _constants
from . import _settings
from . import _db
from . import _io


logger = logging.getLogger(__name__)


def run(dbname: str | None,
        host: str | None,
        port: str | int | None,
        user: str | None,
        schema: str,
        phenotype: str,
        cmd: str) -> None:

    if isinstance(port, int):
        port = str(port)

    args = _settings.SetQueryParameters(dbname, host, port, user, schema, phenotype,
                                        _constants.DEFAULT_CONFIG_FILENAME)

    with (_db.connect(args.dbname, args.host, args.port, args.user, 
                      args.password) as conn,
          _db.data_cur(conn) as data_cur):



        # determine wither schema exists
        if not _db.is_schema_unique(data_cur, args.schema):
            raise ValueError(f"The schema, {args.schema}, is not"
                " uniquely defined in db")
    

        # determine whether schema has required tables
        if not _db.is_table_unique(data_cur, args.schema, _constants.METADATA_TABLENAME):
            raise ValueError(f"The table, {_constants.METADATA_TABLENAME} is not"
                    f" uniquely defined in schema {args.schema}.")


        if not _db.is_table_unique(data_cur, args.schema, _constants.PHENOTYPE_TABLENAME):
            raise ValueError(f"The table, {_constants.PHENOTYPE_TABLENAME} is not"
                    f" uniquely defined in schema {args.schema}.")
        

        # covariate data
        covariate_names = _db.get_covariate_names(data_cur, args.schema, args.phenotype)

        logger.info("Identified covariates: %s",covariate_names)
        colnames, cov_out = _db.get_records(data_cur,
                                            args.schema,
                                            _constants.PHENOTYPE_TABLENAME,  
                                            covariate_names,
                                            _constants.SAMPLE_COLNAME)

        logger.info("Covariate table columns: %s", colnames)

        meta_data = _io.make_output_metadata(args.dbname,
                                               args.schema,
                                               args.phenotype,
                                               cmd)

        _io.write_to_file(args.covariates_file,
                          cov_out,
                          colnames,
                          meta_data = meta_data)


        # phenotype data
        colnames, pheno_out = _db.get_records(data_cur,
                                              args.schema,
                                              _constants.PHENOTYPE_TABLENAME,
                                              [args.phenotype],
                                              _constants.SAMPLE_COLNAME)

        logger.info("Phenotype table columns: %s", colnames)

        _io.write_to_file(args.phenotype_file,
                          pheno_out,
                          colnames,
                          meta_data = meta_data)
