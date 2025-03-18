"""Query data from the Palmer Lab database

By: Robert Vogel
Affiliation: Palmer Lab at UCSD

Acknowledgement:
    This code has been reviewed by Claude, the AI assistant from Anthropic.
    The code was designed and implemented by Robert Vogel, code recommendations
    that were provided by Claude were adapted and implemented by Robert Vogel.
"""

# Query palmer lab database 
#
import sys
import re
import os
import psycopg as pg


from . import _db
from . import _constants



def run(dbname: str | None,
        host: str | None,
        port: str | int | None,
        user: str | None,
        schema: str,
        phenotype: str,
        cmd: str) -> None:

    args = _db.SetQueryParameters(dbname, host, port, user, schema, phenotype,
                           _constants.DEFAULT_CONFIG_FILENAME)

    with (_db.connect(args.dbname, args.host, args.port, args.user, 
                      args.password) as conn,
        conn.cursor() as cur,
        conn.cursor(row_factory=pg.rows.class_row(_db.Metadata)) as cur_mdata,
        conn.cursor(row_factory=pg.rows.class_row(_db.PhenotypeRecord)) as cur_data):



        # determine wither schema exists
        if not _db.is_schema_unique(cur, args.schema):
            raise ValueError(f"The schema, {args.schema}, is not"
                " uniquely defined in db")
    

        # determine whether schema has required tables
        if not _db.is_table_unique(cur, args.schema, _db.METADATA_TABLENAME):
            raise ValueError(f"The table, {_db.METADATA_TABLENAME} is not"
                    f" uniquely defined in schema {args.schema}.")


        if not _db.is_table_unique(cur, args.schema, _db.PHENOTYPE_TABLENAME):
            raise ValueError(f"The table, {_db.PHENOTYPE_TABLENAME} is not"
                    f" uniquely defined in schema {args.schema}.")
            
        

        # use meta data to determine what covariates to use for the 
        # association analysis
        covariate_names = _db.get_covariate_names(cur_mdata, args.schema, args.phenotype)
        column_names = [_db.SAMPLE_COLNAME] + covariate_names

        query = pg.sql.SQL(
                    "SELECT {fields} FROM {schema_name}.{table_name}"
                ).format(
                    fields = pg.sql.SQL(',').join(
                        [
                            pg.sql.Identifier(w) for w in column_names
                        ]
                    ),
                    schema_name = pg.sql.Identifier(args.schema),
                    table_name = pg.sql.Identifier(_db.PHENOTYPE_TABLENAME)
                )

        if ((cov_out := cur.execute(query)) is None
            or cov_out.rowcount == 0):
            raise ValueError("No covariate records found")


        with (open(os.path.join(args.outdir, f"{args.phenotype}_covariates.csv"), "w")
            as fid):
            
            fid.write(_db.make_output_metadata(args.dbname,
                                               args.schema,
                                               args.phenotype,
                                               cmd))

            # write header
            fid.write(f"{_db.OUTPUT_DELIMITER.join(column_names)}\n")

            #write query results
            for w in cov_out:
                fid.write(f"{_db.OUTPUT_DELIMITER.join(w)}\n")

