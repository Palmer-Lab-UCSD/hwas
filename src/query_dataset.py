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
import argparse
import re
import os

import config
import db

import psycopg as pg



def parse_args(args: list[str]) -> argparse.Namespace:

    parser = argparse.ArgumentParser()
    parser.add_argument("--dbname",
            type=str,
            default=config.DEFAULT_DB_NAME,
            help="Name of data base to query from.")
    parser.add_argument("--host",
            type=str,
            default=config.DEFAULT_DB,
            help="Hostname of database")
    parser.add_argument("--user",
            type=str,
            default=config.DEFAULT_DB_USER,
            help="User name for logging into database.")
    parser.add_argument("--port",
            type=str,
            default=config.DEFAULT_DB_PORT,
            help="Port in which to connect to db.")
    parser.add_argument("-p",
            type=str,
            default=None,
            help="Password for logging into database")
    parser.add_argument("-o",
            type=str,
            default=None,
            help="Directory to print queried data output.")

    parser.add_argument("schema_name",
            type=str,
            help="Name of schema to query from.")
    parser.add_argument("phenotype",
            type=str,
            help="Name of phenotype to query.")
    
    return parser.parse_args(args)


def main(args: list[str]) -> None:
    args = parse_args(args)

    connection_name = (f"dbname={args.dbname}"
            f" host={args.host}"
            f" port={args.port}"
            f" user={args.user}")

    if args.p is not None:
        connection_name = f"{connection_name} password={args.p}"


    with (pg.connect(connection_name) as conn,
        conn.cursor() as cur,
        conn.cursor(row_factory=pg.rows.class_row(db.Metadata)) as cur_mdata,
        conn.cursor(row_factory=pg.rows.class_row(db.PhenotypeRecord)) as cur_data):



        # determine wither schema exists
        if not db.is_schema_unique(cur, args.schema_name):
            raise ValueError(f"The schema, {args.schema_name}, is not"
                " uniquely defined in db")
    

        # determine whether schema has required tables
        if not db.is_table_unique(cur, args.schema_name, db.METADATA_TABLENAME):
            raise ValueError(f"The table, {db.METADATA_TABLENAME} is not"
                    f" uniquely defined in schema {args.schema_name}.")


        if not db.is_table_unique(cur, args.schema_name, db.PHENOTYPE_TABLENAME):
            raise ValueError(f"The table, {db.PHENOTYPE_TABLENAME} is not"
                    f" uniquely defined in schema {args.schema_name}.")
            
        

        # use meta data to determine what covariates to use for the 
        # association analysis
        query = pg.sql.SQL(
                        "SELECT * FROM {schema_name}.{table_name}" 
                        " WHERE measure = %s;"
                    ).format(
                        schema_name = pg.sql.Identifier(args.schema_name),
                        table_name = pg.sql.Identifier(db.METADATA_TABLENAME)
                    )

        if ((metadata := cur_mdata.execute(query, (args.phenotype,))) is None 
                or metadata.rowcount != 1):
            raise ValueError((f"Phenotype, {args.phenotype}, is not uniquely defined"
                    f" in {args.schema_name}.{db.METADATA_TABLENAME}."))

        metadata = metadata.fetchone()

        # make sure phenotype is not a covariate
        if config.IS_COVARIATE.match(metadata.trait_covariate) is not None:
            raise ValueError("Input phenotype is a covariate.")


        # Check that specified covariate is a covariate and that a
        # column for that covariate exists in the PHENOTYPE_TABLENAME
        # table.
        covariate_names = metadata.covariates.split(db.COVARIATE_DELIMITER)

        for w in covariate_names:

            if not db.is_covariate(cur, args.schema_name, w):
                raise ValueError(f"The covariate {w} specified in the database"
                                 f" for phenotype {args.phenotype}"
                                 " is not labeled a covariate in the table"
                                 f" {args.schema_name}.{db.METADATA_TABLENAME}.")


        cov_columns = [db.SAMPLE_COLNAME] + covariate_names
        query = pg.sql.SQL(
                    "SELECT {fields} FROM {schema_name}.{table_name}"
                ).format(
                    fields = pg.sql.SQL(',').join(
                        [
                            pg.sql.Identifier(w) for w in cov_columns
                        ]
                    ),
                    schema_name = pg.sql.Identifier(args.schema_name),
                    table_name = pg.sql.Identifier(db.PHENOTYPE_TABLENAME)
                )

        if ((cov_out := cur.execute(query)) is None
            or cov_out.rowcount == 0):
            raise ValueError("No covariate records found")

        with (open(os.path.join(args.o, f"{args.phenotype}_covariates.csv"), "w")
            as fid):

            # write header
            fid.write(f"{db.COVARIATE_DELIMITER.join(covariate_names)}\n")

            #write query results
            for w in cov_out:
                fid.write(f"{db.COVARIATE_DELIMITER.join(w)}\n")

        


if __name__ == "__main__":
    main(sys.argv[1:])
