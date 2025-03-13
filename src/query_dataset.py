# Query palmer lab database 
#
import sys
import argparse
import re
import os

import data_model as dm

import psycopg as pg






covariate_regex = re.compile("covariate")

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--dbname",
            type=str,
            default="PalmerLab_Datasets",
            help="Name of data base to query from.")
    parser.add_argument("--host",
            type=str,
            default="palmerlab-main-database-c2021-08-02.c6sgfwysomht.us-west-2.rds.amazonaws.com",
            help="Hostname of database")
    parser.add_argument("--user",
            type=str,
            default="postgres",
            help="User name for logging into database.")
    parser.add_argument("--port",
            type=str,
            default="5432",
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


def main(args):
    args = parse_args(args)

    connection_name = (f"dbname={args.dbname}"
            f" host={args.host}"
            f" port={args.port}"
            f" user={args.user}")

    if args.p is not None:
        connection_name = f"{connection_name} password={args.p}"


    with (pg.connect(connection_name) as conn,
        conn.cursor() as cur,
        conn.cursor(row_factory=pg.rows.class_row(dm.Metadata)) as cur_mdata,
        conn.cursor(row_factory=pg.rows.class_row(dm.PhenotypeRecord)) as cur_data):



        # determine wither schema exists
        if ((out := cur_mdata.execute("SELECT 1 FROM information_schema.schemata"
            " WHERE schema_name = %s", (args.schema_name,))) is None
            or out.rowcount != 1):
    
            raise ValueError(f"The schema, {args.schema_name}, is not"
                " uniquely defined in db")
    

        # determine whether schema has required tables
        if ((out := cur_mdata.execute(("SELECT 1 FROM"
            " (SELECT * FROM information_schema.tables WHERE table_schema = %s)"
            " WHERE table_name = %s OR table_name = %s"),
            (args.schema_name, dm.METADATA_TABLENAME, dm.PHENOTYPE_TABLENAME))) is None 
            or out.rowcount != 2):

            raise ValueError(f"The schema, {args.schema_name}, does not contain the"
                    f" required tables: '{dm.METADATA_TABLENAME}'"
                    f" and '{dm.PHENOTYPE_TABLENAME}'.")
        

        
        query = pg.sql.SQL(
                    "SELECT * FROM {schema_name}.{table_name} WHERE measure = %s"
                ).format(
                    schema_name = pg.sql.Identifier(args.schema_name),
                    table_name = pg.sql.Identifier(dm.METADATA_TABLENAME)
                )

        if ((out := cur_mdata.execute(query, (args.phenotype,))) is None 
                or out.rowcount != 1):
            raise ValueError((f"Phenotype, {args.phenotype}, is not uniquely defined"
                    f" in {args.schema_name}.description."))

        out = out.fetchone()


        # This code clip is adapted from the example at
        # https://www.psycopg.org/psycopg3/docs/api/sql.html#module-usage
        # It is meant to dynamically make a SQL query with columns names
        # specified by the covariates field from the previous query


        if covariate_regex.match(out.trait_covariate) is not None:
            raise ValueError("Input phenotype is a covariate.")


        # Check that specified covariate is a covariate and that a
        # column for that covariate exists in the PHENOTYPE_TABLENAME
        # table.
        covariate_names = ([dm.SAMPLE_COLNAME] 
                            + out.covariates.split(dm.COVARIATE_DELIMITER))

        for w in covariate_names:
            query = pg.sql.SQL(
                        "SELECT {column_name} IN"
                        " (SELECT column_name FROM information_schema.columns"
                        " WHERE table_schema = %s AND table_name = %s)"
                    ).format(
                        column_name = pg.sql.Literal(w)
                    )

            if ((tmp_out := cur.execute(query,
                            (args.schema_name, dm.PHENOTYPE_TABLENAME))) is None
                or tmp_out.rowcount != 1
                or not tmp_out.fetchone()[0]):

                raise ValueError(f"The covariate, {w}, is not a named column"
                                 f" in the {dm.PHENOTYPE_TABLENAME} table.")
                


        query = pg.sql.SQL(
                    "SELECT {fields} FROM {schema_name}.{table_name}"
                ).format(
                    fields = pg.sql.SQL(',').join(
                        [
                            pg.sql.Identifier(w) for w in covariate_names
                        ]
                    ),
                    schema_name = pg.sql.Identifier(args.schema_name),
                    table_name = pg.sql.Identifier(dm.PHENOTYPE_TABLENAME)
                )

        cov_out = cur.execute(query)

        with (open(os.path.join(args.o, f"{args.phenotype}_covariates.csv"), "w")
            as fid):

            # write header
            fid.write(f"{dm.COVARIATE_DELIMITER.join(covariate_names)}\n")

            #write query results
            for w in cov_out:
                fid.write(f"{dm.COVARIATE_DELIMITER.join(w)}\n")

        


if __name__ == "__main__":
    main(sys.argv[1:])
