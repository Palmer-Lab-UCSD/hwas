# Query palmer lab database 
#
import sys
import argparse
import re
from dataclasses import dataclass

import psycopg as pg



@dataclass
class ColumnMetadata:
    measure: str
    description: str
    trait_covariate: str
    covariates: str


@dataclass
class PhenotypeRecord:
    rfid:str
    measurement:str


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
            required=True,
            help="Password for logging into database")
    parser.add_argument("-o",
            type=str,
            required=True,
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
            f" user={args.user}"
            f" password={args.p}"
            f" host={args.host}"
            f" port={args.port}")

    with (pg.connect(connection_name) as conn,
        conn.cursor(row_factory=pg.rows.class_row(ColumnMetadata)) as cur_mdata,
        conn.cursor(row_factory=pg.rows.class_row(PhenotypeRecord)) as cur_data):



        # determine wither schema exists
        if ((out := cur_mdata.execute("SELECT 1 FROM information_schema.schemata"
            " WHERE schema_name = %s", (args.schema_name,))) is None
            or out.rowcount != 1):
    
            raise ValueError(f"The schema, {args.schema_name}, is not"
                " uniquely defined in db")
    

        # determine whether schema has required tables
        if ((out := cur_mdata.execute(("SELECT 1 FROM"
            " (SELECT * FROM information_schema.tables WHERE table_schema = %s)"
            " WHERE table_name = 'descriptions' OR table_name = 'gwas_phenotypes'"),
            (args.schema_name,))) is None or out.rowcount != 2):

            raise ValueError(f"The schema, {args.schema_name}, does not contain the"
                    " required tables: 'descriptions' and 'gwas_phenotypes'")
        

        
        query = pg.sql.SQL(
                "SELECT * FROM {schema_name}.descriptions WHERE measure = %s"
                ).format(schema_name=pg.sql.Identifier(args.schema_name))

        if ((out := cur_mdata.execute(query, (args.phenotype,))) is None 
                or out.rowcount != 1):
            raise ValueError((f"Phenotype, {args.phenotype}, is not uniquely defined"
                    f" in {args.schema_name}.description."))

        out = out.fetchone()
        print(out)

        if covariate_regex.match(out.trait_covariate) is not None:
            raise ValueError("Input phenotype is a covariate.")


        


if __name__ == "__main__":
    main(sys.argv[1:])
