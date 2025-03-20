"""Data types, tables, and functions for Palmer Lab database

By: Robert Vogel
Affiliation: Palmer Lab at UCSD

Acknowledgement:
    This code has been reviewed by Claude, the AI assistant from Anthropic.
    The code was designed and implemented by Robert Vogel, code recommendations
    that were provided by Claude were adapted and implemented by Robert Vogel.


The Palmer Lab data base is organized, roughly, as follows.
Experimental data produced by a single lab, award, or other
logical unit is organized in a schema.  Each schema contains,
at a minimum, two tables 'descriptions' and 'gwas_phenotypes'.
The 'descriptions' table enumerates all measurements in the 
'gwas_phenotypes' table, and documents which are to be used
as covariates vs. experimental phenotypes.  Moreover, for each
experimental phenotype the set of measured covariates to be
used are documented.

Example

Schema name:
    p50_investigator

Table names:
    descriptions,
    gwas_phenotypes,
    etc.
"""
import os
from collections.abc import Iterable
import datetime
import psycopg as pg

from . import _constants




def is_schema_unique(cur: pg.Cursor,
                     schema_name: str) -> bool:

    out = cur.execute(("SELECT 1 FROM information_schema.schemata"
                        " WHERE schema_name = %s"),
                        (schema_name,))

    if out.rowcount == 1:
        return True

    return False


def is_table_unique(cur: pg.Cursor,
                    schema_name: str,
                    table_name: str) -> bool:

    out = cur.execute(
                "SELECT 1 FROM"
                " (SELECT * FROM information_schema.tables" 
                " WHERE table_schema = %s)"
                " WHERE table_name = %s;",
                (schema_name, table_name)
            )

    if out.rowcount == 1:
        return True

    return False


def is_covariate(cur: pg.Cursor,
                 schema_name: str,
                 measurement: str) -> bool:

    out = cur.execute(
            pg.sql.SQL(
                "SELECT {measurement} IN"
                " (SELECT {data_name_colname} FROM {schema_name}.{table_name}"
                " WHERE {data_type_colname} LIKE %s);"
                ).format(
                    measurement = pg.sql.Literal(measurement),
                    data_name_colname = pg.sql.Identifier(_constants.MEASURE_NAME_COLNAME),
                    schema_name = pg.sql.Identifier(schema_name),
                    table_name = pg.sql.Identifier(_constants.METADATA_TABLENAME),
                    data_type_colname = pg.sql.Identifier(_constants.MEASURE_TYPE_COLNAME)
                ),
                (_constants.COVARIATE_TYPE_TOKEN,)
            )

    if out.rowcount == 1:
        return True

    return False


def get_covariate_names(cur: pg.Cursor,
                        schema: str,
                        phenotype: str) -> list[str]:

        query = pg.sql.SQL(
                        "SELECT * FROM {schema_name}.{table_name}" 
                        " WHERE measure = %s;"
                    ).format(
                        schema_name = pg.sql.Identifier(schema),
                        table_name = pg.sql.Identifier(_constants.METADATA_TABLENAME)
                    )

        if ((metadata := cur.execute(query, (phenotype,))) is None 
                or metadata.rowcount != 1):
            raise ValueError((f"Phenotype, {phenotype}, is not uniquely defined"
                    f" in {schema}.{_constants.METADATA_TABLENAME}."))

        metadata = metadata.fetchone()

        # make sure phenotype is not a covariate
        if _constants.IS_COVARIATE.match(metadata["trait_covariate"]) is not None:
            raise ValueError("Input phenotype is a covariate.")


        # Check that specified covariate is a covariate and that a
        # column for that covariate exists in the PHENOTYPE_TABLENAME
        # table.
        covariate_names = metadata["covariates"].split(_constants.COVARIATE_DELIMITER)

        for w in covariate_names:

            if not is_covariate(cur, schema, w):
                raise ValueError(f"The covariate {w} specified in the database"
                                 f" for phenotype {phenotype}"
                                 " is not labeled a covariate in the table"
                                 f" {schema}.{_constants.METADATA_TABLENAME}.")


        return covariate_names


def get_records(cur: pg.Cursor,
                schema_name: str,
                table_name: str,
                colnames: list[str],
                sample_colname: str | None) -> (Iterable[str], pg.Cursor):

    if not isinstance(sample_colname, str) and sample_colname is not None:
        raise ValueError("{sample_colname} must be either type str or None")

    column_names = colnames

    if sample_colname is not None:
        column_names = [sample_colname] + column_names


    query = pg.sql.SQL(
                "SELECT {fields} FROM {schema_name}.{table_name}"
            ).format(
                fields = pg.sql.SQL(',').join(
                    [
                        pg.sql.Identifier(w) for w in column_names
                    ]
                ),
                schema_name = pg.sql.Identifier(schema_name),
                table_name = pg.sql.Identifier(_constants.PHENOTYPE_TABLENAME)
            )
    
    if ((cov_out := cur.execute(query)) is None
        or cov_out.rowcount == 0):
        raise ValueError("No covariate records found")

    return (column_names, cov_out)



def connect(dbname: str,
            host: str,
            port: str,
            user: str,
            password: str | None) -> pg.Connection:
    """Open a connection with postgres database.

    """
    connection_name = (f"dbname={dbname}"
            f" host={host}"
            f" port={port}"
            f" user={user}")

    if password is not None:
        connection_name = f"{connection_name} password={password}"

    return pg.connect(connection_name) 


def data_cur(conn: pg.Connection) -> pg.Cursor:
    return conn.cursor(row_factory = pg.rows.dict_row)

