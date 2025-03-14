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

from dataclasses import dataclass
from collections.abc import Iterable
import psycopg as pg


PHENOTYPE_TABLENAME = "gwas_phenotypes"
METADATA_TABLENAME = "descriptions"

SAMPLE_COLNAME = "rfid"

COVARIATE_DELIMITER = ','
COVARIATE_TYPE_TOKEN = '%covariate%'

MEASURE_TYPE_COLNAME = 'trait_covariate'
MEASURE_NAME_COLNAME = "measure"


@dataclass
class Metadata:
    measure: str
    description: str
    trait_covariate: str
    covariates: str


@dataclass
class PhenotypeRecord:
    rfid:str
    measurement:str



def is_schema_unique(cur: pg.Cursor,
                     schema_name: str) -> bool:

    out = cur.execute(("SELECT 1 FROM information_schema.schemata"
                        " WHERE schema_name = %s"),
                        (schema_name,))

    if out.rowcount == 0:
        raise ValueError(f"The schema, {schema_name}, is not"
                        " defined in db")

    if out.rowcount > 1:
        return False

    return True


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


    if out.rowcount == 0:
        raise ValueError(f"The table {table_name} is not defined"
                         f" in schema {schema_name}.")

    if out.rowcount > 1:
        return False

    return True


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
                    data_name_colname = pg.sql.Identifier(MEASURE_NAME_COLNAME),
                    schema_name = pg.sql.Identifier(schema_name),
                    table_name = pg.sql.Identifier(METADATA_TABLENAME),
                    data_type_colname = pg.sql.Identifier(MEASURE_TYPE_COLNAME)
                ),
                (COVARIATE_TYPE_TOKEN,)
            )


    if out.rowcount == 0:
        raise ValueError(f"The measurement {measurement} is not defined"
                         f" in {schema_name}.")

    if out.rowcount > 1:
        raise ValueError(f"The phenotype {measurement} is not uniquely defined"
                         f" in {schema_name}.")

    return out.fetchone()[0]


