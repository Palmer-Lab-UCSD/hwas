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

import re
from dataclasses import dataclass
from collections.abc import Iterable
import datetime
import psycopg as pg

from . import _set_parameters
from . import _constants


IS_COVARIATE = re.compile("covariate")
PHENOTYPE_TABLENAME = "gwas_phenotypes"
METADATA_TABLENAME = "descriptions"

SAMPLE_COLNAME = "rfid"

COVARIATE_DELIMITER = ','
COVARIATE_TYPE_TOKEN = '%covariate%'

MEASURE_TYPE_COLNAME = 'trait_covariate'
MEASURE_NAME_COLNAME = "measure"

OUTPUT_DELIMITER = ','


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


class SetQueryParameters(_set_parameters.BaseSetParameter):
    """Determine and set the query parameters

    Query parameters may be specified as command line inputs or
    read from a configuration file.  The values that specified 
    at the command line take precedence over configuration file.
    """

    def __init__(self,
                 dbname: str | None,
                 host: str | None,
                 port: str | None,
                 user: str | None,
                 schema: str | None,
                 phenotype: str | None,
                 config_filename: str | None) -> None:


        super().__init__(config_filename, "query") 

        self.dbname = self.set_parameter("dbname", dbname)
        self.host = self.set_parameter("host", host)
        self.port = self.set_parameter("port", port)
        self.user = self.set_parameter("user", user)
        self.schema = self.set_parameter("schema", schema)
        self.phenotype = self.set_parameter("phenotype", phenotype)

        self.outdir = self.set_parameter("outdir", None)

        password_env_var = self.set_parameter("db_password_env_var", None)

        self.password = None
        if password_env_var in os.environ:
            self.password = os.environ[password_env_var]
        


    def set_parameter(self,
                        name: str,
                        val: str | None) -> str:

        if val is not None:
            return val

        if name in self._config:
            val = self._config[name]
            
        if val is None:
            raise ValueError(f"The parameter {name} is set to None.")

        return val


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


def make_output_metadata(dbname: str,
                         schema: str,
                         phenotype: str,
                         cmd: str) -> str:
    date = datetime.datetime.now(datetime.UTC)
    
    return (f"{_constants.DEFAULT_META_PREFIX}date"
              f"={date.year}-{date.month:02}-{date.day:02}"
              f"-{date.hour:02}:{date.minute:02}:{date.second:02}\n"
        f"{_constants.DEFAULT_META_PREFIX}timezone={date.tzinfo}\n"
        f"{_constants.DEFAULT_META_PREFIX}user={os.environ['USER']}\n"
        f"{_constants.DEFAULT_META_PREFIX}database={dbname}\n"
        f"{_constants.DEFAULT_META_PREFIX}schema={schema}\n"
        f"{_constants.DEFAULT_META_PREFIX}phenotype={phenotype}\n"
        f"{_constants.DEFAULT_META_PREFIX}pipeline_version={_constants.VERSION}\n"
        f"{_constants.DEFAULT_META_PREFIX}input_command={cmd}\n")


def get_covariate_names(cur_meta: pg.Cursor,
                        schema: str,
                        phenotype: str) -> list[str]:

        query = pg.sql.SQL(
                        "SELECT * FROM {schema_name}.{table_name}" 
                        " WHERE measure = %s;"
                    ).format(
                        schema_name = pg.sql.Identifier(schema),
                        table_name = pg.sql.Identifier(METADATA_TABLENAME)
                    )

        if ((metadata := cur_mdata.execute(query, (phenotype,))) is None 
                or metadata.rowcount != 1):
            raise ValueError((f"Phenotype, {args.phenotype}, is not uniquely defined"
                    f" in {args.schema}.{_db.METADATA_TABLENAME}."))

        metadata = metadata.fetchone()

        # make sure phenotype is not a covariate
        if IS_COVARIATE.match(metadata.trait_covariate) is not None:
            raise ValueError("Input phenotype is a covariate.")


        # Check that specified covariate is a covariate and that a
        # column for that covariate exists in the PHENOTYPE_TABLENAME
        # table.
        covariate_names = metadata.covariates.split(_db.COVARIATE_DELIMITER)

        for w in covariate_names:

            if not _db.is_covariate(cur, args.schema, w):
                raise ValueError(f"The covariate {w} specified in the database"
                                 f" for phenotype {args.phenotype}"
                                 " is not labeled a covariate in the table"
                                 f" {args.schema}.{_db.METADATA_TABLENAME}.")


        return covariate_names
