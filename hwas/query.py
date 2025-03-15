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
import datetime
import psycopg as pg


from . import _db
from . import _constants
from . import _set_parameters



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




def run(dbname: str | None,
        host: str | None,
        port: str | int | None,
        user: str | None,
        schema: str,
        phenotype: str,
        ) -> None:

    args = SetQueryParameters(dbname, host, port, user, schema, phenotype,
                           _constants.DEFAULT_CONFIG_FILENAME)

    connection_name = (f"dbname={args.dbname}"
            f" host={args.host}"
            f" port={args.port}"
            f" user={args.user}")

    if args.password is not None:
        connection_name = f"{connection_name} password={args.password}"


    with (pg.connect(connection_name) as conn,
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
        query = pg.sql.SQL(
                        "SELECT * FROM {schema_name}.{table_name}" 
                        " WHERE measure = %s;"
                    ).format(
                        schema_name = pg.sql.Identifier(args.schema),
                        table_name = pg.sql.Identifier(_db.METADATA_TABLENAME)
                    )

        if ((metadata := cur_mdata.execute(query, (args.phenotype,))) is None 
                or metadata.rowcount != 1):
            raise ValueError((f"Phenotype, {args.phenotype}, is not uniquely defined"
                    f" in {args.schema}.{_db.METADATA_TABLENAME}."))

        metadata = metadata.fetchone()

        # make sure phenotype is not a covariate
        if _db.IS_COVARIATE.match(metadata.trait_covariate) is not None:
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


        cov_columns = [_db.SAMPLE_COLNAME] + covariate_names
        query = pg.sql.SQL(
                    "SELECT {fields} FROM {schema_name}.{table_name}"
                ).format(
                    fields = pg.sql.SQL(',').join(
                        [
                            pg.sql.Identifier(w) for w in cov_columns
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

            date = datetime.datetime.now(datetime.UTC)


            fid.write(f"{_constants.DEFAULT_META_PREFIX}date"
                      f"={date.year}-{date.month:02}-{date.day:02}"
                      f"-{date.hour:02}:{date.minute:02}:{date.second:02}\n")
            fid.write(f"{_constants.DEFAULT_META_PREFIX}timezone={date.tzinfo}\n")
            fid.write(f"{_constants.DEFAULT_META_PREFIX}user={os.environ['USER']}\n")
            fid.write(f"{_constants.DEFAULT_META_PREFIX}database={args.dbname}\n")
            fid.write(f"{_constants.DEFAULT_META_PREFIX}schema={args.schema}\n")
            fid.write(f"{_constants.DEFAULT_META_PREFIX}phenotype={args.phenotype}\n")
            fid.write(f"{_constants.DEFAULT_META_PREFIX}pipeline_version={_constants.VERSION}\n")




            # write header
            fid.write(f"{_db.COVARIATE_DELIMITER.join(covariate_names)}\n")

            #write query results
            for w in cov_out:
                fid.write(f"{_db.COVARIATE_DELIMITER.join(w)}\n")

