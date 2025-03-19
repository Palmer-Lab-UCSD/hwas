
import os
import configparser
from . import _constants


class BaseSetParameter:
    def __init__(self,
                 config_filename: str | None,
                 name: str) -> None:

        if config_filename is None:
            config_filename = _constants.DEFAULT_CONFIG_FILENAME

        if not os.path.exists(config_filename):
            raise FileNotFoundError("Required configuration file isn't found."
                                    " Configuration files can be generated by"
                                    " the `init` subcommand.")

        self._config = configparser.ConfigParser(
                            interpolation = configparser.ExtendedInterpolation()
                        )
        self._config.read(config_filename)

        if name in self._config:
            self._config = self._config[name]
        else:
            raise ValueError("The query section in the configuration"
                             " file does not exist.")


class SetQueryParameters(BaseSetParameter):
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
