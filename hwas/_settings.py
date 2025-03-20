
import os
import configparser

from . import _constants


class BaseSetParameter:
    def __init__(self,
                 config_filename: str | None,
                 name: str) -> None:

        self._config_filename = config_filename

        if config_filename is None:
            self._config_filename = _constants.DEFAULT_CONFIG_FILENAME

        if not os.path.exists(self._config_filename):
            raise FileNotFoundError("Required configuration file isn't found."
                                    " Configuration files can be generated by"
                                    " the `init` subcommand.")

        self._config = configparser.ConfigParser(
                            interpolation = configparser.ExtendedInterpolation()
                        )
        self._config.read(self._config_filename)

        if name in self._config:
            self._config_section = self._config[name]
        else:
            raise ValueError("The query section in the configuration"
                             " file does not exist.")

    def set_parameter(self,
                      name: str,
                      val: str | None) -> str:

        if val is not None:
            self._config_section[name] = val
            return val

        if name in self._config_section:
            val = self._config_section[name]
            
        if val is None:
            raise ValueError(f"The parameter {name} is set to None.")

        return val


    def update(self):
        j = 1
        while ((archive_config := f".{self._config_filename}_{j}") is not None
               and os.path.exists(archive_config)):
            j+=1

        os.rename(self._config_filename, archive_config)

        with open(self._config_filename, "w") as fid:
            self._config.write(fid)


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

        self.covariates_file = self.set_parameter(
                        "covariates_file",
                        os.path.join(self._config["common"]["path"],
                                    (f"{self.phenotype}"
                                    f"{_constants.COVARIATE_FILE_SUFFIX}"))
                        )
        self.phenotype_file = self.set_parameter(
                        "phenotype_file",
                        os.path.join(self._config["common"]["path"],
                                           (f"{self.phenotype}"
                                           f"{_constants.PHENOTYPE_FILE_SUFFIX}"))
                        )

        password_env_var = self.set_parameter("env_pw", None)

        self.password = None
        if password_env_var in os.environ:
            self.password = os.environ[password_env_var]

        self.update()




class SetIntersectParameters(BaseSetParameter):
    def __init__(self,
                 vcf: str,
                 config_filename: str | None = None) -> None:

        super().__init__(config_filename, "intersect")

        if not os.path.exists(vcf):
            raise FileNotFoundError(f"VCF {vcf} coult not be found.")


        self.covariates_file = self._config_section["covariates_file"]
        self.phenotype_file = self._config_section["phenotype_file"]
        self.bcftools = self._config_section["bcftools"]
        self.rscript = self._config_section["rscript"]

        self.vcf = self.set_parameter("vcf", vcf)
        self.update()
