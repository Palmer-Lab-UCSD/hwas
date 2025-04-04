"""Tools for managing configuration file

Functions:
    init: load template configuration file and apply environment variables
    merge_cfg: add or replace sections, options and values from one 
                configuration file to another.
    interface: entry point for command line manipulation and reporting of
        config file
"""
import warnings
import re
import os
import typing
import shutil
import configparser

from . import _constants
from . import _templates



class ConfigParser:
    """Load a query data from a configuration file.

    The Python standard library contains a module named `configparser` for
    which the ConfigParser class is defined.  Instantiating a class instance
    takes serveral arguments, e.g. interpolation.  The `hwas` pipepline
    assumes a specific configuration file and options to be set.  This 
    class narrows the use case of ConfigParser for `hwas` and modifies
    the `get` and `set` methods of `configparser.ConfigParser` class.


    Args:
        max_recursion_depth: (optional) the maximum number of elements that
            make a chain of interpolated options.
    """

    def __init__(self, max_recursion_depth: int = 10) -> None:

        if not isinstance(max_recursion_depth,int) or max_recursion_depth < 0:
            raise ValueError("The input maximum_recursion_depth="
                             f" {max_recursion_depth} is not a positive,"
                             " integer as required.")

        self.max_recursion_depth = max_recursion_depth

        self._cfg = configparser.ConfigParser(
                            interpolation = configparser.ExtendedInterpolation(),
                            allow_no_value = True)

    def read(self, filename: str) -> None:
        """Read in data from a configuation file."""
        self._cfg.read(filename)

    def write(self, fileobject: typing.TextIO) -> None:
        self._cfg.write(fileobject)

    def has_option(self, section: str, option: str) -> bool:
        return self._cfg.has_option(section, option)

    def has_section(self, section: str) -> bool:
        return self._cfg.has_section(section)

    def add_section(self, section: str) -> None:
        self._cfg.add_section(section)

    def options(self, section: str) -> list[str]:
        return self._cfg.options(section)

    def sections(self) -> list[str]:
        return self._cfg.sections()

    def get_option_interpolator(self,
                                section: str,
                                option: str) -> tuple[str,str] | None:
        """Reteive the section and option that an interpolator points to.

        Args:
            section: the configuration file section of the interpolator
            option: the configuration file option of the interplator

        Returns:
            tuple: (section, option) that interpolator points to
            None: if value associated with input (section, option) is
                None or not an extended interpolator.

        Raises:
            RuntimeError: incorrect decomposition of extended interpolator

        """
        _section = section
        _option = option
        
        if ((option_val := self._cfg.get(section, option, raw = True)) is None):
            return None

        if ((out := re.match("^\\${([-_\\w]+):?([-_\\w]*)}$", option_val))
            is None):

            return None

        captured_groups = out.groups()

        if len(captured_groups) != 2:
            raise RuntimeError(f"Unexpected decomposition of ({section},"
                               f" {option}), please submit an issue on GitHub"
                               f" at {_constants.GITHUB_URL}")

        if captured_groups.count('') == 1:
            return (_section, captured_groups[0])

        return captured_groups

    def is_interpolation(self,
                         section: str,
                         option: str) -> bool:
        # recall that extended interplations are $opt or ${section:option}
        s = self._cfg.get(section, option, raw = True)

        if s is None:
            return False

        return re.match("^\\${[-_:\\w]+}$", s) is not None

    def _traverse_interpolators(self, section, option) -> tuple[str, str]:

        depth = 0
        while self.is_interpolation(section, option):

            if depth > self.max_recursion_depth:
                raise RecursionError("Interpolation depth exceeds maximum of"
                                     " {self.max_recursion_depth}. To increase"
                                     " depth, set the max_recusion_depth to a"
                                     " larger number.")

            if (out := self.get_option_interpolator(section, option)) is None:
                break

            section, option = out
    
            depth += 1

        return (section, option)

    def section_to_dict(self, section: str) -> dict[str,str]:
        option_dict = {}

        for w in self.options(section):
            option_dict[w] = self.get(section, w)

        return option_dict

    def set(self,
            section: str,
            option: str,
            value: str | None) -> None:
        """Set the value of (section, option) in a loaded configuration."""

        section, option = self._traverse_interpolators(section, option)

        self._cfg.set(section, option, value)

    def get(self,
            section: str,
            option: str,
            raw: bool = False) -> str:
        """Retrieve value set to (section, option)."""

        if raw:
            return self._cfg.get(section, option, raw=raw)

        section, option = self._traverse_interpolators(section, option)

        try:
            return self._cfg.get(section, option)
        except TypeError:
            return self._cfg.get(section, option, raw=True)


class DynamicConfigSection:
    """Simplify option references in configuration file section

    Args:
        section: the section of the configuration file in which
            we want to store option keys and values.
    """
    def __init__(self, section: str) -> None:
        super().__setattr__("name", section)
        super().__setattr__("_dynamic_option_names", [])

    def __iter__(self):
        for w in self._dynamic_option_names:
            yield (w, getattr(self, w))

    def __contains__(self, option: str) -> bool:
        return option in self._dynamic_option_names

    def __getattr__(self, option: str) -> typing.Any:
        if option not in self.__dict__:
            raise AttributeError(f"{option} is not a valid attribute")

        if (val := self.__dict__[option]) == "None":
            return None

        return val

    def __str__(self) -> str:
        s = ""
        for option in self._dynamic_option_names:
            s = f"{s}\n{option} = {self.__getattr__(option)}"
        return s

    def __setattr__(self, option, val):
        if option not in self:
            self._dynamic_option_names.append(option)

        self.__dict__[option] = val

    def update(self, 
               **kwargs: typing.Mapping[str, str | None]) -> None:

        if kwargs is None:
            return 

        for option in self._dynamic_option_names:
            if option in kwargs and kwargs[option] is not None:
                self.__setattr__(option, kwargs[option])

    def is_specification_complete(self) -> bool:
        for option in self._dynamic_option_names:
            if self.__getattr__(option) is None:
                return False

        return True


def get_config_section(config_filename: str,
                       section_name: str) -> DynamicConfigSection:
    """Factory function of DynamicConfigSection class"""

    cfg = ConfigParser()
    cfg.read(config_filename)

    if not cfg.has_section(section_name):
        raise ValueError(f"{section_name} not defined in configuration.")

    pars = DynamicConfigSection(section_name)

    for opt in cfg.options(section_name):
        setattr(pars, opt, cfg.get(section_name, opt))

    return pars


def update_config_section(section_obj: DynamicConfigSection,
                          config_file: str = _constants.FILENAME_CONFIG) -> None:
    cfg = ConfigParser()
    cfg.read(config_file)

    for option, value in section_obj:
        cfg.set(section_obj.name, option, value)

    with open(config_file, "w") as fid:
        cfg.write(fid)


def _load_default_config() -> ConfigParser:

    cfg = ConfigParser()
    if ((config_template := _templates.get_template_filename(_constants.TEMPLATES_CONFIG))
        is None):

        raise FileNotFoundError(f"Template, {config_template}, not found.")

    cfg.read(config_template)

    return cfg


def init() -> ConfigParser:

    cfg = _load_default_config()

    # Set up directories for results
    cfg.set("common", "user",           os.environ.get("USER"))
    cfg.set("query", "dbname",          os.environ.get(_constants.ENV_DB_NAME))
    cfg.set("query", "host",            os.environ.get(_constants.ENV_DB_HOST))
    cfg.set("query", "port",            os.environ.get(_constants.ENV_DB_PORT))
    # cfg.set("query", "db_pw_env",       _constants.ENV_DB_PW)

    cfg.set("slurm", "qos",      os.environ.get(_constants.ENV_SLURM_QOS))
    cfg.set("slurm", "account",  os.environ.get(_constants.ENV_SLURM_ACCOUNT))
    cfg.set("slurm", "partition",os.environ.get(_constants.ENV_SLURM_PARTITION))

    cfg.set("common", "version", _constants.VERSION)
    cfg.set("common", "meta_prefix", _constants.DEFAULT_META_PREFIX)
    cfg.set("common", "header_prefix", _constants.DEFAULT_HEADER_PREFIX)

    
    return cfg
        


def merge_cfg(opt_receiver: ConfigParser,
              opt_donator: ConfigParser) -> None:

    for s in opt_donator.sections():

        if not opt_receiver.has_section(s):
            opt_receiver.add_section(s)

        for option in opt_donator.options(s):

            if (opt_receiver.has_option(s, option)
                and opt_receiver.is_interpolation(s, option)):

                warnings.warn(f"Changing interpolation {s}.{option}"
                              f" = {opt_receiver.get(s, option, raw=True)}"
                              f" value to {opt_donator.get(s, option)}.",
                              Warning)

            opt_receiver.set(s, option, opt_donator.get(s, option))


def interface(option: str | None = None,
              value: str | None = None,
              section: str = "common",
              subcommand: str | None = None) -> None:


    if not os.path.isfile(_constants.FILENAME_CONFIG):
        raise FileNotFoundError("No configuration file in current directory")

    cfg = ConfigParser()
    cfg.read(_constants.FILENAME_CONFIG)

    # set section.option = value
    if (section is not None 
        and option is not None
        and value is not None):

        cfg.set(section, option, value)

        shutil.move(_constants.FILENAME_CONFIG,
                    f".{_constants.FILENAME_CONFIG}")

        with open(_constants.FILENAME_CONFIG, "w") as fid:
            cfg.write(fid)


    # get the value of section.option
    if section is not None and option is not None:
        print(f"{section}.{option} = {cfg.get(section, option)}")
        return

    # enumerate options and sections
    if section is not None:

        print(f"[{section}]")
        max_word = 0

        for opt in cfg.options(section):
            if len(opt) > max_word:
                max_word = len(opt)

        col_num = max_word + 4
        for opt in cfg.options(section):
            end_word = 2 + len(opt)
            num_spaces = col_num - end_word
            s = f"  {opt}"
            for w in range(num_spaces):
                s = f"{s} "
            print(f"{s}{cfg.get(section, opt)}")

        return

    if section is None:
        print("Sections:")
        for sec in cfg.sections():
            print(f"  {sec}")

