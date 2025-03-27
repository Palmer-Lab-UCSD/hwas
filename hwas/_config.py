"""Tools for managing configuration file

Functions:
    init: load template configuration file and apply environment variables
    merge_cfg: add or replace sections, options and values from one 
                configuration file to another.
"""
import warnings
import re
import os
import typing
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

    def write(self, fileobject: typing.TextIO) -> int:
        return self._cfg.write(fileobject)

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

        out = out.groups()

        if len(out) != 2:
            raise RuntimeError(f"Unexpected decomposition of ({section},"
                               f" {option}), please submit an issue on GitHub"
                               f" at {_constants.GITHUB_URL}")

        if out.count('') == 1:
            return (_section, out[0])

        return out

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

    def set(self,
            section: str,
            option: str,
            value: str | None) -> None:
        """Set the value of (section, option) in a loaded configuration."""

        section, option = self._traverse_interpolators(section, option)

        self._cfg.set(section, option, value)

    def get(self,
            section: str,
            option: str) -> str:
        """Retrieve value set to (section, option)."""

        section, option = self._traverse_interpolators(section, option)

        return self._cfg.get(section, option)


class DynamicConfigSection:
    """Simplify option references in configuration file section

    Args:
        section: the section of the configuration file in which
            we want to store option keys and values.

    
    Note:
        This class dynamically sets properties using a classmethod.
        Class methods will carry out changes to all existing and subsequent
        class instances.  To understand how this manifests in this class
        consider the following example. Suppose that I have two class
        instances of DynamicConfigSection
        
          a = DynamicConfigSection("common")
          b = DynamicConfigSection("arbitrary")
        
        and suppose that I set the property for a as follows
        
          a._x = 25
          a._set_property("x")
          print(a.x)
              --> 25
        
        No suprise here.  Now let's consider class instance b
        
          hasattr(b, _x)
              --> False
          'x' in b.__dir__()
              --> True
        
        This is a bit odd.  Despite not defining `x` for instance `b`
        it does exist.  Now suppose we try to query `x` from `b`
        
          print(b.x)
              --> AttributeError, _x not defined
        
        Interesting, the error doesn't say that `x` isn't defined by instead
        the attribute that the property `x` references.
        
          b._x = 10
          print(b.x)
              --> 10
        
        Once we define `_x` attribute, the property defined in instance `a`
        works.  This is because `_x` is defined on the class instance, but
        not the class, meanwhile `x` is defined on the class, not just the
        instance.
    """
    def __init__(self, section: str) -> None:
        self._section: str = section
        self._dynamic_option_names: list[str] = []

    def __iter__(self):
        for w in self._dynamic_option_names:
            yield (w, getattr(self, w))

    @classmethod
    def _set_property(self, option: str) -> None:
        setattr(self,
                option,
                property(lambda self: getattr(self, self._opt_to_attr(option)),
                         lambda self, val: setattr(self, self._opt_to_attr(option), val)))

    @staticmethod
    def _opt_to_attr(option: str) -> str:
        return f"_{option}"

    def add_option(self, option: str, value: str) -> None:
        setattr(self, self._opt_to_attr(option), value)
        self._set_property(option)
        self._dynamic_option_names.append(option)


def get_config_section(config_filename: str,
                       section_name: str) -> DynamicConfigSection:

    cfg = ConfigParser()
    cfg.read(config_filename)

    if not cfg.has_section(section_name):
        raise ValueError(f"{section_name} not defined in configuration.")

    pars = DynamicConfigSection(section_name)

    for opt in cfg.options(section_name):
        pars.add_option(opt, cfg.get(section_name, opt))

    return pars


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
    cfg.set("common", "bin",            os.environ.get(_constants.ENV_BIN))
    cfg.set("query", "dbname",          os.environ.get(_constants.ENV_DB_NAME))
    cfg.set("query", "host",            os.environ.get(_constants.ENV_DB_HOST))
    cfg.set("query", "port",            os.environ.get(_constants.ENV_DB_PORT))
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


