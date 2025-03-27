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

    def __init__(self) -> None:
        self._cfg = configparser.ConfigParser(
                            interpolation = configparser.ExtendedInterpolation(),
                            allow_no_value = True)

    def read(self, object_name: typing.Any) -> None:
        self._cfg.read(object_name)

    def get_option_interpolator(self,
                                section: str,
                                option: str) -> tuple[str] | None:

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

    def set(self,
            section: str,
            option: str,
            value: str | None) -> None:

        while self.is_interpolation(section, option):
            section, option = self.get_option_interpolator(section, option)

        self._cfg.set(section, option, value)

    def get(self,
            section: str,
            option: str) -> str:

        while self.is_interpolation(section, option):
            section, option = self.get_option_interpolator(section, option)

        return self._cfg.get(section, option)


class DynamicConfigSection:
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


