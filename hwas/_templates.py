"""Interface with templates
"""

import string
import pathlib
import os
import importlib.resources

from . import _constants



class HwasTemplate(string.Template):
    delimiter = '@'


def get_template_filename(template_name: str) -> str:

    if ((template_path := importlib.resources
                        .files(_constants.TEMPLATES_PATH)
                        .joinpath(template_name))
        and template_path.is_file()):

        return str(template_path)

    raise FileNotFoundError(str(template_path))


def render(template_path: str, 
           options: dict[str,str]) -> str:

    if (os.path.isfile(template_path) is None):
        raise FileNotFoundError(f"Template {template_path} not found.")

    with open(template_path, 'r') as fid:
        output = HwasTemplate(fid.read())


    for identifier in output.get_identifiers():

        if identifier not in options:
            raise KeyError(f"Template file, {template_path}, identifier is"
                           " not found among the provided option values.")


    return output.substitute(**options)


# def write_templated_string() -> None:
#     hwas_sbatch_template = (importlib.resources
#                             .files(_constants.DEFAULT_TEMPLATE_PATH)
#                             .joinpath(_constants.SBATCH_HWAS_TEMPLATE))
# 
#     if len(hwas_sbatch_template) != 1:
#         raise FileNotFoundError("Can't find template file"
#                                 f" {_constants.SBATCH_HWAS_TEMPLATE}")
# 
# 
#     with (open(os.path.join(path, _constants.SBATCH_SCRIPT_HWAS), "w") as fout,
#         open(importlib.resources
#              .files(_constants.TEMPLATE_PATH)
#              .joinpath(_constants.SBATCH_HWAS_TEMPLATE), "r") as fin):
#             pass
