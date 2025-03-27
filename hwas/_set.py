"""Interface for configuration parameters


Functions:
    list_pars: list all parameters that are and can be set
"""

import os
import configparser


def list_pars():
    
    if not os.path.isfile("config"):
        pass
