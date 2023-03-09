#!/usr/bin/env python3

"""Utilities WIP
"""

import os
import json



#-------# Function definitions #-----------------------------------------------#


def mkconfigs():
    """Creates sample JSON configuration files.
    """

    json_multiplemuts = {
            "mutA" : "path/to/mutA.tsv",
            "mutB" : "path/to/mutB.tsv",
            "mutC" : "path/to/mutC.tsv"
            }

    with open(f"{os.getcwd()}/multiplemuts_sample_config.json", "w") as j:
        j.write(json.dumps(json_multiplemuts, indent=4))


def read_config_json_file(config_file):
    """Reads JSON configuration file.
    """
    # TODO: Read the file even if it has syntax errors.

    with open(config_file, "r") as j:
        config_dictionary = json.loads(j.read())

    return config_dictionary

