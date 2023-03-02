#!/usr/bin/env python3

import os
import json

"""
Utilities like create a directory with a given name or (...)
"""


#-------# Function definitions #-----------------------------------------------#

def mk_new_dir(new_dir_name):
    """
    Create a new directory with a given name. Also returns the newly created
    directory path as a string.
    """

    # Get the current working directory path
    cwd_path = os.getcwd()
    # Create the directory name
    new_dir_path = f"{cwd_path}/{new_dir_name}"

    # If the directory already exists, add an "_n" extension to the name.
    # (n represents an integrer)
    if os.path.isdir(new_dir_path) == False:
        os.mkdir(new_dir_path)
    else:
        accumulator = 1
        new_dir_path = f"{cwd_path}/{new_dir_name}_{accumulator}"
        # While this directory already exists add +1 to the _n extension.
        while os.path.isdir(new_dir_path):
            accumulator += 1
            new_dir_path = f"{cwd_path}/{new_dir_name}_{accumulator}"
        # Actually create the directory if the constructed name doesn't exist.
        os.mkdir(new_dir_path)

    # Returning the path as a string.
    return new_dir_path

def mkconfigs():
    """Creates sample JSON configuration files.
    """

    json_multiplemuts = {
            "dataframe_file_name" : "path/to/df.csv or .xls",
            "fold_change_threshold" : 1.5,
            "p_value_threshold" : 0.05,
            "padj_threshold" : 0.05,
            "wild_type_sample" : "K12",
            "mutant_samples" : ["geneA", "geneB", "geneC"],
            "plot_formats" : ["pdf", "jpg", "svg"],
            "df_to_merge" : {
                "geneA" : "path/to/df_to_merge.csv or .xls"
                }
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

