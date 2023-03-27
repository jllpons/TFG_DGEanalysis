#!/usr/bin/env python3

"""Utilities for dgeappy"""

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
    json_mapgenes = {
                'map' : 'path/to/yourmap.tsv',
                'strains' : {
                    'yourstrain_1' : {
                        'df' : 'path/to/yourstrain_1.tsv',
                        'id_col_in_map' : 'Column name in <map.tsv> for this strain',
                        'id_col_in_df' : 'gene IDs column name in <yourstrain_1.tsv>',
                        },
                    'yourstrain_N' : {
                        'df' : 'path/to/yourstrain_N.tsv',
                        'id_col_in_map' : 'Column name in <map.tsv> for this strain',
                        'id_col_in_df' : 'gene IDs column name in <yourstrain_N.tsv>',
                        },
                    }
                }

    with open(f"{os.getcwd()}/multiplemuts_sample_config.json", "w") as j:
        j.write(json.dumps(json_multiplemuts, indent=4))
    with open(f"{os.getcwd()}/mapgenes_sample_config.json", "w") as j:
        j.write(json.dumps(json_mapgenes, indent=4))


def read_config_json_file(config_file):
    """Reads JSON configuration file.
    """
    # TODO: Read the file even if it has syntax errors.

    with open(config_file, "r") as j:
        config_dictionary = json.loads(j.read())

    return config_dictionary

