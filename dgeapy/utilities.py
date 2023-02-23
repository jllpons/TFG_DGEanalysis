#!/usr/bin/env python3

import os

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
