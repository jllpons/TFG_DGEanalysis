#!/usr/bin/env python3

"""
WIP
"""

import sys
import subprocess

from dgeapy.utilities import mkconfigs

def main():
    description = """
    dgeapy: Script for Differential Gene Expression dataframe analysis.

    usage: dgeapy.py <command>

        -h or --help:   print help message

        dataframe analyisis:
            2muts       WIP
            3muts       analyze a dataframe contaning 3 mutants

        utilities:
            mkconfigs   create config.json file samples
    """
    
    argv_len = len(sys.argv)
    if argv_len == 1:
        print(description)

    if argv_len > 1:
        cmd = sys.argv[1]

        if cmd == "-h" or "--help":
            print(description)

        elif cmd == "3muts":
            subcmd = ["dgeapy_3muts.py"] + sys.argv[2:]
            subprocess.call(subcmd)

        if cmd == "mkconfigs":
            mkconfigs()

if __name__ == "__main__":
    main()

