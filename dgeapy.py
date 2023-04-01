#!/usr/bin/env python3

"""
dgeapy: a script that tries to analyse Differential Gene Expression (DGE)
data at a different levels.
"""

import os
import sys
import subprocess

from dgeapy.utilities import mkconfigs


def main():

    description = """
dgeapy: Script for Differential Gene Expression data analyisis at different levels.

usage: dgeapy.py <command> [options]

    dataframe analyisis:
        multiplemuts    analyze a dataframe contaning 3 mutants

    utilities:
        mapgenes        map geneIDs using a <map.tsv> file
        joindfs         perform a left join on two tables using common column
        mkconfigs       create config.json file samples

    options:
        -h, --help
    """

    arg_len = len(sys.argv)
    if arg_len == 1:
        print(description)

    elif arg_len > 1:
        cmd = sys.argv[1]
        dgeapy_path = os.path.dirname(os.path.realpath(__file__))

        if cmd == "-h" or cmd == "--help":
            print(description)

        elif cmd == "multiplemuts":
            subcmd = ["python", f"{dgeapy_path}/dgeapy_multiplemuts.py",] + sys.argv[2:]
            subprocess.run(subcmd)

        elif cmd == "mapgenes":
            subcmd = ["python", f"{dgeapy_path}/dgeapy_mapgenes.py",] + sys.argv[2:]
            subprocess.run(subcmd)

        elif cmd == "joindfs":
            subcmd = ["python", f"{dgeapy_path}/dgeapy_joindfs.py",] + sys.argv[2:]
            subprocess.run(subcmd)

        elif cmd == "mkconfigs":
            mkconfigs()

        else:
            print(description)


if __name__ == "__main__":
    main()

