#!/usr/bin/env python3

"""License?
"""

import sys
import subprocess

from dgeapy.utilities import mkconfigs


def main():

    description = """
dgeapy: Script for Differential Gene Expression data analyisis from a dataframe.

usage: dgeapy.py <command> [options]

    dataframe analyisis:
        multiplemuts    analyze a dataframe contaning 3 mutants

    utilities:
        mkconfigs       create config.json file samples

    options:
        -h, --help
    """

    arg_len = len(sys.argv)
    if arg_len == 1:
        print(description)

    elif arg_len > 1:
        cmd = sys.argv[1]

        if cmd == "-h" or cmd == "--help":
            print(description)

        elif cmd == "multiplemuts":
            subcmd = ["python3", "dgeapy_multiplemuts.py",] + sys.argv[2:]
            subprocess.call(subcmd)

        elif cmd == "mkconfigs":
            mkconfigs()

        else:
            print(description)


if __name__ == "__main__":
    main()

