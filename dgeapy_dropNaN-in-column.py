#!/usr/bin/env python3

"""
"""

import os
import sys
import argparse

import pandas as pd



def main():

    description = """The script loads the input table and identifies the specified column. It then removes rows that have NaN values in that column using the `notna()` method in pandas.

The modified table, without the rows containing NaN values, is saved as a new TSV and XLSX file with "_dropNaN" appended to the original file name.

    """

    parser = argparse.ArgumentParser(
            description=description, usage='dgeapy.py dropNaN-in-column [args]'
            )

    input = parser.add_argument_group('input')
    input.add_argument(
            '-t', '--table',
            metavar='<table.tsv>',
            type=str,
            help='input table'
                )
    parser.add_argument(
            '-c', '--column',
            metavar='STR',
            type=str,
            help='column where to drop NaN values'
            )

    args = parser.parse_args()

    if not args.table:
        parser.print_help()
        sys.exit("\n** The <table.tsv> file is required **\n")
    if not args.column:
        parser.print_help()
        sys.exit("\n** The column is required **\n")

    table_file = os.path.abspath(args.table)
    if not os.path.isfile(table_file):
        raise FileNotFoundError(f'Could not find file: {table_file}')


    if table_file.endswith('.xlsx'):
        table = pd.read_excel(table_file)
    else:
        table = pd.read_csv(table_file, sep='\t')

    table = table[table[args.column].notna()]


    table.to_csv(
            f'{table_file.split("/")[-1]}_dropNaN.tsv',
            sep='\t',
            index=False
            )
    table.to_excel(
            f'{table_file.split("/")[-1]}_dropNaN.xlsx',
            index=False
            )


if __name__ == "__main__":
    main()



