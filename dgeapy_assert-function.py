#!/usr/bin/env python3

"""
"""

import os
import sys
import argparse
import pdb

import numpy as np
import pandas as pd
from pandas.core.api import notnull

from dgeapy.utilities import read_config_json_file


def main():

    description = """The script reads the JSON file to extract the function names and their associated GO codes and KEGG pathways. It then loads the input table and iterates over each row, assigning functions to entries based on the provided annotations. The annotations can be GO codes, GO ancestors, or KEGG pathways.

The annotated table is saved as a new TSV and XLSX file with "_functions" appended to the original file name.
    """

    parser = argparse.ArgumentParser(
            description=description, usage='dgeapy.py assert-function [args]'
            )

    input = parser.add_argument_group('input')
    input.add_argument(
            '-j', '--json',
            metavar='<config.json>',
            type=str,
            help='JSON file containing function name as key and a list of ' \
                 'all the selected GO codes as value'
                 )
    input.add_argument(
            '-t', '--table',
            metavar='<table.tsv>',
            type=str,
            help='input table'
                )
    parser.add_argument(
            '-s', '--strict',
            action='store_true',
            default=False,
            help='only consider "is_a" GO relationships, on default both ' \
                 '"is_a" and "regulates" are considered'
            )

    args = parser.parse_args()

    if not args.json:
        parser.print_help()
        sys.exit("\n** The <config.json> file is required **\n")
    if not args.table:
        parser.print_help()
        sys.exit("\n** The <table.tsv> file is required **\n")

    json_file = os.path.abspath(args.json)
    table_file = os.path.abspath(args.table)
    if not os.path.isfile(json_file):
        raise FileNotFoundError(f'Could not find file: {json_file}')
    if not os.path.isfile(table_file):
        raise FileNotFoundError(f'Could not find file: {table_file}')

    function_dict = read_config_json_file(json_file)

    if table_file.endswith('.xlsx'):
        table = pd.read_excel(table_file)
    else:
        table = pd.read_csv(table_file, sep='\t')

    go_codes_column = 'Gene Ontology IDs'
    ancestors_column = 'go_ancestors_is_a_and_regulates'
    kegg_pathways_column = 'kegg_pathways'
    if args.strict is True:
        ancestors_column = 'go_ancestors_is_a'

    table.insert(
            loc=(len(table.columns.values.tolist())),
            column='function',
            value=[np.nan] * table.shape[0]
            )

    for index, row in table.iterrows():
        if pd.isnull(row[go_codes_column]) and pd.isnull(row[ancestors_column]) and pd.isnull(row[kegg_pathways_column]):
            table.at[index, 'function'] = 'not annotated'

        if pd.notnull(row[go_codes_column]):
            go_codes = str(row[go_codes_column]).split(';')

            for code in go_codes:
                code = code.strip()
                for key in function_dict:
                    if code in function_dict[key]['GO']:
                        table.at[
                                index,'function'
                                ] = key
                        break
                if pd.notnull(table.at[index, 'function']):
                    break

        if pd.notnull(row[ancestors_column]):
            if pd.isnull(table.at[index, 'function']):
                ancestor_codes = str(row[ancestors_column]).split(';')

                for ancestor in ancestor_codes:
                    for key in function_dict:
                        if ancestor in function_dict[key]['GO']:
                            table.at[
                                    index,'function'
                                    ] = key
                            break
                    if pd.notnull(table.at[index, 'function']):
                        break

        if pd.notnull(row[kegg_pathways_column]):
            if pd.isnull(table.at[index, 'function']):
                kegg_pathways = str(row[kegg_pathways_column])

                for key in function_dict:
                    for pathway in function_dict[key]['KEGG']:
                        if pathway in kegg_pathways:
                            table.at[
                                    index,'function'
                                    ] = key
                            break
                    if pd.notnull(table.at[index, 'function']):
                        break
        if pd.isnull(table.at[index, 'function']):
            table.at[index, 'function'] = 'others'


    table.to_csv(
            f'{table_file.split("/")[-1]}_functions.tsv',
            sep='\t',
            index=False
            )
    table.to_excel(
            f'{table_file.split("/")[-1]}_functions.xlsx',
            index=False
            )


if __name__ == "__main__":
    main()



