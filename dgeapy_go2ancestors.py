#!/usr/bin/env python3

import os
import sys
import argparse

import numpy as np
import pandas as pd
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag


# <http://purl.obolibrary.org/obo/go.obo>
GODAG = GODag(
        # path changed for GitHub
        '../data/gene_onthology_obo/go.obo',
        optional_attrs={'relationship'}
        )

def go2ancestors_is_a(go_code: str) -> list:
    """
    """

    gosubdag = GoSubDag([go_code], GODAG, prt=None)

    return list(gosubdag.rcntobj.go2ancestors[go_code])


def go2ancestors_is_a_and_regulates(go_code: str) -> list:
    """
    """

    optional_relationships = {
        'regulates', 'negatively_regulates', 'positively_regulates'
        }

    gosubdag = GoSubDag(
            [go_code],
            GODAG,
            relationships=optional_relationships,
            prt=None
            )

    return list(gosubdag.rcntobj.go2ancestors[go_code])


def main():

    description = """The script reads the input table, loads the provided GO ontology file, and defines two functions to retrieve ancestors of GO codes. It then iterates over each row in the table, splits the GO codes in the specified column, and retrieves the ancestors for each code. The ancestors are stored in new columns, 'go_ancestors_is_a' and 'go_ancestors_is_a_and_regulates', in the table.
    """

    parser = argparse.ArgumentParser(
            description=description, usage='dgeapy.py go2ancestors [args]'
            )

    input = parser.add_argument_group('input')
    input.add_argument(
            '-t', '--table',
            metavar='<table.tsv>',
            type=str,
            help='input table'
                 )
    input.add_argument(
            '-on',
            metavar='STR',
            type=str,
            help='column on <table.tsv> that contains GO ' \
                 'codes for each gene'
                 )

    args = parser.parse_args()

    if not args.table:
        parser.print_help()
        sys.exit("\n** The <table.tsv> file is required **\n")
    if not args.on:
        parser.print_help()
        sys.exit("\n** Name of the column contains GO codes is required **\n")

    table_file = os.path.abspath(args.table)
    if not os.path.isfile(table_file):
        raise FileNotFoundError(f'Could not find file: {table_file}')

    if table_file.endswith('.xlsx'):
        table = pd.read_excel(table_file)
    else:
        table = pd.read_csv(table_file, sep='\t')

    go_column = args.on

    table_columns = table.columns.values.tolist()

    for i in range(len(table_columns)):
        if table_columns[i] == go_column:
            table.insert(
                    loc=(i+1),
                    column='go_ancestors_is_a',
                    value=[np.nan] * table.shape[0]
                    )
            table.insert(
                    (i+2),
                    'go_ancestors_is_a_and_regulates',
                    value=[np.nan] * table.shape[0]
                    )

    for index, row in table.iterrows():
        if row['Gene Ontology IDs'] is not np.nan:
            go_codes = str(row['Gene Ontology IDs']).split(';')
            go_ancestors_is_a = []
            go_ancestors_is_a_and_regulates = []

            for code in go_codes:
                try:
                    go_ancestors_is_a.extend(go2ancestors_is_a(code.strip()))
                except KeyError: # obsolete GO code
                    pass
                try:
                    go_ancestors_is_a_and_regulates.extend(
                            go2ancestors_is_a_and_regulates(code.strip())
                            )
                except KeyError: # obsolete GO code
                    pass

            table.at[index, 'go_ancestors_is_a'] = ";".join(set(go_ancestors_is_a))
            table.at[index, 'go_ancestors_is_a_and_regulates'] = ";".join(set(go_ancestors_is_a_and_regulates))

    table.to_csv(
            f'{table_file.split("/")[-1]}_go_ancestors.tsv',
            sep='\t',
            index=False
            )
    table.to_excel(
            f'{table_file.split("/")[-1]}_go_ancestors.xlsx',
            index=False
            )


if __name__ == '__main__':
    main()



