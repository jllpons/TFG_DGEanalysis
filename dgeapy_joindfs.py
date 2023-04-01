#!/usr/bin/env python3

"""
Perform a left join of two tables based on a common key column. The first
table is provided as a TSV or XLSX file using the '-t1' or '--table1' option,
and the second table is provided in the same way using the '-t2' or '--table2' option.
The name of the key column in both tables is specified using the '-on' option.

The output consists of two files: the first file contains the joined table,
with all columns from both input tables. The second file contains the rows
from the second table that did not have a match in the first table.
"""

import os
import sys
import argparse

import pandas as pd


def main():

    description = """
    Perform a left join of two tables based on a common key column. The first
    table is provided as a TSV or XLSX file using the '-t1' or '--table1' option,
    and the second table is provided in the same way using the '-t2' or '--table2' option.
    The name of the key column in both tables and is specified using the '-on' option.
    The output consists of two files: the first file contains the joined table,
    with all columns from both input tables. The second file contains the rows
    from the second table that did not have a match in the first table.
    """

    parser = argparse.ArgumentParser(
            description=description, usage='dgeapy.py joindfs [args]'
            )

    input = parser.add_argument_group('input')
    input.add_argument(
            '-t1', '--table1',
            metavar='<table_1.tsv>',
            type=str,
            help='first table, where the join will be applied'
                 )
    input.add_argument(
            '-t2', '--table2',
            metavar='<table_2.tsv>',
            type=str,
            help='second table, matching keys in this table will be joined'
                 )
    input.add_argument(
            '-on',
            metavar='STR',
            type=str,
            help='column on <table_1.tsv> and <table_2.tsv> that will be ' \
                 'taken as key for the join'
                 )
    # TODO:
    #     - Add output options
    #     - Add override option

    args = parser.parse_args()

    if not args.table1:
        parser.print_help()
        sys.exit("\n** The <table_1.tsv> file is required **\n")
    if not args.table2:
        parser.print_help()
        sys.exit("\n** The <table_2.tsv> file is required **\n")
    if not args.on:
        parser.print_help()
        sys.exit("\n** <map.tsv> and <table.tsv> column to use as key is required **\n")

    table_1_file = os.path.abspath(args.table1)
    table_2_file = os.path.abspath(args.table2)
    if not os.path.isfile(table_1_file):
        raise FileNotFoundError(f'Could not find file: {table_1_file}')
    if not os.path.isfile(table_2_file):
        raise FileNotFoundError(f'Could not find file: {table_2_file}')

    if table_1_file.endswith('.xlsx'):
        table_1 = pd.read_excel(table_1_file)
    else:
        table_1 = pd.read_csv(table_1_file, sep='\t')

    if table_2_file.endswith('.xlsx'):
        table_2 = pd.read_excel(table_2_file)
    else:
        table_2 = pd.read_csv(table_2_file, sep='\t')

    key = args.on

    table_1 = table_1.set_index(key)
    table_2 = table_2.set_index(key)

    out_dir = f'{os.getcwd()}/dgeapy_joindfs_output'
    # Crate a directory for the output. If already exists, add _n to the name.
    out_dir_accumulator = 1
    if os.path.isdir(out_dir):
        out_dir_n = f"{out_dir}_{out_dir_accumulator}"
        while os.path.isdir(out_dir_n):
            out_dir_accumulator += 1
            out_dir_n = f"{out_dir}_{out_dir_accumulator}"
        out_dir = out_dir_n
        os.mkdir(out_dir)
    else:
        os.mkdir(out_dir)

    table_1_mapped = table_1.join(table_2, how='left')

    table_1_columns = table_1.columns.values.tolist()
    table_2_columns = table_2.columns.values.tolist()
    table_1_mapped = table_1_mapped.reindex(columns= table_1_columns + table_2_columns)

    table_2_notmapped = table_2[~table_2.index.isin(table_1_mapped.index)]

    table_1_mapped = table_1_mapped.reset_index()
    table_1_mapped = table_1_mapped.loc[:, ~table_1_mapped.columns.str.contains('^Unnamed')]
    table_1_mapped.to_csv(f'{out_dir}/{table_1_file.split("/")[-1]}_joined.tsv', sep='\t', index=False)
    table_1_mapped.to_excel(f'{out_dir}/{table_1_file.split("/")[-1]}_joined.xlsx', index=False)

    table_2_notmapped = table_2_notmapped.reset_index()
    table_2_notmapped = table_2_notmapped.loc[:, ~table_2_notmapped.columns.str.contains('^Unnamed')]
    table_2_notmapped.to_csv(f'{out_dir}/{table_2_file.split("/")[-1]}_notjoined.tsv', sep='\t', index=False)
    table_2_notmapped.to_excel(f'{out_dir}/{table_2_file.split("/")[-1]}_notjoined.xlsx', index=False)

if __name__ == "__main__":
    main()

