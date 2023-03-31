#!/usr/bin/env python3

"""
Perform mapping between a mapping file <map.tsv> and a table file <table.tsv>
based on specified columns. The mapped data is saved to output files in a
directory created by the script.
"""

import os
import sys
import argparse

import numpy as np
import pandas as pd


def main():

    description = """Perform mapping between a mapping file <map.tsv> and
    a table file <table.tsv> based on specified columns. The mapped data is
    saved to output files in a directory created by the script.
    """

    parser = argparse.ArgumentParser(
            description=description, usage='dgeapy.py mapgenes [args]'
            )

    input = parser.add_argument_group('input')
    input.add_argument(
            '-m', '--map',
            metavar='<map.tsv>',
            type=str,
            help='geneIDs column name in <map.ysv> that  will appear as ' \
                 )
    input.add_argument(
            '-c', '--map-columns',
            metavar='STR STR',
            type=str,
            nargs='+',
            default=[],
            help='columns to add, separated by commas and no '
                 )
    input.add_argument(
            '-t', '--table',
            metavar='<table.tsv>',
            type=str,
            help='table'
                 )
    input.add_argument(
            '--on',
            metavar='STR',
            type=str,
            help='column name and values on <map.tsv> and <table.tsv> column'
                 )
    # TODO:
    #     - Add output options
    #     - Add override option

    args = parser.parse_args()

    if not args.map:
        parser.print_help()
        sys.exit("\n** The <map.tsv> file is required **\n")
    if not args.table:
        parser.print_help()
        sys.exit("\n** The <table.tsv> file is required **\n")
    if not args.on:
        parser.print_help()
        sys.exit("\n** <map.tsv> and <table.tsv> column value is required **\n")

    map_file = os.path.abspath(args.map)
    table_file = os.path.abspath(args.table)
    if not os.path.isfile(map_file):
        raise FileNotFoundError(f'Could not find file: {map_file}')
    if not os.path.isfile(table_file):
        raise FileNotFoundError(f'Could not find file: {table_file}')

    map = pd.read_csv(map_file, sep='\t')
    table = pd.read_csv(table_file, sep='\t')

    key = args.on
    if args.map_columns:
        map = map[args.map_columns]

    map['index'] = np.where(map[key].str.endswith('#'), map[key].str[:-1], map[key])
    map = map.set_index('index')
    table = table.set_index(key)

    out_dir = f'{os.getcwd()}/dgeapy_map_output'
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

    table_mapped = table.join(map, how='left')

    table_columns = table.columns.values.tolist()
    map_columns = map.columns.values.tolist()
    table_mapped = table_mapped.reindex(columns= map_columns + table_columns)

    map_notmapped = map[~map.index.isin(table_mapped.index)]

    table_mapped = table_mapped.reset_index()
    table_mapped.to_csv(f'{out_dir}/{table_file.split("/")[-1]}_mapped.tsv', sep='\t', index=False)
    table_mapped.to_excel(f'{out_dir}/{table_file.split("/")[-1]}_mapped.xlsx', index=False)

    map_notmapped = map_notmapped.reset_index()
    map_notmapped = map_notmapped.drop(columns=['index'])
    map_notmapped.to_csv(f'{out_dir}/{map_file.split("/")[-1]}_notmapped.tsv', sep='\t', index=False)
    map_notmapped.to_excel(f'{out_dir}/{map_file.split("/")[-1]}_notmapped.xlsx', index=False)

if __name__ == "__main__":
    main()

