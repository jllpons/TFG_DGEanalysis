#!/usr/bin/env python3

"""
Script for mapping gene IDs between different strains.

Given a .tsv file that serves as a map and a number of dataframes containing 
gene IDs, this script returns the dataframes containing only mapped gene IDs.

Returns:
    Multiple .tsv files that contain the mapped gene IDs for each input 
    dataframe specified in the configuration JSON file.

Dependencies:
    - pandas
    - dgeapy.utilities
"""

import os
import sys
import argparse
from dataclasses import dataclass

import pandas as pd

from dgeapy.utilities import read_config_json_file


@dataclass
class Strain:
    """Store data related to each strain to map"""
    df_filename: str
    df: pd.DataFrame
    id_col_in_map: str
    id_col_in_df: str




def main():

    description = """
    Script for mapping gene IDs between different strains.
    Uses a TSV file that serves as a map and a number of dataframes
    containing geneIDs. Returns the dataframes containing only mapped geneIDs.
    """

    parser = argparse.ArgumentParser(
            description=description, usage='dgeapy.py mapgenes <config.json>'
            )

    parser.add_argument(
            "configuration_json_file",
            metavar="<config.json>",
            nargs="?",
            default="",
            type=str,
            help="path to JSON configuration file",
            )
    parser.add_argument(
            '-k', '--key',
            metavar='STR',
            default='',
            type=str,
            help='this geneIDs will appear as "mapped_geneID" in the '       \
                 'generated dataframes otherwise strain geneIDs present in ' \
                 'map will be used'
                 )
    parser.add_argument(
            '-s', '--stats',
            action='store_true',
            default=False,
            help='generate a <stats.txt> file with the stats about the mapping'
            )
    parser.add_argument(
            '-u', '--unmapped',
            action='store_true',
            default=False,
            help='save to TSV files containing (1) orphan geneIDs and (2) ' \
                 'geneIDs that are mapped but not present in input dfs'
            )
    # TODO:
    #     - Add output options
    #     - Add override option

    args = parser.parse_args()

    # Config file stuff
    if not args.configuration_json_file:
        parser.print_help()
        sys.exit("\n** The JSON configuration file is required **\n")

    config_file = os.path.abspath(args.configuration_json_file)
    if not os.path.isfile(config_file):
        raise FileNotFoundError(f"Could not find file: {config_file}")

    config_dict = read_config_json_file(config_file)

    out_dir = f'{os.getcwd()}/dgeapy_map'
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

    map = pd.read_csv(config_dict['map'], sep='\t')
    # Load data into the dataclass
    strains = []
    for s in config_dict['strains']:
        strains.append(Strain(
                        df_filename=config_dict['strains'][s]['df'],
                        df = pd.read_csv(
                            config_dict['strains'][s]['df'], sep='\t'
                            ),
                        id_col_in_map=config_dict['strains'][s]['id_col_in_map'],
                        id_col_in_df=config_dict['strains'][s]['id_col_in_df'],
                        ))

    # The columns we care about
    strain_cols = [i.id_col_in_map for i in strains]
    if len(args.key):
        strain_cols.append(args.key)

    # Keep only rows in columns we care about that don't contain null values
    map_only_mapped = map[strain_cols].dropna()
    # Only keep rows in columns we care about that contain at
    # least one null value
    map_only_unmapped = map[map[strain_cols].isnull().any(axis=1)]

    # Keep only geneIDs that are present in dfs to map
    map_only_mapped_present_dfs = map_only_mapped
    for s in strains:
        map_only_mapped_present_dfs = map_only_mapped_present_dfs[
                map_only_mapped_present_dfs[s.id_col_in_map].isin(
                    s.df[s.id_col_in_df]
                    )]

    # Keep only mapped geneIDs that are not present in dfs to map
    map_only_mapped_notpresent_dfs = map_only_mapped[
            ~map_only_mapped[strain_cols].isin(
                map_only_mapped_present_dfs[strain_cols]
                )].dropna()

    strain_mapped_genes = []
    for s in strains:
        # New column in df, column name is the same that apper in the map,
        # values are the same that the ones in geneID column
        s.df[s.id_col_in_map] = s.df[s.id_col_in_df]

        if len(args.key):
            k = args.key
            df = pd.merge(
                    map_only_mapped_present_dfs[[s.id_col_in_map, k]],
                    s.df,
                    on=s.id_col_in_map
                    )
            df = df.rename(columns={k : 'mapped_geneID'})

        else:
            df = pd.merge(
                    map_only_mapped_present_dfs[s.id_col_in_map],
                    s.df,
                    on=s.id_col_in_map
                    )
            df = df.rename(columns={s.id_col_in_map : 'mapped_geneID'})

        df = df.set_index('mapped_geneID')
        df = df.drop(columns=s.id_col_in_df)
        df.to_csv(f'{out_dir}/{(s.df_filename.split("/"))[-1]}_mapped.tsv', sep='\t')

        strain_mapped_genes.append(len(df.index))

    if args.unmapped is True:
        map_only_unmapped.to_csv(f'{out_dir}/orphan_geneIDs.tsv', sep='\t')
        map_only_mapped_notpresent_dfs.to_csv(
                f'{out_dir}/mapped_notpresent_in_input_dfs.tsv', sep='\t'
                )

    if args.stats is True:
        strain_specific_text = []

        for i,s in enumerate(strains):
            s_text = f'\t- From {len(s.df.index)} genes in "{s.df_filename}", {strain_mapped_genes[i]} where mapped.\n\t'

            strain_specific_text.append(s_text)

        text = f"""
        dgeapy mapgenes stats:

        - {strain_mapped_genes[0]} genes have been mapped.
        {''.join(strain_specific_text)}

        - From {len(map.index)} genes provided in the map, {len(map_only_unmapped) + len(map_only_mapped_notpresent_dfs)} have not been mapped:
            - {len(map_only_unmapped)} genes did not have at least one geneID in one of the columns specified in the <config.json> file.
            - {len(map_only_mapped_notpresent_dfs.index)} geneIDs where present in all of the columns specified in the input map but could not be found in at least one of the input dataframes.
        """

        with open(f'{out_dir}/dgeapy_mapgenes_stats.txt', 'w') as handle:
            handle.write(text)

if __name__ == "__main__":
    main()

