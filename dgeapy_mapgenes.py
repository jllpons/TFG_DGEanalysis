#!/usr/bin/env python3

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


def mk_mapped_df(
        map: pd.DataFrame,
        strain_df : pd.DataFrame,
        index: str
        ) -> pd.DataFrame:
    """Merges two dataframes based on index."""

    df = pd.merge(map, strain_df, on=index)

    return df


def main():

    description = """
    Script for mapping gene IDs between different strains.
    Given a .tsv file that serves as a map and a number of dataframes
    containing geneIDs. Return the dataframes containing only mapped geneIDs.
    """

    parser = argparse.ArgumentParser(
            description=description,
            usage='dgeapy.py mapgenes <config.json>'
            )

    parser.add_argument(
            "configuration_json_file",
            metavar="<config.json>",
            nargs="?",
            default="",
            type=str,
            help="path to JSON configuration file",
            )
    # TODO:
    #     - Add option to use a column as alias for geneIDs
    #     - Add output options
    #     - Add override option
    #     - Print a list of unmapped genes??

    args = parser.parse_args()

    if not args.configuration_json_file:
        parser.print_help()
        sys.exit("\n** The JSON configuration file is required **\n")
    config_file = os.path.abspath(args.configuration_json_file)
    if not os.path.isfile(config_file):
        raise FileNotFoundError(f"Could not find file: {config_file}")
    # For the moment, read the JSON file
    # TODO: make the funcntion work even if the JSON file has invalid syntax.
    config_dict = read_config_json_file(config_file)

    map = pd.read_csv(config_dict['map'], sep='\t')
    strains = []
    for s in config_dict['strains']:
        strains.append(Strain(
                        df_filename=config_dict['strains'][s]['df'],
                        df = pd.read_csv(
                            config_dict['strains'][s]['df'],
                            sep='\t'
                            ),
                        id_col_in_map=config_dict['strains'][s]['id_col_in_map'],
                        id_col_in_df=config_dict['strains'][s]['id_col_in_df'],
                        ))

    strain_cols = [i.id_col_in_map for i in strains]
    map_wout_unmapped = map[strain_cols].dropna()

    for s in strains:
        map_wout_unmapped = map_wout_unmapped[
                map_wout_unmapped[s.id_col_in_map].isin(s.df[s.id_col_in_df])
                ]

    for s in strains:
        s.df[s.id_col_in_map] = s.df[s.id_col_in_df]
        df = mk_mapped_df(
                map_wout_unmapped[s.id_col_in_map],
                s.df,
                s.id_col_in_map
                )
        df = df.rename(columns={s.id_col_in_map : 'mapped_geneID'})
        df = df.set_index('mapped_geneID')
        df.to_csv(f'{s.df_filename}_mapped.tsv', sep='\t')


if __name__ == "__main__":
    main()

