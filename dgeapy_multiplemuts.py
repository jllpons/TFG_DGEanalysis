#!/usr/bin/env python3

"""
"""

import os
import sys
import argparse
from dataclasses import dataclass

import pandas as pd

import dgeapy


@dataclass
class Sampledata:
    """Store all of the data related to each sample."""
    name: str
    input_df: pd.DataFrame
    df_columns:  dict
    dge_df: pd.DataFrame
    up_df: pd.DataFrame
    down_df: pd.DataFrame


def main(args=None):

    description = "Differential Gene Expression data analysis from multiple " \
                  "dataframes containing mutant vs. wild type experiments."

    parser = argparse.ArgumentParser(
                        description=description,
                        usage="dgeapy.py multiplemuts <config.json>"
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
            "--pvalue",
            nargs=1,
            default=0.05,
            type=float,
            help="p-value threshold, default is 0.05",
            )
    parser.add_argument(
            "--padj",
            nargs=1,
            default=0.05,
            type=float,
            help="adjusted p-value threshold, default is 0.05",
            )
    parser.add_argument(
            "--fc",
            nargs=1,
            default=2.00,
            type=float,
            help="fold change threshold, default is 2.00",
            )
    parser.add_argument(
            "--formats",
            nargs="?",
            default=["png"],
            type=str,
            action="append",
            help="plot formats, defalut is png",
            )
    parser.add_argument(
            "-g",
            action='store_true',
            default=False,
            help="include non-coding transcripts"
            )

    args = parser.parse_args()

    if not args.configuration_json_file:
        parser.print_help()
        sys.exit("\n** The JSON configuration file is required **\n")

    config_file = os.path.abspath(args.configuration_json_file)
    if not os.path.isfile(config_file):
        raise FileNotFoundError(f"Could not find file: {config_file}")

    # For the moment, read the JSON file
    # TODO: make the funcntion work even if the JSON file has invalid syntax.
    config_dictionary = dgeapy.read_config_json_file(config_file)

    FOLD_CHANGE_THRESHOLD = args.fc
    P_VALUE_THRESHOLD = args.pvalue
    PADJ_THRESHOLD = args.padj
    PLOT_FORMATS = args.formats
    DATAFRAMES = config_dictionary
    include_novels = args.g

    if len(DATAFRAMES) != 3:
        sys.exit('\n** dgeapy multiplemuts just supports data from 3' \
                 ' different samples **\n')

    for k in DATAFRAMES:
        if not os.path.isfile(DATAFRAMES[k]):
                raise FileNotFoundError(f"Could not find file: {DATAFRAMES[k]}")

    cwd = os.getcwd()
    output_dir = f"{cwd}/dgeapy_output"

    # Crate a directory for the output. If already exists, add _n to the name.
    output_dir_accumulator = 1
    if os.path.isdir(output_dir):
        output_dir_n = f"{output_dir}_{output_dir_accumulator}"
        while os.path.isdir(output_dir_n):
            output_dir_accumulator += 1
            output_dir_n = f"{output_dir}_{output_dir_accumulator}"
        output_dir = output_dir_n
        os.mkdir(output_dir)
    else:
        os.mkdir(output_dir)

    output_dirs_dict = {
            "df" : f"{output_dir}/dataframes",
            "volcano" : f"{output_dir}/volcano_plots",
            "venn" : f"{output_dir}/venn_diagrams",
            "upset" : f"{output_dir}/upset_plots",
            "sankey" : f"{output_dir}/sankey_diagrams"
            }
    for k in output_dirs_dict:
        os.mkdir(output_dirs_dict[k])

    data = []
    for k in DATAFRAMES:

        sample_df_dir = f'{output_dirs_dict["df"]}/{k}'
        os.mkdir(sample_df_dir)

        df = pd.read_csv(
                    DATAFRAMES[k],
                    sep="\t",
                    na_values=["--", "",]
                    )

        if include_novels is False:
            df = df[~df.gene_id.str.contains("Novel")]
            df = df[~df.gene_id.str.contains("sRNA")]

        df = df.set_index('gene_id')

        dgeapy.add_fold_change_columns(df)
        dgeapy.add_regulation_columns(df)

        # Dictionary containing names of the columnns of:
        # FC, log2FC, p-value and padj values.
        column_names_to_check = dgeapy.column_names_to_check(df)

        # DEG according to FC, pvale and padj values
        dge_df = dgeapy.filter_FC_PVALUE_PADJ(
                            dataframe=df,
                            foldchange_threshold=FOLD_CHANGE_THRESHOLD,
                            foldchange_column_name=column_names_to_check["FoldChange"],
                            p_value_threshold=P_VALUE_THRESHOLD,
                            p_value_column_name=column_names_to_check["pvalue"],
                            padj_threshold=PADJ_THRESHOLD,
                            padj_column_name=column_names_to_check["padj"],
                            )

        # Up and Down regulated genes from DEG
        up_df = dge_df[dge_df[column_names_to_check["Regulation"]] == "Up"]
        donw_df = dge_df[dge_df[column_names_to_check["Regulation"]] == "Down"]

        sample_data = Sampledata(
                        name=k,
                        input_df=df,
                        df_columns=column_names_to_check,
                        dge_df=dge_df,
                        up_df=up_df,
                        down_df=donw_df
                        )
        data.append(sample_data)

        # Saving to csv and excel files.
        sample_data.input_df.to_csv(
                                f'{sample_df_dir}/{sample_data.name}_input.tsv',
                                sep="\t",
                                )
        sample_data.input_df.to_excel(
                                f'{sample_df_dir}/{sample_data.name}_input.xlsx',
                                )
        sample_data.dge_df.to_csv(
                                f'{sample_df_dir}/{sample_data.name}_DEG.tsv',
                                sep="\t",
                                )
        sample_data.dge_df.to_excel(
                                f'{sample_df_dir}/{sample_data.name}_DEG.xlsx',
                                )
        sample_data.up_df.to_csv(
                                f'{sample_df_dir}/{sample_data.name}_UP.tsv',
                                sep="\t",
                                )
        sample_data.up_df.to_excel(
                                f'{sample_df_dir}/{sample_data.name}_UP.xlsx',
                                )
        sample_data.down_df.to_csv(
                                f'{sample_df_dir}/{sample_data.name}_DOWN.tsv',
                                sep="\t",
                                )
        sample_data.dge_df.to_excel(
                                f'{sample_df_dir}/{sample_data.name}_DOWN.xlsx',
                                )

        # Generate a volcano and a count plots
        dgeapy.generate_volcano_plot(
                data=sample_data,
                file_path=output_dirs_dict['volcano'],
                foldchange_threshold=FOLD_CHANGE_THRESHOLD,
                padj_threshold=PADJ_THRESHOLD,
                plot_formats=PLOT_FORMATS,
                )

    # Sankey diagrams
    dgeapy.generate_sankey_diagram(
            data,
            fc_value=FOLD_CHANGE_THRESHOLD,
            padj_value=PADJ_THRESHOLD,
            plot_formats=PLOT_FORMATS,
            path=output_dirs_dict['sankey']
            )

    # For the 3 sets of gene IDs for DEG, UP and DOWN regulated genes:
    #   - Generate 2 venn's diagrams representing the intersections of
    #     gene IDs. One will be defalut, the other will be unweight.
    #   - Generete an upset plot for also representing the intersections
    #   - From the intersections represented, generate  a dataframe for each
    #     containing all of the relevant information and save it to a file.
    dgeapy.mk_venn_upset_and_intersections_dfs(
            data=data,
            plot_formats=PLOT_FORMATS,
            venn_path=output_dirs_dict["venn"],
            upset_path=output_dirs_dict["upset"],
            df_path=output_dirs_dict["df"],
            )

    # Both up/down_regulation_labels are dictionaries conaining
    # the labels for the next venn diagrams we're going to generate.
    # They'll display the actual number of genes considered
    # up and down regulated at the same time.
    if len(data) == 3:
        up_regulation_labels = dgeapy.get_gene_ids_set_for_intersections(
                set1=set(data[0].up_df.index),
                set2=set(data[0].up_df.index),
                set3=set(data[0].up_df.index),
                )
        down_regulation_labels = dgeapy.get_gene_ids_set_for_intersections(
                set1=set(data[0].up_df.index),
                set2=set(data[0].up_df.index),
                set3=set(data[0].up_df.index),
                )
        # Generate the same two diagrams but with the labels
        dgeapy.generate_venn3_diagram_with_regulation_labels(
                mutant1_gene_set=set(data[0].dge_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].dge_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].dge_df.index),
                mutant3_name=data[2].name,
                plot_formats=PLOT_FORMATS,
                up_regulation_labels=up_regulation_labels,
                down_regulation_labels=down_regulation_labels,
                title="Differentially expressed genes",
                file_path=f"{output_dirs_dict['venn']}/venn_DEG_labels",
                )

    # Comparing sets for the 6 possible inverted regulations combinations.
    # Makes a venn diagramm and generates 7 dataframes for the genes of each
    # intersection.
    inverted_reg_venn_dir = f"{output_dirs_dict['venn']}/inverted_regulations"
    inverted_reg_upset_dir = f"{output_dirs_dict['upset']}/inverted_regulations"
    os.mkdir(inverted_reg_venn_dir)
    os.mkdir(inverted_reg_upset_dir)
    dgeapy.get_inverted_regulations_and_mk_venns_and_dataframes(
            data=data,
            plot_formats=PLOT_FORMATS,
            venn_directory_path=inverted_reg_venn_dir,
            upset_directory_path=inverted_reg_upset_dir,
            dataframes_directory_path=output_dirs_dict['df'],
            )

if __name__ == "__main__":
    main()

