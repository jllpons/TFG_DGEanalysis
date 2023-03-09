#!/usr/bin/env python3

"""Lincense¿?
"""

import os
import sys
import argparse

import pandas as pd

import dgeapy

#------------------------------------------------------------------------------#

def main():

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
            help="include genes annotated as novel gene or sRNA"
            )

    args = parser.parse_args()

    if not args.configuration_json_file:
        parser.print_help()
        sys.exit("\n ** The JSON configuration file is required **\n")

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
    MUTANTS = list(config_dictionary.keys())
    DATAFRAMES = config_dictionary
    include_novels = args.g

    for mut in MUTANTS:
        if not os.path.isfile(DATAFRAMES[mut]):
                raise FileNotFoundError(f"Could not find file: {DATAFRAMES[mut]}")

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

    output_df_dir = f"{output_dir}/dataframes"
    os.mkdir(output_df_dir)
    output_volcano_dir = f"{output_dir}/volcano_plots"
    os.mkdir(output_volcano_dir)
    output_venn_dir = f"{output_dir}/venn_diagrams"
    os.mkdir(output_venn_dir)
    output_upset_dir = f"{output_dir}/upset_plots"
    os.mkdir(output_upset_dir)
    output_sankey_dir = f"{output_dir}/sankey_diagrams"
    os.mkdir(output_sankey_dir)

    # Storing df containning DE, Up and Down regulated genes.
    mutant_regulations_dictionary = {}

    for mutant in MUTANTS:

        mutant_df_dir = f"{output_df_dir}/{mutant}"
        os.mkdir(mutant_df_dir)

        mutant_df = pd.read_csv(
                            DATAFRAMES[mutant],
                            sep="\t",
                            na_values=["--", "",]
                            )

        dgeapy.add_fold_change_columns(mutant_df)
        dgeapy.add_regulation_columns(mutant_df)

        if include_novels is False:
            mutant_df = mutant_df[~mutant_df.gene_id.str.contains("Novel")]
            mutant_df = mutant_df[~mutant_df.gene_id.str.contains("sRNA")]

        # Dictionary containing names of the columnns of:
        # FC, log2FC, p-value and padj values.
        column_names_to_check = dgeapy.column_names_to_check(mutant_df)

        # DEG according to FC, pvale and padj values
        mutant_df_deg = dgeapy.filter_FC_PVALUE_PADJ(
                                dataframe=mutant_df,
                                foldchange_threshold=FOLD_CHANGE_THRESHOLD,
                                foldchange_column_name=column_names_to_check["FoldChange"],
                                p_value_threshold=P_VALUE_THRESHOLD,
                                p_value_column_name=column_names_to_check["pvalue"],
                                padj_threshold=PADJ_THRESHOLD,
                                padj_column_name=column_names_to_check["padj"],
                                )

        # Up and Down regulated genes from DEG
        mutant_df_up = mutant_df_deg[mutant_df_deg[
                                        column_names_to_check["Regulation"]
                                                    ] == "Up"]
        mutant_df_down = mutant_df_deg[mutant_df_deg[
                                        column_names_to_check["Regulation"]
                                                    ] == "Down"]

        # Adding dfs of  DE, Up and Down regulated genes
        mutant_regulations_dictionary[mutant] = {
                "DEG" : mutant_df_deg,
                "Up" : mutant_df_up,
                "Down" : mutant_df_down,
                }

        # Saving to csv and excel files.
        mutant_df.to_csv(
                    f"{mutant_df_dir}/{mutant}.tsv",
                    sep="\t",
                    index=False,
                    )
        mutant_df.to_excel(
                    f"{mutant_df_dir}/{mutant}.xlsx",
                    index=False,
                    )
        mutant_df_deg.to_csv(
                        f"{mutant_df_dir}/{mutant}_DEG.tsv",
                        sep="\t",
                        index=False,
                        )
        mutant_df_deg.to_excel(
                        f"{mutant_df_dir}/{mutant}_DEG.xlsx",
                        index=False,
                        )
        mutant_df_up.to_csv(
                        f"{mutant_df_dir}/{mutant}_UP.tsv",
                        sep="\t",
                        index=False,
                        )
        mutant_df_up.to_excel(
                        f"{mutant_df_dir}/{mutant}_UP.xlsx",
                        index=False,
                        )
        mutant_df_down.to_csv(
                        f"{mutant_df_dir}/{mutant}_DOWN.tsv",
                        sep="\t",
                        index=False,
                        )
        mutant_df_down.to_excel(
                        f"{mutant_df_dir}/{mutant}_DOWN.xlsx",
                        index=False,
                        )

        # Generate a volcano and a count plots
        dgeapy.generate_volcano_plot(
                dataframe=mutant_df,
                file_path=f"{output_volcano_dir}/{mutant}",
                foldchange_threshold=FOLD_CHANGE_THRESHOLD,
                log2foldchange_column_name=column_names_to_check["log2FoldChange"],
                padj_threshold=PADJ_THRESHOLD,
                padj_column_name=column_names_to_check["padj"],
                mutant_name=mutant,
                plot_formats=PLOT_FORMATS,
                )

    # Dictionary containing a set of gene IDs for each different
    # mutant gene regulations.
    mut_reg_GeneIdSet_dict = {}
    for mut in mutant_regulations_dictionary:

        mut_reg_GeneIdSet_dict[mut] = {}

        for reg in mutant_regulations_dictionary[mut]:

            mut_reg_GeneIdSet_dict[mut][reg] = set(
                    mutant_regulations_dictionary[mut][reg].gene_id
                    )

    # For the 3 sets of DEG gene IDs:
    #   - Generate 2 venn's diagrams representing the intersections of
    #     gene IDs. One will be defalut, the other will be unweight.
    #   - Generete an upset plot for also representing the intersections
    #   - From the intersections represented, generate  a dataframe for each
    #     containing all of the relevant information and save it to a file.
    dgeapy.mk_venn_upset_and_intersections_dfs(
            mutant1_gene_set=mut_reg_GeneIdSet_dict[MUTANTS[0]]["DEG"],
            mutant1_name=MUTANTS[0],
            mutant2_gene_set=mut_reg_GeneIdSet_dict[MUTANTS[1]]["DEG"],
            mutant2_name=MUTANTS[1],
            mutant3_gene_set=mut_reg_GeneIdSet_dict[MUTANTS[2]]["DEG"],
            mutant3_name=MUTANTS[2],
            plot_formats=PLOT_FORMATS,
            plots_title="Differentially expressed genes.",
            venn_path=f"{output_venn_dir}/venn_DEG",
            upset_path=f"{output_upset_dir}/upset_DEG",
            df_path=output_df_dir,
            sub_dfs=mutant_regulations_dictionary,
            df_filenames="DEG_interesction",
            )

    # Both up/down_regulation_labels are dictionaries conaining
    # the labels for the next venn diagrams we're going to generate.
    # They'll display the actual number of genes considered
    # up and down regulated at the same time.
    up_regulation_labels = dgeapy.get_gene_ids_set_for_intersections(
            set1=mut_reg_GeneIdSet_dict[MUTANTS[0]]["Up"],
            set2=mut_reg_GeneIdSet_dict[MUTANTS[1]]["Up"],
            set3=mut_reg_GeneIdSet_dict[MUTANTS[2]]["Up"],
            )
    down_regulation_labels = dgeapy.get_gene_ids_set_for_intersections(
            set1=mut_reg_GeneIdSet_dict[MUTANTS[0]]["Down"],
            set2=mut_reg_GeneIdSet_dict[MUTANTS[1]]["Down"],
            set3=mut_reg_GeneIdSet_dict[MUTANTS[2]]["Down"],
            )
    # Generate the same two diagrams but with the labels
    dgeapy.generate_venn3_diagram_with_regulation_labels(
            mutant1_gene_set=mut_reg_GeneIdSet_dict[MUTANTS[0]]["DEG"],
            mutant1_name=MUTANTS[0],
            mutant2_gene_set=mut_reg_GeneIdSet_dict[MUTANTS[1]]["DEG"],
            mutant2_name=MUTANTS[1],
            mutant3_gene_set=mut_reg_GeneIdSet_dict[MUTANTS[2]]["DEG"],
            mutant3_name=MUTANTS[2],
            plot_formats=PLOT_FORMATS,
            up_regulation_labels=up_regulation_labels,
            down_regulation_labels=down_regulation_labels,
            title="Differentially expressed genes.",
            file_path=f"{output_venn_dir}/venn_DEG_labels",
            )

    # Upregulated sets:
    dgeapy.mk_venn_upset_and_intersections_dfs(
            mutant1_gene_set=mut_reg_GeneIdSet_dict[MUTANTS[0]]["Up"],
            mutant1_name=MUTANTS[0],
            mutant2_gene_set=mut_reg_GeneIdSet_dict[MUTANTS[1]]["Up"],
            mutant2_name=MUTANTS[1],
            mutant3_gene_set=mut_reg_GeneIdSet_dict[MUTANTS[2]]["Up"],
            mutant3_name=MUTANTS[2],
            plot_formats=PLOT_FORMATS,
            plots_title="Upregulated genes.",
            venn_path=f"{output_venn_dir}/venn_UP",
            upset_path=f"{output_upset_dir}/upset_UP",
            df_path=output_df_dir,
            sub_dfs=mutant_regulations_dictionary,
            df_filenames="UP_interesction",
            )
    # Downregulated sets:
    dgeapy.mk_venn_upset_and_intersections_dfs(
            mutant1_gene_set=mut_reg_GeneIdSet_dict[MUTANTS[0]]["Down"],
            mutant1_name=MUTANTS[0],
            mutant2_gene_set=mut_reg_GeneIdSet_dict[MUTANTS[1]]["Down"],
            mutant2_name=MUTANTS[1],
            mutant3_gene_set=mut_reg_GeneIdSet_dict[MUTANTS[2]]["Down"],
            mutant3_name=MUTANTS[2],
            plot_formats=PLOT_FORMATS,
            plots_title="Downregulated genes.",
            venn_path=f"{output_venn_dir}/venn_DOWN",
            upset_path=f"{output_upset_dir}/upset_DOWN",
            df_path=output_df_dir,
            sub_dfs=mutant_regulations_dictionary,
            df_filenames="DOWN_interesction",
            )

    # Comparing sets for the 6 possible inverted regulations combinations.
    # Makes a venn diagramm and generates 7 dataframes for the genes of each
    # intersection.
    inverted_reg_dir = f"{output_venn_dir}/inverted_regulations"
    os.mkdir(inverted_reg_dir)
    dgeapy.get_inverted_regulations_and_mk_venns_and_dataframes(
            sets_dictionary=mut_reg_GeneIdSet_dict,
            mutants=MUTANTS,
            sub_dfs=mutant_regulations_dictionary,
            plot_formats=PLOT_FORMATS,
            venn_directory_path=inverted_reg_dir,
            dataframes_directory_path=output_df_dir,
            )

if __name__ == "__main__":
    main()

