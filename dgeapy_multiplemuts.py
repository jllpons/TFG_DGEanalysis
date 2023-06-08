#!/usr/bin/env python3

import os
import sys
import argparse
from dataclasses import dataclass

import pandas as pd

import dgeapy


@dataclass
class SampleData:
    """Store data related to each sample."""
    name: str
    input_df: pd.DataFrame
    df_columns:  dict
    dge_df: pd.DataFrame
    up_df: pd.DataFrame
    down_df: pd.DataFrame


def main():

    description = """
    Differential Gene Expression data analysis between 2, 3 and 4
    'mutant vs. wild' like type experiments."""

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
            '--padj',
            metavar="FLOAT",
            default=0.05,
            type=float,
            help="adjusted p-value threshold, default is 0.05",
            )
    parser.add_argument(
            '--fc',
            metavar="FLOAT",
            default=1.50,
            type=float,
            help="fold change threshold, default is 1.50",
            )
    parser.add_argument(
            "--formats",
            metavar="STR,",
            nargs="?",
            default=["png"],
            type=str,
            action="append",
            help="plot formats, defalut is png",
            )
    parser.add_argument(
            '-n', '--non-coding',
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
    PADJ_THRESHOLD = args.padj
    PLOT_FORMATS = args.formats
    DATAFRAMES = config_dictionary
    include_novels = args.non_coding

    if len(DATAFRAMES) not in [1, 2, 3, 4]:
        sys.exit('\n** dgeapy multiplemuts only supports data from 2, 3 or 4' \
                 ' different samples **\n')

    for k in DATAFRAMES:
        if not os.path.isfile(DATAFRAMES[k]):
                raise FileNotFoundError(f"Could not find file: {DATAFRAMES[k]}")

    cwd = os.getcwd()
    output_dir = f"{cwd}/dgeapy_multiplemuts_output"

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

        if DATAFRAMES[k].endswith(".xlsx"):
            df = pd.read_excel(
                        DATAFRAMES[k],
                        na_values=["--", "",]
                        )
        else:
            df = pd.read_csv(
                        DATAFRAMES[k],
                        sep="\t",
                        na_values=["--", "",]
                        )

        dgeapy.add_fold_change_columns(df)

        dgeapy.add_regulation_columns(df)

        column_names = dgeapy.get_column_names(df)

        # TODO: improve this
        for key in column_names:
            df = df.rename(columns={column_names[key] : key})
            column_names[key] = key

        df = df.set_index(column_names['index'])

        if include_novels is False:
            df = df[~df.index.str.contains("Novel")]
            df = df[~df.index.str.contains("sRNA")]

        # Sometime one old locus tag belongs to two new locus tag so we have
        # to remove one
        df = df.loc[~df.index.duplicated(keep='first')]

        # DEG according to FC and padj value3
        dge_df = dgeapy.filter_FC_PADJ(
                            dataframe=df,
                            foldchange_threshold=FOLD_CHANGE_THRESHOLD,
                            foldchange_column_name='FoldChange',
                            padj_threshold=PADJ_THRESHOLD,
                            padj_column_name='padj',
                            )

        # Up and Down regulated genes from DEG
        up_df = dge_df[dge_df[column_names["Regulation"]] == "Up"]
        donw_df = dge_df[dge_df[column_names["Regulation"]] == "Down"]

        sample_data = SampleData(
                        name=k,
                        input_df=df,
                        df_columns=column_names,
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
        sample_data.down_df.to_excel(
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
    if len(data) == 1:
        sys.exit()

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
    if len(data) == 2:
        up_regulation_labels = dgeapy.get_gene_ids_set_for_intersections2(
                set1=set(data[0].up_df.index),
                set2=set(data[1].up_df.index),
                )
        down_regulation_labels = dgeapy.get_gene_ids_set_for_intersections2(
                set1=set(data[0].down_df.index),
                set2=set(data[1].down_df.index),
                )
        # Generate the same two diagrams but with the labels
        dgeapy.generate_venn2_diagram_with_regulation_labels(
                mutant1_gene_set=set(data[0].dge_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].dge_df.index),
                mutant2_name=data[1].name,
                plot_formats=PLOT_FORMATS,
                up_regulation_labels=up_regulation_labels,
                down_regulation_labels=down_regulation_labels,
                title="Differentially expressed genes",
                file_path=f"{output_dirs_dict['venn']}/venn_DEG_labels",
                )

    elif len(data) == 3:
        up_regulation_labels = dgeapy.get_gene_ids_set_for_intersections3(
                set1=set(data[0].up_df.index),
                set2=set(data[1].up_df.index),
                set3=set(data[2].up_df.index),
                )
        down_regulation_labels = dgeapy.get_gene_ids_set_for_intersections3(
                set1=set(data[0].down_df.index),
                set2=set(data[1].down_df.index),
                set3=set(data[2].down_df.index),
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

    elif len(data) == 4:
        up_regulation_labels = dgeapy.get_gene_ids_set_for_intersections4(
                set1=set(data[0].up_df.index),
                set2=set(data[1].up_df.index),
                set3=set(data[2].up_df.index),
                set4=set(data[3].up_df.index),
                )
        down_regulation_labels = dgeapy.get_gene_ids_set_for_intersections4(
                set1=set(data[0].down_df.index),
                set2=set(data[1].down_df.index),
                set3=set(data[2].down_df.index),
                set4=set(data[3].down_df.index),
                )
        # Generate the same two diagrams but with the labels
        dgeapy.generate_venn4_diagram_with_regulation_labels(
                mutant1_gene_set=set(data[0].dge_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].dge_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].dge_df.index),
                mutant3_name=data[2].name,
                mutant4_gene_set=set(data[3].dge_df.index),
                mutant4_name=data[3].name,
                plot_formats=PLOT_FORMATS,
                up_regulation_labels=up_regulation_labels,
                down_regulation_labels=down_regulation_labels,
                title="Differentially expressed genes",
                file_path=f"{output_dirs_dict['venn']}/venn_DEG_labels",
                )

    # Comparing sets for the possible inverted regulations combinations.
    # TODO: compute inverted regulations for 4-sample data.
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

