#!/usr/bin/env python3

import json
import pandas as pd
import dge_analysis

#-------# Script Configuration #-----------------------------------------------#

with open("config.json", 'r') as j:
    config_dict = json.loads(j.read())

# The dataframe file name.
DATAFRAME_FILE_NAME = config_dict["dataframe_file_name"]

# The threshold for considering a gene up or down regulated.
FOLD_CHANGE_THRESHOLD = config_dict["fold_change_threshold"]

# The threshold for considering a p-value.
P_VALUE_THRESHOLD = config_dict["p_value_threshold"]

# The threshold for considering a padj.
PADJ_THRESHOLD = config_dict["padj_threshold"]

# The name of the wild type sample.
WILD_TYPE_SAMPLE = config_dict["wild_type_sample"]

# The name of the mutant samples.
MUTANT_SAMPLES = config_dict["mutant_samples"]

# Formats that the generated plots will have:
PLOT_FORMATS = config_dict["plot_formats"]

# mytypes is a dictionary that we're going to pass to the read_excel function.
# Tells pandas explicitly what type of data each column contains.
mytypes = config_dict["datatypes_dictionary"]

df_to_merge_mutant = config_dict["df_to_merge"]["mutant_name"]
df_to_merge_file_name = config_dict["df_to_merge"]["file_name"]

#-------# Main Function #------------------------------------------------------#

# Don't use the recover function if it has been used previously.
try:
    df_recovered = DATAFRAME_FILE_NAME.replace(".xls", "_RECOVERED.xls")

    df = pd.read_excel(
            # The name of the file to open
            df_recovered,
            # The engine to use with ".xls" files according to the documentation.
            engine="xlrd",
            # What type of data each column contains.
            dtype=mytypes,
            # Take "--" as missing values.
            na_values=["--"],
            )
except:
    # Recovering the intial dataframe, as it seems that the format has some
    # type of error.
    df_recovered = dge_analysis.recover_corrupted_file(filename=DATAFRAME_FILE_NAME)

    df = pd.read_excel(
            # The name of the file to open
            df_recovered,
            # The engine to use with ".xls" files according to the documentation.
            engine="xlrd",
            # What type of data each column contains.
            dtype=mytypes,
            # Take "--" as missing values.
            na_values=["--"],
            )

# Add a "Fold Change" column with values.
dge_analysis.add_fold_change_columns(dataframe=df)
# Add "Regulation" column with values.
dge_analysis.add_regulation_columns(dataframe=df)

# Generate multiple sub-dataframes. One for each mutant.
# The sub-dataframes variable will be a dicctionary where the keys corresponf to
# each mutant from MUTANT_SAMPLES and the correspinding value will be the
# generated df.
sub_dataframes = dge_analysis.generate_sub_dataframes(
        dataframe=df,
        mutants_list=MUTANT_SAMPLES
        )

# Merging another df
df_to_merge = pd.read_excel(
        df_to_merge_file_name,
        engine="xlrd",
        na_values=['--']
        )
dge_analysis.add_fold_change_columns(df_to_merge)
dge_analysis.add_regulation_columns(df_to_merge)
sub_dataframes[df_to_merge_mutant] = df_to_merge


dfs_directory = dge_analysis.mk_new_dir("dataframes")
volcano_directory = dge_analysis.mk_new_dir("volcano_plots")

# Saving three sets of gene_ids for each mutant in a dictionary:
#   - DEG (Up or Down)
#   - Up
#   - Down
mutant_regulation_sets_directory = {}

for mutant in MUTANT_SAMPLES:
    # Save to a file.
    sub_dataframes[mutant].to_csv(
            f"{dfs_directory}/{mutant}.csv",
            index=False
            )

    sub_dataframes[mutant].to_excel(
            f"{dfs_directory}/{mutant}.xlsx",
            index=False,
            )

    # Get the name of the columns that'll use for filtering the sub-dataframes.
    column_names_to_check = dge_analysis.column_names_to_check(
            dataframe=sub_dataframes[mutant],
            )

    # Filter the generate sub-dataframe according to the preestablished
    # thresholds for Fold Change, p-value and p-adjust values. Returns the
    # filtered dataframe.
    sub_dataframe_filtered = dge_analysis.filter_FC_PVALUE_PADJ(
            dataframe=sub_dataframes[mutant],
            foldchange_threshold=FOLD_CHANGE_THRESHOLD,
            p_value_threshold=P_VALUE_THRESHOLD,
            padj_threshold=PADJ_THRESHOLD,
            column_names_to_check=column_names_to_check,
            )

    # Saving the sets of gene_ids
    mutant_regulation_sets_directory[mutant] = {
            "DEG" : set(sub_dataframe_filtered.gene_id),
            "Up" : set(
                sub_dataframe_filtered[
                    sub_dataframe_filtered[f"{mutant}vsWT_Regulation"] == "Up"
                                       ]["gene_id"]
                ),
            "Down" : set(
                sub_dataframe_filtered[
                    sub_dataframe_filtered[f"{mutant}vsWT_Regulation"] == "Down"
                    ]["gene_id"]
                ),
            }

    # Save a ".csv" file for the filtered sub-dataframe
    sub_dataframe_filtered.to_csv(
            f"{dfs_directory}/{mutant}_filtered.csv",
            index=False
            )
    sub_dataframe_filtered.to_excel(
            f"{dfs_directory}/{mutant}_filtered.xlsx",
            index=False,
            )

    # Generate a volcano and a count plots
    dge_analysis.generate_volcano_plot(
            dataframe=sub_dataframes[mutant],
            file_path=f"{volcano_directory}/{mutant}",
            x_axis_values=column_names_to_check["log2FoldChange"],
            y_axis_values=column_names_to_check["padj"],
            foldchange_threshold=FOLD_CHANGE_THRESHOLD,
            padj_threshold=PADJ_THRESHOLD,
            mutant_name=mutant,
            plot_formats=PLOT_FORMATS,
            )

# WIP: may be delated in the future:
# 
#     # Don't use the recover function if it has been already used.
#     try:
#         go_df_recovered = pd.read_excel(
#                 # The name of the file to open
#                 f"{mutant}vsWT_ALL_GOenrich_RECOVERED.xls",
#                 # The engine to use with ".xls" files according to the documentation.
#                 engine="xlrd",
#                 )
#     except:
#         go_df_to_recover = dge_analysis.recover_corrupted_file(
#                 f"{mutant}vsWT_ALL_GOenrich.xls"
#                 )
#         go_df_recovered = pd.read_excel(
#                 # The name of the file to open
#                 go_df_to_recover,
#                 # The engine to use with ".xls" files according to the documentation.
#                 engine="xlrd",
#                 )
#     # Add Gene Ontology annotations for each gene.
#     # Multiple annotations on the same gene will be stacked
#     # and separeded with "/".
#     sub_dataframe_with_GO_info = dge_analysis.add_go_columns(
#             df=sub_dataframes[mutant],
#             go_df=go_df_recovered,
#             )
# 
#     sub_dataframe_with_GO_info.to_csv(
#             f"{mutant_file_name}_with_GO_info.csv",
#             index=False
#             )
# 
#     sub_dataframe_filtered_with_GO_info = dge_analysis.filter_FC_PVALUE_PADJ(
#             dataframe=sub_dataframe_with_GO_info,
#             foldchange_threshold=FOLD_CHANGE_THRESHOLD,
#             p_value_threshold=P_VALUE_THRESHOLD,
#             padj_threshold=PADJ_THRESHOLD,
#             column_names_to_check=column_names_to_check,
#             )
# 
#     sub_dataframe_filtered_with_GO_info.to_csv(
#             f"{mutant_file_name}_with_GO_info_filtered.csv",
#             index=False,
#             )
# 
#     dge_analysis.generate_cooler_volcano_plot(
#             dataframe=sub_dataframes[mutant],
#             file_path=mutant_file_name,
#             x_axis_values=column_names_to_check["log2FoldChange"],
#             y_axis_values=column_names_to_check["padj"],
#             foldchange_threshold=FOLD_CHANGE_THRESHOLD,
#             padj_threshold=PADJ_THRESHOLD,
#             picked_go_term1="GO:0040011",
#             picked_go_term2="not now plis",
#             mutant_name=mutant,
#             plot_formats=PLOT_FORMATS,
#             )
#---> Loop ends here

venn_directory = dge_analysis.mk_new_dir(
        new_dir_name="venn_diagrams"
        )

# Generate 2 venn's diagrams representing all of the differentially expressed
# genes. One will be defalut, the other will be unweight.
dge_analysis.generate_venn3_diagram(
        set1=(mutant_regulation_sets_directory[MUTANT_SAMPLES[0]]["DEG"], MUTANT_SAMPLES[0]),
        set2=(mutant_regulation_sets_directory[MUTANT_SAMPLES[1]]["DEG"], MUTANT_SAMPLES[1]),
        set3=(mutant_regulation_sets_directory[MUTANT_SAMPLES[2]]["DEG"], MUTANT_SAMPLES[2]),
        plot_formats=PLOT_FORMATS,
        title="Differentially expressed genes.",
        file_name=f"{venn_directory}/venn_total_regulation",
        )

# Both up/down_regulation_labels are dictionaries conaining the labels for the
# next venn diagrams we're going to generate. They'll display the actual number
# of genes considered up and down regulated at the same time.
up_regulation_labels = dge_analysis.get_labels_for_venn3_diagram(
        set1=mutant_regulation_sets_directory[MUTANT_SAMPLES[0]]["Up"],
        set2=mutant_regulation_sets_directory[MUTANT_SAMPLES[1]]["Up"],
        set3=mutant_regulation_sets_directory[MUTANT_SAMPLES[2]]["Up"],
        )
down_regulation_labels = dge_analysis.get_labels_for_venn3_diagram(
        set1=mutant_regulation_sets_directory[MUTANT_SAMPLES[0]]["Down"],
        set2=mutant_regulation_sets_directory[MUTANT_SAMPLES[1]]["Down"],
        set3=mutant_regulation_sets_directory[MUTANT_SAMPLES[2]]["Down"],
        )

# Generate the same two diagrams but with the labels
dge_analysis.generate_venn3_diagram_with_regulation_labels(
        set1=(mutant_regulation_sets_directory[MUTANT_SAMPLES[0]]["DEG"], MUTANT_SAMPLES[0]),
        set2=(mutant_regulation_sets_directory[MUTANT_SAMPLES[1]]["DEG"], MUTANT_SAMPLES[1]),
        set3=(mutant_regulation_sets_directory[MUTANT_SAMPLES[2]]["DEG"], MUTANT_SAMPLES[2]),
        up_regulation_labels=up_regulation_labels,
        down_regulation_labels=down_regulation_labels,
        plot_formats=PLOT_FORMATS,
        title="Differentially expressed genes.",
        file_name=f"{venn_directory}/venn_total_regulation_w_labels",
        )

# Generate two more diagrams but only for the up-regulated genes
dge_analysis.generate_venn3_diagram(
        set1=(mutant_regulation_sets_directory[MUTANT_SAMPLES[0]]["Up"], MUTANT_SAMPLES[0]),
        set2=(mutant_regulation_sets_directory[MUTANT_SAMPLES[1]]["Up"], MUTANT_SAMPLES[1]),
        set3=(mutant_regulation_sets_directory[MUTANT_SAMPLES[2]]["Up"], MUTANT_SAMPLES[2]),
        plot_formats=PLOT_FORMATS,
        title="Up-regulated genes.",
        file_name=f"{venn_directory}/venn_up_regulation",
        )

# Generate two more diagrams but only for the down-regulated genes
dge_analysis.generate_venn3_diagram(
        set1=(mutant_regulation_sets_directory[MUTANT_SAMPLES[0]]["Down"], MUTANT_SAMPLES[0]),
        set2=(mutant_regulation_sets_directory[MUTANT_SAMPLES[1]]["Down"], MUTANT_SAMPLES[1]),
        set3=(mutant_regulation_sets_directory[MUTANT_SAMPLES[2]]["Down"], MUTANT_SAMPLES[2]),
        plot_formats=PLOT_FORMATS,
        title="Down-regulated genes.",
        file_name=f"{venn_directory}/venn_down_regulation",
        )

# Comparing sets for inverted regulations between mutants...

# Mutant 1 up, others down.
dge_analysis.generate_venn3_diagram(
        set1=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[0]]["Up"],
            (r"$\uparrow$" + MUTANT_SAMPLES[0]),
            ),
        set2=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[1]]["Down"],
            (r"$\downarrow$" + MUTANT_SAMPLES[1]),
            ),
        set3=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[2]]["Down"],
            (r"$\downarrow$" + MUTANT_SAMPLES[2]),
            ),
        plot_formats=PLOT_FORMATS,
        title="Differentially expressed genes.",
        file_name=f"{venn_directory}/venn_{MUTANT_SAMPLES[0]}_up_others_down",
        )

# Mutants 1 and 2 up, 3 down
dge_analysis.generate_venn3_diagram(
        set1=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[0]]["Up"],
            (r"$\uparrow$" + MUTANT_SAMPLES[0]),
            ),
        set2=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[1]]["Up"],
            (r"$\uparrow$" + MUTANT_SAMPLES[1]),
            ),
        set3=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[2]]["Down"],
            (r"$\downarrow$" + MUTANT_SAMPLES[2]),
            ),
        plot_formats=PLOT_FORMATS,
        title="Differentially expressed genes.",
        file_name=f"{venn_directory}/venn_{MUTANT_SAMPLES[2]}_down_others_up",
        )

# Mutants 1 and 3 up, 2 down
dge_analysis.generate_venn3_diagram(
        set1=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[0]]["Up"],
            (r"$\uparrow$" + MUTANT_SAMPLES[0]),
            ),
        set2=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[1]]["Down"],
            (r"$\downarrow$" + MUTANT_SAMPLES[1]),
            ),
        set3=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[2]]["Up"],
            (r"$\uparrow$" + MUTANT_SAMPLES[2]),
            ),
        plot_formats=PLOT_FORMATS,
        title="Differentially expressed genes.",
        file_name=f"{venn_directory}/venn_{MUTANT_SAMPLES[1]}_down_others_up",
        )

# Mutant 1 down, others up.
dge_analysis.generate_venn3_diagram(
        set1=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[0]]["Down"],
            (r"$\downarrow$" + MUTANT_SAMPLES[0]),
            ),
        set2=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[1]]["Up"],
            (r"$\uparrow$" + MUTANT_SAMPLES[1]),
            ),
        set3=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[2]]["Up"],
            (r"$\uparrow$" + MUTANT_SAMPLES[2]),
            ),
        plot_formats=PLOT_FORMATS,
        title="Differentially expressed genes.",
        file_name=f"{venn_directory}/venn_{MUTANT_SAMPLES[0]}_down_others_up",
        )

# Mutants 1 and 2 down, 3 up
dge_analysis.generate_venn3_diagram(
        set1=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[0]]["Down"],
            (r"$\downarrow$" + MUTANT_SAMPLES[0]),
            ),
        set2=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[1]]["Down"],
            (r"$\downarrow$" + MUTANT_SAMPLES[1]),
            ),
        set3=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[2]]["Up"],
            (r"$\uparrow$" + MUTANT_SAMPLES[2]),
            ),
        plot_formats=PLOT_FORMATS,
        title="Differentially expressed genes.",
        file_name=f"{venn_directory}/venn_{MUTANT_SAMPLES[2]}_up_others_down",
        )

# Mutants 1 and 3 down, 2 up
dge_analysis.generate_venn3_diagram(
        set1=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[0]]["Down"],
            (r"$\downarrow$" + MUTANT_SAMPLES[0]),
            ),
        set2=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[1]]["Up"],
            (r"$\uparrow$" + MUTANT_SAMPLES[1]),
            ),
        set3=(
            mutant_regulation_sets_directory[MUTANT_SAMPLES[2]]["Down"],
            (r"$\downarrow$" + MUTANT_SAMPLES[2]),
            ),
        plot_formats=PLOT_FORMATS,
        title="Differentially expressed genes.",
        file_name=f"{venn_directory}/venn_{MUTANT_SAMPLES[1]}_up_others_down",
        )


