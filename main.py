#!/usr/bin/env python3

import json
import os
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
            na_values=["--", "-", "-//- && -", "",],
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

sub_dataframes_filtered = {}
# Merging another df
df_to_merge = pd.read_excel(
        df_to_merge_file_name,
        engine="xlrd",
        na_values=["--", "-", "-//- && -", "",]
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
    mutant_dir = f"{dfs_directory}/{mutant}"
    os.mkdir(mutant_dir)

    # Save to a file.
    sub_dataframes[mutant].to_csv(
            f"{mutant_dir}/{mutant}.csv",
            index=False
            )

    sub_dataframes[mutant].to_excel(
            f"{mutant_dir}/{mutant}.xlsx",
            index=False,
            )

    # Get the name of the columns that'll use for filtering the sub-dataframes.
    column_names_to_check = dge_analysis.column_names_to_check(
            dataframe=sub_dataframes[mutant],
            )

    # Filter the generate sub-dataframe according to the preestablished
    # thresholds for Fold Change, p-value and p-adjust values. Returns the
    # filtered dataframe.
    sub_dataframes_filtered[mutant] = dge_analysis.filter_FC_PVALUE_PADJ(
            dataframe=sub_dataframes[mutant],
            foldchange_threshold=FOLD_CHANGE_THRESHOLD,
            p_value_threshold=P_VALUE_THRESHOLD,
            padj_threshold=PADJ_THRESHOLD,
            column_names_to_check=column_names_to_check,
            )

    # subframes conaining only up or down-regulated genes.
    sub_frame_up = sub_dataframes_filtered[mutant][
            sub_dataframes_filtered[mutant][f"{mutant}vsWT_Regulation"] == "Up"
            ]
    sub_frame_down = sub_dataframes_filtered[mutant][
            sub_dataframes_filtered[mutant][f"{mutant}vsWT_Regulation"] == "Down"
            ]

    # Adding the sets of gene_ids to the dictionary.
    mutant_regulation_sets_directory[mutant] = {
            "DEG" : set(sub_dataframes_filtered[mutant].gene_id),
            "Up" : set(sub_frame_up.gene_id),
            "Down" : set(sub_frame_down.gene_id),
            }

    # Save filtered sub-dataframes to files.
    sub_dataframes_filtered[mutant].to_csv(
            f"{mutant_dir}/{mutant}_DEG.csv",
            index=False
            )
    sub_dataframes_filtered[mutant].to_excel(
            f"{mutant_dir}/{mutant}_DEG.xlsx",
            index=False,
            )
    sub_frame_up.to_csv(
            f"{mutant_dir}/{mutant}_Up.csv",
            index=False
            )
    sub_frame_up.to_excel(
            f"{mutant_dir}/{mutant}_Up.xlsx",
            index=False,
            )
    sub_frame_down.to_csv(
            f"{mutant_dir}/{mutant}_Down.csv",
            index=False
            )
    sub_frame_down.to_excel(
            f"{mutant_dir}/{mutant}_Down.xlsx",
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
        file_name=f"{venn_directory}/venn_DEG",
        )

dge_analysis.mk_df_for_each_intersection(
        mutant1=MUTANT_SAMPLES[0],
        set1=mutant_regulation_sets_directory[MUTANT_SAMPLES[0]]["DEG"],
        mutant2=MUTANT_SAMPLES[1],
        set2=mutant_regulation_sets_directory[MUTANT_SAMPLES[1]]["DEG"],
        mutant3=MUTANT_SAMPLES[2],
        set3=mutant_regulation_sets_directory[MUTANT_SAMPLES[2]]["DEG"],
        filtered_dataframes=sub_dataframes_filtered,
        dataframes_directory_path=dfs_directory,
        name="DEG_intersections",
        )



# Both up/down_regulation_labels are dictionaries conaining the labels for the
# next venn diagrams we're going to generate. They'll display the actual number
# of genes considered up and down regulated at the same time.
up_regulation_labels = dge_analysis.get_gene_ids_list_for_intersections(
        set1=mutant_regulation_sets_directory[MUTANT_SAMPLES[0]]["Up"],
        set2=mutant_regulation_sets_directory[MUTANT_SAMPLES[1]]["Up"],
        set3=mutant_regulation_sets_directory[MUTANT_SAMPLES[2]]["Up"],
        )
down_regulation_labels = dge_analysis.get_gene_ids_list_for_intersections(
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
        file_name=f"{venn_directory}/venn_DEG_with_labels",
        )

# Generate two more diagrams but only for the up-regulated genes
dge_analysis.generate_venn3_diagram(
        set1=(mutant_regulation_sets_directory[MUTANT_SAMPLES[0]]["Up"], MUTANT_SAMPLES[0]),
        set2=(mutant_regulation_sets_directory[MUTANT_SAMPLES[1]]["Up"], MUTANT_SAMPLES[1]),
        set3=(mutant_regulation_sets_directory[MUTANT_SAMPLES[2]]["Up"], MUTANT_SAMPLES[2]),
        plot_formats=PLOT_FORMATS,
        title="Up-regulated genes.",
        file_name=f"{venn_directory}/venn_Up",
        )
dge_analysis.mk_df_for_each_intersection(
        mutant1=MUTANT_SAMPLES[0],
        set1=mutant_regulation_sets_directory[MUTANT_SAMPLES[0]]["Up"],
        mutant2=MUTANT_SAMPLES[1],
        set2=mutant_regulation_sets_directory[MUTANT_SAMPLES[1]]["Up"],
        mutant3=MUTANT_SAMPLES[2],
        set3=mutant_regulation_sets_directory[MUTANT_SAMPLES[2]]["Up"],
        filtered_dataframes=sub_dataframes_filtered,
        dataframes_directory_path=dfs_directory,
        name="Up_intersections",
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
dge_analysis.mk_df_for_each_intersection(
        mutant1=MUTANT_SAMPLES[0],
        set1=mutant_regulation_sets_directory[MUTANT_SAMPLES[0]]["Down"],
        mutant2=MUTANT_SAMPLES[1],
        set2=mutant_regulation_sets_directory[MUTANT_SAMPLES[1]]["Down"],
        mutant3=MUTANT_SAMPLES[2],
        set3=mutant_regulation_sets_directory[MUTANT_SAMPLES[2]]["Down"],
        filtered_dataframes=sub_dataframes_filtered,
        dataframes_directory_path=dfs_directory,
        name="Down_intersections",
        )

# Comparing sets for the 6 possible inverted regulations combinations.
# Makes a venn diagramm and generates 7 dataframes for the genes of each
# intersection.
inverted_reg_dir = f"{venn_directory}/inverted_regulations"
os.mkdir(inverted_reg_dir)
dge_analysis.get_inverted_regulations_and_mk_venns_and_dataframes(
        sets_dictionary=mutant_regulation_sets_directory,
        mutants=MUTANT_SAMPLES,
        plot_formats=PLOT_FORMATS,
        venn_directory_path=inverted_reg_dir,
        filtered_dataframes=sub_dataframes_filtered,
        dataframes_directory_path=dfs_directory,
        )

