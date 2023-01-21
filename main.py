#!/usr/bin/env python3

import dge_analysis
import pandas as pd
import json

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

for mutant in MUTANT_SAMPLES:
    # Create a directory with the name of the mutant
    mutant_directory = dge_analysis.mk_new_dir(
            new_dir_name=f"{mutant}"
            )

    # Template name for the newly crated files.
    mutant_file_name = f"{mutant_directory}/{mutant}"

    # Save a ".csv" file of the dataframe.
    sub_dataframes[mutant].to_csv(f"{mutant_file_name}.csv")

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

    # Save a ".csv" file for the filtered sub-dataframe
    sub_dataframe_filtered.to_csv(f"{mutant_file_name}_filtered.csv")

    # Generate a volcano and a count plots
    dge_analysis.generate_volcano_plot(
            dataframe=sub_dataframes[mutant],
            file_path=mutant_file_name,
            x_axis_values=column_names_to_check["log2FoldChange"],
            y_axis_values=column_names_to_check["padj"],
            foldchange_threshold=FOLD_CHANGE_THRESHOLD,
            padj_threshold=PADJ_THRESHOLD,
            mutant_name=mutant,
            plot_formats=PLOT_FORMATS,
            )

    # Don't use the recover function if it has been already used.
    try:
        go_df_recovered = pd.read_excel(
                # The name of the file to open
                f"{mutant}vsWT_ALL_GOenrich_RECOVERED.xls",
                # The engine to use with ".xls" files according to the documentation.
                engine="xlrd",
                )
    except:
        go_df_to_recover = dge_analysis.recover_corrupted_file(
                f"{mutant}vsWT_ALL_GOenrich.xls"
                )
        go_df_recovered = pd.read_excel(
                # The name of the file to open
                go_df_to_recover,
                # The engine to use with ".xls" files according to the documentation.
                engine="xlrd",
                )

    # Add Gene Ontology annotations for each gene.
    # Multiple annotations on the same gene will be stacked
    # and separeded with "/".
    sub_dataframe_with_GO_info = dge_analysis.add_go_columns(
            df=sub_dataframes[mutant],
            go_df=go_df_recovered,
            )

    sub_dataframe_with_GO_info.to_csv(f"{mutant_file_name}_with_GO_info.csv")


