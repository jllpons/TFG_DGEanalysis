#!/usr/bin/env python3

import pandas as pd

"""
Dataframe handeling, filtering and sub-dataframe generation form a more complex
and extensive DGE dataframe.
"""


#-------# Function Definitions #-----------------------------------------------#


def generate_sub_dataframes(
        dataframe,
        mutants_list,
        ):
    """
    This function generates a sub-dataframe for each comparison that's been
    done in the given dataframe. In this case, this is done by selecting the
    columns related to each mutant and inserting all of them in their
    respective sub-dataframe. Returns a dictionary where the name of the mutant
    it's the tag and the corresponding dataframem, the value.
    """

    # Getting the column names as a list.
    column_names = dataframe.columns.values.tolist()

    # This is a dictionary where all the sub-dataframes will be stored.
    sub_dataframes = {}

    # For each mutant (comparison made) in the experiment.
    for mutant in mutants_list:
        # This list contains all mutant names but the one
        # selected by the for loop.
        unwanted_mutant_names = [i for i in mutants_list if i is not mutant]
        # Here we'll be storing all the column names that we want to keep.
        wanted_column_names = []

        # For each name in the given dataframe.
        for name in column_names:
            # This will check if the column name does not contain one of the
            # mutant names that are not related to the mutant name selected
            # in the for loop. If it does not contain them, the name is appended
            # to the wanted_column_names list.
            # NOTE: This is not optimal because I'm assuming that the number of
            #       mutants in all of the given dataframes will be 3.
            if (
                    unwanted_mutant_names[0] not in name
                    and
                    unwanted_mutant_names[1] not in name
                    # Here you could add "unwanted_mutant_names[2] if there
                    # are 4 mutants. Or delate the previous line if there are 2.
                ):
                wanted_column_names.append(name)

        # We have the names but we alse need the values.
        new_dataframe_columns = {}
        # For each wanted column name.
        for column_name in wanted_column_names:
            # Add the values with the column name as tag to the dictionary.
            new_dataframe_columns[column_name] = dataframe[
                    column_name
                    ].values.tolist()

        # Generate the new dataframe
        new_dataframe = pd.DataFrame(
                data=new_dataframe_columns
                )
        # Store the datafrane with the mutant as tag.
        sub_dataframes[mutant] = new_dataframe

    return sub_dataframes


def column_names_to_check(dataframe):
    """
    This function checks if "_FoldChange", "pvalue" or "padj" are present in any
    column name of the given dataframe and returns the names as a dictionary.
    """

    # A list of the column names.
    column_names = dataframe.columns.values.tolist()
    # Here we'll be storing the selected column names.
    column_names_to_check = {}

    # For each name in the column names list.
    for name in column_names:
        # If if contains "_FoldChange" add it to the columns_to_check list.
        if "_FoldChange" in name:
            column_names_to_check["FoldChange"] = name
        # Same for the log2 Fold Change Column
        elif "_log2FoldChange" in name:
            column_names_to_check["log2FoldChange"] = name
        # Same for "pvalue".
        elif "pvalue" in name:
            column_names_to_check["pvalue"] = name
        # Also same for the "padj".
        elif "padj" in name:
            column_names_to_check["padj"] = name

    return column_names_to_check


def filter_FC_PVALUE_PADJ(
        dataframe,
        foldchange_threshold,
        p_value_threshold,
        padj_threshold,
        column_names_to_check
        ):
    """
    This function takes a dataframe and returns it with its values filtered
    according to some given values. It needs a list of columns to check in
    a specific order.
    """

    # Assignin the filtered dataframe to a variable.
    filtered_dataframe = dataframe[ # The original dataframe.
            (dataframe[column_names_to_check["FoldChange"]]
                >=
            foldchange_threshold)
            & # "&" equals "and" operator in pandas.
            (dataframe[column_names_to_check["pvalue"]] < p_value_threshold)
            &
            (dataframe[column_names_to_check["padj"]] < padj_threshold)
            ]
    # NOTE: This is not optimal, i should check if elements of the list have
    #       a correct order.
    return filtered_dataframe
