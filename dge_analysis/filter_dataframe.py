#!/usr/bin/env python3

"""
Dataframe handeling, filtering and sub-dataframe generation.
"""


import pandas as pd
from dge_analysis.generate_venn_diagrams import generate_venn3_diagram
from xlrd.formatting import re


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

    return filtered_dataframe

def get_labels_for_venn3_diagram(
        set1,
        set2,
        set3,
        ):
    """
    Computes labels for a venn diagram between 3 sets. Returns a dictionary
    where the values'll be the numeric values for each intersection of the venn
    diagram. Keys are binary numbers that represent each intersection:
        - First bit equals to the first set.
        - Second bit equals to the second set.
        - Third bit equals to the third set.
    """

    labels_dict = {}
    labels_dict["100"] = len([i for i in set1 if i not in set2 and i not in set3])
    labels_dict["110"] = len([i for i in set1 if i in set2 and i not in set3])
    labels_dict["101"] = len([i for i in set1 if i not in set2 and i in set3])
    labels_dict["111"] = len([i for i in set1 if i in set2 and i in set3])
    labels_dict["010"] = len([i for i in set2 if i not in set1 and i not in set3])
    labels_dict["011"] = len([i for i in set2 if i not in set1 and i in set3])
    labels_dict["001"] = len([i for i in set3 if i not in set1 and i not in set2])

    return labels_dict

def get_inverted_regulations_and_mk_venns_and_dataframes(
        sets_dictionary,
        mutants,
        plot_formats,
        venn_directory_path,
        ):
    """
    Work in progress
    """

    for i in range(len(mutants)):
        mutant1 = mutants[i]
        for regulation in sets_dictionary[mutants[i]]:
            if regulation != "DEG":
                mutant1_regulation = sets_dictionary[mutant1][regulation]
                mutants_2_and_3_list = [i for i in mutants if i is not mutant1]
                mutant2 = mutants_2_and_3_list[0]
                mutant3 = mutants_2_and_3_list[1]

                if regulation == "Up":
                    mutant1_label = (r"$\uparrow$" + str(mutant1))
                    mutant2_regulation = sets_dictionary[mutant2]["Down"]
                    mutant2_label = (r"$\downarrow$" + str(mutant2))
                    mutant3_regulation = sets_dictionary[mutant3]["Down"]
                    mutant3_label = (r"$\downarrow$" + str(mutant3))
                    venn_filename = f"{mutant1}_up_others_down"

                elif regulation == "Down":
                    mutant1_label = (r"$\downarrow$" + str(mutant1))
                    mutant2_regulation = sets_dictionary[mutant2]["Up"]
                    mutant2_label = (r"$\uparrow$" + str(mutant2))
                    mutant3_regulation = sets_dictionary[mutant3]["Up"]
                    mutant3_label = (r"$\uparrow$" + str(mutant3))
                    venn_filename = f"{mutant1}_down_others_up"

                reguations = [mutant1_regulation, mutant2_regulation, mutant3_regulation]
                labels = [mutant1_label, mutant2_label, mutant3_label]

                # I'm doing this because i want the mutants to have the same
                # spot and color in the different venn diagrams
                if i == 0:
                    generate_venn3_diagram(
                            set1=(reguations[0], labels[0]),
                            set2=(reguations[1], labels[1]),
                            set3=(reguations[2], labels[2]),
                            plot_formats=plot_formats,
                            title="Differentially expressed genes.",
                            file_name=f"{venn_directory_path}/venn_{venn_filename}",
                            )
                elif i == 1:
                    generate_venn3_diagram(
                            set1=(reguations[1], labels[1]),
                            set2=(reguations[0], labels[0]),
                            set3=(reguations[2], labels[2]),
                            plot_formats=plot_formats,
                            title="Differentially expressed genes.",
                            file_name=f"{venn_directory_path}/venn_{venn_filename}",
                            )
                elif i == 2:
                    generate_venn3_diagram(
                            set1=(reguations[1], labels[1]),
                            set2=(reguations[2], labels[2]),
                            set3=(reguations[0], labels[0]),
                            plot_formats=plot_formats,
                            title="Differentially expressed genes.",
                            file_name=f"{venn_directory_path}/venn_{venn_filename}",
                            )

