#!/usr/bin/env python3

"""
Dataframe handeling, filtering and sub-dataframe generation.
"""


import os
import pandas as pd
from dge_analysis.generate_venn_diagrams import generate_venn3_diagram


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

def get_gene_ids_list_for_intersections(
        set1,
        set2,
        set3,
        ):
    """
    Computes lists contaning gene_ids. Generates one list for each intersection.
    Returns a dictionary where the k eys are binary numbers that represent the
    list for each intersection:
        - First bit equals to the first set.
        - Second bit equals to the second set.
        - Third bit equals to the third set.
    """

    intersections_dict = {}
    intersections_dict["100"] = [i for i in set1 if i not in set2 and i not in set3]
    intersections_dict["110"] = [i for i in set1 if i in set2 and i not in set3]
    intersections_dict["101"] = [i for i in set1 if i not in set2 and i in set3]
    intersections_dict["111"] = [i for i in set1 if i in set2 and i in set3]
    intersections_dict["010"] = [i for i in set2 if i not in set1 and i not in set3]
    intersections_dict["011"] = [i for i in set2 if i not in set1 and i in set3]
    intersections_dict["001"] = [i for i in set3 if i not in set1 and i not in set2]

    return intersections_dict


def sort_df(df):
    """
    """

    column_names = df.columns.values.tolist()

    gene_id = []
    log2FC = []
    fc = []
    regulation = []
    pvalue = []
    padj = []
    gene_info = []
    # There's probably a better way of doing this
    for i in column_names:
        if "gene_id" in i:
            gene_id.append(i)
        elif "_log2FoldChange" in i:
            log2FC.append(i)
        elif "_FoldChange" in i:
            fc.append(i)
        elif "_Regulation" in i:
            regulation.append(i)
        elif "_pvalue" in i:
            pvalue.append(i)
        elif "_padj"in i:
            padj.append(i)
        elif "_chr" in i:
            gene_info.append(i)
        elif "_start"in i:
            gene_info.append(i)
        elif "_end" in i:
            gene_info.append(i)
        elif "_strand" in i:
            gene_info.append(i)
        elif "_length" in i:
            gene_info.append(i)
        elif "_biotype" in i:
            gene_info.append(i)
        elif "_description" in i:
            gene_info.append(i)
        elif "tf_family" in i:
            gene_info.append(i)

    new_colums_names = (
            gene_id
            +
            log2FC
            +
            fc
            +
            regulation
            +
            pvalue
            +
            padj
            +
            gene_info[:8]
            )

    df = df[new_colums_names]
    try:
        df = df.sort_values(
                fc[0],
                ascending=False,
                )
    except:
        pass

    return df


def mk_df_for_each_intersection(
        mutant1,
        set1,
        mutant2,
        set2,
        mutant3,
        set3,
        filtered_dataframes,
        dataframes_directory_path,
        name,
        ):

    new_dir = f"{dataframes_directory_path}/{name}"
    os.mkdir(new_dir)

    intersections_dict = get_gene_ids_list_for_intersections(
            set1=set1,
            set2=set2,
            set3=set3,
            )

    mutant1_df = filtered_dataframes[mutant1][[
        "gene_id",
        f"{mutant1}vsWT_log2FoldChange",
        f"{mutant1}vsWT_FoldChange",
        f"{mutant1}vsWT_Regulation",
        f"{mutant1}vsWT_pvalue",
        f"{mutant1}vsWT_padj",
        "gene_chr",
        "gene_start",
        "gene_end",
        "gene_strand",
        "gene_length",
        "gene_biotype",
        "gene_description",
        "tf_family",
        ]]
    mutant2_df = filtered_dataframes[mutant2][[
        "gene_id",
        f"{mutant2}vsWT_log2FoldChange",
        f"{mutant2}vsWT_FoldChange",
        f"{mutant2}vsWT_Regulation",
        f"{mutant2}vsWT_pvalue",
        f"{mutant2}vsWT_padj",
        "gene_chr",
        "gene_start",
        "gene_end",
        "gene_strand",
        "gene_length",
        "gene_biotype",
        "gene_description",
        "tf_family",
        ]]
    mutant3_df = filtered_dataframes[mutant3][[
        "gene_id",
        f"{mutant3}vsWT_log2FoldChange",
        f"{mutant3}vsWT_FoldChange",
        f"{mutant3}vsWT_Regulation",
        f"{mutant3}vsWT_pvalue",
        f"{mutant3}vsWT_padj",
        "gene_chr",
        "gene_start",
        "gene_end",
        "gene_strand",
        "gene_length",
        "gene_biotype",
        "gene_description",
        "tf_family",
        ]]

    for key in intersections_dict:

        df = pd.DataFrame(data={"gene_id" : intersections_dict[key]})
        df.set_index("gene_id")

        if key[0] == str(1):
            mutant1_df.set_index("gene_id")
            df = pd.merge(
                    df,
                    mutant1_df,
                    on="gene_id",
                    )

        if key[1] == str(1):
            mutant2_df.set_index("gene_id")
            df= pd.merge(
                    df,
                    mutant2_df,
                    on="gene_id",
                    )

        if key[2] == str(1):
            mutant3_df.set_index("gene_id")
            df = pd.merge(
                    df,
                    mutant3_df,
                    on="gene_id",
                    )

        df = sort_df(df)

        df.to_csv(
                f"{new_dir}/{name}_{key}.csv",
                index=False,
                )

        df.to_excel(
                f"{new_dir}/{name}_{key}.xlsx",
                index=False
                )


def get_inverted_regulations_and_mk_venns_and_dataframes(
        sets_dictionary,
        mutants,
        plot_formats,
        venn_directory_path,
        filtered_dataframes,
        dataframes_directory_path,
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
                    name = f"{mutant1}_up_others_down"

                else:
                    mutant1_label = (r"$\downarrow$" + str(mutant1))
                    mutant2_regulation = sets_dictionary[mutant2]["Up"]
                    mutant2_label = (r"$\uparrow$" + str(mutant2))
                    mutant3_regulation = sets_dictionary[mutant3]["Up"]
                    mutant3_label = (r"$\uparrow$" + str(mutant3))
                    name = f"{mutant1}_down_others_up"

                # I'm doing this because i want the mutants to have the same
                # spot and color in the different venn diagrams
                muts = [mutant1, mutant2, mutant3]
                regulations = [mutant1_regulation, mutant2_regulation, mutant3_regulation]
                labels = [mutant1_label, mutant2_label, mutant3_label]
                if i == 0:
                    generate_venn3_diagram(
                            set1=(regulations[0], labels[0]),
                            set2=(regulations[1], labels[1]),
                            set3=(regulations[2], labels[2]),
                            plot_formats=plot_formats,
                            title="Differentially expressed genes.",
                            file_name=f"{venn_directory_path}/venn_{name}",
                            )
                    mk_df_for_each_intersection(
                            mutant1=muts[0],
                            set1=regulations[0],
                            mutant2=muts[1],
                            set2=regulations[1],
                            mutant3=muts[2],
                            set3=regulations[2],
                            filtered_dataframes=filtered_dataframes,
                            dataframes_directory_path=dataframes_directory_path,
                            name=name,
                            )

                elif i == 1:
                    generate_venn3_diagram(
                            set1=(regulations[1], labels[1]),
                            set2=(regulations[0], labels[0]),
                            set3=(regulations[2], labels[2]),
                            plot_formats=plot_formats,
                            title="Differentially expressed genes.",
                            file_name=f"{venn_directory_path}/venn_{name}",
                            )
                    mk_df_for_each_intersection(
                            mutant1=muts[1],
                            set1=regulations[1],
                            mutant2=muts[0],
                            set2=regulations[0],
                            mutant3=muts[2],
                            set3=regulations[2],
                            filtered_dataframes=filtered_dataframes,
                            dataframes_directory_path=dataframes_directory_path,
                            name=name,
                            )


                elif i == 2:
                    generate_venn3_diagram(
                            set1=(regulations[1], labels[1]),
                            set2=(regulations[2], labels[2]),
                            set3=(regulations[0], labels[0]),
                            plot_formats=plot_formats,
                            title="Differentially expressed genes.",
                            file_name=f"{venn_directory_path}/venn_{name}",
                            )
                    mk_df_for_each_intersection(
                            mutant1=muts[1],
                            set1=regulations[1],
                            mutant2=muts[2],
                            set2=regulations[2],
                            mutant3=muts[0],
                            set3=regulations[0],
                            filtered_dataframes=filtered_dataframes,
                            dataframes_directory_path=dataframes_directory_path,
                            name=name,
                            )





