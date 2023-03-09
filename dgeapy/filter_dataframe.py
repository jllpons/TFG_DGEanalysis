#!/usr/bin/env python3

"""Dataframe handeling, filtering and sub-dataframe generation.
"""


import os

import pandas as pd

from .generate_venn_diagrams import generate_venn3_diagram
from .generate_upset_plots import generate_upset_plot

#-------# Function Definitions #-----------------------------------------------#

def generate_sub_dataframes_3muts(
        dataframe,
        mutants_list,
        ):
    """This function generates a sub-dataframe for each comparison that's been
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
            if (
                    unwanted_mutant_names[0] not in name
                    and
                    unwanted_mutant_names[1] not in name
                ):
                wanted_column_names.append(name)

        # We have the names but we alse need the values.
        new_dataframe_columns = {}
        # For each wanted column name.
        for column_name in wanted_column_names:
            # Add the values with the column name as key to the dictionary.
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
    """This function checks if "_FoldChange", "pvalue" or "padj" are present in any
    column name of the given dataframe and returns the names as a dictionary.
    """

    # A list of the column names.
    column_names = dataframe.columns.values.tolist()
    # Here we'll be storing the selected column names.
    column_names_to_check = {}

    # TODO: add more possibilities to each column name
    # For each name in the column names list.
    for name in column_names:
        # If if contains "_FoldChange" add it to the columns_to_check list.
        if "log2FoldChange" in name:
            column_names_to_check["log2FoldChange"] = name
        # Same for the log2 Fold Change Column
        elif "FoldChange" in name:
            column_names_to_check["FoldChange"] = name
        # Same for "pvalue".
        elif "pvalue" in name:
            column_names_to_check["pvalue"] = name
        # Also same for the "padj".
        elif "padj" in name:
            column_names_to_check["padj"] = name
        elif "Regulation" in name:
            column_names_to_check["Regulation"] = name

    return column_names_to_check


def filter_FC_PVALUE_PADJ(
        dataframe,
        foldchange_threshold,
        foldchange_column_name,
        p_value_threshold,
        p_value_column_name,
        padj_threshold,
        padj_column_name,
        ):
    """This function takes a dataframe and returns it with its values filtered
    according to some given values.
    """

    # Assignin the filtered dataframe to a variable.
    filtered_dataframe = dataframe[ # The original dataframe.
            (dataframe[foldchange_column_name]
                >=
            foldchange_threshold)
            & # "&" equals "and" operator in pandas.
            (dataframe[p_value_column_name] < p_value_threshold)
            &
            (dataframe[padj_column_name] < padj_threshold)
            ]

    return filtered_dataframe

def get_gene_ids_set_for_intersections(
        set1,
        set2,
        set3,
        ):
    """Computes lists contaning gene_ids. Generates one list for each intersection.
    Returns a dictionary where the k eys are binary numbers that represent the
    list for each intersection:
        - First bit equals to the first set.
        - Second bit equals to the second set.
        - Third bit equals to the third set.
    """

    intersections_dict = {
            "100" : set(),
            "010" : set(),
            "001" : set(),
            "110" : set(),
            "101" : set(),
            "011" : set(),
            "111" : set()
            }
    for gene in set1:
        if gene in set2 and gene in set3:
            intersections_dict["111"].add(gene)
        elif gene in set2 and gene not in set3:
            intersections_dict["110"].add(gene)
        elif gene not in set2 and gene in set3:
            intersections_dict["101"].add(gene)
        elif gene not in set2 and gene not in set3:
            intersections_dict["100"].add(gene)

    for gene in set2:
        if gene not in set1 and gene in set3:
            intersections_dict["011"].add(gene)
        elif gene not in set1 and gene not in set3:
            intersections_dict["010"].add(gene)

    for gene in set3:
        if gene not in set1 and gene not in set2:
            intersections_dict["001"].add(gene)

    return intersections_dict


def sort_df(df, mutant_names):
    """Takes a df generated in mk_df_for_each_intersection() and sorts the
    columns. Rows are sorted taking the first foldchange column values.
    """

    column_names = df.columns.values.tolist()

    # Column names can be classified into this 7 lists.
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
        elif "log2FoldChange" in i:
            log2FC.append(i)
        elif "FoldChange" in i:
            fc.append(i)
        elif "Regulation" in i:
            regulation.append(i)
        elif "pvalue" in i:
            pvalue.append(i)
        elif "padj"in i:
            padj.append(i)
        elif "chr" in i:
            gene_info.append(i)
        elif "start"in i:
            gene_info.append(i)
        elif "end" in i:
            gene_info.append(i)
        elif "strand" in i:
            gene_info.append(i)
        elif "length" in i:
            gene_info.append(i)
        elif "biotype" in i:
            gene_info.append(i)
        elif "description" in i:
            gene_info.append(i)
        elif "tf_family" in i:
            gene_info.append(i)

    # We have the column names goruped so we can choos the ordrer. Gene info
    # columns have the same values so we only take the first 8.
    new_colums_names = (
            gene_id
            + log2FC
            + fc
            + regulation
            + pvalue
            + padj
            + gene_info[:8]
            )

    df = df[new_colums_names]
    try:
        df = df.sort_values(
                    fc[0],
                    ascending=False,
                    )

        return df

    except:
        return df


def mk_df_for_each_intersection(
        mutant1_gene_set,
        mutant1_name,
        mutant2_gene_set,
        mutant2_name,
        mutant3_gene_set,
        mutant3_name,
        sub_dfs,
        path,
        file_names,
        ):
    """Takes 3 sets of gene IDs and computes the 7 possible intersections
    Creates a dataframe for each intersection that will display all of the
    relevant information for each gene.
    """

    # Storing output files in here
    os.mkdir(f"{path}/{file_names}")

    # Generate a dictionary contaning all of the gene IDs intersections.
    intersections_dict = get_gene_ids_set_for_intersections(
            set1=mutant1_gene_set,
            set2=mutant2_gene_set,
            set3=mutant3_gene_set,
            )

    # Columns that the generated dataframes will contain.
    columns = [
            "gene_id",
            "log2FoldChange",
            "FoldChange",
            "Regulation",
            "pvalue",
            "padj",
            "gene_chr",
            "gene_start",
            "gene_end",
            "gene_strand",
            "gene_length",
            "gene_biotype",
            "gene_description",
            "tf_family",
            ]
    mutant1_df = sub_dfs[mutant1_name]["DEG"][columns]
    mutant1_df = mutant1_df.rename(columns={
                    "log2FoldChange" : f"{mutant1_name}_log2FoldChange",
                    "FoldChange" : f"{mutant1_name}_FoldChange",
                    "padj" : f"{mutant1_name}_padj",
                    "pvalue" : f"{mutant1_name}_pvalue",
                    "Regulation" : f"{mutant1_name}_Regulation",
                    })
    mutant2_df = sub_dfs[mutant2_name]["DEG"][columns]
    mutant2_df = mutant2_df.rename(columns={
                    "log2FoldChange" : f"{mutant2_name}_log2FoldChange",
                    "FoldChange" : f"{mutant2_name}_FoldCange",
                    "padj" : f"{mutant2_name}_padj",
                    "pvalue" : f"{mutant2_name}_pvalue",
                    "Regulation" : f"{mutant2_name}_Regulation",
                    })
    mutant3_df = sub_dfs[mutant3_name]["DEG"][columns]
    mutant3_df = mutant3_df.rename(columns={
                    "log2FoldChange" : f"{mutant3_name}_log2FoldChange",
                    "FoldChange" : f"{mutant3_name}_FoldCange",
                    "padj" : f"{mutant3_name}_padj",
                    "pvalue" : f"{mutant3_name}_pvalue",
                    "Regulation" : f"{mutant3_name}_Regulation",
                    })

    # For each one of the intersections
    for key in intersections_dict:

        # Making a dataframe containing only the gene IDs presents
        # In that intersection. Also setting the gene IDs as index.
        df = pd.DataFrame(data={"gene_id" : list(intersections_dict[key])})
        df.set_index("gene_id")

        # If first bit is 1, there're gene IDs from the 1st mutant DE genes
        # that are present in that intersection. So we take the
        # 1st mutant DE genes dataframe and we add all of the columns
        # showed above for ONLY the gene IDs that are already present
        # in the previosly created df.
        if key[0] == str(1):
            mutant1_df.set_index("gene_id")
            df = pd.merge(
                    df,
                    mutant1_df,
                    on="gene_id",
                    )

        # The second bit being 1 means there're gene IDs form the 2nd mutant
        # DE genes present in that intersection. We do the same as above.
        if key[1] == str(1):
            mutant2_df.set_index("gene_id")
            df= pd.merge(
                    df,
                    mutant2_df,
                    on="gene_id",
                    )

        # The third bit being 1 means there're gene IDs form the 3rd mutant
        # DE genes present in that intersection. We do the same as before.
        if key[2] == str(1):
            mutant3_df.set_index("gene_id")
            df = pd.merge(
                    df,
                    mutant3_df,
                    on="gene_id",
                    )

        # We sort the columns for a cleaner visualization
        df = sort_df(df, (mutant1_name, mutant2_name, mutant3_name))

        df.to_csv(
                f"{path}/{file_names}/{file_names}_{key}.tsv",
                index=False,
                )

        df.to_excel(
                f"{path}/{file_names}/{file_names}_{key}.xlsx",
                index=False
                )

def mk_venn_upset_and_intersections_dfs(
        mutant1_gene_set,
        mutant1_name,
        mutant2_gene_set,
        mutant2_name,
        mutant3_gene_set,
        mutant3_name,
        plot_formats,
        plots_title,
        venn_path,
        upset_path,
        df_path,
        sub_dfs,
        df_filenames
        ):
    """Takes 3 sets of gene IDs and the respective mutant name and generates
    the correspoding venn diagrams, upset plots and a dataframe for each
    one of the 7 intersections.
    """

    generate_venn3_diagram(
        mutant1_gene_set,
        mutant1_name,
        mutant2_gene_set,
        mutant2_name,
        mutant3_gene_set,
        mutant3_name,
        plot_formats,
        plots_title,
        venn_path,
        )
    generate_upset_plot(
        mutant1_gene_set,
        mutant1_name,
        mutant2_gene_set,
        mutant2_name,
        mutant3_gene_set,
        mutant3_name,
        plot_formats,
        plots_title,
        upset_path,
        )
    mk_df_for_each_intersection(
        mutant1_gene_set,
        mutant1_name,
        mutant2_gene_set,
        mutant2_name,
        mutant3_gene_set,
        mutant3_name,
        sub_dfs,
        df_path,
        df_filenames,
        )


def get_inverted_regulations_and_mk_venns_and_dataframes(
        sets_dictionary,
        mutants,
        sub_dfs,
        plot_formats,
        venn_directory_path,
        dataframes_directory_path,
        ):
    """Generates the venn diagramas and the correspoding intersection
    dataframe for each one of the possible inverted regulations combinations
    """

    inverted_regulation_dict = {}
    regulations = ("Up", "Down")

    for mut in mutants:
        other_mutants = [m for m in mutants if m is not mut]

        for reg in regulations:

            inverted_regulation_dict[mut] = {
                    "reg" : reg,
                    "label" : (r"$\uparrow$" + mut)
                    }

            other_reg = "Down"
            inverted_regulation_dict[other_mutants[0]] = {
                "reg" : other_reg,
                "label" : (r"$\downarrow$" + other_mutants[0])
                }
            inverted_regulation_dict[other_mutants[1]] = {
                "reg" : other_reg,
                "label" : (r"$\downarrow$" + other_mutants[1])
                }

            name = f"{mut}_Up_others_Down"

            if reg == "Down":
                inverted_regulation_dict[mut] = {
                        "reg" : reg,
                        "label" : (r"$\downarrow$" + mut)
                        }

                other_reg = "Up"
                inverted_regulation_dict[other_mutants[0]] = {
                    "reg" : other_reg,
                    "label" : (r"$\uparrow$" + other_mutants[0])
                    }
                inverted_regulation_dict[other_mutants[1]] = {
                    "reg" : other_reg,
                    "label" : (r"$\uparrow$" + other_mutants[1])
                    }

                name = f"{mut}_Down_others_Up"

            generate_venn3_diagram(
                    mutant1_gene_set=sets_dictionary[mutants[0]][
                        inverted_regulation_dict[mutants[0]]["reg"]
                        ],
                    mutant1_name=inverted_regulation_dict[mutants[0]]["label"],
                    mutant2_gene_set=sets_dictionary[mutants[1]][
                        inverted_regulation_dict[mutants[1]]["reg"]
                        ],
                    mutant2_name=inverted_regulation_dict[mutants[1]]["label"],
                    mutant3_gene_set=sets_dictionary[mutants[2]][
                        inverted_regulation_dict[mutants[2]]["reg"]
                        ],
                    mutant3_name=inverted_regulation_dict[mutants[2]]["label"],
                    plot_formats=plot_formats,
                    title="Inverted regulations.",
                    path=f"{venn_directory_path}/{name}",
                    )
            mk_df_for_each_intersection(
                mutant1_gene_set=sets_dictionary[mutants[0]][
                    inverted_regulation_dict[mutants[0]]["reg"]
                    ],
                mutant1_name=mutants[0],
                mutant2_gene_set=sets_dictionary[mutants[1]][
                    inverted_regulation_dict[mutants[1]]["reg"]
                    ],
                mutant2_name=mutants[1],
                mutant3_gene_set=sets_dictionary[mutants[2]][
                    inverted_regulation_dict[mutants[2]]["reg"]
                    ],
                mutant3_name=mutants[2],
                sub_dfs=sub_dfs,
                path=dataframes_directory_path,
                file_names=name
                )

