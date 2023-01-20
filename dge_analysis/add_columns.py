#!/usr/bin/env python3

import numpy as np
import pandas as pd


"""
Functions that add columns to a pandas DGE dataframe.
"""


#-------# Function Definitions #-----------------------------------------------#


def add_fold_change_columns(dataframe):
    """
    Adds a "Fold Change" column to the dataframe for each
    "log2 Fold Change" one. Fold Change values are the result of an exponention
    with the number 2 as the base integrer and the corresponding "log2 Fold
    Change" value as the exponent for every "log2 Fold Change" column.
    """

    # We get the column names.
    column_names = dataframe.columns.values.tolist()
    # We set an accumulator so we can keep track of the number of columns 
    # we've been adding.
    accumulator = 1

    # Iterate trough the column names using the index numbers.
    for index_num  in range(len(column_names)):
        # If log2FoldChange" is present in the name of the column:
        if "log2FoldChange" in column_names[index_num]:

            # Saving the column name as a variable.
            log2FoldChange_column_name = str(column_names[index_num])

            # new_column values are the result of exponention of 2 as the base
            # and the corresponding "log2 Fold Change" absolute
            # value as the exponent.
            new_column_values = np.power(
                    2,
                    # The absolute value present in the "log2FoldChange" column.
                    abs(dataframe[log2FoldChange_column_name])
                    )
            # Creating the new column name just by replacing "log2FoldChange".
            # from the previous column with "FoldChange".
            new_column_name = log2FoldChange_column_name.replace(
                "log2FoldChange", 
                "FoldChange",
                )
            # Inserting the new column into the dataframe.
            dataframe.insert(
                    # Location where new_colmn will be inserted.
                    loc=(index_num + accumulator),
                    # Column name
                    column=new_column_name,
                    # Values it will contain
                    value=new_column_values,
                    )
            # Add +1 to the accumulator because we added one more column
            # to the dataframe. By doing it so, the next "FoldChange" column
            # will be placed after the corresponding "log2FoldChange" column.
            accumulator += 1
        else:
            pass


def add_regulation_columns(dataframe):
    """
    Adds a "Regulation" column to the dataframe after each  "Fold Change" one.
    A gene is considered as up-regulated if the "log2 Fold Change" value 
    is positive, and down-regulated if the value is negative. 
    """

    # We get the column names.
    column_names = dataframe.columns.values.tolist()
    # We set an accumulator so we can keep track of the column names we've
    # been adding.
    accumulator = 2

    # Iterate trough the column names using the index numbers.
    for index_num  in range(len(column_names)):
        # If "log2FoldChange" is present in the name of the column:
        if "log2FoldChange" in column_names[index_num]:

            # Saving the column name as a variable.
            log2FoldChange_column_name = str(column_names[index_num])

            # Creating the new column name just by replacing "log2FoldChange"
            # from the previous column with "Regulation".
            new_column_name = log2FoldChange_column_name.replace(
                "log2FoldChange", 
                "Regulation",
                )

            # new_column values are the result of checking if each value form
            # "log2 Fold Change" is positive or negative (greater or less 
            # then zero). We'll fill this list with a for loop.
            new_column_values = []

            # For each value in the "log2FoldChange" column:
            for value in dataframe[log2FoldChange_column_name]:
                # If value is greater than zero, add "Up".
                if value > 0:
                    new_column_values.append("Up")
                # If value is less than zero, add "Down".
                elif value < 0:
                    new_column_values.append("Down")
                # If value is equal to zero, add "Unchanged"
                elif value == 0:
                    new_column_values.append("Unchanged")
                # In any other situation, add "--"
                else:
                    new_column_values.append("--")

            dataframe.insert(
                    # Location where new_colmn will be inserted.
                    loc=(index_num + accumulator),
                    # Columb name
                    column=new_column_name,
                    # Values it will contain
                    value=new_column_values,
                    )

            # Specifiying the new column type as "string".
            dataframe[new_column_name] = dataframe[new_column_name].astype(
                    "string"
                    )

            # Add +1 to the accumulator because we added one more column
            # to the dataframe. By doing it so, the next "Regulation" column
            # will be placed after two positions of the corresponding
            # "log2FoldChange" column.
            accumulator += 1


def create_go_dict(df):
    """
    """
    # TODO

    go_dictionary = {}
    for i in range(len(df)):
        go_category = df["Category"][i]
        go_ID = df["GOID"][i]
        go_descrition = df["Description"][i]

        gene_IDs = df["geneID"][i].split("/")
        for gene in gene_IDs:
            go_dictionary[gene] = {
                    "go_category" : go_category,
                    "go_ID" : go_ID,
                    "go_description" : go_descrition,
                    }

    return go_dictionary



def add_go_columns(
        df,
        go_df,
        ):
    """
    """
    # TODO

    go_dictionary = create_go_dict(go_df)

    df = pd.concat([
            df,
            pd.DataFrame(columns=["GO_category", "GO_ID", "GO_description"])
            ])

    for gene in go_dictionary:
        gene_index = np.where(df["gene_id"] == gene)
        i = int(gene_index[0])
        df.at[i, "GO_category"] = go_dictionary[gene]["go_category"]
        df.at[i, "GO_ID"] = go_dictionary[gene]["go_ID"]
        df.at[i, "GO_description"] = go_dictionary[gene]["go_description"]

    print(df["GO_description"])


