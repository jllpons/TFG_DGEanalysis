#!/usr/bin/env python3

"""Dataframe handeling, filtering and sub-dataframe generation.
"""

import pdb

import os

import pandas as pd

from .venn_diagrams import generate_venn2_diagram, generate_venn3_diagram, generate_venn4_diagram
from .upset_plots import generate_upset_plot


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


def get_column_names(dataframe):
    """Matches column names to keys inn a dictionary
    """

    # A list of the column names.
    column_names = dataframe.columns.values.tolist()
    # Here we'll be storing the selected column names.
    column_names_to_check = {}

    # TODO: add more possibilities to each column name
    # For each name in the column names list.
    for name in column_names:
        if any([x in name for x in ['gene_id', 'Gene', 'mapped_geneID', 'geneID']]):
            column_names_to_check['geneID'] = name
        elif any([x in name for x in ['log2FoldChange', 'log2foldchange']]):
            column_names_to_check["log2FoldChange"] = name
        elif any([x in name for x in ['FoldChange', 'foldchange']]):
            column_names_to_check["FoldChange"] = name
        elif "pvalue" in name:
            column_names_to_check["pvalue"] = name
        elif "padj" in name:
            column_names_to_check["padj"] = name
        elif "Regulation" in name:
            column_names_to_check["Regulation"] = name
        elif any([x in name for x in ['gene_description', 'hypothetical protein', 'Description']]):
            column_names_to_check['Description'] = name
        else:
            # Lots of changes since the initial design
            column_names_to_check[name] = name

    return column_names_to_check


def filter_FC_PADJ(
        dataframe,
        foldchange_threshold,
        foldchange_column_name,
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
            (dataframe[padj_column_name] < padj_threshold)
            ]

    return filtered_dataframe


def get_gene_ids_set_for_intersections2(
        set1,
        set2,
        ):
    """Computes lists contaning gene_ids. Generates one list for each intersection.
    Returns a dictionary where the keys are binary numbers that represent the
    list for each intersection:
        - First bit equals to the first set.
        - Second bit equals to the second set.
    """

    intersections_dict = {
            "10" : set(),
            "01" : set(),
            "11" : set()
            }
    for gene in set1:
        if gene in set2:
            intersections_dict["11"].add(gene)
        elif gene not in set2:
            intersections_dict["10"].add(gene)

    for gene in set2:
        if gene not in set1:
            intersections_dict["01"].add(gene)

    return intersections_dict


def get_gene_ids_set_for_intersections3(
        set1,
        set2,
        set3,
        ):
    """Computes lists contaning gene_ids. Generates one list for each intersection.
    Returns a dictionary where the keys are binary numbers that represent the
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


def get_gene_ids_set_for_intersections4(set1, set2, set3, set4):
    """
    Computes sets containing gene IDs. Generates one set for each intersection.
    Returns a dictionary where the keys are binary numbers that represent the
    list for each intersection:
        - First bit equals to the first set.
        - Second bit equals to the second set.
        - Third bit equals to the third set.
        - Fourth bit equals to the fourth set.
    """
    intersections_dict = {
        "0001": set(),
        "0010": set(),
        "0100": set(),
        "1000": set(),
        "0011": set(),
        "0101": set(),
        "1001": set(),
        "0110": set(),
        "1010": set(),
        "1100": set(),
        "0111": set(),
        "1011": set(),
        "1101": set(),
        "1110": set(),
        "1111": set(),
    }

    for gene in set1:
        if gene in set2 and gene in set3 and gene in set4:
            intersections_dict["1111"].add(gene)
        elif gene in set2 and gene in set3 and gene not in set4:
            intersections_dict["1110"].add(gene)
        elif gene in set2 and gene not in set3 and gene in set4:
            intersections_dict["1101"].add(gene)
        elif gene in set2 and gene not in set3 and gene not in set4:
            intersections_dict["1100"].add(gene)
        elif gene not in set2 and gene in set3 and gene in set4:
            intersections_dict["1011"].add(gene)
        elif gene not in set2 and gene in set3 and gene not in set4:
            intersections_dict["1010"].add(gene)
        elif gene not in set2 and gene not in set3 and gene in set4:
            intersections_dict["1001"].add(gene)
        else:
            intersections_dict["1000"].add(gene)

    for gene in set2:
        if gene not in set1 and gene in set3 and gene in set4:
            intersections_dict['0111'].add(gene)
        elif gene not in set1 and gene in set3 and gene not in set4:
            intersections_dict['0110'].add(gene)
        elif gene not in set1 and gene not in set3 and gene in set4:
            intersections_dict['0101'].add(gene)
        elif gene not in set1 and gene not in set3 and gene not in set4:
            intersections_dict['0100'].add(gene)

    for gene in set3:
        if gene not in set1 and gene not in set2 and gene in set4:
            intersections_dict['0011'].add(gene)
        elif gene not in set1 and gene not in set2 and gene not in set4:
            intersections_dict['0010'].add(gene)

    for gene in set4:
        if gene not in set1 and gene not in set2 and gene not in set3:
            intersections_dict['0001'].add(gene)

    return intersections_dict


def sort_df(df):
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
    go_annotations = []

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
        elif "Description" in i:
            gene_info.append(i)
        elif "tf_family" in i:
            gene_info.append(i)
        elif "GO" in i:
            go_annotations.append(i)
        else:
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
            + gene_info
            + go_annotations
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


def mk_df_for_each_intersection2(
        mutant1_gene_set,
        mutant1_name,
        mutant2_gene_set,
        mutant2_name,
        data,
        path,
        file_names,
        ):
    """Takes 2 sets of gene IDs and computes the 3 possible intersections
    Creates a dataframe for each intersection that will display all of the
    relevant information for each gene.
    """

    # Storing output files in here
    os.mkdir(f"{path}/{file_names}")

    # Generate a dictionary contaning all of the gene IDs intersections.
    intersections_dict = get_gene_ids_set_for_intersections2(
            set1=mutant1_gene_set,
            set2=mutant2_gene_set,
            )

    # Columns that the generated dataframes will contain.
    columns = [
            'log2FoldChange',
            'FoldChange',
            'padj',
            'pvalue',
            'Regulation',
            'Description',
            ]

    mutant1_df = data[0].dge_df[columns]
    mutant2_df = data[1].dge_df[columns]
    for col in columns:
        mutant1_df = mutant1_df.rename(columns={col : f'{data[0].name}_{col}'})
        mutant2_df = mutant2_df.rename(columns={col : f'{data[1].name}_{col}'})

    # For each one of the intersections
    for key in intersections_dict:

        # Making a dataframe containing only the gene IDs presents
        # In that intersection. Also setting the gene IDs as index.
        df = pd.DataFrame(data={'geneID' : list(intersections_dict[key])})
        df = df.set_index('geneID')

        # If first bit is 1, there're gene IDs from the 1st mutant DE genes
        # that are present in that intersection. So we take the
        # 1st mutant DE genes dataframe and we add all of the columns
        # showed above for ONLY the gene IDs that are already present
        # in the previosly created df.
        if key[0] == str(1):
            df = pd.merge(
                    df,
                    mutant1_df,
                    on='geneID',
                    )

        # The second bit being 1 means there're gene IDs form the 2nd mutant
        # DE genes present in that intersection. We do the same as above.
        if key[1] == str(1):
            df= pd.merge(
                    df,
                    mutant2_df,
                    on='geneID',
                    )


        # We sort the columns for a cleaner visualization
        df = sort_df(df)

        df.to_csv(
                f"{path}/{file_names}/{file_names}_{key}.tsv",
                )

        df.to_excel(
                f"{path}/{file_names}/{file_names}_{key}.xlsx",
                )


def mk_df_for_each_intersection3(
        mutant1_gene_set,
        mutant1_name,
        mutant2_gene_set,
        mutant2_name,
        mutant3_gene_set,
        mutant3_name,
        data,
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
    intersections_dict = get_gene_ids_set_for_intersections3(
            set1=mutant1_gene_set,
            set2=mutant2_gene_set,
            set3=mutant3_gene_set,
            )

    # Columns that the generated dataframes will contain.
    columns1 = list(data[0].df_columns.keys())
    columns1.remove('index')
    columns2 = list(data[1].df_columns.keys())
    columns2.remove('index')
    columns3 = list(data[2].df_columns.keys())
    columns3.remove('index')

    mutant1_df = data[0].dge_df[columns1]
    mutant2_df = data[1].dge_df[columns2]
    mutant3_df = data[2].dge_df[columns3]
    for col in columns1:
        mutant1_df = mutant1_df.rename(columns={col : f'{data[0].name}_{col}'})
    for col in columns2:
        mutant2_df = mutant2_df.rename(columns={col : f'{data[1].name}_{col}'})
    for col in columns3:
        mutant3_df = mutant3_df.rename(columns={col : f'{data[2].name}_{col}'})

    # For each one of the intersections
    for key in intersections_dict:

        # Making a dataframe containing only the gene IDs presents
        # In that intersection. Also setting the gene IDs as index.
        df = pd.DataFrame(data={"index" : list(intersections_dict[key])})
        df = df.set_index("index")

        # If first bit is 1, there're gene IDs from the 1st mutant DE genes
        # that are present in that intersection. So we take the
        # 1st mutant DE genes dataframe and we add all of the columns
        # showed above for ONLY the gene IDs that are already present
        # in the previosly created df.
        if key[0] == str(1):
            df = pd.merge(
                    df,
                    mutant1_df,
                    on="index",
                    )

        # The second bit being 1 means there're gene IDs form the 2nd mutant
        # DE genes present in that intersection. We do the same as above.
        if key[1] == str(1):
            df = pd.merge(
                    df,
                    mutant2_df,
                    on="index",
                    )

        # The third bit being 1 means there're gene IDs form the 3rd mutant
        # DE genes present in that intersection. We do the same as before.
        if key[2] == str(1):
            df = pd.merge(
                    df,
                    mutant3_df,
                    on="index",
                    )

        # We sort the columns for a cleaner visualization
        df = sort_df(df)

        df.to_csv(
                f"{path}/{file_names}/{file_names}_{key}.tsv",
                )

        df.to_excel(
                f"{path}/{file_names}/{file_names}_{key}.xlsx",
                )


def mk_df_for_each_intersection4(
        mutant1_gene_set,
        mutant1_name,
        mutant2_gene_set,
        mutant2_name,
        mutant3_gene_set,
        mutant3_name,
        mutant4_gene_set,
        mutant4_name,
        data,
        path,
        file_names,
        ):
    """Takes 4 sets of gene IDs and computes the 16 possible intersections
    Creates a dataframe for each intersection that will display all of the
    relevant information for each gene.
    """

    # Storing output files in here
    os.mkdir(f"{path}/{file_names}")

    # Generate a dictionary contaning all of the gene IDs intersections.
    intersections_dict = get_gene_ids_set_for_intersections4(
            set1=mutant1_gene_set,
            set2=mutant2_gene_set,
            set3=mutant3_gene_set,
            set4=mutant4_gene_set,
            )

    # Columns that the generated dataframes will contain.
    columns1 = list(data[0].df_columns.keys())
    columns1.remove('index')
    columns2 = list(data[1].df_columns.keys())
    columns2.remove('index')
    columns3 = list(data[2].df_columns.keys())
    columns3.remove('index')
    columns4 = list(data[3].df_columns.keys())
    columns4.remove('index')

    mutant1_df = data[0].dge_df[columns1]
    mutant2_df = data[1].dge_df[columns2]
    mutant3_df = data[2].dge_df[columns3]
    mutant4_df = data[3].dge_df[columns4]
    for col in columns1:
        mutant1_df = mutant1_df.rename(columns={col : f'{data[0].name}_{col}'})
    for col in columns2:
        mutant2_df = mutant2_df.rename(columns={col : f'{data[1].name}_{col}'})
    for col in columns3:
        mutant3_df = mutant3_df.rename(columns={col : f'{data[2].name}_{col}'})
    for col in columns4:
        mutant4_df = mutant4_df.rename(columns={col : f'{data[3].name}_{col}'})

    # For each one of the intersections
    for key in intersections_dict:

        # Making a dataframe containing only the gene IDs presents
        # In that intersection. Also setting the gene IDs as index.
        df = pd.DataFrame(data={"index" : list(intersections_dict[key])})
        df = df.set_index("index")

        # If first bit is 1, there're gene IDs from the 1st mutant DE genes
        # that are present in that intersection. So we take the
        # 1st mutant DE genes dataframe and we add all of the columns
        # showed above for ONLY the gene IDs that are already present
        # in the previosly created df.
        if key[0] == str(1):
            df = pd.merge(
                    df,
                    mutant1_df,
                    on="index",
                    )

        # The second bit being 1 means there're gene IDs form the 2nd mutant
        # DE genes present in that intersection. We do the same as above.
        if key[1] == str(1):
            df= pd.merge(
                    df,
                    mutant2_df,
                    on="index",
                    )

        # The third bit being 1 means there're gene IDs form the 3rd mutant
        # DE genes present in that intersection. We do the same as before.
        if key[2] == str(1):
            df = pd.merge(
                    df,
                    mutant3_df,
                    on="index",
                    )

        # The fourth bit being 1 means there're gene IDs form the 4th mutant
        # DE genes present in that intersection. We do the same as before.
        if key[3] == str(1):
            df = pd.merge(
                    df,
                    mutant4_df,
                    on="index",
                    )

        # We sort the columns for a cleaner visualization
        df = sort_df(df)

        df.to_csv(
                f"{path}/{file_names}/{file_names}_{key}.tsv",
                )

        df.to_excel(
                f"{path}/{file_names}/{file_names}_{key}.xlsx",
                )


def mk_venn_upset_and_intersections_dfs(
        data,
        plot_formats,
        venn_path,
        upset_path,
        df_path,
        ):
    """Takes 3 sets of gene IDs and the respective mutant name and generates
    the correspoding venn diagrams, upset plots and a dataframe for each
    one of the 7 intersections.
    """

    if len(data) == 2:

        # Differentially expressed genes
        generate_venn2_diagram(
                mutant1_gene_set=set(data[0].dge_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].dge_df.index),
                mutant2_name=data[1].name,
                plot_formats=plot_formats,
                title='Differentially expressed genes',
                path=f'{venn_path}/venn_DEG',
                )
        generate_upset_plot(
                mutant1_gene_set=set(data[0].dge_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].dge_df.index),
                mutant2_name=data[1].name,
                plot_formats=plot_formats,
                title='Differentially expressed genes',
                path=f'{upset_path}/UpSet_DEG'
                )
        mk_df_for_each_intersection2(
                mutant1_gene_set=set(data[0].dge_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].dge_df.index),
                mutant2_name=data[1].name,
                data=data,
                path=df_path,
                file_names='DEG_intersection',
                )

        # Upregulated genes
        generate_venn2_diagram(
                mutant1_gene_set=set(data[0].up_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].up_df.index),
                mutant2_name=data[1].name,
                plot_formats=plot_formats,
                title='Upregulated genes',
                path=f'{venn_path}/venn_UP',
                )
        generate_upset_plot(
                mutant1_gene_set=set(data[0].up_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].up_df.index),
                mutant2_name=data[1].name,
                plot_formats=plot_formats,
                title='Upregulated genes',
                path=f'{upset_path}/UpSet_UP'
                )
        mk_df_for_each_intersection2(
                mutant1_gene_set=set(data[0].up_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].up_df.index),
                mutant2_name=data[1].name,
                data=data,
                path=df_path,
                file_names='UP_intersection',
                )

        # Downregulated genes
        generate_venn2_diagram(
                mutant1_gene_set=set(data[0].down_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].down_df.index),
                mutant2_name=data[1].name,
                plot_formats=plot_formats,
                title='Downregulated genes',
                path=f'{venn_path}/venn_DOWN',
                )
        generate_upset_plot(
                mutant1_gene_set=set(data[0].down_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].down_df.index),
                mutant2_name=data[1].name,
                plot_formats=plot_formats,
                title='Downregulated genes',
                path=f'{upset_path}/UpSet_DOWN'
                )
        mk_df_for_each_intersection2(
                mutant1_gene_set=set(data[0].down_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].down_df.index),
                mutant2_name=data[1].name,
                data=data,
                path=df_path,
                file_names='DOWN_intersection',
                )

    elif len(data) == 3:

        # Differentially expressed genes
        generate_venn3_diagram(
                mutant1_gene_set=set(data[0].dge_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].dge_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].dge_df.index),
                mutant3_name=data[2].name,
                plot_formats=plot_formats,
                title='Differentially expressed genes',
                path=f'{venn_path}/venn_DEG',
                )
        generate_upset_plot(
                mutant1_gene_set=set(data[0].dge_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].dge_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].dge_df.index),
                mutant3_name=data[2].name,
                plot_formats=plot_formats,
                title='Differentially expressed genes',
                path=f'{upset_path}/UpSet_DEG'
                )
        mk_df_for_each_intersection3(
                mutant1_gene_set=set(data[0].dge_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].dge_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].dge_df.index),
                mutant3_name=data[2].name,
                data=data,
                path=df_path,
                file_names='DEG_intersection',
                )

        # Upregulated genes
        generate_venn3_diagram(
                mutant1_gene_set=set(data[0].up_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].up_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].up_df.index),
                mutant3_name=data[2].name,
                plot_formats=plot_formats,
                title='Upregulated genes',
                path=f'{venn_path}/venn_UP',
                )
        generate_upset_plot(
                mutant1_gene_set=set(data[0].up_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].up_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].up_df.index),
                mutant3_name=data[2].name,
                plot_formats=plot_formats,
                title='Upregulated genes',
                path=f'{upset_path}/UpSet_UP'
                )
        mk_df_for_each_intersection3(
                mutant1_gene_set=set(data[0].up_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].up_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].up_df.index),
                mutant3_name=data[2].name,
                data=data,
                path=df_path,
                file_names='UP_intersection',
                )

        # Downregulated genes
        generate_venn3_diagram(
                mutant1_gene_set=set(data[0].down_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].down_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].down_df.index),
                mutant3_name=data[2].name,
                plot_formats=plot_formats,
                title='Downregulated genes',
                path=f'{venn_path}/venn_DOWN',
                )
        generate_upset_plot(
                mutant1_gene_set=set(data[0].down_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].down_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].down_df.index),
                mutant3_name=data[2].name,
                plot_formats=plot_formats,
                title='Downregulated genes',
                path=f'{upset_path}/UpSet_DOWN'
                )
        mk_df_for_each_intersection3(
                mutant1_gene_set=set(data[0].down_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].down_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].down_df.index),
                mutant3_name=data[2].name,
                data=data,
                path=df_path,
                file_names='DOWN_intersection',
                )

    elif len(data) == 4:

        # Differentially expressed genes
        generate_venn4_diagram(
                mutant1_gene_set=set(data[0].dge_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].dge_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].dge_df.index),
                mutant3_name=data[2].name,
                mutant4_gene_set=set(data[3].dge_df.index),
                mutant4_name=data[3].name,
                plot_formats=plot_formats,
                title='Differentially expressed genes',
                path=f'{venn_path}/venn_DEG',
                )
        generate_upset_plot(
                mutant1_gene_set=set(data[0].dge_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].dge_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].dge_df.index),
                mutant3_name=data[2].name,
                mutant4_gene_set=set(data[3].dge_df.index),
                mutant4_name=data[3].name,
                plot_formats=plot_formats,
                title='Differentially expressed genes',
                path=f'{upset_path}/UpSet_DEG'
                )
        mk_df_for_each_intersection4(
                mutant1_gene_set=set(data[0].dge_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].dge_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].dge_df.index),
                mutant3_name=data[2].name,
                mutant4_gene_set=set(data[3].dge_df.index),
                mutant4_name=data[3].name,
                data=data,
                path=df_path,
                file_names='DEG_intersection',
                )

        # Upregulated genes
        generate_venn4_diagram(
                mutant1_gene_set=set(data[0].up_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].up_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].up_df.index),
                mutant3_name=data[2].name,
                mutant4_gene_set=set(data[3].dge_df.index),
                mutant4_name=data[3].name,
                plot_formats=plot_formats,
                title='Upregulated genes',
                path=f'{venn_path}/venn_UP',
                )
        generate_upset_plot(
                mutant1_gene_set=set(data[0].up_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].up_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].up_df.index),
                mutant3_name=data[2].name,
                mutant4_gene_set=set(data[3].dge_df.index),
                mutant4_name=data[3].name,
                plot_formats=plot_formats,
                title='Upregulated genes',
                path=f'{upset_path}/UpSet_UP'
                )
        mk_df_for_each_intersection4(
                mutant1_gene_set=set(data[0].up_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].up_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].up_df.index),
                mutant3_name=data[2].name,
                mutant4_gene_set=set(data[3].dge_df.index),
                mutant4_name=data[3].name,
                data=data,
                path=df_path,
                file_names='UP_intersection',
                )

        # Downregulated genes
        generate_venn4_diagram(
                mutant1_gene_set=set(data[0].down_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].down_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].down_df.index),
                mutant3_name=data[2].name,
                mutant4_gene_set=set(data[3].dge_df.index),
                mutant4_name=data[3].name,
                plot_formats=plot_formats,
                title='Downregulated genes',
                path=f'{venn_path}/venn_DOWN',
                )
        generate_upset_plot(
                mutant1_gene_set=set(data[0].down_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].down_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].down_df.index),
                mutant3_name=data[2].name,
                mutant4_gene_set=set(data[3].dge_df.index),
                mutant4_name=data[3].name,
                plot_formats=plot_formats,
                title='Downregulated genes',
                path=f'{upset_path}/UpSet_DOWN'
                )
        mk_df_for_each_intersection4(
                mutant1_gene_set=set(data[0].down_df.index),
                mutant1_name=data[0].name,
                mutant2_gene_set=set(data[1].down_df.index),
                mutant2_name=data[1].name,
                mutant3_gene_set=set(data[2].down_df.index),
                mutant3_name=data[2].name,
                mutant4_gene_set=set(data[3].dge_df.index),
                mutant4_name=data[3].name,
                data=data,
                path=df_path,
                file_names='DOWN_intersection',
                )


def get_inverted_regulations_and_mk_venns_and_dataframes(
        data,
        plot_formats,
        venn_directory_path,
        upset_directory_path,
        dataframes_directory_path,
        ):
    """Generates the venn diagramas and the correspoding intersection
    dataframe for each one of the possible inverted regulations combinations
    """

    inverted_regulation_dict = {}

    for d in data:

        inverted_regulation_dict[d.name] = {
                'Up' : {
                    'set' : set(d.up_df.index),
                    'label' : (r'$\uparrow$' + d.name),
                    },
                'Down' : {
                    'set' : set(d.down_df.index),
                    'label' : (r'$\downarrow$' + d.name),
                    },
                }

    regulations = ("Up", "Down")

    for k in inverted_regulation_dict:

        for reg in regulations:

            plot_data_dict = {}

            plot_data_dict[k] = 'Up'

            name = f'{k}_Up_others_Down'

            for notk in inverted_regulation_dict:

                if notk is not k:
                    plot_data_dict[notk] = 'Down'

                    if reg == 'Down':

                        plot_data_dict[k] = 'Down'

                        plot_data_dict[notk] = 'Up'

                        name = f'{k}_Down_others_Up'

            if len(data) == 2:
                generate_venn2_diagram(
                        mutant1_gene_set=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['set'],
                        mutant1_name=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['label'],
                        mutant2_gene_set=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['set'],
                        mutant2_name=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['label'],
                        plot_formats=plot_formats,
                        title="Inverted regulations between mutants",
                        path=f"{venn_directory_path}/{name}",
                        )
                generate_upset_plot(
                        mutant1_gene_set=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['set'],
                        mutant1_name=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['label'],
                        mutant2_gene_set=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['set'],
                        mutant2_name=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['label'],
                        plot_formats=plot_formats,
                        title='Inverted regulations between mutants',
                        path=f'{upset_directory_path}/{name}'
                        )
                mk_df_for_each_intersection2(
                        mutant1_gene_set=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['set'],
                        mutant1_name=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['label'],
                        mutant2_gene_set=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['set'],
                        mutant2_name=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['label'],
                        data=data,
                        path=dataframes_directory_path,
                        file_names=name
                        )

            elif len(data) == 3:

                generate_venn3_diagram(
                        mutant1_gene_set=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['set'],
                        mutant1_name=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['label'],
                        mutant2_gene_set=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['set'],
                        mutant2_name=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['label'],
                        mutant3_gene_set=inverted_regulation_dict[
                            data[2].name][plot_data_dict[data[2].name]]['set'],
                        mutant3_name=inverted_regulation_dict[
                            data[2].name][plot_data_dict[data[2].name]]['label'],
                        plot_formats=plot_formats,
                        title="Inverted regulations between mutants",
                        path=f"{venn_directory_path}/{name}",
                        )
                generate_upset_plot(
                        mutant1_gene_set=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['set'],
                        mutant1_name=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['label'],
                        mutant2_gene_set=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['set'],
                        mutant2_name=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['label'],
                        mutant3_gene_set=inverted_regulation_dict[
                            data[2].name][plot_data_dict[data[2].name]]['set'],
                        mutant3_name=inverted_regulation_dict[
                            data[2].name][plot_data_dict[data[2].name]]['label'],
                        plot_formats=plot_formats,
                        title='Inverted regulations between mutants',
                        path=f'{upset_directory_path}/{name}'
                        )
                mk_df_for_each_intersection3(
                        mutant1_gene_set=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['set'],
                        mutant1_name=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['label'],
                        mutant2_gene_set=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['set'],
                        mutant2_name=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['label'],
                        mutant3_gene_set=inverted_regulation_dict[
                            data[2].name][plot_data_dict[data[2].name]]['set'],
                        mutant3_name=inverted_regulation_dict[
                            data[2].name][plot_data_dict[data[2].name]]['label'],
                        data=data,
                        path=dataframes_directory_path,
                        file_names=name
                        )

            elif len(data) == 4:

                generate_venn4_diagram(
                        mutant1_gene_set=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['set'],
                        mutant1_name=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['label'],
                        mutant2_gene_set=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['set'],
                        mutant2_name=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['label'],
                        mutant3_gene_set=inverted_regulation_dict[
                            data[2].name][plot_data_dict[data[2].name]]['set'],
                        mutant3_name=inverted_regulation_dict[
                            data[2].name][plot_data_dict[data[2].name]]['label'],
                        mutant4_gene_set=inverted_regulation_dict[
                            data[3].name][plot_data_dict[data[3].name]]['set'],
                        mutant4_name=inverted_regulation_dict[
                            data[3].name][plot_data_dict[data[3].name]]['label'],
                        plot_formats=plot_formats,
                        title="Inverted regulations between mutants",
                        path=f"{venn_directory_path}/{name}",
                        )
                generate_upset_plot(
                        mutant1_gene_set=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['set'],
                        mutant1_name=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['label'],
                        mutant2_gene_set=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['set'],
                        mutant2_name=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['label'],
                        mutant3_gene_set=inverted_regulation_dict[
                            data[2].name][plot_data_dict[data[2].name]]['set'],
                        mutant3_name=inverted_regulation_dict[
                            data[2].name][plot_data_dict[data[2].name]]['label'],
                        mutant4_gene_set=inverted_regulation_dict[
                            data[3].name][plot_data_dict[data[3].name]]['set'],
                        mutant4_name=inverted_regulation_dict[
                            data[3].name][plot_data_dict[data[3].name]]['label'],
                        plot_formats=plot_formats,
                        title='Inverted regulations between mutants',
                        path=f'{upset_directory_path}/{name}'
                        )
                mk_df_for_each_intersection4(
                        mutant1_gene_set=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['set'],
                        mutant1_name=inverted_regulation_dict[
                            data[0].name][plot_data_dict[data[0].name]]['label'],
                        mutant2_gene_set=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['set'],
                        mutant2_name=inverted_regulation_dict[
                            data[1].name][plot_data_dict[data[1].name]]['label'],
                        mutant3_gene_set=inverted_regulation_dict[
                            data[2].name][plot_data_dict[data[2].name]]['set'],
                        mutant3_name=inverted_regulation_dict[
                            data[2].name][plot_data_dict[data[2].name]]['label'],
                        mutant4_gene_set=inverted_regulation_dict[
                            data[3].name][plot_data_dict[data[3].name]]['set'],
                        mutant4_name=inverted_regulation_dict[
                            data[3].name][plot_data_dict[data[3].name]]['label'],
                        data=data,
                        path=dataframes_directory_path,
                        file_names=f'{name}'
                        )

