#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

"""
Plot generation for DGE analysis using numpy, pandas, matplotlib and seaborn.
"""


#-------# Function Definitions #-----------------------------------------------#


def generate_regulation_countplot(
        data,
        file_path,
        mutant_name,
        padj_threshold,
        foldchange_threshold,
        plot_formats,
        ):
    """
    Creates a count plot representing the number of genes that are:
        1) Above the p-vale threshold
        2) Under the p-vale and foldchange threshold.
        3) Under the p-vale and above the foldchange threshold.
    """

    df = pd.DataFrame(data, columns=[" "])


    count_plot = sns.countplot(
            x=df[" "],
            palette=("silver", "cornflowerblue", "indianred"),
            )

    title = (
            mutant_name
            +
            f" gene regulation stats for:\np-adj < {padj_threshold}"
            +
            f"\nFold Change > {foldchange_threshold}"
            )

    plt.title(title)
    # The following two lines remove the top and right borders of the plot.
    count_plot.spines["top"].set_visible(False)
    count_plot.spines["right"].set_visible(False)

    fig = count_plot.get_figure()

    for format in plot_formats:
        plot_name = f"{file_path}_regulation-stats.{format}"
        fig.savefig(plot_name, format=format)

    # Clear current figure
    plt.clf()


def generate_volcano_plot(
        dataframe,
        file_path,
        x_axis_values,
        y_axis_values,
        foldchange_threshold,
        padj_threshold,
        mutant_name,
        plot_formats,
        ):
    """
    Generates a volcano plot with log2 Fold change values on the x axis, and
    log10 padj values on the y axis. Colors values according to the
    corresponding preestablished threshold value for each axis.
    It also generates a countplot for a better visualization of the amount of
    genes that have been reported as Up, Down or No sig.
    """

    # Calculating the log2 threshold for FoldChange values
    log2FoldChange_threshold = np.log2(foldchange_threshold)
    # Same for -log10 threshold for padj values
    log10_padj_threshold = -np.log10(padj_threshold)

    # Adding -log10(padj) column to the dataframe
    dataframe["-log10(padj)"] = -np.log10(dataframe[y_axis_values])

    # Creating a colomun wich values will be:
    # -NO if the gene's padj or Foldchange values doesn't surpass the thresholds
    # -DOWN if it's foldchange value is lesser than the negative threshold
    # -UP if it's foldchange value is greater than the positive threshold
    regulation_column = []
    for a, b in zip(dataframe[x_axis_values], dataframe["-log10(padj)"]):

        if b < log10_padj_threshold:
            regulation_column.append("NO")
        elif a >= -log2FoldChange_threshold and a <= log2FoldChange_threshold:
            regulation_column.append("NO")
        elif a < -log2FoldChange_threshold:
            regulation_column.append("DOWN")
        elif a > log2FoldChange_threshold:
            regulation_column.append("UP")
        else:
            regulation_column.append("NO")


    # Replacing UP/DOWN/NO with the same string + the number of instances for
    # each corresponding value in the dataframe.
    # I just think it looks better if it shows the value count
    # in the plot's legend.
    up_count = regulation_column.count("UP")
    down_count = regulation_column.count("DOWN")
    no_count = regulation_column.count("NO")
    up_newname = f"Up ({up_count})"
    down_newname = f"Down ({down_count})"
    no_newname = f"Not sig ({no_count})"

    color_column = []
    for element in regulation_column:
        if element == "UP":
            color_column.append(up_newname)
        elif element == "DOWN":
            color_column.append(down_newname)
        elif element == "NO":
            color_column.append(no_newname)

    # Generate a countplot for a better visualization
    generate_regulation_countplot(
            data=color_column,
            file_path=file_path,
            mutant_name=mutant_name,
            padj_threshold=padj_threshold,
            foldchange_threshold=foldchange_threshold,
            plot_formats=plot_formats,
            )

    # Adding this column to the dataframe
    dataframe.insert(
            loc=1,
            column="color",
            value=color_column,
            )

    # Setting the size
    plt.figure(figsize=(6,7))

    sns.set_palette("muted")

    # Creating the plot
    volcano = sns.scatterplot(
            data=dataframe,
            x=x_axis_values,
            y="-log10(padj)",
            hue="color",
            hue_order=[ down_newname, no_newname, up_newname ],
            palette=("cornflowerblue", "silver", "indianred"),
            )

    # Add lines that will better ilustrate the choosen thresholds:
    # zorder: add to the bottom (all other elements will be added in front)
    # c is for color
    # lw is for line with
    # ls is for line simbol
    volcano.axhline(log10_padj_threshold, zorder=0, c="grey", lw=1, ls="-." )
    volcano.axvline(log2FoldChange_threshold, zorder=0, c="grey",lw=1,ls="-.")
    volcano.axvline(-log2FoldChange_threshold, zorder=0, c="grey",lw=1,ls="-.")
    # The following two lines remove the top and right borders of the plot.
    volcano.spines["top"].set_visible(False)
    volcano.spines["right"].set_visible(False)
    # Creating the legend title. Apparently seaborn/matplotlib uses latex and
    # formated strings have given me problems.
    legend_title = (
            "Adjusted p-value < "
            +
            str(padj_threshold)
            +
            "\n\t-$log_{10}$= "
            +
            str(round(log10_padj_threshold, 3))
            +
            "\nFold Change > "
            +
            str(foldchange_threshold)
            +
            "\n\t$log_{2}$= "
            +
            str(round(log2FoldChange_threshold,3))
            )
    # Modifying legend title. Otherwise it'd be "color"
    volcano.legend(title=legend_title)
    # Add more space to the right
    plt.subplots_adjust(right=0.75)
    # Place the legend into hte bbox_to_anchor coordinates.
    sns.move_legend(
            volcano,
            loc=1,
            bbox_to_anchor=(1.4,0.5),
            frameon=False,
            )

    # Label for x axis
    plt.xlabel("$log_{2}$ Fold change", size=12)
    # Label for y axis
    plt.ylabel("-$log_{10}$ Adjusted p-value", size=12)
    # Title for the plot
    plt.title(f"{mutant_name} vs WT")

    fig = volcano.get_figure()

    # Create the same plot in each specificed format.
    for format in plot_formats:
        plot_name = file_path + "_volcano." + format
        fig.savefig(plot_name, format=format)
    plt.clf()


