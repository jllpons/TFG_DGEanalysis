#!/usr/bin/env python3

"""Plot generation for DGE analysis using numpy, matplotlib and seaborn.
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def generate_regulation_countplot(
        data,
        file_path,
        name,
        padj_threshold,
        foldchange_threshold,
        hue_order,
        plot_formats,
        ):
    """Creates a count plot representing the number of genes that are:
            1) Above the p-vale threshold
            2) Under the p-vale and foldchange threshold.
            3) Under the p-vale and above the foldchange threshold.
    """

    # Setting the size
    plt.figure(figsize=(10,6))
    count_plot = sns.countplot(
                        y=data,
                        palette=("silver", "cornflowerblue", "indianred"),
                        order=hue_order,
                        width=0.4,
                        )

# TODO: show percentages
#     percentage = [42.13, 29.13, 27.98]
#     patches = count_plot.patches
#     for i in range(len(patches)):
#        y = patches[i].get_x() + patches[i].get_width()/2
#        x = patches[i].get_height()+2
#        count_plot.annotate('{:.2f}%'.format(percentage[i]), (x, y), ha='center')

    title = (
            name
            +
            f" gene regulation stats for:\np-adj < {padj_threshold}"
            +
            f"\nFold Change >= {foldchange_threshold}"
            )

    plt.title(title)
    # The following two lines remove the top and right borders of the plot.
    count_plot.spines["top"].set_visible(False)
    count_plot.spines["right"].set_visible(False)

    # Label for x axis
    plt.xlabel("Differentality expressed genes", size=10)

    fig = count_plot.get_figure()

    for format in plot_formats:
        plot_name = f"{file_path}/{name}_regulation-stats.{format}"
        fig.savefig(plot_name, format=format, dpi=300)

        if format == "png":
            plot_name = f"{file_path}/{name}_regulation-stats_transparent-bg.{format}"
            fig.savefig(plot_name, format=format, dpi=300, transparent=True)

    plt.close()


def generate_volcano_plot(
        data,
        file_path,
        foldchange_threshold,
        padj_threshold,
        plot_formats,
        ):
    """Generates a volcano plot with log2 Fold change values on the x axis, and
    log10 padj values on the y axis. Colors values according to the
    corresponding preestablished threshold value for each axis.
    It also generates a countplot for a better visualization of the amount of
    genes that have been reported as Up, Down or No sig.
    """

    df = data.input_df
    log2FC_column = data.df_columns['log2FoldChange']
    padj_column = data.df_columns['padj']

    # Calculating the log2 threshold for FoldChange values
    log2FoldChange_threshold = np.log2(foldchange_threshold)
    # Same for -log10 threshold for padj values
    log10_padj_threshold = -np.log10(padj_threshold)

    # Adding -log10(padj) column to the dataframe
    df["-log10(padj)"] = -np.log10(df[padj_column])

    # Creating a colomun wich values will be:
    # -NO if the gene's padj or Foldchange values doesn't surpass the thresholds
    # -DOWN if it's foldchange value is lesser than the negative threshold
    # -UP if it's foldchange value is greater than the positive threshold
    regulation_column = []
    for a, b in zip(df[log2FC_column], df["-log10(padj)"]):

        if b < log10_padj_threshold:
            regulation_column.append("NO")
        elif a > -log2FoldChange_threshold and a < log2FoldChange_threshold:
            regulation_column.append("NO")
        elif a <= -log2FoldChange_threshold:
            regulation_column.append("DOWN")
        elif a >= log2FoldChange_threshold:
            regulation_column.append("UP")
        else:
            regulation_column.append("NO")

    # Adding this column to the dataframe
    df.insert(loc=len(df.columns), column="color", value=regulation_column,)

    # Replacing UP/DOWN/NO with the same string + the number of instances for
    # each corresponding value in the column.
    # I just think it looks better if it shows the value count
    # in the plot's legend.
    up_count = regulation_column.count("UP")
    down_count = regulation_column.count("DOWN")
    no_count = regulation_column.count("NO")
    up_newname = f"Up ({up_count})"
    down_newname = f"Down ({down_count})"
    no_newname = f"Not sig ({no_count})"

    df["color"] = df["color"].map(
            {"UP":str(up_newname),"DOWN":str(down_newname),"NO":str(no_newname)}
            )

    # Generate a countplot for a better visualization
    generate_regulation_countplot(
            data=df["color"],
            file_path=file_path,
            name=data.name,
            padj_threshold=padj_threshold,
            foldchange_threshold=foldchange_threshold,
            hue_order=[no_newname, down_newname, up_newname],
            plot_formats=plot_formats,
            )


    # Setting the size
    plt.figure(figsize=(6,7))

    # Creating the plot
    volcano = sns.scatterplot(
            data=df,
            s=15,
            linewidth=0.2,
            x=log2FC_column,
            y="-log10(padj)",
            hue="color",
            hue_order=[down_newname, no_newname, up_newname],
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
            "\nFold Change >= "
            +
            str(foldchange_threshold)
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
            #frameon=False,
            )

    # Label for x axis
    plt.xlabel("$log_{2}$ Fold change", size=12)
    # Label for y axis
    plt.ylabel("-$log_{10}$ Adjusted p-value", size=12)

    deletion = r'$\Delta$'
    # Title for the plot
    plt.title('Stationary ' + r'$vs.$' + ' Log phase')

    fig = volcano.get_figure()

    # Create the same plot in each specificed format.
    for format in plot_formats:
        plot_name = f"{file_path}/{data.name}_volcano.{format}"
        fig.savefig(plot_name, format=format, dpi=300)

        if format == "png":
            plot_name = f"{file_path}/{data.name}_volcano_transparent-bg.{format}"
            fig.savefig(plot_name, format=format, dpi=300, transparent=True)

    plt.close()

