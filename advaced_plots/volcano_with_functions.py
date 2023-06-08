import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def mapcolor(data):
    log2FoldChange, log10padj, function = data

    if abs(log2FoldChange) < np.log2(1.5) or log10padj < -np.log10(0.05):
        return 'No significant'
    if function == 'motility':
        return 'Motility'
    elif function == 'biofilm formation':
        return 'Biofilm formation'
    elif function == 'antimicrobial resistance':
        return 'Antimicrobial resistance'

    return 'Significant'

def generate_volcano_plot(
        data,
        foldchange_threshold,
        padj_threshold,
        ):
    """Generates a volcano plot with log2 Fold change values on the x axis, and
    log10 padj values on the y axis. Colors values according to the
    corresponding preestablished threshold value for each axis.
    It also generates a countplot for a better visualization of the amount of
    genes that have been reported as Up, Down or No sig.
    """

    df = data
    log2FC_column = data['log2FoldChange']

    # Calculating the log2 threshold for FoldChange values
    log2FoldChange_threshold = np.log2(foldchange_threshold)
    # Same for -log10 threshold for padj values
    log10_padj_threshold = -np.log10(padj_threshold)

    # Adding -log10(padj) column to the dataframe
    df["-log10(padj)"] = -np.log10(df['padj'])

    df['color'] = df[['log2FoldChange', '-log10(padj)', 'function']].apply(mapcolor ,axis=1)

    # Setting the size

    no_funct_df = df.loc[(df['color'] == 'Significant') | (df['color'] == 'No significant')]
    funct_df = df.loc[(df['color'] != 'Significant') & (df['color'] != 'No significant')]

    cmap = {'No significant' : 'lightgray', 'Significant' : 'Grey'}


    fig, ax = plt.subplots(figsize=(6,6))
    # Creating the plot
    volcano = sns.scatterplot(
            data=no_funct_df,
            x=log2FC_column,
            linewidth=0.2,
            s=20,
            y="-log10(padj)",
            hue="color",
            hue_order=['No significant', 'Significant'],
            palette=cmap,
            ax=ax,
            )
    volcano = sns.scatterplot(
            data=funct_df,
            x=log2FC_column,
            linewidth=0.2,
            s=20,
            y="-log10(padj)",
            hue="color",
            hue_order=['Biofilm formation', 'Antimicrobial resistance', 'Motility'],
            palette=('#FFB000', '#DC267F','#648FFF'),
            ax=ax,
            )

    # Add lines that will better ilustrate the choosen thresholds:
    # zorder: add to the bottom (all other elements will be added in front)
    # c is for color
    # lw is for line with
    # ls is for line simbol
    plt.axhline(log10_padj_threshold, zorder=0, c="grey", lw=1, ls="-." )
    plt.axvline(log2FoldChange_threshold, zorder=0, c="grey",lw=1,ls="-.")
    plt.axvline(-log2FoldChange_threshold, zorder=0, c="grey",lw=1,ls="-.")


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
            "\nFold Change >= |"
            +
            str(foldchange_threshold)
            +
            '|'
            )
    # Modifying legend title. Otherwise it'd be "color"
    plt.legend(title=legend_title)
    # Add more space to the right
    plt.subplots_adjust(right=0.75)
    # Place the legend into hte bbox_to_anchor coordinates.
    sns.move_legend(
            volcano,
            loc=1,
            bbox_to_anchor=(1.4,0.9),
            )
    ax.set_xlim(-9, 9)
    ax.set_ylim(0, 275)

    # Label for x axis
    plt.xlabel("$log_{2}$ Fold change", size=12)
    # Label for y axis
    plt.ylabel("-$log_{10}$ Adjusted p-value", size=12)
    # Title for the plot
    delta = r'$\Delta$'
    vs = r'$\mathit{vs}$'
    plt.title(f'K279a {delta}' + r'$\mathit{rpfF}$' + f' {vs} Wild Type')

    fig = volcano.get_figure()

    # Create the same plot in each specificed format.
    fig.savefig('rpfF-1.png', dpi=300)
    fig.savefig('rpfF-1_nobg.png', dpi=300, transparent=True)

    plt.close()


#df = pd.read_csv('clp/clp.tsv', sep='\t')

df = pd.read_excel('k279 rpfF-1/k279 rpfF-1_input.xlsx')

generate_volcano_plot(
        data=df,
        foldchange_threshold=1.5,
        padj_threshold=0.05
        )

