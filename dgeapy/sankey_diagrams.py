#!/usr/bin/env python3

"""Sankey diagram generation functions. Uses numpy, pandas, matplotlib and
pySankey.
"""

import numpy as np
import matplotlib.pyplot as plt
from pysankey import sankey

def generate_sankey_diagram(
        data,
        fc_value,
        padj_value,
        plot_formats,
        path,
        ):
    """Generate sankey diagrams for gene regulations betwwen 2 mutants,
    form left to right. Will output all of the possible combinations
    """

    log2fc = np.log2(fc_value)

    for a in data:

        a.input_df.loc[a.input_df[
                            a.df_columns['padj']] > padj_value,'Regulation'
                                 ] = "No sig"

        a.input_df.loc[
            (a.input_df[a.df_columns['log2FoldChange']] < log2fc)
            & (a.input_df[a.df_columns['log2FoldChange']] > -log2fc), 'Regulation'
                      ] = "No sig"

        a.input_df.loc[(
            a.input_df['Regulation'] != "No sig")
            & (a.input_df[a.df_columns['log2FoldChange']] >= log2fc),
                'Regulation'
                ] = "Up"

        a.input_df.loc[
                (a.input_df['Regulation'] != "No sig") 
                & (a.input_df[a.df_columns['log2FoldChange']] <= -log2fc),
                'Regulation'
                ] = "Down"

        a.input_df = a.input_df[a.input_df['Regulation'] != 'No sig']

        for b in data:

            if b is not a:

                b.input_df.loc[b.input_df[
                                    b.df_columns['padj']] > padj_value,'Regulation'
                                         ] = "No sig"

                b.input_df.loc[
                    (b.input_df[b.df_columns['log2FoldChange']] < log2fc)
                    & (b.input_df[b.df_columns['log2FoldChange']] > -log2fc), 'Regulation'
                              ] = "No sig"

                b.input_df.loc[(
                    b.input_df['Regulation'] != "No sig")
                    & (b.input_df[b.df_columns['log2FoldChange']] >= log2fc),
                        'Regulation'
                        ] = "Up"

                b.input_df.loc[
                        (b.input_df['Regulation'] != "No sig") 
                        & (b.input_df[b.df_columns['log2FoldChange']] <= -log2fc),
                        'Regulation'
                        ] = "Down"

                b.input_df = b.input_df[b.input_df['Regulation'] != 'No sig']

                df = a.input_df.merge(
                                b.input_df,
                                how='outer',
                                on=a.df_columns['geneID'],
                                suffixes=('_a', '_b')
                                )

                df.loc[df['Regulation_a'].isnull(), 'Regulation_a'] = 'No sig'
                df.loc[df['Regulation_b'].isnull(), 'Regulation_b'] = 'No sig'

                df = df[['Regulation_a', 'Regulation_b']]

                colordict = {
                        'No sig' : "silver",
                        'Down' : "cornflowerblue",
                        'Up' : "indianred",
                        }

                sankey(
                        df['Regulation_a'],
                        df['Regulation_b'],
                        leftLabels=['Down', 'No sig', 'Up',],
                        rightLabels=['Down', 'No sig', 'Up',],
                        colorDict=colordict,
                        fontsize=8,
                        )

                plt.gcf().set_size_inches(6 ,6)

                plt.title(
                    f'Flow of Diffetentially Expressed Genes\n from {a.name} to {b.name}'
                        )

                for format in plot_formats:
                    plt.savefig(f"{path}/{a.name}_vs_{b.name}_sankey.{format}", dpi=300)

                plt.clf()
                plt.close()

# def generate_sankey_diagram(
#         data,
#         fc_value,
#         padj_value,
#         plot_formats,
#         path,
#         ):
#     """Generate sankey diagrams for gene regulations betwwen 2 mutants,
#     form left to right. Will output all of the possible combinations
#     """
# 
#     log2fc = np.log2(fc_value)
# 
#     for a in data:
# 
#         a.input_df.loc[a.input_df[
#                             a.df_columns['padj']] > padj_value,'Regulation'
#                                  ] = "No sig"
# 
#         a.input_df.loc[
#             (a.input_df[a.df_columns['log2FoldChange']] < log2fc)
#             & (a.input_df[a.df_columns['log2FoldChange']] > -log2fc), 'Regulation'
#                       ] = "No sig"
# 
#         a.input_df.loc[(
#             a.input_df['Regulation'] != "No sig")
#             & (a.input_df[a.df_columns['log2FoldChange']] >= log2fc),
#                 'Regulation'
#                 ] = "Up"
# 
#         a.input_df.loc[
#                 (a.input_df['Regulation'] != "No sig") 
#                 & (a.input_df[a.df_columns['log2FoldChange']] <= -log2fc),
#                 'Regulation'
#                 ] = "Down"
# 
#         for b in data:
# 
#             if b is not a:
# 
#                 b.input_df.loc[b.input_df[
#                                     b.df_columns['padj']] > padj_value,'Regulation'
#                                          ] = "No sig"
# 
#                 b.input_df.loc[
#                     (b.input_df[b.df_columns['log2FoldChange']] < log2fc)
#                     & (b.input_df[b.df_columns['log2FoldChange']] > -log2fc), 'Regulation'
#                               ] = "No sig"
# 
#                 b.input_df.loc[(
#                     b.input_df['Regulation'] != "No sig")
#                     & (b.input_df[b.df_columns['log2FoldChange']] >= log2fc),
#                         'Regulation'
#                         ] = "Up"
# 
#                 b.input_df.loc[
#                         (b.input_df['Regulation'] != "No sig") 
#                         & (b.input_df[b.df_columns['log2FoldChange']] <= -log2fc),
#                         'Regulation'
#                         ] = "Down"
# 
#                 df = a.input_df.merge(
#                                 b.input_df,
#                                 how='outer',
#                                 on='gene_id',
#                                 suffixes=('_a', '_b')
#                                 )
# 
#                 df.loc[df['Regulation_a'].isnull(), 'Regulation_a'] = 'No sig'
#                 df.loc[df['Regulation_b'].isnull(), 'Regulation_b'] = 'No sig'
# 
#                 df = df[['Regulation_a', 'Regulation_b']]
# 
#                 a_upcounts = (df['Regulation_a'] == "Up").sum()
#                 a_nosigcounts = (df['Regulation_a'] == "No sig").sum()
#                 a_downcounts = (df['Regulation_a'] == "Down").sum()
# 
#                 a_total = a_upcounts + a_nosigcounts + a_downcounts
#                 a_uppercentage = (100 * (a_upcounts / a_total))
#                 a_nosigpercentage = (100 * (a_nosigcounts / a_total))
#                 a_downpercentage = (100 * (a_downcounts / a_total))
# 
#                 a_uplabel = f"Up\n{a_upcounts}\n({a_uppercentage:.1f}%)"
#                 a_nosiglabel = f"No sig\n{a_nosigcounts}\n({a_nosigpercentage:.1f}%)"
#                 a_downlabel = f"Down\n{a_downcounts}\n({a_downpercentage:.1f}%)"
# 
#                 b_upcounts = (df['Regulation_b'] == "Up").sum()
#                 b_nosigcounts = (df['Regulation_b'] == "No sig").sum()
#                 b_downcounts = (df['Regulation_b'] == "Down").sum()
# 
#                 b_total = b_upcounts + b_nosigcounts + b_downcounts
#                 b_uppercentage = (100 * (b_upcounts / b_total))
#                 b_nosigpercentage = (100 * (b_nosigcounts / b_total))
#                 b_downpercentage = (100 * (b_downcounts / b_total))
# 
#                 b_uplabel = f"Up\n{b_upcounts}\n({b_uppercentage:.1f}%)"
#                 b_nosiglabel = f"No sig\n{b_nosigcounts}\n({b_nosigpercentage:.1f}%)"
#                 b_downlabel = f"Down\n{b_downcounts}\n({b_downpercentage:.1f}%)"
# 
#                 df['Regulation_a'] = df['Regulation_a'].map({
#                     "Up" : a_uplabel,
#                     "No sig" : a_nosiglabel,
#                     "Down" : a_downlabel,
#                     })
#                 df['Regulation_b'] = df['Regulation_b'].map({
#                     "Up" : b_uplabel,
#                     "No sig" : b_nosiglabel,
#                     "Down" : b_downlabel,
#                     })
# 
#                 colordict = {
#                         a_nosiglabel : "silver",
#                         a_downlabel : "cornflowerblue",
#                         a_uplabel : "indianred",
#                         b_nosiglabel : "silver",
#                         b_downlabel : "cornflowerblue",
#                         b_uplabel : "indianred",
#                         }
# 
#                 sankey(
#                         df['Regulation_a'],
#                         df['Regulation_b'],
#                         leftLabels=[a_downlabel, a_nosiglabel, a_uplabel,],
#                         rightLabels=[b_downlabel, b_nosiglabel, b_uplabel,],
#                         colorDict=colordict,
#                         fontsize=8,
#                         )
# 
#                 plt.gcf().set_size_inches(4 ,6)
# 
#                 plt.title(
#                     f'Flow of Diffetentially Expressed Genes\n from {a.name} to {b.name}'
#                         )
# 
#                 for format in plot_formats:
#                     plt.savefig(f"{path}/{a.name}_vs_{b.name}_sankey.{format}", dpi=300)
# 
#                 plt.clf()
#                 plt.close()
# 
# 
