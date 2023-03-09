#!/usr/bin/env python3

"""Sankey diagram generation functions. Uses numpy, pandas, matplotlib and
pySankey.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pysankey import sankey

#------------------------------------------------------------------------------#

def generate_sankey_diagram(
        dataframes,
        mutants,
        fc_column_names,
        padj_column_names,
        fc_value,
        padj_value,
        plot_formats,
        path,
        ):
    """Generate sankey diagrams for gene regulations betwwen 2 mutants,
    form left to right. Will output all of the possible combinations
    """

    log2fc = np.log2(fc_value)

    for a, b, c, m in zip(dataframes, fc_column_names, padj_column_names, mutants):

        a.loc[a[c] > padj_value,"Regulation"] = "No sig"
        a.loc[(a[b] < log2fc) & (a[b] > -log2fc), "Regulation"] = "No sig"
        a.loc[(a["Regulation"] != "No sig")
              & (a[b] >= log2fc), 
                "Regulation"
                ] = "Up"
        a.loc[
                (a["Regulation"] != "No sig") 
                & (a[b] <= -log2fc),
                "Regulation"
                ] = "Down"
        a.loc[a["Regulation"].isnull(), "Regulation"] = "No sig"

        for x, y, z, n in zip(dataframes, fc_column_names, padj_column_names, mutants):

            if x is not a:
                x.loc[x[z] > padj_value,"Regulation"] = "No sig"
                x.loc[(x[y] < log2fc) & (a[y] > -log2fc), "Regulation"] = "No sig"
                x.loc[(x["Regulation"] != "No sig")
                      & (x[y] >= log2fc), 
                        "Regulation"
                        ] = "Up"
                x.loc[
                        (x["Regulation"] != "No sig") 
                        & (x[y] <= -log2fc),
                        "Regulation"
                        ] = "Down"
                x.loc[x["Regulation"].isnull(), "Regulation"] = "No sig"

                df = pd.concat([a["Regulation"], x["Regulation"]],
                               axis=1,
                               keys=["a_Regulation", "x_Regulation"]
                               )

                a_upcounts = (df["a_Regulation"] == "Up").sum()
                a_nosigcounts = (df["a_Regulation"] == "No sig").sum()
                a_downcounts = (df["a_Regulation"] == "Down").sum()

                a_total = a_upcounts + a_nosigcounts + a_downcounts
                a_uppercentage = (100 * (a_upcounts / a_total))
                a_nosigpercentage = (100 * (a_nosigcounts / a_total))
                a_downpercentage = (100 * (a_downcounts / a_total))

                a_uplabel = f"Up\n{a_upcounts}\n({a_uppercentage})"
                a_nosiglabel = f"No sig\n{a_nosigcounts}\n({a_nosigpercentage})"
                a_downlabel = f"Down\n{a_downcounts}\n({a_downpercentage})"

                x_upcounts = (df["x_Regulation"] == "Up").sum()
                x_nosigcounts = (df["x_Regulation"] == "No sig").sum()
                x_downcounts = (df["x_Regulation"] == "Down").sum()

                x_total = x_upcounts + x_nosigcounts + x_downcounts
                x_uppercentage = (100 * (x_upcounts / x_total))
                x_nosigpercentage = (100 * (x_nosigcounts / x_total))
                x_downpercentage = (100 * (x_downcounts / x_total))

                x_uplabel = f"Up\n{x_upcounts}\n({x_uppercentage})"
                x_nosiglabel = f"No sig\n{x_nosigcounts}\n({x_nosigpercentage})"
                x_downlabel = f"Down\n{x_downcounts}\n({x_downpercentage})"

                df["a_Regulation"] = df["a_Regulation"].map({
                    "Up" : a_uplabel,
                    "No sig" : a_nosiglabel,
                    "Down" : a_downlabel,
                    })
                df["x_Regulation"] = df["x_Regulation"].map({
                    "Up" : x_uplabel,
                    "No sig" : x_nosiglabel,
                    "Down" : x_downlabel,
                    })

                colordict = {
                        a_nosiglabel : "silver",
                        a_downlabel : "cornflowerblue",
                        a_uplabel : "indianred",
                        x_nosiglabel : "silver",
                        x_downlabel : "cornflowerblue",
                        x_uplabel : "indianred",
                        }

                sankey(
                        df["a_Regulation"],
                        df["x_Regulation"],
                        leftLabels=[a_downlabel, a_nosiglabel, a_uplabel,],
                        rightLabels=[x_downlabel, x_nosiglabel, x_uplabel,],
                        colorDict=colordict,
                        fontsize=8,
                        )

                plt.gcf().set_size_inches(4 ,6)

                plt.title(
                        "Flow of Diffetentially Expressed Genes:\n" \
                        f"from {m} to {n}."
                        )

                for format in plot_formats:
                    plt.savefig(f"{path}{m}vs{n}sankey.{format}", dpi=300)

                plt.close()






