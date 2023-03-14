#!/usr/bin/env python3

"""Upset plot generation functions. Uses matplotlib and UpSetPlot
"""

import matplotlib.pyplot as plt
from upsetplot import from_contents, UpSet


def generate_upset_plot(
        mutant1_gene_set,
        mutant1_name,
        mutant2_gene_set,
        mutant2_name,
        mutant3_gene_set,
        mutant3_name,
        plot_formats,
        title,
        path,
        ):
    """Generate Upset plot from 3 given sets of gene IDs and their
    respective mutant names.
    """

    plt.style.use(['default'])

    upset_sets = from_contents({
                    mutant1_name : list(mutant1_gene_set),
                    mutant2_name : list(mutant2_gene_set),
                    mutant3_name : list(mutant3_gene_set),
                    })

    upset_plot = UpSet(
                    upset_sets,
                    subset_size="count",
                    show_counts="{:d}",
                    sort_by="cardinality",
                    element_size=55,
                    ).plot()

    plt.suptitle(title)

    for format in plot_formats:
        plt.savefig(f"{path}.{format}", dpi=300)

        if format == "png":
            plt.savefig(
                    f"{path}_transparent-bg.{format}",
                    dpi=300,
                    transparent=True
                    )

    plt.close()

