#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_unweighted, venn3, venn3_unweighted

def generate_venn2_diagram(
        mutant1_gene_set,
        mutant1_name,
        mutant2_gene_set,
        mutant2_name,
        plot_formats,
        title,
        ):
    """Generate a Venn diagram from 2 given sets and their respective names.
    """

    plt.style.use(['default'])

    set_array = [mutant1_gene_set, mutant2_gene_set]

    set_names = [mutant1_name, mutant2_name]

    venn2(set_array[0:2], set_names[0:2])

    plt.title(title)

    for format in plot_formats:
        plt.savefig(f"venn2.{format}", dpi=300, transparent=True)

        #if format == "png":
        #    plt.savefig(f"{path}_transparent-bg.{format}", dpi=300, transparent=True)

    plt.close()

    venn2_unweighted(set_array[0:3], set_names[0:3])

    plt.title(title)

#     for format in plot_formats:
#         plt.savefig(f"venn2_unweighted.{format}", dpi=300)
# 
#         if format == "png":
#             plt.savefig(
#                     f"{path}_unweighted_transparent-bg.{format}",
#                     dpi=300,
#                     transparent=True
#                     )
mutant1deg = pd.read_csv('M30 rpfC-2/M30 rpfC-2_DEG.tsv', sep='\t', index_col='index')
mutant2deg = pd.read_csv('M30 rpfF-2/M30 rpfF-2_DEG.tsv', sep='\t', index_col='index')

delta = r'$\Delta$'
generate_venn2_diagram(
        mutant1_name=f'M30 {delta}' + r'$\mathit{rpfC}$',
        mutant1_gene_set=set(mutant1deg.index.values),
        mutant2_name=f'M30 {delta}' + r'$\mathit{rpfF}$',
        mutant2_gene_set=set(mutant2deg.index.values),
        plot_formats=['png'],
        title='Differentially expressed genes',
        )
