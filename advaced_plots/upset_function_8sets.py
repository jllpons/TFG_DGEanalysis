#!/usr/bin/env python3

import numpy as np
import pandas as pd
from upsetplot import UpSet, from_contents
from matplotlib import pyplot as plt
from matplotlib import cm

mutant1up = pd.read_csv('k279 rpfC-1/k279 rpfC-1_UP.tsv', sep='\t', index_col='index')
mutant1down = pd.read_csv('k279 rpfC-1/k279 rpfC-1_DOWN.tsv', sep='\t', index_col='index')

mutant2up = pd.read_csv('k279 rpfF-1/k279 rpfF-1_UP.tsv', sep='\t', index_col='index')
mutant2down = pd.read_csv('k279 rpfF-1/k279 rpfF-1_DOWN.tsv', sep='\t', index_col='index')

mutant3up = pd.read_csv('M30 rpfC-2/M30 rpfC-2_UP.tsv', sep='\t', index_col='index')
mutant3down = pd.read_csv('M30 rpfC-2/M30 rpfC-2_DOWN.tsv', sep='\t', index_col='index')

mutant4up = pd.read_csv('M30 rpfF-2/M30 rpfF-2_UP.tsv', sep='\t', index_col='index')
mutant4down = pd.read_csv('M30 rpfF-2/M30 rpfF-2_DOWN.tsv', sep='\t', index_col='index')

bool_df = from_contents({
    r'Upregulated in K279a $\Delta$' + r'$\mathit{rpfC}$' : mutant1up.index.values,
    r'Upregulated in K279a $\Delta$' + r'$\mathit{rpfF}$' : mutant2up.index.values,
    r'Upregulated in M30 $\Delta$' + r'$\mathit{rpfC}$' : mutant3up. index.values,
    r'Upregulated in M30 $\Delta$' + r'$\mathit{rpfF}$' : mutant4up. index.values,
    r'Downregulated in K279a $\Delta$' + r'$\mathit{rpfC}$' : mutant1down.index.values,
    r'Downregulated in K279a $\Delta$' + r'$\mathit{rpfF}$' : mutant2down.index.values,
    r'Downregulated in M30 $\Delta$' + r'$\mathit{rpfC}$' : mutant3down. index.values,
    r'Downregulated in M30 $\Delta$' + r'$\mathit{rpfF}$' : list(set(mutant4down.index.values))
    })

functions = mutant1up['function'].to_frame()
functions = pd.concat([functions, mutant1down['function'].to_frame()])
functions = pd.concat([functions, mutant2up['function'].to_frame()])
functions = pd.concat([functions, mutant2down['function'].to_frame()])
functions = pd.concat([functions, mutant3up['function'].to_frame()])
functions = pd.concat([functions, mutant3down['function'].to_frame()])
functions = pd.concat([functions, mutant4up['function'].to_frame()])
functions = pd.concat([functions, mutant4down['function'].to_frame()])

functions = functions.loc[~functions.index.duplicated(keep='first')]

bool_df = bool_df.join(functions, how='inner', on='id')

bool_df['function'] = bool_df['function'].map(
        {'not annotated' : '5. Not annotated',
         'others' : '4. Others',
         'biofilm formation' : '3. Biofilm formation',
         'antimicrobial resistance' : '2. Antimicrobial resistance',
         'motility' : '1. Motility'
         })

upset = UpSet(
        bool_df,
        intersection_plot_elements=0,
        sort_by='cardinality',
        element_size=45,
        show_counts=True,
        sort_categories_by='-input',
        min_subset_size=6,
        )

cmap = {
        '5. Not annotated' : '#080708',
        '4. Others' : 'silver',
        '3. Biofilm formation' : '#FFB000',
        '2. Antimicrobial resistance' : '#DC267F',
        '1. Motility' : '#648FFF',
        }

upset.add_stacked_bars(
        by='function',
        colors=cmap,
        title='Number of genes and function',
        elements=10
        )

upset.plot()

plt.suptitle('Gene distributions between mutants')

plt.savefig('UpSet_functions_and_rest_8sets.png', dpi=300, transparent=True)
