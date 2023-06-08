#!/usr/bin/env python3

import pandas as pd
from upsetplot import UpSet, from_contents
from matplotlib import pyplot as plt

mutant1up = pd.read_csv('rpfC/rpfC_UP.tsv', sep='\t', index_col='index')
mutant1up = mutant1up[mutant1up.index.notnull()]
mutant1down = pd.read_csv('rpfC/rpfC_DOWN.tsv', sep='\t', index_col='index')
mutant1down = mutant1down[mutant1down.index.notnull()]

mutant2up = pd.read_csv('rpfF/rpfF_UP.tsv', sep='\t', index_col='index')
mutant2up = mutant2up[mutant2up.index.notnull()]
mutant2down = pd.read_csv('rpfF/rpfF_DOWN.tsv', sep='\t', index_col='index')
mutant2down = mutant2down[mutant2down.index.notnull()]


bool_df = from_contents({
    (r'Upregulated in K279a $\Delta$' + r'$\mathit{rpfC}$') : list(set(mutant1up.index.values)),
    (r'Upregulated in K279a $\Delta$' + r'$\mathit{rpfF}$') : list(set(mutant2up.index.values)),
    (r'Downregulated in K279a $\Delta$' + r'$\mathit{rpfC}$') : list(set(mutant1down.index.values)),
    (r'Downregulated in K279a $\Delta$' + r'$\mathit{rpfF}$') : list(set(mutant2down.index.values))
    })

for value in list(set(mutant1up.index.values)):
    if value in list((mutant1down.index.values)):
        print(value)

functions = mutant1up['function'].to_frame()
functions = pd.concat([functions, mutant1down['function'].to_frame()])
functions = pd.concat([functions, mutant2up['function'].to_frame()])
functions = pd.concat([functions, mutant2down['function'].to_frame()])

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
        sort_by="cardinality",
        #element_size=45,
        show_counts=True,
        sort_categories_by='-input',
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

plt.suptitle('Nombre de gens i funció')

plt.savefig('UpSet_M30.png', dpi=300, transparent=True)
