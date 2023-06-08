#!/usr/bin/env python3

import numpy as np
import pandas as pd
from upsetplot import UpSet, from_contents
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib_venn import venn3_unweighted

mutant1up = pd.read_csv('clp/clp_UP.tsv', sep='\t', index_col='index')
mutant1down = pd.read_csv('clp/clp_DOWN.tsv', sep='\t', index_col='index')
mutant2up = pd.read_csv('rpfC/rpfC_UP.tsv', sep='\t', index_col='index')
mutant2down = pd.read_csv('rpfC/rpfC_DOWN.tsv', sep='\t', index_col='index')
mutant3up = pd.read_csv('rpfF/rpfF_UP.tsv', sep='\t', index_col='index')
mutant3down = pd.read_csv('rpfF/rpfF_DOWN.tsv', sep='\t', index_col='index')

bool_df = from_contents({
    r'Sobreexpressats a $\Delta$clp' : mutant1up.index.values,
    r'Sobreexpressats a $\Delta$rpfC' : mutant2up.index.values,
    r'Sobreexpressats a $\Delta$rpfF' : mutant3up. index.values,
    r'Subexpressats a $\Delta$clp' : mutant1down.index.values,
    r'Subexpressats a $\Delta$rpfC' : mutant2down.index.values,
    r'Subexpressats a $\Delta$rpfF' : mutant3down. index.values
    })

functions = mutant1up['function'].to_frame()
functions = pd.concat([functions, mutant1down['function'].to_frame()])
functions = pd.concat([functions, mutant2up['function'].to_frame()])
functions = pd.concat([functions, mutant2down['function'].to_frame()])
functions = pd.concat([functions, mutant3up['function'].to_frame()])
functions = pd.concat([functions, mutant3down['function'].to_frame()])

functions = functions.loc[~functions.index.duplicated(keep='first')]

bool_df = bool_df.join(functions, how='inner', on='id')

bool_df['function'] = bool_df['function'].map(
        {'not annotated' : '5. Sense anotació',
         'others' : '4. Altres',
         'biofilm formation' : '3. Formació de biofilm',
         'antimicrobial resistance' : '2. Resistència als antimicrobians',
         'motility' : '1. Motilitat'
         })

upset = UpSet(
        bool_df,
        intersection_plot_elements=0,
        sort_by="cardinality",
        element_size=45,
        show_counts=True,
        sort_categories_by='-input',
        )

cmap = {
        '5. Sense anotació' : '#080708',
        '4. Altres' : 'silver',
        '3. Formació de biofilm' : '#FFB000',
        '2. Resistència als antimicrobians' : '#DC267F',
        '1. Motilitat' : '#648FFF',
        }

upset.add_stacked_bars(
        by='function',
        colors=cmap,
        title='nombre de gens i funció',
        elements=10
        )

upset.plot()

plt.suptitle('Gene distributions between mutants')

plt.savefig('UpSet_functions_and_rest_6sets.png', dpi=300)


