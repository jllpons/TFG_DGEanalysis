#!/usr/bin/env python3

import numpy as np
import pandas as pd
from upsetplot import UpSet, from_contents
from matplotlib import pyplot as plt
from matplotlib import cm

mutant1 = pd.read_csv('clp_DEG.tsv', sep='\t')
mutant1 = mutant1.set_index('index')
mutant2 = pd.read_csv('rpfC_DEG.tsv', sep='\t')
mutant2 = mutant2.set_index('index')
mutant3 = pd.read_csv('rpfF_DEG.tsv', sep='\t')
mutant3 = mutant3.set_index('index')

bool_df = from_contents({
    r'DEG in $\Delta$clp' : mutant1.index.values,
    r'DEG in $\Delta$rpfC' : mutant2.index.values,
    r'DEG in $\Delta$rpfF' : mutant3. index.values
    })

# A man of culture:
mutant1fcname = r'$\Delta$clp $log_{2}$ fold change'
mutant1 = mutant1.rename(columns={'log2FoldChange' : mutant1fcname})
fc1 = mutant1[mutant1fcname]
# bool_df = bool_df.join(fc1, how='outer', on='id')

mutant2fcname = r'$\Delta$rpfC $log_{2}$ fold change'
mutant2 = mutant2.rename(columns={'log2FoldChange' : mutant2fcname})
fc2 = mutant2[mutant2fcname]
# bool_df = bool_df.join(fc2, how='outer', on='id')

upset = UpSet(
        bool_df,
        sort_by="cardinality",
        show_counts=True,
        element_size=55,
        sort_categories_by='-input',
        )

# upset.add_catplot(
#         value=mutant1fcname,
#         kind='strip',
#         color='black',
#         )
# 
upset.plot()


plt.savefig('upset.png', dpi=300)
