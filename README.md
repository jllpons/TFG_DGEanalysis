# DGEAPY: Differential Gene Expression Analysis using Python

**NOTE**: This repository contains all of the scripts I've created for my Bachelor's End of Degree Project. The main motivation is being able to share what I've done and serve as a way for me to gain experience with the git version control system and remote repositories. Thus, **don't expect any fully functional and mantained program**.

---

## Description

dgeapy.py is a script that tries to analyse Differential Gene Expression (DGE) data at a different levels with a focus on data visualisation.

## Dependencies

- For the analysis:
    - [pandas](<https://pypi.org/project/pandas/>): dataframe analysis
    - [numpy](<https://pypi.org/project/numpy/>): computing
- Data visualisation
    - [matplotlib](<https://pypi.org/project/matplotlib/>): low-level manipulations.
    - [seaborn](<https://pypi.org/project/seaborn/>): high-level manipulations.
    - [matplotlib-venn](<https://pypi.org/project/matplotlib-venn/>): venn diagram generation.
    - [pyvenn](<https://pypi.org/project/venn/>): 4-set venn diagram generation.
    - [UpSetPlot](<https://pypi.org/project/UpSetPlot/0.8.0/>): upset plot generation.
    - [pySankey](<https://github.com/Pierre-Sassoulas/pySankey/tree/main>): sankey diagram generation.

## Usage

```
> ./dgeapy.py -h

dgeapy: Script for Differential Gene Expression data analyisis at different levels.

usage: dgeapy.py <command> [options]

    dataframe analyisis:
        multiplemuts    analyze a dataframe contaning 3 mutants

    utilities:
        mapgenes        map geneIDs using a <map.tsv> file
        mkconfigs       create config.json file samples

    options:
        -h, --help
```

For DGE analyisis of multiple samples you can use:

```
> ./dgeapy.py multiplemuts -h
usage: dgeapy.py multiplemuts <config.json>

Differential Gene Expression data analysis between 2, 3 and 4 'mutant vs. wild' like type
experiments.

positional arguments:
  <config.json>     path to JSON configuration file

optional arguments:
  -h, --help        show this help message and exit
  --padj FLOAT      adjusted p-value threshold, default is 0.05
  --fc FLOAT        fold change threshold, default is 1.50
  --formats [STR,]  plot formats, defalut is png
  -n, --non-coding  include non-coding transcripts
```

And if you need to map geneIDs between different strains you can use:

```
> ./dgeapy.py mapgenes -h
usage: dgeapy.py mapgenes <config.json>

Script for mapping gene IDs between different strains. Uses a TSV file that
serves as a map and a number of dataframes containing geneIDs. Returns the
dataframes containing only mapped geneIDs.

positional arguments:
  <config.json>      path to JSON configuration file

optional arguments:
  -h, --help         show this help message and exit
  -k STR, --key STR  geneIDs column name that will appear as "mapped_geneID"
                     in the generated dataframes. Otherwise strain geneIDs will be used
  -s, --stats        generate a <stats.txt> file.
  -u, --unmapped     save two TSV files containing (1) orphan geneIDs and (2) geneIDs
                     that are mapped but not present in input dfs> ./dgeapy.py mapgenes -h
```

