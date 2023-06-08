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
        multiplemuts        analyze a dataframe contaning 3 mutants

    utilities:
        assert-function     assign function to each gene based on a preestablished list of GO codes and KEGG pathways to define each function.
        dropNaN-in-column   drop all NaN values in a specific column
        go2ancestors        obtain all of the GO ancestors from GO codes and store them in a new column
        joindfs             perform a left join on two tables using common column
        mapgenes            map geneIDs using a <map.tsv> file
        mkconfigs           create config.json file samples

    options:
        -h, --help
```

The main DGE analyisis is done with the `multiplemuts`  command:

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

Although it is very WIP, you can check how more advanced plots can be done with specific scripts found in the `advanced_plots` directory.
