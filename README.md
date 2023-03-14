# End of Degree Project: Differential Gene Expression Analysis.

**NOTE**: This repository contains all of the scripts I've created for my Bachelor's End of Degree Project. The main motivation is being able to share what I've done and serve as a way for me to gain experience with the git version control system and remote repositories. Thus, **don't expect any fully functional and mantained program**.

---

## Description

dgeapy.py is a script that tries to analyse Differential Gene Expression (DGE) data at a different levels.

## Dependencies

- For the analysis:
    - [pandas](<https://pypi.org/project/pandas/>): dataframe analysis
    - [numpy](<https://pypi.org/project/numpy/>): computing
- Data visualisation
    - [matplotlib](<https://pypi.org/project/matplotlib/>): low-level manipulations.
    - [seaborn](<https://pypi.org/project/seaborn/>): high-level manipulations.
    - [matplotlib-venn](<https://pypi.org/project/matplotlib-venn/>): venn diagram generation.
    - [UpSetPlot](<https://pypi.org/project/UpSetPlot/0.8.0/>): upset plot generation.
    - [pySankey](<https://github.com/Pierre-Sassoulas/pySankey/tree/main>): sankey diagram generation.

## Usage

```
> ./dgeapy.py -h

dgeapy: Script for Differential Gene Expression data analyisis from a dataframe.

usage: dgeapy.py <command> [options]

    dataframe analyisis:
        multiplemuts    analyze a dataframe contaning 3 mutants

    utilities:
        mkconfigs       create config.json file samples

    options:
        -h, --help
```

For DGE analyisis of multiple samples you can use:

```
> ./dgeapy.py multiplemuts -h

usage: dgeapy.py multiplemuts <config.json>

Differential Gene Expression data analysis from multiple dataframes containing mutant vs. wild type experiments.

positional arguments:
  <config.json>        path to JSON configuration file

optional arguments:
  -h, --help           show this help message and exit
  --pvalue PVALUE      p-value threshold, default is 0.05
  --padj PADJ          adjusted p-value threshold, default is 0.05
  --fc FC              fold change threshold, default is 2.00
  --formats [FORMATS]  plot formats, defalut is png
  -g                   include non-coding transcripts
```

