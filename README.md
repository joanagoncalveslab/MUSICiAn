# MUSICiAn: Genome-wide Identification of Genes Involved in DNA Repair via Control-Free Mutational Spectra Analysis 
[BioArXiv paper](https://doi.org/10.1101/2025.01.27.635038)


==============================

## Abstract

**Motivation:** Understanding the factors involved in DNA double-strand break (DSB) repair is crucial for the development of targeted anti-cancer therapies, yet the roles of many genes remain unclear. Recent studies show that perturbations of certain genes can alter the distribution of sequence-specific mutations left behind after DSB repair. This suggests that genome-wide screening could reveal novel DSB repair factors by identifying genes whose perturbation causes the mutational distribution spectra observed at a given DSB site to deviate significantly from the wild-type. However, designing proper controls for a genome-wide perturbation screen could be challenging. We explore the idea that a genome-wide screen might allow us to forgo the use of traditional non-targeting controls by reframing the analysis as an outlier detection problem, assuming that most genes have minimal influence on DSB repair.

**Results:** We propose MUSICiAn (Mutational Signature Catalogue Analysis), a compositional data analysis method that ranks gene perturbation-specific mutational spectra without controls by measuring deviations from the central tendency in the distributions of all spectra. We show that MUSICiAn can effectively estimate pseudo-controls for the existing Repair-seq dataset, screening 476 genes and 60 non-targeting controls. We further apply MUSICiAn to a genome-wide dataset profiling mutational outcomes induced by CRISPR-Cas9 at three target sites across cells with individual perturbations of 18,406 genes. MUSICiAn successfully recovers known genes, highlights the spliceosome as a lesser-appreciated player in DSB repair, and reveals candidates for further investigation.

This repository contains all the code used to generate process the data and generate the results.

Directory Structure
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    │
    ├── batch              <- Slurm scripts 
    |
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── exploratory    <- Ad-hoc analysis
    │   └── other          <- Miscellanious
    │   └── reports        <- Article figures, etc
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── data           <- Scripts to download or generate data    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │                  predictions
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io


--------