# MUSICiAn: Genome-wide Identification of Genes Involved in DNA Repair via Control-Free Mutational Spectra Analysis 
[BioArXiv paper](https://doi.org/10.1101/2025.01.27.635038)


==============================

## Abstract

**Motivation:** Understanding the factors involved in DNA double-strand break (DSB) repair is crucial for the development of targeted anti-cancer therapies, yet the roles of many genes remain unclear. Recent studies show that perturbations of certain genes can alter the distribution of sequence-specific mutations left behind after DSB repair. This suggests that genome-wide screening could reveal novel DSB repair factors by identifying genes whose perturbation causes the mutational distribution spectra observed at a given DSB site to deviate significantly from the wild-type. However, designing proper controls for a genome-wide perturbation screen could be challenging. We explore the idea that a genome-wide screen might allow us to forgo the use of traditional non-targeting controls by reframing the analysis as an outlier detection problem, assuming that most genes have minimal influence on DSB repair.

**Results:** We propose MUSICiAn (Mutational Signature Catalogue Analysis), a compositional data analysis method that ranks gene perturbation-specific mutational spectra without controls by measuring deviations from the central tendency in the distributions of all spectra. We show that MUSICiAn can effectively estimate pseudo-controls for the existing Repair-seq dataset, screening 476 genes and 60 non-targeting controls. We further apply MUSICiAn to a genome-wide dataset profiling mutational outcomes induced by CRISPR-Cas9 at three target sites across cells with individual perturbations of 18,406 genes. MUSICiAn successfully recovers known genes, highlights the spliceosome as a lesser-appreciated player in DSB repair, and reveals candidates for further investigation.

## Repository details

This repository contains all the code and processed data files used to process the data and generate the results for the paper "[MUSICiAn: Genome-wide Identification of Genes Involved in DNA Repair via Control-Free Mutational Spectra Analysis]((https://doi.org/10.1101/2025.01.27.635038))". The repository assumes the raw FASTQ files have been downloaded and processed via [SIQ](https://pubmed.ncbi.nlm.nih.gov/36071722/) using the parameters `m 2 -c -e 0.05`. See the main article for full details. 

## Running the notebooks

After cloning the repository, please run `pip install -r requirements.txt` from within the top-level directory to install all the necessary packages to run python scripts and jupyter notebooks. The folder `src/data/` contains python scripts for downloading and processing the [SIQ](https://pubmed.ncbi.nlm.nih.gov/36071722/) output into mutational spectra. The Jupyter notebooks under `notebooks/reports` score the mutational spectra and generate the figures used in the article. Please refer to the main article for full details on the algorithm. 


Directory Structure
------------

    ├── LICENSE
    ├── README.md          <- The top-level README for developers using this project.
    │
    ├── batch              <- Slurm scripts 
    │    
    ├── notebooks          <- Jupyter notebooks.
    │   └── exploratory    <- Ad-hoc analysis
    │   └── other          <- Miscellanious
    │   └── reports        <- Article figures, etc
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    └── src                <- Source code for use in this project.
        ├── __init__.py    <- Makes src a Python module
        │
        ├── data           <- Scripts to download, generate and pre-process the SIQ output through to aggregated mutational spectra.
        ├── models         <- Scripts to train models and then use trained models to make
        │                  predictions
        └── visualization  <- Scripts to create exploratory and results oriented visualizations



--------