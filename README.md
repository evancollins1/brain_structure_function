# Mapping the structure-function relationship along macroscale gradients in the human brain

<p align="center">
  <img src="data/misc/cover_wide.png" alt="alt text" width="800px" align="middle"/>
</p>

*Image: Yale Brain Atlas shown left. Generalized WM tractography from HCP shown right. Evan Collins / Yale*

## Description

This repository contains all data processing code and analysis discussed in the paper "Mapping the structure-function relationship along macroscale gradients in the human brain". It is currently in revision. 

_____

# Table of contents

- [Abstract](#abstract)
- [Repository structure](#repository-structure)
- [System requirements](#system-requirements)
- [Installation guide](#installation-guide)
- [Demo and Instructions](#demo-and-instructions)

_____

# Abstract

Functional coactivation between human brain regions is partly explained by white matter connections; however, how the structure-function relationship varies by function remains unclear. Here, we reference large data repositories to compute maps of structure-function correspondence across hundreds of specific functions and brain regions. We use natural language processing to accurately predict structure-function correspondence for specific functions and to identify macroscale gradients across the brain that correlate with structure-function correspondence as well as cortical thickness. Our findings suggest structure-function correspondence unfolds along a sensory-fugal organizational axis, with higher correspondence in primary sensory and motor cortex for perceptual and motor functions, and lower correspondence in association cortex for cognitive functions. Our study bridges neuroscience and natural language to describe how structure-function coupling varies by region and function in the brain, offering insight into the diversity and evolution of neural network properties.

# Repository structure

The repository is structured as follows:

- `/data`: contains all data needed to run analyses. Note due to licensing reasons, the user will need to download [LIWC2015](https://www.liwc.app/) and move it to the `/nlp` subfolder.
    - `/atlas`: contains all data pertinent to the Yale Brain Atlas (696-parcel-version) coordinates and indices.
    
    - `/atlas_690`: contains all data pertinent to the Yale Brain Atlas (690-parcel-version) coordinates and indices. This atlas version without corpus callosum parcels was used for our cortical thickness analysis.
     
    - `/brain_map`: contains brain images generated from the `plot_brain_map` function in `/scripts/utils.R`; this folder can be ignored for all practical purposes.
    
    - `/connectivity`: contains all data pertinent to functional and structural connectivities.
    
        - `/additional_sc_metrics`: contains the files generated for our data from the approach adapted from Esfahlani et al. (2022) to survey many additional structural connectivity measures.
        
    - `/cortical_thickness`: contains cortical thickness data in Yale Brain Atlas space for 200 healthy individuals from Human Connectome Project (HCP).
    
    - `/misc`: contains miscellaneous data files, including demographics and data processing files saved along the way of running chunks in `analyze_structure_function.Rmd`.
    
    - `/metaanalysis`: contains dataframe of cleaned functional term list as well as separate folders for all other files pertinent to Neurosynth or NeuroQuery.
        
        - `/neuroquery`: contains functional activation dataframe for NeuroQuery.
        
        - `/neurosynth`: contains functional activation, embeddings, abstract text, and other dataframes for Neurosynth.  
    
    - `/nlp`: contains data pertinent to natural language processing (NLP) used in this study.
    
    - `/rsfMRI`: contains rsfMRI connectivity matrices in Yale Brain Atlas space for 34 healthy subjects.
    
    - `/tractography`: contains the connectivity matrices for white matter (WM) streamline count and WM length in Yale Brain Atlas space for 1,065 healthy individuals from HCP. Initial data acquired from [Yeh (2022)](https://brain.labsolver.org/hcp_ya.html). Please acknowledge the source to Yeh (2022) and WU-Minn HCP.
    
- `/figures`: contains the figures included in the paper.

    - `/subplots`: contains the individual images used to compile the figures.

- `/scripts`: contains the R, Python, and MATLAB scripts used to process data and generate figures.

    - `/sc_metrics_functions`: contains MATLAB functions from the [GitHub folder](https://github.com/brain-networks/local_scfc/tree/main/fcn) for Esfahlani et al. (2022) to survey many stuctural connectivity measures.

# System requirements

## Hardware requirements

To run the analysis code in this repository, it is recommended to have a computer with enough RAM (> 8 GB) to support the in-memory operations. 

## Software requirements

This code has been implemented with `R` version 4.3.1. Other `R` versions will likely work. See below for specific `R` packages required. To run `.Rmd` files, RStudio is recommended.

A small subset of analysis, namely, spectral clustering and diffusion mapping, processing NeuroQuery data, creation of a structure-function network plot, processing NeuroQuery data, creating ML models for SF prediction was conducted using `Python` version 3.10.12. See below for creating an appropriate conda environment. To run `.ipynb` files, Visual Studio Code is recommended.

Lastly, structural connectivity metric optimization was carried out in `MATLAB` version 2022b using code from [Esfahlani et al. (2022)](https://www.nature.com/articles/s41467-022-29770-y). Other `MATLAB` versions will likely work.

### OS requirements

This code has been tested on the following systems, although it should work generally for other OS types (e.g., Windows) with potentially minor required changes to package versions.

- macOS: Ventura 13.2.1

- Linus: Ubuntu 22.04

# Installation guide

## Install repository

```
$ git clone https://github.com/evancollins1/brain_structure_function.git
```

## Install R package dependencies

From the `R` terminal, install the following packages:

```
install.packages(c("dplyr", "plotly", "ggplot2", "Rmisc", "png", "grid", "ggpubr", "reactable", "reactablefmtr", "readobj", "svMisc", "httr", "easyPubMed", "word2vec", "reticulate", "mclust", "tidyverse", "ggseg3d", "lsa", "corrplot", "ggbreak", "patchwork", "Rtsne", "ggrepel", "ggcharts", "igraph", "quanteda", "MASS", "factoextra", "ggtext", "rmatio", "leaps", "car", "igraph"))
```

The `.Rmd` files function will all the packages in their versions as they appear on CRAN on October 1, 2023. The specific package versions are as follows:

```
dplyr=1.1.2
plotly=4.10.2
ggplot2=3.4.3
Rmisc=1.5.1
png=0.1.8
grid=4.3.1
ggpubr=0.6.0
reactable=0.4.4
reactablefmtr=2.0.0
readobj=0.4.1
svMisc=1.2.3
httr=1.4.6
easyPubMed=2.13
word2vec=0.3.4
reticulate=1.30
mclust=6.0.0
tidyverse=2.0.0
ggseg3d=1.6.3
lsa=0.73.3
corrplot=0.92
ggbreak=0.1.2
patchwork=1.1.2
Rtsne=0.16
ggrepel=0.9.3
ggcharts=0.2.1
igraph=1.5.0
quanteda=3.3.1
MASS=7.3.60
factoextra=1.0.7
ggtext=0.1.2
rmatio=0.19.0
leaps=3.1
car=3.1.2
igraph=1.5.0
```

Moreover, one function in `utils.R` requires the use of the `orca` package to export static images. `orca` version 1.3.1 was used for this study. Installation instructions for `orca` can be viewed [here](https://github.com/plotly/orca). If Node.js is installed, one way to download `orca` is with `npm` as follows:

```
$ npm install -g electron@6.1.4 orca
```

## Install Python modules

In the `brain_structure_function` main folder, create a conda environnment (`sf_python`) with the modules specified in `environment.yml`.

```
$ conda env create -f environment.yml
$ conda activate sf_python
$ conda install jupyter
```

## Install MATLAB

Using `MATLAB` code from Esfahlani et al. (2022), additional structural connectivity metrics were computed. This script was implemented in `MATLAB` version 2022b. Instructions for downloading `MATLAB` can be found [here](https://www.mathworks.com/products/matlab.html).

# Demo and Instructions

This repository contains all the data, data processing code, and analysis code to replicate all findings and figures shown in the paper. There is one exception, due to licensing reasons, the [LIWC2015](https://www.liwc.app/) dataset is not available in this repository; the user will need to download and place it into the `/data/nlp` folder as `LIWC2015_Dictionary.dic`.

This repository also contains the processed data resulting from the data processing code and figures generated from the analysis code; thus, it is not necessary to run each script in the order we initially executed them.

Nevertheless, the following list details the precise order (and purpose) of executing the scripts in the `/scripts` folder to generate all figures included in the paper. Note that each subsequent script as shown in the list below often relies on the processed data generated from the prior scripts. Again, all processed data is already included in the repository.

1. `process_neurosynth_data.Rmd` was run to process functional activation data from Neurosynth and transform it into Yale Brain Atlas space.

2. `process_neuroquery_data.ipynb` was run to process functional activation data from NeuroQuery and transform it into Yale Brain Atlas space.

2. `compute_word_embeddings.Rmd` was run to compute word embeddings for the functional terms of Neurosynth.

3. `compute_connectivities.Rmd` was run to compute functional and structural connectivities.

4. `compute_additional_sc_metrics.m` was to compute additional structural connectivity metrics adapted from [Esfahlani et al. (2022)](https://www.nature.com/articles/s41467-022-29770-y).

4. `analyze_structure_function.Rmd` was run to conduct most analyses and generate most figures. Code chunks are labeled by the figure number they generate. Note that sections of this extensive markdown are marked for the user to refer to another script to generate processed data or a figure. As the user generates the figures in order, they will encounter the need to run additional scripts in this order:

    i. `generate_network_plot.ipynb` was used to generate Figure 2A.
  
    ii. `spectral_clustering_and_diffusion_mapping.ipynb` was used to generate spectral clustering dataframes for Figures 2B-C, S1A, S3.
    
    iii. `generate_ml_models_for_sf_by_term.ipynb` was used to create and evaluate the support vector machine and feed-forward neural network featured in Figure 4F.
    

The expected output of this demo is the generation of all figures included in the paper. Running all code should take ~6 hours.
