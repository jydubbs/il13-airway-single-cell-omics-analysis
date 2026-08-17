# Integrating cell-type-resolved proteomics and scRNA-seq
## Chronic IL-13 airway epithelial remodelling

Chronic IL-13 drives airway epithelial remodelling in type 2-high asthma, including increased mucus secretion and loss of ciliary function. Although these responses have been studied extensively using transcriptomics, how closely these responses at the RNA level are reflected at the protein level has not been investigated. This study integrated low-input, cell-type-resolved DIA proteomics of FACS-sorted airway epithelial cells with a published single-cell RNA-seq dataset to investigate RNA-protein concordance after chronic IL-13 treatment in airway epithelial remodelling.

## Overview

This repository contains the code used to integrate low-input DIA proteomics with scRNA-seq to investigate RNA-protein concordance following chronic IL-13 treatment in goblet, club and multiciliated airway epithelial cells.

## Analysis workflow

1. Proteomics preprocessing and differential abundance
2. scRNA-seq preprocessing and cell-type annotation
3. Differential expression analysis
4. RNA-protein matching by HGNC gene symbol
5. RNA-protein concordance classification
6. GO Biological Process enrichment
7. Protein missingness analysis
8. Candidate prioritisation
9. Cross-dataset validation

## Repository structure

analysis/
    ...
data/
    ...
results/
    ...

## Data

### Proteomics
Unpublished data, input data not provided.

### Jackson et al. (2020)
GEO: GSE145013

### GSE229202
GEO: GSE229202

## Requirements

Python ...
R ...

## Citation

...

## Licence

...
