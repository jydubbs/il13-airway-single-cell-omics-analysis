# Integrating cell-type-resolved proteomics and scRNA-seq
## Chronic IL-13 airway epithelial remodelling

Chronic IL-13 drives airway epithelial remodelling in type 2-high asthma, including increased mucus secretion and loss of ciliary function. Although these responses have been studied extensively using transcriptomics, how closely these responses at the RNA level are reflected at the protein level has not been investigated. This study integrated low-input, cell-type-specific DIA proteomics of FACS-sorted airway epithelial cells with a published single-cell RNA-seq dataset to investigate RNA-protein concordance after chronic IL-13 treatment in airway epithelial remodelling.

This repository contains the analysis code associated with the MSc Bioinformatics research project.

---

## Overview

The main analysis workflow:
1. Process the low-input DIA-MS proteomics data and perform differential protein abundance analysis.
2. Process and analyse the transcriptomic datasets.
3. Match protein groups to gene symbols and integrate RNA and protein responses.
4. Classify genes according to RNA–protein concordance.
5. Assess whether protein-only responses could be influenced by proteomics missingness or imputation.
6. Perform Gene Ontology pathway enrichment for the different response categories.
7. Prioritise genes of potential biological interest and compare protein-only genes with airway disease-associated GWAS loci.

The primary RNA–protein integration uses the Jackson et al. (2020) scRNA-seq dataset. GSE229202 was analysed separately as an independent transcriptomic reference for comparison of selected findings.

---

## Repository structure

## Repository structure

```text

.
├── 01_proteomics_processing/
│   ├── 01_basal_proteomics_preprocessing.R
│   ├── 02_basal_proteomics_detailed_imputed_values.R
│   ├── 03_goblet_proteomics_preprocessing.R
│   ├── 04_goblet_proteomics_detailed_imputed_values.R
│   ├── 05_multiciliated_proteomics_preprocessing.R
│   ├── 06_multiciliated_proteomics_detailed_imputed_values.R
│   ├── 07_club_proteomics_preprocessing.R
│   ├── 08_club_proteomics_detailed_imputed_values.R
│   └── 09_proteomics_qc_analysis.R
│
├── 02_scrnaseq_processing/
│   ├── 01_jackson_deg_deseq2.R
│   └── 02_GSE229202_processing.ipynb
│
├── 03_rna_protein_integration/
│   ├── 01_jackson_integration.ipynb
│   ├── 02_GSE229202_integration.ipynb
│   └── 03_missingness_imputation_check.ipynb
│
├── 04_pathway_enrichment/
│   ├── 01_GO_overrepresentation_global_protein.R
│   └── 02_GO_overrepresentation.R
│
├── 05_candidate_prioritisation/
│   └── 01_gwas_correlation.ipynb
│
└── README.md

```

## Data sources

### Proteomics
Unpublished data, input data not provided.

### Jackson et al. (2020)
Publicly available scRNA-seq dataset
GEO: GSE145013

### GSE229202
Publicly available scRNA-seq dataset
GEO: GSE229202

## Requirements

### Python

- Python 3
- scanpy
- pandas
- numpy
- scipy
- matplotlib
- seaborn

### R

- R
- DESeq2
- clusterProfiler
- org.Hs.eg.db
- ggplot2
- dplyr
- writexl

