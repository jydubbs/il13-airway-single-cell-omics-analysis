# Integrating cell-type-resolved proteomics and scRNA-seq
## Chronic IL-13 airway epithelial remodelling

Chronic IL-13 drives airway epithelial remodelling in type 2-high asthma, including increased mucus secretion and loss of ciliary function. Although these responses have been studied extensively using transcriptomics, how closely these responses at the RNA level are reflected at the protein level has not been investigated.

This study integrated low-input, cell-type-specific DIA-MS proteomics of FACS-sorted airway epithelial cells with published single-cell RNA-seq data to investigate RNA–protein concordance after chronic IL-13 treatment in airway epithelial remodelling. The primary integration used the Jackson et al. (2020) dataset, with GSE229202 used as an second independent transcriptomic reference for selected findings (e.g. CLCA2, ORAI1, STIM1, ITPR2, ITPR3).

This repository contains the analysis code associated with the QMUL MSc Bioinformatics research project (2025-2026).

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

---


## Scripts

#### 01 — Proteomics processing

* 01_basal_proteomics_preprocessing.R : Processes basal stem cell proteomics data, including log2 transformation, protein-specific missing-value imputation, differential abundance testing and missingness classification.
* 02_basal_proteomics_detailed_imputed_values.R : Generates detailed observed and imputed replicate-level values for assessing potential imputation-driven results.
* 03_goblet_proteomics_preprocessing.R : Performs the equivalent proteomics processing and differential abundance analysis for goblet cells.
* 04_goblet_proteomics_detailed_imputed_values.R : Generates detailed observed and imputed values for goblet cell proteins.
* 05_multiciliated_proteomics_preprocessing.R : Performs proteomics processing and differential abundance analysis for multiciliated cells.
* 06_multiciliated_proteomics_detailed_imputed_values.R : Generates detailed observed and imputed values for multiciliated cell proteins.
* 07_club_proteomics_preprocessing.R : Performs proteomics processing and differential abundance analysis for club cells.
* 08_club_proteomics_detailed_imputed_values.R : Generates detailed observed and imputed values for club cell proteins.
* 09_proteomics_qc_analysis.R : Performs proteomics quality-control analyses and visualisation across cell types and treatment conditions.

#### 02 — scRNA-seq processing

* 01_jackson_deg_deseq2.R : Performs donor-level pseudobulk differential-expression analysis of the Jackson et al. (2020) scRNA-seq dataset using DESeq2.
* 02_GSE229202_processing.ipynb : Processes the independent GSE229202 scRNA-seq dataset, including quality control, cell-type identification and Wilcoxon differential-expression analysis.

#### 03 — RNA–protein integration

* 01_jackson_integration.ipynb : Main RNA–protein integration pipeline. Resolves ambiguous protein-to-gene mappings, integrates proteomics with Jackson et al. (2020) RNA measurements and differential-expression results, and assigns RNA–protein concordance categories.
* 02_GSE229202_integration.ipynb : Applies the same RNA–protein integration and concordance framework using GSE229202 as an independent transcriptomic reference.
* 03_missingness_imputation_check.ipynb : Examines proteomics missingness patterns across RNA–protein concordance categories to assess whether protein-only results may be influenced by missing value imputation.

#### 04 — Pathway enrichment

* 01_GO_overrepresentation_global_protein.R : Performs GO Biological Process over-representation analysis of the overall differential protein response.
* 02_GO_overrepresentation.R : Performs GO Biological Process over-representation analysis separately for RNA–protein concordance categories and cell types.

#### 05 — Candidate prioritisation

* 01_gwas_correlation.ipynb : Compares protein-only genes with genes associated with airway diseases (asthma, COPD and lung-function) GWAS loci.

---

## Data sources

### Proteomics
Unpublished data, input data not provided.

### Jackson et al. (2020)
Publicly available scRNA-seq dataset  
- GEO: [GSE145013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145013)

### GSE229202
Publicly available scRNA-seq dataset  
- GEO: [GSE229202](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229202)

---

## Packages

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

