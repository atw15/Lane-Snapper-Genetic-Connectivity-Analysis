# Lane Snapper Genetic Connectivity Analysis

This repository contains R code and data analysis workflow for studying the genetic connectivity of Lane Snapper (Lutjanus synagris) populations across the Gulf of Mexico, southeastern coast of the United States, and Belize. The analysis includes data filtering, population structure analysis, Hardy-Weinberg equilibrium testing, and gene flow estimation.

## Table of Contents

- [Introduction](#introduction)
- [Data Preparation and Filtering](#data-preparation-and-filtering)
- [Function Definitions](#function-definitions)
- [Data Segregation and HWE Testing](#data-segregation-and-hwe-testing)
- [Bayesian Modeling with JAGS](#bayesian-modeling-with-jags)
- [Model Comparison](#model-comparison)
- [Visualization](#visualization)
- [Gene Flow and Admixture Analysis](#gene-flow-and-admixture-analysis)
- [Acknowledgements](#acknowledgements)
- [References](#references)

## Introduction

This project aims to investigate the genetic diversity and connectivity among Lane Snapper populations using genotyping by random amplicon sequencing (GRAS-Di). We focus on analyzing population structure, Hardy-Weinberg equilibrium, and gene flow using filtered SNP datasets.

## Data Preparation and Filtering

### Initial Filtering

1. **Start with 777,138 sites.**
2. **dDocent default settings filtered to 15,392 sites.**
3. **Three-step filter applied:**
    - Variants successfully genotyped in 60% of individuals.
    - Minimum quality score of 30.
    - Minor allele count of 3.
    - **Result:** 6552 out of 15,392 sites retained.

### Filtering Paths

#### Path 1

- Remove individuals with >60% missing data.
- **Result:** 5 individuals removed, 6552 sites retained.

#### Path 2

- Apply minimum mean genotype depth filter (minDP 3).
- Remove individuals with >60% missing data.
- **Result:** 8 individuals removed, 6552 sites retained.

### Further Filtering

- Retest Hardy-Weinberg equilibrium (HWE) for filtered populations.
- **Result:** Approximately 1949 SNPs retained.

## Function Definitions

### HWE Testing and Basic Stats

- Functions to convert data formats (genind to hierfstat) and calculate basic statistics.
- Hardy-Weinberg equilibrium testing for each population.
- Combination of HWE test results with basic statistics for further filtering.

## Data Segregation and HWE Testing

- Segregate populations into individual genlight objects.
- Run HWE tests on each population.
- Combine HWE test results with basic statistics to identify loci with significant deviation from HWE.
- Retest HWE for filtered populations and check for consistency.

## Bayesian Modeling with JAGS

- Bayesian models used to re-estimate growth parameters, focusing on models preferred in preliminary evaluations.
- Truncated normal probability density function applied for residual variability in length.
- Models fit using a Gibbs Sampler implemented in JAGS.
- Deviance Information Criterion (DIC) used for model selection.

## Model Comparison

- Compare models using DIC values to determine the best-fitting model.

## Visualization

- Admixture bar plots and pie charts to visualize genetic structure and mean admixture proportions for each site.
- Basemap with pie charts for geographic visualization of genetic diversity.

## Gene Flow and Admixture Analysis

- Calculate pairwise Fst values and gene flow metrics.
- Conduct admixture analysis using LEA package and plot ancestry matrices.
- Evaluate population structure using cross-entropy and select optimal number of clusters (K).

## Acknowledgements

This project was developed under the guidance of Dr. Elizabeth Babcock. Her expertise and support were invaluable in completing this analysis.


## References

- **Enoki, H., et al. (2018).** GRAS-Di: Genotyping by random amplicon sequencing, direct. *Nature Communications*.
- **Hosoya, S., et al. (2019).** Improving RADseq by integrating random amplicon sequencing with GRAS-Di. *BMC Genomics*.
- **Ogle, D. H. (2016).** *Introductory Fisheries Analyses with R*. CRC Press.
- **Lunn, D., et al. (2012).** The BUGS project: Evolution, critique and future directions. *Statistics in Medicine*.
- **Kamvar, Z. N., et al. (2014).** Poppr: An R package for genetic analysis of populations with clonal, partially clonal, and/or sexual reproduction. *Molecular Ecology Resources*.
- **Paradis, E. (2010).** pegas: An R package for population genetics with an integrated-modular approach. *Bioinformatics*.
- **Pritchard, J. K., et al. (2000).** Inference of population structure using multilocus genotype data. *Genetics*.
- **Pearman, W. S., et al. (2022).** Revisiting the role of Hardy-Weinberg equilibrium in population genomics and its applications. *Molecular Ecology*.
- **Hudson, R. R., et al. (1992).** Gene flow and the geographic structure of natural populations. *Science*.
- **R Core Team (2022).** R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.
- **Frichot, E., et al. (2015).** LEA: An R package for landscape and ecological association studies. *Methods in Ecology and Evolution*.
- **Glaubitz, J. C., et al. (2014).** TASSEL: Software to analyze diversity among SNPs in populations. *Molecular Ecology Resources*.
- **Gruber, B., et al. (2018).** dartR: Importing and analysing SNP and silicodart data generated by DArT P/L. *Molecular Ecology Resources*.
- **Jombart, T., et al. (2011).** adegenet 1.3-1: New tools for the analysis of genome-wide SNP data. *Bioinformatics*.

---

This repository provides a comprehensive workflow for analyzing genetic connectivity of Lane Snapper, from data filtering and HWE testing to advanced modeling and visualization. The annotated code offers clear explanations for each step, making it easy to follow and reproduce the analysis.
