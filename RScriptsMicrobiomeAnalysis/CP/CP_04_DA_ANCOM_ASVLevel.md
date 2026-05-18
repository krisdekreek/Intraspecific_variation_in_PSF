---
title: "CP_04_differential_abundance_ANCOM - Less filtered data and outliers removed - ASV level"
author: "Kris de Kreek"
date: "2025-12-09"
output: 
  html_document:
    toc: true
    keep_md: true
editor_options: 
  chunk_output_type: console
---

This script is adapted from a script form Melissa Uribe Acosta. Extra input comes from a script of Pedro Beschore da Costa.   
   
   
Tutorials   
- [Heatmap for final plot](http://rstudio-pubs-static.s3.amazonaws.com/288398_185f2889a5f641c6b9aa7b14fa15b634.html)
- [Github page Ancom with FAQ including an explanation about what intercepts are](https://github.com/FrederickHuangLin/ANCOMBC)   
- [Info on Ancom](https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html)   
- [And more info on Ancom](https://microbiome.github.io/course_2022_turku/differential-abundance-analysis.html#ancom-bc)   
   
   
# 4.0 load libraries, improve memory use and subset objects for the desired treatment comparisons
### Load libraries

``` r
R.version$version.string # prints R version
```

```
## [1] "R version 4.5.1 (2025-06-13 ucrt)"
```

``` r
library(phyloseq)
packageVersion("phyloseq")
```

```
## [1] '1.52.0'
```

``` r
library(metamisc) #for phyloseq_sep_variable (and more?)
packageVersion("metamisc")
```

```
## [1] '0.4.0'
```

``` r
library(DESeq2)
packageVersion("DESeq2")
```

```
## [1] '1.48.2'
```

``` r
library(zinbwave)
packageVersion("zinbwave")
```

```
## [1] '1.30.0'
```

``` r
library(ANCOMBC)
packageVersion("ANCOMBC")
```

```
## [1] '2.10.1'
```

``` r
#library(rmeta) #rusher321/rmeta https://github.com/rusher321/rmeta/
#packageVersion("rmeta") 
library(dplyr)
packageVersion("dplyr")
```

```
## [1] '1.1.4'
```

# 4.2 ANCOM-BC ASV level
## 4.2.1 ANCOM including structural zeroes detection
### Accession
#### Run ANCOM

``` r
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_unnormalized_bac_ps_ForDA_Acc.RData")

ancomWZ_Bac_Acc <- lapply(CP_unnormalized_bac_ps_ForDA_Acc, function (x)
                           ancombc(data = x,
                             #tax_level = "Family",
                             formula = "Cat_treatment",
                             group = "Cat_treatment",
                             p_adj_method = "fdr",
                             prv_cut = 0.1, #default, taxa that are in less than the fraction present are discarded
                             lib_cut = 0, # threshold for filtering samples based on library size, 0 has no filtering
                             struc_zero = TRUE,
                             neg_lb = TRUE,
                             tol = 1e-5,
                             max_iter = 100,
                             conserve = TRUE,
                             alpha = 0.05,
                             global = FALSE))
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 11, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## Loading required package: foreach
```

```
## Loading required package: rngtools
```

```
## ANCOM-BC primary results ...
```

```
## Merge the information of structural zeros ... 
## Note that taxa with structural zeros will have 0 p/q-values and SEs
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 11, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## Merge the information of structural zeros ... 
## Note that taxa with structural zeros will have 0 p/q-values and SEs
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 11, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## Merge the information of structural zeros ... 
## Note that taxa with structural zeros will have 0 p/q-values and SEs
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 11, Mb = 11
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## Merge the information of structural zeros ... 
## Note that taxa with structural zeros will have 0 p/q-values and SEs
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 12, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## Merge the information of structural zeros ... 
## Note that taxa with structural zeros will have 0 p/q-values and SEs
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 11, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## Merge the information of structural zeros ... 
## Note that taxa with structural zeros will have 0 p/q-values and SEs
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 10, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## Merge the information of structural zeros ... 
## Note that taxa with structural zeros will have 0 p/q-values and SEs
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 12, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## Merge the information of structural zeros ... 
## Note that taxa with structural zeros will have 0 p/q-values and SEs
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 12, Mb = 11
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## Merge the information of structural zeros ... 
## Note that taxa with structural zeros will have 0 p/q-values and SEs
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 12, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## Merge the information of structural zeros ... 
## Note that taxa with structural zeros will have 0 p/q-values and SEs
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 12, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## Merge the information of structural zeros ... 
## Note that taxa with structural zeros will have 0 p/q-values and SEs
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 10, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## Merge the information of structural zeros ... 
## Note that taxa with structural zeros will have 0 p/q-values and SEs
```

``` r
save(ancomWZ_Bac_Acc, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/ancomWZ_Bac_Acc.RData")
```

#### Organize and filter results

``` r
ancomWZ_res_df_Acc_Bac <- lapply(ancomWZ_Bac_Acc, function (x)
  data.frame(
  Species = x$res$lfc$taxon, # ASVs
  lfc = unlist(x$res$lfc["Cat_treatmentMb"]),
  se = unlist(x$res$se["Cat_treatmentMb"]),
  W = unlist(x$res$W["Cat_treatmentMb"]), # test statistic
  p_val = unlist(x$res$p_val["Cat_treatmentMb"]),
  q_val = unlist(x$res$q_val["Cat_treatmentMb"]), # adjusted p-values
  diff_abn = unlist(x$res$diff_abn["Cat_treatmentMb"]))) # true if q value is less than alpha

fdr_ancomWZ_Acc_Bac <- lapply(ancomWZ_res_df_Acc_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.05)))
fdr_ancomWZ_Acc_Bac_0.001 <- lapply(ancomWZ_res_df_Acc_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.001)))
fdr_ancomWZ_Acc_Bac_0.0001 <- lapply(ancomWZ_res_df_Acc_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.0001)))

head(fdr_ancomWZ_Acc_Bac$VL)
```

```
##                     Species        lfc         se          W        p_val
## Cat_treatmentMb68   bASV_72 -0.4323658 0.08996851 -4.8057455 1.541759e-06
## Cat_treatmentMb131 bASV_138  0.3929632 0.00000000  0.7046867 0.000000e+00
## Cat_treatmentMb142 bASV_149  0.2991132 0.00000000  0.6792483 0.000000e+00
## Cat_treatmentMb158 bASV_165  0.5512705 0.00000000  1.0276822 0.000000e+00
## Cat_treatmentMb176 bASV_184  1.4275005 0.00000000  2.1869220 0.000000e+00
## Cat_treatmentMb190 bASV_198  0.2178036 0.00000000  0.8032644 0.000000e+00
##                          q_val diff_abn
## Cat_treatmentMb68  0.003312785     TRUE
## Cat_treatmentMb131 0.000000000     TRUE
## Cat_treatmentMb142 0.000000000     TRUE
## Cat_treatmentMb158 0.000000000     TRUE
## Cat_treatmentMb176 0.000000000     TRUE
## Cat_treatmentMb190 0.000000000     TRUE
```

``` r
# Number of differential ASV at different alphas
length(fdr_ancomWZ_Acc_Bac$VL$q_val)
```

```
## [1] 1118
```

``` r
length(fdr_ancomWZ_Acc_Bac_0.001$VL$q_val)
```

```
## [1] 1112
```

``` r
length(fdr_ancomWZ_Acc_Bac_0.0001$VL$q_val)
```

```
## [1] 1112
```

``` r
save(fdr_ancomWZ_Acc_Bac, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/fdr_ancomWZ_Acc_Bac.RData")
```

``` r
# Acc_Bac_ancomNoZ_ASVs <- lapply(fdr_ancomWZ_Acc_Bac, function(x) as.vector(x[,'Species']))
# fdr_ancomWZ_Acc_Bac$VL[,'Species']
# fdr_ancomWZ_Acc_Bac$VL$Species
# 
# head(ancomWZ_Bac_Acc$VL$res$diff_abn, 20)
# summary(ancomWZ_Bac_Acc$VL$res$diff_abn)
# summary(ancomWZ_Bac_Acc$VL$res$q_val)
# length(ancomWZ_Bac_Acc$VL$res$q_val$`(Intercept)`)
# 
# head(fdr_ancomWZ_Acc_Bac$VL, 20)
# summary(fdr_ancomWZ_Acc_Bac$VL)
# 
# length(ancomWZ_Bac_Acc$VL$res$diff_abn)
# length(ancomWZ_res_df_Acc_Bac$VL$q_val)

# Clean environment
rm(list = ls())
```

#### Run ANCOM

``` r
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_unnormalized_bac_ps_ForDA_Dom.RData")

ancomWZ_Bac_Dom <- lapply(CP_unnormalized_bac_ps_ForDA_Dom, function (x)
                           ancombc(data = x,
                             #tax_level = "Family",
                             formula = "Accession + Cat_treatment",
                             group = "Cat_treatment",
                             p_adj_method = "fdr",
                             prv_cut = 0.1, #default, taxa that are in less than the fraction present are discarded
                             lib_cut = 0, # threshold for filtering samples based on library size, 0 has no filtering
                             struc_zero = TRUE,
                             neg_lb = TRUE,
                             tol = 1e-5,
                             max_iter = 100,
                             conserve = TRUE,
                             alpha = 0.05,
                             global = FALSE))
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Accession, Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 56, Mb = 59
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## Merge the information of structural zeros ... 
## Note that taxa with structural zeros will have 0 p/q-values and SEs
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Accession, Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 79, Mb = 83
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## Merge the information of structural zeros ... 
## Note that taxa with structural zeros will have 0 p/q-values and SEs
```

``` r
save(ancomWZ_Bac_Dom, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/ancomWZ_Bac_Dom.RData")
```

#### Organize and filter results

``` r
ancomWZ_res_df_Dom_Bac <- lapply(ancomWZ_Bac_Dom, function (x)
  data.frame(
  Species = x$res$lfc$taxon, # ASVs
  lfc = unlist(x$res$lfc["Cat_treatmentMb"]),
  se = unlist(x$res$se["Cat_treatmentMb"]),
  W = unlist(x$res$W["Cat_treatmentMb"]), # test statistic
  p_val = unlist(x$res$p_val["Cat_treatmentMb"]),
  q_val = unlist(x$res$q_val["Cat_treatmentMb"]), # adjusted p-values
  diff_abn = unlist(x$res$diff_abn["Cat_treatmentMb"]))) # true if q value is less than alpha

fdr_ancomWZ_Dom_Bac <- lapply(ancomWZ_res_df_Dom_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.05)))
fdr_ancomWZ_Dom_Bac_0.001 <- lapply(ancomWZ_res_df_Dom_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.001)))
fdr_ancomWZ_Dom_Bac_0.0001 <- lapply(ancomWZ_res_df_Dom_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.0001)))

head(fdr_ancomWZ_Dom_Bac$VL)
```

```
## NULL
```

``` r
# Number of differential ASV at different alphas
length(fdr_ancomWZ_Dom_Bac$VL$q_val)
```

```
## [1] 0
```

``` r
length(fdr_ancomWZ_Dom_Bac_0.001$VL$q_val)
```

```
## [1] 0
```

``` r
length(fdr_ancomWZ_Dom_Bac_0.0001$VL$q_val)
```

```
## [1] 0
```

``` r
save(fdr_ancomWZ_Dom_Bac, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/fdr_ancomWZ_Dom_Bac.RData")
```

``` r
rm(list = ls())
```

### All data
#### Run Ancom

``` r
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_unnormalized_bac_ps_ForDA.RData")

ancomWZ_Bac_All <- ancombc(data = CP_unnormalized_bac_ps_ForDA,
                           #tax_level = "Family",
                           formula = "Accession + Cat_treatment",
                           group = "Cat_treatment",
                           p_adj_method = "fdr",
                           prv_cut = 0.1, #default, taxa that are in less than the fraction present are discarded
                           lib_cut = 0, # threshold for filtering samples based on library size, 0 has no filtering
                           struc_zero = TRUE,
                           neg_lb = TRUE,
                           tol = 1e-5,
                           max_iter = 100,
                           conserve = TRUE,
                           alpha = 0.05,
                           global = FALSE)
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Accession, Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 135, Mb = 142
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## Merge the information of structural zeros ... 
## Note that taxa with structural zeros will have 0 p/q-values and SEs
```

``` r
save(ancomWZ_Bac_All, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/ancomWZ_Bac_All.RData")
```

#### Organize and filter results

``` r
ancomWZ_res_df_All_Bac <- data.frame(
    Species = ancomWZ_Bac_All$res$lfc$taxon, # ASVs
    lfc = unlist(ancomWZ_Bac_All$res$lfc["Cat_treatmentMb"]),
    se = unlist(ancomWZ_Bac_All$res$se["Cat_treatmentMb"]),
    W = unlist(ancomWZ_Bac_All$res$W["Cat_treatmentMb"]), # test statistic
    p_val = unlist(ancomWZ_Bac_All$res$p_val["Cat_treatmentMb"]),
    q_val = unlist(ancomWZ_Bac_All$res$q_val["Cat_treatmentMb"]), # adjusted p-values
    diff_abn = unlist(ancomWZ_Bac_All$res$diff_abn["Cat_treatmentMb"]))

fdr_ancomWZ_All_Bac <- ancomWZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.05)
fdr_ancomWZ_All_Bac_0.001 <- ancomWZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.001)
fdr_ancomWZ_All_Bac_0.0001 <- ancomWZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.0001)

head(fdr_ancomWZ_All_Bac)
```

```
##                     Species        lfc        se          W        p_val
## Cat_treatmentMb259 bASV_269 -1.5634104 0.1420518 -11.005918 3.578511e-28
## Cat_treatmentMb277 bASV_287 -1.0437351 0.1614311  -6.465514 1.009549e-10
## Cat_treatmentMb290 bASV_300 -0.5573383 0.1464165  -3.806528 1.409315e-04
## Cat_treatmentMb346 bASV_359  0.4800426 0.1388645   3.456914 5.463991e-04
## Cat_treatmentMb352 bASV_365 -0.5502962 0.1373861  -4.005473 6.189348e-05
## Cat_treatmentMb373 bASV_387  0.7202894 0.1347347   5.345982 8.992820e-08
##                           q_val diff_abn
## Cat_treatmentMb259 1.184487e-24     TRUE
## Cat_treatmentMb277 8.354015e-08     TRUE
## Cat_treatmentMb290 1.608563e-02     TRUE
## Cat_treatmentMb346 4.888057e-02     TRUE
## Cat_treatmentMb352 8.907279e-03     TRUE
## Cat_treatmentMb373 2.480520e-05     TRUE
```

``` r
# Number of differential ASV at different alphas
length(fdr_ancomWZ_All_Bac$q_val)
```

```
## [1] 37
```

``` r
length(fdr_ancomWZ_All_Bac_0.001$q_val)
```

```
## [1] 19
```

``` r
length(fdr_ancomWZ_All_Bac_0.0001$q_val)
```

```
## [1] 14
```

``` r
save(fdr_ancomWZ_All_Bac, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/fdr_ancomWZ_All_Bac.RData")
```

``` r
# Clean environment
rm(list = ls())
```


## 4.2.1 ANCOM WITHOUT structural zeroes detection
### Accession
#### Run ANCOM

``` r
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_unnormalized_bac_ps_ForDA_Acc.RData")

ancomNoZ_Bac_Acc <- lapply(CP_unnormalized_bac_ps_ForDA_Acc, function (x)
                           ancombc(data = x,
                             #tax_level = "Family",
                             formula = "Cat_treatment",
                             group = "Cat_treatment",
                             p_adj_method = "fdr",
                             prv_cut = 0.1, #default, taxa that are in less than the fraction present are discarded
                             lib_cut = 0, # threshold for filtering samples based on library size, 0 has no filtering
                             struc_zero = FALSE,
                             neg_lb = TRUE,
                             tol = 1e-5,
                             max_iter = 100,
                             conserve = TRUE,
                             alpha = 0.05,
                             global = FALSE))
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 11, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 11, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 11, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 11, Mb = 11
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 12, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 11, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 10, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 12, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 12, Mb = 11
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 12, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 12, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 10, Mb = 12
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

``` r
save(ancomNoZ_Bac_Acc, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/ancomNoZ_Bac_Acc.RData")
```

#### Organize and filter results

``` r
ancomNoZ_res_df_Acc_Bac <- lapply(ancomNoZ_Bac_Acc, function (x)
    data.frame(
      Species = x$res$lfc$taxon, # ASVs
      lfc = unlist(x$res$lfc["Cat_treatmentMb"]),
      se = unlist(x$res$se["Cat_treatmentMb"]),
      W = unlist(x$res$W["Cat_treatmentMb"]), # test statistic
      p_val = unlist(x$res$p_val["Cat_treatmentMb"]),
      q_val = unlist(x$res$q_val["Cat_treatmentMb"]), # adjusted p-values
      diff_abn = unlist(x$res$diff_abn["Cat_treatmentMb"]))) # true if q value is less than alpha

fdr_ancomNoZ_Acc_Bac <- lapply(ancomNoZ_res_df_Acc_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.05)))
fdr_ancomNoZ_Acc_Bac_0.001 <- lapply(ancomNoZ_res_df_Acc_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.001)))
fdr_ancomNoZ_Acc_Bac_0.0001 <- lapply(ancomNoZ_res_df_Acc_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.0001)))

head(fdr_ancomNoZ_Acc_Bac$VL)
```

```
##                     Species        lfc         se         W        p_val
## Cat_treatmentMb68   bASV_72 -0.4323658 0.08996851 -4.805745 1.541759e-06
## Cat_treatmentMb304 bASV_315  1.7483926 0.39669991  4.407343 1.046464e-05
## Cat_treatmentMb550 bASV_581  1.5402534 0.33306479  4.624486 3.755286e-06
## Cat_treatmentMb568 bASV_601 -1.5633761 0.36292834 -4.307672 1.649821e-05
## Cat_treatmentMb741 bASV_792  1.4523420 0.34746647  4.179805 2.917594e-05
## Cat_treatmentMb924 bASV_996  1.4751350 0.36342576  4.058972 4.928920e-05
##                          q_val diff_abn
## Cat_treatmentMb68  0.003312785     TRUE
## Cat_treatmentMb304 0.009130400     TRUE
## Cat_treatmentMb550 0.004368649     TRUE
## Cat_treatmentMb568 0.011515749     TRUE
## Cat_treatmentMb741 0.016970675     TRUE
## Cat_treatmentMb924 0.024574187     TRUE
```

``` r
# Number of differential ASV at different alphas
length(fdr_ancomNoZ_Acc_Bac$VL$q_val)
```

```
## [1] 7
```

``` r
length(fdr_ancomNoZ_Acc_Bac_0.001$VL$q_val)
```

```
## [1] 0
```

``` r
length(fdr_ancomNoZ_Acc_Bac_0.0001$VL$q_val)
```

```
## [1] 0
```

``` r
save(fdr_ancomNoZ_Acc_Bac, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/fdr_ancomNoZ_Acc_Bac.RData")

# Clean environment
rm(list = ls())
```

### Domestication
#### Run ANCOM

``` r
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_unnormalized_bac_ps_ForDA_Dom.RData")

ancomNoZ_Bac_Dom <- lapply(CP_unnormalized_bac_ps_ForDA_Dom, function (x)
                           ancombc(data = x,
                             #tax_level = "Family",
                             formula = "Accession + Cat_treatment",
                             group = "Cat_treatment",
                             p_adj_method = "fdr",
                             prv_cut = 0.1, #default, taxa that are in less than the fraction present are discarded
                             lib_cut = 0, # threshold for filtering samples based on library size, 0 has no filtering
                             struc_zero = FALSE,
                             neg_lb = TRUE,
                             tol = 1e-5,
                             max_iter = 100,
                             conserve = TRUE,
                             alpha = 0.05,
                             global = FALSE))
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Accession, Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 56, Mb = 59
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Accession, Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 79, Mb = 83
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

``` r
save(ancomNoZ_Bac_Dom, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/ancomNoZ_Bac_Dom.RData")
```

#### Organize and filter results

``` r
ancomNoZ_res_df_Dom_Bac <- lapply(ancomNoZ_Bac_Dom, function (x)
    data.frame(
      Species = x$res$lfc$taxon, # ASVs
      lfc = unlist(x$res$lfc["Cat_treatmentMb"]),
      se = unlist(x$res$se["Cat_treatmentMb"]),
      W = unlist(x$res$W["Cat_treatmentMb"]), # test statistic
      p_val = unlist(x$res$p_val["Cat_treatmentMb"]),
      q_val = unlist(x$res$q_val["Cat_treatmentMb"]), # adjusted p-values
      diff_abn = unlist(x$res$diff_abn["Cat_treatmentMb"]))) # true if q value is less than alpha

fdr_ancomNoZ_Dom_Bac <- lapply(ancomNoZ_res_df_Dom_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.05)))
fdr_ancomNoZ_Dom_Bac_0.001 <- lapply(ancomNoZ_res_df_Dom_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.001)))
fdr_ancomNoZ_Dom_Bac_0.0001 <- lapply(ancomNoZ_res_df_Dom_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.0001)))

head(fdr_ancomNoZ_Dom_Bac$VL)
```

```
## NULL
```

``` r
# Number of differential ASV at different alphas
length(fdr_ancomNoZ_Dom_Bac$VL$q_val)
```

```
## [1] 0
```

``` r
length(fdr_ancomNoZ_Dom_Bac_0.001$VL$q_val)
```

```
## [1] 0
```

``` r
length(fdr_ancomNoZ_Dom_Bac_0.0001$VL$q_val)
```

```
## [1] 0
```

``` r
save(fdr_ancomNoZ_Dom_Bac, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/fdr_ancomNoZ_Dom_Bac.RData")

# Clean environment
rm(list = ls())
```


### All data
#### Run Ancom

``` r
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_unnormalized_bac_ps_ForDA.RData")

ancomNoZ_Bac_All <- ancombc(data = CP_unnormalized_bac_ps_ForDA,
                           #tax_level = "Family",
                           formula = "Accession + Cat_treatment",
                           group = "Cat_treatment",
                           p_adj_method = "fdr",
                           prv_cut = 0.1, #default, taxa that are in less than the fraction present are discarded
                           lib_cut = 0, # threshold for filtering samples based on library size, 0 has no filtering
                           struc_zero = FALSE,
                           neg_lb = TRUE,
                           tol = 1e-5,
                           max_iter = 100,
                           conserve = TRUE,
                           alpha = 0.05,
                           global = FALSE)
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: Accession, Cat_treatment
```

```
## The available variables in the sample metadata: Phase_PSF, SampleNr, Accession, Domestication, Cat_treatment, Soil_conditioning, TreatmentNr, PlantNr, Batch, Table, Block_n, Rows, Cols, plots, InductionNr_CP, InductionNr_FP, Caterpillar_biomass_Sample, Leaf_Damage_Sample, GSL_leaf_analyis, C_N_Determination, Leave_samples_left_over_after_GSL_analysis, Rhizosphere_Sample, Rhizoplane_Sample, Root_GSL_Sample, Root_samples_left_over_after_GSL_analysis, Root_biomass, Bulk_Soil_Sample, Top_Soil_Sample, Caterpillar_Microbe_Sample, Phylosphere_Sample, Glycerol_Stock, DNA_rhizosphere, DNA_rhizoplane, DNA_top_soil, DNA_Caterpillar, Note, L_GSL_Leaf_weight_mg, L_Total_GSL, L_Total_Aliphatic, L_Total_Indol, L_Progoitrin, L_Glucoraphanin, L_Sinigrin, L_Gluconapin, L_Glucoerucin, L_Glucobrassicin, L_Methoxy_glucobrassicin, L_Gluconasturtiin, L_Neoglucobrassicin, L_X3MSOP, L_X5MSOP, L_X4.Pent, L_X4OHI3M, L_Note_GSL, R_GSL_Leaf_weight_mg, R_Total_GSL, R_Total_Aliphatic, R_Total_Indol, R_Progoitrin, R_Glucoraphanin, R_Sinigrin, R_Gluconapin, R_Glucoerucin, R_Glucobrassicin, R_Methoxy_glucobrassicin, R_Neoglucobrassicin, R_Gluconasturtiin, R_X3MSOP, R_X5MSOP, R_X4.Pent, R_X4OHI3M, R_Note_GSL, L_N_content, L_C_content, L_C_N_ratio, L_Note_C_N, Mean_Weight_Cat, Note_Cat_Weight, InductionNr, Survival, Dead, Survival_fraction, Root_Weight, Note_RootBiomass, library_sizes_prefiltering, Mitochondria_reads, Plastid_reads, Host_DNA_n_reads, Host_DNA_contamination_pct, library_size, is.neg
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 2
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## The sample size per group is: Co = 135, Mb = 142
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## ANCOM-BC primary results ...
```

``` r
save(ancomNoZ_Bac_All, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/ancomNoZ_Bac_All.RData")
```

#### Organize and filter results

``` r
ancomNoZ_res_df_All_Bac <- data.frame(
    Species = ancomNoZ_Bac_All$res$lfc$taxon, # ASVs
    lfc = unlist(ancomNoZ_Bac_All$res$lfc["Cat_treatmentMb"]),
    se = unlist(ancomNoZ_Bac_All$res$se["Cat_treatmentMb"]),
    W = unlist(ancomNoZ_Bac_All$res$W["Cat_treatmentMb"]), # test statistic
    p_val = unlist(ancomNoZ_Bac_All$res$p_val["Cat_treatmentMb"]),
    q_val = unlist(ancomNoZ_Bac_All$res$q_val["Cat_treatmentMb"]), # adjusted p-values
    diff_abn = unlist(ancomNoZ_Bac_All$res$diff_abn["Cat_treatmentMb"]))

fdr_ancomNoZ_All_Bac <- ancomNoZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.05)
fdr_ancomNoZ_All_Bac_0.001 <- ancomNoZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.001)
fdr_ancomNoZ_All_Bac_0.0001 <- ancomNoZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.0001)

head(fdr_ancomNoZ_All_Bac)
```

```
##                     Species        lfc        se          W        p_val
## Cat_treatmentMb259 bASV_269 -1.5634104 0.1420518 -11.005918 3.578511e-28
## Cat_treatmentMb277 bASV_287 -1.0437351 0.1614311  -6.465514 1.009549e-10
## Cat_treatmentMb290 bASV_300 -0.5573383 0.1464165  -3.806528 1.409315e-04
## Cat_treatmentMb346 bASV_359  0.4800426 0.1388645   3.456914 5.463991e-04
## Cat_treatmentMb352 bASV_365 -0.5502962 0.1373861  -4.005473 6.189348e-05
## Cat_treatmentMb373 bASV_387  0.7202894 0.1347347   5.345982 8.992820e-08
##                           q_val diff_abn
## Cat_treatmentMb259 1.184487e-24     TRUE
## Cat_treatmentMb277 8.354015e-08     TRUE
## Cat_treatmentMb290 1.608563e-02     TRUE
## Cat_treatmentMb346 4.888057e-02     TRUE
## Cat_treatmentMb352 8.907279e-03     TRUE
## Cat_treatmentMb373 2.480520e-05     TRUE
```

``` r
# Number of differential ASV at different alphas
length(fdr_ancomNoZ_All_Bac$q_val)
```

```
## [1] 37
```

``` r
length(fdr_ancomNoZ_All_Bac_0.001$q_val)
```

```
## [1] 19
```

``` r
length(fdr_ancomNoZ_All_Bac_0.0001$q_val)
```

```
## [1] 14
```

``` r
save(fdr_ancomNoZ_All_Bac, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/fdr_ancomNoZ_All_Bac.RData")

# Clean environment
rm(list = ls())
```
