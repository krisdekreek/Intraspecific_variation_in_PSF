---
title: "FP_05_differential_abundance_ANCOM - Less filtered data and outliers removed - ASV level"
author: "Kris de Kreek"
date: "2026-02-13"
output: 
  html_document:
    toc: true
    keep_md: true
editor_options: 
  chunk_output_type: console
---

This script is adapted from a script form Melissa Uribe Acosta. Extra input comes from a script of Pedro Beschore da Costa.   
   
   
Tutorials   
- [Heatmap for final plot](http://rstudio-pubs-static.s3.amazonaws.com/288398_185f2889a5f641c6b9aa7b14fa15b635.html)
- [Github page Ancom with FAQ including an explanation about what intercepts are](https://github.com/FrederickHuangLin/ANCOMBC)   
- [Info on Ancom](https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html)   
- [And more info on Ancom](https://microbiome.github.io/course_2022_turku/differential-abundance-analysis.html#ancom-bc)   
   
   
# 5.0 load libraries, improve memory use and subset objects for the desired treatment comparisons
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

# 5.2 ANCOM-BC ASV level
## 5.2.1 ANCOM including structural zeroes detection
### Accession
#### Run ANCOM

``` r
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_unnormalized_bac_ps_ForDA_clean_Acc.RData")

# Remove Batch 1
FP_unnormalized_bac_ps_ForDA_clean_Acc <- lapply(FP_unnormalized_bac_ps_ForDA_clean_Acc, function(x) subset_samples(x, Batch != 1))

ancomWZ_Bac_Acc <- lapply(FP_unnormalized_bac_ps_ForDA_clean_Acc, function (x)
                           ancombc(data = x,
                             #tax_level = "Family",
                             formula = "Soil_conditioning",
                             group = "Soil_conditioning",
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
## The specified variables in the formula: Soil_conditioning
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
## The sample size per group is: Co = 6, Mb = 5
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
## The specified variables in the formula: Soil_conditioning
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
## The sample size per group is: Co = 6, Mb = 6
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
## The specified variables in the formula: Soil_conditioning
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
## The sample size per group is: Co = 4, Mb = 3
```

```
## Warning: Small sample size detected for the following group(s): 
## Co, Mb
## Variance estimation would be unstable when the sample size is < 5 per group
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
## The specified variables in the formula: Soil_conditioning
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
## The sample size per group is: Co = 5, Mb = 6
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
## The specified variables in the formula: Soil_conditioning
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
## The sample size per group is: Co = 6, Mb = 6
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
save(ancomWZ_Bac_Acc, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ASVLevel/ancomWZ_Bac_Acc.RData")
```

#### Organize and filter results

``` r
ancomWZ_res_df_Acc_Bac <- lapply(ancomWZ_Bac_Acc, function (x)
  data.frame(
  Species = x$res$lfc$taxon, # ASVs
  lfc = unlist(x$res$lfc[3]),
  se = unlist(x$res$se[3]),
  W = unlist(x$res$W[3]), # test statistic
  p_val = unlist(x$res$p_val[3]),
  q_val = unlist(x$res$q_val[3]), # adjusted p-values
  diff_abn = unlist(x$res$diff_abn[3]))) # true if q value is less than alpha

fdr_ancomWZ_Acc_Bac <- lapply(ancomWZ_res_df_Acc_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.05)))
fdr_ancomWZ_Acc_Bac_0.001 <- lapply(ancomWZ_res_df_Acc_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.001)))
fdr_ancomWZ_Acc_Bac_0.0001 <- lapply(ancomWZ_res_df_Acc_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.0001)))

head(fdr_ancomWZ_Acc_Bac$VL)
```

```
##                       Species       lfc         se        W        p_val
## Soil_conditioningMb6   bASV_7 1.0426517 0.10462714 9.965404 2.159940e-23
## Soil_conditioningMb9  bASV_10 0.3892605 0.11432480 3.404865 6.619687e-04
## Soil_conditioningMb16 bASV_17 0.4797463 0.12068267 3.975271 7.029924e-05
## Soil_conditioningMb21 bASV_22 0.4960659 0.08522070 5.820955 5.851233e-09
## Soil_conditioningMb26 bASV_27 0.7074373 0.07166162 9.871914 5.510298e-23
## Soil_conditioningMb38 bASV_39 0.5186245 0.12356570 4.197156 2.702876e-05
##                              q_val diff_abn
## Soil_conditioningMb6  5.697922e-21     TRUE
## Soil_conditioningMb9  2.328364e-02     TRUE
## Soil_conditioningMb16 3.838147e-03     TRUE
## Soil_conditioningMb21 9.079738e-07     TRUE
## Soil_conditioningMb26 1.321470e-20     TRUE
## Soil_conditioningMb38 1.739070e-03     TRUE
```

``` r
# Number of differential ASV at different alphas
length(fdr_ancomWZ_Acc_Bac$VL$q_val)
```

```
## [1] 660
```

``` r
length(fdr_ancomWZ_Acc_Bac_0.001$VL$q_val)
```

```
## [1] 618
```

``` r
length(fdr_ancomWZ_Acc_Bac_0.0001$VL$q_val)
```

```
## [1] 614
```

``` r
save(fdr_ancomWZ_Acc_Bac, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ASVLevel/fdr_ancomWZ_Acc_Bac.RData")
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

### All data
#### Run Ancom

``` r
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_unnormalized_bac_ps_ForDA_clean.RData")

# Remove Batch 1
FP_unnormalized_bac_ps_ForDA_clean <- subset_samples(FP_unnormalized_bac_ps_ForDA_clean, Batch != 1)

ancomWZ_Bac_All <- ancombc(data = FP_unnormalized_bac_ps_ForDA_clean,
                           #tax_level = "Family",
                           formula = "Accession + Soil_conditioning",
                           group = "Soil_conditioning",
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
## The specified variables in the formula: Accession, Soil_conditioning
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
## The sample size per group is: Co = 27, Mb = 26
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
save(ancomWZ_Bac_All, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ASVLevel/ancomWZ_Bac_All.RData")
```

#### Organize and filter results

``` r
ancomWZ_res_df_All_Bac <- data.frame(
    Species = ancomWZ_Bac_All$res$lfc$taxon, # ASVs
    lfc = unlist(ancomWZ_Bac_All$res$lfc[7]),
    se = unlist(ancomWZ_Bac_All$res$se[7]),
    W = unlist(ancomWZ_Bac_All$res$W[7]), # test statistic
    p_val = unlist(ancomWZ_Bac_All$res$p_val[7]),
    q_val = unlist(ancomWZ_Bac_All$res$q_val[7]), # adjusted p-values
    diff_abn = unlist(ancomWZ_Bac_All$res$diff_abn[7]))

fdr_ancomWZ_All_Bac <- ancomWZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.05)
fdr_ancomWZ_All_Bac_0.001 <- ancomWZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.001)
fdr_ancomWZ_All_Bac_0.0001 <- ancomWZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.0001)

head(fdr_ancomWZ_All_Bac)
```

```
##                         Species         lfc       se          W        p_val
## Soil_conditioningMb137 bASV_141 -1.98898115 0.420815 -4.7264974 2.284257e-06
## Soil_conditioningMb234 bASV_240  0.54534846 0.000000  1.9001147 0.000000e+00
## Soil_conditioningMb264 bASV_270 -0.05254552 0.000000 -0.1758851 0.000000e+00
## Soil_conditioningMb275 bASV_281  0.56214427 0.000000  2.4049740 0.000000e+00
## Soil_conditioningMb280 bASV_286 -0.51437004 0.000000 -1.5982200 0.000000e+00
## Soil_conditioningMb313 bASV_320 -0.16926823 0.000000 -0.7620455 0.000000e+00
##                              q_val diff_abn
## Soil_conditioningMb137 0.004303495     TRUE
## Soil_conditioningMb234 0.000000000     TRUE
## Soil_conditioningMb264 0.000000000     TRUE
## Soil_conditioningMb275 0.000000000     TRUE
## Soil_conditioningMb280 0.000000000     TRUE
## Soil_conditioningMb313 0.000000000     TRUE
```

``` r
# Number of differential ASV at different alphas
length(fdr_ancomWZ_All_Bac$q_val)
```

```
## [1] 879
```

``` r
length(fdr_ancomWZ_All_Bac_0.001$q_val)
```

```
## [1] 878
```

``` r
length(fdr_ancomWZ_All_Bac_0.0001$q_val)
```

```
## [1] 878
```

``` r
save(fdr_ancomWZ_All_Bac, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ASVLevel/fdr_ancomWZ_All_Bac.RData")
```

``` r
# Clean environment
rm(list = ls())
```


## 5.2.1 ANCOM WITHOUT structural zeroes detection
### Accession
#### Run ANCOM

``` r
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_unnormalized_bac_ps_ForDA_clean_Acc.RData")

# Remove Batch 1
FP_unnormalized_bac_ps_ForDA_clean_Acc <- lapply(FP_unnormalized_bac_ps_ForDA_clean_Acc, function(x) subset_samples(x, Batch != 1))

ancomNoZ_Bac_Acc <- lapply(FP_unnormalized_bac_ps_ForDA_clean_Acc, function (x)
                           ancombc(data = x,
                             #tax_level = "Family",
                             formula = "Soil_conditioning",
                             group = "Soil_conditioning",
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
## The specified variables in the formula: Soil_conditioning
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
## The sample size per group is: Co = 6, Mb = 5
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
## The specified variables in the formula: Soil_conditioning
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
## The sample size per group is: Co = 6, Mb = 6
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
## The specified variables in the formula: Soil_conditioning
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
## The sample size per group is: Co = 4, Mb = 3
```

```
## Warning: Small sample size detected for the following group(s): 
## Co, Mb
## Variance estimation would be unstable when the sample size is < 5 per group
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
## The specified variables in the formula: Soil_conditioning
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
## The sample size per group is: Co = 5, Mb = 6
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
## The specified variables in the formula: Soil_conditioning
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
## The sample size per group is: Co = 6, Mb = 6
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
save(ancomNoZ_Bac_Acc, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ASVLevel/ancomNoZ_Bac_Acc.RData")
```

#### Organize and filter results

``` r
ancomNoZ_res_df_Acc_Bac <- lapply(ancomNoZ_Bac_Acc, function (x)
    data.frame(
      Species = x$res$lfc$taxon, # ASVs
      lfc = unlist(x$res$lfc[3]),
      se = unlist(x$res$se[3]),
      W = unlist(x$res$W[3]), # test statistic
      p_val = unlist(x$res$p_val[3]),
      q_val = unlist(x$res$q_val[3]), # adjusted p-values
      diff_abn = unlist(x$res$diff_abn[3]))) # true if q value is less than alpha

fdr_ancomNoZ_Acc_Bac <- lapply(ancomNoZ_res_df_Acc_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.05)))
fdr_ancomNoZ_Acc_Bac_0.001 <- lapply(ancomNoZ_res_df_Acc_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.001)))
fdr_ancomNoZ_Acc_Bac_0.0001 <- lapply(ancomNoZ_res_df_Acc_Bac, function(x)
  (x %>% dplyr::filter(q_val < 0.0001)))

head(fdr_ancomNoZ_Acc_Bac$VL)
```

```
##                       Species       lfc         se        W        p_val
## Soil_conditioningMb6   bASV_7 1.0426517 0.10462714 9.965404 2.159940e-23
## Soil_conditioningMb9  bASV_10 0.3892605 0.11432480 3.404865 6.619687e-04
## Soil_conditioningMb16 bASV_17 0.4797463 0.12068267 3.975271 7.029924e-05
## Soil_conditioningMb21 bASV_22 0.4960659 0.08522070 5.820955 5.851233e-09
## Soil_conditioningMb26 bASV_27 0.7074373 0.07166162 9.871914 5.510298e-23
## Soil_conditioningMb38 bASV_39 0.5186245 0.12356570 4.197156 2.702876e-05
##                              q_val diff_abn
## Soil_conditioningMb6  5.697922e-21     TRUE
## Soil_conditioningMb9  2.328364e-02     TRUE
## Soil_conditioningMb16 3.838147e-03     TRUE
## Soil_conditioningMb21 9.079738e-07     TRUE
## Soil_conditioningMb26 1.321470e-20     TRUE
## Soil_conditioningMb38 1.739070e-03     TRUE
```

``` r
# Number of differential ASV at different alphas
length(fdr_ancomNoZ_Acc_Bac$VL$q_val)
```

```
## [1] 114
```

``` r
length(fdr_ancomNoZ_Acc_Bac_0.001$VL$q_val)
```

```
## [1] 37
```

``` r
length(fdr_ancomNoZ_Acc_Bac_0.0001$VL$q_val)
```

```
## [1] 20
```

``` r
save(fdr_ancomNoZ_Acc_Bac, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ASVLevel/fdr_ancomNoZ_Acc_Bac.RData")

# Clean environment
rm(list = ls())
```

### All data
#### Run Ancom

``` r
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_unnormalized_bac_ps_ForDA_clean.RData")

# Remove Batch 1
FP_unnormalized_bac_ps_ForDA_clean <- subset_samples(FP_unnormalized_bac_ps_ForDA_clean, Batch != 1)

ancomNoZ_Bac_All <- ancombc(data = FP_unnormalized_bac_ps_ForDA_clean,
                           #tax_level = "Family",
                           formula = "Accession + Soil_conditioning",
                           group = "Soil_conditioning",
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
## The specified variables in the formula: Accession, Soil_conditioning
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
## The sample size per group is: Co = 27, Mb = 26
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
save(ancomNoZ_Bac_All, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ASVLevel/ancomNoZ_Bac_All.RData")
```

#### Organize and filter results

``` r
ancomNoZ_res_df_All_Bac <- data.frame(
    Species = ancomNoZ_Bac_All$res$lfc$taxon, # ASVs
    lfc = unlist(ancomNoZ_Bac_All$res$lfc[7]),
    se = unlist(ancomNoZ_Bac_All$res$se[7]),
    W = unlist(ancomNoZ_Bac_All$res$W[7]), # test statistic
    p_val = unlist(ancomNoZ_Bac_All$res$p_val[7]),
    q_val = unlist(ancomNoZ_Bac_All$res$q_val[7]), # adjusted p-values
    diff_abn = unlist(ancomNoZ_Bac_All$res$diff_abn[7]))

fdr_ancomNoZ_All_Bac <- ancomNoZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.05)
fdr_ancomNoZ_All_Bac_0.001 <- ancomNoZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.001)
fdr_ancomNoZ_All_Bac_0.0001 <- ancomNoZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.0001)

head(fdr_ancomNoZ_All_Bac)
```

```
##                           Species       lfc        se         W        p_val
## Soil_conditioningMb137   bASV_141 -1.988981 0.4208150 -4.726497 2.284257e-06
## Soil_conditioningMb1001 bASV_1164  1.653547 0.3741318  4.419692 9.884144e-06
## Soil_conditioningMb1527 bASV_2195  1.360677 0.2907397  4.680051 2.868041e-06
##                               q_val diff_abn
## Soil_conditioningMb137  0.004303495     TRUE
## Soil_conditioningMb1001 0.009887439     TRUE
## Soil_conditioningMb1527 0.004303495     TRUE
```

``` r
# Number of differential ASV at different alphas
length(fdr_ancomNoZ_All_Bac$q_val)
```

```
## [1] 3
```

``` r
length(fdr_ancomNoZ_All_Bac_0.001$q_val)
```

```
## [1] 0
```

``` r
length(fdr_ancomNoZ_All_Bac_0.0001$q_val)
```

```
## [1] 0
```

``` r
save(fdr_ancomNoZ_All_Bac, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ASVLevel/fdr_ancomNoZ_All_Bac.RData")

# Clean environment
rm(list = ls())
```

<!-- # Holm -->
<!-- # 5.2 ANCOM-BC ASV level -->
<!-- ## 5.2.1 ANCOM including structural zeroes detection -->
<!-- ### Accession -->
<!-- #### Run ANCOM -->
<!-- ```{r} -->
<!-- load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ASVLevel/Acc_unnorm_bac_ps_l_Ancom.RData") -->

<!-- ancomWZ_Bac_Acc <- lapply(Acc_unnorm_bac_ps_l_Ancom, function (x) -->
<!--                            ancombc(data = x, -->
<!--                              #tax_level = "Family", -->
<!--                              formula = "Soil_conditioning", -->
<!--                              group = "Soil_conditioning", -->
<!--                              p_adj_method = "holm", -->
<!--                              prv_cut = 0.1, #default, taxa that are in less than the fraction present are discarded -->
<!--                              lib_cut = 0, # threshold for filtering samples based on library size, 0 has no filtering -->
<!--                              struc_zero = TRUE, -->
<!--                              neg_lb = TRUE, -->
<!--                              tol = 1e-5, -->
<!--                              max_iter = 100, -->
<!--                              conserve = TRUE, -->
<!--                              alpha = 0.05, -->
<!--                              global = FALSE)) -->

<!-- save(ancomWZ_Bac_Acc, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_LessFilteredDataHolm/ancomWZ_Bac_Acc.RData") -->
<!-- ``` -->

<!-- #### Organize and filter results -->
<!-- ```{r} -->
<!-- ancomWZ_res_df_Acc_Bac <- lapply(ancomWZ_Bac_Acc, function (x) -->
<!--   data.frame( -->
<!--   Species = x$res$lfc$taxon, # ASVs -->
<!--   lfc = unlist(x$res$lfc[3]), -->
<!--   se = unlist(x$res$se[3]), -->
<!--   W = unlist(x$res$W[3]), # test statistic -->
<!--   p_val = unlist(x$res$p_val[3]), -->
<!--   q_val = unlist(x$res$q_val[3]), # adjusted p-values -->
<!--   diff_abn = unlist(x$res$diff_abn[3]))) # true if q value is less than alpha -->

<!-- fdr_ancomWZ_Acc_Bac <- lapply(ancomWZ_res_df_Acc_Bac, function(x) -->
<!--   (x %>% dplyr::filter(q_val < 0.05))) -->
<!-- fdr_ancomWZ_Acc_Bac_0.001 <- lapply(ancomWZ_res_df_Acc_Bac, function(x) -->
<!--   (x %>% dplyr::filter(q_val < 0.001))) -->
<!-- fdr_ancomWZ_Acc_Bac_0.0001 <- lapply(ancomWZ_res_df_Acc_Bac, function(x) -->
<!--   (x %>% dplyr::filter(q_val < 0.0001))) -->

<!-- head(fdr_ancomWZ_Acc_Bac$VL) -->

<!-- # Number of differential ASV at different alphas -->
<!-- length(fdr_ancomWZ_Acc_Bac$VL$q_val) -->
<!-- length(fdr_ancomWZ_Acc_Bac_0.001$VL$q_val) -->
<!-- length(fdr_ancomWZ_Acc_Bac_0.0001$VL$q_val) -->

<!-- save(fdr_ancomWZ_Acc_Bac, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_LessFilteredDataHolm/fdr_ancomWZ_Acc_Bac.RData") -->
<!-- ``` -->
<!-- ```{r} -->
<!-- # Acc_Bac_ancomNoZ_ASVs <- lapply(fdr_ancomWZ_Acc_Bac, function(x) as.vector(x[,'Species'])) -->
<!-- # fdr_ancomWZ_Acc_Bac$VL[,'Species'] -->
<!-- # fdr_ancomWZ_Acc_Bac$VL$Species -->
<!-- #  -->
<!-- # head(ancomWZ_Bac_Acc$VL$res$diff_abn, 20) -->
<!-- # summary(ancomWZ_Bac_Acc$VL$res$diff_abn) -->
<!-- # summary(ancomWZ_Bac_Acc$VL$res$q_val) -->
<!-- # length(ancomWZ_Bac_Acc$VL$res$q_val$`(Intercept)`) -->
<!-- #  -->
<!-- # head(fdr_ancomWZ_Acc_Bac$VL, 20) -->
<!-- # summary(fdr_ancomWZ_Acc_Bac$VL) -->
<!-- #  -->
<!-- # length(ancomWZ_Bac_Acc$VL$res$diff_abn) -->
<!-- # length(ancomWZ_res_df_Acc_Bac$VL$q_val) -->

<!-- # Clean environment -->
<!-- rm(list = ls()) -->
<!-- ``` -->

<!-- ### All data -->
<!-- #### Run Ancom -->
<!-- ```{r} -->
<!-- load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ASVLevel/FP_unnormalized_bac_ps_Ancom_noS.RData") -->

<!-- ancomWZ_Bac_All <- ancombc(data = FP_unnormalized_bac_ps_Ancom_noS, -->
<!--                            #tax_level = "Family", -->
<!--                            formula = "Accession + Soil_conditioning", -->
<!--                            group = "Soil_conditioning", -->
<!--                            p_adj_method = "holm", -->
<!--                            prv_cut = 0.1, #default, taxa that are in less than the fraction present are discarded -->
<!--                            lib_cut = 0, # threshold for filtering samples based on library size, 0 has no filtering -->
<!--                            struc_zero = TRUE, -->
<!--                            neg_lb = TRUE, -->
<!--                            tol = 1e-5, -->
<!--                            max_iter = 100, -->
<!--                            conserve = TRUE, -->
<!--                            alpha = 0.05, -->
<!--                            global = FALSE) -->

<!-- save(ancomWZ_Bac_All, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_LessFilteredDataHolm/ancomWZ_Bac_All.RData") -->
<!-- ``` -->

<!-- #### Organize and filter results -->
<!-- ```{r} -->
<!-- ancomWZ_res_df_All_Bac <- data.frame( -->
<!--     Species = ancomWZ_Bac_All$res$lfc$taxon, # ASVs -->
<!--     lfc = unlist(ancomWZ_Bac_All$res$lfc[7]), -->
<!--     se = unlist(ancomWZ_Bac_All$res$se[7]), -->
<!--     W = unlist(ancomWZ_Bac_All$res$W[7]), # test statistic -->
<!--     p_val = unlist(ancomWZ_Bac_All$res$p_val[7]), -->
<!--     q_val = unlist(ancomWZ_Bac_All$res$q_val[7]), # adjusted p-values -->
<!--     diff_abn = unlist(ancomWZ_Bac_All$res$diff_abn[7])) -->

<!-- fdr_ancomWZ_All_Bac <- ancomWZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.05) -->
<!-- fdr_ancomWZ_All_Bac_0.001 <- ancomWZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.001) -->
<!-- fdr_ancomWZ_All_Bac_0.0001 <- ancomWZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.0001) -->

<!-- head(fdr_ancomWZ_All_Bac) -->

<!-- # Number of differential ASV at different alphas -->
<!-- length(fdr_ancomWZ_All_Bac$q_val) -->
<!-- length(fdr_ancomWZ_All_Bac_0.001$q_val) -->
<!-- length(fdr_ancomWZ_All_Bac_0.0001$q_val) -->

<!-- save(fdr_ancomWZ_All_Bac, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_LessFilteredDataHolm/fdr_ancomWZ_All_Bac.RData") -->
<!-- ``` -->
<!-- ```{r} -->
<!-- # Clean environment -->
<!-- rm(list = ls()) -->
<!-- ``` -->


<!-- ## 5.2.1 ANCOM WITHOUT structural zeroes detection -->
<!-- ### Accession -->
<!-- #### Run ANCOM -->
<!-- ```{r} -->
<!-- load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ASVLevel/Acc_unnorm_bac_ps_l_Ancom.RData") -->

<!-- ancomNoZ_Bac_Acc <- lapply(Acc_unnorm_bac_ps_l_Ancom, function (x) -->
<!--                            ancombc(data = x, -->
<!--                              #tax_level = "Family", -->
<!--                              formula = "Soil_conditioning", -->
<!--                              group = "Soil_conditioning", -->
<!--                              p_adj_method = "holm", -->
<!--                              prv_cut = 0.1, #default, taxa that are in less than the fraction present are discarded -->
<!--                              lib_cut = 0, # threshold for filtering samples based on library size, 0 has no filtering -->
<!--                              struc_zero = FALSE, -->
<!--                              neg_lb = TRUE, -->
<!--                              tol = 1e-5, -->
<!--                              max_iter = 100, -->
<!--                              conserve = TRUE, -->
<!--                              alpha = 0.05, -->
<!--                              global = FALSE)) -->

<!-- save(ancomNoZ_Bac_Acc, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_LessFilteredDataHolm/ancomNoZ_Bac_Acc.RData") -->
<!-- ``` -->

<!-- #### Organize and filter results -->
<!-- ```{r} -->
<!-- ancomNoZ_res_df_Acc_Bac <- lapply(ancomNoZ_Bac_Acc, function (x) -->
<!--     data.frame( -->
<!--       Species = x$res$lfc$taxon, # ASVs -->
<!--       lfc = unlist(x$res$lfc[3]), -->
<!--       se = unlist(x$res$se[3]), -->
<!--       W = unlist(x$res$W[3]), # test statistic -->
<!--       p_val = unlist(x$res$p_val[3]), -->
<!--       q_val = unlist(x$res$q_val[3]), # adjusted p-values -->
<!--       diff_abn = unlist(x$res$diff_abn[3]))) # true if q value is less than alpha -->

<!-- fdr_ancomNoZ_Acc_Bac <- lapply(ancomNoZ_res_df_Acc_Bac, function(x) -->
<!--   (x %>% dplyr::filter(q_val < 0.05))) -->
<!-- fdr_ancomNoZ_Acc_Bac_0.001 <- lapply(ancomNoZ_res_df_Acc_Bac, function(x) -->
<!--   (x %>% dplyr::filter(q_val < 0.001))) -->
<!-- fdr_ancomNoZ_Acc_Bac_0.0001 <- lapply(ancomNoZ_res_df_Acc_Bac, function(x) -->
<!--   (x %>% dplyr::filter(q_val < 0.0001))) -->

<!-- head(fdr_ancomNoZ_Acc_Bac$VL) -->

<!-- # Number of differential ASV at different alphas -->
<!-- length(fdr_ancomNoZ_Acc_Bac$VL$q_val) -->
<!-- length(fdr_ancomNoZ_Acc_Bac_0.001$VL$q_val) -->
<!-- length(fdr_ancomNoZ_Acc_Bac_0.0001$VL$q_val) -->

<!-- save(fdr_ancomNoZ_Acc_Bac, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_LessFilteredDataHolm/fdr_ancomNoZ_Acc_Bac.RData") -->

<!-- # Clean environment -->
<!-- rm(list = ls()) -->
<!-- ``` -->

<!-- ### All data -->
<!-- #### Run Ancom -->
<!-- ```{r} -->
<!-- load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ASVLevel/FP_unnormalized_bac_ps_Ancom_noS.RData") -->

<!-- ancomNoZ_Bac_All <- ancombc(data = FP_unnormalized_bac_ps_Ancom_noS, -->
<!--                            #tax_level = "Family", -->
<!--                            formula = "Accession + Soil_conditioning", -->
<!--                            group = "Soil_conditioning", -->
<!--                            p_adj_method = "holm", -->
<!--                            prv_cut = 0.1, #default, taxa that are in less than the fraction present are discarded -->
<!--                            lib_cut = 0, # threshold for filtering samples based on library size, 0 has no filtering -->
<!--                            struc_zero = FALSE, -->
<!--                            neg_lb = TRUE, -->
<!--                            tol = 1e-5, -->
<!--                            max_iter = 100, -->
<!--                            conserve = TRUE, -->
<!--                            alpha = 0.05, -->
<!--                            global = FALSE) -->

<!-- save(ancomNoZ_Bac_All, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_LessFilteredDataHolm/ancomNoZ_Bac_All.RData") -->
<!-- ``` -->

<!-- #### Organize and filter results -->
<!-- ```{r} -->
<!-- ancomNoZ_res_df_All_Bac <- data.frame( -->
<!--     Species = ancomNoZ_Bac_All$res$lfc$taxon, # ASVs -->
<!--     lfc = unlist(ancomNoZ_Bac_All$res$lfc[7]), -->
<!--     se = unlist(ancomNoZ_Bac_All$res$se[7]), -->
<!--     W = unlist(ancomNoZ_Bac_All$res$W[7]), # test statistic -->
<!--     p_val = unlist(ancomNoZ_Bac_All$res$p_val[7]), -->
<!--     q_val = unlist(ancomNoZ_Bac_All$res$q_val[7]), # adjusted p-values -->
<!--     diff_abn = unlist(ancomNoZ_Bac_All$res$diff_abn[7])) -->

<!-- fdr_ancomNoZ_All_Bac <- ancomNoZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.05) -->
<!-- fdr_ancomNoZ_All_Bac_0.001 <- ancomNoZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.001) -->
<!-- fdr_ancomNoZ_All_Bac_0.0001 <- ancomNoZ_res_df_All_Bac %>% dplyr::filter(q_val < 0.0001) -->

<!-- head(fdr_ancomNoZ_All_Bac) -->

<!-- # Number of differential ASV at different alphas -->
<!-- length(fdr_ancomNoZ_All_Bac$q_val) -->
<!-- length(fdr_ancomNoZ_All_Bac_0.001$q_val) -->
<!-- length(fdr_ancomNoZ_All_Bac_0.0001$q_val) -->

<!-- save(fdr_ancomNoZ_All_Bac, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_LessFilteredDataHolm/fdr_ancomNoZ_All_Bac.RData") -->

<!-- # Clean environment -->
<!-- rm(list = ls()) -->
<!-- ``` -->
