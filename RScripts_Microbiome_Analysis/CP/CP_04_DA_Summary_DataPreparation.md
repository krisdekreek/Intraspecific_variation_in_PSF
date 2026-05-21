---
title: "CP_04_differential_abundance_Summary_Results - Data preparation"
author: "Kris de Kreek"
date: "2026-02-11"
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
library(DESeq2)
packageVersion("DESeq2")
```

```
## [1] '1.48.2'
```

``` r
library(stringr)
packageVersion("stringr")
```

```
## [1] '1.5.2'
```

``` r
library(dplyr)
packageVersion("dplyr")
```

```
## [1] '1.1.4'
```

# 4.3 Create lists of DA ASVs
## 4.3.1 Accession ASV level
### Load DA ASVs

``` r
## Load outputs DESeq
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ASVLevel/sigtab_stddds2_Wald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ASVLevel/sigtab_stddds2_LRT_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ASVLevel/sigtab_zinbWald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ASVLevel/sigtab_zinbLRT_Acc_Bac.RData")

## Load output Ancom
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/fdr_ancomWZ_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/fdr_ancomNoZ_Acc_Bac.RData")
```

### Extract DA ASVs

``` r
## From each data frame, extract the ASV IDs caught by each test
Acc_Bac_DESeqLRT_ASVs <- lapply(sigtab_stddds2_LRT_Acc_Bac, function(x) rownames(x))
Acc_Bac_DESeqWald_ASVs <- lapply(sigtab_stddds2_Wald_Acc_Bac, function(x) rownames(x))
Acc_Bac_zinbLRT_ASVs <- lapply(sigtab_zinbLRT_Acc_Bac, function(x) rownames(x))
Acc_Bac_zinbWald_ASVs <- lapply(sigtab_zinbWald_Acc_Bac, function(x) rownames(x))
Acc_Bac_ancomWZ_ASVs <- lapply(fdr_ancomWZ_Acc_Bac, function(x) as.vector(x[,'Species']))
Acc_Bac_ancomNoZ_ASVs <- lapply(fdr_ancomNoZ_Acc_Bac, function(x) as.vector(x[,'Species']))

# Within each plant concatenate the identified asvs
Acc_Bac_Total_DA_ASV <- mapply(function(s,t,w,x,y,z){
  input <- c(s,t,w,x,y,z)},
  s = Acc_Bac_DESeqLRT_ASVs,
  t = Acc_Bac_DESeqWald_ASVs,
  w = Acc_Bac_zinbLRT_ASVs,
  x = Acc_Bac_zinbWald_ASVs,
  y = Acc_Bac_ancomWZ_ASVs,
  z = Acc_Bac_ancomNoZ_ASVs,
  SIMPLIFY = FALSE)
```

### Total number of ASVs

``` r
# Total number of DA ASVs
lapply(Acc_Bac_Total_DA_ASV, function(x) length(x))
```

```
## $OH
## [1] 1156
## 
## $DD
## [1] 895
## 
## $HE
## [1] 1066
## 
## $KI
## [1] 1147
## 
## $VL
## [1] 1129
## 
## $CD
## [1] 911
## 
## $RI
## [1] 1035
## 
## $KT
## [1] 1330
## 
## $MC
## [1] 1166
## 
## $HM
## [1] 1167
## 
## $IT1
## [1] 1095
## 
## $GO1
## [1] 975
```

``` r
save(Acc_Bac_Total_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_Total_DA_ASV.RData")

# Remove duplicated asvs names from that list
Acc_Bac_NoDup_DA_ASV <- lapply(Acc_Bac_Total_DA_ASV, function(x)
  unique(x))

# Number of DA ASVs without duplicates
lapply(Acc_Bac_NoDup_DA_ASV, function(x) length(x))
```

```
## $OH
## [1] 1144
## 
## $DD
## [1] 894
## 
## $HE
## [1] 1065
## 
## $KI
## [1] 1142
## 
## $VL
## [1] 1120
## 
## $CD
## [1] 908
## 
## $RI
## [1] 1032
## 
## $KT
## [1] 1322
## 
## $MC
## [1] 1164
## 
## $HM
## [1] 1166
## 
## $IT1
## [1] 1092
## 
## $GO1
## [1] 974
```

``` r
save(Acc_Bac_NoDup_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_NoDup_DA_ASV.RData")

# Extracting taxonomy
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_unnormalized_bac_ps.RData")

# Extract taxonomy of DA ASVs
CP_unnormalized_bac_ps_Tax <- as.data.frame(tax_table(CP_unnormalized_bac_ps))
Acc_Bac_NoDup_DA_ASV_Tax <- lapply(Acc_Bac_NoDup_DA_ASV, function(x) CP_unnormalized_bac_ps_Tax[x, ])

save(Acc_Bac_NoDup_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects//CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_NoDup_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
Acc_Bac_NoDup_DA_ASV_Tax_Table <- lapply(Acc_Bac_NoDup_DA_ASV_Tax, function(df) {
  lapply(df[ , 2:6], function(col) table(col, useNA = "ifany"))
})
Acc_Bac_NoDup_DA_ASV_Tax_Table <- lapply(Acc_Bac_NoDup_DA_ASV_Tax_Table, function(ranklist) {
                                        lapply(ranklist, as.data.frame)})

# Examples
Acc_Bac_NoDup_DA_ASV_Tax_Table$HM$Class
```

```
##                            col Freq
## 1            c__Acidimicrobiia   32
## 2            c__Acidobacteriae   30
## 3            c__Actinobacteria  147
## 4                       c__AD3    2
## 5       c__Alphaproteobacteria   92
## 6              c__Anaerolineae    2
## 7                   c__Bacilli  159
## 8               c__bacteriap25    2
## 9           c__Bacteriovoracia    1
## 10              c__Bacteroidia   28
## 11 c__BD2-11_terrestrial_group    1
## 12          c__Bdellovibrionia    3
## 13           c__Blastocatellia   15
## 14             c__Chloroflexia   14
## 15               c__Clostridia    9
## 16       c__Desulfitobacteriia    2
## 17         c__Desulfotomaculia    2
## 18      c__Gammaproteobacteria  170
## 19           c__Gemmatimonadia   55
## 20              c__Gitt-GS-136    2
## 21          c__Gracilibacteria    1
## 22               c__Holophagae   13
## 23           c__Incertae_Sedis    2
## 24             c__JG30-KF-CM66    4
## 25                   c__KD4-96   12
## 26          c__Ktedonobacteria   24
## 27             c__Limnochordia    6
## 28           c__Longimicrobiia    3
## 29                c__MB-A2-108    7
## 30         c__Methylomirabilia    7
## 31               c__Myxococcia    3
## 32              c__Nitrospiria    2
## 33                    c__OLB14    1
## 34              c__Oligoflexia    2
## 35               c__Polyangiia   29
## 36  c__S0134_terrestrial_group    3
## 37          c__Saccharimonadia    1
## 38              c__Subgroup_22    1
## 39              c__Subgroup_25    1
## 40               c__Subgroup_5    1
## 41               c__Sumerlaeia    1
## 42          c__Symbiobacteriia    2
## 43         c__Syntrophomonadia    1
## 44        c__Thermaerobacteria    1
## 45      c__Thermoanaerobaculia    1
## 46          c__Thermoleophilia  170
## 47                     c__TK10    2
## 48         c__Vampirivibrionia    1
## 49         c__Verrucomicrobiia   13
## 50         c__Vicinamibacteria   77
## 51                        <NA>    6
```

``` r
Acc_Bac_NoDup_DA_ASV_Tax_Table$VL$Class
```

```
##                            col Freq
## 1            c__Acidimicrobiia   20
## 2            c__Acidobacteriae   26
## 3            c__Actinobacteria  125
## 4                       c__AD3    6
## 5       c__Alphaproteobacteria   93
## 6              c__Anaerolineae    1
## 7                   c__Bacilli  111
## 8               c__bacteriap25    1
## 9           c__Bacteriovoracia    1
## 10              c__Bacteroidia   26
## 11 c__BD2-11_terrestrial_group    1
## 12          c__Bdellovibrionia    2
## 13           c__Blastocatellia   20
## 14             c__Chloroflexia   13
## 15               c__Clostridia    6
## 16               c__Deinococci    1
## 17            c__Fibrobacteria    2
## 18           c__Fimbriimonadia    1
## 19      c__Gammaproteobacteria  174
## 20           c__Gemmatimonadia   57
## 21            c__Halanaerobiia    1
## 22               c__Holophagae   11
## 23           c__Incertae_Sedis    2
## 24             c__JG30-KF-CM66    4
## 25                   c__KD4-96   25
## 26          c__Ktedonobacteria   28
## 27             c__Limnochordia    6
## 28           c__Longimicrobiia    2
## 29                c__MB-A2-108    9
## 30         c__Methylomirabilia    5
## 31               c__Myxococcia    4
## 32              c__Nitrospiria    4
## 33              c__Oligoflexia    2
## 34                   c__P2-11E    1
## 35            c__Phycisphaerae    1
## 36               c__Polyangiia   19
## 37  c__S0134_terrestrial_group    3
## 38               c__Subgroup_5    1
## 39            c__Sulfobacillia    1
## 40               c__Sumerlaeia    1
## 41          c__Symbiobacteriia    4
## 42          c__Thermacetogenia    1
## 43      c__Thermoanaerobaculia    4
## 44          c__Thermoleophilia  164
## 45                     c__TK10    3
## 46         c__Verrucomicrobiia   21
## 47         c__Vicinamibacteria  100
## 48                        <NA>    6
```

``` r
save(Acc_Bac_NoDup_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_NoDup_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring twice

``` r
# Filter to ASVs detected by at least 2 tests
Acc_Bac_NumbOcc_DA_ASV <- lapply(Acc_Bac_Total_DA_ASV, function(x) table(unlist(x)))
Acc_Bac_TwoTimes_DA_ASV <- lapply(Acc_Bac_NumbOcc_DA_ASV, function(x) names(x[x >= 2]))
Acc_Bac_TwoTimes_DA_ASV
```

```
## $OH
## [1] "bASV_1224" "bASV_1272" "bASV_165"  "bASV_1760" "bASV_387"  "bASV_586" 
## [7] "bASV_724"  "bASV_868"  "bASV_911" 
## 
## $DD
## [1] "bASV_269"
## 
## $HE
## [1] "bASV_2629"
## 
## $KI
## [1] "bASV_1376" "bASV_1842" "bASV_4526" "bASV_601"  "bASV_656" 
## 
## $VL
## [1] "bASV_1147" "bASV_3119" "bASV_315"  "bASV_581"  "bASV_601"  "bASV_72"  
## [7] "bASV_792"  "bASV_996" 
## 
## $CD
## [1] "bASV_2350" "bASV_2677" "bASV_724" 
## 
## $RI
## [1] "bASV_1310" "bASV_1656"
## 
## $KT
## [1] "bASV_1272" "bASV_1314" "bASV_1413" "bASV_1760" "bASV_2121" "bASV_2416"
## [7] "bASV_724"  "bASV_911" 
## 
## $MC
## [1] "bASV_1258" "bASV_502" 
## 
## $HM
## [1] "bASV_2978"
## 
## $IT1
## [1] "bASV_1822" "bASV_2045" "bASV_3764"
## 
## $GO1
## [1] "bASV_733"
```

``` r
# Number of DA ASVs occurring in two tests
lapply(Acc_Bac_TwoTimes_DA_ASV, function(x) length(x))
```

```
## $OH
## [1] 9
## 
## $DD
## [1] 1
## 
## $HE
## [1] 1
## 
## $KI
## [1] 5
## 
## $VL
## [1] 8
## 
## $CD
## [1] 3
## 
## $RI
## [1] 2
## 
## $KT
## [1] 8
## 
## $MC
## [1] 2
## 
## $HM
## [1] 1
## 
## $IT1
## [1] 3
## 
## $GO1
## [1] 1
```

``` r
save(Acc_Bac_TwoTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_TwoTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
Acc_Bac_TwoTimes_DA_ASV_Tax <- lapply(Acc_Bac_TwoTimes_DA_ASV, function(x) CP_unnormalized_bac_ps_Tax[x, ])

save(Acc_Bac_TwoTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_TwoTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
Acc_Bac_TwoTimes_DA_ASV_Tax_Table <- lapply(Acc_Bac_TwoTimes_DA_ASV_Tax, function(df) {
  lapply(df[ , 2:6], function(col) table(col, useNA = "ifany"))
})
Acc_Bac_TwoTimes_DA_ASV_Tax_Table <- lapply(Acc_Bac_TwoTimes_DA_ASV_Tax_Table, function(ranklist) {
                                            lapply(ranklist, as.data.frame)})

# Examples
Acc_Bac_TwoTimes_DA_ASV_Tax_Table$HM$Class
```

```
##                  col Freq
## 1 c__Thermoleophilia    1
```

``` r
Acc_Bac_TwoTimes_DA_ASV_Tax_Table$VL$Class
```

```
##                      col Freq
## 1      c__Acidobacteriae    1
## 2      c__Actinobacteria    1
## 3 c__Alphaproteobacteria    1
## 4 c__Gammaproteobacteria    3
## 5              c__KD4-96    1
## 6     c__Thermoleophilia    1
```

``` r
save(Acc_Bac_TwoTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_TwoTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring three times

``` r
# Filter to ASVs detected by at least 3 tests
Acc_Bac_NumbOcc_DA_ASV <- lapply(Acc_Bac_Total_DA_ASV, function(x) table(unlist(x)))
Acc_Bac_ThreeTimes_DA_ASV <- lapply(Acc_Bac_NumbOcc_DA_ASV, function(x) names(x[x >= 3]))
Acc_Bac_ThreeTimes_DA_ASV
```

```
## $OH
## [1] "bASV_724" "bASV_911"
## 
## $DD
## character(0)
## 
## $HE
## character(0)
## 
## $KI
## character(0)
## 
## $VL
## [1] "bASV_72"
## 
## $CD
## character(0)
## 
## $RI
## [1] "bASV_1310"
## 
## $KT
## character(0)
## 
## $MC
## character(0)
## 
## $HM
## character(0)
## 
## $IT1
## character(0)
## 
## $GO1
## character(0)
```

``` r
# Number of DA ASVs occurring in three tests
lapply(Acc_Bac_ThreeTimes_DA_ASV, function(x) length(x))
```

```
## $OH
## [1] 2
## 
## $DD
## [1] 0
## 
## $HE
## [1] 0
## 
## $KI
## [1] 0
## 
## $VL
## [1] 1
## 
## $CD
## [1] 0
## 
## $RI
## [1] 1
## 
## $KT
## [1] 0
## 
## $MC
## [1] 0
## 
## $HM
## [1] 0
## 
## $IT1
## [1] 0
## 
## $GO1
## [1] 0
```

``` r
save(Acc_Bac_ThreeTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_ThreeTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
Acc_Bac_ThreeTimes_DA_ASV_Tax <- lapply(Acc_Bac_ThreeTimes_DA_ASV, function(x) CP_unnormalized_bac_ps_Tax[x, ])

save(Acc_Bac_ThreeTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_ThreeTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
Acc_Bac_ThreeTimes_DA_ASV_Tax_Table <- lapply(Acc_Bac_ThreeTimes_DA_ASV_Tax, function(df) {
  lapply(df[ , 2:6], function(col) table(col, useNA = "ifany"))
})
Acc_Bac_ThreeTimes_DA_ASV_Tax_Table <- lapply(Acc_Bac_ThreeTimes_DA_ASV_Tax_Table, function(ranklist) {
                                            lapply(ranklist, as.data.frame)})

# Examples
Acc_Bac_ThreeTimes_DA_ASV_Tax_Table$HM$Class
```

```
## [1] Freq
## <0 rows> (or 0-length row.names)
```

``` r
Acc_Bac_ThreeTimes_DA_ASV_Tax_Table$VL$Class
```

```
##                 col Freq
## 1 c__Acidobacteriae    1
```

``` r
save(Acc_Bac_ThreeTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_ThreeTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring four times

``` r
# Filter to ASVs detected by at least 4 tests
Acc_Bac_FourTimes_DA_ASV <- lapply(Acc_Bac_NumbOcc_DA_ASV, function(x) names(x[x >= 4]))
Acc_Bac_FourTimes_DA_ASV
```

```
## $OH
## [1] "bASV_724"
## 
## $DD
## character(0)
## 
## $HE
## character(0)
## 
## $KI
## character(0)
## 
## $VL
## character(0)
## 
## $CD
## character(0)
## 
## $RI
## character(0)
## 
## $KT
## character(0)
## 
## $MC
## character(0)
## 
## $HM
## character(0)
## 
## $IT1
## character(0)
## 
## $GO1
## character(0)
```

``` r
# Number of DA ASVs occurring in four tests
lapply(Acc_Bac_FourTimes_DA_ASV, function(x) length(x))
```

```
## $OH
## [1] 1
## 
## $DD
## [1] 0
## 
## $HE
## [1] 0
## 
## $KI
## [1] 0
## 
## $VL
## [1] 0
## 
## $CD
## [1] 0
## 
## $RI
## [1] 0
## 
## $KT
## [1] 0
## 
## $MC
## [1] 0
## 
## $HM
## [1] 0
## 
## $IT1
## [1] 0
## 
## $GO1
## [1] 0
```

``` r
save(Acc_Bac_FourTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_FourTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
Acc_Bac_FourTimes_DA_ASV_Tax <- lapply(Acc_Bac_FourTimes_DA_ASV, function(x) CP_unnormalized_bac_ps_Tax[x, ])

save(Acc_Bac_FourTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_FourTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
Acc_Bac_FourTimes_DA_ASV_Tax_Table <- lapply(Acc_Bac_FourTimes_DA_ASV_Tax, function(df) {
  lapply(df[ , 2:6], function(col) table(col, useNA = "ifany"))
})
Acc_Bac_FourTimes_DA_ASV_Tax_Table <- lapply(Acc_Bac_FourTimes_DA_ASV_Tax_Table, function(ranklist) {
                                            lapply(ranklist, as.data.frame)})

# Examples
Acc_Bac_FourTimes_DA_ASV_Tax_Table$HM$Class
```

```
## [1] Freq
## <0 rows> (or 0-length row.names)
```

``` r
Acc_Bac_FourTimes_DA_ASV_Tax_Table$VL$Class
```

```
## [1] Freq
## <0 rows> (or 0-length row.names)
```

``` r
save(Acc_Bac_FourTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_FourTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring five times

``` r
# Filter to ASVs detected by at least 5 tests
Acc_Bac_FiveTimes_DA_ASV <- lapply(Acc_Bac_NumbOcc_DA_ASV, function(x) names(x[x >= 5]))
Acc_Bac_FiveTimes_DA_ASV
```

```
## $OH
## character(0)
## 
## $DD
## character(0)
## 
## $HE
## character(0)
## 
## $KI
## character(0)
## 
## $VL
## character(0)
## 
## $CD
## character(0)
## 
## $RI
## character(0)
## 
## $KT
## character(0)
## 
## $MC
## character(0)
## 
## $HM
## character(0)
## 
## $IT1
## character(0)
## 
## $GO1
## character(0)
```

``` r
# Number of DA ASVs occurring in Five tests
lapply(Acc_Bac_FiveTimes_DA_ASV, function(x) length(x))
```

```
## $OH
## [1] 0
## 
## $DD
## [1] 0
## 
## $HE
## [1] 0
## 
## $KI
## [1] 0
## 
## $VL
## [1] 0
## 
## $CD
## [1] 0
## 
## $RI
## [1] 0
## 
## $KT
## [1] 0
## 
## $MC
## [1] 0
## 
## $HM
## [1] 0
## 
## $IT1
## [1] 0
## 
## $GO1
## [1] 0
```

``` r
save(Acc_Bac_FiveTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_FiveTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
Acc_Bac_FiveTimes_DA_ASV_Tax <- lapply(Acc_Bac_FiveTimes_DA_ASV, function(x) CP_unnormalized_bac_ps_Tax[x, ])

save(Acc_Bac_FiveTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_FiveTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
Acc_Bac_FiveTimes_DA_ASV_Tax_Table <- lapply(Acc_Bac_FiveTimes_DA_ASV_Tax, function(df) {
  lapply(df[ , 2:6], function(col) table(col, useNA = "ifany"))
})
Acc_Bac_FiveTimes_DA_ASV_Tax_Table <- lapply(Acc_Bac_FiveTimes_DA_ASV_Tax_Table, function(ranklist) {
                                            lapply(ranklist, as.data.frame)})

# Examples
Acc_Bac_FiveTimes_DA_ASV_Tax_Table$HM$Class
```

```
## [1] Freq
## <0 rows> (or 0-length row.names)
```

``` r
Acc_Bac_FiveTimes_DA_ASV_Tax_Table$VL$Class
```

```
## [1] Freq
## <0 rows> (or 0-length row.names)
```

``` r
save(Acc_Bac_FiveTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_FiveTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring six times

``` r
# Filter to ASVs detected by at least 6 tests
Acc_Bac_SixTimes_DA_ASV <- lapply(Acc_Bac_NumbOcc_DA_ASV, function(x) names(x[x >= 6]))
Acc_Bac_SixTimes_DA_ASV
```

```
## $OH
## character(0)
## 
## $DD
## character(0)
## 
## $HE
## character(0)
## 
## $KI
## character(0)
## 
## $VL
## character(0)
## 
## $CD
## character(0)
## 
## $RI
## character(0)
## 
## $KT
## character(0)
## 
## $MC
## character(0)
## 
## $HM
## character(0)
## 
## $IT1
## character(0)
## 
## $GO1
## character(0)
```

``` r
# Number of DA ASVs occurring in Six tests
lapply(Acc_Bac_SixTimes_DA_ASV, function(x) length(x))
```

```
## $OH
## [1] 0
## 
## $DD
## [1] 0
## 
## $HE
## [1] 0
## 
## $KI
## [1] 0
## 
## $VL
## [1] 0
## 
## $CD
## [1] 0
## 
## $RI
## [1] 0
## 
## $KT
## [1] 0
## 
## $MC
## [1] 0
## 
## $HM
## [1] 0
## 
## $IT1
## [1] 0
## 
## $GO1
## [1] 0
```

``` r
save(Acc_Bac_SixTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_SixTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
Acc_Bac_SixTimes_DA_ASV_Tax <- lapply(Acc_Bac_SixTimes_DA_ASV, function(x) CP_unnormalized_bac_ps_Tax[x, ])

save(Acc_Bac_SixTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_SixTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
Acc_Bac_SixTimes_DA_ASV_Tax_Table <- lapply(Acc_Bac_SixTimes_DA_ASV_Tax, function(df) {
  lapply(df[ , 2:6], function(col) table(col, useNA = "ifany"))
})
Acc_Bac_SixTimes_DA_ASV_Tax_Table <- lapply(Acc_Bac_SixTimes_DA_ASV_Tax_Table, function(ranklist) {
                                            lapply(ranklist, as.data.frame)})

# Print output
Acc_Bac_SixTimes_DA_ASV_Tax_Table
```

```
## $OH
## $OH$Phylum
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $OH$Class
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $OH$Order
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $OH$Family
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $OH$Genus
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## 
## $DD
## $DD$Phylum
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $DD$Class
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $DD$Order
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $DD$Family
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $DD$Genus
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## 
## $HE
## $HE$Phylum
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $HE$Class
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $HE$Order
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $HE$Family
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $HE$Genus
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## 
## $KI
## $KI$Phylum
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $KI$Class
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $KI$Order
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $KI$Family
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $KI$Genus
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## 
## $VL
## $VL$Phylum
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $VL$Class
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $VL$Order
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $VL$Family
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $VL$Genus
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## 
## $CD
## $CD$Phylum
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $CD$Class
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $CD$Order
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $CD$Family
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $CD$Genus
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## 
## $RI
## $RI$Phylum
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $RI$Class
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $RI$Order
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $RI$Family
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $RI$Genus
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## 
## $KT
## $KT$Phylum
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $KT$Class
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $KT$Order
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $KT$Family
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $KT$Genus
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## 
## $MC
## $MC$Phylum
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $MC$Class
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $MC$Order
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $MC$Family
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $MC$Genus
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## 
## $HM
## $HM$Phylum
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $HM$Class
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $HM$Order
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $HM$Family
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $HM$Genus
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## 
## $IT1
## $IT1$Phylum
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $IT1$Class
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $IT1$Order
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $IT1$Family
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $IT1$Genus
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## 
## $GO1
## $GO1$Phylum
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $GO1$Class
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $GO1$Order
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $GO1$Family
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $GO1$Genus
## [1] Freq
## <0 rows> (or 0-length row.names)
```

``` r
save(Acc_Bac_SixTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_SixTimes_DA_ASV_Tax_Table.RData")
```

``` r
# Clean environment
rm(list = ls())
```


## 4.3.2 Accession Family level
### Load DA Family

``` r
## Load outputs DESeq
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_FamilyLevel/sigtab_stddds2_Wald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_FamilyLevel/sigtab_stddds2_LRT_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_FamilyLevel/sigtab_zinbWald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_FamilyLevel/sigtab_zinbLRT_Acc_Bac.RData")

## Load output Ancom
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_FamilyLevel/fdr_ancomWZ_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_FamilyLevel/fdr_ancomNoZ_Acc_Bac.RData")

## Load Phyloseq object
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_unnormalized_bac_ps_ForDA.RData")
```

### Extract DA Family
#### Extracting data from tests

``` r
# From each data frame, extract the Family or ASV IDs caught by each test
Acc_Bac_DESeqLRT_ASV_Family <- lapply(sigtab_stddds2_LRT_Acc_Bac, function(x) rownames(x))
Acc_Bac_DESeqWald_ASV_Family <- lapply(sigtab_stddds2_Wald_Acc_Bac, function(x) rownames(x))
Acc_Bac_zinbLRT_ASV_Family <- lapply(sigtab_zinbLRT_Acc_Bac, function(x) rownames(x))
Acc_Bac_zinbWald_ASV_Family <- lapply(sigtab_zinbWald_Acc_Bac, function(x) rownames(x))
Acc_Bac_ancomWZ_Family <- lapply(fdr_ancomWZ_Acc_Bac, function(x) as.vector(x[,'Species']))
Acc_Bac_ancomNoZ_Family <- lapply(fdr_ancomNoZ_Acc_Bac, function(x) as.vector(x[,'Species']))
```

#### Merging tests and transposing list

``` r
Acc_DESeq_DA_ASV_Family_l <- list(DESeqLRT = Acc_Bac_DESeqLRT_ASV_Family, DESeqWald = Acc_Bac_DESeqWald_ASV_Family, zinbLRT = Acc_Bac_zinbLRT_ASV_Family, zinbWald = Acc_Bac_zinbWald_ASV_Family)

# transpose list
Acc <- names(Acc_DESeq_DA_ASV_Family_l[[1]])
Acc_DESeq_DA_ASV_Family_l_t <- vector("list", length(Acc))
names(Acc_DESeq_DA_ASV_Family_l_t) <- Acc

for (cd in Acc) {
  Acc_DESeq_DA_ASV_Family_l_t[[cd]] <- lapply(Acc_DESeq_DA_ASV_Family_l, function(m) m[[cd]])
}

## DA Family Ancom
Acc_Ancom_DA_Family_l <- list(ancomWZ = Acc_Bac_ancomWZ_Family, ancomNoZ = Acc_Bac_ancomNoZ_Family)

# transpose list
Acc_Ancom_DA_Family_l_t <- vector("list", length(Acc))
names(Acc_Ancom_DA_Family_l_t) <- Acc

for (cd in Acc) {
  Acc_Ancom_DA_Family_l_t[[cd]] <- lapply(Acc_Ancom_DA_Family_l, function(m) m[[cd]])
}
```

### Extract taxomomy
#### DESeq

``` r
tax <- as.data.frame(tax_table(CP_unnormalized_bac_ps_ForDA))
tax$ASV <- rownames(tax)

Acc_NoDup_DESeq_DA_ASV_Family_Tax <- lapply(Acc_DESeq_DA_ASV_Family_l_t, function(x) {
  lapply(x, function(y) {
    unique(tax[tax$ASV %in% y, c("Family", "Order", "Class", "Phylum")])
    })})
```

#### Acnom
##### Fix issues with Unknown and Unclassified taxa

``` r
# fix_collapsed_taxonomy <- function(df) {
#   for (i in seq_len(nrow(df))) {
#     
#     # PHYLUM LEVEL FIX
#     if (!is.na(df$Phylum[i]) && grepl("_Unknown|_Unclassified", df$Phylum[i])) {
#       parts <- strsplit(sub("^p__", "", df$Phylum[i]), "_")[[1]]
#       df$Phylum[i] <- paste0("p__", parts[1])
#       lower <- parts[-1]
#       if (length(lower) >= 1) df$Class[i]  <- paste0("c__", lower[1])
#       if (length(lower) >= 2) df$Order[i]  <- paste0("o__", lower[2])
#       if (length(lower) >= 3) df$Family[i] <- paste0("f__", lower[3])
#     }
#     
#     # CLASS LEVEL FIX
#     if (!is.na(df$Class[i]) && grepl("_Unknown|_Unclassified", df$Class[i])) {
#       parts <- strsplit(sub("^c__", "", df$Class[i]), "_")[[1]]
#       df$Class[i] <- paste0("c__", parts[1])
#       lower <- parts[-1]
#       if (length(lower) >= 1) df$Order[i]  <- paste0("o__", lower[1])
#       if (length(lower) >= 2) df$Family[i] <- paste0("f__", lower[2])
#     }
#     
#     # ORDER LEVEL FIX
#     if (!is.na(df$Order[i]) && grepl("_Unclassified", df$Order[i])) {
#       parts <- strsplit(sub("^o__", "", df$Order[i]), "_")[[1]]
#       df$Order[i] <- paste0("o__", parts[1])
#       lower <- parts[-1]
#       if (length(lower) >= 1) df$Family[i] <- paste0("f__", lower[1])
#     }
#   }
#   return(df)
# }
```

##### Extract taxonomy

``` r
#lapply(Acc_Ancom_DA_Family_l, function(x) length(x))

Acc_NoDup_Ancom_DA_Family_tax <- lapply(Acc_Ancom_DA_Family_l_t, function(y) {
  lapply(y, function(x) {
    # Take unique names in case of duplicates
    DA_tax <- unique(x)
    
    # Split data between long and short taxonomic names
    full_tax_strings   <- DA_tax[grepl("k__", DA_tax)]
    family_only_strings <- DA_tax[!grepl("k__", DA_tax)]
    
    # Extract taxomomy from long names
    tax_full <- tibble(Taxon = full_tax_strings) %>%
      mutate(
        Phylum = str_extract(Taxon, "p__.*?(?=_c__|$)"),
        Class  = str_extract(Taxon, "c__.*?(?=_o__|$)"),
        Order  = str_extract(Taxon, "o__.*?(?=_f__|$)"),
        Family = str_extract(Taxon, "f__.*?(?=_g__|$)")
      ) %>%
      select(Phylum, Class, Order, Family)
    
    # Extract taxonomy from ps object for short names
    tax_ps <- as.data.frame(tax_table(CP_unnormalized_bac_ps_ForDA))
    #family_only <- sub("^f__", "", family_only_strings)
    
    tax_family_map <- unique(tax_ps[tax_ps$Family %in% family_only_strings, c("Family", "Order", "Class", "Phylum")])
    #nrow(tax_family_map)
    
    # Remove left over Incertae Sedis from short names object
    tax_family_map <- tax_family_map %>%
      filter(Family != "f__Incertae_Sedis")
   
    # Combine long and short names and replace NAs with unkown and unclassified
    tax_combined <- as.data.frame(bind_rows(tax_full, tax_family_map))
    #tax_combined <- fix_collapsed_taxonomy(tax_combined)
})})

#lapply(Acc_NoDup_Ancom_DA_Family_tax, function(x) nrow(x))
```

#### Merge Ancom and DESeq

``` r
Acc_NoDup_DESeq_DA_ASV_Family_Tax_m <- lapply(Acc_NoDup_DESeq_DA_ASV_Family_Tax, function(x) do.call(rbind, x))
Acc_NoDup_Ancom_DA_Family_tax_m <- lapply(Acc_NoDup_Ancom_DA_Family_tax, function(x) do.call(rbind, x))

# Merge and filter on unique family
Acc_Total_Merged_DA_Family <- mapply(function(y,z){
  input <- rbind(y,z)},
  y = Acc_NoDup_DESeq_DA_ASV_Family_Tax_m,
  z = Acc_NoDup_Ancom_DA_Family_tax_m,
  SIMPLIFY = FALSE)

lapply(Acc_Total_Merged_DA_Family, function(x) nrow(x))
```

```
## $OH
## [1] 40
## 
## $DD
## [1] 34
## 
## $HE
## [1] 39
## 
## $KI
## [1] 48
## 
## $VL
## [1] 40
## 
## $CD
## [1] 51
## 
## $RI
## [1] 31
## 
## $KT
## [1] 36
## 
## $MC
## [1] 44
## 
## $HM
## [1] 37
## 
## $IT1
## [1] 44
## 
## $GO1
## [1] 58
```

``` r
lapply(Acc_Total_Merged_DA_Family, function(x) table(x$Family))
```

```
## $OH
## 
##      f__Actinospicaceae       f__Alcaligenaceae     f__Caloramatoraceae 
##                       1                       1                       1 
##  f__Cryptosporangiaceae f__Desulfitibacteraceae            f__Family_XI 
##                       1                       1                       1 
##    f__Fimbriimonadaceae    f__Flavobacteriaceae          f__Gemmataceae 
##                       1                       1                       1 
##     f__Haloplasmataceae    f__Heliobacteriaceae        f__Holophagaceae 
##                       1                       1                       1 
##       f__Incertae_Sedis     f__Micavibrionaceae        f__Moraxellaceae 
##                      13                       1                       1 
##     f__Natronincolaceae   f__Obscuribacteraceae          f__Opitutaceae 
##                       1                       1                       1 
##     f__Proteiniboraceae    f__Rickettsiellaceae   f__Saccharimonadaceae 
##                       1                       1                       1 
## f__Sedimentibacteraceae     f__Sulfobacillaceae    f__Terrimicrobiaceae 
##                       1                       1                       1 
##         f__Unclassified 
##                       4 
## 
## $DD
## 
##                 f__B1-7BS  f__Caldicoprobacteraceae       f__Caloramatoraceae 
##                         1                         1                         1 
##    f__Christensenellaceae   f__CPla-3_termite_group              f__Family_XI 
##                         1                         1                         1 
##          f__Holophagaceae         f__Incertae_Sedis            f__Kaistiaceae 
##                         1                        10                         1 
##       f__Micavibrionaceae            f__Opitutaceae          f__Pirellulaceae 
##                         1                         1                         1 
##         f__Rhodocyclaceae        f__Rhodothermaceae      f__Rickettsiellaceae 
##                         1                         1                         1 
##   f__Sedimentibacteraceae       f__Sulfobacillaceae          f__Sumerlaeaceae 
##                         1                         1                         1 
##    f__Thermacetogeniaceae           f__Unclassified                   f__WX65 
##                         1                         4                         1 
## f__Xiphinematobacteraceae 
##                         1 
## 
## $HE
## 
##               f__A0839     f__Actinospicaceae     f__Anaerolineaceae 
##                      1                      1                      1 
##      f__Caldilineaceae    f__Caloramatoraceae   f__Catenulisporaceae 
##                      1                      1                      1 
##     f__Chloroflexaceae   f__Dethiobacteraceae          f__Dongiaceae 
##                      1                      1                      1 
##          f__env.OPS_17   f__Fimbriimonadaceae   f__Flavobacteriaceae 
##                      1                      1                      1 
##       f__Holophagaceae       f__Holosporaceae      f__Incertae_Sedis 
##                      1                      1                     12 
##      f__Legionellaceae    f__Magnetospiraceae   f__Methylomonadaceae 
##                      1                      1                      1 
##    f__Micavibrionaceae       f__Moraxellaceae         f__Opitutaceae 
##                      1                      1                      1 
##       f__Paracoccaceae    f__Proteiniboraceae      f__Rhodocyclaceae 
##                      1                      1                      1 
##  f__Saccharimonadaceae       f__Sumerlaeaceae f__Syntrophomonadaceae 
##                      1                      1                      1 
##        f__Unclassified 
##                      1 
## 
## $KI
## 
##              f__27F-1492R        f__Actinospicaceae        f__Anaerolineaceae 
##                         1                         1                         1 
##                  f__BSV26      f__Crocinitomicaceae     f__Ethanoligenenaceae 
##                         1                         1                         1 
##      f__Flavobacteriaceae            f__Gemmataceae     f__Gracilibacteraceae 
##                         1                         1                         1 
##      f__Heliobacteriaceae         f__Incertae_Sedis         f__Isosphaeraceae 
##                         1                        17                         1 
##       f__Magnetospiraceae        f__Marinococcaceae       f__Micavibrionaceae 
##                         1                         1                         1 
##   f__NS11-12_marine_group            f__Opitutaceae       f__Oscillospiraceae 
##                         1                         1                         2 
##   f__Paracaedibacteraceae       f__Proteiniboraceae         f__Rhodocyclaceae 
##                         1                         1                         1 
##      f__Rhodospirillaceae      f__Rickettsiellaceae          f__Sporomusaceae 
##                         1                         1                         1 
##          f__Sumerlaeaceae      f__Terrimicrobiaceae           f__Unclassified 
##                         1                         1                         3 
##         f__Vermiphilaceae f__Xiphinematobacteraceae 
##                         1                         1 
## 
## $VL
## 
##                   f__Actinospicaceae                    f__Caldilineaceae 
##                                    1                                    1 
##                 f__Chitinimonadaceae                    f__Deinococcaceae 
##                                    1                                    1 
## f__Desulfotomaculales_Incertae_Sedis                        f__Dongiaceae 
##                                    1                                    1 
##                     f__Holophagaceae                     f__Holosporaceae 
##                                    1                                    1 
##                    f__Incertae_Sedis                            f__KD3-93 
##                                   14                                    1 
##               f__Magnetospirillaceae                 f__Methylomonadaceae 
##                                    1                                    1 
##                  f__Micavibrionaceae                    f__Oligoflexaceae 
##                                    1                                    1 
##              f__Propionibacteriaceae                    f__Rhodocyclaceae 
##                                    1                                    1 
##                 f__Rhodospirillaceae                 f__Rickettsiellaceae 
##                                    1                                    1 
##                 f__Sanguibacteraceae                            f__SM2D12 
##                                    1                                    1 
##               f__Syntrophomonadaceae                   f__Thermincolaceae 
##                                    1                                    1 
##                      f__Unclassified                f__Vulgatibacteraceae 
##                                    3                                    1 
##                             f__WWH38 
##                                    1 
## 
## $CD
## 
##                             f__37-13                             f__A0839 
##                                    1                                    1 
##                               f__A4b                   f__Actinospicaceae 
##                                    1                                    1 
##                           f__AKIW781                    f__Alcaligenaceae 
##                                    1                                    1 
##             f__Caldicoprobacteraceae                    f__Caldilineaceae 
##                                    1                                    1 
##                  f__Caloramatoraceae                 f__Catenulisporaceae 
##                                    1                                    1 
##                     f__Cytophagaceae                    f__Deinococcaceae 
##                                    1                                    1 
##             f__Desulfitobacteriaceae f__Desulfotomaculales_Incertae_Sedis 
##                                    1                                    1 
##                        f__env.OPS_17                 f__Flavobacteriaceae 
##                                    1                                    1 
##                     f__Holosporaceae                    f__Incertae_Sedis 
##                                    1                                   14 
##                              f__LWQ8               f__Magnetospirillaceae 
##                                    1                                    1 
##                   f__Marinococcaceae                  f__Micavibrionaceae 
##                                    1                                    1 
##                   f__Microscillaceae                   f__Microtrichaceae 
##                                    1                                    1 
##                f__Obscuribacteraceae                       f__Opitutaceae 
##                                    1                                    1 
##                  f__Pelotomaculaceae                  f__Pseudomonadaceae 
##                                    1                                    1 
##                    f__Rhodocyclaceae                  f__Sulfobacillaceae 
##                                    1                                    1 
##                     f__Sumerlaeaceae               f__Syntrophomonadaceae 
##                                    1                                    1 
##                 f__Tepidisphaeraceae               f__Thermacetogeniaceae 
##                                    1                                    1 
##                      f__Unclassified                f__Vulgatibacteraceae 
##                                    3                                    1 
## 
## $RI
## 
##                f__A0839         f__Amb-16S-1323               f__B1-7BS 
##                       1                       1                       1 
##      f__Chloroflexaceae f__CPla-3_termite_group    f__Dethiobacteraceae 
##                       1                       1                       1 
##          f__Elsteraceae     f__Haloplasmataceae        f__Holophagaceae 
##                       1                       1                       1 
##        f__Holosporaceae       f__Incertae_Sedis  f__Magnetospirillaceae 
##                       1                      14                       1 
##        f__Moraxellaceae   f__Obscuribacteraceae          f__Opitutaceae 
##                       1                       1                       1 
## f__Paracaedibacteraceae    f__Rickettsiellaceae      f__Spirosomataceae 
##                       1                       1                       1 
## 
## $KT
## 
##       f__Alcaligenaceae         f__Amb-16S-1323      f__Anaerolineaceae 
##                       1                       1                       1 
##  f__Christensenellaceae        f__Cytophagaceae f__Desulfitibacteraceae 
##                       1                       1                       1 
##       f__Geobacteraceae   f__Gracilibacteraceae    f__Heliobacteriaceae 
##                       1                       1                       1 
##        f__Holophagaceae       f__Incertae_Sedis       f__Isosphaeraceae 
##                       1                      12                       1 
##    f__Methylomonadaceae       f__Oligoflexaceae          f__Opitutaceae 
##                       1                       1                       1 
##   f__Proteinivoracaceae     f__Pseudomonadaceae    f__Rickettsiellaceae 
##                       1                       1                       1 
##   f__Saccharimonadaceae         f__Unclassified          f__Woeseiaceae 
##                       1                       5                       1 
## 
## $MC
## 
##                   f__Actinospicaceae                    f__Alcaligenaceae 
##                                    1                                    1 
##                 f__Defluviicoccaceae f__Desulfotomaculales_Incertae_Sedis 
##                                    1                                    1 
##                        f__Dongiaceae                 f__Fimbriimonadaceae 
##                                    1                                    1 
##                 f__Flavobacteriaceae                    f__Geobacteraceae 
##                                    1                                    1 
##                 f__Heliobacteriaceae                     f__Holophagaceae 
##                                    1                                    1 
##                    f__Incertae_Sedis                  f__Magnetospiraceae 
##                                   11                                    1 
##                     f__Moraxellaceae                       f__Opitutaceae 
##                                    1                                    2 
##              f__Paracaedibacteraceae                   f__Pedosphaeraceae 
##                                    1                                    1 
##                  f__Pelotomaculaceae              f__Propionibacteriaceae 
##                                    1                                    1 
##                  f__Proteiniboraceae                  f__Pseudomonadaceae 
##                                    1                                    1 
##                    f__Rhodocyclaceae                 f__Rickettsiellaceae 
##                                    1                                    1 
##                f__Saccharimonadaceae                              f__SRB2 
##                                    1                                    1 
##               f__Syntrophomonadaceae                 f__Terrimicrobiaceae 
##                                    1                                    1 
##                      f__Unclassified                f__Vulgatibacteraceae 
##                                    4                                    1 
##                             f__WWH38            f__Xiphinematobacteraceae 
##                                    1                                    1 
## 
## $HM
## 
##            f__27F-1492R      f__Actinospicaceae       f__Alcaligenaceae 
##                       1                       1                       1 
##      f__Chloroflexaceae    f__Dethiobacteraceae           f__Dongiaceae 
##                       1                       1                       1 
##   f__Gracilibacteraceae     f__Haloplasmataceae    f__Heliobacteriaceae 
##                       1                       1                       1 
##       f__Incertae_Sedis          f__Kaistiaceae               f__KD3-93 
##                       9                       1                       1 
##                 f__LWQ8     f__Magnetospiraceae f__Paracaedibacteraceae 
##                       1                       1                       1 
##    f__Parachlamydiaceae       f__Peptococcaceae f__Sedimentibacteraceae 
##                       1                       1                       1 
##      f__Spirosomataceae                 f__SRB2        f__Sumerlaeaceae 
##                       1                       1                       2 
## f__Thermaerobacteraceae      f__Thermincolaceae         f__Unclassified 
##                       1                       1                       4 
##   f__Vulgatibacteraceae 
##                       1 
## 
## $IT1
## 
##                         f__27F-1492R                             f__A0839 
##                                    1                                    1 
##                               f__A4b                               f__AB1 
##                                    1                                    1 
##                   f__Anaerolineaceae                  f__Caloramatoraceae 
##                                    1                                    1 
##                   f__Chloroflexaceae               f__Christensenellaceae 
##                                    1                                    1 
## f__Desulfotomaculales_Incertae_Sedis                       f__Elsteraceae 
##                                    1                                    1 
##                 f__Fimbriimonadaceae                    f__Incertae_Sedis 
##                                    1                                   18 
##                    f__Isosphaeraceae                   f__Marinococcaceae 
##                                    1                                    1 
##                f__Obscuribacteraceae                    f__Oligoflexaceae 
##                                    1                                    1 
##              f__Paracaedibacteraceae                  f__Pelotomaculaceae 
##                                    1                                    1 
##              f__Propionibacteriaceae                  f__Pseudomonadaceae 
##                                    1                                    1 
##                f__Saccharimonadaceae                   f__Spirosomataceae 
##                                    1                                    1 
##                              f__SRB2                 f__Terrimicrobiaceae 
##                                    1                                    1 
##                      f__Unclassified 
##                                    3 
## 
## $GO1
## 
##                             f__A0839                               f__A4b 
##                                    1                                    1 
##                   f__Actinospicaceae                            f__B1-7BS 
##                                    1                                    1 
##                  f__Caloramatoraceae                 f__Catenulisporaceae 
##                                    1                                    1 
##                 f__Chitinimonadaceae                 f__Crocinitomicaceae 
##                                    1                                    1 
##                     f__Cytophagaceae               f__Desulfotomaculaceae 
##                                    1                                    1 
## f__Desulfotomaculales_Incertae_Sedis                        f__Dongiaceae 
##                                    1                                    1 
##                       f__Gemmataceae                  f__Haloplasmataceae 
##                                    1                                    1 
##                     f__Holophagaceae                     f__Holosporaceae 
##                                    1                                    1 
##                    f__Incertae_Sedis                       f__Kaistiaceae 
##                                   19                                    1 
##                  f__Leptolyngbyaceae                  f__Magnetospiraceae 
##                                    1                                    1 
##                   f__Microscillaceae                    f__Oligoflexaceae 
##                                    1                                    1 
##              f__Paracaedibacteraceae                  f__Pelotomaculaceae 
##                                    1                                    1 
##                     f__Pirellulaceae                  f__Pseudomonadaceae 
##                                    1                                    1 
##                    f__Rhodocyclaceae                 f__Rhodospirillaceae 
##                                    1                                    1 
##                 f__Rickettsiellaceae              f__Sedimentibacteraceae 
##                                    1                                    1 
##                   f__Spirosomataceae                  f__Sulfobacillaceae 
##                                    1                                    1 
##                 f__Terrimicrobiaceae               f__Thermacetogeniaceae 
##                                    1                                    1 
##                           f__UCG-010                      f__Unclassified 
##                                    1                                    4 
##                f__Vulgatibacteraceae 
##                                    1
```

``` r
save(Acc_Total_Merged_DA_Family, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_FamilyLevel/Acc_Total_Merged_DA_Family.RData")
```

### Remove duplicat Family

``` r
Acc_NoDup_Merged_DA_Family <- lapply(Acc_Total_Merged_DA_Family, function(x) unique(x))
lapply(Acc_NoDup_Merged_DA_Family, function(x) nrow(x))
```

```
## $OH
## [1] 40
## 
## $DD
## [1] 34
## 
## $HE
## [1] 39
## 
## $KI
## [1] 47
## 
## $VL
## [1] 40
## 
## $CD
## [1] 51
## 
## $RI
## [1] 31
## 
## $KT
## [1] 36
## 
## $MC
## [1] 43
## 
## $HM
## [1] 36
## 
## $IT1
## [1] 44
## 
## $GO1
## [1] 57
```

``` r
lapply(Acc_NoDup_Merged_DA_Family, function(x) table(x$Family))
```

```
## $OH
## 
##      f__Actinospicaceae       f__Alcaligenaceae     f__Caloramatoraceae 
##                       1                       1                       1 
##  f__Cryptosporangiaceae f__Desulfitibacteraceae            f__Family_XI 
##                       1                       1                       1 
##    f__Fimbriimonadaceae    f__Flavobacteriaceae          f__Gemmataceae 
##                       1                       1                       1 
##     f__Haloplasmataceae    f__Heliobacteriaceae        f__Holophagaceae 
##                       1                       1                       1 
##       f__Incertae_Sedis     f__Micavibrionaceae        f__Moraxellaceae 
##                      13                       1                       1 
##     f__Natronincolaceae   f__Obscuribacteraceae          f__Opitutaceae 
##                       1                       1                       1 
##     f__Proteiniboraceae    f__Rickettsiellaceae   f__Saccharimonadaceae 
##                       1                       1                       1 
## f__Sedimentibacteraceae     f__Sulfobacillaceae    f__Terrimicrobiaceae 
##                       1                       1                       1 
##         f__Unclassified 
##                       4 
## 
## $DD
## 
##                 f__B1-7BS  f__Caldicoprobacteraceae       f__Caloramatoraceae 
##                         1                         1                         1 
##    f__Christensenellaceae   f__CPla-3_termite_group              f__Family_XI 
##                         1                         1                         1 
##          f__Holophagaceae         f__Incertae_Sedis            f__Kaistiaceae 
##                         1                        10                         1 
##       f__Micavibrionaceae            f__Opitutaceae          f__Pirellulaceae 
##                         1                         1                         1 
##         f__Rhodocyclaceae        f__Rhodothermaceae      f__Rickettsiellaceae 
##                         1                         1                         1 
##   f__Sedimentibacteraceae       f__Sulfobacillaceae          f__Sumerlaeaceae 
##                         1                         1                         1 
##    f__Thermacetogeniaceae           f__Unclassified                   f__WX65 
##                         1                         4                         1 
## f__Xiphinematobacteraceae 
##                         1 
## 
## $HE
## 
##               f__A0839     f__Actinospicaceae     f__Anaerolineaceae 
##                      1                      1                      1 
##      f__Caldilineaceae    f__Caloramatoraceae   f__Catenulisporaceae 
##                      1                      1                      1 
##     f__Chloroflexaceae   f__Dethiobacteraceae          f__Dongiaceae 
##                      1                      1                      1 
##          f__env.OPS_17   f__Fimbriimonadaceae   f__Flavobacteriaceae 
##                      1                      1                      1 
##       f__Holophagaceae       f__Holosporaceae      f__Incertae_Sedis 
##                      1                      1                     12 
##      f__Legionellaceae    f__Magnetospiraceae   f__Methylomonadaceae 
##                      1                      1                      1 
##    f__Micavibrionaceae       f__Moraxellaceae         f__Opitutaceae 
##                      1                      1                      1 
##       f__Paracoccaceae    f__Proteiniboraceae      f__Rhodocyclaceae 
##                      1                      1                      1 
##  f__Saccharimonadaceae       f__Sumerlaeaceae f__Syntrophomonadaceae 
##                      1                      1                      1 
##        f__Unclassified 
##                      1 
## 
## $KI
## 
##              f__27F-1492R        f__Actinospicaceae        f__Anaerolineaceae 
##                         1                         1                         1 
##                  f__BSV26      f__Crocinitomicaceae     f__Ethanoligenenaceae 
##                         1                         1                         1 
##      f__Flavobacteriaceae            f__Gemmataceae     f__Gracilibacteraceae 
##                         1                         1                         1 
##      f__Heliobacteriaceae         f__Incertae_Sedis         f__Isosphaeraceae 
##                         1                        17                         1 
##       f__Magnetospiraceae        f__Marinococcaceae       f__Micavibrionaceae 
##                         1                         1                         1 
##   f__NS11-12_marine_group            f__Opitutaceae       f__Oscillospiraceae 
##                         1                         1                         1 
##   f__Paracaedibacteraceae       f__Proteiniboraceae         f__Rhodocyclaceae 
##                         1                         1                         1 
##      f__Rhodospirillaceae      f__Rickettsiellaceae          f__Sporomusaceae 
##                         1                         1                         1 
##          f__Sumerlaeaceae      f__Terrimicrobiaceae           f__Unclassified 
##                         1                         1                         3 
##         f__Vermiphilaceae f__Xiphinematobacteraceae 
##                         1                         1 
## 
## $VL
## 
##                   f__Actinospicaceae                    f__Caldilineaceae 
##                                    1                                    1 
##                 f__Chitinimonadaceae                    f__Deinococcaceae 
##                                    1                                    1 
## f__Desulfotomaculales_Incertae_Sedis                        f__Dongiaceae 
##                                    1                                    1 
##                     f__Holophagaceae                     f__Holosporaceae 
##                                    1                                    1 
##                    f__Incertae_Sedis                            f__KD3-93 
##                                   14                                    1 
##               f__Magnetospirillaceae                 f__Methylomonadaceae 
##                                    1                                    1 
##                  f__Micavibrionaceae                    f__Oligoflexaceae 
##                                    1                                    1 
##              f__Propionibacteriaceae                    f__Rhodocyclaceae 
##                                    1                                    1 
##                 f__Rhodospirillaceae                 f__Rickettsiellaceae 
##                                    1                                    1 
##                 f__Sanguibacteraceae                            f__SM2D12 
##                                    1                                    1 
##               f__Syntrophomonadaceae                   f__Thermincolaceae 
##                                    1                                    1 
##                      f__Unclassified                f__Vulgatibacteraceae 
##                                    3                                    1 
##                             f__WWH38 
##                                    1 
## 
## $CD
## 
##                             f__37-13                             f__A0839 
##                                    1                                    1 
##                               f__A4b                   f__Actinospicaceae 
##                                    1                                    1 
##                           f__AKIW781                    f__Alcaligenaceae 
##                                    1                                    1 
##             f__Caldicoprobacteraceae                    f__Caldilineaceae 
##                                    1                                    1 
##                  f__Caloramatoraceae                 f__Catenulisporaceae 
##                                    1                                    1 
##                     f__Cytophagaceae                    f__Deinococcaceae 
##                                    1                                    1 
##             f__Desulfitobacteriaceae f__Desulfotomaculales_Incertae_Sedis 
##                                    1                                    1 
##                        f__env.OPS_17                 f__Flavobacteriaceae 
##                                    1                                    1 
##                     f__Holosporaceae                    f__Incertae_Sedis 
##                                    1                                   14 
##                              f__LWQ8               f__Magnetospirillaceae 
##                                    1                                    1 
##                   f__Marinococcaceae                  f__Micavibrionaceae 
##                                    1                                    1 
##                   f__Microscillaceae                   f__Microtrichaceae 
##                                    1                                    1 
##                f__Obscuribacteraceae                       f__Opitutaceae 
##                                    1                                    1 
##                  f__Pelotomaculaceae                  f__Pseudomonadaceae 
##                                    1                                    1 
##                    f__Rhodocyclaceae                  f__Sulfobacillaceae 
##                                    1                                    1 
##                     f__Sumerlaeaceae               f__Syntrophomonadaceae 
##                                    1                                    1 
##                 f__Tepidisphaeraceae               f__Thermacetogeniaceae 
##                                    1                                    1 
##                      f__Unclassified                f__Vulgatibacteraceae 
##                                    3                                    1 
## 
## $RI
## 
##                f__A0839         f__Amb-16S-1323               f__B1-7BS 
##                       1                       1                       1 
##      f__Chloroflexaceae f__CPla-3_termite_group    f__Dethiobacteraceae 
##                       1                       1                       1 
##          f__Elsteraceae     f__Haloplasmataceae        f__Holophagaceae 
##                       1                       1                       1 
##        f__Holosporaceae       f__Incertae_Sedis  f__Magnetospirillaceae 
##                       1                      14                       1 
##        f__Moraxellaceae   f__Obscuribacteraceae          f__Opitutaceae 
##                       1                       1                       1 
## f__Paracaedibacteraceae    f__Rickettsiellaceae      f__Spirosomataceae 
##                       1                       1                       1 
## 
## $KT
## 
##       f__Alcaligenaceae         f__Amb-16S-1323      f__Anaerolineaceae 
##                       1                       1                       1 
##  f__Christensenellaceae        f__Cytophagaceae f__Desulfitibacteraceae 
##                       1                       1                       1 
##       f__Geobacteraceae   f__Gracilibacteraceae    f__Heliobacteriaceae 
##                       1                       1                       1 
##        f__Holophagaceae       f__Incertae_Sedis       f__Isosphaeraceae 
##                       1                      12                       1 
##    f__Methylomonadaceae       f__Oligoflexaceae          f__Opitutaceae 
##                       1                       1                       1 
##   f__Proteinivoracaceae     f__Pseudomonadaceae    f__Rickettsiellaceae 
##                       1                       1                       1 
##   f__Saccharimonadaceae         f__Unclassified          f__Woeseiaceae 
##                       1                       5                       1 
## 
## $MC
## 
##                   f__Actinospicaceae                    f__Alcaligenaceae 
##                                    1                                    1 
##                 f__Defluviicoccaceae f__Desulfotomaculales_Incertae_Sedis 
##                                    1                                    1 
##                        f__Dongiaceae                 f__Fimbriimonadaceae 
##                                    1                                    1 
##                 f__Flavobacteriaceae                    f__Geobacteraceae 
##                                    1                                    1 
##                 f__Heliobacteriaceae                     f__Holophagaceae 
##                                    1                                    1 
##                    f__Incertae_Sedis                  f__Magnetospiraceae 
##                                   11                                    1 
##                     f__Moraxellaceae                       f__Opitutaceae 
##                                    1                                    1 
##              f__Paracaedibacteraceae                   f__Pedosphaeraceae 
##                                    1                                    1 
##                  f__Pelotomaculaceae              f__Propionibacteriaceae 
##                                    1                                    1 
##                  f__Proteiniboraceae                  f__Pseudomonadaceae 
##                                    1                                    1 
##                    f__Rhodocyclaceae                 f__Rickettsiellaceae 
##                                    1                                    1 
##                f__Saccharimonadaceae                              f__SRB2 
##                                    1                                    1 
##               f__Syntrophomonadaceae                 f__Terrimicrobiaceae 
##                                    1                                    1 
##                      f__Unclassified                f__Vulgatibacteraceae 
##                                    4                                    1 
##                             f__WWH38            f__Xiphinematobacteraceae 
##                                    1                                    1 
## 
## $HM
## 
##            f__27F-1492R      f__Actinospicaceae       f__Alcaligenaceae 
##                       1                       1                       1 
##      f__Chloroflexaceae    f__Dethiobacteraceae           f__Dongiaceae 
##                       1                       1                       1 
##   f__Gracilibacteraceae     f__Haloplasmataceae    f__Heliobacteriaceae 
##                       1                       1                       1 
##       f__Incertae_Sedis          f__Kaistiaceae               f__KD3-93 
##                       9                       1                       1 
##                 f__LWQ8     f__Magnetospiraceae f__Paracaedibacteraceae 
##                       1                       1                       1 
##    f__Parachlamydiaceae       f__Peptococcaceae f__Sedimentibacteraceae 
##                       1                       1                       1 
##      f__Spirosomataceae                 f__SRB2        f__Sumerlaeaceae 
##                       1                       1                       1 
## f__Thermaerobacteraceae      f__Thermincolaceae         f__Unclassified 
##                       1                       1                       4 
##   f__Vulgatibacteraceae 
##                       1 
## 
## $IT1
## 
##                         f__27F-1492R                             f__A0839 
##                                    1                                    1 
##                               f__A4b                               f__AB1 
##                                    1                                    1 
##                   f__Anaerolineaceae                  f__Caloramatoraceae 
##                                    1                                    1 
##                   f__Chloroflexaceae               f__Christensenellaceae 
##                                    1                                    1 
## f__Desulfotomaculales_Incertae_Sedis                       f__Elsteraceae 
##                                    1                                    1 
##                 f__Fimbriimonadaceae                    f__Incertae_Sedis 
##                                    1                                   18 
##                    f__Isosphaeraceae                   f__Marinococcaceae 
##                                    1                                    1 
##                f__Obscuribacteraceae                    f__Oligoflexaceae 
##                                    1                                    1 
##              f__Paracaedibacteraceae                  f__Pelotomaculaceae 
##                                    1                                    1 
##              f__Propionibacteriaceae                  f__Pseudomonadaceae 
##                                    1                                    1 
##                f__Saccharimonadaceae                   f__Spirosomataceae 
##                                    1                                    1 
##                              f__SRB2                 f__Terrimicrobiaceae 
##                                    1                                    1 
##                      f__Unclassified 
##                                    3 
## 
## $GO1
## 
##                             f__A0839                               f__A4b 
##                                    1                                    1 
##                   f__Actinospicaceae                            f__B1-7BS 
##                                    1                                    1 
##                  f__Caloramatoraceae                 f__Catenulisporaceae 
##                                    1                                    1 
##                 f__Chitinimonadaceae                 f__Crocinitomicaceae 
##                                    1                                    1 
##                     f__Cytophagaceae               f__Desulfotomaculaceae 
##                                    1                                    1 
## f__Desulfotomaculales_Incertae_Sedis                        f__Dongiaceae 
##                                    1                                    1 
##                       f__Gemmataceae                  f__Haloplasmataceae 
##                                    1                                    1 
##                     f__Holophagaceae                     f__Holosporaceae 
##                                    1                                    1 
##                    f__Incertae_Sedis                       f__Kaistiaceae 
##                                   18                                    1 
##                  f__Leptolyngbyaceae                  f__Magnetospiraceae 
##                                    1                                    1 
##                   f__Microscillaceae                    f__Oligoflexaceae 
##                                    1                                    1 
##              f__Paracaedibacteraceae                  f__Pelotomaculaceae 
##                                    1                                    1 
##                     f__Pirellulaceae                  f__Pseudomonadaceae 
##                                    1                                    1 
##                    f__Rhodocyclaceae                 f__Rhodospirillaceae 
##                                    1                                    1 
##                 f__Rickettsiellaceae              f__Sedimentibacteraceae 
##                                    1                                    1 
##                   f__Spirosomataceae                  f__Sulfobacillaceae 
##                                    1                                    1 
##                 f__Terrimicrobiaceae               f__Thermacetogeniaceae 
##                                    1                                    1 
##                           f__UCG-010                      f__Unclassified 
##                                    1                                    4 
##                f__Vulgatibacteraceae 
##                                    1
```

``` r
save(Acc_NoDup_Merged_DA_Family, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_FamilyLevel/Acc_NoDup_Merged_DA_Family.RData")
```

### DA Family occurring twice

``` r
# Filter to Family detected by at least 2 tests
Acc_TwoTimes_Merged_DA_Family <- lapply(Acc_Total_Merged_DA_Family, function(x) {
    dup <- x[duplicated(x) | duplicated(x, fromLast = TRUE), , drop = FALSE]
    unique(dup)
  })

# Number of DA Family occurring in two tests
lapply(Acc_TwoTimes_Merged_DA_Family, function(x) nrow(x))
```

```
## $OH
## [1] 0
## 
## $DD
## [1] 0
## 
## $HE
## [1] 0
## 
## $KI
## [1] 1
## 
## $VL
## [1] 0
## 
## $CD
## [1] 0
## 
## $RI
## [1] 0
## 
## $KT
## [1] 0
## 
## $MC
## [1] 1
## 
## $HM
## [1] 1
## 
## $IT1
## [1] 0
## 
## $GO1
## [1] 1
```

``` r
save(Acc_TwoTimes_Merged_DA_Family, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_FamilyLevel/Acc_TwoTimes_Merged_DA_Family.RData")
```

``` r
# Clean environment
rm(list = ls())
```


## 4.3.3 Accession Order level
### Load DA Order

``` r
## Load outputs DESeq
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_OrderLevel/sigtab_stddds2_Wald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_OrderLevel/sigtab_stddds2_LRT_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_OrderLevel/sigtab_zinbWald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_OrderLevel/sigtab_zinbLRT_Acc_Bac.RData")

## Load output Ancom
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_OrderLevel/fdr_ancomWZ_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_OrderLevel/fdr_ancomNoZ_Acc_Bac.RData")

## Load Phyloseq object
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_unnormalized_bac_ps_ForDA.RData")
```

### Extract DA Order
#### Extracting data from tests

``` r
# From each data frame, extract the Order or ASV IDs caught by each test
Acc_Bac_DESeqLRT_ASV_Order <- lapply(sigtab_stddds2_LRT_Acc_Bac, function(x) rownames(x))
Acc_Bac_DESeqWald_ASV_Order <- lapply(sigtab_stddds2_Wald_Acc_Bac, function(x) rownames(x))
Acc_Bac_zinbLRT_ASV_Order <- lapply(sigtab_zinbLRT_Acc_Bac, function(x) rownames(x))
Acc_Bac_zinbWald_ASV_Order <- lapply(sigtab_zinbWald_Acc_Bac, function(x) rownames(x))
Acc_Bac_ancomWZ_Order <- lapply(fdr_ancomWZ_Acc_Bac, function(x) as.vector(x[,'Species']))
Acc_Bac_ancomNoZ_Order <- lapply(fdr_ancomNoZ_Acc_Bac, function(x) as.vector(x[,'Species']))
```

#### Merging tests and transposing list

``` r
Acc_DESeq_DA_ASV_Order_l <- list(DESeqLRT = Acc_Bac_DESeqLRT_ASV_Order, DESeqWald = Acc_Bac_DESeqWald_ASV_Order, zinbLRT = Acc_Bac_zinbLRT_ASV_Order, zinbWald = Acc_Bac_zinbWald_ASV_Order)

# transpose list
Acc <- names(Acc_DESeq_DA_ASV_Order_l[[1]])
Acc_DESeq_DA_ASV_Order_l_t <- vector("list", length(Acc))
names(Acc_DESeq_DA_ASV_Order_l_t) <- Acc

for (cd in Acc) {
  Acc_DESeq_DA_ASV_Order_l_t[[cd]] <- lapply(Acc_DESeq_DA_ASV_Order_l, function(m) m[[cd]])
}

## DA Order Ancom
Acc_Ancom_DA_Order_l <- list(ancomWZ = Acc_Bac_ancomWZ_Order, ancomNoZ = Acc_Bac_ancomNoZ_Order)

# transpose list
Acc_Ancom_DA_Order_l_t <- vector("list", length(Acc))
names(Acc_Ancom_DA_Order_l_t) <- Acc

for (cd in Acc) {
  Acc_Ancom_DA_Order_l_t[[cd]] <- lapply(Acc_Ancom_DA_Order_l, function(m) m[[cd]])
}
```

### Extract taxomomy
#### DESeq

``` r
tax <- as.data.frame(tax_table(CP_unnormalized_bac_ps_ForDA))
tax$ASV <- rownames(tax)

Acc_NoDup_DESeq_DA_ASV_Order_Tax <- lapply(Acc_DESeq_DA_ASV_Order_l_t, function(x) {
  lapply(x, function(y) {
    unique(tax[tax$ASV %in% y, c("Order", "Class", "Phylum")])
    })})
```

#### Acnom

``` r
#lapply(Acc_Ancom_DA_Order_l, function(x) length(x))

Acc_NoDup_Ancom_DA_Order_tax <- lapply(Acc_Ancom_DA_Order_l_t, function(y) {
  lapply(y, function(x) {
    # Take unique names in case of duplicates
    DA_tax <- unique(x)
    
    # Split data between long and short taxonomic names
    full_tax_strings   <- DA_tax[grepl("k__", DA_tax)]
    family_only_strings <- DA_tax[!grepl("k__", DA_tax)]
    
    # Extract taxomomy from long names
    tax_full <- tibble(Taxon = full_tax_strings) %>%
      mutate(
        Phylum = str_extract(Taxon, "p__.*?(?=_c__|$)"),
        Class  = str_extract(Taxon, "c__.*?(?=_o__|$)"),
        Order  = str_extract(Taxon, "o__.*?(?=_f__|$)")
      ) %>%
      select(Phylum, Class, Order)

    # Extract taxonomy from ps object for short names
    tax_ps <- as.data.frame(tax_table(CP_unnormalized_bac_ps_ForDA))
    #family_only <- sub("^f__", "", family_only_strings)
    
    tax_family_map <- unique(tax_ps[tax_ps$Order %in% family_only_strings, c("Order", "Class", "Phylum")])
    #nrow(tax_family_map)
    
    # Remove left over Incertae Sedis from short names object
    tax_family_map <- tax_family_map %>%
      filter(Order != "f__Incertae_Sedis")
   
    # Combine long and short names
    tax_combined <- as.data.frame(bind_rows(tax_full, tax_family_map))
})})

#lapply(Acc_NoDup_Ancom_DA_Order_tax, function(x) nrow(x))
```

#### Merge Ancom and DESeq

``` r
Acc_NoDup_DESeq_DA_ASV_Order_Tax_m <- lapply(Acc_NoDup_DESeq_DA_ASV_Order_Tax, function(x) do.call(rbind, x))
Acc_NoDup_Ancom_DA_Order_tax_m <- lapply(Acc_NoDup_Ancom_DA_Order_tax, function(x) do.call(rbind, x))

# Merge and filter on unique family
Acc_Total_Merged_DA_Order <- mapply(function(y,z){
  input <- rbind(y,z)},
  y = Acc_NoDup_DESeq_DA_ASV_Order_Tax_m,
  z = Acc_NoDup_Ancom_DA_Order_tax_m,
  SIMPLIFY = FALSE)

lapply(Acc_Total_Merged_DA_Order, function(x) nrow(x))
```

```
## $OH
## [1] 24
## 
## $DD
## [1] 26
## 
## $HE
## [1] 29
## 
## $KI
## [1] 31
## 
## $VL
## [1] 32
## 
## $CD
## [1] 27
## 
## $RI
## [1] 20
## 
## $KT
## [1] 25
## 
## $MC
## [1] 27
## 
## $HM
## [1] 23
## 
## $IT1
## [1] 25
## 
## $GO1
## [1] 32
```

``` r
lapply(Acc_Total_Merged_DA_Order, function(x) table(x$Order))
```

```
## $OH
## 
##     o__Actinomarinales             o__AKIW659             o__Blfdi19 
##                      1                      1                      1 
## o__Desulfitibacterales    o__Fimbriimonadales    o__Flavobacteriales 
##                      1                      1                      1 
##          o__Gemmatales     o__Haloplasmatales        o__Holophagales 
##                      1                      1                      1 
##      o__Incertae_Sedis          o__Lineage_IV     o__Micavibrionales 
##                      8                      1                      1 
##   o__Obscuribacterales          o__Opitutales    o__Rickettsiellales 
##                      1                      1                      1 
##     o__Sulfobacillales        o__Unclassified 
##                      1                      1 
## 
## $DD
## 
##                o__11-24              o__AKIW659             o__B10-SB3A 
##                       1                       1                       1 
##           o__Babeliales o__Caldicoprobacterales   o__Christensenellales 
##                       1                       1                       1 
##     o__Defluviicoccales         o__Holophagales       o__Incertae_Sedis 
##                       1                       1                       6 
##      o__Micavibrionales           o__Opitutales         o__Pirellulales 
##                       1                       1                       1 
##       o__Rhodothermales     o__Rickettsiellales      o__Sulfobacillales 
##                       1                       1                       1 
##         o__Sumerlaeales   o__Thermacetogeniales         o__Unclassified 
##                       1                       1                       4 
## 
## $HE
## 
##            o__AKIW659     o__Anaerolineales           o__B10-SB3A 
##                     1                     1                     1 
##         o__Babeliales      o__Caldilineales       o__Cytophagales 
##                     1                     1                     3 
##   o__Dethiobacterales          o__Dongiales   o__Fimbriimonadales 
##                     1                     1                     1 
##       o__Holophagales       o__Holosporales     o__Incertae_Sedis 
##                     1                     1                     7 
##      o__Legionellales    o__Methylococcales    o__Micavibrionales 
##                     1                     1                     1 
##         o__Opitutales              o__PeM15    o__Rhodobacterales 
##                     1                     1                     1 
##       o__Sumerlaeales o__Syntrophomonadales       o__Unclassified 
##                     1                     1                     1 
## 
## $KI
## 
##                o__Actinomarinales                        o__AKIW659 
##                                 1                                 1 
##                 o__Anaerolineales                      o__B12-WMSP1 
##                                 1                                 1 
##               o__Defluviicoccales                         o__DTU014 
##                                 1                                 1 
##                     o__Gemmatales                 o__Incertae_Sedis 
##                                 1                                 7 
##                  o__Isosphaerales                o__Kapabacteriales 
##                                 1                                 1 
##                   o__Kryptoniales                     o__Lineage_IV 
##                                 1                                 1 
##                o__Methylococcales                o__Micavibrionales 
##                                 1                                 1 
##                     o__Opitutales                o__Oscillospirales 
##                                 1                                 3 
##            o__Paracaedibacterales               o__Rickettsiellales 
##                                 1                                 1 
##                    o__Subgroup_13                   o__Sumerlaeales 
##                                 1                                 1 
##                   o__Unclassified o__Veillonellales-Selenomonadales 
##                                 2                                 1 
## 
## $VL
## 
##  o__Aggregatilineales            o__AKIW659           o__B10-SB3A 
##                     1                     1                     1 
##         o__Babeliales      o__Caldilineales   o__Defluviicoccales 
##                     1                     1                     1 
##      o__Deinococcales          o__Dongiales             o__DS-100 
##                     1                     1                     1 
##             o__DTU014       o__Holophagales       o__Holosporales 
##                     1                     1                     1 
##     o__Incertae_Sedis         o__Lineage_IV    o__Methylococcales 
##                     6                     1                     1 
##    o__Micavibrionales           o__MSB-4B10      o__Oligoflexales 
##                     1                     1                     1 
##             o__PLTA13      o__Rickettsiales   o__Rickettsiellales 
##                     1                     1                     1 
## o__Syntrophomonadales     o__Thermincolales       o__Unclassified 
##                     1                     1                     3 
## o__Vampirovibrionales 
##                     1 
## 
## $CD
## 
##    o__Aggregatilineales             o__B10-SB3A o__Caldicoprobacterales 
##                       1                       1                       1 
##        o__Caldilineales        o__Deinococcales         o__Holosporales 
##                       1                       1                       1 
##       o__Incertae_Sedis        o__Kallotenuales      o__Kapabacteriales 
##                       7                       1                       1 
##      o__Micavibrionales             o__MSB-4B10    o__Obscuribacterales 
##                       1                       1                       1 
##           o__Opitutales                o__PeM15      o__Sulfobacillales 
##                       1                       1                       1 
##         o__Sumerlaeales   o__Syntrophomonadales   o__Thermacetogeniales 
##                       1                       1                       1 
##         o__Unclassified 
##                       3 
## 
## $RI
## 
##   o__Aggregatilineales             o__AKIW659    o__Dethiobacterales 
##                      1                      1                      1 
##              o__DTU014           o__Elev-1554     o__Haloplasmatales 
##                      1                      1                      1 
##        o__Holophagales        o__Holosporales      o__Incertae_Sedis 
##                      1                      1                      4 
##     o__Kapabacteriales               o__MBA03   o__Obscuribacterales 
##                      1                      1                      1 
##          o__Opitutales o__Paracaedibacterales               o__PeM15 
##                      1                      1                      1 
##    o__Rickettsiellales  o__Vampirovibrionales 
##                      1                      1 
## 
## $KT
## 
##   o__Aggregatilineales             o__AKIW659      o__Anaerolineales 
##                      1                      1                      1 
##  o__Christensenellales o__Desulfitibacterales              o__DS-100 
##                      1                      1                      1 
##       o__Geobacterales        o__Holophagales      o__Incertae_Sedis 
##                      1                      1                      7 
##       o__Isosphaerales       o__Oligoflexales          o__Opitutales 
##                      1                      1                      1 
##   o__Proteinivoracales    o__Rickettsiellales              o__SJA-15 
##                      1                      1                      1 
##         o__Subgroup_13        o__Unclassified 
##                      1                      3 
## 
## $MC
## 
##                  o__11-24        o__Actinomarinales             o__Babeliales 
##                         1                         1                         1 
##              o__Dongiales       o__Fimbriimonadales       o__Flavobacteriales 
##                         1                         1                         1 
##          o__Geobacterales           o__Holophagales         o__Incertae_Sedis 
##                         1                         1                         7 
##        o__Kapabacteriales             o__Opitutales    o__Paracaedibacterales 
##                         1                         2                         1 
##         o__Pedosphaerales       o__Rickettsiellales     o__Syntrophomonadales 
##                         1                         1                         1 
## o__Thermoanaerobacterales           o__Unclassified     o__Vampirovibrionales 
##                         1                         3                         1 
## 
## $HM
## 
##               o__B10-SB3A           o__Chlamydiales       o__Defluviicoccales 
##                         1                         1                         1 
##       o__Dethiobacterales              o__Dongiales                 o__DS-100 
##                         1                         1                         1 
##        o__Haloplasmatales         o__Incertae_Sedis        o__Kapabacteriales 
##                         1                         5                         1 
##    o__Paracaedibacterales          o__Peptococcales            o__Subgroup_13 
##                         1                         1                         1 
##           o__Sumerlaeales    o__Thermaerobacterales         o__Thermincolales 
##                         2                         1                         1 
## o__Thermoanaerobacterales           o__Unclassified     o__Vampirovibrionales 
##                         1                         1                         1 
## 
## $IT1
## 
##                    o__AKIW659             o__Anaerolineales 
##                             1                             1 
##                  o__B12-WMSP1                    o__Blfdi19 
##                             1                             1 
##         o__Christensenellales o__Clostridia_vadinBB60_group 
##                             1                             1 
##           o__Defluviicoccales                     o__DS-100 
##                             1                             1 
##                     o__DTU014           o__Fimbriimonadales 
##                             1                             1 
##           o__Flavobacteriales             o__Incertae_Sedis 
##                             1                             6 
##              o__Isosphaerales            o__Kapabacteriales 
##                             1                             1 
##          o__Obscuribacterales              o__Oligoflexales 
##                             1                             1 
##        o__Paracaedibacterales     o__Thermoanaerobacterales 
##                             1                             1 
##               o__Unclassified 
##                             2 
## 
## $GO1
## 
##               o__11-24            o__B10-SB3A        o__Chlamydiales 
##                      1                      1                      1 
##           o__Dongiales          o__Gemmatales     o__Haloplasmatales 
##                      1                      1                      1 
##        o__Holophagales        o__Holosporales      o__Incertae_Sedis 
##                      1                      1                      7 
##     o__Kapabacteriales     o__Leptolyngbyales          o__Lineage_IV 
##                      1                      1                      1 
##             o__mle1-27            o__MSB-4B10       o__Oligoflexales 
##                      2                      1                      1 
## o__Paracaedibacterales        o__Pirellulales       o__Rickettsiales 
##                      1                      1                      1 
##    o__Rickettsiellales              o__SJA-15     o__Sulfobacillales 
##                      1                      1                      1 
##  o__Thermacetogeniales        o__Unclassified 
##                      1                      3
```

``` r
save(Acc_Total_Merged_DA_Order, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_OrderLevel/Acc_Total_Merged_DA_Order.RData")
```

### Remove duplicat Order

``` r
Acc_NoDup_Merged_DA_Order <- lapply(Acc_Total_Merged_DA_Order, function(x) unique(x))
lapply(Acc_NoDup_Merged_DA_Order, function(x) nrow(x))
```

```
## $OH
## [1] 24
## 
## $DD
## [1] 26
## 
## $HE
## [1] 27
## 
## $KI
## [1] 29
## 
## $VL
## [1] 32
## 
## $CD
## [1] 27
## 
## $RI
## [1] 20
## 
## $KT
## [1] 25
## 
## $MC
## [1] 26
## 
## $HM
## [1] 22
## 
## $IT1
## [1] 25
## 
## $GO1
## [1] 31
```

``` r
lapply(Acc_NoDup_Merged_DA_Order, function(x) table(x$Order))
```

```
## $OH
## 
##     o__Actinomarinales             o__AKIW659             o__Blfdi19 
##                      1                      1                      1 
## o__Desulfitibacterales    o__Fimbriimonadales    o__Flavobacteriales 
##                      1                      1                      1 
##          o__Gemmatales     o__Haloplasmatales        o__Holophagales 
##                      1                      1                      1 
##      o__Incertae_Sedis          o__Lineage_IV     o__Micavibrionales 
##                      8                      1                      1 
##   o__Obscuribacterales          o__Opitutales    o__Rickettsiellales 
##                      1                      1                      1 
##     o__Sulfobacillales        o__Unclassified 
##                      1                      1 
## 
## $DD
## 
##                o__11-24              o__AKIW659             o__B10-SB3A 
##                       1                       1                       1 
##           o__Babeliales o__Caldicoprobacterales   o__Christensenellales 
##                       1                       1                       1 
##     o__Defluviicoccales         o__Holophagales       o__Incertae_Sedis 
##                       1                       1                       6 
##      o__Micavibrionales           o__Opitutales         o__Pirellulales 
##                       1                       1                       1 
##       o__Rhodothermales     o__Rickettsiellales      o__Sulfobacillales 
##                       1                       1                       1 
##         o__Sumerlaeales   o__Thermacetogeniales         o__Unclassified 
##                       1                       1                       4 
## 
## $HE
## 
##            o__AKIW659     o__Anaerolineales           o__B10-SB3A 
##                     1                     1                     1 
##         o__Babeliales      o__Caldilineales       o__Cytophagales 
##                     1                     1                     1 
##   o__Dethiobacterales          o__Dongiales   o__Fimbriimonadales 
##                     1                     1                     1 
##       o__Holophagales       o__Holosporales     o__Incertae_Sedis 
##                     1                     1                     7 
##      o__Legionellales    o__Methylococcales    o__Micavibrionales 
##                     1                     1                     1 
##         o__Opitutales              o__PeM15    o__Rhodobacterales 
##                     1                     1                     1 
##       o__Sumerlaeales o__Syntrophomonadales       o__Unclassified 
##                     1                     1                     1 
## 
## $KI
## 
##                o__Actinomarinales                        o__AKIW659 
##                                 1                                 1 
##                 o__Anaerolineales                      o__B12-WMSP1 
##                                 1                                 1 
##               o__Defluviicoccales                         o__DTU014 
##                                 1                                 1 
##                     o__Gemmatales                 o__Incertae_Sedis 
##                                 1                                 7 
##                  o__Isosphaerales                o__Kapabacteriales 
##                                 1                                 1 
##                   o__Kryptoniales                     o__Lineage_IV 
##                                 1                                 1 
##                o__Methylococcales                o__Micavibrionales 
##                                 1                                 1 
##                     o__Opitutales                o__Oscillospirales 
##                                 1                                 1 
##            o__Paracaedibacterales               o__Rickettsiellales 
##                                 1                                 1 
##                    o__Subgroup_13                   o__Sumerlaeales 
##                                 1                                 1 
##                   o__Unclassified o__Veillonellales-Selenomonadales 
##                                 2                                 1 
## 
## $VL
## 
##  o__Aggregatilineales            o__AKIW659           o__B10-SB3A 
##                     1                     1                     1 
##         o__Babeliales      o__Caldilineales   o__Defluviicoccales 
##                     1                     1                     1 
##      o__Deinococcales          o__Dongiales             o__DS-100 
##                     1                     1                     1 
##             o__DTU014       o__Holophagales       o__Holosporales 
##                     1                     1                     1 
##     o__Incertae_Sedis         o__Lineage_IV    o__Methylococcales 
##                     6                     1                     1 
##    o__Micavibrionales           o__MSB-4B10      o__Oligoflexales 
##                     1                     1                     1 
##             o__PLTA13      o__Rickettsiales   o__Rickettsiellales 
##                     1                     1                     1 
## o__Syntrophomonadales     o__Thermincolales       o__Unclassified 
##                     1                     1                     3 
## o__Vampirovibrionales 
##                     1 
## 
## $CD
## 
##    o__Aggregatilineales             o__B10-SB3A o__Caldicoprobacterales 
##                       1                       1                       1 
##        o__Caldilineales        o__Deinococcales         o__Holosporales 
##                       1                       1                       1 
##       o__Incertae_Sedis        o__Kallotenuales      o__Kapabacteriales 
##                       7                       1                       1 
##      o__Micavibrionales             o__MSB-4B10    o__Obscuribacterales 
##                       1                       1                       1 
##           o__Opitutales                o__PeM15      o__Sulfobacillales 
##                       1                       1                       1 
##         o__Sumerlaeales   o__Syntrophomonadales   o__Thermacetogeniales 
##                       1                       1                       1 
##         o__Unclassified 
##                       3 
## 
## $RI
## 
##   o__Aggregatilineales             o__AKIW659    o__Dethiobacterales 
##                      1                      1                      1 
##              o__DTU014           o__Elev-1554     o__Haloplasmatales 
##                      1                      1                      1 
##        o__Holophagales        o__Holosporales      o__Incertae_Sedis 
##                      1                      1                      4 
##     o__Kapabacteriales               o__MBA03   o__Obscuribacterales 
##                      1                      1                      1 
##          o__Opitutales o__Paracaedibacterales               o__PeM15 
##                      1                      1                      1 
##    o__Rickettsiellales  o__Vampirovibrionales 
##                      1                      1 
## 
## $KT
## 
##   o__Aggregatilineales             o__AKIW659      o__Anaerolineales 
##                      1                      1                      1 
##  o__Christensenellales o__Desulfitibacterales              o__DS-100 
##                      1                      1                      1 
##       o__Geobacterales        o__Holophagales      o__Incertae_Sedis 
##                      1                      1                      7 
##       o__Isosphaerales       o__Oligoflexales          o__Opitutales 
##                      1                      1                      1 
##   o__Proteinivoracales    o__Rickettsiellales              o__SJA-15 
##                      1                      1                      1 
##         o__Subgroup_13        o__Unclassified 
##                      1                      3 
## 
## $MC
## 
##                  o__11-24        o__Actinomarinales             o__Babeliales 
##                         1                         1                         1 
##              o__Dongiales       o__Fimbriimonadales       o__Flavobacteriales 
##                         1                         1                         1 
##          o__Geobacterales           o__Holophagales         o__Incertae_Sedis 
##                         1                         1                         7 
##        o__Kapabacteriales             o__Opitutales    o__Paracaedibacterales 
##                         1                         1                         1 
##         o__Pedosphaerales       o__Rickettsiellales     o__Syntrophomonadales 
##                         1                         1                         1 
## o__Thermoanaerobacterales           o__Unclassified     o__Vampirovibrionales 
##                         1                         3                         1 
## 
## $HM
## 
##               o__B10-SB3A           o__Chlamydiales       o__Defluviicoccales 
##                         1                         1                         1 
##       o__Dethiobacterales              o__Dongiales                 o__DS-100 
##                         1                         1                         1 
##        o__Haloplasmatales         o__Incertae_Sedis        o__Kapabacteriales 
##                         1                         5                         1 
##    o__Paracaedibacterales          o__Peptococcales            o__Subgroup_13 
##                         1                         1                         1 
##           o__Sumerlaeales    o__Thermaerobacterales         o__Thermincolales 
##                         1                         1                         1 
## o__Thermoanaerobacterales           o__Unclassified     o__Vampirovibrionales 
##                         1                         1                         1 
## 
## $IT1
## 
##                    o__AKIW659             o__Anaerolineales 
##                             1                             1 
##                  o__B12-WMSP1                    o__Blfdi19 
##                             1                             1 
##         o__Christensenellales o__Clostridia_vadinBB60_group 
##                             1                             1 
##           o__Defluviicoccales                     o__DS-100 
##                             1                             1 
##                     o__DTU014           o__Fimbriimonadales 
##                             1                             1 
##           o__Flavobacteriales             o__Incertae_Sedis 
##                             1                             6 
##              o__Isosphaerales            o__Kapabacteriales 
##                             1                             1 
##          o__Obscuribacterales              o__Oligoflexales 
##                             1                             1 
##        o__Paracaedibacterales     o__Thermoanaerobacterales 
##                             1                             1 
##               o__Unclassified 
##                             2 
## 
## $GO1
## 
##               o__11-24            o__B10-SB3A        o__Chlamydiales 
##                      1                      1                      1 
##           o__Dongiales          o__Gemmatales     o__Haloplasmatales 
##                      1                      1                      1 
##        o__Holophagales        o__Holosporales      o__Incertae_Sedis 
##                      1                      1                      7 
##     o__Kapabacteriales     o__Leptolyngbyales          o__Lineage_IV 
##                      1                      1                      1 
##             o__mle1-27            o__MSB-4B10       o__Oligoflexales 
##                      1                      1                      1 
## o__Paracaedibacterales        o__Pirellulales       o__Rickettsiales 
##                      1                      1                      1 
##    o__Rickettsiellales              o__SJA-15     o__Sulfobacillales 
##                      1                      1                      1 
##  o__Thermacetogeniales        o__Unclassified 
##                      1                      3
```

``` r
save(Acc_NoDup_Merged_DA_Order, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_OrderLevel/Acc_NoDup_Merged_DA_Order.RData")
```

### DA Order occurring twice

``` r
# Filter to Order detected by at least 2 tests
Acc_TwoTimes_Merged_DA_Order <- lapply(Acc_Total_Merged_DA_Order, function(x) {
    dup <- x[duplicated(x) | duplicated(x, fromLast = TRUE), , drop = FALSE]
    unique(dup)
  })

# Number of DA Order occurring in two tests
lapply(Acc_TwoTimes_Merged_DA_Order, function(x) nrow(x))
```

```
## $OH
## [1] 0
## 
## $DD
## [1] 0
## 
## $HE
## [1] 1
## 
## $KI
## [1] 1
## 
## $VL
## [1] 0
## 
## $CD
## [1] 0
## 
## $RI
## [1] 0
## 
## $KT
## [1] 0
## 
## $MC
## [1] 1
## 
## $HM
## [1] 1
## 
## $IT1
## [1] 0
## 
## $GO1
## [1] 1
```

``` r
save(Acc_TwoTimes_Merged_DA_Order, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_OrderLevel/Acc_TwoTimes_Merged_DA_Order.RData")
```

``` r
# Clean environment
rm(list = ls())
```

## 4.3.4 Accession Class level
### Load DA Class

``` r
## Load outputs DESeq
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ClassLevel/sigtab_stddds2_Wald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ClassLevel/sigtab_stddds2_LRT_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ClassLevel/sigtab_zinbWald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ClassLevel/sigtab_zinbLRT_Acc_Bac.RData")

## Load output Ancom
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ClassLevel/fdr_ancomWZ_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ClassLevel/fdr_ancomNoZ_Acc_Bac.RData")

## Load Phyloseq object
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_unnormalized_bac_ps_ForDA.RData")
```

### Extract DA Class
#### Extracting data from tests

``` r
# From each data frame, extract the Class or ASV IDs caught by each test
Acc_Bac_DESeqLRT_ASV_Class <- lapply(sigtab_stddds2_LRT_Acc_Bac, function(x) rownames(x))
Acc_Bac_DESeqWald_ASV_Class <- lapply(sigtab_stddds2_Wald_Acc_Bac, function(x) rownames(x))
Acc_Bac_zinbLRT_ASV_Class <- lapply(sigtab_zinbLRT_Acc_Bac, function(x) rownames(x))
Acc_Bac_zinbWald_ASV_Class <- lapply(sigtab_zinbWald_Acc_Bac, function(x) rownames(x))
Acc_Bac_ancomWZ_Class <- lapply(fdr_ancomWZ_Acc_Bac, function(x) as.vector(x[,'Species']))
Acc_Bac_ancomNoZ_Class <- lapply(fdr_ancomNoZ_Acc_Bac, function(x) as.vector(x[,'Species']))
```

#### Merging tests and transposing list

``` r
Acc_DESeq_DA_ASV_Class_l <- list(DESeqLRT = Acc_Bac_DESeqLRT_ASV_Class, DESeqWald = Acc_Bac_DESeqWald_ASV_Class, zinbLRT = Acc_Bac_zinbLRT_ASV_Class, zinbWald = Acc_Bac_zinbWald_ASV_Class)

# transpose list
Acc <- names(Acc_DESeq_DA_ASV_Class_l[[1]])
Acc_DESeq_DA_ASV_Class_l_t <- vector("list", length(Acc))
names(Acc_DESeq_DA_ASV_Class_l_t) <- Acc

for (cd in Acc) {
  Acc_DESeq_DA_ASV_Class_l_t[[cd]] <- lapply(Acc_DESeq_DA_ASV_Class_l, function(m) m[[cd]])
}

## DA Class Ancom
Acc_Ancom_DA_Class_l <- list(ancomWZ = Acc_Bac_ancomWZ_Class, ancomNoZ = Acc_Bac_ancomNoZ_Class)

# transpose list
Acc_Ancom_DA_Class_l_t <- vector("list", length(Acc))
names(Acc_Ancom_DA_Class_l_t) <- Acc

for (cd in Acc) {
  Acc_Ancom_DA_Class_l_t[[cd]] <- lapply(Acc_Ancom_DA_Class_l, function(m) m[[cd]])
}
```

### Extract taxomomy
#### DESeq

``` r
tax <- as.data.frame(tax_table(CP_unnormalized_bac_ps_ForDA))
tax$ASV <- rownames(tax)

Acc_NoDup_DESeq_DA_ASV_Class_Tax <- lapply(Acc_DESeq_DA_ASV_Class_l_t, function(x) {
  lapply(x, function(y) {
    unique(tax[tax$ASV %in% y, c("Class", "Phylum")])
    })})
```

#### Acnom

``` r
#lapply(Acc_Ancom_DA_Class_l, function(x) length(x))

Acc_NoDup_Ancom_DA_Class_tax <- lapply(Acc_Ancom_DA_Class_l_t, function(y) {
  lapply(y, function(x) {
    # Take unique names in case of duplicates
    DA_tax <- unique(x)
    
    # Split data between long and short taxonomic names
    full_tax_strings   <- DA_tax[grepl("k__", DA_tax)]
    family_only_strings <- DA_tax[!grepl("k__", DA_tax)]
    
    tax_full <- tibble(Taxon = full_tax_strings) %>%
      mutate(
        Phylum = str_extract(Taxon, "p__.*?(?=_c__|$)"),
        Class  = str_extract(Taxon, "c__.*?(?=_o__|$)")
      ) %>%
      select(Phylum, Class)
    
    # Extract taxonomy from ps object for short names
    tax_ps <- as.data.frame(tax_table(CP_unnormalized_bac_ps_ForDA))
    #family_only <- sub("^f__", "", family_only_strings)
    
    tax_family_map <- unique(tax_ps[tax_ps$Class %in% family_only_strings, c("Class", "Phylum")])
    #nrow(tax_family_map)
    
    # Remove left over Incertae Sedis from short names object
    tax_family_map <- tax_family_map %>%
      filter(Class != "f__Incertae_Sedis")
   
    # Combine long and short names
    tax_combined <- as.data.frame(bind_rows(tax_full, tax_family_map))
})})

#lapply(Acc_NoDup_Ancom_DA_Class_tax, function(x) nrow(x))
```

#### Merge Ancom and DESeq

``` r
Acc_NoDup_DESeq_DA_ASV_Class_Tax_m <- lapply(Acc_NoDup_DESeq_DA_ASV_Class_Tax, function(x) do.call(rbind, x))
Acc_NoDup_Ancom_DA_Class_tax_m <- lapply(Acc_NoDup_Ancom_DA_Class_tax, function(x) do.call(rbind, x))

# Merge and filter on unique family
Acc_Total_Merged_DA_Class <- mapply(function(y,z){
  input <- rbind(y,z)},
  y = Acc_NoDup_DESeq_DA_ASV_Class_Tax_m,
  z = Acc_NoDup_Ancom_DA_Class_tax_m,
  SIMPLIFY = FALSE)

lapply(Acc_Total_Merged_DA_Class, function(x) nrow(x))
```

```
## $OH
## [1] 12
## 
## $DD
## [1] 12
## 
## $HE
## [1] 14
## 
## $KI
## [1] 13
## 
## $VL
## [1] 12
## 
## $CD
## [1] 15
## 
## $RI
## [1] 6
## 
## $KT
## [1] 7
## 
## $MC
## [1] 14
## 
## $HM
## [1] 15
## 
## $IT1
## [1] 9
## 
## $GO1
## [1] 14
```

``` r
lapply(Acc_Total_Merged_DA_Class, function(x) table(x$Class))
```

```
## $OH
## 
##    c__Elusimicrobia   c__Fimbriimonadia   c__Incertae_Sedis        c__Moorellia 
##                   1                   1                   4                   1 
##    c__Parcubacteria      c__Subgroup_22      c__Subgroup_25    c__Sulfobacillia 
##                   1                   1                   1                   1 
## c__Vampirivibrionia 
##                   1 
## 
## $DD
## 
##         c__Babeliae      c__bacteriap25   c__Incertae_Sedis            c__OLB14 
##                   1                   1                   3                   1 
##     c__Rhodothermia           c__SHA-26    c__Sulfobacillia       c__Sumerlaeia 
##                   1                   1                   1                   1 
##  c__Thermacetogenia c__Vampirivibrionia 
##                   1                   1 
## 
## $HE
## 
##                 c__Babeliae c__BD2-11_terrestrial_group 
##                           1                           1 
##           c__Dethiobacteria           c__Eremiobacteria 
##                           1                           1 
##           c__Fimbriimonadia          c__Gracilibacteria 
##                           1                           1 
##           c__Incertae_Sedis                   c__SHA-26 
##                           1                           1 
##              c__Subgroup_22               c__Subgroup_5 
##                           1                           1 
##               c__Sumerlaeia         c__Syntrophomonadia 
##                           1                           1 
##             c__Unclassified         c__Vampirivibrionia 
##                           1                           1 
## 
## $KI
## 
##  c__Elusimicrobia c__Eremiobacteria c__Incertae_Sedis   c__Kapabacteria 
##                 1                 1                 4                 1 
##      c__Kryptonia  c__Negativicutes         c__P2-11E c__Planctomycetes 
##                 1                 1                 1                 1 
##         c__SHA-26     c__Sumerlaeia 
##                 1                 1 
## 
## $VL
## 
##         c__Babeliae       c__Deinococci c__Desulfuromonadia    c__Elusimicrobia 
##                   1                   1                   1                   1 
##   c__Incertae_Sedis           c__P2-11E       c__Subgroup_5 c__Syntrophomonadia 
##                   4                   1                   1                   1 
## c__Vampirivibrionia 
##                   1 
## 
## $CD
## 
##       c__Deinococci   c__Incertae_Sedis     c__Kapabacteria            c__OLB14 
##                   1                   2                   1                   1 
##           c__SHA-26      c__Subgroup_22      c__Subgroup_25       c__Subgroup_5 
##                   1                   1                   1                   1 
##    c__Sulfobacillia       c__Sumerlaeia c__Syntrophomonadia  c__Thermacetogenia 
##                   1                   1                   1                   1 
##     c__Unclassified c__Vampirivibrionia 
##                   1                   1 
## 
## $RI
## 
## c__Dethiobacteria c__Incertae_Sedis   c__Kapabacteria         c__SHA-26 
##                 1                 2                 1                 1 
##    c__Subgroup_25 
##                 1 
## 
## $KT
## 
## c__Gracilibacteria  c__Incertae_Sedis       c__Moorellia           c__OM190 
##                  1                  2                  1                  1 
##          c__P2-11E   c__Parcubacteria 
##                  1                  1 
## 
## $MC
## 
##             c__Babeliae     c__Desulfuromonadia       c__Eremiobacteria 
##                       1                       1                       1 
##       c__Fimbriimonadia       c__Incertae_Sedis         c__Kapabacteria 
##                       1                       3                       1 
##                c__OM190               c__P2-11E        c__Parcubacteria 
##                       1                       1                       1 
##     c__Syntrophomonadia c__Thermoanaerobacteria         c__Unclassified 
##                       1                       1                       1 
## 
## $HM
## 
## c__BD2-11_terrestrial_group               c__Chlamydiia 
##                           1                           1 
##         c__Desulfuromonadia           c__Dethiobacteria 
##                           1                           1 
##           c__Incertae_Sedis             c__Kapabacteria 
##                           4                           1 
##            c__Parcubacteria                   c__SHA-26 
##                           1                           1 
##               c__Sumerlaeia        c__Thermaerobacteria 
##                           2                           1 
##     c__Thermoanaerobacteria 
##                           1 
## 
## $IT1
## 
##       c__Eremiobacteria       c__Fimbriimonadia      c__Gracilibacteria 
##                       1                       1                       1 
##       c__Incertae_Sedis         c__Kapabacteria               c__P2-11E 
##                       3                       1                       1 
## c__Thermoanaerobacteria 
##                       1 
## 
## $GO1
## 
##      c__Chlamydiia   c__Elusimicrobia c__Gracilibacteria  c__Incertae_Sedis 
##                  1                  1                  1                  2 
##    c__Kapabacteria          c__P2-11E          c__SHA-26     c__Subgroup_22 
##                  1                  1                  1                  1 
##     c__Subgroup_25      c__Subgroup_5   c__Sulfobacillia c__Thermacetogenia 
##                  1                  1                  1                  1 
##    c__Unclassified 
##                  1
```

``` r
save(Acc_Total_Merged_DA_Class, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ClassLevel/Acc_Total_Merged_DA_Class.RData")
```
 
### Remove duplicat Class

``` r
Acc_NoDup_Merged_DA_Class <- lapply(Acc_Total_Merged_DA_Class, function(x) unique(x))
lapply(Acc_NoDup_Merged_DA_Class, function(x) nrow(x))
```

```
## $OH
## [1] 12
## 
## $DD
## [1] 12
## 
## $HE
## [1] 14
## 
## $KI
## [1] 13
## 
## $VL
## [1] 12
## 
## $CD
## [1] 15
## 
## $RI
## [1] 6
## 
## $KT
## [1] 7
## 
## $MC
## [1] 14
## 
## $HM
## [1] 14
## 
## $IT1
## [1] 9
## 
## $GO1
## [1] 14
```

``` r
lapply(Acc_NoDup_Merged_DA_Class, function(x) table(x$Class))
```

```
## $OH
## 
##    c__Elusimicrobia   c__Fimbriimonadia   c__Incertae_Sedis        c__Moorellia 
##                   1                   1                   4                   1 
##    c__Parcubacteria      c__Subgroup_22      c__Subgroup_25    c__Sulfobacillia 
##                   1                   1                   1                   1 
## c__Vampirivibrionia 
##                   1 
## 
## $DD
## 
##         c__Babeliae      c__bacteriap25   c__Incertae_Sedis            c__OLB14 
##                   1                   1                   3                   1 
##     c__Rhodothermia           c__SHA-26    c__Sulfobacillia       c__Sumerlaeia 
##                   1                   1                   1                   1 
##  c__Thermacetogenia c__Vampirivibrionia 
##                   1                   1 
## 
## $HE
## 
##                 c__Babeliae c__BD2-11_terrestrial_group 
##                           1                           1 
##           c__Dethiobacteria           c__Eremiobacteria 
##                           1                           1 
##           c__Fimbriimonadia          c__Gracilibacteria 
##                           1                           1 
##           c__Incertae_Sedis                   c__SHA-26 
##                           1                           1 
##              c__Subgroup_22               c__Subgroup_5 
##                           1                           1 
##               c__Sumerlaeia         c__Syntrophomonadia 
##                           1                           1 
##             c__Unclassified         c__Vampirivibrionia 
##                           1                           1 
## 
## $KI
## 
##  c__Elusimicrobia c__Eremiobacteria c__Incertae_Sedis   c__Kapabacteria 
##                 1                 1                 4                 1 
##      c__Kryptonia  c__Negativicutes         c__P2-11E c__Planctomycetes 
##                 1                 1                 1                 1 
##         c__SHA-26     c__Sumerlaeia 
##                 1                 1 
## 
## $VL
## 
##         c__Babeliae       c__Deinococci c__Desulfuromonadia    c__Elusimicrobia 
##                   1                   1                   1                   1 
##   c__Incertae_Sedis           c__P2-11E       c__Subgroup_5 c__Syntrophomonadia 
##                   4                   1                   1                   1 
## c__Vampirivibrionia 
##                   1 
## 
## $CD
## 
##       c__Deinococci   c__Incertae_Sedis     c__Kapabacteria            c__OLB14 
##                   1                   2                   1                   1 
##           c__SHA-26      c__Subgroup_22      c__Subgroup_25       c__Subgroup_5 
##                   1                   1                   1                   1 
##    c__Sulfobacillia       c__Sumerlaeia c__Syntrophomonadia  c__Thermacetogenia 
##                   1                   1                   1                   1 
##     c__Unclassified c__Vampirivibrionia 
##                   1                   1 
## 
## $RI
## 
## c__Dethiobacteria c__Incertae_Sedis   c__Kapabacteria         c__SHA-26 
##                 1                 2                 1                 1 
##    c__Subgroup_25 
##                 1 
## 
## $KT
## 
## c__Gracilibacteria  c__Incertae_Sedis       c__Moorellia           c__OM190 
##                  1                  2                  1                  1 
##          c__P2-11E   c__Parcubacteria 
##                  1                  1 
## 
## $MC
## 
##             c__Babeliae     c__Desulfuromonadia       c__Eremiobacteria 
##                       1                       1                       1 
##       c__Fimbriimonadia       c__Incertae_Sedis         c__Kapabacteria 
##                       1                       3                       1 
##                c__OM190               c__P2-11E        c__Parcubacteria 
##                       1                       1                       1 
##     c__Syntrophomonadia c__Thermoanaerobacteria         c__Unclassified 
##                       1                       1                       1 
## 
## $HM
## 
## c__BD2-11_terrestrial_group               c__Chlamydiia 
##                           1                           1 
##         c__Desulfuromonadia           c__Dethiobacteria 
##                           1                           1 
##           c__Incertae_Sedis             c__Kapabacteria 
##                           4                           1 
##            c__Parcubacteria                   c__SHA-26 
##                           1                           1 
##               c__Sumerlaeia        c__Thermaerobacteria 
##                           1                           1 
##     c__Thermoanaerobacteria 
##                           1 
## 
## $IT1
## 
##       c__Eremiobacteria       c__Fimbriimonadia      c__Gracilibacteria 
##                       1                       1                       1 
##       c__Incertae_Sedis         c__Kapabacteria               c__P2-11E 
##                       3                       1                       1 
## c__Thermoanaerobacteria 
##                       1 
## 
## $GO1
## 
##      c__Chlamydiia   c__Elusimicrobia c__Gracilibacteria  c__Incertae_Sedis 
##                  1                  1                  1                  2 
##    c__Kapabacteria          c__P2-11E          c__SHA-26     c__Subgroup_22 
##                  1                  1                  1                  1 
##     c__Subgroup_25      c__Subgroup_5   c__Sulfobacillia c__Thermacetogenia 
##                  1                  1                  1                  1 
##    c__Unclassified 
##                  1
```

``` r
save(Acc_NoDup_Merged_DA_Class, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ClassLevel/Acc_NoDup_Merged_DA_Class.RData")
```

### DA Class occurring twice

``` r
# Filter to Class detected by at least 2 tests
Acc_TwoTimes_Merged_DA_Class <- lapply(Acc_Total_Merged_DA_Class, function(x) {
    dup <- x[duplicated(x) | duplicated(x, fromLast = TRUE), , drop = FALSE]
    unique(dup)
  })

# Number of DA Class occurring in two tests
lapply(Acc_TwoTimes_Merged_DA_Class, function(x) nrow(x))
```

```
## $OH
## [1] 0
## 
## $DD
## [1] 0
## 
## $HE
## [1] 0
## 
## $KI
## [1] 0
## 
## $VL
## [1] 0
## 
## $CD
## [1] 0
## 
## $RI
## [1] 0
## 
## $KT
## [1] 0
## 
## $MC
## [1] 0
## 
## $HM
## [1] 1
## 
## $IT1
## [1] 0
## 
## $GO1
## [1] 0
```

``` r
save(Acc_TwoTimes_Merged_DA_Class, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ClassLevel/Acc_TwoTimes_Merged_DA_Class.RData")
```

``` r
# Clean environment
rm(list = ls())
```


## 4.3.5 Accession Phylum level
### Load DA Phylum

``` r
## Load outputs DESeq
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_PhylumLevel/sigtab_stddds2_Wald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_PhylumLevel/sigtab_stddds2_LRT_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_PhylumLevel/sigtab_zinbWald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_PhylumLevel/sigtab_zinbLRT_Acc_Bac.RData")

## Load output Ancom
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_PhylumLevel/fdr_ancomWZ_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_PhylumLevel/fdr_ancomNoZ_Acc_Bac.RData")

## Load Phyloseq object
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_unnormalized_bac_ps_ForDA.RData")
```

### Extract DA Phylum
#### Extracting data from tests

``` r
# From each data frame, extract the Phylum or ASV IDs caught by each test
Acc_Bac_DESeqLRT_ASV_Phylum <- lapply(sigtab_stddds2_LRT_Acc_Bac, function(x) rownames(x))
Acc_Bac_DESeqWald_ASV_Phylum <- lapply(sigtab_stddds2_Wald_Acc_Bac, function(x) rownames(x))
Acc_Bac_zinbLRT_ASV_Phylum <- lapply(sigtab_zinbLRT_Acc_Bac, function(x) rownames(x))
Acc_Bac_zinbWald_ASV_Phylum <- lapply(sigtab_zinbWald_Acc_Bac, function(x) rownames(x))
Acc_Bac_ancomWZ_Phylum <- lapply(fdr_ancomWZ_Acc_Bac, function(x) as.vector(x[,'Species']))
Acc_Bac_ancomNoZ_Phylum <- lapply(fdr_ancomNoZ_Acc_Bac, function(x) as.vector(x[,'Species']))
```

#### Merging tests and transposing list

``` r
Acc_DESeq_DA_ASV_Phylum_l <- list(DESeqLRT = Acc_Bac_DESeqLRT_ASV_Phylum, DESeqWald = Acc_Bac_DESeqWald_ASV_Phylum, zinbLRT = Acc_Bac_zinbLRT_ASV_Phylum, zinbWald = Acc_Bac_zinbWald_ASV_Phylum)

# transpose list
Acc <- names(Acc_DESeq_DA_ASV_Phylum_l[[1]])
Acc_DESeq_DA_ASV_Phylum_l_t <- vector("list", length(Acc))
names(Acc_DESeq_DA_ASV_Phylum_l_t) <- Acc

for (cd in Acc) {
  Acc_DESeq_DA_ASV_Phylum_l_t[[cd]] <- lapply(Acc_DESeq_DA_ASV_Phylum_l, function(m) m[[cd]])
}

## DA Phylum Ancom
Acc_Ancom_DA_Phylum_l <- list(ancomWZ = Acc_Bac_ancomWZ_Phylum, ancomNoZ = Acc_Bac_ancomNoZ_Phylum)

# transpose list
Acc_Ancom_DA_Phylum_l_t <- vector("list", length(Acc))
names(Acc_Ancom_DA_Phylum_l_t) <- Acc

for (cd in Acc) {
  Acc_Ancom_DA_Phylum_l_t[[cd]] <- lapply(Acc_Ancom_DA_Phylum_l, function(m) m[[cd]])
}
```

### Extract taxomomy
#### DESeq

``` r
tax <- as.data.frame(tax_table(CP_unnormalized_bac_ps_ForDA))
tax$ASV <- rownames(tax)

Acc_NoDup_DESeq_DA_ASV_Phylum_Tax <- lapply(Acc_DESeq_DA_ASV_Phylum_l_t, function(x) {
  lapply(x, function(y) {
    unique(tax[tax$ASV %in% y, c("Phylum")])
    })})
```

#### Acnom

``` r
#lapply(Acc_Ancom_DA_Phylum_l, function(x) length(x))

Acc_NoDup_Ancom_DA_Phylum_tax <- lapply(Acc_Ancom_DA_Phylum_l_t, function(y) {
  lapply(y, function(x) {
    # Take unique names in case of duplicates
    DA_tax <- unique(x)
    
    # Split data between long and short taxonomic names
    full_tax_strings   <- DA_tax[grepl("k__", DA_tax)]
    family_only_strings <- DA_tax[!grepl("k__", DA_tax)]
    
    # Extract taxomomy from long names
    tax_full <- tibble(Taxon = full_tax_strings) %>%
      mutate(
        Phylum = str_extract(Taxon, "p__.*?(?=_c__|$)"),
        Class  = str_extract(Taxon, "c__.*?(?=_o__|$)")
      ) %>%
      select(Phylum, Class)
    
    # Extract taxonomy from ps object for short names
    tax_ps <- as.data.frame(tax_table(CP_unnormalized_bac_ps_ForDA))
    #family_only <- sub("^f__", "", family_only_strings)
    
    tax_family_map <- unique(tax_ps[tax_ps$Phylum %in% family_only_strings, c("Class", "Phylum")])
    #nrow(tax_family_map)
    
    # Remove left over Incertae Sedis from short names object
    tax_family_map <- tax_family_map %>%
      filter(Phylum != "f__Incertae_Sedis")
   
    # Combine long and short names
    tax_combined <- as.data.frame(bind_rows(tax_full, tax_family_map))
    tax_combined <- tax_combined$Phylum
})})

#lapply(Acc_NoDup_Ancom_DA_Phylum_tax, function(x) nrow(x))
```

#### Merge Ancom and DESeq

``` r
Acc_NoDup_DESeq_DA_ASV_Phylum_Tax_m <- lapply(Acc_NoDup_DESeq_DA_ASV_Phylum_Tax, function(x) unlist(x, use.names = FALSE))
Acc_NoDup_Ancom_DA_Phylum_tax_m <- lapply(Acc_NoDup_Ancom_DA_Phylum_tax, function(x) unlist(x, use.names = FALSE))

# Merge and filter on unique family
Acc_Total_Merged_DA_Phylum <- mapply(function(y,z){
  input <- c(y,z)},
  y = Acc_NoDup_DESeq_DA_ASV_Phylum_Tax_m,
  z = Acc_NoDup_Ancom_DA_Phylum_tax_m,
  SIMPLIFY = FALSE)

lapply(Acc_Total_Merged_DA_Phylum, function(x) length(x))
```

```
## $OH
## [1] 6
## 
## $DD
## [1] 6
## 
## $HE
## [1] 7
## 
## $KI
## [1] 6
## 
## $VL
## [1] 5
## 
## $CD
## [1] 4
## 
## $RI
## [1] 2
## 
## $KT
## [1] 1
## 
## $MC
## [1] 17
## 
## $HM
## [1] 6
## 
## $IT1
## [1] 7
## 
## $GO1
## [1] 7
```

``` r
lapply(Acc_Total_Merged_DA_Phylum, function(x) table(x))
```

```
## $OH
## x
##              p__Elusimicrobiota                       p__MBNT15 
##                               3                               1 
##                        p__NB1-j p__SAR324_clade(Marine_group_B) 
##                               1                               1 
## 
## $DD
## x
##                 p__Dependentiae                       p__MBNT15 
##                               1                               1 
##                        p__NB1-j                p__Rhodothermota 
##                               1                               1 
## p__SAR324_clade(Marine_group_B)                  p__Sumerlaeota 
##                               1                               1 
## 
## $HE
## x
##             p__Armatimonadota p__Candidatus_Eremiobacterota 
##                             3                             1 
##               p__Dependentiae                p__Sumerlaeota 
##                             1                             1 
##                        p__WS2 
##                             1 
## 
## $KI
## x
## p__Candidatus_Eremiobacterota    p__Candidatus_Kapabacteria 
##                             1                             1 
##       p__Candidatus_Kryptonia                     p__MBNT15 
##                             1                             1 
##                      p__NB1-j                p__Sumerlaeota 
##                             1                             1 
## 
## $VL
## x
## p__Acidobacteriota   p__Chloroflexota    p__Deinococcota    p__Dependentiae 
##                  1                  1                  1                  1 
##           p__NB1-j 
##                  1 
## 
## $CD
## x
## p__Candidatus_Kapabacteria            p__Deinococcota 
##                          1                          1 
##             p__Sumerlaeota                     p__WS2 
##                          1                          1 
## 
## $RI
## x
## p__Candidatus_Kapabacteria                   p__NB1-j 
##                          1                          1 
## 
## $KT
## x
## p__WS2 
##      1 
## 
## $MC
## x
##   p__Candidatus_Eremiobacterota      p__Candidatus_Kapabacteria 
##                               1                               1 
##                 p__Dependentiae                        p__NB1-j 
##                               1                               1 
##              p__Planctomycetota p__SAR324_clade(Marine_group_B) 
##                              12                               1 
## 
## $HM
## x
## p__Candidatus_Kapabacteria             p__Chlamydiota 
##                          1                          1 
##                   p__NB1-j             p__Sumerlaeota 
##                          1                          2 
##                     p__WS2 
##                          1 
## 
## $IT1
## x
##             p__Armatimonadota p__Candidatus_Eremiobacterota 
##                             3                             1 
##    p__Candidatus_Kapabacteria                     p__MBNT15 
##                             1                             1 
##                      p__NB1-j 
##                             1 
## 
## $GO1
## x
## p__Candidatus_Kapabacteria             p__Chlamydiota 
##                          1                          1 
##         p__Elusimicrobiota                   p__NB1-j 
##                          3                          1 
##                     p__WS2 
##                          1
```

``` r
save(Acc_Total_Merged_DA_Phylum, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_PhylumLevel/Acc_Total_Merged_DA_Phylum.RData")
```

### Remove duplicat Phylum

``` r
Acc_NoDup_Merged_DA_Phylum <- lapply(Acc_Total_Merged_DA_Phylum, function(x) unique(x))
lapply(Acc_NoDup_Merged_DA_Phylum, function(x) length(x))
```

```
## $OH
## [1] 4
## 
## $DD
## [1] 6
## 
## $HE
## [1] 5
## 
## $KI
## [1] 6
## 
## $VL
## [1] 5
## 
## $CD
## [1] 4
## 
## $RI
## [1] 2
## 
## $KT
## [1] 1
## 
## $MC
## [1] 6
## 
## $HM
## [1] 5
## 
## $IT1
## [1] 5
## 
## $GO1
## [1] 5
```

``` r
lapply(Acc_NoDup_Merged_DA_Phylum, function(x) table(x))
```

```
## $OH
## x
##              p__Elusimicrobiota                       p__MBNT15 
##                               1                               1 
##                        p__NB1-j p__SAR324_clade(Marine_group_B) 
##                               1                               1 
## 
## $DD
## x
##                 p__Dependentiae                       p__MBNT15 
##                               1                               1 
##                        p__NB1-j                p__Rhodothermota 
##                               1                               1 
## p__SAR324_clade(Marine_group_B)                  p__Sumerlaeota 
##                               1                               1 
## 
## $HE
## x
##             p__Armatimonadota p__Candidatus_Eremiobacterota 
##                             1                             1 
##               p__Dependentiae                p__Sumerlaeota 
##                             1                             1 
##                        p__WS2 
##                             1 
## 
## $KI
## x
## p__Candidatus_Eremiobacterota    p__Candidatus_Kapabacteria 
##                             1                             1 
##       p__Candidatus_Kryptonia                     p__MBNT15 
##                             1                             1 
##                      p__NB1-j                p__Sumerlaeota 
##                             1                             1 
## 
## $VL
## x
## p__Acidobacteriota   p__Chloroflexota    p__Deinococcota    p__Dependentiae 
##                  1                  1                  1                  1 
##           p__NB1-j 
##                  1 
## 
## $CD
## x
## p__Candidatus_Kapabacteria            p__Deinococcota 
##                          1                          1 
##             p__Sumerlaeota                     p__WS2 
##                          1                          1 
## 
## $RI
## x
## p__Candidatus_Kapabacteria                   p__NB1-j 
##                          1                          1 
## 
## $KT
## x
## p__WS2 
##      1 
## 
## $MC
## x
##   p__Candidatus_Eremiobacterota      p__Candidatus_Kapabacteria 
##                               1                               1 
##                 p__Dependentiae                        p__NB1-j 
##                               1                               1 
##              p__Planctomycetota p__SAR324_clade(Marine_group_B) 
##                               1                               1 
## 
## $HM
## x
## p__Candidatus_Kapabacteria             p__Chlamydiota 
##                          1                          1 
##                   p__NB1-j             p__Sumerlaeota 
##                          1                          1 
##                     p__WS2 
##                          1 
## 
## $IT1
## x
##             p__Armatimonadota p__Candidatus_Eremiobacterota 
##                             1                             1 
##    p__Candidatus_Kapabacteria                     p__MBNT15 
##                             1                             1 
##                      p__NB1-j 
##                             1 
## 
## $GO1
## x
## p__Candidatus_Kapabacteria             p__Chlamydiota 
##                          1                          1 
##         p__Elusimicrobiota                   p__NB1-j 
##                          1                          1 
##                     p__WS2 
##                          1
```

``` r
save(Acc_NoDup_Merged_DA_Phylum, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_PhylumLevel/Acc_NoDup_Merged_DA_Phylum.RData")
```

### DA Phylum occurring twice

``` r
# Filter to Phylum detected by at least 2 tests
Acc_TwoTimes_Merged_DA_Phylum <- lapply(Acc_Total_Merged_DA_Phylum, function(x) {
    dup <- x[duplicated(x) | duplicated(x, fromLast = TRUE), drop = FALSE]
    unique(dup)
  })

# Number of DA Phylum occurring in two tests
lapply(Acc_TwoTimes_Merged_DA_Phylum, function(x) length(x))
```

```
## $OH
## [1] 1
## 
## $DD
## [1] 0
## 
## $HE
## [1] 1
## 
## $KI
## [1] 0
## 
## $VL
## [1] 0
## 
## $CD
## [1] 0
## 
## $RI
## [1] 0
## 
## $KT
## [1] 0
## 
## $MC
## [1] 1
## 
## $HM
## [1] 1
## 
## $IT1
## [1] 1
## 
## $GO1
## [1] 1
```

``` r
save(Acc_TwoTimes_Merged_DA_Phylum, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_PhylumLevel/Acc_TwoTimes_Merged_DA_Phylum.RData")
```

``` r
# Clean environment
rm(list = ls())
```


## 4.3.4 Domestication ASV level
### Load DA ASVs

``` r
## Load outputs DESeq
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ASVLevel/sigtab_stddds2_Wald_Dom_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ASVLevel/sigtab_stddds2_LRT_Dom_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ASVLevel/sigtab_zinbWald_Dom_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ASVLevel/sigtab_zinbLRT_Dom_Bac.RData")

## Load output Ancom
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/fdr_ancomWZ_Dom_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/fdr_ancomNoZ_Dom_Bac.RData")
```

### Extract DA ASVs

``` r
## From each data frame, extract the ASV IDs caught by each test
Dom_Bac_DESeqLRT_ASVs <- lapply(sigtab_stddds2_LRT_Dom_Bac, function(x) rownames(x))
Dom_Bac_DESeqWald_ASVs <- lapply(sigtab_stddds2_Wald_Dom_Bac, function(x) rownames(x))
Dom_Bac_zinbLRT_ASVs <- lapply(sigtab_zinbLRT_Dom_Bac, function(x) rownames(x))
Dom_Bac_zinbWald_ASVs <- lapply(sigtab_zinbWald_Dom_Bac, function(x) rownames(x))
Dom_Bac_ancomWZ_ASVs <- lapply(fdr_ancomWZ_Dom_Bac, function(x) as.vector(x[,'Species']))
Dom_Bac_ancomNoZ_ASVs <- lapply(fdr_ancomNoZ_Dom_Bac, function(x) as.vector(x[,'Species']))

# Within each plant concatenate the identified asvs
Dom_Bac_Total_DA_ASV <- mapply(function(s,t,w,x,y,z){
  input <- c(s,t,w,x,y,z)},
  s = Dom_Bac_DESeqLRT_ASVs,
  t = Dom_Bac_DESeqWald_ASVs,
  w = Dom_Bac_zinbLRT_ASVs,
  x = Dom_Bac_zinbWald_ASVs,
  y = Dom_Bac_ancomWZ_ASVs,
  z = Dom_Bac_ancomNoZ_ASVs,
  SIMPLIFY = FALSE)
```

### Total number of ASVs

``` r
# Total number of DA ASVs
lapply(Dom_Bac_Total_DA_ASV, function(x) length(x))
```

```
## $Wild
## [1] 81
## 
## $Cultivated
## [1] 48
```

``` r
save(Dom_Bac_Total_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_Total_DA_ASV.RData")

# Remove duplicated asvs names from that list
Dom_Bac_NoDup_DA_ASV <- lapply(Dom_Bac_Total_DA_ASV, function(x)
  unique(x))

# Number of DA ASVs without duplicates
lapply(Dom_Bac_NoDup_DA_ASV, function(x) length(x))
```

```
## $Wild
## [1] 72
## 
## $Cultivated
## [1] 27
```

``` r
save(Dom_Bac_NoDup_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_NoDup_DA_ASV.RData")

# Extracting taxonomy
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_unnormalized_bac_ps.RData")

# Extract taxonomy of DA ASVs
CP_unnormalized_bac_ps_Tax <- as.data.frame(tax_table(CP_unnormalized_bac_ps))
Dom_Bac_NoDup_DA_ASV_Tax <- lapply(Dom_Bac_NoDup_DA_ASV, function(x) CP_unnormalized_bac_ps_Tax[x, ])

save(Dom_Bac_NoDup_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects//CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_NoDup_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
Dom_Bac_NoDup_DA_ASV_Tax_Table <- lapply(Dom_Bac_NoDup_DA_ASV_Tax, function(df) {
  lapply(df[ , 2:6], function(col) table(col, useNA = "ifany"))
})
Dom_Bac_NoDup_DA_ASV_Tax_Table <- lapply(Dom_Bac_NoDup_DA_ASV_Tax_Table, function(ranklist) {
                                        lapply(ranklist, as.data.frame)})

# Examples
Dom_Bac_NoDup_DA_ASV_Tax_Table$HM$Class
```

```
## NULL
```

``` r
Dom_Bac_NoDup_DA_ASV_Tax_Table$VL$Class
```

```
## NULL
```

``` r
save(Dom_Bac_NoDup_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_NoDup_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring twice

``` r
# Filter to ASVs detected by at least 2 tests
Dom_Bac_NumbOcc_DA_ASV <- lapply(Dom_Bac_Total_DA_ASV, function(x) table(unlist(x)))
Dom_Bac_TwoTimes_DA_ASV <- lapply(Dom_Bac_NumbOcc_DA_ASV, function(x) names(x[x >= 2]))
Dom_Bac_TwoTimes_DA_ASV
```

```
## $Wild
## [1] "bASV_1025" "bASV_1389" "bASV_269"  "bASV_287"  "bASV_581"  "bASV_601" 
## [7] "bASV_656"  "bASV_724"  "bASV_963" 
## 
## $Cultivated
##  [1] "bASV_1083" "bASV_1272" "bASV_1298" "bASV_1490" "bASV_1701" "bASV_1760"
##  [7] "bASV_2032" "bASV_2402" "bASV_2416" "bASV_2677" "bASV_269"  "bASV_287" 
## [13] "bASV_387"  "bASV_586"  "bASV_601"  "bASV_683"  "bASV_724"  "bASV_878" 
## [19] "bASV_911"  "bASV_938"  "bASV_963"
```

``` r
# Number of DA ASVs occurring in two tests
lapply(Dom_Bac_TwoTimes_DA_ASV, function(x) length(x))
```

```
## $Wild
## [1] 9
## 
## $Cultivated
## [1] 21
```

``` r
save(Dom_Bac_TwoTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_TwoTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
Dom_Bac_TwoTimes_DA_ASV_Tax <- lapply(Dom_Bac_TwoTimes_DA_ASV, function(x) CP_unnormalized_bac_ps_Tax[x, ])

save(Dom_Bac_TwoTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_TwoTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
Dom_Bac_TwoTimes_DA_ASV_Tax_Table <- lapply(Dom_Bac_TwoTimes_DA_ASV_Tax, function(df) {
  lapply(df[ , 2:6], function(col) table(col, useNA = "ifany"))
})
Dom_Bac_TwoTimes_DA_ASV_Tax_Table <- lapply(Dom_Bac_TwoTimes_DA_ASV_Tax_Table, function(ranklist) {
                                            lapply(ranklist, as.data.frame)})

# Examples
Dom_Bac_TwoTimes_DA_ASV_Tax_Table$HM$Class
```

```
## NULL
```

``` r
Dom_Bac_TwoTimes_DA_ASV_Tax_Table$VL$Class
```

```
## NULL
```

``` r
save(Dom_Bac_TwoTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_TwoTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring three times

``` r
# Filter to ASVs detected by at least 3 tests
Dom_Bac_NumbOcc_DA_ASV <- lapply(Dom_Bac_Total_DA_ASV, function(x) table(unlist(x)))
Dom_Bac_ThreeTimes_DA_ASV <- lapply(Dom_Bac_NumbOcc_DA_ASV, function(x) names(x[x >= 3]))
Dom_Bac_ThreeTimes_DA_ASV
```

```
## $Wild
## character(0)
## 
## $Cultivated
## character(0)
```

``` r
# Number of DA ASVs occurring in three tests
lapply(Dom_Bac_ThreeTimes_DA_ASV, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_ThreeTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_ThreeTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
Dom_Bac_ThreeTimes_DA_ASV_Tax <- lapply(Dom_Bac_ThreeTimes_DA_ASV, function(x) CP_unnormalized_bac_ps_Tax[x, ])

save(Dom_Bac_ThreeTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_ThreeTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
Dom_Bac_ThreeTimes_DA_ASV_Tax_Table <- lapply(Dom_Bac_ThreeTimes_DA_ASV_Tax, function(df) {
  lapply(df[ , 2:6], function(col) table(col, useNA = "ifany"))
})
Dom_Bac_ThreeTimes_DA_ASV_Tax_Table <- lapply(Dom_Bac_ThreeTimes_DA_ASV_Tax_Table, function(ranklist) {
                                            lapply(ranklist, as.data.frame)})

# Examples
Dom_Bac_ThreeTimes_DA_ASV_Tax_Table$HM$Class
```

```
## NULL
```

``` r
Dom_Bac_ThreeTimes_DA_ASV_Tax_Table$VL$Class
```

```
## NULL
```

``` r
save(Dom_Bac_ThreeTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_ThreeTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring four times

``` r
# Filter to ASVs detected by at least 4 tests
Dom_Bac_FourTimes_DA_ASV <- lapply(Dom_Bac_NumbOcc_DA_ASV, function(x) names(x[x >= 4]))
Dom_Bac_FourTimes_DA_ASV
```

```
## $Wild
## character(0)
## 
## $Cultivated
## character(0)
```

``` r
# Number of DA ASVs occurring in four tests
lapply(Dom_Bac_FourTimes_DA_ASV, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_FourTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_FourTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
Dom_Bac_FourTimes_DA_ASV_Tax <- lapply(Dom_Bac_FourTimes_DA_ASV, function(x) CP_unnormalized_bac_ps_Tax[x, ])

save(Dom_Bac_FourTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_FourTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
Dom_Bac_FourTimes_DA_ASV_Tax_Table <- lapply(Dom_Bac_FourTimes_DA_ASV_Tax, function(df) {
  lapply(df[ , 2:6], function(col) table(col, useNA = "ifany"))
})
Dom_Bac_FourTimes_DA_ASV_Tax_Table <- lapply(Dom_Bac_FourTimes_DA_ASV_Tax_Table, function(ranklist) {
                                            lapply(ranklist, as.data.frame)})

# Examples
Dom_Bac_FourTimes_DA_ASV_Tax_Table$HM$Class
```

```
## NULL
```

``` r
Dom_Bac_FourTimes_DA_ASV_Tax_Table$VL$Class
```

```
## NULL
```

``` r
save(Dom_Bac_FourTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_FourTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring five times

``` r
# Filter to ASVs detected by at least 5 tests
Dom_Bac_FiveTimes_DA_ASV <- lapply(Dom_Bac_NumbOcc_DA_ASV, function(x) names(x[x >= 5]))
Dom_Bac_FiveTimes_DA_ASV
```

```
## $Wild
## character(0)
## 
## $Cultivated
## character(0)
```

``` r
# Number of DA ASVs occurring in Five tests
lapply(Dom_Bac_FiveTimes_DA_ASV, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_FiveTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_FiveTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
Dom_Bac_FiveTimes_DA_ASV_Tax <- lapply(Dom_Bac_FiveTimes_DA_ASV, function(x) CP_unnormalized_bac_ps_Tax[x, ])

save(Dom_Bac_FiveTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_FiveTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
Dom_Bac_FiveTimes_DA_ASV_Tax_Table <- lapply(Dom_Bac_FiveTimes_DA_ASV_Tax, function(df) {
  lapply(df[ , 2:6], function(col) table(col, useNA = "ifany"))
})
Dom_Bac_FiveTimes_DA_ASV_Tax_Table <- lapply(Dom_Bac_FiveTimes_DA_ASV_Tax_Table, function(ranklist) {
                                            lapply(ranklist, as.data.frame)})

# Examples
Dom_Bac_FiveTimes_DA_ASV_Tax_Table$HM$Class
```

```
## NULL
```

``` r
Dom_Bac_FiveTimes_DA_ASV_Tax_Table$VL$Class
```

```
## NULL
```

``` r
save(Dom_Bac_FiveTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_FiveTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring six times

``` r
# Filter to ASVs detected by at least 6 tests
Dom_Bac_SixTimes_DA_ASV <- lapply(Dom_Bac_NumbOcc_DA_ASV, function(x) names(x[x >= 6]))
Dom_Bac_SixTimes_DA_ASV
```

```
## $Wild
## character(0)
## 
## $Cultivated
## character(0)
```

``` r
# Number of DA ASVs occurring in Six tests
lapply(Dom_Bac_SixTimes_DA_ASV, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_SixTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_SixTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
Dom_Bac_SixTimes_DA_ASV_Tax <- lapply(Dom_Bac_SixTimes_DA_ASV, function(x) CP_unnormalized_bac_ps_Tax[x, ])

save(Dom_Bac_SixTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_SixTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
Dom_Bac_SixTimes_DA_ASV_Tax_Table <- lapply(Dom_Bac_SixTimes_DA_ASV_Tax, function(df) {
  lapply(df[ , 2:6], function(col) table(col, useNA = "ifany"))
})
Dom_Bac_SixTimes_DA_ASV_Tax_Table <- lapply(Dom_Bac_SixTimes_DA_ASV_Tax_Table, function(ranklist) {
                                            lapply(ranklist, as.data.frame)})

# Print output
Dom_Bac_SixTimes_DA_ASV_Tax_Table
```

```
## $Wild
## $Wild$Phylum
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $Wild$Class
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $Wild$Order
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $Wild$Family
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $Wild$Genus
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## 
## $Cultivated
## $Cultivated$Phylum
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $Cultivated$Class
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $Cultivated$Order
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $Cultivated$Family
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $Cultivated$Genus
## [1] Freq
## <0 rows> (or 0-length row.names)
```

``` r
save(Dom_Bac_SixTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/Dom_Bac_SixTimes_DA_ASV_Tax_Table.RData")
```

``` r
# Clean environment
rm(list = ls())
```


## 4.3.5 Domestication Class level
### Load DA Class

``` r
## Load outputs DESeq
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ClassLevel/sigtab_stddds2_Wald_Dom_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ClassLevel/sigtab_stddds2_LRT_Dom_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ClassLevel/sigtab_zinbWald_Dom_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ClassLevel/sigtab_zinbLRT_Dom_Bac.RData")

## Load output Ancom
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ClassLevel/fdr_ancomWZ_Dom_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ClassLevel/fdr_ancomNoZ_Dom_Bac.RData")

## Load Phyloseq object
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_unnormalized_bac_ps_ForDA.RData")
```

### Extract DA Class

``` r
# From each data frame, extract the Class or ASV IDs caught by each test
Dom_Bac_DESeqLRT_ASV <- lapply(sigtab_stddds2_LRT_Dom_Bac, function(x) rownames(x))
Dom_Bac_DESeqWald_ASV <- lapply(sigtab_stddds2_Wald_Dom_Bac, function(x) rownames(x))
Dom_Bac_zinbLRT_ASV <- lapply(sigtab_zinbLRT_Dom_Bac, function(x) rownames(x))
Dom_Bac_zinbWald_ASV <- lapply(sigtab_zinbWald_Dom_Bac, function(x) rownames(x))
Dom_Bac_ancomWZ_Class <- lapply(fdr_ancomWZ_Dom_Bac, function(x) as.vector(x[,'Species']))
Dom_Bac_ancomNoZ_Class <- lapply(fdr_ancomNoZ_Dom_Bac, function(x) as.vector(x[,'Species']))

# Turn ASV IDS into classes for DESeq results, if no ASVs are present, it will return nothing
Dom_Bac_DESeqLRT_Class <- lapply(Dom_Bac_DESeqLRT_ASV, function(x) {
  if (length(x) == 0) {
    return(character(0))}
  out <- tax_table(CP_unnormalized_bac_ps_ForDA)[x, "Class", drop = FALSE]
  return(out)
})
Dom_Bac_DESeqWald_Class <- lapply(Dom_Bac_DESeqWald_ASV, function(x) {
  if (length(x) == 0) {
    return(character(0))}
  out <- tax_table(CP_unnormalized_bac_ps_ForDA)[x, "Class", drop = FALSE]
  return(out)
})
Dom_Bac_zinbLRT_Class <- lapply(Dom_Bac_zinbLRT_ASV, function(x) {
  if (length(x) == 0) {
    return(character(0))}
  out <- tax_table(CP_unnormalized_bac_ps_ForDA)[x, "Class", drop = FALSE]
  return(out)
})
Dom_Bac_zinbWald_Class <- lapply(Dom_Bac_zinbWald_ASV, function(x) {
  if (length(x) == 0) {
    return(character(0))}
  out <- tax_table(CP_unnormalized_bac_ps_ForDA)[x, "Class", drop = FALSE]
  return(out)
})

# Replace taxa which include kingdom and phylum
Dom_Bac_ancomWZ_Class <- lapply(Dom_Bac_ancomWZ_Class, function(x) sub(".*c__", "c__", x))
Dom_Bac_ancomNoZ_Class <- lapply(Dom_Bac_ancomNoZ_Class, function(x) sub(".*c__", "c__", x))

# Within each plant concatenate the identified asvs
Dom_Bac_Total_DA_Class <- mapply(function(s,t,w,x,y,z){
  input <- c(s,t,w,x,y,z)},
  s = Dom_Bac_DESeqLRT_Class,
  t = Dom_Bac_DESeqWald_Class,
  w = Dom_Bac_zinbLRT_Class,
  x = Dom_Bac_zinbWald_Class,
  y = Dom_Bac_ancomWZ_Class,
  z = Dom_Bac_ancomNoZ_Class,
  SIMPLIFY = FALSE)
```

### Process DA Class data

``` r
# Total number of DA Class
lapply(Dom_Bac_Total_DA_Class, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_Total_DA_Class, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ClassLevel/Dom_Bac_Total_DA_Class.RData")

# Remove duplicated asvs names from that list
Dom_Bac_NoDup_DA_Class <- lapply(Dom_Bac_Total_DA_Class, function(x)
  unique(x))

# Number of DA Class without duplicates
lapply(Dom_Bac_NoDup_DA_Class, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_NoDup_DA_Class, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ClassLevel/Dom_Bac_NoDup_DA_Class.RData")
```

### DA Class occurring twice

``` r
# Filter to Class detected by at least 2 tests
Dom_Bac_NumbOcc_DA_Class <- lapply(Dom_Bac_Total_DA_Class, function(x) table(unlist(x)))
Dom_Bac_TwoTimes_DA_Class <- lapply(Dom_Bac_NumbOcc_DA_Class, function(x) names(x[x >= 2]))
Dom_Bac_TwoTimes_DA_Class
```

```
## $Wild
## NULL
## 
## $Cultivated
## NULL
```

``` r
# Number of DA Class occurring in two tests
lapply(Dom_Bac_TwoTimes_DA_Class, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_TwoTimes_DA_Class, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ClassLevel/Dom_Bac_TwoTimes_DA_Class.RData")
```

### DA Class occurring three times

``` r
# Filter to Class detected by at least 3 tests
Dom_Bac_ThreeTimes_DA_Class <- lapply(Dom_Bac_NumbOcc_DA_Class, function(x) names(x[x >= 3]))
Dom_Bac_ThreeTimes_DA_Class
```

```
## $Wild
## NULL
## 
## $Cultivated
## NULL
```

``` r
# Number of DA Class occurring in three tests
lapply(Dom_Bac_ThreeTimes_DA_Class, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_ThreeTimes_DA_Class, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ClassLevel/Dom_Bac_ThreeTimes_DA_Class.RData")
```

### DA Class occurring four times

``` r
# Filter to Class detected by at least 4 tests
Dom_Bac_FourTimes_DA_Class <- lapply(Dom_Bac_NumbOcc_DA_Class, function(x) names(x[x >= 4]))
Dom_Bac_FourTimes_DA_Class
```

```
## $Wild
## NULL
## 
## $Cultivated
## NULL
```

``` r
# Number of DA Class occurring in four tests
lapply(Dom_Bac_FourTimes_DA_Class, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_FourTimes_DA_Class, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ClassLevel/Dom_Bac_FourTimes_DA_Class.RData")
```

### DA Class occurring five times

``` r
# Filter to Class detected by at least 5 tests
Dom_Bac_FiveTimes_DA_Class <- lapply(Dom_Bac_NumbOcc_DA_Class, function(x) names(x[x >= 5]))
Dom_Bac_FiveTimes_DA_Class
```

```
## $Wild
## NULL
## 
## $Cultivated
## NULL
```

``` r
# Number of DA Class occurring in Five tests
lapply(Dom_Bac_FiveTimes_DA_Class, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_FiveTimes_DA_Class, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ClassLevel/Dom_Bac_FiveTimes_DA_Class.RData")
```

### DA Class occurring six times

``` r
# Filter to Class detected by at least 6 tests
Dom_Bac_SixTimes_DA_Class <- lapply(Dom_Bac_NumbOcc_DA_Class, function(x) names(x[x >= 6]))
Dom_Bac_SixTimes_DA_Class
```

```
## $Wild
## NULL
## 
## $Cultivated
## NULL
```

``` r
# Number of DA Class occurring in Six tests
lapply(Dom_Bac_SixTimes_DA_Class, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_SixTimes_DA_Class, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ClassLevel/Dom_Bac_SixTimes_DA_Class.RData")
```

``` r
# Clean environment
rm(list = ls())
```


## 4.3.6 Domestication Phylum level
### Load DA Phylum

``` r
## Load outputs DESeq
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_PhylumLevel/sigtab_stddds2_Wald_Dom_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_PhylumLevel/sigtab_stddds2_LRT_Dom_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_PhylumLevel/sigtab_zinbWald_Dom_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_PhylumLevel/sigtab_zinbLRT_Dom_Bac.RData")

## Load output Ancom
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_PhylumLevel/fdr_ancomWZ_Dom_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_PhylumLevel/fdr_ancomNoZ_Dom_Bac.RData")

## Load Phyloseq object
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_unnormalized_bac_ps_ForDA.RData")
```

### Extract DA Phylum

``` r
# From each data frame, extract the Phylum or ASV IDs caught by each test
Dom_Bac_DESeqLRT_ASV <- lapply(sigtab_stddds2_LRT_Dom_Bac, function(x) rownames(x))
Dom_Bac_DESeqWald_ASV <- lapply(sigtab_stddds2_Wald_Dom_Bac, function(x) rownames(x))
Dom_Bac_zinbLRT_ASV <- lapply(sigtab_zinbLRT_Dom_Bac, function(x) rownames(x))
Dom_Bac_zinbWald_ASV <- lapply(sigtab_zinbWald_Dom_Bac, function(x) rownames(x))
Dom_Bac_ancomWZ_Phylum <- lapply(fdr_ancomWZ_Dom_Bac, function(x) as.vector(x[,'Species']))
Dom_Bac_ancomNoZ_Phylum <- lapply(fdr_ancomNoZ_Dom_Bac, function(x) as.vector(x[,'Species']))

# Turn ASV IDS into classes for DESeq results, if no ASVs are present, it will return nothing
Dom_Bac_DESeqLRT_Phylum <- lapply(Dom_Bac_DESeqLRT_ASV, function(x) {
  if (length(x) == 0) {
    return(character(0))}
  out <- tax_table(CP_unnormalized_bac_ps_ForDA)[x, "Phylum", drop = FALSE]
  return(out)
})
Dom_Bac_DESeqWald_Phylum <- lapply(Dom_Bac_DESeqWald_ASV, function(x) {
  if (length(x) == 0) {
    return(character(0))}
  out <- tax_table(CP_unnormalized_bac_ps_ForDA)[x, "Phylum", drop = FALSE]
  return(out)
})
Dom_Bac_zinbLRT_Phylum <- lapply(Dom_Bac_zinbLRT_ASV, function(x) {
  if (length(x) == 0) {
    return(character(0))}
  out <- tax_table(CP_unnormalized_bac_ps_ForDA)[x, "Phylum", drop = FALSE]
  return(out)
})
Dom_Bac_zinbWald_Phylum <- lapply(Dom_Bac_zinbWald_ASV, function(x) {
  if (length(x) == 0) {
    return(character(0))}
  out <- tax_table(CP_unnormalized_bac_ps_ForDA)[x, "Phylum", drop = FALSE]
  return(out)
})

# Within each plant concatenate the identified asvs
Dom_Bac_Total_DA_Phylum <- mapply(function(s,t,w,x,y,z){
  input <- c(s,t,w,x,y,z)},
  s = Dom_Bac_DESeqLRT_Phylum,
  t = Dom_Bac_DESeqWald_Phylum,
  w = Dom_Bac_zinbLRT_Phylum,
  x = Dom_Bac_zinbWald_Phylum,
  y = Dom_Bac_ancomWZ_Phylum,
  z = Dom_Bac_ancomNoZ_Phylum,
  SIMPLIFY = FALSE)
```

### Process DA Phylum data

``` r
# Total number of DA Phylum
lapply(Dom_Bac_Total_DA_Phylum, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_Total_DA_Phylum, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_PhylumLevel/Dom_Bac_Total_DA_Phylum.RData")

# Remove duplicated asvs names from that list
Dom_Bac_NoDup_DA_Phylum <- lapply(Dom_Bac_Total_DA_Phylum, function(x)
  unique(x))

# Number of DA Phylum without duplicates
lapply(Dom_Bac_NoDup_DA_Phylum, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_NoDup_DA_Phylum, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_PhylumLevel/Dom_Bac_NoDup_DA_Phylum.RData")
```

### DA Phylum occurring twice

``` r
# Filter to Phylum detected by at least 2 tests
Dom_Bac_NumbOcc_DA_Phylum <- lapply(Dom_Bac_Total_DA_Phylum, function(x) table(unlist(x)))
Dom_Bac_TwoTimes_DA_Phylum <- lapply(Dom_Bac_NumbOcc_DA_Phylum, function(x) names(x[x >= 2]))
Dom_Bac_TwoTimes_DA_Phylum
```

```
## $Wild
## NULL
## 
## $Cultivated
## NULL
```

``` r
# Number of DA Phylum occurring in two tests
lapply(Dom_Bac_TwoTimes_DA_Phylum, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_TwoTimes_DA_Phylum, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_PhylumLevel/Dom_Bac_TwoTimes_DA_Phylum.RData")
```

### DA Phylum occurring three times

``` r
# Filter to Phylum detected by at least 3 tests
Dom_Bac_ThreeTimes_DA_Phylum <- lapply(Dom_Bac_NumbOcc_DA_Phylum, function(x) names(x[x >= 3]))
Dom_Bac_ThreeTimes_DA_Phylum
```

```
## $Wild
## NULL
## 
## $Cultivated
## NULL
```

``` r
# Number of DA Phylum occurring in three tests
lapply(Dom_Bac_ThreeTimes_DA_Phylum, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_ThreeTimes_DA_Phylum, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_PhylumLevel/Dom_Bac_ThreeTimes_DA_Phylum.RData")
```

### DA Phylum occurring four times

``` r
# Filter to Phylum detected by at least 4 tests
Dom_Bac_FourTimes_DA_Phylum <- lapply(Dom_Bac_NumbOcc_DA_Phylum, function(x) names(x[x >= 4]))
Dom_Bac_FourTimes_DA_Phylum
```

```
## $Wild
## NULL
## 
## $Cultivated
## NULL
```

``` r
# Number of DA Phylum occurring in four tests
lapply(Dom_Bac_FourTimes_DA_Phylum, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_FourTimes_DA_Phylum, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_PhylumLevel/Dom_Bac_FourTimes_DA_Phylum.RData")
```

### DA Phylum occurring five times

``` r
# Filter to Phylum detected by at least 5 tests
Dom_Bac_FiveTimes_DA_Phylum <- lapply(Dom_Bac_NumbOcc_DA_Phylum, function(x) names(x[x >= 5]))
Dom_Bac_FiveTimes_DA_Phylum
```

```
## $Wild
## NULL
## 
## $Cultivated
## NULL
```

``` r
# Number of DA Phylum occurring in Five tests
lapply(Dom_Bac_FiveTimes_DA_Phylum, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_FiveTimes_DA_Phylum, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_PhylumLevel/Dom_Bac_FiveTimes_DA_Phylum.RData")
```

### DA Phylum occurring six times

``` r
# Filter to Phylum detected by at least 6 tests
Dom_Bac_SixTimes_DA_Phylum <- lapply(Dom_Bac_NumbOcc_DA_Phylum, function(x) names(x[x >= 6]))
Dom_Bac_SixTimes_DA_Phylum
```

```
## $Wild
## NULL
## 
## $Cultivated
## NULL
```

``` r
# Number of DA Phylum occurring in Six tests
lapply(Dom_Bac_SixTimes_DA_Phylum, function(x) length(x))
```

```
## $Wild
## [1] 0
## 
## $Cultivated
## [1] 0
```

``` r
save(Dom_Bac_SixTimes_DA_Phylum, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_PhylumLevel/Dom_Bac_SixTimes_DA_Phylum.RData")
```

``` r
# Clean environment
rm(list = ls())
```


## 4.3.7 All data ASV level
### Load DA ASVs

``` r
## Load outputs DESeq
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ASVLevel/sigtab_stddds2_Wald_All_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ASVLevel/sigtab_stddds2_LRT_All_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ASVLevel/sigtab_zinbWald_All_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_DESeq_ASVLevel/sigtab_zinbLRT_All_Bac.RData")

## Load output Ancom
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/fdr_ancomWZ_All_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_Ancom_ASVLevel/fdr_ancomNoZ_All_Bac.RData")
```

### Extract DA ASVs

``` r
## From each data frame, extract the ASV IDs caught by each test
All_Bac_DESeqLRT_ASVs <- rownames(sigtab_stddds2_LRT_All_Bac)
All_Bac_DESeqWald_ASVs <- rownames(sigtab_stddds2_Wald_All_Bac)
All_Bac_zinbLRT_ASVs <- rownames(sigtab_zinbLRT_All_Bac)
All_Bac_zinbWald_ASVs <- rownames(sigtab_zinbWald_All_Bac)
All_Bac_ancomWZ_ASVs <- as.vector(fdr_ancomWZ_All_Bac[ ,'Species'])
All_Bac_ancomNoZ_ASVs <- as.vector(fdr_ancomNoZ_All_Bac[ ,'Species'])

All_Bac_Total_DA_ASV <- c(All_Bac_DESeqLRT_ASVs, All_Bac_DESeqWald_ASVs, All_Bac_zinbLRT_ASVs, All_Bac_zinbWald_ASVs, All_Bac_ancomWZ_ASVs, All_Bac_ancomNoZ_ASVs)
```

### Total number of ASVs

``` r
# Total number of DA ASVs
length(All_Bac_Total_DA_ASV)
```

```
## [1] 79
```

``` r
save(All_Bac_Total_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_Total_DA_ASV.RData")

# Remove duplicated asvs names from that list
All_Bac_NoDup_DA_ASV <- unique(All_Bac_Total_DA_ASV)

# Number of DA ASVs without duplicates
length(All_Bac_NoDup_DA_ASV)
```

```
## [1] 39
```

``` r
save(All_Bac_NoDup_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_NoDup_DA_ASV.RData")

# Extracting taxonomy
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_unnormalized_bac_ps.RData")

# Extract taxonomy of DA ASVs
CP_unnormalized_bac_ps_Tax <- as.data.frame(tax_table(CP_unnormalized_bac_ps))
All_Bac_NoDup_DA_ASV_Tax <- CP_unnormalized_bac_ps_Tax[All_Bac_NoDup_DA_ASV, ]

save(All_Bac_NoDup_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects//CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_NoDup_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
All_Bac_NoDup_DA_ASV_Tax_Table <- lapply(All_Bac_NoDup_DA_ASV_Tax[ , 2:6], function(df) table(df))
All_Bac_NoDup_DA_ASV_Tax_Table <- lapply(All_Bac_NoDup_DA_ASV_Tax_Table, as.data.frame)

# Examples
All_Bac_NoDup_DA_ASV_Tax_Table$Class
```

```
##                        df Freq
## 1       c__Actinobacteria    1
## 2  c__Alphaproteobacteria    2
## 3              c__Bacilli    8
## 4       c__Blastocatellia    1
## 5  c__Gammaproteobacteria    1
## 6       c__Gemmatimonadia    5
## 7               c__KD4-96    3
## 8  c__Thermoanaerobaculia    1
## 9      c__Thermoleophilia   15
## 10    c__Vicinamibacteria    2
```

``` r
save(All_Bac_NoDup_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_NoDup_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring twice

``` r
# Filter to ASVs detected by at least 2 tests
All_Bac_NumbOcc_DA_ASV <- table(unlist(All_Bac_Total_DA_ASV))
All_Bac_TwoTimes_DA_ASV <- names(All_Bac_NumbOcc_DA_ASV[All_Bac_NumbOcc_DA_ASV >= 2])
All_Bac_TwoTimes_DA_ASV
```

```
##  [1] "bASV_1025" "bASV_1083" "bASV_1224" "bASV_1272" "bASV_1298" "bASV_1389"
##  [7] "bASV_1490" "bASV_1701" "bASV_1760" "bASV_1850" "bASV_2073" "bASV_2402"
## [13] "bASV_2416" "bASV_2417" "bASV_2629" "bASV_2677" "bASV_269"  "bASV_2712"
## [19] "bASV_287"  "bASV_300"  "bASV_3025" "bASV_3143" "bASV_359"  "bASV_365" 
## [25] "bASV_387"  "bASV_4949" "bASV_586"  "bASV_601"  "bASV_656"  "bASV_683" 
## [31] "bASV_724"  "bASV_751"  "bASV_878"  "bASV_911"  "bASV_938"  "bASV_963" 
## [37] "bASV_979"
```

``` r
# Number of DA ASVs occurring in two tests
length(All_Bac_TwoTimes_DA_ASV)
```

```
## [1] 37
```

``` r
save(All_Bac_TwoTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_TwoTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
All_Bac_TwoTimes_DA_ASV_Tax <- CP_unnormalized_bac_ps_Tax[All_Bac_TwoTimes_DA_ASV, ]

save(All_Bac_TwoTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects//CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_TwoTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
All_Bac_TwoTimes_DA_ASV_Tax_Table <- lapply(All_Bac_TwoTimes_DA_ASV_Tax[ , 2:6], function(df) table(df))
All_Bac_TwoTimes_DA_ASV_Tax_Table <- lapply(All_Bac_TwoTimes_DA_ASV_Tax_Table, as.data.frame)

# Examples
All_Bac_TwoTimes_DA_ASV_Tax_Table$Class
```

```
##                       df Freq
## 1      c__Actinobacteria    1
## 2 c__Alphaproteobacteria    2
## 3             c__Bacilli    8
## 4      c__Blastocatellia    1
## 5      c__Gemmatimonadia    5
## 6              c__KD4-96    3
## 7     c__Thermoleophilia   15
## 8    c__Vicinamibacteria    2
```

``` r
save(All_Bac_TwoTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_TwoTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring three times

``` r
# Filter to ASVs detected by at least 3 tests
All_Bac_ThreeTimes_DA_ASV <- names(All_Bac_NumbOcc_DA_ASV[All_Bac_NumbOcc_DA_ASV >= 3])
All_Bac_ThreeTimes_DA_ASV
```

```
## [1] "bASV_269" "bASV_601"
```

``` r
# Number of DA ASVs occurring in three tests
length(All_Bac_ThreeTimes_DA_ASV)
```

```
## [1] 2
```

``` r
save(All_Bac_ThreeTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_ThreeTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
All_Bac_ThreeTimes_DA_ASV_Tax <- CP_unnormalized_bac_ps_Tax[All_Bac_ThreeTimes_DA_ASV, ]

save(All_Bac_ThreeTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects//CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_ThreeTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
All_Bac_ThreeTimes_DA_ASV_Tax_Table <- lapply(All_Bac_ThreeTimes_DA_ASV_Tax[ , 2:6], function(df) table(df))
All_Bac_ThreeTimes_DA_ASV_Tax_Table <- lapply(All_Bac_ThreeTimes_DA_ASV_Tax_Table, as.data.frame)

# Examples
All_Bac_ThreeTimes_DA_ASV_Tax_Table
```

```
## $Phylum
##                  df Freq
## 1 p__Actinomycetota    2
## 
## $Class
##                   df Freq
## 1 c__Thermoleophilia    2
## 
## $Order
##              df Freq
## 1 o__Gaiellales    2
## 
## $Family
##                  df Freq
## 1 f__Incertae_Sedis    2
## 
## $Genus
##                  df Freq
## 1 g__Incertae_Sedis    2
```

``` r
save(All_Bac_ThreeTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_ThreeTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring four times

``` r
# Filter to ASVs detected by at least 3 tests
All_Bac_FourTimes_DA_ASV <- names(All_Bac_NumbOcc_DA_ASV[All_Bac_NumbOcc_DA_ASV >= 4])
All_Bac_FourTimes_DA_ASV
```

```
## [1] "bASV_601"
```

``` r
# Number of DA ASVs occurring in three tests
length(All_Bac_FourTimes_DA_ASV)
```

```
## [1] 1
```

``` r
save(All_Bac_FourTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_FourTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
All_Bac_FourTimes_DA_ASV_Tax <- CP_unnormalized_bac_ps_Tax[All_Bac_FourTimes_DA_ASV, ]

save(All_Bac_FourTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects//CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_FourTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
All_Bac_FourTimes_DA_ASV_Tax_Table <- lapply(All_Bac_FourTimes_DA_ASV_Tax[ , 2:6], function(df) table(df))
All_Bac_FourTimes_DA_ASV_Tax_Table <- lapply(All_Bac_FourTimes_DA_ASV_Tax_Table, as.data.frame)

# Examples
All_Bac_FourTimes_DA_ASV_Tax_Table
```

```
## $Phylum
##                  df Freq
## 1 p__Actinomycetota    1
## 
## $Class
##                   df Freq
## 1 c__Thermoleophilia    1
## 
## $Order
##              df Freq
## 1 o__Gaiellales    1
## 
## $Family
##                  df Freq
## 1 f__Incertae_Sedis    1
## 
## $Genus
##                  df Freq
## 1 g__Incertae_Sedis    1
```

``` r
save(All_Bac_FourTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_FourTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring five times

``` r
# Filter to ASVs detected by at least 3 tests
All_Bac_FiveTimes_DA_ASV <- names(All_Bac_NumbOcc_DA_ASV[All_Bac_NumbOcc_DA_ASV >= 5])
All_Bac_FiveTimes_DA_ASV
```

```
## character(0)
```

``` r
# Number of DA ASVs occurring in three tests
length(All_Bac_FiveTimes_DA_ASV)
```

```
## [1] 0
```

``` r
save(All_Bac_FiveTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_FiveTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
All_Bac_FiveTimes_DA_ASV_Tax <- CP_unnormalized_bac_ps_Tax[All_Bac_FiveTimes_DA_ASV, ]

save(All_Bac_FiveTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects//CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_FiveTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
All_Bac_FiveTimes_DA_ASV_Tax_Table <- lapply(All_Bac_FiveTimes_DA_ASV_Tax[ , 2:6], function(df) table(df))
All_Bac_FiveTimes_DA_ASV_Tax_Table <- lapply(All_Bac_FiveTimes_DA_ASV_Tax_Table, as.data.frame)

# Examples
All_Bac_FiveTimes_DA_ASV_Tax_Table
```

```
## $Phylum
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $Class
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $Order
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $Family
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $Genus
## [1] Freq
## <0 rows> (or 0-length row.names)
```

``` r
save(All_Bac_FiveTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_FiveTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring six times

``` r
# Filter to ASVs detected by at least 3 tests
All_Bac_SixTimes_DA_ASV <- names(All_Bac_NumbOcc_DA_ASV[All_Bac_NumbOcc_DA_ASV >= 6])
All_Bac_SixTimes_DA_ASV
```

```
## character(0)
```

``` r
# Number of DA ASVs occurring in three tests
length(All_Bac_SixTimes_DA_ASV)
```

```
## [1] 0
```

``` r
save(All_Bac_SixTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_SixTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
All_Bac_SixTimes_DA_ASV_Tax <- CP_unnormalized_bac_ps_Tax[All_Bac_SixTimes_DA_ASV, ]

save(All_Bac_SixTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects//CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_SixTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
All_Bac_SixTimes_DA_ASV_Tax_Table <- lapply(All_Bac_SixTimes_DA_ASV_Tax[ , 2:6], function(df) table(df))
All_Bac_SixTimes_DA_ASV_Tax_Table <- lapply(All_Bac_SixTimes_DA_ASV_Tax_Table, as.data.frame)

# Examples
All_Bac_SixTimes_DA_ASV_Tax_Table
```

```
## $Phylum
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $Class
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $Order
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $Family
## [1] Freq
## <0 rows> (or 0-length row.names)
## 
## $Genus
## [1] Freq
## <0 rows> (or 0-length row.names)
```

``` r
save(All_Bac_SixTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/CP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_SixTimes_DA_ASV_Tax_Table.RData")
```

``` r
# Clean environment
rm(list = ls())
```


## 4.3.8 All data Class and Phylum level
No DA classes and phyla detected when analysing whole data.
