---
title: "FP_05_differential_abundance_Summary_Results - Data preperation - Only batch 2"
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


# 5.3 Create lists of DA ASVs
## 5.3.1 Accession ASV level
### Load DA ASVs

``` r
## Load outputs DESeq
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_ASVLevel/sigtab_stddds2_Wald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_ASVLevel/sigtab_stddds2_LRT_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_ASVLevel/sigtab_zinbWald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_ASVLevel/sigtab_zinbLRT_Acc_Bac.RData")

## Load output Ancom
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ASVLevel/fdr_ancomWZ_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ASVLevel/fdr_ancomNoZ_Acc_Bac.RData")
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
## $VL
## [1] 822
## 
## $CD
## [1] 761
## 
## $RI
## [1] 925
## 
## $HM
## [1] 1341
## 
## $GO1
## [1] 717
```

``` r
save(Acc_Bac_Total_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_Total_DA_ASV.RData")

# Remove duplicated asvs names from that list
Acc_Bac_NoDup_DA_ASV <- lapply(Acc_Bac_Total_DA_ASV, function(x)
  unique(x))

# Number of DA ASVs without duplicates
lapply(Acc_Bac_NoDup_DA_ASV, function(x) length(x))
```

```
## $VL
## [1] 667
## 
## $CD
## [1] 700
## 
## $RI
## [1] 858
## 
## $HM
## [1] 919
## 
## $GO1
## [1] 661
```

``` r
save(Acc_Bac_NoDup_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_NoDup_DA_ASV.RData")

# Extracting taxonomy
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_unnormalized_bac_ps.RData")

# Extract taxonomy of DA ASVs
FP_unnormalized_bac_ps_Tax <- as.data.frame(tax_table(FP_unnormalized_bac_ps))
Acc_Bac_NoDup_DA_ASV_Tax <- lapply(Acc_Bac_NoDup_DA_ASV, function(x) FP_unnormalized_bac_ps_Tax[x, ])

save(Acc_Bac_NoDup_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects//FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_NoDup_DA_ASV_Tax.RData")

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
##                           col Freq
## 1           c__Acidimicrobiia   20
## 2           c__Acidobacteriae   19
## 3           c__Actinobacteria  112
## 4      c__Alphaproteobacteria  155
## 5                 c__Babeliae    1
## 6                  c__Bacilli   54
## 7              c__Bacteroidia   28
## 8          c__Bdellovibrionia    3
## 9           c__Blastocatellia    5
## 10            c__Chloroflexia   12
## 11              c__Clostridia    4
## 12              c__Deinococci    1
## 13        c__Desulfovibrionia    1
## 14        c__Desulfuromonadia    2
## 15           c__Fibrobacteria    1
## 16     c__Gammaproteobacteria  231
## 17          c__Gemmatimonadia   78
## 18             c__Gitt-GS-136    2
## 19         c__Gracilibacteria    1
## 20          c__Incertae_Sedis    3
## 21            c__JG30-KF-CM66    2
## 22                  c__KD4-96   21
## 23         c__Ktedonobacteria   27
## 24            c__Limnochordia    1
## 25               c__MB-A2-108    3
## 26              c__Myxococcia    1
## 27             c__Nitrospiria    2
## 28                   c__OLB14    1
## 29             c__Oligoflexia    8
## 30           c__Parcubacteria    2
## 31          c__Planctomycetes    1
## 32              c__Polyangiia    9
## 33 c__S0134_terrestrial_group    1
## 34         c__Saccharimonadia    2
## 35         c__Symbiobacteriia    1
## 36     c__Thermoanaerobaculia    1
## 37         c__Thermoleophilia   58
## 38                    c__TK10    4
## 39        c__Vampirivibrionia    1
## 40        c__Verrucomicrobiia   14
## 41        c__Vicinamibacteria   21
## 42                       <NA>    5
```

``` r
Acc_Bac_NoDup_DA_ASV_Tax_Table$VL$Class
```

```
##                           col Freq
## 1           c__Acidimicrobiia   13
## 2           c__Acidobacteriae   14
## 3           c__Actinobacteria   87
## 4      c__Alphaproteobacteria  107
## 5                  c__Bacilli   57
## 6              c__Bacteroidia   15
## 7          c__Bdellovibrionia    1
## 8           c__Blastocatellia    4
## 9             c__Chloroflexia   10
## 10              c__Clostridia    4
## 11           c__Fibrobacteria    1
## 12     c__Gammaproteobacteria  172
## 13          c__Gemmatimonadia   53
## 14             c__Gitt-GS-136    2
## 15           c__Halanaerobiia    1
## 16              c__Holophagae    3
## 17          c__Incertae_Sedis    2
## 18            c__JG30-KF-CM66    1
## 19                  c__KD4-96    7
## 20         c__Ktedonobacteria   10
## 21          c__Longimicrobiia    1
## 22               c__MB-A2-108    1
## 23             c__Nitrospiria    1
## 24                   c__OLB14    1
## 25             c__Oligoflexia    1
## 26          c__Planctomycetes    1
## 27              c__Polyangiia    3
## 28 c__S0134_terrestrial_group    1
## 29         c__Saccharimonadia    1
## 30         c__Symbiobacteriia    1
## 31         c__Thermoleophilia   60
## 32                    c__TK10    2
## 33        c__Verrucomicrobiia    3
## 34        c__Vicinamibacteria   21
## 35                       <NA>    5
```

``` r
save(Acc_Bac_NoDup_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_NoDup_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring twice

``` r
# Filter to ASVs detected by at least 2 tests
Acc_Bac_NumbOcc_DA_ASV <- lapply(Acc_Bac_Total_DA_ASV, function(x) table(unlist(x)))
Acc_Bac_TwoTimes_DA_ASV <- lapply(Acc_Bac_NumbOcc_DA_ASV, function(x) names(x[x >= 2]))
Acc_Bac_TwoTimes_DA_ASV
```

```
## $VL
##   [1] "bASV_10"    "bASV_1003"  "bASV_10085" "bASV_1025"  "bASV_1028" 
##   [6] "bASV_103"   "bASV_1067"  "bASV_1080"  "bASV_109"   "bASV_11178"
##  [11] "bASV_11225" "bASV_1130"  "bASV_1153"  "bASV_1155"  "bASV_1157" 
##  [16] "bASV_12161" "bASV_1224"  "bASV_12339" "bASV_124"   "bASV_126"  
##  [21] "bASV_1285"  "bASV_1329"  "bASV_1331"  "bASV_1345"  "bASV_1352" 
##  [26] "bASV_1401"  "bASV_14493" "bASV_1456"  "bASV_1539"  "bASV_154"  
##  [31] "bASV_16090" "bASV_166"   "bASV_17"    "bASV_1728"  "bASV_178"  
##  [36] "bASV_187"   "bASV_1962"  "bASV_199"   "bASV_2012"  "bASV_2027" 
##  [41] "bASV_2120"  "bASV_214"   "bASV_22"    "bASV_2296"  "bASV_2391" 
##  [46] "bASV_2399"  "bASV_2459"  "bASV_2473"  "bASV_2487"  "bASV_24934"
##  [51] "bASV_2528"  "bASV_256"   "bASV_2578"  "bASV_27"    "bASV_2701" 
##  [56] "bASV_2709"  "bASV_2724"  "bASV_2770"  "bASV_290"   "bASV_2902" 
##  [61] "bASV_303"   "bASV_315"   "bASV_319"   "bASV_3204"  "bASV_3207" 
##  [66] "bASV_3208"  "bASV_3276"  "bASV_329"   "bASV_3526"  "bASV_3753" 
##  [71] "bASV_3806"  "bASV_385"   "bASV_386"   "bASV_3873"  "bASV_39"   
##  [76] "bASV_4029"  "bASV_417"   "bASV_4207"  "bASV_430"   "bASV_4392" 
##  [81] "bASV_47"    "bASV_474"   "bASV_4782"  "bASV_515"   "bASV_534"  
##  [86] "bASV_541"   "bASV_5565"  "bASV_559"   "bASV_560"   "bASV_5682" 
##  [91] "bASV_57"    "bASV_5805"  "bASV_5974"  "bASV_60"    "bASV_63"   
##  [96] "bASV_656"   "bASV_664"   "bASV_6708"  "bASV_7"     "bASV_732"  
## [101] "bASV_827"   "bASV_8335"  "bASV_85"    "bASV_8506"  "bASV_888"  
## [106] "bASV_8895"  "bASV_8917"  "bASV_899"   "bASV_91"    "bASV_910"  
## [111] "bASV_939"   "bASV_9579"  "bASV_968"   "bASV_969"  
## 
## $CD
##  [1] "bASV_1003" "bASV_1055" "bASV_1321" "bASV_1371" "bASV_141"  "bASV_1519"
##  [7] "bASV_1591" "bASV_1652" "bASV_1994" "bASV_2073" "bASV_214"  "bASV_2146"
## [13] "bASV_2290" "bASV_2445" "bASV_2747" "bASV_2777" "bASV_288"  "bASV_302" 
## [19] "bASV_3075" "bASV_3190" "bASV_340"  "bASV_345"  "bASV_426"  "bASV_444" 
## [25] "bASV_446"  "bASV_496"  "bASV_5157" "bASV_528"  "bASV_5707" "bASV_578" 
## [31] "bASV_594"  "bASV_61"   "bASV_650"  "bASV_660"  "bASV_729"  "bASV_900" 
## [37] "bASV_917" 
## 
## $RI
##  [1] "bASV_1061"  "bASV_10663" "bASV_1164"  "bASV_136"   "bASV_1400" 
##  [6] "bASV_143"   "bASV_1463"  "bASV_1525"  "bASV_1562"  "bASV_1732" 
## [11] "bASV_1734"  "bASV_174"   "bASV_1800"  "bASV_1839"  "bASV_185"  
## [16] "bASV_1853"  "bASV_1858"  "bASV_190"   "bASV_1956"  "bASV_1965" 
## [21] "bASV_2016"  "bASV_2195"  "bASV_2316"  "bASV_2322"  "bASV_246"  
## [26] "bASV_2608"  "bASV_2960"  "bASV_3012"  "bASV_3147"  "bASV_319"  
## [31] "bASV_3347"  "bASV_3691"  "bASV_3749"  "bASV_3784"  "bASV_3980" 
## [36] "bASV_4078"  "bASV_430"   "bASV_438"   "bASV_450"   "bASV_4811" 
## [41] "bASV_5341"  "bASV_5489"  "bASV_557"   "bASV_572"   "bASV_6336" 
## [46] "bASV_634"   "bASV_6384"  "bASV_6577"  "bASV_6708"  "bASV_6736" 
## [51] "bASV_6872"  "bASV_7325"  "bASV_7603"  "bASV_7803"  "bASV_848"  
## [56] "bASV_8504"  "bASV_90"    "bASV_9406"  "bASV_980"  
## 
## $HM
##   [1] "bASV_1"     "bASV_1002"  "bASV_1003"  "bASV_1034"  "bASV_1062" 
##   [6] "bASV_107"   "bASV_1071"  "bASV_1077"  "bASV_10830" "bASV_109"  
##  [11] "bASV_1092"  "bASV_1107"  "bASV_1108"  "bASV_1114"  "bASV_11148"
##  [16] "bASV_1139"  "bASV_1153"  "bASV_1163"  "bASV_1164"  "bASV_11668"
##  [21] "bASV_1177"  "bASV_1202"  "bASV_1231"  "bASV_1259"  "bASV_126"  
##  [26] "bASV_12754" "bASV_1289"  "bASV_129"   "bASV_1293"  "bASV_13076"
##  [31] "bASV_13121" "bASV_1315"  "bASV_1325"  "bASV_1345"  "bASV_1354" 
##  [36] "bASV_1366"  "bASV_1371"  "bASV_1385"  "bASV_1405"  "bASV_141"  
##  [41] "bASV_1428"  "bASV_1446"  "bASV_1463"  "bASV_148"   "bASV_149"  
##  [46] "bASV_1499"  "bASV_15018" "bASV_1515"  "bASV_1520"  "bASV_1522" 
##  [51] "bASV_1538"  "bASV_1539"  "bASV_1544"  "bASV_1564"  "bASV_1575" 
##  [56] "bASV_1584"  "bASV_1593"  "bASV_1602"  "bASV_1603"  "bASV_1618" 
##  [61] "bASV_1665"  "bASV_1730"  "bASV_1742"  "bASV_176"   "bASV_1762" 
##  [66] "bASV_178"   "bASV_1780"  "bASV_179"   "bASV_1801"  "bASV_182"  
##  [71] "bASV_183"   "bASV_1853"  "bASV_187"   "bASV_1899"  "bASV_190"  
##  [76] "bASV_1908"  "bASV_19120" "bASV_1933"  "bASV_1934"  "bASV_1940" 
##  [81] "bASV_1962"  "bASV_1963"  "bASV_1969"  "bASV_1988"  "bASV_202"  
##  [86] "bASV_203"   "bASV_2062"  "bASV_2105"  "bASV_2132"  "bASV_2145" 
##  [91] "bASV_2163"  "bASV_2188"  "bASV_2195"  "bASV_2218"  "bASV_2258" 
##  [96] "bASV_227"   "bASV_2277"  "bASV_228"   "bASV_2290"  "bASV_22948"
## [101] "bASV_2299"  "bASV_234"   "bASV_2352"  "bASV_2358"  "bASV_2368" 
## [106] "bASV_2378"  "bASV_2387"  "bASV_2393"  "bASV_2425"  "bASV_243"  
## [111] "bASV_2502"  "bASV_254"   "bASV_2540"  "bASV_2565"  "bASV_2577" 
## [116] "bASV_2650"  "bASV_2655"  "bASV_268"   "bASV_27"    "bASV_2717" 
## [121] "bASV_2727"  "bASV_2737"  "bASV_2767"  "bASV_2770"  "bASV_279"  
## [126] "bASV_2861"  "bASV_2867"  "bASV_288"   "bASV_292"   "bASV_2949" 
## [131] "bASV_2955"  "bASV_2973"  "bASV_3"     "bASV_3012"  "bASV_3046" 
## [136] "bASV_309"   "bASV_3138"  "bASV_3158"  "bASV_326"   "bASV_330"  
## [141] "bASV_3348"  "bASV_335"   "bASV_339"   "bASV_3393"  "bASV_3403" 
## [146] "bASV_3408"  "bASV_345"   "bASV_3469"  "bASV_3504"  "bASV_3526" 
## [151] "bASV_3533"  "bASV_3537"  "bASV_3546"  "bASV_3547"  "bASV_356"  
## [156] "bASV_363"   "bASV_3660"  "bASV_3726"  "bASV_376"   "bASV_3768" 
## [161] "bASV_378"   "bASV_382"   "bASV_3823"  "bASV_385"   "bASV_3862" 
## [166] "bASV_3920"  "bASV_394"   "bASV_3962"  "bASV_397"   "bASV_4036" 
## [171] "bASV_4079"  "bASV_410"   "bASV_411"   "bASV_4116"  "bASV_412"  
## [176] "bASV_4149"  "bASV_416"   "bASV_42"    "bASV_423"   "bASV_430"  
## [181] "bASV_4314"  "bASV_4395"  "bASV_4415"  "bASV_45"    "bASV_4528" 
## [186] "bASV_4569"  "bASV_464"   "bASV_4655"  "bASV_4664"  "bASV_474"  
## [191] "bASV_4762"  "bASV_477"   "bASV_4802"  "bASV_483"   "bASV_4871" 
## [196] "bASV_4872"  "bASV_491"   "bASV_4912"  "bASV_496"   "bASV_499"  
## [201] "bASV_5"     "bASV_5015"  "bASV_5040"  "bASV_510"   "bASV_511"  
## [206] "bASV_515"   "bASV_5307"  "bASV_5339"  "bASV_541"   "bASV_542"  
## [211] "bASV_5456"  "bASV_5497"  "bASV_5529"  "bASV_557"   "bASV_561"  
## [216] "bASV_564"   "bASV_566"   "bASV_5682"  "bASV_5762"  "bASV_5791" 
## [221] "bASV_5795"  "bASV_5864"  "bASV_5893"  "bASV_617"   "bASV_62"   
## [226] "bASV_6343"  "bASV_6380"  "bASV_65"    "bASV_650"   "bASV_6504" 
## [231] "bASV_664"   "bASV_670"   "bASV_673"   "bASV_69"    "bASV_697"  
## [236] "bASV_7"     "bASV_713"   "bASV_7274"  "bASV_73"    "bASV_7439" 
## [241] "bASV_753"   "bASV_812"   "bASV_8134"  "bASV_8136"  "bASV_8148" 
## [246] "bASV_8223"  "bASV_824"   "bASV_831"   "bASV_840"   "bASV_8409" 
## [251] "bASV_855"   "bASV_88"    "bASV_8819"  "bASV_910"   "bASV_916"  
## [256] "bASV_928"   "bASV_930"   "bASV_9493"  "bASV_9605"  "bASV_961"  
## [261] "bASV_974"   "bASV_9867"  "bASV_9884"  "bASV_991"   "bASV_995"  
## 
## $GO1
##  [1] "bASV_1003"  "bASV_1060"  "bASV_1146"  "bASV_11729" "bASV_1196" 
##  [6] "bASV_1203"  "bASV_1355"  "bASV_14"    "bASV_15"    "bASV_160"  
## [11] "bASV_1621"  "bASV_17"    "bASV_176"   "bASV_177"   "bASV_1798" 
## [16] "bASV_192"   "bASV_2152"  "bASV_2377"  "bASV_25"    "bASV_2540" 
## [21] "bASV_26"    "bASV_2612"  "bASV_27"    "bASV_271"   "bASV_2979" 
## [26] "bASV_3249"  "bASV_340"   "bASV_3744"  "bASV_385"   "bASV_386"  
## [31] "bASV_41"    "bASV_426"   "bASV_46"    "bASV_477"   "bASV_5196" 
## [36] "bASV_549"   "bASV_60"    "bASV_6205"  "bASV_65"    "bASV_650"  
## [41] "bASV_66"    "bASV_79"    "bASV_794"   "bASV_9540"  "bASV_9867"
```

``` r
# Number of DA ASVs occurring in two tests
lapply(Acc_Bac_TwoTimes_DA_ASV, function(x) length(x))
```

```
## $VL
## [1] 114
## 
## $CD
## [1] 37
## 
## $RI
## [1] 59
## 
## $HM
## [1] 265
## 
## $GO1
## [1] 45
```

``` r
save(Acc_Bac_TwoTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_TwoTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
Acc_Bac_TwoTimes_DA_ASV_Tax <- lapply(Acc_Bac_TwoTimes_DA_ASV, function(x) FP_unnormalized_bac_ps_Tax[x, ])

save(Acc_Bac_TwoTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_TwoTimes_DA_ASV_Tax.RData")

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
##                           col Freq
## 1           c__Acidimicrobiia   11
## 2           c__Acidobacteriae    9
## 3           c__Actinobacteria   36
## 4      c__Alphaproteobacteria   60
## 5                  c__Bacilli    3
## 6              c__Bacteroidia   11
## 7          c__Bdellovibrionia    1
## 8           c__Blastocatellia    1
## 9             c__Chloroflexia    4
## 10              c__Deinococci    1
## 11        c__Desulfovibrionia    1
## 12        c__Desulfuromonadia    1
## 13     c__Gammaproteobacteria   62
## 14          c__Gemmatimonadia   30
## 15          c__Incertae_Sedis    1
## 16            c__JG30-KF-CM66    1
## 17         c__Ktedonobacteria    3
## 18              c__Myxococcia    1
## 19             c__Oligoflexia    2
## 20          c__Planctomycetes    1
## 21              c__Polyangiia    1
## 22 c__S0134_terrestrial_group    1
## 23         c__Thermoleophilia    7
## 24                    c__TK10    3
## 25        c__Vampirivibrionia    1
## 26        c__Verrucomicrobiia    3
## 27        c__Vicinamibacteria    8
## 28                       <NA>    1
```

``` r
Acc_Bac_TwoTimes_DA_ASV_Tax_Table$VL$Class
```

```
##                           col Freq
## 1           c__Actinobacteria   20
## 2      c__Alphaproteobacteria   16
## 3                  c__Bacilli   10
## 4              c__Bacteroidia    2
## 5             c__Chloroflexia    3
## 6               c__Clostridia    2
## 7      c__Gammaproteobacteria   20
## 8           c__Gemmatimonadia   18
## 9              c__Oligoflexia    1
## 10          c__Planctomycetes    1
## 11 c__S0134_terrestrial_group    1
## 12         c__Thermoleophilia   14
## 13                    c__TK10    1
## 14        c__Verrucomicrobiia    1
## 15        c__Vicinamibacteria    3
## 16                       <NA>    1
```

``` r
save(Acc_Bac_TwoTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_TwoTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring three times

``` r
# Filter to ASVs detected by at least 3 tests
Acc_Bac_NumbOcc_DA_ASV <- lapply(Acc_Bac_Total_DA_ASV, function(x) table(unlist(x)))
Acc_Bac_ThreeTimes_DA_ASV <- lapply(Acc_Bac_NumbOcc_DA_ASV, function(x) names(x[x >= 3]))
Acc_Bac_ThreeTimes_DA_ASV
```

```
## $VL
##  [1] "bASV_1028" "bASV_1153" "bASV_126"  "bASV_1345" "bASV_1539" "bASV_178" 
##  [7] "bASV_1962" "bASV_2120" "bASV_2578" "bASV_27"   "bASV_290"  "bASV_303" 
## [13] "bASV_329"  "bASV_385"  "bASV_386"  "bASV_417"  "bASV_7"    "bASV_732" 
## [19] "bASV_888"  "bASV_968" 
## 
## $CD
##  [1] "bASV_1055" "bASV_141"  "bASV_1652" "bASV_214"  "bASV_302"  "bASV_340" 
##  [7] "bASV_345"  "bASV_426"  "bASV_444"  "bASV_729"  "bASV_900" 
## 
## $RI
## [1] "bASV_136"  "bASV_1525" "bASV_1956" "bASV_90"   "bASV_980" 
## 
## $HM
##  [1] "bASV_1"    "bASV_1002" "bASV_109"  "bASV_1107" "bASV_1108" "bASV_1163"
##  [7] "bASV_1164" "bASV_1231" "bASV_126"  "bASV_1325" "bASV_1354" "bASV_1366"
## [13] "bASV_1446" "bASV_148"  "bASV_1522" "bASV_1618" "bASV_1665" "bASV_178" 
## [19] "bASV_182"  "bASV_183"  "bASV_190"  "bASV_1934" "bASV_202"  "bASV_2195"
## [25] "bASV_2258" "bASV_228"  "bASV_2290" "bASV_234"  "bASV_2378" "bASV_243" 
## [31] "bASV_2655" "bASV_268"  "bASV_288"  "bASV_3"    "bASV_326"  "bASV_3348"
## [37] "bASV_345"  "bASV_363"  "bASV_376"  "bASV_3920" "bASV_394"  "bASV_3962"
## [43] "bASV_42"   "bASV_423"  "bASV_4395" "bASV_45"   "bASV_483"  "bASV_491" 
## [49] "bASV_496"  "bASV_5"    "bASV_511"  "bASV_515"  "bASV_541"  "bASV_566" 
## [55] "bASV_5795" "bASV_62"   "bASV_6504" "bASV_664"  "bASV_670"  "bASV_697" 
## [61] "bASV_7"    "bASV_73"   "bASV_753"  "bASV_831"  "bASV_928"  "bASV_930" 
## [67] "bASV_995" 
## 
## $GO1
## [1] "bASV_1003" "bASV_1060" "bASV_2612" "bASV_385"  "bASV_426"  "bASV_46"  
## [7] "bASV_477"
```

``` r
# Number of DA ASVs occurring in three tests
lapply(Acc_Bac_ThreeTimes_DA_ASV, function(x) length(x))
```

```
## $VL
## [1] 20
## 
## $CD
## [1] 11
## 
## $RI
## [1] 5
## 
## $HM
## [1] 67
## 
## $GO1
## [1] 7
```

``` r
save(Acc_Bac_ThreeTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_ThreeTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
Acc_Bac_ThreeTimes_DA_ASV_Tax <- lapply(Acc_Bac_ThreeTimes_DA_ASV, function(x) FP_unnormalized_bac_ps_Tax[x, ])

save(Acc_Bac_ThreeTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_ThreeTimes_DA_ASV_Tax.RData")

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
##                       col Freq
## 1       c__Acidimicrobiia    2
## 2       c__Actinobacteria   16
## 3  c__Alphaproteobacteria   18
## 4              c__Bacilli    1
## 5          c__Bacteroidia    1
## 6      c__Bdellovibrionia    1
## 7  c__Gammaproteobacteria   12
## 8       c__Gemmatimonadia   12
## 9      c__Ktedonobacteria    1
## 10    c__Verrucomicrobiia    2
## 11    c__Vicinamibacteria    1
```

``` r
Acc_Bac_ThreeTimes_DA_ASV_Tax_Table$VL$Class
```

```
##                      col Freq
## 1      c__Actinobacteria    7
## 2 c__Alphaproteobacteria    4
## 3        c__Chloroflexia    2
## 4 c__Gammaproteobacteria    3
## 5      c__Gemmatimonadia    2
## 6      c__Planctomycetes    1
## 7     c__Thermoleophilia    1
```

``` r
save(Acc_Bac_ThreeTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_ThreeTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring four times

``` r
# Filter to ASVs detected by at least 4 tests
Acc_Bac_FourTimes_DA_ASV <- lapply(Acc_Bac_NumbOcc_DA_ASV, function(x) names(x[x >= 4]))
Acc_Bac_FourTimes_DA_ASV
```

```
## $VL
##  [1] "bASV_1028" "bASV_1153" "bASV_126"  "bASV_1345" "bASV_1539" "bASV_178" 
##  [7] "bASV_1962" "bASV_2120" "bASV_2578" "bASV_27"   "bASV_329"  "bASV_385" 
## [13] "bASV_417"  "bASV_7"    "bASV_888" 
## 
## $CD
## [1] "bASV_1055" "bASV_141"  "bASV_1652" "bASV_214"  "bASV_302"  "bASV_340" 
## [7] "bASV_426"  "bASV_444"  "bASV_900" 
## 
## $RI
## [1] "bASV_1525" "bASV_90"   "bASV_980" 
## 
## $HM
##  [1] "bASV_1"    "bASV_1002" "bASV_109"  "bASV_1107" "bASV_1108" "bASV_1163"
##  [7] "bASV_1164" "bASV_1231" "bASV_126"  "bASV_1325" "bASV_1354" "bASV_1366"
## [13] "bASV_1446" "bASV_178"  "bASV_182"  "bASV_183"  "bASV_1934" "bASV_2195"
## [19] "bASV_2290" "bASV_234"  "bASV_243"  "bASV_2655" "bASV_268"  "bASV_288" 
## [25] "bASV_3"    "bASV_326"  "bASV_3348" "bASV_363"  "bASV_394"  "bASV_42"  
## [31] "bASV_423"  "bASV_45"   "bASV_483"  "bASV_496"  "bASV_515"  "bASV_541" 
## [37] "bASV_566"  "bASV_5795" "bASV_62"   "bASV_6504" "bASV_664"  "bASV_7"   
## [43] "bASV_753"  "bASV_831"  "bASV_928"  "bASV_995" 
## 
## $GO1
## [1] "bASV_1003" "bASV_426"  "bASV_477"
```

``` r
# Number of DA ASVs occurring in four tests
lapply(Acc_Bac_FourTimes_DA_ASV, function(x) length(x))
```

```
## $VL
## [1] 15
## 
## $CD
## [1] 9
## 
## $RI
## [1] 3
## 
## $HM
## [1] 46
## 
## $GO1
## [1] 3
```

``` r
save(Acc_Bac_FourTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_FourTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
Acc_Bac_FourTimes_DA_ASV_Tax <- lapply(Acc_Bac_FourTimes_DA_ASV, function(x) FP_unnormalized_bac_ps_Tax[x, ])

save(Acc_Bac_FourTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_FourTimes_DA_ASV_Tax.RData")

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
##                      col Freq
## 1      c__Actinobacteria   13
## 2 c__Alphaproteobacteria   11
## 3             c__Bacilli    1
## 4         c__Bacteroidia    1
## 5     c__Bdellovibrionia    1
## 6 c__Gammaproteobacteria   10
## 7      c__Gemmatimonadia    8
## 8    c__Verrucomicrobiia    1
```

``` r
Acc_Bac_FourTimes_DA_ASV_Tax_Table$VL$Class
```

```
##                      col Freq
## 1      c__Actinobacteria    7
## 2 c__Alphaproteobacteria    4
## 3        c__Chloroflexia    1
## 4 c__Gammaproteobacteria    1
## 5      c__Gemmatimonadia    1
## 6      c__Planctomycetes    1
```

``` r
save(Acc_Bac_FourTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_FourTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring five times

``` r
# Filter to ASVs detected by at least 5 tests
Acc_Bac_FiveTimes_DA_ASV <- lapply(Acc_Bac_NumbOcc_DA_ASV, function(x) names(x[x >= 5]))
Acc_Bac_FiveTimes_DA_ASV
```

```
## $VL
## [1] "bASV_126" "bASV_178" "bASV_27"  "bASV_329" "bASV_7"  
## 
## $CD
## [1] "bASV_214" "bASV_302" "bASV_426" "bASV_444"
## 
## $RI
## character(0)
## 
## $HM
##  [1] "bASV_1"    "bASV_109"  "bASV_1107" "bASV_1108" "bASV_126"  "bASV_178" 
##  [7] "bASV_182"  "bASV_183"  "bASV_234"  "bASV_243"  "bASV_268"  "bASV_3"   
## [13] "bASV_326"  "bASV_363"  "bASV_394"  "bASV_42"   "bASV_483"  "bASV_515" 
## [19] "bASV_541"  "bASV_566"  "bASV_62"   "bASV_7"    "bASV_831"  "bASV_928" 
## [25] "bASV_995" 
## 
## $GO1
## [1] "bASV_426"
```

``` r
# Number of DA ASVs occurring in Five tests
lapply(Acc_Bac_FiveTimes_DA_ASV, function(x) length(x))
```

```
## $VL
## [1] 5
## 
## $CD
## [1] 4
## 
## $RI
## [1] 0
## 
## $HM
## [1] 25
## 
## $GO1
## [1] 1
```

``` r
save(Acc_Bac_FiveTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_FiveTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
Acc_Bac_FiveTimes_DA_ASV_Tax <- lapply(Acc_Bac_FiveTimes_DA_ASV, function(x) FP_unnormalized_bac_ps_Tax[x, ])

save(Acc_Bac_FiveTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_FiveTimes_DA_ASV_Tax.RData")

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
##                      col Freq
## 1      c__Actinobacteria    9
## 2 c__Alphaproteobacteria    7
## 3 c__Gammaproteobacteria    5
## 4      c__Gemmatimonadia    4
```

``` r
Acc_Bac_FiveTimes_DA_ASV_Tax_Table$VL$Class
```

```
##                      col Freq
## 1      c__Actinobacteria    4
## 2 c__Alphaproteobacteria    1
```

``` r
save(Acc_Bac_FiveTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_FiveTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring six times

``` r
# Filter to ASVs detected by at least 6 tests
Acc_Bac_SixTimes_DA_ASV <- lapply(Acc_Bac_NumbOcc_DA_ASV, function(x) names(x[x >= 6]))
Acc_Bac_SixTimes_DA_ASV
```

```
## $VL
## [1] "bASV_126"
## 
## $CD
## character(0)
## 
## $RI
## character(0)
## 
## $HM
##  [1] "bASV_109" "bASV_126" "bASV_178" "bASV_182" "bASV_183" "bASV_234"
##  [7] "bASV_243" "bASV_268" "bASV_326" "bASV_363" "bASV_394" "bASV_42" 
## [13] "bASV_483" "bASV_515" "bASV_541" "bASV_566" "bASV_7"   "bASV_831"
## [19] "bASV_928"
## 
## $GO1
## character(0)
```

``` r
# Number of DA ASVs occurring in Six tests
lapply(Acc_Bac_SixTimes_DA_ASV, function(x) length(x))
```

```
## $VL
## [1] 1
## 
## $CD
## [1] 0
## 
## $RI
## [1] 0
## 
## $HM
## [1] 19
## 
## $GO1
## [1] 0
```

``` r
save(Acc_Bac_SixTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_SixTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
Acc_Bac_SixTimes_DA_ASV_Tax <- lapply(Acc_Bac_SixTimes_DA_ASV, function(x) FP_unnormalized_bac_ps_Tax[x, ])

save(Acc_Bac_SixTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_SixTimes_DA_ASV_Tax.RData")

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
## $VL
## $VL$Phylum
##                 col Freq
## 1 p__Actinomycetota    1
## 
## $VL$Class
##                 col Freq
## 1 c__Actinobacteria    1
## 
## $VL$Order
##                      col Freq
## 1 o__Streptosporangiales    1
## 
## $VL$Family
##    col Freq
## 1 <NA>    1
## 
## $VL$Genus
##    col Freq
## 1 <NA>    1
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
## $HM
## $HM$Phylum
##                  col Freq
## 1  p__Actinomycetota    7
## 2 p__Gemmatimonadota    2
## 3  p__Pseudomonadota   10
## 
## $HM$Class
##                      col Freq
## 1      c__Actinobacteria    7
## 2 c__Alphaproteobacteria    7
## 3 c__Gammaproteobacteria    3
## 4      c__Gemmatimonadia    2
## 
## $HM$Order
##                      col Freq
## 1      o__Azospirillales    1
## 2     o__Burkholderiales    2
## 3    o__Gemmatimonadales    2
## 4    o__Hyphomicrobiales    6
## 5    o__Kitasatosporales    1
## 6      o__Lysobacterales    1
## 7   o__Pseudonocardiales    1
## 8 o__Streptosporangiales    5
## 
## $HM$Family
##                       col Freq
## 1      f__Azospirillaceae    1
## 2     f__Beijerinckiaceae    3
## 3     f__Burkholderiaceae    1
## 4    f__Gemmatimonadaceae    2
## 5     f__Oxalobacteraceae    1
## 6   f__Pseudonocardiaceae    1
## 7         f__Rhizobiaceae    3
## 8   f__Rhodanobacteraceae    1
## 9    f__Streptomycetaceae    1
## 10 f__Thermomonosporaceae    2
## 11                   <NA>    3
## 
## $HM$Genus
##                                              col Freq
## 1                             g__Actinoallomurus    1
## 2                                g__Actinomadura    1
## 3                                g__Azospirillum    1
## 4  g__Burkholderia-Caballeronia-Paraburkholderia    1
## 5                                g__Gemmatimonas    2
## 6                                    g__Massilia    1
## 7                            g__Methylobacterium    1
## 8                               g__Methylorosula    1
## 9                              g__Pseudonocardia    1
## 10                                  g__Rhizobium    1
## 11                                 g__Tahibacter    1
## 12                                          <NA>    7
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
save(Acc_Bac_SixTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_SixTimes_DA_ASV_Tax_Table.RData")
```

``` r
# Clean environment
rm(list = ls())
```

## 5.3.2 Accession Family level
### Load DA Family

``` r
## Load outputs DESeq
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_FamilyLevel/sigtab_stddds2_Wald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_FamilyLevel/sigtab_stddds2_LRT_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_FamilyLevel/sigtab_zinbWald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_FamilyLevel/sigtab_zinbLRT_Acc_Bac.RData")

## Load output Ancom
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_FamilyLevel/fdr_ancomWZ_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_FamilyLevel/fdr_ancomNoZ_Acc_Bac.RData")

## Load Phyloseq object
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_unnormalized_bac_ps_ForDA_clean.RData")
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
tax <- as.data.frame(tax_table(FP_unnormalized_bac_ps_ForDA_clean))
tax$ASV <- rownames(tax)

Acc_NoDup_DESeq_DA_ASV_Family_Tax <- lapply(Acc_DESeq_DA_ASV_Family_l_t, function(x) {
  lapply(x, function(y) {
    unique(tax[tax$ASV %in% y, c("Family", "Order", "Class", "Phylum")])
    })})
```

#### Acnom

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
      ) %>% select(Phylum, Class, Order, Family)
    
    # Extract taxonomy from ps object for short names
    tax_ps <- as.data.frame(tax_table(FP_unnormalized_bac_ps_ForDA_clean))
    #family_only <- sub("^f__", "", family_only_strings)
    
    tax_family_map <- unique(tax_ps[tax_ps$Family %in% family_only_strings, c("Family", "Order", "Class", "Phylum")])
    #nrow(tax_family_map)
    
    # Remove left over Incertae Sedis from short names object
    tax_family_map <- tax_family_map %>%
      filter(Family != "f__Incertae_Sedis")
   
    # Combine long and short names
    tax_combined <- as.data.frame(bind_rows(tax_full, tax_family_map))
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
## $VL
## [1] 66
## 
## $CD
## [1] 57
## 
## $RI
## [1] 70
## 
## $HM
## [1] 248
## 
## $GO1
## [1] 71
```

``` r
lapply(Acc_Total_Merged_DA_Family, function(x) table(x$Family))
```

```
## $VL
## 
##                              f__A0839                     f__Alcaligenaceae 
##                                     1                                     1 
##                             f__B1-7BS                            f__BIrii41 
##                                     1                                     1 
##              f__Caldalkalibacillaceae                     f__Caldilineaceae 
##                                     1                                     1 
##                f__Desulfotomaculaceae                        f__Elsteraceae 
##                                     2                                     1 
##                 f__Enterobacteriaceae                f__Erysipelotrichaceae 
##                                     1                                     1 
##                   f__Fibrobacteraceae                        f__Frankiaceae 
##                                     1                                     2 
##                    f__Geminicoccaceae                   f__Halanaerobiaceae 
##                                     1                                     2 
##                 f__Herpetosiphonaceae                     f__Incertae_Sedis 
##                                     1                                     7 
##                     f__Isosphaeraceae                       f__JG30-KF-CM45 
##                                     2                                     3 
##                               f__LWQ8                f__Magnetospirillaceae 
##                                     2                                     1 
##                        f__Nevskiaceae                      f__Paracoccaceae 
##                                     1                                     1 
##              f__Peptostreptococcaceae                     f__Rhodocyclaceae 
##                                     1                                     1 
##                  f__Rhodomicrobiaceae                    f__Sandaracinaceae 
##                                     1                                     1 
##                  f__Sphingomonadaceae                    f__Spirosomataceae 
##                                     3                                     1 
##                f__Steroidobacteraceae f__Streptosporangiales_Incertae_Sedis 
##                                     1                                     2 
##                      f__Sumerlaeaceae             f__Thermoanaerobaculaceae 
##                                     1                                     1 
##                 f__Thermomicrobiaceae                f__Thermomonosporaceae 
##                                     5                                     4 
##                       f__Unclassified 
##                                     9 
## 
## $CD
## 
##       f__Anaerolineaceae   f__Aneurinibacillaceae                f__B1-7BS 
##                        1                        1                        1 
##      f__Cellvibrionaceae        f__Comamonadaceae f__Desulfitobacteriaceae 
##                        5                        1                        1 
##   f__Desulfotomaculaceae      f__Fibrobacteraceae           f__Gemmataceae 
##                        1                        1                        1 
##       f__Hyphomonadaceae        f__Incertae_Sedis         f__Inquilinaceae 
##                        1                       12                        1 
##           f__JG30-KF-AS9    f__Ktedonobacteraceae   f__Magnetospirillaceae 
##                        5                        6                        1 
##        f__Micrococcaceae      f__Oscillospiraceae  f__Paracaedibacteraceae 
##                        5                        1                        1 
##      f__Phaselicystaceae      f__Pseudomonadaceae        f__Rhodocyclaceae 
##                        1                        1                        1 
##                f__SM2D12       f__Spirochaetaceae         f__Sumerlaeaceae 
##                        1                        2                        1 
##     f__Terrimicrobiaceae          f__Unclassified                 f__WWH38 
##                        1                        2                        1 
## 
## $RI
## 
##  f__Caldalkalibacillaceae      f__Cellulomonadaceae       f__Cellvibrionaceae 
##                         1                         2                         1 
##     f__Chromobacteriaceae   f__CPla-3_termite_group    f__Desulfotomaculaceae 
##                         1                         1                         1 
##    f__Desulfovibrionaceae             f__env.OPS_17    f__Erysipelotrichaceae 
##                         1                         1                         1 
##       f__Fibrobacteraceae      f__Fimbriimonadaceae          f__Garciellaceae 
##                         1                         1                         1 
##        f__Geminicoccaceae        f__Glycomycetaceae     f__Herpetosiphonaceae 
##                         1                         1                         1 
##         f__Incertae_Sedis          f__Inquilinaceae         f__Isosphaeraceae 
##                        13                         1                         1 
##            f__Kaistiaceae        f__Koribacteraceae         f__Legionellaceae 
##                         1                         1                         2 
##        f__Microscillaceae          f__Myxococcaceae     f__Obscuribacteraceae 
##                         1                         1                         2 
##         f__Oligoflexaceae       f__Oscillospiraceae      f__Parachlamydiaceae 
##                         1                         2                         1 
##          f__Paracoccaceae          f__Polyangiaceae     f__Rhodanobacteraceae 
##                         1                         1                         2 
##         f__Rhodocyclaceae      f__Rickettsiellaceae         f__Rubritaleaceae 
##                         3                         1                         1 
##        f__Spirosomataceae    f__Steroidobacteraceae       f__Sulfobacillaceae 
##                         1                         1                         1 
##          f__Sumerlaeaceae      f__Symbiobacteraceae      f__Terrimicrobiaceae 
##                         2                         1                         1 
##    f__Thermacetogeniaceae f__Thermoanaerobaculaceae           f__Unclassified 
##                         1                         1                         8 
##               f__URHD0088 
##                         1 
## 
## $HM
## 
##                   f__Acetobacteraceae                  f__Acidimicrobiaceae 
##                                     6                                     6 
##     f__Acidobacteriaceae_(Subgroup_1)                            f__AKIW781 
##                                     2                                     2 
##                     f__Alcaligenaceae              f__Anaeromyxobacteraceae 
##                                     3                                     2 
##                    f__Azospirillaceae                             f__B1-7BS 
##                                     6                                     1 
##                   f__Beijerinckiaceae                            f__BIrii41 
##                                     6                                     2 
##                    f__Bryobacteraceae                   f__Burkholderiaceae 
##                                     3                                     2 
##              f__Caldalkalibacillaceae                   f__Caulobacteraceae 
##                                     1                                     6 
##                  f__Cellulomonadaceae                   f__Chitinophagaceae 
##                                     3                                     5 
##                f__Christensenellaceae                     f__Clostridiaceae 
##                                     1                                     2 
##                       f__Coxiellaceae  f__Desulfotomaculales_Incertae_Sedis 
##                                     1                                     1 
##                f__Desulfovibrionaceae                 f__Enterobacteriaceae 
##                                     2                                     1 
##                   f__Fibrobacteraceae                  f__Flavobacteriaceae 
##                                     1                                     1 
##                     f__Geobacteraceae                f__Geodermatophilaceae 
##                                     1                                     5 
##                    f__Glycomycetaceae                   f__Halanaerobiaceae 
##                                     1                                     1 
##                 f__Hydrogenophilaceae                  f__Hymenobacteraceae 
##                                     2                                     1 
##                    f__Hyphomonadaceae                     f__Incertae_Sedis 
##                                     1                                    43 
##                      f__Inquilinaceae                     f__Isosphaeraceae 
##                                     2                                     4 
##                         f__KF-JG30-B3                    f__Lysobacteraceae 
##                                     1                                     6 
##                   f__Magnetospiraceae                f__Magnetospirillaceae 
##                                     2                                     4 
##                   f__Micavibrionaceae                     f__Micrococcaceae 
##                                     1                                     6 
##                     f__Micropepsaceae                        f__Opitutaceae 
##                                     5                                     1 
##                   f__Oscillospiraceae                   f__Paenibacillaceae 
##                                     1                                     3 
##                   f__Phaselicystaceae                      f__Polyangiaceae 
##                                     1                                     3 
##           f__Pseudobdellovibrionaceae                   f__Pseudomonadaceae 
##                                     3                                     4 
##                 f__Pseudonocardiaceae                     f__Reyranellaceae 
##                                     5                                     2 
##                 f__Rhodanobacteraceae                     f__Rhodocyclaceae 
##                                     6                                     5 
##                  f__Rickettsiellaceae                     f__Roseiflexaceae 
##                                     1                                     3 
##                     f__Rubritaleaceae                    f__Sandaracinaceae 
##                                     5                                     1 
##                    f__Solibacteraceae               f__Solirubrobacteraceae 
##                                     6                                     2 
##                    f__Spirosomataceae                f__Steroidobacteraceae 
##                                     1                                     1 
## f__Streptosporangiales_Incertae_Sedis                  f__Symbiobacteraceae 
##                                     1                                     1 
##                  f__Terrimicrobiaceae                f__Thermacetogeniaceae 
##                                     1                                     1 
##                 f__Thermomicrobiaceae                f__Thermomonosporaceae 
##                                     3                                     6 
##                       f__Unclassified                           f__URHD0088 
##                                    13                                     4 
##                f__Verrucomicrobiaceae                              f__WWH38 
##                                     5                                     1 
##                  f__Xanthobacteraceae             f__Xiphinematobacteraceae 
##                                     2                                     4 
## 
## $GO1
## 
##             f__AKIW781            f__AKYG1722             f__BIrii41 
##                      1                      1                      1 
##      f__Caldilineaceae       f__Chlamydiaceae        f__Coxiellaceae 
##                      1                      1                      1 
## f__Desulfovibrionaceae          f__env.OPS_17    f__Fibrobacteraceae 
##                      1                      1                      1 
##   f__Fimbriimonadaceae       f__Holosporaceae     f__Hyphomonadaceae 
##                      1                      1                      1 
##      f__Incertae_Sedis  f__Intrasporangiaceae        f__JG30-KF-CM45 
##                     21                      2                      2 
##         f__Kaistiaceae  f__Ktedonobacteraceae    f__Micavibrionaceae 
##                      2                      2                      1 
##   f__Microbacteriaceae       f__Moraxellaceae     f__Nocardioidaceae 
##                      1                      1                      5 
##    f__Oscillospiraceae       f__Paracoccaceae     f__Parvibaculaceae 
##                      1                      2                      1 
##   f__Rhodomicrobiaceae     f__Sandaracinaceae f__Sphingobacteriaceae 
##                      1                      1                      2 
##     f__Sporichthyaceae f__Steroidobacteraceae   f__Terrimicrobiaceae 
##                      3                      1                      1 
##  f__Thermomicrobiaceae        f__Unclassified f__Vampirovibrionaceae 
##                      2                      4                      1 
##   f__Xanthobacteraceae 
##                      2
```

``` r
save(Acc_Total_Merged_DA_Family, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_FamilyLevel/Acc_Total_Merged_DA_Family.RData")
```

### Remove duplicat Family

``` r
Acc_NoDup_Merged_DA_Family <- lapply(Acc_Total_Merged_DA_Family, function(x) unique(x))
lapply(Acc_NoDup_Merged_DA_Family, function(x) nrow(x))
```

```
## $VL
## [1] 44
## 
## $CD
## [1] 37
## 
## $RI
## [1] 58
## 
## $HM
## [1] 98
## 
## $GO1
## [1] 54
```

``` r
lapply(Acc_NoDup_Merged_DA_Family, function(x) table(x$Family))
```

```
## $VL
## 
##                              f__A0839                     f__Alcaligenaceae 
##                                     1                                     1 
##                             f__B1-7BS                            f__BIrii41 
##                                     1                                     1 
##              f__Caldalkalibacillaceae                     f__Caldilineaceae 
##                                     1                                     1 
##                f__Desulfotomaculaceae                        f__Elsteraceae 
##                                     1                                     1 
##                 f__Enterobacteriaceae                f__Erysipelotrichaceae 
##                                     1                                     1 
##                   f__Fibrobacteraceae                        f__Frankiaceae 
##                                     1                                     1 
##                    f__Geminicoccaceae                   f__Halanaerobiaceae 
##                                     1                                     1 
##                 f__Herpetosiphonaceae                     f__Incertae_Sedis 
##                                     1                                     7 
##                     f__Isosphaeraceae                       f__JG30-KF-CM45 
##                                     1                                     1 
##                               f__LWQ8                f__Magnetospirillaceae 
##                                     1                                     1 
##                        f__Nevskiaceae                      f__Paracoccaceae 
##                                     1                                     1 
##              f__Peptostreptococcaceae                     f__Rhodocyclaceae 
##                                     1                                     1 
##                  f__Rhodomicrobiaceae                    f__Sandaracinaceae 
##                                     1                                     1 
##                  f__Sphingomonadaceae                    f__Spirosomataceae 
##                                     1                                     1 
##                f__Steroidobacteraceae f__Streptosporangiales_Incertae_Sedis 
##                                     1                                     1 
##                      f__Sumerlaeaceae             f__Thermoanaerobaculaceae 
##                                     1                                     1 
##                 f__Thermomicrobiaceae                f__Thermomonosporaceae 
##                                     1                                     1 
##                       f__Unclassified 
##                                     4 
## 
## $CD
## 
##       f__Anaerolineaceae   f__Aneurinibacillaceae                f__B1-7BS 
##                        1                        1                        1 
##      f__Cellvibrionaceae        f__Comamonadaceae f__Desulfitobacteriaceae 
##                        1                        1                        1 
##   f__Desulfotomaculaceae      f__Fibrobacteraceae           f__Gemmataceae 
##                        1                        1                        1 
##       f__Hyphomonadaceae        f__Incertae_Sedis         f__Inquilinaceae 
##                        1                       10                        1 
##           f__JG30-KF-AS9    f__Ktedonobacteraceae   f__Magnetospirillaceae 
##                        1                        1                        1 
##        f__Micrococcaceae      f__Oscillospiraceae  f__Paracaedibacteraceae 
##                        1                        1                        1 
##      f__Phaselicystaceae      f__Pseudomonadaceae        f__Rhodocyclaceae 
##                        1                        1                        1 
##                f__SM2D12       f__Spirochaetaceae         f__Sumerlaeaceae 
##                        1                        1                        1 
##     f__Terrimicrobiaceae          f__Unclassified                 f__WWH38 
##                        1                        2                        1 
## 
## $RI
## 
##  f__Caldalkalibacillaceae      f__Cellulomonadaceae       f__Cellvibrionaceae 
##                         1                         1                         1 
##     f__Chromobacteriaceae   f__CPla-3_termite_group    f__Desulfotomaculaceae 
##                         1                         1                         1 
##    f__Desulfovibrionaceae             f__env.OPS_17    f__Erysipelotrichaceae 
##                         1                         1                         1 
##       f__Fibrobacteraceae      f__Fimbriimonadaceae          f__Garciellaceae 
##                         1                         1                         1 
##        f__Geminicoccaceae        f__Glycomycetaceae     f__Herpetosiphonaceae 
##                         1                         1                         1 
##         f__Incertae_Sedis          f__Inquilinaceae         f__Isosphaeraceae 
##                        11                         1                         1 
##            f__Kaistiaceae        f__Koribacteraceae         f__Legionellaceae 
##                         1                         1                         1 
##        f__Microscillaceae          f__Myxococcaceae     f__Obscuribacteraceae 
##                         1                         1                         1 
##         f__Oligoflexaceae       f__Oscillospiraceae      f__Parachlamydiaceae 
##                         1                         1                         1 
##          f__Paracoccaceae          f__Polyangiaceae     f__Rhodanobacteraceae 
##                         1                         1                         1 
##         f__Rhodocyclaceae      f__Rickettsiellaceae         f__Rubritaleaceae 
##                         1                         1                         1 
##        f__Spirosomataceae    f__Steroidobacteraceae       f__Sulfobacillaceae 
##                         1                         1                         1 
##          f__Sumerlaeaceae      f__Symbiobacteraceae      f__Terrimicrobiaceae 
##                         1                         1                         1 
##    f__Thermacetogeniaceae f__Thermoanaerobaculaceae           f__Unclassified 
##                         1                         1                         6 
##               f__URHD0088 
##                         1 
## 
## $HM
## 
##                   f__Acetobacteraceae                  f__Acidimicrobiaceae 
##                                     1                                     1 
##     f__Acidobacteriaceae_(Subgroup_1)                            f__AKIW781 
##                                     1                                     1 
##                     f__Alcaligenaceae              f__Anaeromyxobacteraceae 
##                                     1                                     1 
##                    f__Azospirillaceae                             f__B1-7BS 
##                                     1                                     1 
##                   f__Beijerinckiaceae                            f__BIrii41 
##                                     1                                     1 
##                    f__Bryobacteraceae                   f__Burkholderiaceae 
##                                     1                                     1 
##              f__Caldalkalibacillaceae                   f__Caulobacteraceae 
##                                     1                                     1 
##                  f__Cellulomonadaceae                   f__Chitinophagaceae 
##                                     1                                     1 
##                f__Christensenellaceae                     f__Clostridiaceae 
##                                     1                                     1 
##                       f__Coxiellaceae  f__Desulfotomaculales_Incertae_Sedis 
##                                     1                                     1 
##                f__Desulfovibrionaceae                 f__Enterobacteriaceae 
##                                     1                                     1 
##                   f__Fibrobacteraceae                  f__Flavobacteriaceae 
##                                     1                                     1 
##                     f__Geobacteraceae                f__Geodermatophilaceae 
##                                     1                                     1 
##                    f__Glycomycetaceae                   f__Halanaerobiaceae 
##                                     1                                     1 
##                 f__Hydrogenophilaceae                  f__Hymenobacteraceae 
##                                     1                                     1 
##                    f__Hyphomonadaceae                     f__Incertae_Sedis 
##                                     1                                    22 
##                      f__Inquilinaceae                     f__Isosphaeraceae 
##                                     1                                     1 
##                         f__KF-JG30-B3                    f__Lysobacteraceae 
##                                     1                                     1 
##                   f__Magnetospiraceae                f__Magnetospirillaceae 
##                                     1                                     1 
##                   f__Micavibrionaceae                     f__Micrococcaceae 
##                                     1                                     1 
##                     f__Micropepsaceae                        f__Opitutaceae 
##                                     1                                     1 
##                   f__Oscillospiraceae                   f__Paenibacillaceae 
##                                     1                                     1 
##                   f__Phaselicystaceae                      f__Polyangiaceae 
##                                     1                                     1 
##           f__Pseudobdellovibrionaceae                   f__Pseudomonadaceae 
##                                     1                                     1 
##                 f__Pseudonocardiaceae                     f__Reyranellaceae 
##                                     1                                     1 
##                 f__Rhodanobacteraceae                     f__Rhodocyclaceae 
##                                     1                                     1 
##                  f__Rickettsiellaceae                     f__Roseiflexaceae 
##                                     1                                     1 
##                     f__Rubritaleaceae                    f__Sandaracinaceae 
##                                     1                                     1 
##                    f__Solibacteraceae               f__Solirubrobacteraceae 
##                                     1                                     1 
##                    f__Spirosomataceae                f__Steroidobacteraceae 
##                                     1                                     1 
## f__Streptosporangiales_Incertae_Sedis                  f__Symbiobacteraceae 
##                                     1                                     1 
##                  f__Terrimicrobiaceae                f__Thermacetogeniaceae 
##                                     1                                     1 
##                 f__Thermomicrobiaceae                f__Thermomonosporaceae 
##                                     1                                     1 
##                       f__Unclassified                           f__URHD0088 
##                                     6                                     1 
##                f__Verrucomicrobiaceae                              f__WWH38 
##                                     1                                     1 
##                  f__Xanthobacteraceae             f__Xiphinematobacteraceae 
##                                     1                                     1 
## 
## $GO1
## 
##             f__AKIW781            f__AKYG1722             f__BIrii41 
##                      1                      1                      1 
##      f__Caldilineaceae       f__Chlamydiaceae        f__Coxiellaceae 
##                      1                      1                      1 
## f__Desulfovibrionaceae          f__env.OPS_17    f__Fibrobacteraceae 
##                      1                      1                      1 
##   f__Fimbriimonadaceae       f__Holosporaceae     f__Hyphomonadaceae 
##                      1                      1                      1 
##      f__Incertae_Sedis  f__Intrasporangiaceae        f__JG30-KF-CM45 
##                     18                      1                      1 
##         f__Kaistiaceae  f__Ktedonobacteraceae    f__Micavibrionaceae 
##                      1                      1                      1 
##   f__Microbacteriaceae       f__Moraxellaceae     f__Nocardioidaceae 
##                      1                      1                      1 
##    f__Oscillospiraceae       f__Paracoccaceae     f__Parvibaculaceae 
##                      1                      1                      1 
##   f__Rhodomicrobiaceae     f__Sandaracinaceae f__Sphingobacteriaceae 
##                      1                      1                      1 
##     f__Sporichthyaceae f__Steroidobacteraceae   f__Terrimicrobiaceae 
##                      1                      1                      1 
##  f__Thermomicrobiaceae        f__Unclassified f__Vampirovibrionaceae 
##                      1                      4                      1 
##   f__Xanthobacteraceae 
##                      1
```

``` r
save(Acc_NoDup_Merged_DA_Family, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_FamilyLevel/Acc_NoDup_Merged_DA_Family.RData")
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
## $VL
## [1] 11
## 
## $CD
## [1] 7
## 
## $RI
## [1] 11
## 
## $HM
## [1] 54
## 
## $GO1
## [1] 13
```

``` r
save(Acc_TwoTimes_Merged_DA_Family, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_FamilyLevel/Acc_TwoTimes_Merged_DA_Family.RData")
```

``` r
# Clean environment
rm(list = ls())
```


## 5.3.3 Accession Order level
### Load DA Order

``` r
## Load outputs DESeq
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_OrderLevel/sigtab_stddds2_Wald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_OrderLevel/sigtab_stddds2_LRT_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_OrderLevel/sigtab_zinbWald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_OrderLevel/sigtab_zinbLRT_Acc_Bac.RData")

## Load output Ancom
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_OrderLevel/fdr_ancomWZ_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_OrderLevel/fdr_ancomNoZ_Acc_Bac.RData")

## Load Phyloseq object
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_unnormalized_bac_ps_ForDA_clean.RData")
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
tax <- as.data.frame(tax_table(FP_unnormalized_bac_ps_ForDA_clean))
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
      ) %>% select(Phylum, Class, Order)
    
    # Extract taxonomy from ps object for short names
    tax_ps <- as.data.frame(tax_table(FP_unnormalized_bac_ps_ForDA_clean))
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
## $VL
## [1] 41
## 
## $CD
## [1] 29
## 
## $RI
## [1] 53
## 
## $HM
## [1] 153
## 
## $GO1
## [1] 44
```

``` r
lapply(Acc_Total_Merged_DA_Order, function(x) table(x$Order))
```

```
## $VL
## 
##                             o__Blfdi19                o__Caldalkalibacillales 
##                                      1                                      1 
##                       o__Caldilineales                  o__Desulfotomaculales 
##                                      1                                      2 
##                    o__Enterobacterales                  o__Erysipelotrichales 
##                                      1                                      1 
##                     o__Fibrobacterales                     o__Halanaerobiales 
##                                      1                                      2 
##                      o__Incertae_Sedis                       o__Isosphaerales 
##                                      5                                      3 
## o__Peptostreptococcales-Tissierellales                     o__Rhodobacterales 
##                                      1                                      1 
##                    o__Salinisphaerales                    o__Sphingomonadales 
##                                      1                                      2 
##                  o__Steroidobacterales                 o__Streptosporangiales 
##                                      1                                      6 
##                        o__Sumerlaeales               o__Thermoanaerobaculales 
##                                      1                                      1 
##                   o__Thermomicrobiales                        o__Tistrellales 
##                                      6                                      1 
##                        o__Unclassified 
##                                      2 
## 
## $CD
## 
##               o__Anaerolineales           o__Aneurinibacillales 
##                               1                               1 
##                     o__B10-SB3A                      o__Blfdi19 
##                               2                               1 
## o__Candidatus_Zambryskibacteria                        o__CCD24 
##                               1                               1 
##                 o__Chlamydiales         o__Desulfitobacteriales 
##                               1                               1 
##              o__Fibrobacterales                   o__Gemmatales 
##                               1                               1 
##               o__Incertae_Sedis            o__Ktedonobacterales 
##                               4                               6 
##              o__Oscillospirales          o__Paracaedibacterales 
##                               1                               1 
##             o__Rhodospirillales               o__Spirochaetales 
##                               2                               2 
##                 o__Sumerlaeales                 o__Unclassified 
##                               1                               1 
## 
## $RI
## 
##            o__Babeliales               o__Blfdi19  o__Caldalkalibacillales 
##                        2                        2                        1 
##          o__Chlamydiales    o__Desulfotomaculales    o__Desulfovibrionales 
##                        1                        1                        1 
##            o__Elsterales    o__Erysipelotrichales         o__Eubacteriales 
##                        2                        1                        1 
##       o__Fibrobacterales      o__Fimbriimonadales        o__Glycomycetales 
##                        1                        1                        1 
##        o__Incertae_Sedis         o__Isosphaerales         o__Legionellales 
##                        5                        1                        2 
##     o__Obscuribacterales         o__Oligoflexales       o__Oscillospirales 
##                        2                        2                        2 
##       o__Paenibacillales       o__Rhodobacterales      o__Rhodospirillales 
##                        2                        1                        2 
##      o__Rickettsiellales       o__Rokubacteriales    o__Steroidobacterales 
##                        1                        2                        1 
##           o__Subgroup_17            o__Subgroup_2            o__Subgroup_7 
##                        1                        1                        2 
##       o__Sulfobacillales          o__Sumerlaeales     o__Symbiobacteriales 
##                        1                        2                        1 
##    o__Thermacetogeniales o__Thermoanaerobaculales          o__Tistrellales 
##                        1                        2                        1 
##          o__Unclassified 
##                        3 
## 
## $HM
## 
##      o__Acetobacterales     o__Acidimicrobiales       o__Azospirillales 
##                       6                       6                       6 
##             o__B10-SB3A           o__Bacillales       o__Bryobacterales 
##                       1                       2                       2 
## o__Caldalkalibacillales      o__Caulobacterales                o__CCD24 
##                       1                       6                       4 
##      o__Chitinophagales         o__Chlamydiales       o__Chloroflexales 
##                       5                       1                       3 
##   o__Christensenellales        o__Clostridiales          o__Coxiellales 
##                       1                       2                       2 
##   o__Desulfovibrionales     o__Enterobacterales      o__Fibrobacterales 
##                       2                       1                       1 
##        o__Geobacterales       o__Glycomycetales      o__Halanaerobiales 
##                       1                       1                       1 
##       o__Incertae_Sedis        o__Isosphaerales        o__Kallotenuales 
##                      33                       3                       2 
##      o__Kapabacteriales           o__Lineage_IV       o__Lysobacterales 
##                       1                       1                       6 
##        o__Micropepsales           o__Opitutales      o__Oscillospirales 
##                       3                       1                       1 
##      o__Paenibacillales         o__Polyangiales      o__Pseudomonadales 
##                       3                       6                       5 
##    o__Pseudonocardiales     o__Rhodospirillales     o__Rickettsiellales 
##                       3                       1                       1 
##    o__Saccharimonadales       o__Solibacterales   o__Steroidobacterales 
##                       2                       6                       1 
##  o__Streptosporangiales          o__Subgroup_17           o__Subgroup_2 
##                       6                       1                       1 
##           o__Subgroup_7    o__Symbiobacteriales   o__Thermacetogeniales 
##                       1                       2                       1 
##         o__Unclassified   o__Verrucomicrobiales 
##                       1                       6 
## 
## $GO1
## 
##                  o__B10-SB3A                   o__Blfdi19 
##                            1                            1 
##             o__Caldilineales o__Candidatus_Kaiserbacteria 
##                            1                            2 
##              o__Chlamydiales               o__Coxiellales 
##                            1                            1 
##          o__Defluviicoccales        o__Desulfovibrionales 
##                            1                            1 
##                 o__FFCH16263           o__Fibrobacterales 
##                            1                            1 
##          o__Fimbriimonadales              o__Holosporales 
##                            1                            1 
##            o__Incertae_Sedis             o__Kallotenuales 
##                            9                            1 
##                o__Lineage_IV            o__Lysobacterales 
##                            1                            2 
##           o__Oscillospirales            o__Parvibaculales 
##                            1                            1 
##       o__Propionibacteriales           o__Rhodobacterales 
##                            5                            2 
##             o__Rickettsiales                      o__S085 
##                            1                            1 
##        o__Steroidobacterales         o__Thermomicrobiales 
##                            1                            2 
##              o__Unclassified        o__Vampirovibrionales 
##                            2                            1 
##             o__Zavarziniales 
##                            1
```

``` r
save(Acc_Total_Merged_DA_Order, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_OrderLevel/Acc_Total_Merged_DA_Order.RData")
```

### Remove duplicat Order

``` r
Acc_NoDup_Merged_DA_Order <- lapply(Acc_Total_Merged_DA_Order, function(x) unique(x))
lapply(Acc_NoDup_Merged_DA_Order, function(x) nrow(x))
```

```
## $VL
## [1] 26
## 
## $CD
## [1] 21
## 
## $RI
## [1] 40
## 
## $HM
## [1] 59
## 
## $GO1
## [1] 34
```

``` r
lapply(Acc_NoDup_Merged_DA_Order, function(x) table(x$Order))
```

```
## $VL
## 
##                             o__Blfdi19                o__Caldalkalibacillales 
##                                      1                                      1 
##                       o__Caldilineales                  o__Desulfotomaculales 
##                                      1                                      1 
##                    o__Enterobacterales                  o__Erysipelotrichales 
##                                      1                                      1 
##                     o__Fibrobacterales                     o__Halanaerobiales 
##                                      1                                      1 
##                      o__Incertae_Sedis                       o__Isosphaerales 
##                                      5                                      1 
## o__Peptostreptococcales-Tissierellales                     o__Rhodobacterales 
##                                      1                                      1 
##                    o__Salinisphaerales                    o__Sphingomonadales 
##                                      1                                      1 
##                  o__Steroidobacterales                 o__Streptosporangiales 
##                                      1                                      1 
##                        o__Sumerlaeales               o__Thermoanaerobaculales 
##                                      1                                      1 
##                   o__Thermomicrobiales                        o__Tistrellales 
##                                      1                                      1 
##                        o__Unclassified 
##                                      2 
## 
## $CD
## 
##               o__Anaerolineales           o__Aneurinibacillales 
##                               1                               1 
##                     o__B10-SB3A                      o__Blfdi19 
##                               1                               1 
## o__Candidatus_Zambryskibacteria                        o__CCD24 
##                               1                               1 
##                 o__Chlamydiales         o__Desulfitobacteriales 
##                               1                               1 
##              o__Fibrobacterales                   o__Gemmatales 
##                               1                               1 
##               o__Incertae_Sedis            o__Ktedonobacterales 
##                               4                               1 
##              o__Oscillospirales          o__Paracaedibacterales 
##                               1                               1 
##             o__Rhodospirillales               o__Spirochaetales 
##                               1                               1 
##                 o__Sumerlaeales                 o__Unclassified 
##                               1                               1 
## 
## $RI
## 
##            o__Babeliales               o__Blfdi19  o__Caldalkalibacillales 
##                        1                        1                        1 
##          o__Chlamydiales    o__Desulfotomaculales    o__Desulfovibrionales 
##                        1                        1                        1 
##            o__Elsterales    o__Erysipelotrichales         o__Eubacteriales 
##                        1                        1                        1 
##       o__Fibrobacterales      o__Fimbriimonadales        o__Glycomycetales 
##                        1                        1                        1 
##        o__Incertae_Sedis         o__Isosphaerales         o__Legionellales 
##                        5                        1                        1 
##     o__Obscuribacterales         o__Oligoflexales       o__Oscillospirales 
##                        1                        1                        1 
##       o__Paenibacillales       o__Rhodobacterales      o__Rhodospirillales 
##                        1                        1                        1 
##      o__Rickettsiellales       o__Rokubacteriales    o__Steroidobacterales 
##                        1                        1                        1 
##           o__Subgroup_17            o__Subgroup_2            o__Subgroup_7 
##                        1                        1                        1 
##       o__Sulfobacillales          o__Sumerlaeales     o__Symbiobacteriales 
##                        1                        1                        1 
##    o__Thermacetogeniales o__Thermoanaerobaculales          o__Tistrellales 
##                        1                        1                        1 
##          o__Unclassified 
##                        3 
## 
## $HM
## 
##      o__Acetobacterales     o__Acidimicrobiales       o__Azospirillales 
##                       1                       1                       1 
##             o__B10-SB3A           o__Bacillales       o__Bryobacterales 
##                       1                       1                       1 
## o__Caldalkalibacillales      o__Caulobacterales                o__CCD24 
##                       1                       1                       1 
##      o__Chitinophagales         o__Chlamydiales       o__Chloroflexales 
##                       1                       1                       1 
##   o__Christensenellales        o__Clostridiales          o__Coxiellales 
##                       1                       1                       1 
##   o__Desulfovibrionales     o__Enterobacterales      o__Fibrobacterales 
##                       1                       1                       1 
##        o__Geobacterales       o__Glycomycetales      o__Halanaerobiales 
##                       1                       1                       1 
##       o__Incertae_Sedis        o__Isosphaerales        o__Kallotenuales 
##                      13                       1                       1 
##      o__Kapabacteriales           o__Lineage_IV       o__Lysobacterales 
##                       1                       1                       1 
##        o__Micropepsales           o__Opitutales      o__Oscillospirales 
##                       1                       1                       1 
##      o__Paenibacillales         o__Polyangiales      o__Pseudomonadales 
##                       1                       1                       1 
##    o__Pseudonocardiales     o__Rhodospirillales     o__Rickettsiellales 
##                       1                       1                       1 
##    o__Saccharimonadales       o__Solibacterales   o__Steroidobacterales 
##                       1                       1                       1 
##  o__Streptosporangiales          o__Subgroup_17           o__Subgroup_2 
##                       1                       1                       1 
##           o__Subgroup_7    o__Symbiobacteriales   o__Thermacetogeniales 
##                       1                       1                       1 
##         o__Unclassified   o__Verrucomicrobiales 
##                       1                       1 
## 
## $GO1
## 
##                  o__B10-SB3A                   o__Blfdi19 
##                            1                            1 
##             o__Caldilineales o__Candidatus_Kaiserbacteria 
##                            1                            1 
##              o__Chlamydiales               o__Coxiellales 
##                            1                            1 
##          o__Defluviicoccales        o__Desulfovibrionales 
##                            1                            1 
##                 o__FFCH16263           o__Fibrobacterales 
##                            1                            1 
##          o__Fimbriimonadales              o__Holosporales 
##                            1                            1 
##            o__Incertae_Sedis             o__Kallotenuales 
##                            7                            1 
##                o__Lineage_IV            o__Lysobacterales 
##                            1                            1 
##           o__Oscillospirales            o__Parvibaculales 
##                            1                            1 
##       o__Propionibacteriales           o__Rhodobacterales 
##                            1                            1 
##             o__Rickettsiales                      o__S085 
##                            1                            1 
##        o__Steroidobacterales         o__Thermomicrobiales 
##                            1                            1 
##              o__Unclassified        o__Vampirovibrionales 
##                            2                            1 
##             o__Zavarziniales 
##                            1
```

``` r
save(Acc_NoDup_Merged_DA_Order, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_OrderLevel/Acc_NoDup_Merged_DA_Order.RData")
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
## $VL
## [1] 6
## 
## $CD
## [1] 4
## 
## $RI
## [1] 13
## 
## $HM
## [1] 32
## 
## $GO1
## [1] 7
```

``` r
save(Acc_TwoTimes_Merged_DA_Order, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_OrderLevel/Acc_TwoTimes_Merged_DA_Order.RData")
```

``` r
# Clean environment
rm(list = ls())
```

## 5.3.4 Accession Class level
### Load DA Class

``` r
## Load outputs DESeq
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_ClassLevel/sigtab_stddds2_Wald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_ClassLevel/sigtab_stddds2_LRT_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_ClassLevel/sigtab_zinbWald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_ClassLevel/sigtab_zinbLRT_Acc_Bac.RData")

## Load output Ancom
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ClassLevel/fdr_ancomWZ_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ClassLevel/fdr_ancomNoZ_Acc_Bac.RData")

## Load Phyloseq object
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_unnormalized_bac_ps_ForDA_clean.RData")
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
tax <- as.data.frame(tax_table(FP_unnormalized_bac_ps_ForDA_clean))
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
    
    # Extract taxomomy from long names
    tax_full <- tibble(Taxon = full_tax_strings) %>%
      mutate(
        Phylum = str_extract(Taxon, "p__.*?(?=_c__|$)"),
        Class  = str_extract(Taxon, "c__.*?(?=_o__|$)")
      ) %>% select(Phylum, Class)
    
    # Extract taxonomy from ps object for short names
    tax_ps <- as.data.frame(tax_table(FP_unnormalized_bac_ps_ForDA_clean))
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
## $VL
## [1] 19
## 
## $CD
## [1] 16
## 
## $RI
## [1] 24
## 
## $HM
## [1] 78
## 
## $GO1
## [1] 12
```

``` r
lapply(Acc_Total_Merged_DA_Class, function(x) table(x$Class))
```

```
## $VL
## 
##        c__Chloroflexia    c__Desulfotomaculia       c__Fibrobacteria 
##                      5                      2                      1 
##       c__Halanaerobiia      c__Incertae_Sedis               c__OLB14 
##                      2                      2                      1 
##      c__Planctomycetes   c__Sericytochromatia          c__Sumerlaeia 
##                      2                      1                      1 
## c__Thermoanaerobaculia        c__Unclassified 
##                      1                      1 
## 
## $CD
## 
##         c__Chlamydiia c__Desulfitobacteriia     c__Eremiobacteria 
##                     1                     1                     1 
##      c__Fibrobacteria     c__Incertae_Sedis    c__Ktedonobacteria 
##                     1                     2                     6 
##       c__Spirochaetia         c__Sumerlaeia       c__Unclassified 
##                     2                     1                     1 
## 
## $RI
## 
##                     c__AD3                c__Babeliae 
##                          1                          2 
##              c__Chlamydiia        c__Desulfotomaculia 
##                          1                          1 
##        c__Desulfovibrionia          c__Eremiobacteria 
##                          1                          1 
##           c__Fibrobacteria          c__Fimbriimonadia 
##                          1                          1 
##                 c__Kazania        c__Methylomirabilia 
##                          1                          1 
##          c__Planctomycetes              c__Polyangiia 
##                          1                          2 
## c__S0134_terrestrial_group           c__Sulfobacillia 
##                          1                          1 
##              c__Sumerlaeia         c__Symbiobacteriia 
##                          2                          1 
##         c__Thermacetogenia     c__Thermoanaerobaculia 
##                          1                          1 
##            c__Unclassified        c__Vampirivibrionia 
##                          1                          2 
## 
## $HM
## 
##          c__Acidimicrobiia          c__Acidobacteriae 
##                          5                          5 
##          c__Actinobacteria                     c__AD3 
##                          6                          1 
##                 c__Bacilli             c__Bacteroidia 
##                          3                          1 
##         c__Bdellovibrionia              c__Chlamydiia 
##                          2                          1 
##              c__Clostridia        c__Desulfovibrionia 
##                          6                          4 
##        c__Desulfuromonadia           c__Elusimicrobia 
##                          1                          1 
##           c__Fibrobacteria             c__Gitt-GS-136 
##                          1                          3 
##         c__Gracilibacteria           c__Halanaerobiia 
##                          1                          1 
##          c__Incertae_Sedis            c__Kapabacteria 
##                          7                          1 
##                 c__Kazania                   c__OLB14 
##                          1                          1 
##           c__Parcubacteria          c__Planctomycetes 
##                          1                          1 
##              c__Polyangiia c__S0134_terrestrial_group 
##                          3                          5 
##         c__Saccharimonadia         c__Symbiobacteriia 
##                          3                          2 
##         c__Thermacetogenia                    c__TK10 
##                          1                          6 
##            c__Unclassified        c__Verrucomicrobiia 
##                          1                          3 
## 
## $GO1
## 
##             c__Anaerolineae c__BD2-11_terrestrial_group 
##                           1                           1 
##               c__Chlamydiia          c__Dehalococcoidia 
##                           1                           1 
##         c__Desulfovibrionia            c__Elusimicrobia 
##                           1                           1 
##            c__Fibrobacteria           c__Fimbriimonadia 
##                           1                           1 
##           c__Incertae_Sedis                  c__Kazania 
##                           2                           1 
##                    c__OLB14 
##                           1
```

``` r
save(Acc_Total_Merged_DA_Class, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ClassLevel/Acc_Total_Merged_DA_Class.RData")
```

### Remove duplicat Class

``` r
Acc_NoDup_Merged_DA_Class <- lapply(Acc_Total_Merged_DA_Class, function(x) unique(x))
lapply(Acc_NoDup_Merged_DA_Class, function(x) nrow(x))
```

```
## $VL
## [1] 12
## 
## $CD
## [1] 10
## 
## $RI
## [1] 20
## 
## $HM
## [1] 32
## 
## $GO1
## [1] 12
```

``` r
lapply(Acc_NoDup_Merged_DA_Class, function(x) table(x$Class))
```

```
## $VL
## 
##        c__Chloroflexia    c__Desulfotomaculia       c__Fibrobacteria 
##                      1                      1                      1 
##       c__Halanaerobiia      c__Incertae_Sedis               c__OLB14 
##                      1                      2                      1 
##      c__Planctomycetes   c__Sericytochromatia          c__Sumerlaeia 
##                      1                      1                      1 
## c__Thermoanaerobaculia        c__Unclassified 
##                      1                      1 
## 
## $CD
## 
##         c__Chlamydiia c__Desulfitobacteriia     c__Eremiobacteria 
##                     1                     1                     1 
##      c__Fibrobacteria     c__Incertae_Sedis    c__Ktedonobacteria 
##                     1                     2                     1 
##       c__Spirochaetia         c__Sumerlaeia       c__Unclassified 
##                     1                     1                     1 
## 
## $RI
## 
##                     c__AD3                c__Babeliae 
##                          1                          1 
##              c__Chlamydiia        c__Desulfotomaculia 
##                          1                          1 
##        c__Desulfovibrionia          c__Eremiobacteria 
##                          1                          1 
##           c__Fibrobacteria          c__Fimbriimonadia 
##                          1                          1 
##                 c__Kazania        c__Methylomirabilia 
##                          1                          1 
##          c__Planctomycetes              c__Polyangiia 
##                          1                          1 
## c__S0134_terrestrial_group           c__Sulfobacillia 
##                          1                          1 
##              c__Sumerlaeia         c__Symbiobacteriia 
##                          1                          1 
##         c__Thermacetogenia     c__Thermoanaerobaculia 
##                          1                          1 
##            c__Unclassified        c__Vampirivibrionia 
##                          1                          1 
## 
## $HM
## 
##          c__Acidimicrobiia          c__Acidobacteriae 
##                          1                          1 
##          c__Actinobacteria                     c__AD3 
##                          1                          1 
##                 c__Bacilli             c__Bacteroidia 
##                          1                          1 
##         c__Bdellovibrionia              c__Chlamydiia 
##                          1                          1 
##              c__Clostridia        c__Desulfovibrionia 
##                          1                          1 
##        c__Desulfuromonadia           c__Elusimicrobia 
##                          1                          1 
##           c__Fibrobacteria             c__Gitt-GS-136 
##                          1                          1 
##         c__Gracilibacteria           c__Halanaerobiia 
##                          1                          1 
##          c__Incertae_Sedis            c__Kapabacteria 
##                          3                          1 
##                 c__Kazania                   c__OLB14 
##                          1                          1 
##           c__Parcubacteria          c__Planctomycetes 
##                          1                          1 
##              c__Polyangiia c__S0134_terrestrial_group 
##                          1                          1 
##         c__Saccharimonadia         c__Symbiobacteriia 
##                          1                          1 
##         c__Thermacetogenia                    c__TK10 
##                          1                          1 
##            c__Unclassified        c__Verrucomicrobiia 
##                          1                          1 
## 
## $GO1
## 
##             c__Anaerolineae c__BD2-11_terrestrial_group 
##                           1                           1 
##               c__Chlamydiia          c__Dehalococcoidia 
##                           1                           1 
##         c__Desulfovibrionia            c__Elusimicrobia 
##                           1                           1 
##            c__Fibrobacteria           c__Fimbriimonadia 
##                           1                           1 
##           c__Incertae_Sedis                  c__Kazania 
##                           2                           1 
##                    c__OLB14 
##                           1
```

``` r
save(Acc_NoDup_Merged_DA_Class, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ClassLevel/Acc_NoDup_Merged_DA_Class.RData")
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
## $VL
## [1] 4
## 
## $CD
## [1] 2
## 
## $RI
## [1] 4
## 
## $HM
## [1] 16
## 
## $GO1
## [1] 0
```

``` r
save(Acc_TwoTimes_Merged_DA_Class, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ClassLevel/Acc_TwoTimes_Merged_DA_Class.RData")
```

``` r
# Clean environment
rm(list = ls())
```


## 5.3.5 Accession Phylum level
### Load DA Phylum

``` r
## Load outputs DESeq
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_PhylumLevel/sigtab_stddds2_Wald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_PhylumLevel/sigtab_stddds2_LRT_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_PhylumLevel/sigtab_zinbWald_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_PhylumLevel/sigtab_zinbLRT_Acc_Bac.RData")

## Load output Ancom
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_PhylumLevel/fdr_ancomWZ_Acc_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_PhylumLevel/fdr_ancomNoZ_Acc_Bac.RData")

## Load Phyloseq object
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_unnormalized_bac_ps_ForDA_clean.RData")
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
tax <- as.data.frame(tax_table(FP_unnormalized_bac_ps_ForDA_clean))
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
      ) %>% select(Phylum, Class)
    
    # Extract taxonomy from ps object for short names
    tax_ps <- as.data.frame(tax_table(FP_unnormalized_bac_ps_ForDA_clean))
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
## $VL
## [1] 11
## 
## $CD
## [1] 45
## 
## $RI
## [1] 17
## 
## $HM
## [1] 65
## 
## $GO1
## [1] 5
```

``` r
lapply(Acc_Total_Merged_DA_Phylum, function(x) table(x))
```

```
## $VL
## x
## p__Bdellovibrionota   p__Fibrobacterota p__Halanaerobiaeota      p__Sumerlaeota 
##                   6                   1                   2                   1 
##              p__WS2 
##                   1 
## 
## $CD
## x
##               p__Armatimonadota             p__Bdellovibrionota 
##                               2                               1 
##   p__Candidatus_Eremiobacterota                  p__Chlamydiota 
##                               1                               1 
##                p__Chloroflexota               p__Fibrobacterota 
##                              30                               1 
##                 p__Nitrospirota                      p__RCP2-54 
##                               2                               1 
## p__SAR324_clade(Marine_group_B)                p__Spirochaetota 
##                               1                               4 
##                  p__Sumerlaeota 
##                               1 
## 
## $RI
## x
##             p__Armatimonadota p__Candidatus_Eremiobacterota 
##                             2                             1 
##                p__Chlamydiota               p__Dependentiae 
##                             1                             2 
##             p__Fibrobacterota          p__Methylomirabilota 
##                             1                             2 
##                p__Myxococcota                p__Sumerlaeota 
##                             6                             2 
## 
## $HM
## x
##         p__Acidobacteriota          p__Actinomycetota 
##                         19                         10 
##               p__Bacillota            p__Bacteroidota 
##                          2                          2 
## p__Candidatus_Kapabacteria             p__Chlamydiota 
##                          1                          1 
##         p__Elusimicrobiota          p__Fibrobacterota 
##                          3                          1 
##        p__Halanaerobiaeota         p__Patescibacteria 
##                          1                         13 
##                 p__RCP2-54 p__Thermodesulfobacteriota 
##                          3                          6 
##            p__Unclassified       p__Verrucomicrobiota 
##                          1                          2 
## 
## $GO1
## x
##               p__Armatimonadota                  p__Chlamydiota 
##                               2                               1 
##               p__Fibrobacterota p__SAR324_clade(Marine_group_B) 
##                               1                               1
```

``` r
save(Acc_Total_Merged_DA_Phylum, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_PhylumLevel/Acc_Total_Merged_DA_Phylum.RData")
```

### Remove duplicat Phylum

``` r
Acc_NoDup_Merged_DA_Phylum <- lapply(Acc_Total_Merged_DA_Phylum, function(x) unique(x))
lapply(Acc_NoDup_Merged_DA_Phylum, function(x) length(x))
```

```
## $VL
## [1] 5
## 
## $CD
## [1] 11
## 
## $RI
## [1] 8
## 
## $HM
## [1] 14
## 
## $GO1
## [1] 4
```

``` r
lapply(Acc_NoDup_Merged_DA_Phylum, function(x) table(x))
```

```
## $VL
## x
## p__Bdellovibrionota   p__Fibrobacterota p__Halanaerobiaeota      p__Sumerlaeota 
##                   1                   1                   1                   1 
##              p__WS2 
##                   1 
## 
## $CD
## x
##               p__Armatimonadota             p__Bdellovibrionota 
##                               1                               1 
##   p__Candidatus_Eremiobacterota                  p__Chlamydiota 
##                               1                               1 
##                p__Chloroflexota               p__Fibrobacterota 
##                               1                               1 
##                 p__Nitrospirota                      p__RCP2-54 
##                               1                               1 
## p__SAR324_clade(Marine_group_B)                p__Spirochaetota 
##                               1                               1 
##                  p__Sumerlaeota 
##                               1 
## 
## $RI
## x
##             p__Armatimonadota p__Candidatus_Eremiobacterota 
##                             1                             1 
##                p__Chlamydiota               p__Dependentiae 
##                             1                             1 
##             p__Fibrobacterota          p__Methylomirabilota 
##                             1                             1 
##                p__Myxococcota                p__Sumerlaeota 
##                             1                             1 
## 
## $HM
## x
##         p__Acidobacteriota          p__Actinomycetota 
##                          1                          1 
##               p__Bacillota            p__Bacteroidota 
##                          1                          1 
## p__Candidatus_Kapabacteria             p__Chlamydiota 
##                          1                          1 
##         p__Elusimicrobiota          p__Fibrobacterota 
##                          1                          1 
##        p__Halanaerobiaeota         p__Patescibacteria 
##                          1                          1 
##                 p__RCP2-54 p__Thermodesulfobacteriota 
##                          1                          1 
##            p__Unclassified       p__Verrucomicrobiota 
##                          1                          1 
## 
## $GO1
## x
##               p__Armatimonadota                  p__Chlamydiota 
##                               1                               1 
##               p__Fibrobacterota p__SAR324_clade(Marine_group_B) 
##                               1                               1
```

``` r
save(Acc_NoDup_Merged_DA_Phylum, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_PhylumLevel/Acc_NoDup_Merged_DA_Phylum.RData")
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
## $VL
## [1] 2
## 
## $CD
## [1] 4
## 
## $RI
## [1] 5
## 
## $HM
## [1] 9
## 
## $GO1
## [1] 1
```

``` r
save(Acc_TwoTimes_Merged_DA_Phylum, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_PhylumLevel/Acc_TwoTimes_Merged_DA_Phylum.RData")
```

``` r
# Clean environment
rm(list = ls())
```


## 5.3.6 All data ASV level
### Load DA ASVs

``` r
## Load outputs DESeq
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_ASVLevel/sigtab_stddds2_Wald_All_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_ASVLevel/sigtab_stddds2_LRT_All_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_ASVLevel/sigtab_zinbWald_All_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_DESeq_ASVLevel/sigtab_zinbLRT_All_Bac.RData")

## Load output Ancom
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ASVLevel/fdr_ancomWZ_All_Bac.RData")
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_Ancom_ASVLevel/fdr_ancomNoZ_All_Bac.RData")
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
## [1] 891
```

``` r
save(All_Bac_Total_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_Total_DA_ASV.RData")

# Remove duplicated asvs names from that list
All_Bac_NoDup_DA_ASV <- unique(All_Bac_Total_DA_ASV)

# Number of DA ASVs without duplicates
length(All_Bac_NoDup_DA_ASV)
```

```
## [1] 887
```

``` r
save(All_Bac_NoDup_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_NoDup_DA_ASV.RData")

# Extracting taxonomy
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_unnormalized_bac_ps.RData")

# Extract taxonomy of DA ASVs
FP_unnormalized_bac_ps_Tax <- as.data.frame(tax_table(FP_unnormalized_bac_ps))
All_Bac_NoDup_DA_ASV_Tax <- FP_unnormalized_bac_ps_Tax[All_Bac_NoDup_DA_ASV, ]

save(All_Bac_NoDup_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects//FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_NoDup_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
All_Bac_NoDup_DA_ASV_Tax_Table <- lapply(All_Bac_NoDup_DA_ASV_Tax[ , 2:6], function(df) table(df))
All_Bac_NoDup_DA_ASV_Tax_Table <- lapply(All_Bac_NoDup_DA_ASV_Tax_Table, as.data.frame)

# Examples
All_Bac_NoDup_DA_ASV_Tax_Table$Class
```

```
##                            df Freq
## 1           c__Acidimicrobiia   18
## 2           c__Acidobacteriae   20
## 3           c__Actinobacteria   68
## 4      c__Alphaproteobacteria   91
## 5             c__Anaerolineae    1
## 6                  c__Bacilli   87
## 7          c__Bacteriovoracia    1
## 8              c__Bacteroidia   41
## 9          c__Bdellovibrionia    4
## 10          c__Berkelbacteria    1
## 11          c__Blastocatellia   12
## 12            c__Chloroflexia    7
## 13              c__Clostridia    3
## 14         c__Dehalococcoidia    1
## 15        c__Desulfotomaculia    1
## 16        c__Desulfuromonadia    2
## 17          c__Fimbriimonadia    1
## 18     c__Gammaproteobacteria  227
## 19          c__Gemmatimonadia   50
## 20             c__Gitt-GS-136    2
## 21              c__Holophagae    3
## 22          c__Incertae_Sedis    3
## 23            c__JG30-KF-CM66    2
## 24            c__Kapabacteria    2
## 25                 c__Kazania    2
## 26                  c__KD4-96    8
## 27         c__Ktedonobacteria   25
## 28            c__Limnochordia    7
## 29               c__MB-A2-108    4
## 30              c__Myxococcia    3
## 31             c__Oligoflexia   14
## 32           c__Parcubacteria    5
## 33           c__Phycisphaerae    8
## 34              c__Polyangiia   15
## 35 c__S0134_terrestrial_group    2
## 36         c__Saccharimonadia    6
## 37       c__Sericytochromatia    2
## 38              c__Subgroup_5    1
## 39         c__Symbiobacteriia    2
## 40     c__Thermoanaerobaculia    1
## 41         c__Thermoleophilia   65
## 42                    c__TK10    3
## 43        c__Verrucomicrobiia   19
## 44        c__Vicinamibacteria   45
```

``` r
save(All_Bac_NoDup_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_NoDup_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring twice

``` r
# Filter to ASVs detected by at least 2 tests
All_Bac_NumbOcc_DA_ASV <- table(unlist(All_Bac_Total_DA_ASV))
All_Bac_TwoTimes_DA_ASV <- names(All_Bac_NumbOcc_DA_ASV[All_Bac_NumbOcc_DA_ASV >= 2])
All_Bac_TwoTimes_DA_ASV
```

```
## [1] "bASV_1164" "bASV_141"  "bASV_2195" "bASV_458"
```

``` r
# Number of DA ASVs occurring in two tests
length(All_Bac_TwoTimes_DA_ASV)
```

```
## [1] 4
```

``` r
save(All_Bac_TwoTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_TwoTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
All_Bac_TwoTimes_DA_ASV_Tax <- FP_unnormalized_bac_ps_Tax[All_Bac_TwoTimes_DA_ASV, ]

save(All_Bac_TwoTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects//FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_TwoTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
All_Bac_TwoTimes_DA_ASV_Tax_Table <- lapply(All_Bac_TwoTimes_DA_ASV_Tax[ , 2:6], function(df) table(df))
All_Bac_TwoTimes_DA_ASV_Tax_Table <- lapply(All_Bac_TwoTimes_DA_ASV_Tax_Table, as.data.frame)

# Examples
All_Bac_TwoTimes_DA_ASV_Tax_Table$Class
```

```
##                       df Freq
## 1      c__Acidimicrobiia    1
## 2      c__Actinobacteria    1
## 3 c__Alphaproteobacteria    1
## 4             c__Bacilli    1
```

``` r
save(All_Bac_TwoTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_TwoTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring three times

``` r
# Filter to ASVs detected by at least 3 tests
All_Bac_ThreeTimes_DA_ASV <- names(All_Bac_NumbOcc_DA_ASV[All_Bac_NumbOcc_DA_ASV >= 3])
All_Bac_ThreeTimes_DA_ASV
```

```
## character(0)
```

``` r
# Number of DA ASVs occurring in three tests
length(All_Bac_ThreeTimes_DA_ASV)
```

```
## [1] 0
```

``` r
save(All_Bac_ThreeTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_ThreeTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
All_Bac_ThreeTimes_DA_ASV_Tax <- FP_unnormalized_bac_ps_Tax[All_Bac_ThreeTimes_DA_ASV, ]

save(All_Bac_ThreeTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects//FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_ThreeTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
All_Bac_ThreeTimes_DA_ASV_Tax_Table <- lapply(All_Bac_ThreeTimes_DA_ASV_Tax[ , 2:6], function(df) table(df))
All_Bac_ThreeTimes_DA_ASV_Tax_Table <- lapply(All_Bac_ThreeTimes_DA_ASV_Tax_Table, as.data.frame)

# Examples
All_Bac_ThreeTimes_DA_ASV_Tax_Table
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
save(All_Bac_ThreeTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_ThreeTimes_DA_ASV_Tax_Table.RData")
```

### DA ASVs occuring four times

``` r
# Filter to ASVs detected by at least 3 tests
All_Bac_FourTimes_DA_ASV <- names(All_Bac_NumbOcc_DA_ASV[All_Bac_NumbOcc_DA_ASV >= 4])
All_Bac_FourTimes_DA_ASV
```

```
## character(0)
```

``` r
# Number of DA ASVs occurring in three tests
length(All_Bac_FourTimes_DA_ASV)
```

```
## [1] 0
```

``` r
save(All_Bac_FourTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_FourTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
All_Bac_FourTimes_DA_ASV_Tax <- FP_unnormalized_bac_ps_Tax[All_Bac_FourTimes_DA_ASV, ]

save(All_Bac_FourTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects//FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_FourTimes_DA_ASV_Tax.RData")

# Make a table of occurrence at each taxonomic level
All_Bac_FourTimes_DA_ASV_Tax_Table <- lapply(All_Bac_FourTimes_DA_ASV_Tax[ , 2:6], function(df) table(df))
All_Bac_FourTimes_DA_ASV_Tax_Table <- lapply(All_Bac_FourTimes_DA_ASV_Tax_Table, as.data.frame)

# Examples
All_Bac_FourTimes_DA_ASV_Tax_Table
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
save(All_Bac_FourTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_FourTimes_DA_ASV_Tax_Table.RData")
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
save(All_Bac_FiveTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_FiveTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
All_Bac_FiveTimes_DA_ASV_Tax <- FP_unnormalized_bac_ps_Tax[All_Bac_FiveTimes_DA_ASV, ]

save(All_Bac_FiveTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects//FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_FiveTimes_DA_ASV_Tax.RData")

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
save(All_Bac_FiveTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_FiveTimes_DA_ASV_Tax_Table.RData")
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
save(All_Bac_SixTimes_DA_ASV, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_SixTimes_DA_ASV.RData")

# Extract taxonomy of DA ASVs
All_Bac_SixTimes_DA_ASV_Tax <- FP_unnormalized_bac_ps_Tax[All_Bac_SixTimes_DA_ASV, ]

save(All_Bac_SixTimes_DA_ASV_Tax, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects//FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_SixTimes_DA_ASV_Tax.RData")

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
save(All_Bac_SixTimes_DA_ASV_Tax_Table, file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/All_Bac_SixTimes_DA_ASV_Tax_Table.RData")
```

``` r
# Clean environment
rm(list = ls())
```


## 5.3.5 All data Class and Phylum level
No DA classes and phyla detected when analysing whole data.
