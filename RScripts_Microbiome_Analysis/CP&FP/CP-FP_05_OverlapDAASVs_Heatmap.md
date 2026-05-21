---
title: "CP_FP_05 - Heatmaps of ASVs that are DA in CP and in FP"
author: "Kris de Kreek"
date: "2026-05-21"
output: 
  html_document:
    toc: true
    keep_md: true
editor_options: 
  chunk_output_type: console
---

   
Tutorials   
- [Heatmap for final plot](http://rstudio-pubs-static.s3.amazonaws.com/288398_185f2889a5f641c6b9aa7b14fa15b634.html)
- [Excamples of pheatmap package](https://davetang.github.io/muse/pheatmap.html)   
   
   
# 5.0 load libraries
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
library(pheatmap)
packageVersion("pheatmap")
```

```
## [1] '1.0.13'
```


# 5.1 Heatmap overlapping ASVs Acc
## 5.1.1 Load data

``` r
# load phyloseq object CP-FP together
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/Data/R_objects/CP_FP_together/unnormalized_bac_ps.RData")

# load DA ASVs CP for Acc
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/Data/R_objects/CP_DA_SummaryFiles_ASVLevel/Acc_Bac_TwoTimes_DA_ASV.RData")
Acc_Bac_TwoTimes_DA_ASV_CP <- Acc_Bac_TwoTimes_DA_ASV
rm(Acc_Bac_TwoTimes_DA_ASV)

# load DA ASVs FP for Acc
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/Data/R_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_TwoTimes_DA_ASV.RData")
Acc_Bac_TwoTimes_DA_ASV_FP <- Acc_Bac_TwoTimes_DA_ASV
rm(Acc_Bac_TwoTimes_DA_ASV)

# ASVs overlapping in DA tests CP and FP for Acc
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/Data/R_objects/CP_FP_OverlapDA/Acc_overlap_DA.RData")

# load logfold change values CP
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/Data/R_objects/CP_DA_Ancom_ASVLevel/ancomWZ_Bac_Acc.RData")
ancomWZ_Bac_Acc_CP <- ancomWZ_Bac_Acc
rm(ancomWZ_Bac_Acc)

# load logfold change values FP
load(file = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/Data/R_objects/FP_DA_Ancom_ASVLevel/ancomWZ_Bac_Acc.RData")
ancomWZ_Bac_Acc_FP <- ancomWZ_Bac_Acc
rm(ancomWZ_Bac_Acc)
```

## 5.1.2 Extract logfold change
### CP

``` r
# Filter logfold change for DA ASVs per accession
lfc_DA_specific_CP <- lapply(names(ancomWZ_Bac_Acc_CP), function(name) {
  da_asvs <- intersect(Acc_overlap_DA, Acc_Bac_TwoTimes_DA_ASV_CP[[name]]) #Only ASVs that are shared between CP and FP and were DA in CP
  lfc_table <- ancomWZ_Bac_Acc_CP[[name]]$res$lfc
  rownames(lfc_table) <- lfc_table$taxon
  lfc_table <- lfc_table[lfc_table$taxon %in% da_asvs, "Cat_treatmentMb", drop = FALSE]
  lfc_table$ASV <- rownames(lfc_table)
  colnames(lfc_table)[1] <- paste0(name) #, "_CP"
  lfc_table
})

names(lfc_DA_specific_CP) <- names(lfc_DA_specific_CP)

# merge all accessions together
AllAcc_lfc_CP <- Reduce(function(x, y) merge(x, y, by = "ASV", all = TRUE), lfc_DA_specific_CP)
rownames(AllAcc_lfc_CP) <- AllAcc_lfc_CP$ASV
AllAcc_lfc_CP
```

```
##                 ASV        OH DD HE        KI       VL CD RI KT MC HM IT1 GO1
## bASV_1224 bASV_1224 -1.254837 NA NA        NA       NA NA NA NA NA NA  NA  NA
## bASV_315   bASV_315        NA NA NA        NA 1.748393 NA NA NA NA NA  NA  NA
## bASV_656   bASV_656        NA NA NA -1.719751       NA NA NA NA NA NA  NA  NA
```

### FP

``` r
# Filter logfold change for DA ASVs per accession
lfc_DA_specific_FP <- lapply(names(ancomWZ_Bac_Acc_FP), function(name) {
  da_asvs <- intersect(Acc_overlap_DA, Acc_Bac_TwoTimes_DA_ASV_FP[[name]]) #Only ASVs that are shared between CP and FP and were DA in FP
  lfc_table <- ancomWZ_Bac_Acc_FP[[name]]$res$lfc
  rownames(lfc_table) <- lfc_table$taxon
  lfc_table <- lfc_table[lfc_table$taxon %in% da_asvs, "Soil_conditioningMb", drop = FALSE]
  lfc_table$ASV <- rownames(lfc_table)
  colnames(lfc_table)[1] <- paste0(" ", name)
  lfc_table
})

names(lfc_DA_specific_FP) <- names(lfc_DA_specific_FP)

# merge all accessions together
AllAcc_lfc_FP <- Reduce(function(x, y) merge(x, y, by = "ASV", all = TRUE), lfc_DA_specific_FP)
rownames(AllAcc_lfc_FP) <- AllAcc_lfc_FP$ASV

# Remove RI
AllAcc_lfc_FP <- AllAcc_lfc_FP[ , -4]
AllAcc_lfc_FP
```

```
##                 ASV         VL  CD  HM  GO1
## bASV_1224 bASV_1224  1.9549544  NA  NA   NA
## bASV_315   bASV_315  0.4535477  NA  NA   NA
## bASV_656   bASV_656 -1.3437874  NA  NA   NA
```

### Merge CP and FP

``` r
AllAcc_lfc <- merge(AllAcc_lfc_CP, AllAcc_lfc_FP, by = "ASV")
rownames(AllAcc_lfc) <- AllAcc_lfc$ASV
```

### Adding taxonomy

``` r
# Only keep DA in ps object
ps <- prune_taxa(Acc_overlap_DA, unnormalized_bac_ps)

# Phyloseq object to data frame with relative abundance, meta data and taxonomy. Each row is a ASV-sample combination
df_clean <- psmelt(ps)
head(df_clean[, c(1:9, 95:101)], 12)
```

```
##          OTU Sample Abundance Phase_PSF SampleNr Accession Domestication
## 782 bASV_315   F538       127        FP      538        RI    Cultivated
## 780 bASV_315   F479       108        FP      479        VL          Wild
## 777 bASV_315   F478        91        FP      478        VL          Wild
## 568 bASV_315   F010        90        FP       10        CD    Cultivated
## 771 bASV_315   F476        90        FP      476        VL          Wild
## 795 bASV_315   F526        89        FP      526        RI    Cultivated
## 794 bASV_315   FE98        88        FP      E98        RI    Cultivated
## 586 bASV_315   F020        87        FP       20        CD    Cultivated
## 769 bASV_315   F394        83        FP      394        RI    Cultivated
## 779 bASV_315   F536        81        FP      536        RI    Cultivated
## 796 bASV_315   F457        77        FP      457        VL          Wild
## 775 bASV_315   F396        76        FP      396        RI    Cultivated
##     Cat_treatment Soil_conditioning is.neg     Kingdom            Phylum
## 782            Co                Co  FALSE k__Bacteria p__Pseudomonadota
## 780            Mb                Mb  FALSE k__Bacteria p__Pseudomonadota
## 777            Mb                Mb  FALSE k__Bacteria p__Pseudomonadota
## 568            Mb                Co  FALSE k__Bacteria p__Pseudomonadota
## 771            Mb                Mb  FALSE k__Bacteria p__Pseudomonadota
## 795            Co                Mb  FALSE k__Bacteria p__Pseudomonadota
## 794            Mb                Co  FALSE k__Bacteria p__Pseudomonadota
## 586            Mb                Co  FALSE k__Bacteria p__Pseudomonadota
## 769            Mb                Mb  FALSE k__Bacteria p__Pseudomonadota
## 779            Co                Co  FALSE k__Bacteria p__Pseudomonadota
## 796            Mb                Co  FALSE k__Bacteria p__Pseudomonadota
## 775            Mb                Mb  FALSE k__Bacteria p__Pseudomonadota
##                      Class              Order            Family          Genus
## 782 c__Gammaproteobacteria o__Burkholderiales f__Comamonadaceae g__Ramlibacter
## 780 c__Gammaproteobacteria o__Burkholderiales f__Comamonadaceae g__Ramlibacter
## 777 c__Gammaproteobacteria o__Burkholderiales f__Comamonadaceae g__Ramlibacter
## 568 c__Gammaproteobacteria o__Burkholderiales f__Comamonadaceae g__Ramlibacter
## 771 c__Gammaproteobacteria o__Burkholderiales f__Comamonadaceae g__Ramlibacter
## 795 c__Gammaproteobacteria o__Burkholderiales f__Comamonadaceae g__Ramlibacter
## 794 c__Gammaproteobacteria o__Burkholderiales f__Comamonadaceae g__Ramlibacter
## 586 c__Gammaproteobacteria o__Burkholderiales f__Comamonadaceae g__Ramlibacter
## 769 c__Gammaproteobacteria o__Burkholderiales f__Comamonadaceae g__Ramlibacter
## 779 c__Gammaproteobacteria o__Burkholderiales f__Comamonadaceae g__Ramlibacter
## 796 c__Gammaproteobacteria o__Burkholderiales f__Comamonadaceae g__Ramlibacter
## 775 c__Gammaproteobacteria o__Burkholderiales f__Comamonadaceae g__Ramlibacter
```

``` r
# Extract unique ASV–Class-Phylum mapping
tax <- df_clean[!duplicated(df_clean$OTU), c("OTU", "Genus", "Family", "Order", "Class", "Phylum")]

# Row names must match heat map rows
rownames(tax) <- tax$OTU
tax$OTU <- NULL
head(tax)
```

```
##                       Genus            Family              Order
## bASV_315     g__Ramlibacter f__Comamonadaceae o__Burkholderiales
## bASV_656  g__Incertae_Sedis f__Incertae_Sedis      o__Gaiellales
## bASV_1224 g__Incertae_Sedis f__Incertae_Sedis      o__Gaiellales
##                            Class            Phylum
## bASV_315  c__Gammaproteobacteria p__Pseudomonadota
## bASV_656      c__Thermoleophilia p__Actinomycetota
## bASV_1224     c__Thermoleophilia p__Actinomycetota
```

``` r
# Add genus to matrix
AllAcc_lfc2 <- AllAcc_lfc
tax$ASV <- rownames(tax)
AllAcc_lfc2 <- merge(tax, AllAcc_lfc2, ID = "ASV")

# Add higher taxonomic level if genus is NA or Incertae Sedis
taxon <- AllAcc_lfc2$Genus
taxon[is.na(taxon)] <- AllAcc_lfc2$Family[is.na(taxon)]
taxon[is.na(taxon)] <- AllAcc_lfc2$Order[is.na(taxon)]

taxon[taxon == "g__Incertae_Sedis"] <- AllAcc_lfc2$Family[taxon == "g__Incertae_Sedis"]
taxon[taxon == "f__Incertae_Sedis"] <- AllAcc_lfc2$Order[taxon == "f__Incertae_Sedis"]
rownames(AllAcc_lfc2) <- paste(AllAcc_lfc2$ASV, taxon, sep = "_")

AllAcc_lfc2$Phylum <- sub("p__", "", AllAcc_lfc2$Phylum)
AllAcc_lfc2$Class <- sub("c__", "", AllAcc_lfc2$Class)

# Subset data columns
tax <- AllAcc_lfc2[1:6]
AllAcc_lfc2 <- AllAcc_lfc2[7:22]
head(AllAcc_lfc2)
```

```
##                                OH DD HE        KI       VL CD RI KT MC HM IT1
## bASV_1224_o__Gaiellales -1.254837 NA NA        NA       NA NA NA NA NA NA  NA
## bASV_315_g__Ramlibacter        NA NA NA        NA 1.748393 NA NA NA NA NA  NA
## bASV_656_o__Gaiellales         NA NA NA -1.719751       NA NA NA NA NA NA  NA
##                         GO1         VL  CD  HM  GO1
## bASV_1224_o__Gaiellales  NA  1.9549544  NA  NA   NA
## bASV_315_g__Ramlibacter  NA  0.4535477  NA  NA   NA
## bASV_656_o__Gaiellales   NA -1.3437874  NA  NA   NA
```

``` r
# Make color pallet for class
colors15 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
              "#882255", "#44AA99", "#117733", "#332288", "#AA4499", "#DDCC77") # "#000000", "#999999", 

# Define colors
classes2 <- unique(tax$Class)
class_colors2 <- setNames(colors15[1:length(classes2)], classes2)
class_colors2 <- list(Class = class_colors2)

phylums2 <- unique(tax$Phylum)
phylums_colors2 <- setNames(colors15[1:length(phylums2)], phylums2)
phylums_colors2 <- list(Phylum = phylums_colors2)

colors <- c(phylums_colors2, class_colors2)

# Row labels of genus in italic 
rn <- rownames(AllAcc_lfc2)

## Create expression labels
row_labels <- sapply(rn, function(x) {
  tax_part <- sub("^[^_]+_[^_]+_", "", x) # Extract taxonomy
  if (grepl("^g__", tax_part)) { # grab all names with genus
    prefix <- sub("(.*g__).*", "\\1", x) # keep part that does not need to be in italic
    genus <- sub("^g__", "", tax_part) # extract genus
    bquote(.(prefix) * italic(.(genus))) # combine genus and part not in italic
  } else {
    x # change nothing for higher taxonomic levels
}})

## Convert to expression vector
row_labels <- as.expression(row_labels)
# Turn data frame into a matrix
AllAcc_lfc$ASV <- NULL
AllAcc_lfc <- as.matrix(AllAcc_lfc)
AllAcc_lfc2 <- as.matrix(AllAcc_lfc2)

# replace NAs by 0s
AllAcc_lfc_clu <- AllAcc_lfc
AllAcc_lfc2_clu <- AllAcc_lfc2
AllAcc_lfc_clu[is.na(AllAcc_lfc_clu)] <- 0
AllAcc_lfc2_clu[is.na(AllAcc_lfc2_clu)] <- 0
```

### Heatmap

``` r
max_abs <- max(abs(AllAcc_lfc2), na.rm = TRUE)
pheatmap(AllAcc_lfc2,
          clustering_distance_rows = dist(AllAcc_lfc2_clu),
          clustering_distance_cols = dist(t(AllAcc_lfc2_clu)),
         #cluster_cols = FALSE, 
         #cluster_rows = FALSE, 
          color = colorRampPalette(c("blue", "white", "red"))(50),
          breaks = seq(-max_abs, max_abs, length.out = 51),
          gaps_col = c(12),
          #main = "Differentially abundant taxa (log2FC by Family)",
          fontsize_row = 10,
          fontsize_col = 14,
          na_col = "grey90",
          angle_col = 90, 
          annotation_row = tax["Phylum"],
          annotation_colors = colors,
          labels_row = row_labels)
```

![](CP-FP_05_OverlapDAASVs_Heatmap_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

``` r
max_abs <- max(abs(AllAcc_lfc2), na.rm = TRUE)
p <- pheatmap(AllAcc_lfc2,
          #clustering_distance_rows = dist(AllAcc_lfc2_clu),
          clustering_distance_cols = dist(t(AllAcc_lfc2_clu)),
          cluster_cols = FALSE, 
         #cluster_rows = FALSE, 
          color = colorRampPalette(c("blue", "white", "red"))(50),
          breaks = seq(-max_abs, max_abs, length.out = 51),
          gaps_col = c(12),
          #main = "Differentially abundant taxa (log2FC by Family)",
          fontsize_row = 10,
          fontsize_col = 14,
          na_col = "grey90",
          angle_col = 90, 
          annotation_row = tax["Class"],
          annotation_colors = colors,
          labels_row = row_labels)
print(p)
```

![](CP-FP_05_OverlapDAASVs_Heatmap_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

``` r
files <- c("CP_FP_DAHeatmap_Overlap_ASV.svg", "CP_FP_DAHeatmap_Overlap_ASV.png")

mapply(function(x){
  file_path <- file.path("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs", x)

  if (grepl("svg$", x)) {
    svg(file_path, width = 24/2.54, height = 3.2/2.54)
  } else {
    png(file_path, width = 24, height = 3.2, units = "cm", res = 300)
  }
  grid::grid.draw(p$gtable)
  dev.off()
}, files)
```

```
## CP_FP_DAHeatmap_Overlap_ASV.svg.png CP_FP_DAHeatmap_Overlap_ASV.png.png 
##                                   2                                   2
```

``` r
pheatmap(AllAcc_lfc2,
          clustering_distance_rows = dist(AllAcc_lfc2_clu),
          clustering_distance_cols = dist(t(AllAcc_lfc2_clu)),
         #cluster_cols = FALSE, 
         #cluster_rows = FALSE, 
          color = colorRampPalette(c("blue", "white", "red"))(50),
          breaks = seq(-max_abs, max_abs, length.out = 51),
          #main = "Differentially abundant taxa (log2FC by Family)",
          fontsize_row = 10,
          fontsize_col = 14,
          na_col = "grey90",
          angle_col = 90, 
          annotation_row = tax[5:6],
          annotation_colors = colors,
          labels_row = row_labels)
```

![](CP-FP_05_OverlapDAASVs_Heatmap_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

``` r
pheatmap(AllAcc_lfc2,
          clustering_distance_rows = dist(AllAcc_lfc2_clu),
          #clustering_distance_cols = dist(t(AllAcc_lfc2_clu)),
         cluster_cols = FALSE, 
         #cluster_rows = FALSE, 
          color = colorRampPalette(c("blue", "white", "red"))(50),
          breaks = seq(-max_abs, max_abs, length.out = 51),
          #main = "Differentially abundant taxa (log2FC by Family)",
          fontsize_row = 10,
          fontsize_col = 14,
          na_col = "grey90",
          angle_col = 90, 
          annotation_row = tax[5:6],
          annotation_colors = colors,
          labels_row = row_labels)
```

![](CP-FP_05_OverlapDAASVs_Heatmap_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

``` r
pheatmap(AllAcc_lfc2,
          clustering_distance_rows = dist(AllAcc_lfc2_clu),
          clustering_distance_cols = dist(t(AllAcc_lfc2_clu)),
         #cluster_cols = FALSE, 
         #cluster_rows = FALSE, 
          color = colorRampPalette(c("blue", "white", "red"))(50),
          breaks = seq(-max_abs, max_abs, length.out = 51),
          #main = "Differentially abundant taxa (log2FC by Family)",
          fontsize_row = 10,
          fontsize_col = 14,
          na_col = "grey90",
          angle_col = 90, 
          annotation_row = tax[3:6],
          annotation_colors = colors,
          labels_row = row_labels)
```

![](CP-FP_05_OverlapDAASVs_Heatmap_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

``` r
pheatmap(AllAcc_lfc2,
          clustering_distance_rows = dist(AllAcc_lfc2_clu),
          #clustering_distance_cols = dist(t(AllAcc_lfc2_clu)),
          cluster_cols = FALSE, 
         #cluster_rows = FALSE, 
          color = colorRampPalette(c("blue", "white", "red"))(50),
          breaks = seq(-max_abs, max_abs, length.out = 51),
          gaps_col = c(12),
          #main = "Differentially abundant taxa (log2FC by Family)",
          fontsize_row = 10,
          fontsize_col = 14,
          na_col = "grey90",
          angle_col = 90, 
          annotation_row = tax[5:6],
          annotation_colors = colors,
          labels_row = row_labels)
```

![](CP-FP_05_OverlapDAASVs_Heatmap_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

``` r
pheatmap(AllAcc_lfc2,
          clustering_distance_rows = dist(AllAcc_lfc2_clu),
          #clustering_distance_cols = dist(t(AllAcc_lfc2_clu)),
          cluster_cols = FALSE, 
         #cluster_rows = FALSE, 
          color = colorRampPalette(c("blue", "white", "red"))(50),
          breaks = seq(-max_abs, max_abs, length.out = 51),
          gaps_col = c(12),
          #main = "Differentially abundant taxa (log2FC by Family)",
          fontsize_row = 10,
          fontsize_col = 14,
          na_col = "grey90",
          angle_col = 90, 
          annotation_row = tax[3:6],
          annotation_colors = colors,
          labels_row = row_labels)
```

![](CP-FP_05_OverlapDAASVs_Heatmap_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

``` r
# Clean environment
rm(list = ls())
```
