---
title: "FP_05_ Relative abundance per sample - Heatmaps"
author: "Kris de Kreek"
date: "2026-03-18"
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
- [Excamples of pheatmap package](https://davetang.github.io/muse/pheatmap.html)   
   
   
# 4.0 Preperation
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
library(pheatmap)
packageVersion("pheatmap")
```

```
## [1] '1.0.13'
```

``` r
library(tidyr)
packageVersion("tidyr")
```

```
## [1] '1.3.1'
```

``` r
# library(metamisc) #for phyloseq_sep_variable (and more?)
# packageVersion("metamisc")
# library(zinbwave)
# packageVersion("zinbwave")
# library(UpSetR)
# packageVersion("UpSetR")
# library(tuple)
# packageVersion("tuple")
# library(ggpubr)
# packageVersion("ggpubr")
# library(stringr)
# packageVersion("stringr")
# library(ggVennDiagram)
# packageVersion("ggVennDiagram")
```


# 4.10 Heatmap relative abundance per sample
## 4.10.1 Accessions seperate
### Heatmap At ASV level DA ASVs occuring twice
#### Load data

``` r
# load phyloseq object
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_unnormalized_bac_ps_ForDA.RData")

# load DA ASVs
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/TwelveAccessionExperiment/MicrobiomeAnalysis/Data/Phyloseq_objects/FP_DA_SummaryFiles/SummaryFiles_ASVLevel/Acc_Bac_TwoTimes_DA_ASV.RData")
```

#### Remove RI and batch 1

``` r
FP_unnormalized_bac_ps_ForDA <- subset_samples(FP_unnormalized_bac_ps_ForDA, Accession != "RI")
FP_unnormalized_bac_ps_ForDA <- subset_samples(FP_unnormalized_bac_ps_ForDA, Batch != 1)
```

#### Prepare data

``` r
# Put all unique DA ASVs in one vector
Acc_Bac_TwoTimes_DA_ASV$RI <- NULL
Acc_Bac_TwoTimes_DA_ASV_comb <- unique(unlist(Acc_Bac_TwoTimes_DA_ASV, use.names = FALSE))

# Calculate relative abundance
FP_ps <- transform_sample_counts(FP_unnormalized_bac_ps_ForDA, function(x) x / sum(x))

# Only keep DA in ps object
FP_ps <- prune_taxa(Acc_Bac_TwoTimes_DA_ASV_comb, FP_ps)

# Phyloseq object to data frame with relative abundance, meta data and taxonomy. Each row is a ASV-sample combination
df <- psmelt(FP_ps)
head(df[, c(1:9)], 12)
```

```
##          OTU Sample  Abundance Phase_PSF SampleNr Accession Domestication
## 16920 bASV_7   F475 0.08862542        FP      475        VL          Wild
## 16892 bASV_7   F478 0.07961704        FP      478        VL          Wild
## 16891 bASV_7   F479 0.07325433        FP      479        VL          Wild
## 16913 bASV_7   F176 0.07283707        FP      176        HM    Cultivated
## 16883 bASV_7   F175 0.06932283        FP      175        HM    Cultivated
## 16903 bASV_7   F476 0.06616079        FP      476        VL          Wild
## 16902 bASV_7   F477 0.05984501        FP      477        VL          Wild
## 16896 bASV_7   F116 0.05406350        FP      116       GO1    Cultivated
## 16925 bASV_7   F177 0.05290581        FP      177        HM    Cultivated
## 16909 bASV_7   F095 0.05268142        FP       95       GO1    Cultivated
## 16921 bASV_7   F117 0.05239759        FP      117       GO1    Cultivated
## 16884 bASV_7   F179 0.05215760        FP      179        HM    Cultivated
##       Cat_treatment Soil_conditioning
## 16920            Mb                Mb
## 16892            Mb                Mb
## 16891            Mb                Mb
## 16913            Mb                Co
## 16883            Mb                Co
## 16903            Mb                Mb
## 16902            Mb                Mb
## 16896            Mb                Mb
## 16925            Mb                Co
## 16909            Mb                Co
## 16921            Mb                Mb
## 16884            Mb                Co
```

``` r
# unique accessions
Acc <- names(Acc_Bac_TwoTimes_DA_ASV)
```

#### Matrix per accession

``` r
# Subset data per accession and DA ASVs
Rel_Ab_l <- list()
for (i in seq_along(Acc)) {
  Acc_sub <- Acc[i]
  sub_df <- df[df$Accession == Acc_sub & df$OTU %in% Acc_Bac_TwoTimes_DA_ASV[[Acc_sub]], c(1:3, 9, 96:100)]
  sub_df <- sub_df[order(sub_df$Soil_conditioning, sub_df$Sample), ]
  sub_df$Sample <- paste(sub_df$Sample, sub_df$Soil_conditioning, sep = "_")
    #colnames(sub_df)[3] <- Acc_sub 
  Rel_Ab_l[[i]] <- sub_df
}

names(Rel_Ab_l) <- Acc

# Number of Co Samples
Numb_Co <- lapply(Rel_Ab_l, function(x) length(x[x$Soil_conditioning == "Co", "Soil_conditioning"])/length(unique(x$OTU)))

# Put each sample in seperate column
Rel_Ab_l_wide <- lapply(Rel_Ab_l, function(df) {
  df$Soil_conditioning <- NULL
  pivot_wider(df, names_from  = Sample, values_from = Abundance, values_fill = 0)
})

# Make taxonomy rownames
Rel_Ab_l_wide <- lapply(Rel_Ab_l_wide, function(df) {
  df <- as.data.frame(df)
  taxon <- df$Genus
  taxon[is.na(taxon)] <- df$Family[is.na(taxon)]
  taxon[is.na(taxon)] <- df$Order[is.na(taxon)]
  taxon[is.na(taxon)] <- df$Class[is.na(taxon)]
  taxon[is.na(taxon)] <- df$Phylum[is.na(taxon)]
  
  taxon[taxon == "g__Incertae_Sedis"] <- df$Family[taxon == "g__Incertae_Sedis"]
  taxon[taxon == "f__Incertae_Sedis"] <- df$Order[taxon == "f__Incertae_Sedis"]
  taxon[taxon == "o__Incertae_Sedis"] <- df$Class[taxon == "o__Incertae_Sedis"]
  taxon[taxon == "c__Incertae_Sedis"] <- df$Phylum[taxon == "c__Incertae_Sedis"]
  
  taxon[taxon == "g__Unclassified"] <- df$Family[taxon == "g__Unclassified"]
  taxon[taxon == "f__Unclassified"] <- df$Order[taxon == "f__Unclassified"]
  taxon[taxon == "o__Unclassified"] <- df$Class[taxon == "o__Unclassified"]
  taxon[taxon == "c__Unclassified"] <- df$Phylum[taxon == "c__Unclassified"]
  
  rownames(df) <- paste(df$OTU, taxon, sep = "_")
  #df$taxon <- make.unique(paste(taxon, sep = "_"))
  
  df$Phylum <- sub("p__", "", df$Phylum)
  df$Class <- sub("c__", "", df$Class)
  return(df)
})

# Extract taxonomy df
Rel_Ab_l_wide_tax <- lapply(Rel_Ab_l_wide, function(df){
  df[, 1:6]
})

Rel_Ab_l_wide2 <- Rel_Ab_l_wide

# Extract abundance
Rel_Ab_l_wide_m <- lapply(Rel_Ab_l_wide, function(df) {
  df <- df[, 7:ncol(df)]
  return(df)
})
```

#### Define colors class annotation

``` r
# Make color pallet for class
colors35 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#332288", "#CC79A7", 
              "#D55E00", "#0072B2", "#000000", "#999933",
              "#117733", "#88CCEE", "#44AA99", "#DDCC77", "#CC6677", "#882255", "#AA4499",
              "#661100", "#6699CC", "#AA4466", "#4477AA", "#EE6677", "#BBBBBB", 
              "#2E3440", "#8FBCBB", "#BF616A", "#A3BE8C", "#EBCB8B", "#B48EAD", 
              "#5E81AC", "#D08770", "#3B4252", "#7B9E87", "#999999")

phylum_class_colors <- lapply(names(Rel_Ab_l_wide_tax), function(x) {
  classes <- unique(Rel_Ab_l_wide_tax[[x]]$Class)
  class_colors <- setNames(colors35[1:length(classes)], classes)
  class_colors <- list(Class = class_colors)
  
  phylums <- unique(Rel_Ab_l_wide_tax[[x]]$Phylum)
  phylums_colors <- setNames(colors35[1:length(phylums)], phylums)
  phylums_colors <- list(Phylum = phylums_colors)
  
  phylum_class_colors <- c(phylums_colors, class_colors)
  return(phylum_class_colors)
})

names(phylum_class_colors) <- names(Rel_Ab_l_wide_tax)

# Turn data frame into a matrix
Rel_Ab_l_wide_m <- lapply(Rel_Ab_l_wide_m, function(df) {
  as.matrix(df)
})

# replace NAs by 0s
Rel_Ab_l_wide_m_clu <- lapply(Rel_Ab_l_wide_m, function(df) {
  df[is.na(df)] <- 0
  return(df)
})
```


#### Print plots VL

``` r
# breaks in headmap
max_abs <- max(abs(Rel_Ab_l_wide_m$VL), na.rm = TRUE)
breaks <- unique(c(seq(0, 0.001*max_abs, length.out = 40), # More breaks in this range for sensitivity
                   seq(0.001*max_abs, 0.01*max_abs, length.out = 40),
                   seq(0.01*max_abs, 0.1*max_abs, length.out = 200),
                   seq(0.1*max_abs, max_abs, length.out = 15)))
#breaks <- c(0, 10^seq(log10(max_abs * 1e-4), log10(max_abs), length.out = 100))

# Headmap
pheatmap(Rel_Ab_l_wide_m$VL,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        main = "VL",
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$VL),
        na_col = "grey90",
        angle_col = 90,
        #legend_breaks = c(0.001*round(max_abs, 2), 0.01*round(max_abs, 2), round(max_abs, 2))
        )
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

``` r
# Headmap with class annotation
pheatmap(Rel_Ab_l_wide_m$VL,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        main = "VL",
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$VL),
        na_col = "grey90",
        angle_col = 90,
        #legend_breaks = c(0.001*round(max_abs, 2), 0.01*round(max_abs, 2), round(max_abs, 2))
        annotation_row = Rel_Ab_l_wide_tax$VL["Class"],
        annotation_colors = phylum_class_colors$VL)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

``` r
# Headmap without labels
pheatmap(Rel_Ab_l_wide_m$VL,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        #main = "VL",
        show_colnames = FALSE,
        show_rownames = FALSE,
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$VL),
        na_col = "grey90",
        angle_col = 90,
        annotation_row = Rel_Ab_l_wide_tax$VL["Class"],
        annotation_colors = phylum_class_colors$VL)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

#### Print plots CD

``` r
# breaks in headmap
max_abs <- max(abs(Rel_Ab_l_wide_m$CD), na.rm = TRUE)
breaks <- unique(c(seq(0, 0.5*max_abs, length.out = 30), # More breaks in this range for sensitivity
                   seq(0.5*max_abs, max_abs, length.out = 20)))

# Headmap
pheatmap(Rel_Ab_l_wide_m$CD,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        main = "CD",
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$CD),
        na_col = "grey90",
        angle_col = 90)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

``` r
# Headmap
pheatmap(Rel_Ab_l_wide_m$CD,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        main = "CD",
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$CD),
        na_col = "grey90",
        angle_col = 90,
        annotation_row = Rel_Ab_l_wide_tax$CD["Class"],
        annotation_colors = phylum_class_colors$CD)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

``` r
pheatmap(Rel_Ab_l_wide_m$CD,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        #main = "CD",
        show_colnames = FALSE,
        show_rownames = FALSE,
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$CD),
        na_col = "grey90",
        angle_col = 90,
        annotation_row = Rel_Ab_l_wide_tax$CD["Class"],
        annotation_colors = phylum_class_colors$CD)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

#### Print plots HM

``` r
# breaks in headmap
max_abs <- max(abs(Rel_Ab_l_wide_m$HM), na.rm = TRUE)
breaks <- unique(c(seq(0, 0.01*max_abs, length.out = 80), # More breaks in this range for sensitivity
                   seq(0.01*max_abs, 0.1*max_abs, length.out = 200),
                   seq(0.1*max_abs, max_abs, length.out = 15)))

# Headmap
pheatmap(Rel_Ab_l_wide_m$HM,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        main = "HM",
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$HM),
        na_col = "grey90",
        angle_col = 90)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

``` r
# Headmap
pheatmap(Rel_Ab_l_wide_m$HM,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        main = "HM",
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$HM),
        na_col = "grey90",
        angle_col = 90,
        annotation_row = Rel_Ab_l_wide_tax$HM["Class"],
        annotation_colors = phylum_class_colors$HM)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

``` r
# Headmap
pheatmap(Rel_Ab_l_wide_m$HM,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        #main = "HM",
        show_colnames = FALSE,
        show_rownames = FALSE,
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$HM),
        na_col = "grey90",
        angle_col = 90,
        annotation_row = Rel_Ab_l_wide_tax$HM["Class"],
        annotation_colors = phylum_class_colors$HM)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

#### Print plots GO1

``` r
# breaks in headmap
max_abs <- max(abs(Rel_Ab_l_wide_m$GO1), na.rm = TRUE)
breaks <- unique(c(seq(0, 0.02*max_abs, length.out = 120), # More breaks in this range for sensitivity
                   seq(0.02*max_abs, 0.5*max_abs, length.out = 300),
                   seq(0.5*max_abs, max_abs, length.out = 30)))

# Headmap
pheatmap(Rel_Ab_l_wide_m$GO1,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        main = "GO1",
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$GO1),
        na_col = "grey90",
        angle_col = 90)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

``` r
# Headmap
pheatmap(Rel_Ab_l_wide_m$GO1,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        main = "GO1",
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$GO1),
        na_col = "grey90",
        angle_col = 90,
        annotation_row = Rel_Ab_l_wide_tax$GO1["Class"],
        annotation_colors = phylum_class_colors$GO1)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

``` r
# Headmap
pheatmap(Rel_Ab_l_wide_m$GO1,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        #main = "GO1",
        show_colnames = FALSE,
        show_rownames = FALSE,
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$GO1),
        na_col = "grey90",
        angle_col = 90,
        annotation_row = Rel_Ab_l_wide_tax$GO1["Class"],
        annotation_colors = phylum_class_colors$GO1)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-18-1.png)<!-- -->


### Heatmaps with top 10 classes

``` r
Rel_Ab_l_wide3 <- lapply(Rel_Ab_l_wide2, function(df) {
  # Order by Class
  tab_cla <- as.data.frame(table(df$Class))
  tab_cla <- tab_cla[order(tab_cla$Freq, decreasing = TRUE), ]
  
  # select top 7 classes
  top10_cla <- as.character(tab_cla$Var1[1:10])
  
  # replace other classes "Other"
  df2 <- df
  df2$Class <- ifelse(df2$Class %in% top10_cla,
                         df2$Class,
                         "Other")
  
  # put other last for plotting
  if ("Other" %in% df2$Class) {
  df2$Class <- factor(df2$Class, levels = c(top10_cla, "Other"))
} else {
  df2$Class <- factor(df2$Class, levels = top10_cla)
}
  return(df2)
})

# Extract taxonomy df
Rel_Ab_l_wide_tax3 <- lapply(Rel_Ab_l_wide3, function(df){
  df[, 1:6]
})

# Extract abundance
Rel_Ab_l_wide_m3 <- lapply(Rel_Ab_l_wide3, function(df) {
  df <- df[, 7:ncol(df)]
  return(df)
})
```

#### Define colors class annotation

``` r
phylum_class_colors3 <- lapply(names(Rel_Ab_l_wide_tax3), function(x) {
  classes <- levels(Rel_Ab_l_wide_tax3[[x]]$Class)
  
  # Check if "Other" is present
  has_other <- "Other" %in% classes
  
  if (has_other) {
    class_colors <- setNames(colors35[1:(length(classes)-1)], classes[classes != "Other"])
    class_colors <- c(class_colors, Other = "#BBBBBB")
  } else {
    class_colors <- setNames(colors35[1:length(classes)], classes)
  }
  
  class_colors <- list(Class = class_colors)
  
  # Phylum colors (unchanged)
  phylums <- unique(Rel_Ab_l_wide_tax3[[x]]$Phylum)
  phylums_colors <- setNames(colors35[1:length(phylums)], phylums)
  phylums_colors <- list(Phylum = phylums_colors)
  
  phylum_class_colors <- c(phylums_colors, class_colors)
  return(phylum_class_colors)
})

names(phylum_class_colors3) <- names(Rel_Ab_l_wide_tax3)

# Turn data frame into a matrix
Rel_Ab_l_wide_m3 <- lapply(Rel_Ab_l_wide_m3, function(df) {
  as.matrix(df)
})

# replace NAs by 0s
Rel_Ab_l_wide_m3_clu <- lapply(Rel_Ab_l_wide_m3, function(df) {
  df[is.na(df)] <- 0
  return(df)
})
```

#### Print plots VL

``` r
# breaks in headmap
max_abs <- max(abs(Rel_Ab_l_wide_m3$VL), na.rm = TRUE)
breaks <- unique(c(seq(0, 0.001*max_abs, length.out = 40), # More breaks in this range for sensitivity
                   seq(0.001*max_abs, 0.01*max_abs, length.out = 40),
                   seq(0.01*max_abs, 0.1*max_abs, length.out = 200),
                   seq(0.1*max_abs, max_abs, length.out = 15)))
```

``` r
# Headmap with class annotation
pheatmap(Rel_Ab_l_wide_m3$VL,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        main = "VL",
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$VL),
        na_col = "grey90",
        angle_col = 90,
        #legend_breaks = c(0.001*round(max_abs, 2), 0.01*round(max_abs, 2), round(max_abs, 2))
        annotation_row = Rel_Ab_l_wide_tax3$VL["Class"],
        annotation_colors = phylum_class_colors3$VL)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

``` r
# Headmap without labels
p <- pheatmap(Rel_Ab_l_wide_m3$VL,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        #main = "VL",
        show_colnames = FALSE,
        show_rownames = FALSE,
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$VL),
        na_col = "grey90",
        angle_col = 90,
        annotation_row = Rel_Ab_l_wide_tax3$VL["Class"],
        annotation_colors = phylum_class_colors3$VL)
print(p)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

``` r
files <- c("FP_DAHeatmap_ASV_persample_VL.svg", "FP_DAHeatmap_ASV_persample_VL.png")

mapply(function(x){
  file_path <- file.path("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs", x)

  if (grepl("svg$", x)) {
    svg(file_path, width = 17/2.54, height = 18/2.54)
  } else {
    png(file_path, width = 17, height = 18, units = "cm", res = 300)
  }
  grid::grid.draw(p$gtable)
  dev.off()
}, files)
```

```
## FP_DAHeatmap_ASV_persample_VL.svg.png FP_DAHeatmap_ASV_persample_VL.png.png 
##                                     2                                     2
```

#### Print plots CD

``` r
# breaks in headmap
max_abs <- max(abs(Rel_Ab_l_wide_m3$CD), na.rm = TRUE)
breaks <- unique(c(seq(0, 0.5*max_abs, length.out = 30), # More breaks in this range for sensitivity
                   seq(0.5*max_abs, max_abs, length.out = 20)))
```

``` r
# Headmap
pheatmap(Rel_Ab_l_wide_m3$CD,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        main = "CD",
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$CD),
        na_col = "grey90",
        angle_col = 90,
        annotation_row = Rel_Ab_l_wide_tax3$CD["Class"],
        annotation_colors = phylum_class_colors3$CD)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

``` r
q <- pheatmap(Rel_Ab_l_wide_m3$CD,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        #main = "CD",
        show_colnames = FALSE,
        show_rownames = FALSE,
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$CD),
        na_col = "grey90",
        angle_col = 90,
        annotation_row = Rel_Ab_l_wide_tax3$CD["Class"],
        annotation_colors = phylum_class_colors3$CD)
print(q)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-26-1.png)<!-- -->

``` r
files <- c("FP_DAHeatmap_ASV_persample_CD.svg", "FP_DAHeatmap_ASV_persample_CD.png")

mapply(function(x){
  file_path <- file.path("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs", x)

  if (grepl("svg$", x)) {
    svg(file_path, width = 18/2.54, height = 7/2.54)
  } else {
    png(file_path, width = 18, height = 7, units = "cm", res = 300)
  }
  grid::grid.draw(q$gtable)
  dev.off()
}, files)
```

```
## FP_DAHeatmap_ASV_persample_CD.svg.png FP_DAHeatmap_ASV_persample_CD.png.png 
##                                     2                                     2
```

#### Print plots HM

``` r
# breaks in headmap
max_abs <- max(abs(Rel_Ab_l_wide_m3$HM), na.rm = TRUE)
breaks <- unique(c(seq(0, 0.01*max_abs, length.out = 80), # More breaks in this range for sensitivity
                   seq(0.01*max_abs, 0.1*max_abs, length.out = 200),
                   seq(0.1*max_abs, max_abs, length.out = 15)))
```

``` r
# Headmap
pheatmap(Rel_Ab_l_wide_m3$HM,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        main = "HM",
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$HM),
        na_col = "grey90",
        angle_col = 90,
        annotation_row = Rel_Ab_l_wide_tax3$HM["Class"],
        annotation_colors = phylum_class_colors3$HM)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

``` r
# Headmap
r <- pheatmap(Rel_Ab_l_wide_m3$HM,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        #main = "HM",
        show_colnames = FALSE,
        show_rownames = FALSE,
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$HM),
        na_col = "grey90",
        angle_col = 90,
        annotation_row = Rel_Ab_l_wide_tax3$HM["Class"],
        annotation_colors = phylum_class_colors3$HM)
print(r)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-29-1.png)<!-- -->

``` r
files <- c("FP_DAHeatmap_ASV_persample_HM.svg", "FP_DAHeatmap_ASV_persample_HM.png")

mapply(function(x){
  file_path <- file.path("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs", x)

  if (grepl("svg$", x)) {
    svg(file_path, width = 16/2.54, height = 18/2.54)
  } else {
    png(file_path, width = 16, height = 18, units = "cm", res = 300)
  }
  grid::grid.draw(r$gtable)
  dev.off()
}, files)
```

```
## FP_DAHeatmap_ASV_persample_HM.svg.png FP_DAHeatmap_ASV_persample_HM.png.png 
##                                     2                                     2
```

#### Print plots GO1

``` r
# breaks in headmap
max_abs <- max(abs(Rel_Ab_l_wide_m3$GO1), na.rm = TRUE)
breaks <- unique(c(seq(0, 0.02*max_abs, length.out = 120), # More breaks in this range for sensitivity
                   seq(0.02*max_abs, 0.5*max_abs, length.out = 300),
                   seq(0.5*max_abs, max_abs, length.out = 30)))
```

``` r
# Headmap
pheatmap(Rel_Ab_l_wide_m3$GO1,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        main = "GO1",
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$GO1),
        na_col = "grey90",
        angle_col = 90,
        annotation_row = Rel_Ab_l_wide_tax3$GO1["Class"],
        annotation_colors = phylum_class_colors3$GO1)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-31-1.png)<!-- -->

``` r
# Headmap
s <- pheatmap(Rel_Ab_l_wide_m3$GO1,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red", "darkred", "black"))(length(breaks) - 1),
        breaks = breaks,
        cellwidth = 20,
        #main = "GO1",
        show_colnames = FALSE,
        show_rownames = FALSE,
        fontsize_row = 7,
        fontsize_col = 14,
        gaps_col = c(Numb_Co$GO1),
        na_col = "grey90",
        angle_col = 90,
        annotation_row = Rel_Ab_l_wide_tax3$GO1["Class"],
        annotation_colors = phylum_class_colors3$GO1)

print(s)
```

![](FP_05_Heatmaps_RelativeAbund_PerSample_files/figure-html/unnamed-chunk-32-1.png)<!-- -->

``` r
files <- c("FP_DAHeatmap_ASV_persample_GO1.svg", "FP_DAHeatmap_ASV_persample_GO1.png")

mapply(function(x){
  file_path <- file.path("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs", x)

  if (grepl("svg$", x)) {
    svg(file_path, width = 18/2.54, height = 7/2.54)
  } else {
    png(file_path, width = 18, height = 7, units = "cm", res = 300)
  }
  grid::grid.draw(s$gtable)
  dev.off()
}, files)
```

```
## FP_DAHeatmap_ASV_persample_GO1.svg.png FP_DAHeatmap_ASV_persample_GO1.png.png 
##                                      2                                      2
```
