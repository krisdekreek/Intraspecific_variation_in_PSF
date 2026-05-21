---
title: "FP_02_beta_diversity - Outliers removed - Only Batch 2 without RI"
author: "Kris de Kreek"
date: "2026-05-21"
output: 
  html_document:
    toc: true
    keep_md: true
---

# 2.0 - Basic beta diversity analysis
This script should be run after script "01_loading_and_pre_processing.rmd" because it generates all the global environment objects and also loads the libraries

On this script we will evaluate beta diversity - the differences in the microbial community compositions across samples. We will make ordination plots and multivariate tests

This script is originally made by Pedro Beschoren da Costa for the [MeJA_Pilot](https://github.com/PedroBeschoren/MeJA_Pilot) and modified to fit my data.

### Load libraries

``` r
library(purrr)
packageVersion("purrr")
```

```
## [1] '1.1.0'
```

``` r
library(phyloseq)
packageVersion("phyloseq")
```

```
## [1] '1.52.0'
```

``` r
library(ggplot2)
packageVersion("ggplot2")
```

```
## [1] '4.0.0'
```

``` r
library(vegan)
packageVersion("vegan")
```

```
## [1] '2.7.2'
```

``` r
library(tibble) #needed for function rownames_to_columns
packageVersion("tibble")
```

```
## [1] '3.3.0'
```

``` r
library(factoextra) #for heatmap function fviz_dist
packageVersion("factoextra")
```

```
## [1] '1.0.7'
```

``` r
library(patchwork) # for wrap_plot
packageVersion("patchwork")
```

```
## [1] '1.3.2'
```

``` r
library(cowplot) # for extracting legend
packageVersion("cowplot")
```

```
## [1] '1.2.0'
```

``` r
library(ggh4x) #for coloring strip facet wrap
packageVersion("ggh4x")
```

```
## [1] '0.3.1'
```


# 2.1 - Beta Diversity plots
Beta diversity plots are the beating heart or microbiome analysis. Here you will be able to visually tell if communities differ according to treatment or not. It can be a very long topic but here I only use one option.

## 2.1.1 - Loading data

``` r
load("C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/Data/R_objects/FP_CSS_bac_ps.RData")
```

### Remove outliers

``` r
Outliers <- c("F178")
samples_to_keep <- setdiff(sample_names(FP_CSS_bac_ps), Outliers) # only keep samples that are not in both objects
FP_CSS_bac_ps <- prune_samples(samples_to_keep, FP_CSS_bac_ps) #only keep samples that are in samples_to_keep
```

### Remove batch 1 and RI

``` r
FP_CSS_bac_ps <- subset_samples(FP_CSS_bac_ps, Batch != 1)
FP_CSS_bac_ps <- subset_samples(FP_CSS_bac_ps, Accession != "RI")
levels(sample_data(FP_CSS_bac_ps)$Accession)
```

```
## [1] "VL"  "CD"  "HM"  "GO1"
```

### Remove unplanted CTRL
Unplanted CTRL samples appears to from the CP.

``` r
FP_CSS_bac_ps <- subset_samples(FP_CSS_bac_ps, Accession != "Un")
```

### Remove Co cat treatment

``` r
FP_CSS_bac_ps <- subset_samples(FP_CSS_bac_ps, Cat_treatment != "Co")
```

### Phyloseq object per accession

``` r
Acc <- c("VL", "CD", "HM", "GO1")
filtered_ps <- list()

filtered_ps <- lapply(Acc, function(x) {
  samples_to_keep <- sample_names(FP_CSS_bac_ps)[sample_data(FP_CSS_bac_ps)$Accession == x]
  prune_samples(samples_to_keep, FP_CSS_bac_ps)
})
names(filtered_ps) <- Acc  # Name each list element with the corresponding accession
```


## 2.1.3 PCoA plots at ASV level
### PCoA plot without unplanted CTRL

``` r
# Calculate Bray-Curtis distance matrix
bray_dist <- distance(FP_CSS_bac_ps, method = "bray")

# Run classical MDS with cmdscale()
mds_result <- cmdscale(bray_dist, k = 2, eig = TRUE)

# Extract coordinates
mds_coords <- as.data.frame(mds_result$points)
colnames(mds_coords) <- c("MDS1", "MDS2")

# Add sample IDs as a column
mds_coords$SampleID <- rownames(mds_coords)

# Extract sample metadata
metadata_df <- data.frame(sample_data(FP_CSS_bac_ps))
metadata_df$SampleID <- rownames(metadata_df)

# Merge ordination coordinates with metadata
plot_df <- merge(mds_coords, metadata_df, by = "SampleID")

# Eigenvalues percentage explained variance
perc_expl1 <- 100 * (mds_result$eig / sum(mds_result$eig))
perc_expl1[1:4] # percentage explained first four axes
```

```
## [1] 8.712538 5.257216 4.758322 3.795875
```

``` r
barplot(perc_expl1[1:10], names = paste ('PCoA', 1:10), las = 3, ylab = 'eigenvalues %')
```

![](FP_02_beta_diversity_OutliersRemoved_Batch2_final_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

#### Accessions colored 

``` r
colors4 <- c("#0072B2", "#117733", "#D55E00", "#999999") #"#009E73", "#CC79A7", "#E69F00", "#000000", "#F0E442", "#56B4E9", "#117733"

p <- ggplot(plot_df, aes(x = MDS1, y = MDS2, color = Accession, shape = Soil_conditioning)) +
  geom_point(size = 1.5, alpha = 1) +
  stat_ellipse(aes(group = Accession),
               show.legend = FALSE) +
  scale_color_manual(values = colors4) + 
  scale_shape_manual(name = "Soil treatment", 
                     values = c(Co = 16, Mb = 15),
                     labels = c(Co = "Control", Mb = expression(italic("M. brassicae")))) +
  xlab(paste0("PCoA1 (", round(perc_expl1[1], 1), " %)")) +
  ylab(paste0("PCoA2 (", round(perc_expl1[2], 1), " %)")) +
  ggtitle("PERMANOVA Accession: p = 0.005;\nSoil: p = 0.005, Interaction: p = 0.005") +
  theme(panel.background = element_rect(fill = "white"), #theme of the graph
        panel.border = element_blank(), 
        axis.line = element_line(), #show line for x and y axis
        plot.title = element_text(size = 10, hjust = 0),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 14),
        axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_text(vjust = 2.5), #, hjust = -0
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        #strip.text = element_text(size = 13),
        plot.margin = margin(5.5, 12, 5.5, 5.5, "points"),
        legend.position = "none")
print(p)
```

```
## Warning: The following aesthetics were dropped during statistical transformation: shape.
## ℹ This can happen when ggplot fails to infer the correct grouping structure in
##   the data.
## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical
##   variable into a factor?
```

![](FP_02_beta_diversity_OutliersRemoved_Batch2_final_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

``` r
mapply(function(x)
  ggsave(plot = last_plot(), 
         filename = x,
         path = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs",
         scale = 2.2,
         width = 4.3,
         height = 4.45,
         units = "cm",
         dpi = 300),
  x = c("FP_BetaDiv_Acc_Soil.svg", "FP_BetaDiv_Acc_Soil.png"))
```

```
## Warning: The following aesthetics were dropped during statistical transformation: shape.
## ℹ This can happen when ggplot fails to infer the correct grouping structure in
##   the data.
## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical
##   variable into a factor?
## The following aesthetics were dropped during statistical transformation: shape.
## ℹ This can happen when ggplot fails to infer the correct grouping structure in
##   the data.
## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical
##   variable into a factor?
```

```
##                                                                                                                     FP_BetaDiv_Acc_Soil.svg 
## "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs/FP_BetaDiv_Acc_Soil.svg" 
##                                                                                                                     FP_BetaDiv_Acc_Soil.png 
## "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs/FP_BetaDiv_Acc_Soil.png"
```

``` r
# Legend
legend_acc <- get_legend(p + theme(legend.position = "right"))
```

```
## Warning: The following aesthetics were dropped during statistical transformation: shape.
## ℹ This can happen when ggplot fails to infer the correct grouping structure in
##   the data.
## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical
##   variable into a factor?
```

``` r
ggdraw(legend_acc)
```

![](FP_02_beta_diversity_OutliersRemoved_Batch2_final_files/figure-html/unnamed-chunk-9-2.png)<!-- -->

``` r
mapply(function(x)
  ggsave(plot = last_plot(), 
         filename = x,
         path = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs",
         scale = 2.2,
         width = 1.8,
         height = 2.9,
         units = "cm",
         dpi = 300),
  x = c("FP_BetaDiv_Acc_Soil_legend.svg", "FP_BetaDiv_Acc_Soil_legend.png"))
```

```
##                                                                                                                     FP_BetaDiv_Acc_Soil_legend.svg 
## "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs/FP_BetaDiv_Acc_Soil_legend.svg" 
##                                                                                                                     FP_BetaDiv_Acc_Soil_legend.png 
## "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs/FP_BetaDiv_Acc_Soil_legend.png"
```

#### Accession colored without elipse

``` r
ggplot(plot_df, aes(x = MDS1, y = MDS2, color = Accession, shape = Soil_conditioning)) +
  geom_point(size = 1.5, alpha = 1) +
  xlab(paste0("PCoA1 (", round(perc_expl1[1], 1), " %)")) +
  ylab(paste0("PCoA2 (", round(perc_expl1[2], 1), " %)")) +
  ggtitle("PERMANOVA Accession: p = 0.001") +
  theme(panel.background = element_rect(fill = "white"), #theme of the graph
        panel.border = element_blank(),
        axis.line = element_line(), #show line for x and y axis
        plot.title = element_text(size = 10, hjust = 0),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_text(vjust = 2.5), #, hjust = -0
        strip.text = element_text(size = 13),
        plot.margin = margin(5.5, 12, 5.5, 5.5, "points"),
        legend.position = "right")
```

![](FP_02_beta_diversity_OutliersRemoved_Batch2_final_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

#### Accession colored with elipse

``` r
ggplot(plot_df, aes(x = MDS1, y = MDS2, color = Accession)) +
  geom_point(size = 1.5, alpha = 1) +
  stat_ellipse() +
  xlab(paste0("PCoA1 (", round(perc_expl1[1], 1), " %)")) +
  ylab(paste0("PCoA2 (", round(perc_expl1[2], 1), " %)")) +
  ggtitle("PERMANOVA Accession: p = 0.001") +
  theme(panel.background = element_rect(fill = "white"), #theme of the graph
        panel.border = element_blank(),
        axis.line = element_line(), #show line for x and y axis
        plot.title = element_text(size = 10, hjust = 0),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_text(vjust = 2.5), #, hjust = -0
        strip.text = element_text(size = 13),
        plot.margin = margin(5.5, 12, 5.5, 5.5, "points"),
        legend.position = "right")
```

![](FP_02_beta_diversity_OutliersRemoved_Batch2_final_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

#### Soil conditioning colored

``` r
q <- ggplot(plot_df, aes(x = MDS1, y = MDS2, color = Soil_conditioning, shape = Soil_conditioning)) +
  geom_point(size = 1.5, alpha = 1) +
  stat_ellipse() +
  scale_color_manual(name = "Soil treatment", 
                     values = c("#882255", "#009E73"),
                     labels = c(Co = "Control", Mb = expression(italic("M. brassicae")))) +
  scale_shape_manual(name = "Soil treatment", 
                     values = c(Co = 16, Mb = 15),
                     labels = c(Co = "Control", Mb = expression(italic("M. brassicae")))) +
  xlab(paste0("PCoA1 (", round(perc_expl1[1], 1), " %)")) +
  ylab(paste0("PCoA2 (", round(perc_expl1[2], 1), " %)")) +
  ggtitle("PERMANOVA Soil conditioning: p = 0.005") +
  theme(panel.background = element_rect(fill = "white"), #theme of the graph
        panel.border = element_blank(),
        axis.line = element_line(), #show line for x and y axis
        plot.title = element_text(size = 10, hjust = 0),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_text(vjust = 2.5), #, hjust = -0
        strip.text = element_text(size = 13),
        plot.margin = margin(5.5, 12, 5.5, 5.5, "points"),
        legend.position = "right")

# Legend
legend_acc <- get_legend(q + theme(legend.position = "right"))
ggdraw(legend_acc)
```

![](FP_02_beta_diversity_OutliersRemoved_Batch2_final_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

``` r
mapply(function(x)
  ggsave(plot = last_plot(), 
         filename = x,
         path = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs",
         scale = 2.2,
         width = 1.4,
         height = 1,
         units = "cm",
         dpi = 300),
  x = c("FP_BetaDiv_Soil_legend.svg", "FP_BetaDiv_Soil_legend.png"))
```

```
##                                                                                                                     FP_BetaDiv_Soil_legend.svg 
## "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs/FP_BetaDiv_Soil_legend.svg" 
##                                                                                                                     FP_BetaDiv_Soil_legend.png 
## "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs/FP_BetaDiv_Soil_legend.png"
```

#### Soil conditioning coloured with elipse

``` r
ggplot(plot_df, aes(x = MDS1, y = MDS2, color = Soil_conditioning)) +
  geom_point(size = 1.5, alpha = 1) +
  stat_ellipse() +
  scale_color_manual(values = c("#882255", "#009E73")) +
  xlab(paste0("PCoA1 (", round(perc_expl1[1], 1), " %)")) +
  ylab(paste0("PCoA2 (", round(perc_expl1[2], 1), " %)")) +
  ggtitle("PERMANOVA Soil conditioning: p = 0.006") +
  theme(panel.background = element_rect(fill = "white"), #theme of the graph
        panel.border = element_blank(),
        axis.line = element_line(), #show line for x and y axis
        plot.title = element_text(size = 10, hjust = 0),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_text(vjust = 2.5), #, hjust = -0
        strip.text = element_text(size = 13),
        plot.margin = margin(5.5, 12, 5.5, 5.5, "points"),
        legend.position = "right")
```

![](FP_02_beta_diversity_OutliersRemoved_Batch2_final_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

#### PCoA labeling outliers

``` r
ggplot(plot_df, aes(x = MDS1, y = MDS2, color = Accession, shape = Soil_conditioning)) +
  geom_point(size = 1.5, alpha = 1) +
  geom_text(label = plot_df$SampleID, hjust = 0, vjust = 0, size = 2.5) +
  xlab(paste(round(perc_expl1[1], 1), "%")) +
  ylab(paste(round(perc_expl1[2], 1), "%")) +
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        legend.position = "right")
```

![](FP_02_beta_diversity_OutliersRemoved_Batch2_final_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

``` r
rm(mds_coords, mds_result, metadata_df, plot_df, bray_dist, perc_expl1)
```


## 2.1.7 - PCoA per accession
### Define the function 

``` r
stats <- c("0.010", "0.020", "0.025", "0.020")
stats <- lapply(stats, function(x) {
  paste0("PERMANOVA Soil: p = ", x)
})

color_fw <- c(rep("#0072B2", 1), rep("#D55E00", 3))
tcolor_fw <- c(rep("white", 1), rep("black", 3))
tsize_fw <- rep(18, 4)


MDS_listing <- function(physeq_list) {
  plot_list <- lapply(physeq_list, function(x) {
    # Calculate Bray-Curtis distance matrix
    bray_dist <- phyloseq::distance(x, method = "bray")
    
    # Run classical MDS with cmdscale()
    mds_result <- cmdscale(bray_dist, k = 2, eig = TRUE)
    
    # Extract coordinates
    mds_coords <- as.data.frame(mds_result$points)
    colnames(mds_coords) <- c("MDS1", "MDS2")
    
    # Add sample IDs as a column
    mds_coords$SampleID <- rownames(mds_coords)
    
    # Extract sample metadata
    metadata_df <- data.frame(sample_data(x))
    metadata_df$SampleID <- rownames(metadata_df)
    
    # Percentage explained variance
    perc_expl1 <- 100 * (mds_result$eig / sum(mds_result$eig))
    
    # Merge ordination coordinates with metadata
    plot_df <- merge(mds_coords, metadata_df, by = "SampleID")
  
    # Ordination plots
    p <- ggplot(plot_df, aes(x = MDS1, y = MDS2, color = Soil_conditioning, shape = Soil_conditioning)) +
      geom_point(size = 1.5, alpha = 1) +
      scale_color_manual(name = "Soil treatment", 
                         values = c("#882255", "#009E73"),
                         labels = c(Co = "Control", Mb = expression(italic("M. brassicae")))) +
      scale_shape_manual(name = "Soil treatment", 
                         values = c(Co = 16, Mb = 15),
                         labels = c(Co = "Control", Mb = expression(italic("M. brassicae")))) +
      xlab(paste0("PCoA1 (", round(perc_expl1[1], 1), " %)")) +
      ylab(paste0("PCoA2 (", round(perc_expl1[2], 1), " %)")) +
      theme(panel.background = element_rect(fill = "white"), #theme of the graph
          panel.border = element_blank(),
          axis.line = element_line(), #show line for x and y axis
          plot.title = element_text(size = 18, hjust = 0),
          plot.subtitle = element_text(size = 14),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 18),
          axis.title.x = element_text(vjust = -0.5),
          axis.title.y = element_text(vjust = 2.5), #, hjust = -0
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          #strip.text = element_text(size = 13),
          #plot.margin = margin(5.5, 12, 5.5, 5.5, "points"),
          plot.margin = margin(10, 20, 10, 20, "points"), # top . bottom .
          legend.position = "none")
    return(p)
  })
  
  ## Add titles using names of the original list
  titles <- names(plot_list)
  
plot_list_titled <- mapply(function(p, title, subtitle) {
    p + labs(title = title,
             subtitle = subtitle)
  }, p = plot_list, title = titles, subtitle = stats, SIMPLIFY = FALSE)

  plot_list_titled_elips <- mapply(function(p) {
    p + stat_ellipse()
  }, p = plot_list_titled, SIMPLIFY = FALSE)
  
  plot_list_batch <- mapply(function(p) {
    p + aes(shape = Batch)
  }, p = plot_list_titled, SIMPLIFY = FALSE)
  
  plot_list_batch_elips <- mapply(function(p) {
    p + aes(shape = Batch) + stat_ellipse()
  }, p = plot_list_titled, SIMPLIFY = FALSE)
    
  plot_list_batch_elips_legend <- mapply(function(p) {
    p + aes(shape = Batch) + stat_ellipse() + theme(legend.position = "right")
  }, p = plot_list_titled, SIMPLIFY = FALSE)
  
  # Facet wrap
  plot_list_facetwrap <- mapply(function(p, color_fw, tcolor_fw, tsize_fw) {
     strip0 <- strip_themed(background_x = elem_list_rect(fill = color_fw),
                            text_x = elem_list_text(colour = tcolor_fw, size = tsize_fw))
     p + facet_wrap2(~ Accession, #facet wrap2 needed to color strip
                    strip = strip0) +
       labs(title = NULL)
  }, p = plot_list_titled_elips, color_fw = color_fw, tcolor_fw = tcolor_fw, tsize_fw = tsize_fw, SIMPLIFY = FALSE)

  
  return(list(Plot_list = plot_list_titled, Plot_elips = plot_list_titled_elips, Plot_Batch = plot_list_batch, Plot_Batch_elipse = plot_list_batch_elips, Plot_Batch_elipse_legend = plot_list_batch_elips_legend, plot_list_facetwrap = plot_list_facetwrap))
}
```

### Run function and plot

``` r
# let's run our custom function
set.seed(5235)
MDS_Accessions_split <- MDS_listing(filtered_ps)
```

### Print plot
#### Four columns

``` r
# let's run our custom function
wrap_plots(MDS_Accessions_split$Plot_elips, ncol = 4) &
  theme(panel.spacing = unit(8, "lines")) #Allow larger margins ggplot
```

![](FP_02_beta_diversity_OutliersRemoved_Batch2_final_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

``` r
mapply(function(x)
  ggsave(plot = last_plot(), 
         filename = x,
         path = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs",
         scale = 2.2,
         width = 16,
         height = 4.3,
         units = "cm",
         dpi = 300),
  x = c("FP_BetaDiv_AccSep_Cat4.svg", "FP_BetaDiv_AccSep_Cat4.png"))
```

```
##                                                                                                                     FP_BetaDiv_AccSep_Cat4.svg 
## "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs/FP_BetaDiv_AccSep_Cat4.svg" 
##                                                                                                                     FP_BetaDiv_AccSep_Cat4.png 
## "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs/FP_BetaDiv_AccSep_Cat4.png"
```

#### Two columns

``` r
wrap_plots(MDS_Accessions_split$Plot_elips, ncol = 2) &
  theme(panel.spacing = unit(8, "lines")) #Allow larger margins ggplot
```

![](FP_02_beta_diversity_OutliersRemoved_Batch2_final_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

``` r
mapply(function(x)
  ggsave(plot = last_plot(), 
         filename = x,
         path = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs",
         scale = 2.2,
         width = 8,
         height = 8.4,
         units = "cm",
         dpi = 300),
  x = c("FP_BetaDiv_AccSep_Cat2.svg", "FP_BetaDiv_AccSep_Cat2.png"))
```

```
##                                                                                                                     FP_BetaDiv_AccSep_Cat2.svg 
## "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs/FP_BetaDiv_AccSep_Cat2.svg" 
##                                                                                                                     FP_BetaDiv_AccSep_Cat2.png 
## "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs/FP_BetaDiv_AccSep_Cat2.png"
```

#### Two columns with facet wrap

``` r
wrap_plots(MDS_Accessions_split$plot_list_facetwrap, ncol = 2) &
  theme(panel.spacing = unit(8, "lines")) #Allow larger margins ggplot
```

![](FP_02_beta_diversity_OutliersRemoved_Batch2_final_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

``` r
mapply(function(x)
  ggsave(plot = last_plot(), 
         filename = x,
         path = "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs",
         scale = 2.2,
         width = 8,
         height = 8.4,
         units = "cm",
         dpi = 300),
  x = c("FP_BetaDiv_AccSep_Cat2_FW.svg", "FP_BetaDiv_AccSep_Cat2_FW.png"))
```

```
##                                                                                                                     FP_BetaDiv_AccSep_Cat2_FW.svg 
## "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs/FP_BetaDiv_AccSep_Cat2_FW.svg" 
##                                                                                                                     FP_BetaDiv_AccSep_Cat2_FW.png 
## "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/Result/Final graphs/FP_BetaDiv_AccSep_Cat2_FW.png"
```


## 2.1.10 - PCoA per accession outliers labeled
### Define the function 

``` r
MDS_listing <- function(physeq_list) {
  plot_list <- lapply(physeq_list, function(x) {
    # Calculate Bray-Curtis distance matrix
    bray_dist <- distance(x, method = "bray")
    
    # Run classical MDS with cmdscale()
    mds_result <- cmdscale(bray_dist, k = 2, eig = TRUE)
    
    # Percentage explained variance per axis
    perc_expl1 <- 100 * (mds_result$eig / sum(mds_result$eig))
    
    # Extract coordinates
    mds_coords <- as.data.frame(mds_result$points)
    colnames(mds_coords) <- c("MDS1", "MDS2")
    
    # Add sample IDs as a column
    mds_coords$SampleID <- rownames(mds_coords)
    
    # Extract sample metadata
    metadata_df <- data.frame(sample_data(x))
    metadata_df$SampleID <- rownames(metadata_df)
    
    # Merge ordination coordinates with metadata
    plot_df <- merge(mds_coords, metadata_df, by = "SampleID")
  
    # Ordination plots
    p <- ggplot(plot_df, aes(x = MDS1, y = MDS2, color = Soil_conditioning)) +
      geom_point(size = 1.5, alpha = 1) +
      geom_text(label = plot_df$SampleID, hjust = 0.5, vjust = 0, size = 2.5) +
      scale_color_manual(values = c("#0072B2", "#D55E00")) + #"#56B4E9", "#009E73"
      xlab(paste0("PCoA1 (", round(perc_expl1[1], 1), " %)")) +
      ylab(paste0("PCoA2 (", round(perc_expl1[2], 1), " %)")) +
      theme(panel.background = element_rect(fill = "white"), #theme of the graph
          panel.border = element_blank(), #no border around the graph (for border, use: element_rect(color = "gray30", fill = NA))
          axis.line = element_line(), #show line for x and y axis
          plot.title = element_text(size = 14, hjust = 0),
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 14),
          axis.title.x = element_text(vjust = -0.5),
          axis.title.y = element_text(vjust = 2.5), #, hjust = -0
          strip.text = element_text(size = 13),
          plot.margin = margin(5.5, 12, 5.5, 5.5, "points"),
          legend.position = "none")
    return(p)
  })
  
  ## Add titles using names of the original list
  titles <- names(plot_list)
  
   plot_list_titled <- mapply(function(p, title) {
    p + ggtitle(title)
  }, p = plot_list, title = titles, SIMPLIFY = FALSE)
  
  plot_list_titled_nolegend <- mapply(function(p) {
    p + theme(legend.position = "none")
  }, p = plot_list_titled, SIMPLIFY = FALSE)
  
  return(list(Plot_list = plot_list_titled, Plot_list_nolegend = plot_list_titled_nolegend))
}
```

### Run function

``` r
# let's run our custom function
set.seed(5235)
MDS_Accessions_split <- MDS_listing(filtered_ps)

wrap_plots(MDS_Accessions_split$Plot_list, ncol = 5)
```

![](FP_02_beta_diversity_OutliersRemoved_Batch2_final_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

``` r
wrap_plots(MDS_Accessions_split$Plot_list_nolegend, ncol = 5)
```

![](FP_02_beta_diversity_OutliersRemoved_Batch2_final_files/figure-html/unnamed-chunk-21-2.png)<!-- -->


``` r
rm(MDS_Accessions_split)
```


# 2.3 - P values for beta diversity
A permutation anova will tell if the differences in the microbial community structure are significant or not. They will essentially help you separate the data clouds of your ordination with confidence levels. It will compare centroids and spread of the different groups. It is based on the distance matrix and does not assume normality of the data (non-parametric). Be careful when beta dispersion is different for different groups. It may give misleading low p-values when centroids are apart. [More info](https://archetypalecology.wordpress.com/2018/02/21/permutational-multivariate-analysis-of-variance-permanova-in-r-preliminary/).

## 2.3.1 - Run PERMANOVA at ASV level by terms
You will need to run, test and check several different models and data slices to have final insight into the data set you are evaluating. get used with testing multiple models!

### Without cat treatment Co

``` r
metadata_Bac <- as(sample_data(FP_CSS_bac_ps),"data.frame")

set.seed(5235)
dist <- phyloseq::distance(t(otu_table(FP_CSS_bac_ps)), method = "bray")

set.seed(5235)
Acc_per <- adonis2(dist ~ Accession*Soil_conditioning, data = metadata_Bac, method = "bray", by = "terms", permutations = how(blocks = metadata_Bac$Block_g))
Acc_per
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Blocks:  metadata_Bac$Block_g 
## Permutation: free
## Number of permutations: 199
## 
## adonis2(formula = dist ~ Accession * Soil_conditioning, data = metadata_Bac, permutations = how(blocks = metadata_Bac$Block_g), method = "bray", by = "terms")
##                             Df SumOfSqs      R2      F Pr(>F)   
## Accession                    3   0.7117 0.13589 2.3324  0.005 **
## Soil_conditioning            1   0.1554 0.02968 1.5282  0.005 **
## Accession:Soil_conditioning  3   0.5051 0.09644 1.6553  0.005 **
## Residual                    38   3.8649 0.73799                 
## Total                       45   5.2371 1.00000                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
Permanova_Acc_FP <- as.data.frame(Acc_per)[1:3, ]
write.csv(Permanova_Acc_FP, "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/Data/Statistical_output/FP/Permanova_Acc_FP.csv")
```


### PERMANOVA per accession

``` r
Perm_per_acc <- lapply(filtered_ps, function(x) {
  metadata_Bac <- as(sample_data(x),"data.frame")
  
  set.seed(5235)
  dist <- phyloseq::distance(t(otu_table(x)), method = "bray")
  
  set.seed(5235)
  perm_soil <- adonis2(dist ~ Soil_conditioning, data = metadata_Bac, method = "bray", by = "terms", permutations = how(blocks = metadata_Bac$Block_g))
  
  return(Soil_conditioning = perm_soil)
})
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

``` r
names(Perm_per_acc) <- names(filtered_ps)

Perm_per_acc
```

```
## $VL
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Blocks:  metadata_Bac$Block_g 
## Permutation: free
## Number of permutations: 199
## 
## adonis2(formula = dist ~ Soil_conditioning, data = metadata_Bac, permutations = how(blocks = metadata_Bac$Block_g), method = "bray", by = "terms")
##                   Df SumOfSqs      R2      F Pr(>F)   
## Soil_conditioning  1  0.14758 0.13296 1.3802   0.01 **
## Residual           9  0.96236 0.86704                 
## Total             10  1.10994 1.00000                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## $CD
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Blocks:  metadata_Bac$Block_g 
## Permutation: free
## Number of permutations: 199
## 
## adonis2(formula = dist ~ Soil_conditioning, data = metadata_Bac, permutations = how(blocks = metadata_Bac$Block_g), method = "bray", by = "terms")
##                   Df SumOfSqs      R2      F Pr(>F)  
## Soil_conditioning  1  0.15508 0.12821 1.4707   0.02 *
## Residual          10  1.05453 0.87179                
## Total             11  1.20961 1.00000                
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## $HM
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Blocks:  metadata_Bac$Block_g 
## Permutation: free
## Number of permutations: 863
## 
## adonis2(formula = dist ~ Soil_conditioning, data = metadata_Bac, permutations = how(blocks = metadata_Bac$Block_g), method = "bray", by = "terms")
##                   Df SumOfSqs      R2      F Pr(>F)  
## Soil_conditioning  1  0.23701 0.20673 2.3455  0.025 *
## Residual           9  0.90947 0.79327                
## Total             10  1.14648 1.00000                
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## $GO1
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Blocks:  metadata_Bac$Block_g 
## Permutation: free
## Number of permutations: 199
## 
## adonis2(formula = dist ~ Soil_conditioning, data = metadata_Bac, permutations = how(blocks = metadata_Bac$Block_g), method = "bray", by = "terms")
##                   Df SumOfSqs      R2      F Pr(>F)  
## Soil_conditioning  1  0.12083 0.11405 1.2874   0.02 *
## Residual          10  0.93857 0.88595                
## Total             11  1.05940 1.00000                
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
df_acc <- lapply(names(Perm_per_acc), function(x){
  df <- as.data.frame(Perm_per_acc[[x]])[1, ]
  df$Accession <- rep(x, 1)
  df$Treatment <- rownames(Perm_per_acc[[x]])[1]
  return(df)
})
names(df_acc) <- names(Perm_per_acc)
Permanova_perAcc_FP <- do.call(rbind, df_acc)[c(6,7,1:5)]
Permanova_perAcc_FP
```

```
##     Accession         Treatment Df  SumOfSqs        R2        F Pr(>F)
## VL         VL Soil_conditioning  1 0.1475818 0.1329636 1.380187  0.010
## CD         CD Soil_conditioning  1 0.1550844 0.1282103 1.470656  0.020
## HM         HM Soil_conditioning  1 0.2370145 0.2067321 2.345474  0.025
## GO1       GO1 Soil_conditioning  1 0.1208280 0.1140535 1.287363  0.020
```

``` r
write.csv(Permanova_perAcc_FP, "C:/Users/kreek001/OneDrive - Wageningen University & Research/Chapters/Experimental Chapter 3/RScripts/GitHub/Data/Statistical_output/FP/Permanova_perAcc_FP.csv")
```

## 2.3.5 pairwise PERMANOVA
Now that we know the fixed factor effects and interactions, let's make pairwise comparisons.
we will use 2 different functions - EcolUtils::adonis.pair and pairwiseAdonis::pairwise.adonis2.

### EcolUtils::adonis.pair - with lists and p adjust, no blocks.
This is not a very good package to test pairwise comparisons, but there is no better one.
There is a minimum value for a p-value when using specific number of permutations. This can result in identical p-values when having a small data set. Increasing the number of permutations can solve the problem. [See more info](https://github.com/vegandevs/vegan/issues/339).

### Function adonis.pairs with adonis2
Check this: https://rdrr.io/github/GuillemSalazar/EcolUtils/man/adonis.pair.html 


``` r
adonis.pair2 <- function(dist.mat, Factor, nper = 999, p.adjust.m = "BH") {
  comb.fact <- combn(levels(Factor), 2)
  p.val <- numeric(ncol(comb.fact))
  R2 <- numeric(ncol(comb.fact))

  for (i in 1:ncol(comb.fact)) {
    sub.ind <- Factor %in% comb.fact[, i]
    sub.dist <- as.dist(as.matrix(dist.mat)[sub.ind, sub.ind])
    sub.factor <- droplevels(Factor[sub.ind])

    adonis.res <- vegan::adonis2(sub.dist ~ sub.factor, permutations = nper)

    p.val[i] <- adonis.res$`Pr(>F)`[1]
    R2[i] <- adonis.res$R2[1]
  }

  comp <- paste(comb.fact[1, ], "vs", comb.fact[2, ])
  out <- data.frame(comp = comp, p.value = p.val, R2 = R2,
                    p.adjusted = p.adjust(p.val, method = p.adjust.m))
  return(out)
}
```

### Run pairwise comparisons

``` r
set.seed(5235)
dist <- phyloseq::distance(t(otu_table(FP_CSS_bac_ps)), method = "bray")

adonis.pair2(dist.mat = dist,
                        Factor = as.factor(as(phyloseq::sample_data(FP_CSS_bac_ps),"data.frame")$Accession),
                        nper = 9999)
```

```
##        comp p.value         R2 p.adjusted
## 1  VL vs CD  0.0001 0.10431644    0.00012
## 2  VL vs HM  0.0001 0.08789819    0.00012
## 3 VL vs GO1  0.0001 0.13150323    0.00012
## 4  CD vs HM  0.0064 0.07016572    0.00640
## 5 CD vs GO1  0.0001 0.07748724    0.00012
## 6 HM vs GO1  0.0001 0.09814504    0.00012
```

# 2.4 Check library size 
## Original library size

``` r
# check library size
a <- sample_data(FP_CSS_bac_ps)[ , c("library_size", "Accession")]
a$Sample <- rownames(a)
a <- a[order(a$library_size), ]
head(a, 20)
```

```
##      library_size Accession Sample
## F019        17146        CD   F019
## F200        17396        HM   F200
## F175        19278        HM   F175
## F098        19281       GO1   F098
## F037        19605        CD   F037
## F478        20116        VL   F478
## F479        20388        VL   F479
## F460        20638        VL   F460
## F456        20754        VL   F456
## F096        20767       GO1   F096
## F457        20845        VL   F457
## F477        20865        VL   F477
## F099        20965       GO1   F099
## F199        21369        HM   F199
## F018        21569        CD   F018
## F015        21624        CD   F015
## F035        21820        CD   F035
## F017        21839        CD   F017
## F039        21916        CD   F039
## F036        21956        CD   F036
```

``` r
a[order(a$Accession), 1:2]
```

```
##      library_size Accession
## F478        20116        VL
## F479        20388        VL
## F460        20638        VL
## F456        20754        VL
## F457        20845        VL
## F477        20865        VL
## F475        21963        VL
## F458        22498        VL
## F459        23742        VL
## F455        24097        VL
## F476        24334        VL
## F019        17146        CD
## F037        19605        CD
## F018        21569        CD
## F015        21624        CD
## F035        21820        CD
## F017        21839        CD
## F039        21916        CD
## F036        21956        CD
## F016        22693        CD
## F038        24015        CD
## F040        27933        CD
## F020        33821        CD
## F200        17396        HM
## F175        19278        HM
## F199        21369        HM
## F195        22858        HM
## F193        23189        HM
## F194        23590        HM
## F179        23949        HM
## F176        24313        HM
## F197        24592        HM
## F177        24934        HM
## F180        24989        HM
## F098        19281       GO1
## F096        20767       GO1
## F099        20965       GO1
## F097        22647       GO1
## F116        22851       GO1
## F095        23208       GO1
## F120        23244       GO1
## F119        23633       GO1
## F118        23987       GO1
## F100        24174       GO1
## F117        25135       GO1
## F115        25238       GO1
```

``` r
summary(a)
```

```
##   library_size   Accession    Sample         
##  Min.   :17146   VL :11    Length:46         
##  1st Qu.:20890   CD :12    Class :character  
##  Median :22670   HM :11    Mode  :character  
##  Mean   :22647   GO1:12                      
##  3rd Qu.:24008                               
##  Max.   :33821
```

``` r
table(a$Accession)
```

```
## 
##  VL  CD  HM GO1 
##  11  12  11  12
```

``` r
aggregate(a$library_size, by = list(a$Accession), FUN = mean)
```

```
##   Group.1        x
## 1      VL 21840.00
## 2      CD 22994.75
## 3      HM 22768.82
## 4     GO1 22927.50
```
