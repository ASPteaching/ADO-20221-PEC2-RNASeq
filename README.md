# Infographics Data visualization on Bioinformatics.

## Introduction

This repository contains a simplified analysis of RNA-Seq data with R and Bioconductor. The objective is to show a workflow in a simple way that can serve as the basis for an infographic that will be made with the [GRBio](http://grbio.upc.edu) group.

Some graphics than you can see in the [infographics](https://www.canva.com/design/DAFDMgmnNBw/ArYoC2fGU1OAQYMg_X0jAA/edit?utm_content=DAFDMgmnNBw&utm_campaign=designshare&utm_medium=link2&utm_source=sharebutton) are below.

## Quality control &  Exploration

Boxplots and other graphics help to check data distributions.
Ideally, one might expect that samples tend to be more similar within groups than between groups. Distinct techniques such as PCA or Hierarchical clustering are used to check this assumption.

### Boxplot of log counts

![image](figures/boxplot.png)

### MDS plot

The library [emojifont](https://cran.r-project.org/web/packages/emojifont/vignettes/emojifont.html) was used for this plot

![image](figures/MDS_emotis.png)

### Dendogram

![image](figures/dendograma.png)


## Data analysis

A statistical test allows selecting significant differentially expressed genes. The volcano plot shows statistical versus biological significance. 

Heatmap displays the expressions of selected genes in a grid (genes in rows & samples in columns). The color scale reflects the intensity of gene expression in each sample. On the margins, dendrograms group genes or samples based on the similarity of their gene expression pattern. This is useful for identifying genes that are commonly regulated.

### Volcano plot

![image](figures/vp2.png)


### Heatmap

![image](figures/pheatmap_transp.png)

## Biological interpretation

Genes are annotated in different knowledge databases by terms or categories describing their biological role. 
The distribution of annotations of selected genes is compared with the distribution of the same annotations in the genome. This allows determining which biological processes might be associated with our gene list.
We end up linking those differentially expressed genes that are included in the most represented biological categories using a network plot (bottom left-hand figure)


### Dotplot

![image](figures/dotplot_transp.png)

### CNET plot

![image](figures/cnetplot_transp.png)

