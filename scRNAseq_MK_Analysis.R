' 
April 1 2019
------------
MK Single Cell RNAseq

Resources
---------
https://satijalab.org/seurat/pbmc3k_tutorial.html
https://davetang.org/muse/2017/08/01/getting-started-seurat/
'
install.packages("VGAM")
install.packages("devtools")
remotes::install_github("HenrikBengtsson/seurat", ref="feature/hdf5r-optional")

library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)

setwd("/Volumes/Lesleys_Projects/Morrell_Lab/MK_Single_Cell_RNA/data/html_files/deliv_Morrell_031319_cellRanger/Sample_Lung/outs/filtered_feature_bc_matrix")
data_bm.data <- Read10X(data.dir = "/Volumes/Lesleys_Projects/Morrell_Lab/MK_Single_Cell_RNA/data/html_files/deliv_Morrell_031319_cellRanger/Sample_Lung/outs/filtered_feature_bc_matrix/")
# Change 'features.tsv' to 'genes.tsv'
class(data_bm.data)
#attr(,"package")
# 31,053 genes and 1,985 single cells
dim(data_bm.data)
# check out the first six genes and cells
data_bm.data[1:6,1:6]
MK_genes <- data_bm[c("Itga2b", "Pf4", "Vwf"),]
MK_genes
# summary of total expression per single cell
summary(colSums(data_bm.data))
# check how many genes have at least one transcript in each cell
at_least_one <- apply(data_bm.data, 2, function(x) sum(x>0))
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")

hist(colSums(data_bm.data),
     breaks = 100, main = "Expression sum per cell",
     xlab = "Sum expression")
'
# LMC TODO: 
- How to find median number of detected genes within each single cell using histogram  
- Also find median sum of expression among single cells
'

'
Select Cells for further analysis
- Keep genes that are express in 3 or more cells
- Keep cells that contain at least 200 genes
- Talk to Craig, this cutoff can change
'

# Find number of genes detected in at least 3 or more cells
# Roughly 1/2 of the genes are expressed in 3 or more cells [FALSE:16890 TRUE:14163 ]
tmp <- apply(data_bm.data, 1, function(x) sum(x>0))
table(tmp>=3)

bm_cells <- CreateSeuratObject(raw.data = data_bm.data,
                               min.cells = 1,
                               min.genes = 10,
                               project = "10X_BM_Cells")
class(bm_cells)
bm_cells
slotNames(data_bm.data)

bm_cells2 <- CreateSeuratObject(raw.data = data_bm.data,
                                min.cells = 3,
                                min.genes = 200,
                                project = "10X_BM_Cells")
class(bm_cells2)
bm_cells2
slotNames(data_bm.data)

' 
Remove mitochondrial RNA
- It looks there are no mitochondial transcripts
'

mito.genes <- grep(pattern = "^MT-", x = rownames(x = bm_cells@data), value = TRUE)
length(mito.genes)
percent.mito <- Matrix::colSums(bm_cells@raw.data[mito.genes, ]) / Matrix::colSums(bm_cells@raw.data)

# Add mitochondrial RNA statistics for later violin plot
bm_cells <- AddMetaData(object = bm_cells,
                        metadata = percent.mito,
                        col.name = "percent.mito")

# plot number of genes, UMIs, and % mitochondria
# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
VlnPlot(object = bm_cells,
        features.plot = c("nGene", "nUMI", "percent.mito"),
        nCol = 3)

# NOTE TODO: Theres no mitochondrial RNA, you may not want to run that section/plot
# This plot is showing number genes captured vs number of unique molecules
par(mfrow = c(1, 1))
#GenePlot(object = bm_cells, gene1 = "nUMI", gene2 = "percent.mito", pch.use = '.')
GenePlot(object = bm_cells, gene1 = "nUMI", gene2 = "nGene", pch.use = '.')

'
Filter Cells
- Remove cells that contain less than 200genes
- Remove cells that have unique gene counts over 5000 
- 202 cells were filtered out because they had more than 5000 genes
- Check with Craig on these numbers
'

table(bm_cells@meta.data$percent.mito < 0.05 & bm_cells@meta.data$nGene<5000)
# Filter cells using the FilterCells() function
bm_cells <- FilterCells(object = bm_cells,
                        subset.names = c("nGene", "percent.mito"),
                        low.thresholds = c(200, -Inf),
                        high.thresholds = c(5000, 0.05))
bm_cells

'
Normalize Data
- This step is done so that the values from each cell can be compared with other cells
- Gene expression values for each cell are normalised by the total expression in each cell
- By default, global-scaling normalization method “LogNormalize” is used
- Normalizes the gene expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
'
hist(colSums(bm_cells@data),
     breaks = 100,
     main = "Total expression before normalization",
     xlab = "Sum of expression")

bm_cells <- NormalizeData(object = bm_cells,
                          normalization.method = "LogNormalize",
                          scale.factor = 1e4)

hist(colSums(bm_cells@data),
     breaks = 100,
     main = "Total expression after normalization",
     xlab = "Sum of expression")

'
Detect variable genes across single cells
- take each gene, get the average expression and variance of the gene across the 1618 cells
- categorise genes into bins (default is 20) based on their expression and variance
- finally normalise the variance in each bin. This was the same approach in Macosko et al. 
- parameters used below are typical settings for UMI data that is normalised to a total of 10,000 molecules and will identify around 2,000 variable genes.
- TODO: explore different parameters
'
bm_cells@var.genes
logical(0)
bm_cells <- FindVariableGenes(object = bm_cells,
                              mean.function = ExpMean,
                              dispersion.function = LogVMR,
                              x.low.cutoff = 0.0125,
                              x.high.cutoff = 3,
                              y.cutoff = 0.5)
# vector of variable genes
head(bm_cells@var.genes)
# number of variable genes
length(bm_cells@var.genes)
# store mean and variance for each gene
head(bm_cells@hvg.info)

'
Scale Data and Remove Unwanted Variation
----------------------------------------
- single cell dataset likely contains ‘uninteresting’ sources of variation
- This could include: technical noise, batch effects, or even biological sources of variation (cell cycle stage)
- Removing these signals could improve downstream dimensionality reduction and clustering
- scaled z-scored residuals of these models are stored in the scale.data slot, and are used for dimensionality reduction and clustering
Seurat v2.0 implements regression as part of the data scaling process 
-Seurat constructs linear models to predict gene expression based on user-defined variables to help remove unwanted sources of variation.
- The following can be regressed out:
cell-cell variation in gene expression driven by batch (if applicable)
cell alignment rate (as provided by Drop-seq tools for Drop-seq data)
the number of detected molecules, and mitochondrial gene expression
For cycling cells, we can also learn a ‘cell-cycle’ score
- FYI: RegressOut function has been deprecated, and replaced with the vars.to.regress argument in ScaleData.
- Ref: Buettner et al, NBT, 2015 --> this group showed how correcting for these confounding factors improved their downstream analysis
'

bm_cells@scale.data
# build linear model using nUMI and percent.mito
bm_cells <- ScaleData(object = bm_cells,
                      vars.to.regress = c("nUMI", "percent.mito"))
class(bm_cells@scale.data)
bm_cells@scale.data[1:6, 1:6]
'
PCA
___
# RunPCA() function performs the PCA on genes in the @var.genes slot by default and this can be changed using the pc.genes parameter
- ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
- used to identify markers
# that are strongly correlated with cellular heterogeneity, but may not have
# passed through variable gene selection
-The PCHeatmap() function produces a heatmap based on the PCA; by default the function uses the first principal component and plots 30 genes across the number of cells specified in cells
'
bm_cells <- RunPCA(object = bm_cells,
                   pc.genes = bm_cells@var.genes,
                   do.print = TRUE,
                   pcs.print = 1:5,
                   genes.print = 5)

PrintPCAParams(bm_cells)
PrintPCA(object = bm_cells, pcs.print = 1:2, genes.print = 5, use.full = FALSE)
VizPCA(object = bm_cells, pcs.use = 1:2)
PCAPlot(object = bm_cells, dim.1 = 1, dim.2 = 2)
bm_cells <- ProjectPCA(object = bm_cells, do.print = FALSE)
# HeatMap
# Plot 500 Cells
PCHeatmap(object = bm_cells,
          pc.use = 1,
          cells.use = 500,
          do.balanced = TRUE,
          label.columns = FALSE)
#plot all cells
PCHeatmap(object = bm_cells,
          pc.use = 1,
          do.balanced = TRUE,
          label.columns = FALSE)

'
Determine statistically significant princial components
-------------------------------------------------------
- determine how many principal components to use in downstream analyses 
- each PC essentially representing a ‘metagene’ that combines information across a correlated gene set
2 Methods Used"

JackStraw() 
- test based on a random null model, which is time-consuming for large datasets due to the resampling and may not return a clear cutoff and the other is a commonly used heuristic
- End result is a p-value for each gene’s association with each principal component
- this step compares the distribution of p-values for each principal component with a uniform distribution (dashed line)
- “Significant” principal components will show a strong enrichment of genes with low p-values (solid curve above the dashed line)

PCElbowPlot()
- Examine the standard deviations of the principle components
'
# Method 1
bm_cells <- JackStraw(object = bm_cells, 
                      num.replicate = 100, 
                      display.progress = FALSE)
JackStrawPlot(object = bm_cells, PCs = 1:12)
# Method 2
PCElbowPlot(object = bm_cells)

'
Cluster Cells
------------
- cells are embedded in a graph structure (e.g. a K-nearest neighbour (KNN) graph) with edges drawn between cells with similar gene expression patterns
- attempt to partition this graph into highly interconnected “quasi-cliques” or “communities.” As in PhenoGraph, we first construct a KNN graph based on the Euclidean distance in PCA space
- refine the edge weights between any two cells based on the shared overlap in their local neighbourhoods (Jaccard distance). 
- TODO: consider using more than 10 principal components
'
bm_cells <- FindClusters(object = bm_cells,
                         reduction.type = "pca",
                         dims.use = 1:10,
                         resolution = 0.6,
                         print.output = 0,
                         save.SNN = TRUE)
PrintFindClustersParams(object = bm_cells)
bm_cells <- RunTSNE(object = bm_cells,
                    dims.use = 1:10,
                    do.fast = TRUE)

TSNEPlot(object = bm_cells, do.label = TRUE)

# find all markers of cluster 0
# Possibly monocytes: Ccr2
cluster0.markers <- FindMarkers(object = bm_cells,
                                ident.1 = 0,
                                min.pct = 0.25)
head(cluster0.markers)

# find all markers of cluster 1
# May be neutrophils: LCN2, NGP, ORM1
cluster1.markers <- FindMarkers(object = bm_cells,
                                ident.1 = 1,
                                min.pct = 0.25)
head(cluster1.markers)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object = bm_cells,
                                ident.1 = 5,
                                ident.2 = c(0,3),
                                min.pct = 0.25)
head(cluster5.markers)
# find all markers distinguishing cluster 0 from remaining clusters
cluster0.markers <- FindMarkers(object = bm_cells,
                                ident.1 = 0,
                                ident.2 = c(1,2,3,4,5,6,7,8,9,10),
                                min.pct = 0.25)
head(cluster0.markers)

# find all markers distinguishing cluster 1 from remaining clusters
cluster1.markers <- FindMarkers(object = bm_cells,
                                 ident.1 = 1,
                                 ident.2 = c(0,2,3,4,5,6,7,8,9),
                                 min.pct = 0.25)
head(cluster1.markers)

# find all markers distinguishing cluster 2 from remaining clusters
cluster2.markers <- FindMarkers(object = bm_cells,
                                ident.1 = 2,
                                ident.2 = c(2,1,3,4,5,6,7,8,9),
                                min.pct = 0.25)
head(cluster2.markers)


# find markers for every cluster compared to all remaining cells, report
# only the positive ones
bm_cells.markers <- FindAllMarkers(object = bm_cells,
                                   only.pos = TRUE,
                                   min.pct = 0.25,
                                   thresh.use = 0.25)
head(bm_cells.markers)
bm_cells.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

'
# Seurat has several tests for differential expression which can be set with the test.use parameter

FindMarkers() function:

ROC test
t-test
LRT test based on zero-inflated data
LRT test based on tobit-censoring models

'
levels(bm_cells@ident)
table(bm_cells@ident)

my_bimod <- FindMarkers(object = bm_cells,
                        ident.1 = 1,
                        thresh.use = 0.25,
                        test.use = "bimod",
                        only.pos = TRUE)

my_roc <- FindMarkers(object = bm_cells,
                      ident.1 = 1,
                      thresh.use = 0.25,
                      test.use = "roc",
                      only.pos = TRUE)

my_t <- FindMarkers(object = bm_cells,
                    ident.1 = 1,
                    thresh.use = 0.25,
                    test.use = "t",
                    only.pos = TRUE)

my_tobit <- FindMarkers(object = bm_cells,
                        ident.1 = 1,
                        thresh.use = 0.25,
                        test.use = "tobit",
                        only.pos = TRUE)

# identical set of genes
dim(my_bimod)
dim(my_roc)
dim(my_t)
dim(my_tobit)

# the rankings of the genes are quite similar between the methods
my_gene <- row.names(my_bimod)
a <- 1:length(my_gene)
b <- match(my_gene, row.names(my_roc))
c <- match(my_gene, row.names(my_t))
d <- match(my_gene, row.names(my_tobit))

# bimod vs. bimod
cor(a, a, method = "spearman")

# bimod vs. roc
cor(a, b, method = "spearman")

# bimod vs. t
cor(a, c, method = "spearman")

# bimod vs. tobit
cor(a, d, method = "spearman")

par(mfrow=c(2,2))
barplot(a, main = 'bimod')
barplot(b, main = 'roc')
barplot(c, main = 't')
barplot(d, main = 'tobit')

VlnPlot(object = bm_cells, features.plot = c("Matk", "Vwf"))
VlnPlot(object = bm_cells, features.plot = c("Matk", "Pf4"))
head(my_roc)
head(my_tobit)
VlnPlot(object = bm_cells, features.plot = c("Gata1", "Itga2b"))
VlnPlot(object = bm_cells, features.plot = c("Matk", "Itga2b"))

VlnPlot(object = bm_cells,
        features.plot = c("Lcn2", "Pf4"),
        use.raw = TRUE,
        y.log = TRUE)

FeaturePlot(object = bm_cells,
            features.plot = c("Pf4", "Matk", "Itga2b", "Vwf", "Cxcr4", "F8a", "Cd36", "Anxa1", "Vwf"),
            cols.use = c("grey", "blue"),
            reduction.use = "tsne")

head(bm_cells.markers)

# Differential Gene Expression
MK_markers <- FindMarkers(object = bm_cells, ident.1 = 10, ident.2 = NULL, only.pos = TRUE)
# view results
head(x = MK_markers, 15)








bm_cells.markers %>%
  group_by(cluster) %>%
  top_n(10, avg_logFC) -> top10
head(top10)

DoHeatmap(object = bm_cells,
          genes.use = top10$gene,
          slim.col.label = TRUE,
          remove.key = TRUE)

FeaturePlot(object = bm_cells,
            features.plot = c("IL7R", "CD14", "LYZ", "MS4A1", "CD8A", "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP"),
            cols.use = c("grey", "blue"),
            reduction.use = "tsne")
