# Analysis part 2

# Normalisation and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression
#
# Key packages used in this analysis were sctransform_0.2.1 and Seurat_3.1.5
#
# Data for this analysis were generated from two samples (WT = wild type, MT = mutant) that were sequenced together.
#
# sctransform package was chosen because of very fast computational approach and capability to remove technical noise while preserving biological heterogeneity of the data.
#
# Load libraries
library(scater)
library(sctransform)
library(Seurat)

# sce object was transformed into Seurat object. During this transformaion were preserved only cells that express at least 200 genes. Only features that are expressed in at least 3 cells were preserved.
counts <- as.matrix(sce@assays$data$counts)
rownames(counts) <- rownames(sce)

seu <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)

# adding meta data (type WT or MT) to the seurat object
type_1 <- as.data.frame(t(sce$type))
type_2 <- type[rownames(type_1) %in% rownames(seu@meta.data),]
type <- as.data.frame(type_2)
seu <- AddMetaData(seu, type)

# normalisation with sctransform
seu <- SCTransform(object = seu, verbose = FALSE)

