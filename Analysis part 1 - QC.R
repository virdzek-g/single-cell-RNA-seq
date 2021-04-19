# Analysis part 1

# Quality control of single cell RNA seq data.
#
# This part of the analysis heavily relied on scater_1.12.2 R package.
#
# Data for this analysis were generated from two samples (WT = wild type, MT = mutant) that were sequenced together.
#
# These steps were adapted from: Lun ATL, McCarthy DJ and Marioni JC. A step-by-step workfrlow for low-level analysis of single-cell RNA-seq data with Bioconductor F1000Research 2016, 5:2122 (doi:10.12688/f1000research.9501.2)
#
# To determine appropriate number fo cells and genes to keep, I performed  visualisation of outliers and low-quality cells. I advice visualising data before and after quality control to determine right parameters for your data.
#
#
#
# Create single cell experiment object
data <- dplyr::inner_join(WT, MT, by='genes')
counts <- as.matrix(data[-1,-1])
library(scater)
sce <- SingleCellExperiment(list(counts=counts))
rownames(sce) <- rownames(counts)
sce
#dim: 12611 3200  

# Addition of meta data to the sce object. type corresponds either to WT or MT
type <- data[1,]
type <- type[,-1]
colData(sce)$type <- type

# Compute metrics
sce <- calculateQCMetrics(sce)

# Keep only genes which are expressed in 1 or more cells
keep_feature <- rowSums(counts(sce) > 0) > 0
sce <- sce[keep_feature,]
sce # 12611 3200

# Library size was used as a cell quality control measure. Library size represents total sum of counts in the cell across all features. Here were identified cells with log-library size more than 1 median absolute deviations (MADs) above and below median of the log-library size as a low quality cells. MADs was chosen according to quality control metrics shown in ‘Visualization of QC results’ section
libsize.drop <- isOutlier(sce$total_counts, nmads=1, type="both", log=TRUE)

# Cells with high number of genes with zero counts is likely to have a poor quality. Cells with a log-trnasformed number of expressed genes more than 1 MADs above and below the median were therefore removed
feature.drop <- isOutlier(sce$total_features_by_counts, nmads=1, type="both", log=TRUE)

sce <- sce[,!(libsize.drop | feature.drop)]

# Information about removed and kept cells
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), Remaining=ncol(sce))

#ByLibSize ByFeature Remaining
#  1120      1122      2034

# Keep only genes with average count equal or higher than 0.5. Low abundance genes below this threshold would not contribute with enough information for downstream analysis and could interfere with statistical exploration of the data
ave.counts <- calcAverage(sce)
keep <- ave.counts >= 0.01
summary(keep)
sce  <- sce[keep,]
#Mode   FALSE    TRUE 
#logical   9269    8894 


### Visualisation of QC results

#Examining the most expressed features
plotHighestExprs(sce, exprs_values = "counts")

# Frequency of gene expression as a function of mean
plotExprsFreqMean(sce)

# PCA
library(ggplot2)
sce <- runPCA(sce, use_coldata = TRUE, detect_outliers = TRUE)
dims <- sce@reducedDims$PCA_coldata
batch <- as.data.frame(sce$batch)
batch <- as.data.frame(t(batch))
rownames(dims) <- batch$batch
pca <- princomp(dims)
groups <- factor(rownames(pca$scores))
plot <- ggplot2:: qplot(Comp.1, Comp.2, data=pca, colour=groups, xlab="Dimension 1", ylab=="Dimension 2")

# Frequency of library size and number of expressed genes
par(mfrow=c(1,2))
hist(sce$total_counts, breaks = 50, col="lightblue", cex.main=0.7, cex.lab=0.7, xlab='Library size', main=NULL)
hist(sce$total_features_by_counts, breaks = 50, col="lightgreen", cex.main=0.7, cex.lab=0.7, xlab='Number of expressed genes', main=NULL)

# Histogram of log average counts for all genes
hist(log10(ave.counts), breaks=100, main="", col="grey80", xlab=expression(Log[10]~"average count"))



