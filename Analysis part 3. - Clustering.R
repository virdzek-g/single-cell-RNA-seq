# Analysis part 3

# Linear dimensionality reduction with PCA, graph-based clustering based on identified principal components and random forrest regression
#
# This analysis is done on filtered and normalised data


# Load the libraries
library(Seurat) # Seurat_3.1.5
library(randomForest) # randomForest_4.6-14
library(Rtsne) # Rtsne_0.15

# Linear dimensional reduction
# PCA scores are used for cell clustering.
seu <- RunPCA(object = seu, verbose = FALSE, seed.use=30)
# Primary source of heterogeneity can be explored by visualizing cells and genes with most extreme PCA scores in the first dimension.
DimHeatmap (seu, dims = 1, cells = 500, balanced = TRUE)

# Calculation of k-NN and construction of SNN graph
seu <- FindNeighbors(object = seu, dims = 1:10, verbose = FALSE, k.param=20)

# Identification of clusters of cells
seu <- FindClusters(object = seu, verbose = FALSE, resolution=0.5)

# Dimensionality reduction on selected features
seu <- RunTSNE(object=seu, verbose=FALSE, dims.use = 1:10)
TSNEPlot(seu)

# create an object with all the relevant data that will be used to run the random forest
seu10 <- RunTSNE(seu,dims = 1:10, verbose = FALSE, n.components = 10)
tsne10 <- seu10@reductions$tsne@cell.embeddings
type <- seu$type
count <- seu$nCount_RNA
clust <- seu@meta.data$seurat_clusters
pca10 <- seu@reductions$pca@cell.embeddings
datRF <- data.frame(clust, tsne10, type, pca10, count)
datRF$clust <- as.factor(datRF$clust)
datRF$type <- as.factor(datRF$type)

# Random forest regression analysis
set.seed(101)
# Search for the optimal hypoparameters for the RF
t <- tuneRF(datRF[,-1],datRF[,1], stepFactor=1.5, plot=TRUE, ntreeTry=500, trace=TRUE, improve=0.001)

# Random forest algorithm was implemented to calculate the proximity of each pair of cells in a decision tree, with the aim to improve the low-dimensionality projection of the data.
rf <- randomForest(x=datRF[,-1], y=datRF[,1], data=datRF, ntree=3000, mtry=7, importance=TRUE, proximity=TRUE, do.trace=F)

result <- data.frame(cbind(datRF$clust, rf$predicted))
names(result) <- c("Original", "Predicted")

# Proximity calculated for each pair of cells based on proximity of cells in a decision tree.
prox <- data.frame(rf$proximity)
# Calculation of squared distance in a Euclidian space of dimension not greater than the number cells
d <- sqrt(1-prox)

# Create a low dimensional embedding of distances calculated by random forest
tsne_res <- Rtsne(d, dims = 2, perplexity = 95, theta = 0.0, check_duplicates = FALSE, pca = F, max_iter = 800, verbose = TRUE, is_distance = T)

# Extract matrix containing the new representations for the objects
z <- tsne_res$Y
# Initial exploratory analysis after RF
plot(z, col=result$Original, pch=19, main="Predicted")

# Create a new seurat object and overwrite cell embedding with the RF analysis output
# This step is performed to make the plot more esthetical
seu.new <- seu
seu.new@reductions$tsne@cell.embeddings <- z
colnames(seu.new@reductions$tsne@cell.embeddings) <- c('tSNE_1', 'tSNE_2')

# Add RF output of predicted cluster assignment as metadata 
new.clusters<- as.matrix(result$Predicted)
rownames(new.clusters) <- rownames(result)
seu.new <- AddMetaData(seu.new, new.clusters, col.name='new.clusters')
