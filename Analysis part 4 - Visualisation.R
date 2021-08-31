# Analysis part 4

# Visualisation of scRNA-seq data
#
# This analysis is done on filtered, normalised and clustered data after dimensionality reduction


# Load the libraries
library(Seurat)
library(dplyr)

# seu is the seurat object containing normalised data.

# Visualise gene expression on a dimensional reduction plot
FeaturePlot(seu, features=c('Foxp3', 'Cd8a', 'Cd4', 'Pdcd1', 'Tcf7', 'Sell', 'Il7r', 'Rpl30', 'Rps2', 'Cxcr3', 'Gzmb', 'Gzma', 'B2m'), reduction='tsne', cols = c("#052A6E","#4577D4", "#FFFFFF","#FF4040" ,"#A60000"))

# Visualisation of gene expression across clusters or treatments.The size of the dot encodes the percentage of cells within a class, while the color encodes the average expression level across all cells within a class.
DotPlot(seu, dot.scale = 8, features=c('Ezh2', 'Pcna','Top2a','Mki67', 
'Eomes','Id2','Id3','Tbx21',
'Gzma', 'Gzmb', 'Cd44','Klrg1','S1pr5','Cx3cr1',
'Sell','Il7r', 'Tcf7'
)) + RotatedAxis() + scale_color_gradient2(low="#052A6E", mid="#FFFFFF", high = "#A60000")

##################################################
##################################################

# Adding gene modules of relevant signalling pathways
# Calculate the average expression levels of each group of genes on single cell level, subtracted by the aggregated expression of control feature sets. All analyzed features are binned based on averaged expression, and the control features are randomly selected from each bin.

# adding modules AP_1, rb, mm, and act for different signalling pathways
AP_1 <- list(c('Klf6','Cd69','Klf2','Id2','Zfp36','Fos','Fosb','Jund','Junb','Jun','Ccl5','Ccl4'))
seu <- AddModuleScore(object=seu, features=AP_1, ctrl=5, name = "AP_1")

rb <- list(c('Rps8','Rplp0','Rpl15','Rpl28','Rpl7','Rpl11','Rps7','Rps27a','Rps20','Rps2','Rps5','Rpl9','Rpl9','Rpl35'))
seu.new <- AddModuleScore(object=seu.new, features=rb, ctrl=5, name = "rb")

mm <- list(c('Il7r','Tcf7','Sell','Bcl2', 'Cxcr3','Id3','Ccr7'))
seu.new <- AddModuleScore(object=seu.new, features=mm, ctrl=3, name = "mm")

act <- list(c('Ccl4','Itga4','Ahnak','S100a4','Klrk1','Id2','S100a10','Il18r1','Lgals1','Ccl5','Ccr2','Klrc1','S100a11','Gzma','Il18rap','Srgn','Gzmb','Ifng','Hmgb2'))
seu.new <- AddModuleScore(object=seu.new, features=act, ctrl=5, name = "act")


# Visualising modules
names(x = seu[[]])
p1 <- DimPlot(seu, group.by='new.clusters', reduction='tsne')
p2 <- FeaturePlot(seu.new, features='AP_11', cols = c("#052A6E","#4577D4", "#FF4040" ,"#A60000"), reduction='tsne')
p3 <- FeaturePlot(seu.new, features='rb1', cols = c("#052A6E","#4577D4", "#FF4040" ,"#A60000"), reduction='tsne')
p4 <- FeaturePlot(seu.new, features='act1', cols = c("#052A6E","#4577D4", "#FF4040" ,"#A60000"), reduction='tsne')
p5 <- FeaturePlot(seu.new, features='mm1', cols = c("#052A6E","#4577D4", "#FF4040" ,"#A60000"), reduction='tsne')
p1 | p2 | p3 | p4 | p5

##################################################
##################################################

# Heatmap on cells grouped in each individual cluster done on scRNA-seq data
library(Seurat)
library(gplots)

# seu is the seurat object containing normalised data in @assay SCT.
# Function AverageExpression() calculates average expression of each gene in each identity class
avgexp = AverageExpression(seu, return.seurat = T, add.ident = 'seurat_clusters')
sct <- avgexp@assays$SCT

# genes of interest
features <- c('Prf1','Gzmb','Gzma','Cd244','Klrc3','Klra9','Klrb1c','Klrc1','Klre1','Klrg1','Sell','Il7r','Cxcr5','Cxcr3','Slamf6','Ccr7','Itgax','Ccl4','Ccl9','Cxcr1','S1pr5','Eomes','Myc','Id3','Tcf7','Runx1','Bhlhe40','Zeb2','Prdm1')

# extraction of expression data for gene of interest
data <- sct[rownames(sct) %in% features,]
data <- as.matrix(data)

# heatmap
heatmap.2(data, scale='column',Rowv=F,trace="none",col =  colorRampPalette(c("darkblue","white","darkred"))(100))

##################################################
##################################################

# GSEA - results of this analysis can be used as an input for GSEA program

markers <- FindMarkers(seu,group.by='class', ident.1 = 'MT', ident.2='WT', logfc.threshold=0.0001)

dif_ex_clust <- data.frame(genes= rownames(markers), logFC = markers$avg_logFC)
rownames(dif_ex_clust) <- toupper(dif_ex_clust$genes)
write.csv(dif_ex_clust, file='dif_ex_clust.csv')

