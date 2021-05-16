# Analysis part 4

# Visualisation of scRNA-seq data
#
# This analysis is done on filtered, normalised and clustered data after dimensionality reduction


# Load the libraries
library(Seurat)
library(dplyr)

TSNEPlot(seu.new, group.by='new.clusters')
TSNEPlot(seu.new, group.by='type')

# features argument is a vactor of genes of interest
DoHeatmap(seu.new, group.by='new.clusters', features=features, assay = 'SCT')+
 scale_fill_gradientn(colors = c("#02114E", "#FFFFFF", "#CE060F"))

DotPlot(seu.new, features=features, dot.scale = 8, group.by='final.clusters') + 
RotatedAxis() + scale_color_gradient2(low="#052A6E", mid="#FFFFFF", high = "#A60000")

VlnPlot(seu.new, features=features, n=9)), idents=c(3,7), ncol=3, cols=c('#34D54E','#FF2FA6'))

RidgePlot(seu.new, features=feautures, group.by='new.clusters')









