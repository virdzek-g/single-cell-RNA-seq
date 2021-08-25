# Construction of single cell trajectories using Monocle 2 R package on Seurat 3 object
# Analysis done on raw read counts from single cell RNA seq experiment.
# This analysis is adapted from the Monocle 2 vignette

# raw counts data converted into sparse matrix object
## rows = genes, columns = cells after filtration
data <- as(as.matrix(counts), 'sparseMatrix')

# object containing phenotype data corresponding to each cell
## rows = colnames(data)
pd <- new('AnnotatedDataFrame', data = seu@meta.data)

# object containing features data corresponding to each gene
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
# rows = rownames(data)
fd <- new('AnnotatedDataFrame', data = fData)

# create object with the counts data and metadata
mon_obj <- newCellDataSet(data, phenoData = pd,featureData = fd,expressionFamily = negbinomial.size())

# pre-calculate data for differential expression analysis
mon_obj <- estimateSizeFactors(mon_obj)
mon_obj <- estimateDispersions(mon_obj)

# Counts how many cells express each feature in a CellDataSet object that are detectably expressed above a minimum threshold. Counts the number of genes above this threshold detectable in each cell.
mon_obj <- detectGenes(mon_obj, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(mon_obj),
    num_cells_expressed >= 10))
# expressed_genes matrix holds the identifiers for genes expressed in at least 50 cells of the data set.

# Log-transform each value in the expression matrix.
log_m <- log(exprs(mon_obj[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(log_m))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "density", data = melted_dens_df) +
stat_function(fun = dnorm, size = 0.5, color = 'red') +
xlab("Standardized log(counts)") +
ylab("Density")

# Clustering cells without marker genes. Groups cells together according to global expression profile.
## Filter genes based on average expression level, and select genes that are variable across cells.
disp_table <- dispersionTable(mon_obj)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)

# Set list of gene ids to be used for ordering in the monocle object
mon_obj <- setOrderingFilter(mon_obj, unsup_clustering_genes$gene_id)
plot_ordering_genes(mon_obj)

# Cell clustering
mon_obj <- reduceDimension(mon_obj, max_components = 2, num_dim = 10,
                reduction_method = 'tSNE', verbose = T)
mon_obj <- clusterCells(mon_obj, num_clusters = 2)

# Compare groups of cells to find differentially expressed genes, controlling for treatments (new_class)
## Find all genes that are differentially expressed between clusters identified in the upstream analysis done with Seurat 3 package
## Output is a list of differentially expressed genes
diff_test_res <- differentialGeneTest(mon_obj[expressed_genes,],
              fullModelFormulaStr = "~new_class")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.1))
length(ordering_genes)

# set the gene list into the monocle object
mon_obj <- setOrderingFilter(mon_obj, ordering_genes)
plot_ordering_genes(mon_obj)

# Reduce data dimensionality
mon_obj <- reduceDimension(mon_obj, max_components = 2,
    method = 'DDRTree')

# Order cells in pseudotime along a trajectory    
mon_obj <- orderCells(mon_obj)

# Plot trajectory and colour cells according to different clustering resolutions and pseudotime
p1 <- plot_cell_trajectory(mon_obj, color_by = "SCT_snn_res.1")
p2 <- plot_cell_trajectory(mon_obj, color_by = "SCT_snn_res.0.8")
p3 <- plot_cell_trajectory(mon_obj, color_by = "Pseudotime")
p1 | p2 | p3

