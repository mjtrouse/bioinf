#Run seurat pipeline on CF data

library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
#install.packages("R.utils")
library(R.utils)

#unzip gz file
#CF_DATA<-Matrix::readMM("CF_Data/matrix.mtx.gz")

#new way to load matrix as mtx file from ucsc instructions(did not work)
#mat<-fread("CF_Data/matrix.mtx.gz")
#meta <- read.table("CF_Data/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
#genes = mat[,1][[1]]
#genes = gsub(".+[|]", "", genes)
#mat = data.frame(mat[,-1], row.names=genes)
#cf.1 <- CreateSeuratObject(counts = mat, project = "CFsputum", meta.data=meta)

#try seurat way
mat2<-Seurat::ReadMtx("CF_Data/matrix.mtx.gz",cells="CF_Data/barcodes.tsv.gz",
                      features="CF_Data/genes.tsv.gz", feature.column = 1)

meta = read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
# Initialize the Seurat object with the raw (non-normalized data).
cf.1 <- CreateSeuratObject(counts = mat2, project = "cf", min.cells = 3, min.features = 200,
                           meta.data=meta)
cf.1

#percentage of mitochondrial genes(high perc. = low quality cells)
cf.1[["percent.mt"]] <- PercentageFeatureSet(cf.1, pattern = "^MT-")

#subset data with unique features between 200 and 2500
#and less than 5 percent mitochondrial genes
cf.1 <- subset(cf.1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#normalize data
cf.1 <- NormalizeData(cf.1, normalization.method = "LogNormalize", scale.factor = 10000)

#feature selection(find highly variable features)
cf.1 <- FindVariableFeatures(cf.1, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cf.1), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(cf.1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scale to mean 0 and variance 1
all.genes <- rownames(cf.1)
cf.1 <- ScaleData(cf.1, features = all.genes)

#PCA dimensionaily reduction on subset of scaled variable features
cf.1 <- RunPCA(cf.1, features = VariableFeatures(object = cf.1))

#visualize reduced dimensions(only 1st and 2nd dimensions)
VizDimLoadings(cf.1, dims = 1:2, reduction = "pca")

#dimension heat map(first dimension- primary sources for heterogen.)
#cells and features ordered by PCA scores
DimHeatmap(cf.1, dims = 1, cells = 500, balanced = TRUE)

#Determine "dimensionality" of dataset
#how many PC to include?
#random permutation of 1% of data and rerun PCA
#identify PC's that enrich low p-value features
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
cf.1 <- JackStraw(cf.1, dims=50, num.replicate = 100)
cf.1 <- ScoreJackStraw(cf.1, dims = 1:50)

#plot JackStraw compare dist of p-values for each PC 
#to uniform dist
JackStrawPlot(cf.1, dims=1:50) #drops off around 26 but some after sign.


#look at p-values for all dims
#JS(object = cf.1[['pca']], slot = 'empirical')

#alternatively, look at elbow plot
#where is the turn in values? 
ElbowPlot(cf.1, ndims = 25)  #drop between 10-15 here?

#Identifying true dimensionality uncertain
#best to overestimate 
#best to repeat downstream an. with different numbers of PC's

#######Cluster the cells
#graph-based approach (K-nearest neighbor graph, partition)
#KNN based on euclidean distance in PCA space
#refine edge weights between 2 cells based on shared overlap in local hood
cf.1 <- FindNeighbors(cf.1, dims = 1:20) 

#cluster with modularity optimization techniques(Louvain or SLM)
#iteratively group cells together
#resolution sets "granularity"
#higher res=more clusters (0.4-1.2 good for 3K cells)
cf.1 <- FindClusters(cf.1, resolution = 1)

#find clusters using 
#Idents(cf.1)

#Look at cluster ID's of first 5 cells
head(Idents(cf.1), 5) #different results here than tut bc of smaller PCA repeat

##########Run non-linear dimensional reduction(UMAP/tSNE)
#try to place similar cells together in a low-dimensional space
#suggested use fo same PC's as clustering analysis
# If you haven't installed UMAP, you can do so via 
#reticulate::py_install(packages = 'umap-learn')

#UMAP-Unif.ManifoldAppx&projection
cf.1 <- RunUMAP(cf.1, dims = 1:12)

#visualize UMAP
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(cf.1, reduction = "umap", label=T)

#Run tSNE and plot
cf.1<- RunTSNE(cf.1, reduction="pca", dims=1:12)

TSNEPlot(cf.1, reduction="tsne", label=T)

unique(Idents(cf.1))

#save object up to now
#saveRDS(pbmc, file = "pbmc_tutorial.rds")

################################################
###############################################
#Find cluster biomarkers
#(differently expressed features)

#ID's "positive and negative" markers
#can automate all clusters but
#can test groups against all others or 
#against another group

#find all markers of cluster 2
cluster1.markers <- FindMarkers(cf.1, ident.1=1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 3 from clusters 0 and 2
cluster3.markers <- FindMarkers(cf.1, ident.1 = 3, ident.2 = c(0, 2), min.pct = 0.25)
head(cluster3.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
cf1.markers <- FindAllMarkers(cf.1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cf1.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#ROC returns classification power for indiv.markers(0-1)
cluster0.markers <- FindMarkers(cf.1, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster0.markers

#violin plot for raw counts of specific features 
# you can plot raw counts as well
VlnPlot(cf.1, features = c("H3F3B"), slot = "counts", log = TRUE)

#visualize features on a tSNE(t-dist stoch neigh. embedd.) plot
FeaturePlot(cf.1, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

#heatmap for cells and features grouped by cluster
cf1.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(cf.1, features = top10$gene) + NoLegend()

