#REPO code for CF_seurat analysis

library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(R.utils)
library(ggplot2)


set.seed(412)

##Initial pre-processing and analysis based on paper methods
#load gene matrix
mat<-Seurat::ReadMtx("CF_Data/matrix.mtx.gz",cells="CF_Data/barcodes.tsv.gz",
                      features="CF_Data/genes.tsv.gz", feature.column = 1)
#load metadata
meta = read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)

# Initialize the Seurat object with the raw (non-normalized data).
#New seurat object for paper analysis
cf.2 <- CreateSeuratObject(counts = mat, project = "cf", min.cells = 3, min.features = 200,
                           meta.data=meta)

#get percent mitochondrial genes
cf.2[["percent.mt"]] <- PercentageFeatureSet(cf.2, pattern = "^MT-")

#subset data with unique features between less than 2250 & less than 20% mitochon.
cf.2 <- subset(cf.2, subset = nFeature_RNA>200 & nFeature_RNA < 2250 & percent.mt <20)

#normalize data
cf.2 <- NormalizeData(cf.2, normalization.method = "LogNormalize", scale.factor = 10000)

#feature selection(find highly variable features)
cf.2 <- FindVariableFeatures(cf.2, selection.method = "vst", nfeatures = 3000)

#scale to mean 0 and variance 1
all.genes <- rownames(cf.2)
cf.2<-ScaleData(cf.2, features=all.genes, verbose=F)

#regress out percentage number of mitochondrial genes
cf.2 <-ScaleData(cf.2, vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt"), do.scale=F, do.center=F, verbose=F)

#PCA dimensionaily reduction on subset of scaled variable features
cf.2 <- RunPCA(cf.2, features = VariableFeatures(object = cf.2))

#visualize reduced dimensions
VizDimLoadings(cf.2, dims = 1:10, reduction = "pca")

#dimension heat map(first dimension- primary sources for heterogen.)
#cells and features ordered by PCA scores
DimHeatmap(cf.2, dims = 1:20, cells = 500, balanced = TRUE)

#Determine "dimensionality" of dataset (paper used 12)
ElbowPlot(cf.2, ndims=20, reduction="pca")

#######Cluster the cells
#nearest neighbors
cf.2 <- FindNeighbors(cf.2, dims = 1:12)

#cluster with Louvain (no resolution size outlined in methods)
cf.2 <- FindClusters(cf.2, resolution = 0.5)

##########Run non-linear dimensional reduction(UMAP/tSNE)
#try to place similar cells together in a low-dimensional space
#suggested use fo same PC's as clustering analysis

#UMAP-Unif.ManifoldAppx&projection
cf.2 <- RunUMAP(cf.2, dims = 1:12)

cf.2<-RunTSNE(cf.2, reduction="pca", dims = 1:12)

#visualize UMAP/tsne
DimPlot(cf.2, reduction = "umap", label=T)

TSNEPlot(cf.2, reduction="tsne", label=T)

#########Get cluster markers
#find only positive markers with over 25% of cells in cluster expressing
cf2.markers <- FindAllMarkers(cf.2, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.25)
CF2.MARKERS<-as.data.frame(cf2.markers %>%
                                  group_by(cluster) %>%
                                  slice_max(n = 20, order_by = avg_log2FC)
)

################################################################################
################################################################################
#Subsetting CF.2 based on markers for Monocytes(Mono), Alveolar Macrophages(AM),
# & Monocyte derived Macrophages(MoMO)

#look for markers, remove clusters, redo downstream analysis, repeat

###############################################
##Subset 1

#CF.2 MARKERS dataframe used to manually label clusters based on highly
#expressed genes from literature

#common macrophage markers
VlnPlot(cf.2, features=c("FNIP2", "ALCAM"))

#remove neutrophil clusters and epitheleal cells
cf.2.small.1<-subset(cf.2, idents=c(9,7,5,4,3))   #3983 samples, 35451 features

#redo downstream with cf2.subset1
cf.2.small.1 <- FindVariableFeatures(cf.2.small.1, selection.method = "vst", nfeatures = 3000)

#scale and center data, regress out variables
cf.2.small.1<-ScaleData(cf.2.small.1, 
                        features=VariableFeatures(cf.2.small.1), verbose=F)
cf.2.small.1 <-ScaleData(cf.2.small.1, 
                         vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt"), 
                         do.scale=F, do.center=F, verbose=F)

#subset1 PCA dimensionaily reduction on subset of scaled variable features
cf.2.small.1 <- RunPCA(cf.2.small.1,
                       features = VariableFeatures(object = cf.2.small.1))

#visualize pc's
DimHeatmap(cf.2.small.1, dims = 1:20,cells=500, balanced=T, reduction = "pca")

#Determine "dimensionality" of subset1
ElbowPlot(cf.2.small.1, ndims=20, reduction="pca")

#######Cluster the cells
#nearest neighbors
cf.2.small.1 <- FindNeighbors(cf.2.small.1, dims = 1:12)

#cluster with Louvain (no resolution size outlined in methods)
cf.2.small.1 <- FindClusters(cf.2.small.1, resolution = 0.5)

##########Run non-linear dimensional reduction
#subset1 tsne
cf.2.small.1<-RunTSNE(cf.2.small.1, reduction="pca", dims=1:12)

#visualize
DimPlot(cf.2.small.1, reduction="tsne", label=T)

#find clusters with monocyte/macro upreg
VlnPlot(cf.2.small.1, features=c("SERPINB9", "TIMP1", "JARID2", #monocyte top
                                 "APOE","STMN1", "ANXA1",  #AM top
                                 "FNIP2", "ALCAM", "GPR183"      #M top
))

#compare cluster 2 expressions to 5 and 7
cluster2<-FindMarkers(cf.2.small.1, ident.1 =2, 
                      ident.2= c(5,7),logfc.threshold = 0.25,
                      test.use = "roc", only.pos = TRUE)

####################################################################
#SUBSET 2

#remove clusters with low expression of above 
cf.2.small.2<-subset(cf.2.small.1, idents=c(2,4,5,6,7,8)) 
#35451 features, 1998 samples

#redo downstream with cf2.subset1
cf.2.small.2 <- FindVariableFeatures(cf.2.small.2, selection.method = "vst", nfeatures = 3000)

#scale to mean 0 and variance 1 & regress out percentage number of mitochondrial genes
cf.2.small.2<-ScaleData(cf.2.small.2, 
                        features=VariableFeatures(cf.2.small.2), verbose=F)
cf.2.small.2 <-ScaleData(cf.2.small.2, 
                         vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt"), 
                         do.scale=F, do.center=F, verbose=F)


#subset1 PCA dimensionaily reduction on subset 2 of scaled variable features
cf.2.small.2 <- RunPCA(cf.2.small.2, 
                       features = VariableFeatures(object = cf.2.small.2))

#Determine "dimensionality" of subset2
ElbowPlot(cf.2.small.2, ndims=20, reduction="pca")

#######Cluster the cells (subset2)
#nearest neighbors 
cf.2.small.2 <- FindNeighbors(cf.2.small.2, dims = 1:15) #higher dim number

#cluster with Louvain 
cf.2.small.2 <- FindClusters(cf.2.small.2, resolution = 0.5)

##########Run non-linear dimensional reduction
#subset2 tsne
cf.2.small.2<-RunTSNE(cf.2.small.2, reduction="pca", dims=1:15)

#visualize
DimPlot(cf.2.small.2, reduction="tsne", label=T)

#Get features from this subset
cf2.small.2.markers <- FindAllMarkers(cf.2.small.2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CF2.SMALL.2.MARKERS<-as.data.frame(cf2.small.2.markers %>%
                                     group_by(cluster) %>%
                                     slice_max(n = 10, order_by = avg_log2FC)
)

write.csv(CF2.SMALL.2.MARKERS, "CF_Data/small.2.markers(2).csv")
#########################################################################
##subset 3

#find clusters with monocyte/macro upreg
VlnPlot(cf.2.small.2, features=c("SERPINB9", "TIMP1", "JARID2", #monocyte top
                                 "APOE","STMN1", "ANXA1",  #AM top
                                 "FNIP2", "ALCAM", "GPR183"      #M top
))

#remove clusters based on top 10 markers(cf.2.small.2 df) and vln plot
cf.2.small.3<-subset(cf.2.small.2, idents=c(0,1,2,3,7))
#1698 samples

#redo downstream with cf2.subset3
cf.2.small.3 <- FindVariableFeatures(cf.2.small.3, selection.method = "vst", nfeatures = 3000)

#scale to mean 0 and variance 1 & regress out percentage number of mitochondrial genes
cf.2.small.3<-ScaleData(cf.2.small.3, features=VariableFeatures(cf.2.small.3), verbose=F)
cf.2.small.3 <-ScaleData(cf.2.small.3, 
                         vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt"), 
                         do.scale=F, do.center=F, verbose=F)

#subset1 PCA dimensionaily reduction on subset 3 of scaled variable features
cf.2.small.3 <- RunPCA(cf.2.small.3, 
                       features = VariableFeatures(object = cf.2.small.3))

#visualize pc dimensions
DimHeatmap(cf.2.small.3, dims=1:20, reduction="pca", cells=500,balanced = T)

#Determine "dimensionality" of subset1
ElbowPlot(cf.2.small.3, ndims=20, reduction="pca")

#######Cluster the cells (subset2)
#nearest neighbors 
cf.2.small.3 <- FindNeighbors(cf.2.small.3, dims = 1:15) #higher dimension number

#cluster with Louvain
cf.2.small.3 <- FindClusters(cf.2.small.3, resolution = 0.5)

##########Run non-linear dimensional reduction(UMAP/tSNE)
#subset1 tsne
cf.2.small.3<-RunTSNE(cf.2.small.3, reduction="pca", dims=1:15)

#visualize
DimPlot(cf.2.small.3, reduction="tsne", label=T)

#Get features from this subset
cf2.small.3.markers <- FindAllMarkers(cf.2.small.3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CF2.SMALL.3.MARKERS<-as.data.frame(cf2.small.3.markers %>%
                                     group_by(cluster) %>%
                                     slice_max(n = 10, order_by = avg_log2FC)
)

write.csv(CF2.SMALL.3.MARKERS, "CF_Data/cf.2.small.3.markers.csv")

small.3.topgenes<-cf2.small.3.markers %>%
  top_n(n = 10, wt = -(p_val_adj)) 

########################
#Rename clusters for subset 3(based on literature & top 10 genes)
new.cluster.ids <- c("Mono1", "Mono2", "AM1", "Mono3","AM2", "Mult")
names(new.cluster.ids) <- levels(cf.2.small.3)
cf.2.small.3 <- RenameIdents(cf.2.small.3, new.cluster.ids)

#replot tsne with new names
DimPlot(cf.2.small.3, reduction="tsne", label=T)

#compare cf vs. control subset 3
VlnPlot(cf.2.small.3, 
        features = small.3.topgenes$gene , 
        split.by = "disease.ident",split.plot=T,log = T)+
  theme(legend.position = "right")

###########################################################
#look at known necroptosis and senescence markers in subset 3
#necrop: RIP1, RIP2, MLKL
#senesc: CDKN2A, CDKN1A,TP53

VlnPlot(cf.2.small.3, 
        features=c("RIPK1", "RIPK2", "MLKL", "CDKN2A", "CDKN1A", "TP53"), 
        log=T)
#three Monocyte clusters have markers present

#compare ripk2 and cdkn1a cf/control that are present
VlnPlot(cf.2.small.3, features=c("RIPK2", "CDKN1A"), split.by="disease.ident",
        split.plot = T, log=T)+
  theme(legend.position = "right")

#############################################################
#subset 4
#remove cluster with multiple cell types

cf.2.small.4<-subset(cf.2.small.3, idents=c("Mono1", "Mono2", "Mono3", 
                                            "AM1","AM2")) #1634 samples

#redo downstream with cf2.subset4
cf.2.small.4<- FindVariableFeatures(cf.2.small.4, selection.method = "vst", nfeatures = 3000)

#scale to mean 0 and variance 1 & regress out percentage number of mitochondrial genes

cf.2.small.4<-ScaleData(cf.2.small.4, 
                        features=VariableFeatures(cf.2.small.4), verbose=F)
cf.2.small.4 <-ScaleData(cf.2.small.4, 
                         vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt"), 
                         do.scale=F, do.center=F, verbose=F)


#subset4 PCA dimensionaily reduction on subset 4 of scaled variable features
cf.2.small.4 <- RunPCA(cf.2.small.4, features = VariableFeatures(object = cf.2.small.4))

DimHeatmap(cf.2.small.4, dims=1:20, cells=500, balanced=T, reduction = "pca")

#Determine "dimensionality" of subset4
ElbowPlot(cf.2.small.4, ndims=20, reduction="pca")

#######Cluster the cells (subset4)
#nearest neighbors 
cf.2.small.4 <- FindNeighbors(cf.2.small.4, dims = 1:12) 

#cluster with Louvain
cf.2.small.4 <- FindClusters(cf.2.small.4, resolution = 0.8)

##########Run non-linear dimensional reduction(UMAP/tSNE)
#subset4 tsne
cf.2.small.4<-RunTSNE(cf.2.small.4, reduction="pca", dims=1:12)

#visualize
DimPlot(cf.2.small.4, reduction="tsne", label=T)

DimPlot(cf.2.small.4, reduction="pca", label=T)

DimHeatmap(cf.2.small.4, reduction = "pca", dims = 1:12,cells = 500, balanced=T)

#Get features from this subset
cf2.small.4.markers <- FindAllMarkers(cf.2.small.4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CF2.SMALL.4.MARKERS<-as.data.frame(cf2.small.4.markers %>%
                                     group_by(cluster) %>%
                                     slice_max(n = 20, order_by = avg_log2FC)
)

write.csv(CF2.SMALL.4.MARKERS, "CF_Data/cf.2.small.4.markers.csv")

small.4.topgenes<-cf2.small.4.markers %>%
  top_n(n = 10, wt = -(p_val_adj)) 

#compare cf vs. control for top 10 genes
VlnPlot(cf.2.small.4, 
        features = small.4.topgenes$gene , 
        split.by = "disease.ident",split.plot=T,log = T)+
  theme(legend.position = "right")

########################
#Rename clusters for subset 4(literature and top 20 DF)
new.cluster.ids <- c("Mono1", "AM1", "Mono2/MoMO", "Mono3","Mono4", "AM2")
names(new.cluster.ids) <- levels(cf.2.small.4)
cf.2.small.4 <- RenameIdents(cf.2.small.4, new.cluster.ids)

#replot tsne with new names
DimPlot(cf.2.small.4, reduction="tsne",shape.by = "disease.ident",label=T)


#look at known necroptosis and senescence markers in subset 4
#necrop: RIP1, RIP2, MLKL
#senesc: CDKN2A, CDKN1A,TP53

VlnPlot(cf.2.small.4, 
        features=c("RIPK1", "RIPK2", "MLKL", "CDKN2A", "CDKN1A", "TP53"), 
        log=T)
#all monocytes have expression upreg

#compare ripk2 and cdkn1a cf/control that are present
VlnPlot(cf.2.small.4, features=c("RIPK2", "CDKN1A"), split.by="disease.ident",
        split.plot = T, log=T)+
  theme(legend.position = "right")

saveRDS(cf.2.small.4, file="CF_Data/cf.2.small.4.rds")

######################################################
#try different downstream dim=8 (no change)
{
  #nearest neighbors 
  cf.2.small.4.8 <- FindNeighbors(cf.2.small.4, dims = 1:8) 
  
  #cluster with Louvain
  cf.2.small.4.8 <- FindClusters(cf.2.small.4.8, resolution = 0.8)
  
  #subset4 tsne
  cf.2.small.4.8<-RunTSNE(cf.2.small.4.8, reduction="pca", dims=1:8)
  #visualize
  DimPlot(cf.2.small.4.8, reduction="tsne", label=T)
  
  DimHeatmap(cf.2.small.4.8, reduction = "pca", dims = 1:8,
             cells = 500, balanced=T)
}
#######################################################
###Markers from Yong
#no MHCII, cd64, CD11b, ly6c, CD169
#look at all cf.2
VlnPlot(cf.2, 
        features=c("RIPK3" ,"ADGRE1", "CD68", 
                   "MERTK", "CD14"), log=T)

#look at subset
VlnPlot(cf.2.small.4, 
        features=c("RIPK3" ,"ADGRE1", "CD68", 
                   "MERTK", "CD14"), log=T)

