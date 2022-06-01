folders=list()
c1=c("F","N","Q")
c2=c("1","2")

ii=1
for( i in c1)
for( j in  c2)

{ tmp=paste(i,j,sep="")
 folders[[ii]]=tmp
 ii=ii+1
}
folders=unlist(folders)


library(Seurat)
library(ggplot2)
datList=list()
for(i in 1:length(folders))datList[[i]] =Read10X(paste("cellrange_out",folders[i],sep="/"))
sceList=list()
for (i in 1:length(folders))
{
 tmp=CreateSeuratObject(counts = datList[[i]], project =folders[i],min.cells = 1, min.features = 200)
 tmp[["percent.mt"]]  = PercentageFeatureSet(tmp, pattern = "^MT-")
 sceList[[i]]=tmp
}
save(sceList,file="sceList.rda")

library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(patchwork)
dir.create('cluster6a1')

folders=list()
c1=c("F","N","Q")
c2=c("1","2")
ii=1
for( i in c1)
for( j in  c2)
{ tmp=paste(i,j,sep="")
	folders[[ii]]=tmp
	ii=ii+1
}
folders=unlist(folders)
folders=paste(folders,"_sc",sep="")

load("sceList.rda")

sce.f <- merge( sceList[[ 1 ]],y=c(sceList[[ 2 ]]), add.cell.ids = folders[1:2], project = "beeF")
sce.n <- merge( sceList[[ 3 ]],y=c(sceList[[ 4 ]]), add.cell.ids = folders[3:4], project = "beeN")
sce.q <- merge( sceList[[ 5 ]],y=c(sceList[[ 6 ]]), add.cell.ids = folders[5:6], project = "beeQ")

f.Q3Count <- quantile(sce.f$nCount_RNA)[4]
f.Q1Count <- quantile(sce.f$nCount_RNA)[2]
f.upper <- f.Q3Count + 1.5 * (f.Q3Count - f.Q1Count)
print(f.upper)
sce.f <- subset(sce.f, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt <15  & nCount_RNA< f.upper & nCount_RNA > 500)

n.Q3Count <- quantile(sce.n$nCount_RNA)[4]
n.Q1Count <- quantile(sce.n$nCount_RNA)[2]
n.upper <- n.Q3Count + 1.5 * (n.Q3Count - n.Q1Count)
print(n.upper)
sce.n <- subset(sce.n, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt <15  & nCount_RNA< n.upper & nCount_RNA > 500)

q.Q3Count <- quantile(sce.q$nCount_RNA)[4]
q.Q1Count <- quantile(sce.q$nCount_RNA)[2]
q.upper <- q.Q3Count + 1.5 * (q.Q3Count - q.Q1Count)
print(q.upper)
sce.q <- subset(sce.q, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt <15  & nCount_RNA< q.upper & nCount_RNA > 500)

sce.f$f <- "F"
sce.f <- NormalizeData(sce.f, verbose = FALSE)
sce.f <- FindVariableFeatures(sce.f, selection.method = "vst", nfeatures = 3000)

sce.n$n <- "N"
sce.n <- NormalizeData(sce.n, verbose = FALSE)
sce.n <- FindVariableFeatures(sce.n, selection.method = "vst", nfeatures = 3000)

sce.q$q <- "Q"
sce.q <- NormalizeData(sce.q, verbose = FALSE)
sce.q <- FindVariableFeatures(sce.q, selection.method = "vst", nfeatures = 3000)

bee.anchors <- FindIntegrationAnchors(object.list = list(sce.f, sce.n, sce.q), dims = 1:30)
sce.big <- IntegrateData(anchorset = bee.anchors, dims = 1:30)

# Run the standard workflow for visualization and clustering
sce.big <- ScaleData(sce.big, verbose = FALSE)
sce.big <- RunPCA(sce.big, npcs = 50, verbose = FALSE)

DimPlot(object = sce.big, reduction = "pca")
ElbowPlot(sce.big)
#t-SNE and Clustering
sce.big <- RunUMAP(sce.big, reduction = "pca", dims = 1:30)
sce.big <- RunTSNE(sce.big, reduction = "pca", dims = 1:30)
sce.big <- FindNeighbors(sce.big, reduction = "pca", dims = 1:30)
sce.big <- FindClusters(sce.big, resolution = 1)
sce.big <- FindVariableFeatures(sce.big, selection.method = "vst", nfeatures = 3000)

print(dim(sce.big))
print("after filter")

scRNA1=sce.big
print(dim(scRNA1))
print(table(scRNA1@meta.data$orig.ident))

metadata=scRNA1@meta.data
tt=metadata$orig.ident
tt=gsub(pattern = "1", replace ="", x = tt)
tt=gsub(pattern = "2", replace ="", x = tt)
tt1=tt
tt1[tt1=="F"]="FN"
tt1[tt1=="N"]="FN"

metadata=cbind(metadata,type=tt,type1=tt1)
scRNA1@meta.data=metadata

print(table(scRNA1@meta.data$type))
print(head(scRNA1@meta.data))

#tSNE
#group_by_cluster
plot1 = DimPlot(scRNA1, reduction = "tsne", label=T,size=1.2) 
ggsave("cluster6a1/0.1.tSNE.png", plot = plot1, width = 12, height = 8)
#group_by_sampl
plot2 = DimPlot(scRNA1, reduction = "tsne", group.by='orig.ident') 
ggsave("cluster6a1/0.2.tSNE_sample.png", plot = plot2, width = 8, height = 6)

#combinate
plotc <- plot1+plot2
ggsave("cluster6a1/0.3.tSNE_cluster_sample.png", plot = plotc, width = 14, height = 4)

#UMAP
embed_umap <- Embeddings(scRNA1, 'umap')
write.csv(embed_umap,'cluster6a1/embed_umap.csv') 
#group_by_cluster2
plot3 = DimPlot(scRNA1, reduction = "umap", label=T) 
ggsave("cluster6a1/1.1.UMAP.png", plot = plot3, width = 12, height = 8)
#group_by_sample
plot4 = DimPlot(scRNA1, reduction = "umap", group.by='orig.ident')
ggsave("cluster6a1/1.2.UMAP_sample.png", plot = plot4, width = 8, height = 6)
#combinate
plotc <- plot3+plot4
ggsave("cluster6a1/1.3.UMAP_cluster_sample.png", plot = plotc, width = 14, height = 4)

plotc <- plot2+plot4+ plot_layout(guides = 'collect')
ggsave("cluster6a1/2.0.tSNE_UMAP.png", plot = plotc, width = 15, height = 5)

save(scRNA1,file="cluster6a1/scRNA1_final.rda")
# Identify the 10 most highly variable genes

plot3 = DimPlot(scRNA1, reduction = "umap", label=T) 
ggsave("cluster6a1/3.UMAP.png", plot = plot3, width = 12, height = 8)
#group_by_sample
plot3 = DimPlot(scRNA1, reduction = "umap", label=T,split.by ="type",size=0.8)
ggsave("cluster6a1/3.UMAP_cellType_fq.png", plot = plot3, width = 15, height = 5)

plot3 = DimPlot(scRNA1, reduction = "umap", label=T,split.by ="type1",size=0.8)
ggsave("cluster6a1/3.UMAP_cellType1_fq.png", plot = plot3, width = 15, height = 8)

plot3 = DimPlot(scRNA1, reduction = "tsne", label=T)
ggsave("cluster6a1/3.tsne.png", plot = plot3, width = 12, height = 8)
#group_by_sample

plot4 = DimPlot(scRNA1, reduction = "tsne", split.by='type',size=0.8)
ggsave("cluster6a1/3.tsne_cellType_fq.png", plot = plot4, width = 15, height = 6)
plot4 = DimPlot(scRNA1, reduction = "tsne", split.by='type1',size=0.8)
ggsave("cluster6a1/3.tsne_cellType1_fq.png", plot = plot4, width = 15, height = 8)

top20 <- head(VariableFeatures(scRNA1), 20)
pdf("cluster6a1/0.0.top20genes.pdf",width=8,height=6)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scRNA1)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
dev.off()

pdf("cluster6a1/0.2.visDimPCs_.pdf",width=5,height=10)
VizDimLoadings(scRNA1)
dev.off()

print(scRNA1[["pca"]], dims = 1:5, nfeatures = 5)
pdf("cluster6a1/0.3.visDimPCs_.pdf",width=12,height=12)
DimHeatmap(scRNA1, dims = 1:10, cells = 500, balanced = TRUE)
dev.off()


scRNA1<-BuildClusterTree(scRNA1)
Tool(object = scRNA1, slot = 'BuildClusterTree')
pdf("cluster6a1/0.4.clusterTree.pdf",width=9,height=6)
PlotClusterTree(scRNA1)
dev.off()

#find all markers of cluster 
sce.markers_bim <- FindAllMarkers(scRNA1, only.pos = TRUE, min.pct = 0.1, logfc.threshold = .5,return.thresh=0.01,test.use="bimod")
save(sce.markers_bim,file="cluster6a1/sce.marker_bim.rda")

