library(dplyr)
library(Seurat)
library(garnett)
library(monocle)
library(splitstackshape)
library(celldex)
library(scater)
library(scRNAseq)
library(Matrix)

#classification packages
library(infercnv)
library(SingleR)
library(scPred)

#VITAL ORGAN SINGLE CELL DATA

###############################################################
                           #HEART#
###############################################################
setwd('D:/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/HEART/')
setwd('~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/HEART/')

counts=read.csv('GSE109816_normal_heart_umi_matrix.csv')
row.names(counts)=counts$X
counts$X=NULL
anno=read.table('GSE109816_normal_heart_cell_cluster_info.txt',header=T)
table(anno$CellType)

anno[anno$CellType=='CM','cell.name']='Cardiomyocytes'
anno[anno$CellType=='EC','cell.name']='Endothelial'
anno[anno$CellType=='FB','cell.name']='Fibroblasts'
anno[anno$CellType=='MP','cell.name']='Macrophages'
anno[anno$CellType=='SMC','cell.name']="Smooth muscle"

row.names(anno)=anno$ID

healthy_heart <- CreateSeuratObject(counts = counts, project = "hcc10k", min.cells = 3, min.features = 10,meta.data = anno)
healthy_heart[["percent.mt"]] <- PercentageFeatureSet(healthy_heart, pattern = "^MT-")
VlnPlot(healthy_heart, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 0)
healthy_heart <- subset(healthy_heart, subset = nFeature_RNA > 0 & percent.mt < 50 )
healthy_heart <- NormalizeData(healthy_heart, normalization.method = "LogNormalize", scale.factor = 10000)
healthy_heart <- FindVariableFeatures(healthy_heart, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(healthy_heart)
healthy_heart <- ScaleData(healthy_heart, features = all.genes)
healthy_heart <- FindVariableFeatures(healthy_heart)
healthy_heart <- RunPCA(healthy_heart, features = VariableFeatures(object = healthy_heart))

VizDimLoadings(healthy_heart, dims = 1:2, reduction = "pca")
healthy_heart <- FindNeighbors(healthy_heart, dims = 1:30,verbose=TRUE)

set.seed(123)
healthy_heart <- FindClusters(healthy_heart, resolution = 0.8)
head(Idents(healthy_heart), 5)
healthy_heart <- RunUMAP(healthy_heart, dims = 1:30)
ElbowPlot(healthy_heart,reduction='umap')

pdf('heart_vital.pdf')
{
  DimPlot(healthy_heart, reduction = "umap",label=T)+ggtitle('UMAP of healthy heart tissues')
  DimPlot(healthy_heart, reduction = "umap",label=T,group.by = 'cell.name')+ggtitle('UMAP of healthy lung tissues')
  
  
  favorable_markers=c('GPC3','PIGY','TMEM150B','MUC13','TREM2','SLC51B','MELK','DTL','TGM3','GPR158')
  favorable_markers=favorable_markers[favorable_markers%in%row.names(counts)]
  
  FeaturePlot(healthy_heart, features = favorable_markers)
  
  VlnPlot(healthy_heart, features = favorable_markers,combine=F,group.by = 'cell.name')
}
dev.off()

saveRDS(healthy_heart,file="~/Desktop/BISPECIFIC_CODE_EUGENE/BSABS_SCRNA/healthy_heart.rds")

###############################################################
                             #LUNG#
###############################################################
load('GSE130148_raw_counts.RData')

cell_clusters=read.table('GSE130148_barcodes_cell_types.txt',header=T,sep='\t')
row.names(cell_clusters)=cell_clusters$cell.barcode

healthy_lung <- CreateSeuratObject(counts = raw_counts, project = "hcc10k", min.cells = 3, min.features = 10,meta.data = cell_clusters)
healthy_lung[["percent.mt"]] <- PercentageFeatureSet(healthy_lung, pattern = "^MT-")
VlnPlot(healthy_lung, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 0)
healthy_lung <- subset(healthy_lung, subset = nFeature_RNA > 0 & percent.mt < 50 )
healthy_lung <- NormalizeData(healthy_lung, normalization.method = "LogNormalize", scale.factor = 10000)
healthy_lung <- FindVariableFeatures(healthy_lung, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(healthy_lung)
healthy_lung <- ScaleData(healthy_lung, features = all.genes)
healthy_lung <- FindVariableFeatures(healthy_lung)
healthy_lung <- RunPCA(healthy_lung, features = VariableFeatures(object = healthy_lung))

VizDimLoadings(healthy_lung, dims = 1:2, reduction = "pca")
healthy_lung <- FindNeighbors(healthy_lung, dims = 1:30,verbose=TRUE)

set.seed(123)
healthy_lung <- FindClusters(healthy_lung, resolution = 0.8)
head(Idents(healthy_lung), 5)
healthy_lung <- RunUMAP(healthy_lung, dims = 1:30)
ElbowPlot(healthy_lung,reduction='umap')

pdf('lung_vital.pdf')
{
DimPlot(healthy_lung, reduction = "umap",label=T)+ggtitle('UMAP of healthy lung tissues')
DimPlot(healthy_lung, reduction = "umap",label=T,group.by = 'celltype')+ggtitle('UMAP of healthy lung tissues')

favorable_markers=c('GPC3','PIGY','TMEM150B','MUC13','TREM2','SLC51B','MELK','DTL','TGM3','GPR158')
favorable_markers=favorable_markers[favorable_markers%in%row.names(raw_counts)]

FeaturePlot(healthy_lung, features = favorable_markers)

VlnPlot(healthy_lung, features = 'HOPX',combine=F,group.by = 'celltype')
}

dev.off()

saveRDS(healthy_lung,file="~/Desktop/BISPECIFIC_CODE_EUGENE/BSABS_SCRNA/healthy_lung.rds")

###############################################################
                               #KIDNEY#
###############################################################

kidney_reference=read.csv('kidney_ref.csv')
table(kidney_reference$cluster)
head(kidney_reference)

counts2=read.csv('~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/KIDNEY/GSM5224978_counts_A1_S1.csv')
counts2=read.csv('D:/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/KIDNEY/GSM5224978_counts_A1_S1.csv')

counts2=counts2[-c(1:2),]
row.names(counts2)=counts2$X
counts2$X=NULL

kidney_scrna <- CreateSeuratObject(counts = counts2, project = "hcc10k", min.cells = 0, min.features = 0)
kidney_scrna[["percent.mt"]] <- PercentageFeatureSet(kidney_scrna, pattern = "^MT-")
VlnPlot(kidney_scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 0)
kidney_scrna <- subset(kidney_scrna, subset = nFeature_RNA > 0 & percent.mt < 50 )
kidney_scrna <- NormalizeData(kidney_scrna, normalization.method = "LogNormalize", scale.factor = 10000)
kidney_scrna <- FindVariableFeatures(kidney_scrna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(kidney_scrna)
kidney_scrna <- ScaleData(kidney_scrna, features = all.genes)
kidney_scrna <- FindVariableFeatures(kidney_scrna)
kidney_scrna <- RunPCA(kidney_scrna, features = VariableFeatures(object = kidney_scrna))
VizDimLoadings(kidney_scrna, dims = 1:2, reduction = "pca")
kidney_scrna <- FindNeighbors(kidney_scrna, dims = 1:30,verbose=TRUE)

set.seed(123)

kidney_scrna <- FindClusters(kidney_scrna, resolution = 0.8)
head(Idents(kidney_scrna), 5)
kidney_scrna <- RunUMAP(kidney_scrna, dims = 1:30)
ElbowPlot(kidney_scrna,reduction='umap')

DimPlot(kidney_scrna, reduction = "umap",label=T)+ggtitle('UMAP of healthy kidney tissues')

##############################define clusters
CombinePlots(plots = list(
  FeaturePlot(kidney_scrna, features = 'LINC01983',label=T), 
  FeaturePlot(kidney_scrna, features = 'CRABP1',label=T), 
  FeaturePlot(kidney_scrna, features = 'APOC1',label=T) 
))
#Glomerulus=4
CombinePlots(plots = list(
  FeaturePlot(kidney_scrna, features = 'NPHS1',label=T),
  FeaturePlot(kidney_scrna, features = 'TYRO3',label=T),
  FeaturePlot(kidney_scrna, features = 'ATP10A',label=T)
))
#Interstitium =8
CombinePlots(plots = list(
  FeaturePlot(kidney_scrna, features = 'MYH11',label=T),
  FeaturePlot(kidney_scrna, features = 'RERGL',label=T),
  FeaturePlot(kidney_scrna, features = 'MCAM',label=T)
))
#PT pure =2
CombinePlots(plots = list(
  FeaturePlot(kidney_scrna, features = 'LINC01874',label=T),
  FeaturePlot(kidney_scrna, features = 'ALB',label=T),
  FeaturePlot(kidney_scrna, features = 'MCAM',label=T)
))
#TAL pure=7,0
CombinePlots(plots = list(
  FeaturePlot(kidney_scrna, features = 'ANKRD2',label=T),
  FeaturePlot(kidney_scrna, features = 'TDGF1',label=T),
  FeaturePlot(kidney_scrna, features = 'PLAU',label=T),
  FeaturePlot(kidney_scrna, features = 'GNG7',label=T),
  FeaturePlot(kidney_scrna, features = 'PRR15L',label=T),
  FeaturePlot(kidney_scrna, features = 'HSPB7',label=T)
))

new.cluster.ids <- c("TAL pure","TAL/PT","PT pure","DCT-CNT","Glomerulus","CD","DCT","TAL pure","Interstitium")
names(new.cluster.ids) <- levels(kidney_scrna)
kidney_scrna <- RenameIdents(kidney_scrna, new.cluster.ids)
DimPlot(kidney_scrna, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()+ggtitle('Labeled healthy cell populations')
kidney_scrna$cell_type=kidney_scrna@active.ident

saveRDS(kidney_scrna,file='kidney_scrna_2.rds')

favorable_markers=c('GPC3','PIGY','TMEM150B','MUC13','TREM2','SLC51B','MELK','DTL','TGM3','GPR158')
favorable_markers=favorable_markers[favorable_markers%in%row.names(counts2)]

FeaturePlot(kidney_scrna, features = favorable_markers,label=T)

VlnPlot(kidney_scrna, features = favorable_markers,combine=F)


##############################################################################################
                                   #BRAIN#
#############################################################################################
healthy_brain=readH5AD('habib17.processed.h5ad')
brain_meta=data.frame(cellid=colnames(healthy_brain),celltype=healthy_brain$CellType)
row.names(brain_meta)=brain_meta$cellid
healthy_brain <- CreateSeuratObject(counts=as.matrix(assay(healthy_brain)),meta=brain_meta)
healthy_brain[["percent.mt"]] <- PercentageFeatureSet(healthy_brain, pattern = "^MT-")
VlnPlot(healthy_brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 0)
healthy_brain <- subset(healthy_brain, subset = nFeature_RNA > 0 & percent.mt < 50 )
healthy_brain <- NormalizeData(healthy_brain, normalization.method = "LogNormalize", scale.factor = 10000)
healthy_brain <- FindVariableFeatures(healthy_brain, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(healthy_brain)
healthy_brain <- ScaleData(healthy_brain, features = all.genes)
healthy_brain <- FindVariableFeatures(healthy_brain)
healthy_brain <- RunPCA(healthy_brain, features = VariableFeatures(object = healthy_brain))
VizDimLoadings(healthy_brain, dims = 1:2, reduction = "pca")
healthy_brain <- FindNeighbors(healthy_brain, dims = 1:30,verbose=TRUE)

set.seed(123)
healthy_brain <- FindClusters(healthy_brain, resolution = 0.8)
head(Idents(healthy_brain), 5)
healthy_brain <- RunUMAP(healthy_brain, dims = 1:30)
ElbowPlot(healthy_brain,reduction='umap')

FeaturePlot(healthy_brain, features = favorable_markers)

VlnPlot(healthy_brain, features = favorable_markers,combine=F,group.by = 'cellType')

DimPlot(healthy_brain, reduction = "umap",label=T)+ggtitle('UMAP of healthy brain tissues')

saveRDS(healthy_brain,file="~/Desktop/BISPECIFIC_CODE_EUGENE/BSABS_SCRNA/healthy_brain.rds")


###############################################################
                              #LIVER#

###############################################################
 
healthy_liver=read.table('GSE124395_Normalhumanlivercellatlasdata.txt')
liver_clusters=read.table('GSE124395_clusterpartition.txt')
liver_clusters$cell.type='Unknown'
liver_clusters[liver_clusters$sct.cpart%in%c(5,1,19,3,28,12,18),'cell.type']='NK, NKT, T cells'
liver_clusters[liver_clusters$sct.cpart%in%c(20,9,13,32,23),'cell.type']='Liver sinusoidal end'
liver_clusters[liver_clusters$sct.cpart%in%c(29,33,35,26),'cell.type']='Macrovascular end'
liver_clusters[liver_clusters$sct.cpart%in%c(36,21,9),'cell.type']='Unknown'
liver_clusters[liver_clusters$sct.cpart%in%c(11,17,30,14,8),'cell.type']='Hepatocytes'
liver_clusters[liver_clusters$sct.cpart%in%c(22,38,16,34,37),'cell.type']='B cells'
liver_clusters[liver_clusters$sct.cpart%in%c(25,31,6,23),'cell.type']='Kupfer cells'
liver_clusters[liver_clusters$sct.cpart%in%c(4,7,24,39,11),'cell.type']='EPCAM+'

table(liver_clusters$cell.type)
unique(subset(liver_clusters,cell.type=='A')$sct.cpart)

healthy_liver <- CreateSeuratObject(counts = healthy_liver, project = "hcc10k", min.cells = 0, min.features = 0,meta.data = liver_clusters)
healthy_liver[["percent.mt"]] <- PercentageFeatureSet(healthy_liver, pattern = "^MT-")
VlnPlot(healthy_liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 0)
healthy_liver <- subset(healthy_liver, subset = nFeature_RNA > 0 & percent.mt < 50 )
healthy_liver <- NormalizeData(healthy_liver, normalization.method = "LogNormalize", scale.factor = 10000)
healthy_liver <- FindVariableFeatures(healthy_liver, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(healthy_liver)
healthy_liver <- ScaleData(healthy_liver, features = all.genes)
healthy_liver <- FindVariableFeatures(healthy_liver)
healthy_liver <- RunPCA(healthy_liver, features = VariableFeatures(object = healthy_liver))
VizDimLoadings(healthy_liver, dims = 1:2, reduction = "pca")
healthy_liver <- FindNeighbors(healthy_liver, dims = 1:30,verbose=TRUE)

set.seed(123)
healthy_liver <- FindClusters(healthy_liver, resolution = 0.8)
head(Idents(healthy_liver), 5)
healthy_liver <- RunUMAP(healthy_liver, dims = 1:30)
ElbowPlot(healthy_liver,reduction='umap')

saveRDS(healthy_liver,file='healthy_liver.rds')

pdf('liver_vital.pdf')
CombinePlots(plots = list(
  DimPlot(healthy_liver, reduction = "umap",label=T,group.by = 'cell.type')+ggtitle('UMAP of healthy liver tissues'),
  DimPlot(healthy_liver, reduction = "umap",label=T)+ggtitle('UMAP of healthy liver tissues')))


CombinePlots(plots = list(
  FeaturePlot(healthy_liver, features = 'LINC01983',label=T), 
  FeaturePlot(healthy_liver, features = 'CRABP1',label=T), 
  FeaturePlot(healthy_liver, features = 'APOC1',label=T) 
))


favorable_markers=c('GPC3','PIGY','TMEM150B','MUC13','TREM2','SLC51B','MELK','DTL','TGM3','GPR158')

FeaturePlot(healthy_liver, features = favorable_markers_short)
FeaturePlot(healthy_liver, features = favorable_markers)

VlnPlot(healthy_liver, features = favorable_markers,combine=F,group.by = 'cell.type')
dev.off()

###########################################################################################################################################










