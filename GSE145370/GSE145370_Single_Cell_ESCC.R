###Single cell esophagus cancer
library(Seurat)
library(dplyr)
#sce <- Read10X('./GSM4317409/filtered_feature_bc_matrix/')
#sce
createOb <- function(sample){
  #rawcount <- read.csv(paste0(sample,'_hs_RSEC_MolsPerCell.csv'),comment.char="#",row.names=1)
  sampleDir <- paste0('./',sample,'/filtered_feature_bc_matrix/')
  print(sampleDir)
  dat1 <- Read10X(sampleDir)
  samplename <- sample
  
  
  nCoV.seurat <- CreateSeuratObject(counts = dat1, project = samplename, min.cells = 3, min.features = 200)
  nCoV.seurat[['percent.mito']] <- PercentageFeatureSet(nCoV.seurat, pattern = "^MT-")
  pdf(paste0(samplename,'feature_ncount_percentmito.pdf'),width=15)
  print(VlnPlot(object = nCoV.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3))
  dev.off()
  nCoV.seurat <- subset(nCoV.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mito < 20)
  return(nCoV.seurat)
}
samples=dir('./',pattern = "GSM")

sceList=lapply(samples,createOb)
options(future.globals.maxSize = 4000 * 1024^2)

for (i in 1:length(sceList)) {
  sceList[[i]] <- NormalizeData(sceList[[i]], verbose = FALSE)
  sceList[[i]] <- FindVariableFeatures(sceList[[i]], selection.method = "vst",nfeatures = 2000, verbose = FALSE)
}
names(sceList) <- samples
sceList.anchors <- FindIntegrationAnchors(object.list = sceList, dims = 1:30)
sceList.integrated <- IntegrateData(anchorset = sceList.anchors, dims = 1:30)
library(data.table)
saveRDS(sceList.integrated,file="sceList.integrated_step1_v1.RDS")



sceList.integrated <- readRDS('./sceList.integrated_step1_v1.RDS')
sceList.integrated <- ScaleData(sceList.integrated, verbose = FALSE)
sceList.integrated <- RunPCA(sceList.integrated, npcs = 50, verbose = FALSE)

#COVIDRunner.integrated <- readRDS('./sceList.integrated_step2_40PCs.RDS')
COVIDRunner.integrated <- sceList.integrated
pdf('PCA_VizDimLoading.pdf')
VizDimLoadings(COVIDRunner.integrated, dims = 1:2, reduction = "pca")
dev.off()
pdf('PCA_Dim_plot.pdf')
DimPlot(COVIDRunner.integrated, reduction = "pca")
dev.off()
pdf('PCA_DimHeatmap.pdf')
DimHeatmap(COVIDRunner.integrated, dims = 1:5, cells = 500, balanced = TRUE)
dev.off()

COVIDRunner.integrated <- JackStraw(COVIDRunner.integrated, num.replicate = 100)
COVIDRunner.integrated <- ScoreJackStraw(COVIDRunner.integrated, dims = 1:20)
pdf('JackStrawPlot.pdf')
JackStrawPlot(COVIDRunner.integrated, dims = 1:15)
dev.off()
pdf('ElbowPlot.pdf')
ElbowPlot(COVIDRunner.integrated)
dev.off()



COVIDRunner.integrated <- FindNeighbors(COVIDRunner.integrated, dims = 1:50)
COVIDRunner.integrated <- FindClusters(COVIDRunner.integrated, resolution = 1)
COVIDRunner.integrated <- RunUMAP(COVIDRunner.integrated, dims = 1:50)

pdf('Umap_plot_50PCs.pdf')
DimPlot(COVIDRunner.integrated, reduction = "umap",label=T)
dev.off()

saveRDS(COVIDRunner.integrated,file="sceList.integrated_step2_50PCs.RDS")

COVIDRunner.integrated <- readRDS("sceList.integrated_step2_50PCs.RDS")
sce.markers <- FindAllMarkers(COVIDRunner.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(sce.markers,file="50PCs_sce.markers.RDS")

write.table(as.data.frame(sce.markers),file="40PCs_clusters_markers.xls",sep="\t",quote = F)
top10 <- sce.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(COVIDRunner.integrated, features = top10$gene) + NoLegend()
MonacoImmuneData <- readRDS("/home/myang/soft/software/singleCell/SingleR_Ref/MonacoImmuneData.rds")
sce_for_SingleR <- GetAssayData(COVIDRunner.integrated, slot="data")
clusters <- COVIDRunner.integrated@meta.data$seurat_clusters
pred.hesc <- SingleR(test = sce_for_SingleR, ref =MonacoImmuneData,labels =MonacoImmuneData$label.main,method = "cluster", 
                     clusters = clusters)



library(data.table)
pheno <- fread('GSE145370_phen_data.csv')
pheno[1:5,1:5]
loc <- match(COVIDRunner.integrated@meta.data$orig.ident,pheno$geo_accession)
head(loc)
COVIDRunner.integrated@meta.data[,'Group'] = pheno$`tissue type:ch1`[loc]

cellannot <- fread('cell_annotation.xls')
head(cellannot)
loc <- match(COVIDRunner.integrated@meta.data$seurat_clusters,cellannot$V1)
COVIDRunner.integrated@meta.data[,'celltype'] = cellannot$V2[loc]
Idents(COVIDRunner.integrated) = "celltype"


###Figure2A
library(ggplot2)
library(ggsci)
p = DimPlot(COVIDRunner.integrated,split.by="Group",label=T) + scale_color_d3()
ggsave(p,filename="UMAP_plot_Group.pdf",width=10)

###Figure2B
p = VlnPlot(COVIDRunner.integrated,features=c('TREM2','SPP1','APOE','C1QC','C1QB','C1QA'),stack=T,flip=T)+scale_fill_d3()
ggsave(p,filename="Figure2B_VlnPlot.pdf")

###Figure2D
p = FeaturePlot(COVIDRunner.integrated,features=c('TREM2'),cols=c('gray','red'),label=T,raster=FALSE)
ggsave(p,filename="TREM2_featurePlot.pdf")

###Figure2E
MDSC.sce = subset(x = COVIDRunner.integrated, idents = "MDSC")

p = VlnPlot(MDSC.sce,features=c('TREM2','SPP1','APOE','C1QC','C1QB','C1QA'),stack=T,flip=T,split.by="Group") +scale_fill_d3()
ggsave(p,filename="MDSC_Vlnplot_tumor_vs_normal.pdf")


###Figure2F
plotdat <- COVIDRunner.integrated@meta.data
p=ggplot(data=plotdat, mapping=aes(x=Group,fill=celltype))+geom_bar(stat="count",width=0.5,position='fill')+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA)) +scale_fill_d3(palette='category20')

###Figure2G
sampleCellcount =as.data.frame(table(plotdat$orig.ident,plotdat$celltype))

samples = levels(factor(plotdat$orig.ident))
sampleCellpercent = sampleCellcount[sampleCellcount$Var1==samples[1],]
sampleCellpercent[,"Percentage"] = sampleCellpercent$Freq/sum(sampleCellpercent$Freq)

for(i in samples[2:length(samples)]){
cellper = sampleCellcount[sampleCellcount$Var1==i,]
cellper[,"Percentage"] = cellper$Freq/sum(cellper$Freq)
sampleCellpercent = rbind(sampleCellpercent,cellper)
}
colnames(sampleCellpercent) <- c("Sample","Celltype","Count","Percentage")
library(ggpubr)
library(ggsci)
groupdf <- plotdat[!duplicated(plotdat$orig.ident),]
groupdf <- groupdf[,c(1,7)]
loc <- match(sampleCellpercent$Sample,groupdf$orig.ident)
sampleCellpercent[,'group'] = groupdf$Group[loc]
MDSC = sampleCellpercent[sampleCellpercent$Celltype=="MDSC",]
p <- ggboxplot(MDSC, x = "group", y = "Percentage",color = "group",add = "jitter", shape = "group",bxp.errorbar = TRUE)+scale_color_aaas()+theme(panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA))
my_comparisons <- list( levels(factor(sampleCellpercent$group)))

p1=p + stat_compare_means(comparisons = my_comparisons)

ggsave(p1,filename="MDSC_cell_percent.pdf",width=5,height=5)


### supplementary talbe S3
sce.markers <- FindAllMarkers(COVIDRunner.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
 write.table(as.data.frame(sce.markers),file="cell_type_clusters_markers.xls",sep="\t",quote = F)
 
### supplementary talbe S4
 markers <- FindMarkers(COVIDRunner.integrated, ident.1 = "tumors", group.by = 'group', subset.ident = "MDSC")
 write.table(as.data.frame(MDSC.markers),file="MDSC_tumors_vs_Normal_markers.xls",sep="\t",quote = F)







