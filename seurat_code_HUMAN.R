Sys.setenv(LANGUAGE = "en") #显示英文报错信息
gc()
memory.limit(9999999999)
set.seed(123)
rm(list = ls())  
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)

setwd("D:/analysis/forpublication/human_abc/GSE149938")

data<-read.csv("umi_matrix.csv")
head(data[1:10,1:4])
YDL <- CreateSeuratObject(counts = t(data),min.cells=10,min.genes=200,project = "human_abc")#11559 features across 3637 samples
###展示基因及线粒体百分比（这里将其进行标记并统计其分布频率，"nFeature_RNA"为基因数，"nCount_RNA"为UMI数，"percent.mt"为线粒体占比）
VlnPlot(object = YDL, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#对数据进行标准化
##表达量数据标准化,LogNormalize的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000
YDL <- NormalizeData(object = YDL, normalization.method = "LogNormalize", scale.factor = 10000)
#提取那些在细胞间变异系数较大的基因
##鉴定表达高变基因(2000个）,用于下游分析,如PCA；
YDL <- FindVariableFeatures(object = YDL, selection.method = "vst", nfeatures = 2000)
YDL <- ScaleData(object = YDL, features = rownames(YDL))
#线性降维（PCA）,默认用高变基因集,但也可通过features参数自己指定；
YDL=RunPCA(object= YDL,npcs = 20,pc.genes=VariableFeatures(object = YDL))     #PCA分析
ElbowPlot(YDL)#选择top20个PC
pcSelect=15
YDL <- FindNeighbors(object = YDL, dims = 1:pcSelect)                #计算邻接距离
##接着优化模型,resolution参数决定下游聚类分析得到的分群数,对于3K左右的细胞,设为0.4-1.2 能得到较好的结果(官方说明)；如果数据量增大,该参数也应该适当增大。
YDL <- FindClusters(object = YDL, resolution = 0.3)                  #对细胞分组,优化标准模块化
##使用Idents（）函数可查看不同细胞的分群；
head(Idents(YDL), 5)

##这里采用基于TSNE的聚类方法。
YDL <- RunTSNE(object = YDL, dims = 1:pcSelect,check_duplicates = FALSE)                      #TSNE聚类
pdf(file="TSNE_BM.pdf",width=6.5,height=6)
TSNEPlot(object = YDL, pt.size = 2, label = TRUE)    #TSNE可视化

#另一个可视化的方法
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
dev.off()
write.table(YDL$seurat_clusters,file="tsneCluster_BM.txt",quote=F,sep="\t",col.names=F)


##这里采用基于umap的聚类方法
#这里采用基于图论的聚类方法
YDL <- RunUMAP(object = YDL, dims = 1:pcSelect)                      #umap聚类
pdf(file="umap__BM.pdf",width=6.5,height=6)
UMAPPlot(object = YDL, pt.size = 1.5, label = TRUE)    #umap可视化
#另一个可视化的方法
DimPlot(object=YDL,label = TRUE,reduction="umap")
##用DimPlot()函数绘制散点图,reduction = "tsne",指定绘制类型；如果不指定,默认先从搜索 umap,然后 tsne, 再然后 pca；也可以直接使用这3个函数PCAPlot()、TSNEPlot()、UMAPPlot()； cols,pt.size分别调整分组颜色和点的大小；
dev.off()

##细胞周期归类

YDL<- CellCycleScoring(object = YDL, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

head(x = YDL@meta.data)
pdf(file="CellCycle_BM.pdf",width=6.5,height=6)
DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,group.by="Phase",pt.size = 1.5)
dev.off()
#定义细胞类型，根据barcode定义。

YDL[["celltype"]] <- YDL@meta.data$orig.ident

head(colnames(YDL),30)
DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="celltype",pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="celltype",pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,group.by="celltype",pt.size = 1.5)
DimPlot(YDL,reduction = "umap",label = TRUE,group.by="celltype",pt.size = 1.5,split.by = "orig.ident")


#选择特定类型的细胞
Idents(YDL) <- YDL@meta.data$orig.ident
#YDL <- subset(YDL, idents = 'Redblood')
#save(YDL, file = 'YDL.Redblood.Rdata')
Idents(YDL) <- YDL@meta.data$seurat_clusters
#另一个可视化的方法
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)


YDL.markers <- FindAllMarkers(YDL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'roc')
# install.packages("magrittr") # package installations are only needed the first time you use it
# install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run evYDL time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
YDL.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
##存储marker

write.csv(YDL.markers,file="allmarker_celltpye.csv")

#绘制分cluster的热图
top10 <- YDL.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#绘制marker在各个cluster的热图
pdf(file="tsneHeatmap.pdf",width=12,height=36)
DoHeatmap(object = YDL, features = top10$gene) + NoLegend()
DoHeatmap(subset(YDL, downsample = 100), features = top10$gene, size = 3)+ NoLegend()
dev.off()





ery<-subset(YDL, cells= rownames(YDL@meta.data[YDL@meta.data$orig.ident=="ery",]))
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(ery,reduction = "tsne",label = TRUE,pt.size = 1.5)

##鉴定表达高变基因(2000个）,用于下游分析,如PCA；
ery <- FindVariableFeatures(object = ery, selection.method = "vst", nfeatures = 2000)
ery <- ScaleData(object = ery, features = rownames(ery))
#线性降维（PCA）,默认用高变基因集,但也可通过features参数自己指定；
ery=RunPCA(object= ery,npcs = 20,pc.genes=VariableFeatures(object = ery))     #PCA分析
ElbowPlot(ery)#选择top20个PC
pcSelect=15
ery <- FindNeighbors(object = ery, dims = 1:pcSelect)     
##接着优化模型,resolution参数决定下游聚类分析得到的分群数,对于3K左右的细胞,设为0.4-1.2 能得到较好的结果(官方说明)；如果数据量增大,该参数也应该适当增大。
ery<- FindClusters(object = ery, resolution = 0.3)                  #对细胞分组,优化标准模块化
##使用Idents（）函数可查看不同细胞的分群；
head(Idents(ery), 5)

##这里采用基于TSNE的聚类方法。
ery<- RunTSNE(object = ery, dims = 1:pcSelect,check_duplicates = FALSE)                      #TSNE聚类
pdf(file="TSNE_rty.pdf",width=6.5,height=6)
TSNEPlot(object = ery, pt.size = 2, label = TRUE)    #TSNE可视化

#另一个可视化的方法
DimPlot(ery,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(ery,reduction = "tsne",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
dev.off()
write.table(ery$seurat_clusters,file="tsneCluster_ery.txt",quote=F,sep="\t",col.names=F)


##这里采用基于umap的聚类方法
#这里采用基于图论的聚类方法
ery <- RunUMAP(object = ery, dims = 1:pcSelect)                      #umap聚类
pdf(file="umap__ery.pdf",width=6.5,height=6)
UMAPPlot(object = ery, pt.size = 1.5, label = TRUE)    #umap可视化
#另一个可视化的方法
DimPlot(object=ery,label = TRUE,reduction="umap")
##用DimPlot()函数绘制散点图,reduction = "tsne",指定绘制类型；如果不指定,默认先从搜索 umap,然后 tsne, 再然后 pca；也可以直接使用这3个函数PCAPlot()、TSNEPlot()、UMAPPlot()； cols,pt.size分别调整分组颜色和点的大小；
dev.off()

##细胞周期归类

ery<- CellCycleScoring(object = ery, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

head(x = ery@meta.data)
pdf(file="CellCycle_ery.pdf",width=6.5,height=6)
DimPlot(ery,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5)
DimPlot(ery,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5,split.by = "orig.ident")
DimPlot(ery,reduction = "umap",label = TRUE,group.by="Phase",pt.size = 1.5)
dev.off()


ery.markers <- FindAllMarkers(ery, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'roc')
# install.packages("magrittr") # package installations are only needed the first time you use it
# install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run evYDL time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
ery.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
##存储marker

write.csv(ery.markers,file="allmarker_celltpye_ery.csv")

#绘制分cluster的热图
top10 <- ery.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#绘制marker在各个cluster的热图
pdf(file="tsneHeatmap_ery.pdf",width=12,height=36)
DoHeatmap(object = ery, features = top10$gene) + NoLegend()
DoHeatmap(subset(ery, downsample = 100), features = top10$gene, size = 3)+ NoLegend()
dev.off()


FeaturePlot(object = ery, reduction = "tsne",pt.size = 1.5,features = c("PTPRC"))#actin
FeaturePlot(object = ery, reduction = "tsne",pt.size = 1.5,features = c("S100A8"))#actin
FeaturePlot(object = ery, reduction = "tsne",pt.size = 1.5,features = c("CD79A"))#actin
FeaturePlot(object = ery, reduction = "tsne",pt.size = 1.5,features = c("HBB","HBB-Y"))#actin
FeaturePlot(object = ery, reduction = "tsne",pt.size = 1.5,features = c("GYPA"))#actin


FeaturePlot(object = ery, reduction = "tsne",pt.size = 1.5,features = c("TFRC"))#actin
FeaturePlot(object = ery, reduction = "tsne",pt.size = 1.5,features = c("HLA-DRA"))#actin
FeaturePlot(object = ery, reduction = "tsne",pt.size = 1.5,features = c("CD74"))#actin
FeaturePlot(object = ery, reduction = "tsne",pt.size = 1.5,features = c("VCAN"))#actin
FeaturePlot(object = ery, reduction = "tsne",pt.size = 1.5,features = c("GATA1"))#actin
FeaturePlot(object = ery, reduction = "tsne",pt.size = 1.5,features = c("S100A9"))#actin

FeaturePlot(object = ery, reduction = "tsne",pt.size = 1.5,features = c("XPO7"))#actin
saveRDS(ery,"ery.rds")



library(monocle)
#准备monocle分析需要的文件
monocle.matrix=as.matrix(ery@assays$RNA@data)
monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
write.table(monocle.matrix,file="monocleMatrix.txt",quote=F,sep="\t",row.names=F)
monocle.sample=as.matrix(ery@meta.data)
monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
write.table(monocle.sample,file="monocleSample.txt",quote=F,sep="\t",row.names=F)
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
write.table(monocle.geneAnn,file="monocleGene.txt",quote=F,sep="\t",row.names=F)
write.table(ery.markers,file="monocleMarkers.txt",sep="\t",row.names=F,quote=F)


#设置工作目录
monocle.matrix=read.table("monocleMatrix.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.sample=read.table("monocleSample.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.geneAnn=read.table("monocleGene.txt",sep="\t",header=T,row.names=1,check.names=F)
marker=read.table("monocleMarkers.txt",sep="\t",header=T,check.names=F)

#将Seurat结果转换为monocle需要的细胞矩阵，细胞注释表和基因注释表表
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)


#给其中一列数据重命名
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])
# saveRDS(cds,"cds.rds")
# rm(list = ls())  
# marker<-read.csv("allmarker.csv")
# cds<-readRDS("cds.rds")
#伪时间分析流程
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, marker$gene)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,reduction_method = 'DDRTree')
cds <- orderCells(cds)
pdf(file="cluster.trajectory_SINGLET_ery.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "Cluster")
plot_cell_trajectory(cds,color_by="Cluster")+facet_wrap(~Cluster,nrow=3,ncol = 3)
eplot_cell_trajectory(cds,color_by="orig.ident")+facet_wrap(~orig.ident,nrow=2,ncol = 5)
plot_cell_trajectory(cds,color_by="Cluster")+facet_wrap(~orig.ident,nrow=2,ncol = 3)

plot_cell_trajectory(cds, color_by = "ery@active.ident")
plot_cell_trajectory(cds, color_by="Pseudotime", show_backbone=FALSE)
# 可以很明显看到细胞的发育轨迹 
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "State") +facet_wrap(~State, nrow = 1)
dev.off()

saveRDS(cds,"cds_SINGLET_ery.rds")
