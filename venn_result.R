setwd("D:/课题总结/TGM2转录组分析/count_wjx/王佳鑫-figure")
#程涛的ery的regulon列表

rasMat<-read.table("D:/analysis/forpublication/human_abc/GSE149938/ery/int/AUCell.txt")
chengtao_ery<-sub("_.*", "", colnames(rasMat))
chengtao_ery<-unique(chengtao_ery)

write.csv(chengtao_ery,"chengtao_ery.csv")



ery<-readRDS("D:/analysis/forpublication/human_abc/GSE149938/ery.rds")
#平均表达值
YDL.AVERAGE<-AverageExpression(object = ery,return.seurat=F)
write.csv(YDL.AVERAGE,file = "AVERAGE_ery_程涛.csv")

pdf("REGULON_EXPRESSION_IKZF1_程涛.pdf",height = 5,width = 5)
DimPlot(ery, reduction = "tsne", label = TRUE, pt.size = 1.5)
p1<-FeaturePlot(object = ery, reduction = "tsne",pt.size = 1.5,features = "IKZF1",cols = c("gray", "red"),order = T)#actin
print(p1)
dev.off()



# #徐湘民的regulon列表
# rasMat<-read.table("D:/analysis/forpublication/human_bm_ubc/int/AUCell.txt")
# PNAS_ery<-sub("_.*", "", colnames(rasMat))
# PNAS_ery<-unique(PNAS_ery)
# 


#NCB的ery的regulon列表:小鼠
rasMat<-read.table("D:/课题总结/TGM2转录组分析/count_wjx/NCB/AUCell.txt")
NCB_ery<-sub("_.*", "", colnames(rasMat))
NCB_ery<-unique(NCB_ery)
NCB_ery<-toupper(NCB_ery)
write.csv(NCB_ery,"NCB_ery.csv")

ery<-readRDS("D:/analysis/bone_marrow_niche/ery.rds")
current.cluster.ids <- c(0, 1, 2, 3)
new.cluster.ids <- c(
  "Erythroblasts",
  "Erythroblasts",
  "Erythroid megakaryocyte","Erythroid progenitor")
names(new.cluster.ids) <- levels(ery)
ery <- RenameIdents(ery, new.cluster.ids)
DimPlot(ery, reduction = "tsne", label = TRUE, pt.size = 1.5)
#平均表达值
YDL.AVERAGE<-AverageExpression(object = ery,return.seurat=F)
write.csv(YDL.AVERAGE,file = "AVERAGE_ery_NCB.csv")


pdf("REGULON_EXPRESSION_IKZF1_NCB.pdf",height = 5,width = 5)
DimPlot(ery, reduction = "tsne", label = TRUE, pt.size = 1.5)
p1<-FeaturePlot(object = ery, reduction = "tsne",pt.size = 1.5,features = "Ikzf1",cols = c("gray", "red"),order = T)#actin
print(p1)

dev.off()


library (ggvenn)
x<-list("chengtao_ery"=chengtao_ery,
        "NCB_ery"=NCB_ery)
ggvenn(x)

intersect<-intersect(chengtao_ery,NCB_ery)
write.csv(intersect,"intersect.csv")


x<-list("chengtao_ery"=chengtao_ery,
        "PNAS_ery"=PNAS_ery,
        "NCB_ery"=NCB_ery)
ggvenn(x)

intersect<-intersect(chengtao_ery,PNAS_ery)
intersect<-intersect(intersect,NCB_ery)

library(data.table)
merged_list <- rbindlist(chengtao_ery,PNAS_ery,fill = TRUE)
data<-cbind(chengtao_ery,PNAS_ery)
write.table(as.matrix(x),"三个样本的regulon.txt")









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
ery<-readRDS("ery.rds")
DimPlot(ery,reduction = "tsne",label = TRUE,pt.size = 1.5)


ery<-readRDS("D:/analysis/bone_marrow_niche/ery.rds")
DimPlot(ery,reduction = "tsne",label = TRUE,pt.size = 1.5)
