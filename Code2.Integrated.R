# Title     :code2 数据整合分群
# Objective : pSS数据处理
# Created by: LJK
# Created on: 2022/5/16
###################################################################

rm(list = ls())
gc()

library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
setwd("~/data/project/pSS/jidiao/afternames/clusters/All")
source("~/data/project/pSS/jidiao/afternames/clusters/All/color.r")

#组的颜色是
color.group <- fetch_color(2,"tsne", "set2")
#样品的颜色
color.sample <- fetch_color(6,"tsne", "set3")

scRNA <- readRDS( "./data/scRNA_QC.Rdata")

scRNAlist<- SplitObject(scRNA, split.by = "orig.ident")

for (i in 1:length(scRNAlist)) {
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst")
}

library(future) 
options(future.globals.maxSize = 500 * 1024^3) #500g内存，根据服务器来选
plan(multisession, workers = 10) #开启多核运算 (10个核)
#plan('sequential') #终止多核运算 
##以VariableFeatures为基础寻找锚点，运行时间较长
scRNA<- FindIntegrationAnchors(object.list = scRNAlist, dims = 1:50)
##利用锚点整合数据，运行时间较长
scRNA<- IntegrateData(anchorset =scRNA, dims = 1:50)
DefaultAssay(scRNA) <- "integrated"
#ElbowPlot(scRNA, ndims = 50)
pc.num=1:30
# 数据标准化
  scRNA <- ScaleData(scRNA, features = rownames(scRNA), verbose = FALSE)
#数据降维分群
  scRNA <- RunPCA(scRNA, verbose = FALSE) %>% 
  RunTSNE(reduction="pca", dims=pc.num) %>% 
  RunUMAP(reduction="pca", dims=pc.num) %>%
  FindNeighbors(reduction = "pca", dims=pc.num) %>% 
  FindClusters(resolution=0.5)

scRNA$group <- NA
scRNA$group[which(str_detect(scRNA$orig.ident, "^Control_"))] <- "HC"
scRNA$group[which(str_detect(scRNA$orig.ident, "^pSS_"))] <- "pSS"
table(scRNA$group)

mytheme<-theme_bw()+  
  theme(plot.title = element_blank())+
  theme(panel.grid =element_blank()) + ## 删去网格线
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))

dir.create("./data/Integrate")
dir.create("./Fig/Integrate")
##tSNE图展示结果
color.cluster <- fetch_color(27,"tsne", "set1")
p1 <- DimPlot(scRNA, reduction = "tsne", cols= color.sample, label = F,
              group.by = "orig.ident")  + mytheme
p2 <- DimPlot(scRNA, reduction = "tsne", cols= color.cluster, label = F) + mytheme+
  theme(legend.key.size = unit(10, "pt"))+ 
  guides(
    color = guide_legend(
      ncol = 3,
      override.aes = list(size = 4)))#+ guides(color=guide_legend(title = "Clusters"))
pc1 = p1|p2
pc1 
ggsave("./Fig/Integrate/tSNE_overview.png", pc1, width = 8, height = 3)
ggsave("./Fig/Integrate/tSNE_overview.pdf", pc1, width = 8, height = 3)
##UMAP图展示结果
p3 <- DimPlot(scRNA, reduction = "umap", cols= color.sample, label = F,
              group.by = "orig.ident")  + mytheme
p4 <- DimPlot(scRNA, reduction = "umap", cols= color.cluster, label = F) + mytheme+
  theme(legend.key.size = unit(8, "pt"))+ 
  guides(
    color = guide_legend(
      ncol = 3,
      override.aes = list(size = 4)))

pc2 = p3|p4
pc2
ggsave("./Fig/Integrate/UMAP_overview.png", pc2, width = 8, height = 3)
ggsave("./Fig/Integrate/UMAP_overview.pdf", pc2, width = 8, height = 3)

##查看批次效应
p5 <- DimPlot(scRNA, reduction = "tsne", group.by = "orig.ident", cols = color.sample) + mytheme
p6 <- DimPlot(scRNA, reduction = "umap", group.by = "orig.ident", cols = color.sample) + mytheme
  
pc3 = p5 + p6 + plot_layout(guides = "collect")
pc3 
ggsave("./Fig/Integrate/batch_overview.png", pc3, width = 8, height = 4)
ggsave("./Fig/Integrate/batch_overview.pdf", pc3, width = 8, height = 4)

##查看细胞周期
#细胞周期评分绘制小提琴图
theme.set2 = theme_bw()+
  theme(axis.title.x=element_blank())+
  theme(panel.grid =element_blank()) + ## 删去网格线
  theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))

plot.featrures = c("S.Score", "G2M.Score")
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by="seurat_clusters",
                       pt.size = 0,
                       cols = color.cluster,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=1)    
violin
ggsave("./Fig/Integrate/CellCycle_violin.png", plot = violin, width = 8, height = 6)  
ggsave("./Fig/Integrate/CellCycle_violin.pdf", plot = violin, width = 8, height = 6)  

#绘制tsne和umap图
p7 <- DimPlot(scRNA, reduction = "tsne", group.by = "Phase")+ mytheme
p8<-  DimPlot(scRNA, reduction = "umap", group.by = "Phase")+ mytheme
pc4= p7 + p8 + plot_layout(guides = "collect")
pc4
ggsave("./Fig/Integrate/CellCycle_dimplot.png", pc4, width = 8, height = 4)
ggsave("./Fig/Integrate/CellCycle_dimplot.pdf", pc4, width = 8, height = 4)

saveRDS(scRNA, file = "./data/scRNA_Integrate.rds")
