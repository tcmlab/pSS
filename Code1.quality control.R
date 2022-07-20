## Title     :1.数据整合质控
# Objective : pSS数据
# Created by: LJK
# Created on: 2022/5/13
##========================
##=================================================================##
##==================第一节：创建Seurat对象并质控===================##
##=================================================================##

setwd("~/data/project/pSS/jidiao/afternames/clusters/All")
if(T){
  
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggplot2)
  rm(list=ls())
  gc()
  ##===创建seurat对象列表===##
  ##设置文件目录与样本名称
  dir <- dir("./rawdata/")
  dir <- paste0("rawdata/", dir)
  
  #查看文件顺序
  dir                         
  #按文件顺序给样本命名，名称不要以数字开头，中间不能有空格 
  samples_name = c('Control_1', 'Control_2', 'Control_3', 'pSS_1', 'pSS_2','pSS_3')
  
  ##使用循环命令批量创建seurat对象
  scRNAlist <- list()
  for(i in 1:length(dir)){
    counts <- Read10X(data.dir = dir[i])
    #不设置min.cells过滤基因会导致CellCycleScoring报错：
    #Insufficient data values to produce 24 bins.  
    scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i],
                                         min.cells=3, min.features = 200)
    #给细胞barcode加个前缀，防止合并后barcode重名
    scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])   
    #计算线粒体基因比例，小鼠用第二行命令
    if(T){    
      scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-") 
      #scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-")
    }
    #计算核糖体基因比例，小鼠用第二行命令
    if(T){
      scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RB[SL]")
      #scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^Rb[sl]")
    }
  }
  
  #给列表命名并保存数据
  names(scRNAlist) <- samples_name
  #save(scRNAlist, file = "scRNAlist0.Rdata")
  saveRDS(scRNAlist, file = "./data/scRNAlist0.rds")
  
  
  ##===数据质控===##
  dir.create("./data/QC")
  #合并数据
  scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
  scRNA$proj <- rep("10x", ncol(scRNA))
  #scRNAlist <- SplitObject(scRNA, split.by = "orig.ident")
  head(scRNA@meta.data)
  
  ##绘制质控小提琴图
  #设置可能用到的主题
  theme.set1 = theme(axis.title.x=element_blank(), 
                     axis.text.x=element_blank(), 
                     axis.ticks.x=element_blank())
  theme.set2 = theme_bw()+
               theme(axis.title.x=element_blank())+
               theme(panel.grid =element_blank()) + ## 删去网格线
               theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
  
  theme.set3 = theme_bw()+
    theme(panel.grid =element_blank()) + ## 删去网格线
    theme(plot.title = element_blank())+
    theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
  
  source("~/data/project/pSS/jidiao/afternames/clusters/All/color.r")
  color.cluster <- fetch_color(8,"tsne", "set1")
  #组的颜色是
  color.group <- fetch_color(2,"tsne", "set2")
  #样品的颜色
  color.sample <- fetch_color(6,"tsne", "set3")
  
  
  #设置绘图元素
  plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb")

  #质控前小提琴图
  plots = list()
  for(i in seq_along(plot.featrures)){
    plots[[i]] = VlnPlot(scRNA, group.by='orig.ident', pt.size = 0, cols = color.sample,
                         features = plot.featrures[i]) +theme.set2+RotatedAxis() + NoLegend()
    }
  violin <- wrap_plots(plots = plots, nrow=1)    
  ggsave("Fig/QC/vlnplot_before_qc.pdf", plot = violin, width = 10, height = 3) 
  ggsave("Fig/QC/vlnplot_before_qc.png", plot = violin, width = 10, height = 3)  
  
  
  ##设置质控标准
  minGene=200
  maxGene=3000
  CountRNA=10000
  pctMT=10
  pctRB=0.5
  
  ##数据质控并绘制小提琴图
  scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < 
                    maxGene & nCount_RNA < CountRNA & percent.mt < pctMT & percent.rb < pctRB)
  plots = list()
  for(i in seq_along(plot.featrures)){
                                 plots[[i]] = VlnPlot(scRNA, group.by='orig.ident', pt.size = 0, cols = color.sample,
                                 features = plot.featrures[i]) +theme.set2+RotatedAxis() + NoLegend()
                                }
  violin <- wrap_plots(plots = plots, nrow=1)     
  ggsave("Fig/QC/vlnplot_after_qc.pdf", plot = violin, width = 10, height = 3) 
  ggsave("Fig/QC/vlnplot_after_qc.png", plot = violin, width = 10, height = 3)
  
  
  ##===细胞周期评分===##
  #默认情况下，Seurat使用global-scaling的归一化方法，称为“LogNormalize”，这种方法是利用总的表达量对每个细胞里的基因表达值进行归一化，乘以一个scale factor（默认值是10000），再用log转换一下。归一化后的数据存放在pbmc[["RNA"]]@data里。
  #Seurat使用FindVariableFeatures函数鉴定高可变基因，默认情况下，会返回2,000个高可变基因用于下游的分析，如PCA等。
  #ScaleData默认对之前鉴定到的2000个高可变基因进行标准化，也可以通过vars.to.regress参数指定其他的变量对数据进行标准化，表达矩阵进行scaling后，其结果存储在pbmc[["RNA"]]@scale.data中
  
  scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures() %>%
    ScaleData(features = rownames(scRNA))
  
  ##人源样本的细胞周期评分
  if(T){
    g2m_genes <- cc.genes$g2m.genes
    g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(scRNA))
    s_genes <- cc.genes$s.genes    
    s_genes <- CaseMatch(search=s_genes, match=rownames(scRNA))
    scRNA <- CellCycleScoring(scRNA, g2m.features=g2m_genes, s.features=s_genes)
    tmp <- RunPCA(scRNA, features = c(g2m_genes, s_genes), verbose = F)
    p <- DimPlot(tmp, reduction = "pca", group.by = "orig.ident", cols = color.sample)+ 
         theme.set3+theme(plot.title = element_blank())
    ggsave("Fig/QC/CellCycle_pca.png", p, width = 8, height = 6)
    ggsave("Fig/QC/CellCycle_pca.pdf", p, width = 8, height = 6)
    rm(tmp)
  }
  
 
  ##===线粒体和核糖体影响作图===##
  ##查看线粒体对样本的影响，可直接用于小鼠
  mt.genes <- grep("^MT-", rownames(scRNA), value=T, ignore.case=T)
  tmp <- RunPCA(scRNA, features = mt.genes, verbose = F)
  p <- DimPlot(tmp, reduction = "pca", group.by = "orig.ident",cols = color.sample)+
       theme.set3+theme(plot.title = element_blank())
  ggsave("Fig/QC/mito_pca.png", p, width = 8, height = 6)
  ggsave("Fig/QC/mito_pca.pdf", p, width = 8, height = 6)
  rm(tmp)
  
  ##查看核糖体对样本的影响，可直接用于小鼠
  rb.genes <- grep("^RP[SL]", rownames(scRNA), value=T, ignore.case=T)
  tmp <- RunPCA(scRNA, features = rb.genes, verbose = F)
  p <- DimPlot(tmp, reduction = "pca", group.by = "orig.ident" ,cols = color.sample)+
       theme.set3+theme(plot.title = element_blank())
  ggsave("Fig/QC/ribo_pca.png", p, width = 8, height = 6)
  ggsave("Fig/QC/ribo_pca.pdf", p, width = 8, height = 6)
  rm(tmp)
  
  ##保存质控后的结果
  #saveRDS(scRNA, "scRNA.rds")
  save(scRNA, file = "./data/scRNA_QC.Rdata")
}

