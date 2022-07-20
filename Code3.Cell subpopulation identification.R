# Title     :code3 细胞亚群鉴定
# Objective : pSS数据处理
# Created by: LJK
# Created on: 2022/5/16

#############
#细胞亚群鉴定

library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
setwd("~/data/project/pSS/jidiao/afternames/clusters/All")
source("~/data/project/pSS/jidiao/afternames/clusters/All/color.r")

#细胞亚群的颜色是
color.cluster <- fetch_color(5,"tsne", "set1")
#组的颜色是
color.group <- fetch_color(2,"tsne", "set2")
#样品的颜色
color.sample <- fetch_color(6,"tsne", "set3")


marker1<-read.table('./all-marker.txt',header = T)
##保存结果

scRNA <- readRDS(file = "./data/scRNA_Integrate.rds")
Idents(scRNA)<'seurat_clusters'

checkgene<- c('CD79A','CD79B',"MS4A1", 
              "S100A8","S100A9","CST3",'CD68',
              'CD14','FCGR3A','MS4A7',
              'CD1C','FCER1A','CLEC10A',
              'GNLY',"NKG7", 'CST7', 'FCER1G',
              'CD3D','CD3E','IL7R','CD4', 
              'CD8A','CD8B', 'FOXP3', 'CCR7',
              'PPBP', 'TUBB1','PF4')

p13<-DoHeatmap(object = subset(scRNA, downsample= 300),
               features =checkgene, 
               slot = 'scale.data',
               size=4, label = T )
p13
table(scRNA$seurat_clusters)
{
  celltype=data.frame(ClusterID=c(0:26),
                      celltype= 'cells') 
  celltype[celltype$ClusterID %in% c( 0,2,3,7,8,12,13,16,25),2]='T_cells'  
  celltype[celltype$ClusterID %in% c( 1,6, 9,18, 19,21,23,26 ),2]='Myeloid'
  celltype[celltype$ClusterID %in% c( 4,5,17,20),2]='NK_cells'
  celltype[celltype$ClusterID %in% c(19, 22 ),2]='Doublets'
  celltype[celltype$ClusterID %in% c(10,11,24 ),2]='B_cells'
  celltype[celltype$ClusterID %in% c(14, 15),2]='Platelets'
  celltype
  
  scRNA@meta.data$celltype = "NA"
  for(i in 1:nrow(celltype)){
    scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
  
  table(scRNA@meta.data$celltype)
}

Idents(scRNA ) <- "celltype" 
scRNA <- subset(scRNA, idents = "Doublets", invert = TRUE)
scRNA$celltype <- as.factor(as.character(scRNA$celltype)) 
saveRDS(scRNA, "./data/scRNA_celltype.rds") 
#
DefaultAssay(scRNA) <- 'RNA'
scRNA <- readRDS( "./data/scRNA_celltype.rds") 
ScaleData(scRNA, features = rownames(scRNA), verbose = FALSE)
## 在寻找差异基因之前，把默认的assay切换为RNA

ClusterMarker <- FindAllMarkers(scRNA, assay = "RNA", slot = "data", only.pos = T)

save(ClusterMarker, file = "./data/ClusterMarker.rda")

ClusterMarker <- ClusterMarker[,c(7,1:6)]
write.csv(ClusterMarker,'./data/ClusterMarkers.csv', row.names=F)

### 提取差异显著的marker genes
top = 30   #可根据需要调整
TopMarkers1 <- ClusterMarker %>% filter(p_val_adj == 0) %>% group_by(cluster) %>% 
  top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>%
  top_n(n = top, wt = avg_log2FC)
TopMarkers <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(TopMarkers,'./data/TopMarkers.csv', row.names=F)

### 提取没有核糖体的Markers
ClusterMarker_noRibo <- ClusterMarker[!grepl("^RP[SL]", ClusterMarker$gene, ignore.case = F),]
ClusterMarker_noRibo_noMito <- ClusterMarker_noRibo[!grepl("^MT-", ClusterMarker_noRibo$gene, ignore.case = F),]
top = 30   #可根据需要调整
TopMarkers1 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
ClusterMarker_noRibo_noMito <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(ClusterMarker_noRibo_noMito,'./data/TopMarkers_noRibo_noMito.csv', row.names=F)

mytheme<-theme_bw()+  
  theme(plot.title = element_blank())+
  theme(panel.grid =element_blank()) + ## 删去网格线
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))


p1=DimPlot(scRNA, 
           reduction = 'umap', 
           group.by = 'orig.ident',
           cols =  color.sample,
           label = F, pt.size = 0.5) + 
  mytheme +
  theme(legend.position = c(0.85, 0.2))+
  theme(legend.key.size = unit(6, "pt"),
        legend.key = element_rect(
          #color = "red", # 框线色
          fill = "transparent",colour = NA))+ 
  guides(
    color = guide_legend(
      ncol= 1,
      override.aes = list(size =4)
    )
  )

p2=DimPlot(scRNA, 
           reduction = 'umap', 
           group.by = 'celltype',
           cols =  color.cluster,
           label = F, pt.size = 0.5) + mytheme +
  theme(legend.position = c(0.85, 0.2),
        legend.key = element_blank())+
  theme(legend.key.size = unit(8, "pt"),
        legend.key = element_rect(
          #color = "red", # 框线色
          fill = "transparent",colour = NA))+ 
  guides(
    color = guide_legend(
      ncol= 1,
      override.aes = list(size = 4)
    )
  )
x=p1|p2
x
ggsave("./Fig/Celltype/Dimplot.pdf", x, width = 8, height = 4)
ggsave("./Fig/Celltype/Dimplot.png", x, width = 8, height = 4)

checkgene<- c('CD79A','CD79B',"MS4A1",
              "S100A8","S100A9","CST3",
              'GNLY',"NKG7", 'CST7',
              'PPBP', 'TUBB1','PF4',
              'CD3D','CD3G','IL7R'
)
checkgene <- CaseMatch(checkgene, rownames(scRNA))
checkgene <- as.character(checkgene)

features<- factor(checkgene)
My_levels <- c( "B cells","Myeloid cells", "NK cells", 'Platelets', "T cells")
Idents(scRNA) <- factor(Idents(scRNA), levels= My_levels)
Idents(scRNA) <-'celltype'

p3<-DoHeatmap(object = subset(scRNA, downsample=1000),
              group.colors=color.cluster,
              features =checkgene,
              slot = 'scale.data',
              size=4, 
              label =F)
p3

ggsave("./Fig/Celltype/Markers_heatmap_1.pdf",plot = p3, width =6, height = 4)
ggsave("./Fig/Celltype/Markers_heatmap_1.png",plot = p3, width =6, height = 4)

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
p5<-DoHeatmap(object = subset(scRNA, downsample=1000),
              features =checkgene,
              group.colors=color.cluster,
              label = F, #label = F不在热图的上方标注细胞类型，
              angle = 0,
              slot = "scale.data",#slot = "scale.data"使用scale之后的矩阵画图，默认就是这个
              disp.min=-2.5, 
              size = 4, 
              disp.max=2.5)+
  theme(axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8))+ 
  scale_fill_gradientn(colours = rev(mapal))

p5

ggsave("./Fig/Celltype/Markers_2.png",plot = p5, width =8, height = 6)
ggsave("./Fig/Celltype/Markers_2.pdf",plot = p5, width =8, height = 6)


top = 3   #可根据需要调整
TopMarkers1 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
p6<-DoHeatmap(object = subset(scRNA, downsample=1000),
              features =TopMarkers1$gene,
              group.colors=color.cluster,
              label = F, #label = F不在热图的上方标注细胞类型，
              angle = 0,
              slot = "scale.data",#slot = "scale.data"使用scale之后的矩阵画图，默认就是这个
              disp.min=-2.5, 
              size = 4, 
              disp.max=2.5)+
  theme(axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8))+ 
  scale_fill_gradientn(colours = rev(mapal))

p6
#########################
p7 <- VlnPlot(scRNA,  
              features =checkgene %>% rev(),
              #group.colors=color.cluster,
              pt.size = 0, 
              group.by = 'celltype',stack = T)+
  NoLegend() 
              
p7
ggsave("./Fig/Celltype/Markers_VlnPlot.pdf",plot = p7, width =8, height = 6)

library(reshape2)
vln.df=as.data.frame(scRNA[["RNA"]]@data[checkgene,])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")

scRNA$CB<- rownames(as.data.frame(scRNA@active.ident))
anno=scRNA@meta.data[,c("CB","celltype")]
vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = checkgene) #为了控制画图的基因顺序

p8 <-vln.df%>%ggplot(aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+
  scale_fill_brewer(palette = "Set3", direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
ggsave("./Fig/Celltype/Markers_VlnPlot_2.pdf",plot = p8, width =8, height = 8)
###############
## 好看的气泡图基于好看的配色方案
library(RColorBrewer) 
library(viridis)
library(wesanderson)
#install.packages('wesanderson')
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pie(rep(6,n), col=sample(color, n))
col_vector
col_vector =c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest1"), wes_palette("Cavalcanti1"), wes_palette("GrandBudapest2"), wes_palette("FantasticFox1"))
pal <- wes_palette("Zissou1", 10, type = "continuous")
pal2 <- wes_palette("Zissou1", 5, type = "continuous")
pal[3:10]


# dotplot
Idents(scRNA) <- "celltype"
scRNA$celltype<- factor(x = scRNA$celltype, levels = c("Myeloid","T_cells", "B_cells", "NK_cells", 'Platelets'))
table(scRNA$celltype)
library(ggplot2)
final.markers<- c("S100A8","S100A9","CST3",
                  'CD3D','CD3G','IL7R',
                  'CD79A','CD79B',"MS4A1",
                  'GNLY',"NKG7", 'CST7',
                  'PPBP', 'TUBB1','PF4'
                             
)
##  RotatedAxis() scale_colour_gradientn都是比较重要的
 dottheme<-theme_bw()+  
  theme(plot.title = element_blank())+
  theme(panel.grid =element_blank()) + ## 删去网格线
  theme(axis.title.x =  element_blank(), axis.title.y =  element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))

 p9 <- DotPlot(scRNA,dot.scale = 8,
        features = final.markers %>%rev() ) + coord_flip()+
  RotatedAxis() +  dottheme + 
  theme(axis.text.x = element_text(angle = 45, face="bold",color = 'black', hjust=1), axis.text.y = element_text(color = 'black',face="bold"))+ 
  scale_colour_gradientn(colours = rev(mapal))+ 
  theme(legend.position="right")

ggsave("./Fig/Celltype/Markers_DotPlot_1.pdf",plot = p9, width =6, height = 5)

p10 <- DotPlot(scRNA,dot.scale = 8,
              features = final.markers %>%rev() ) + coord_flip()+
  RotatedAxis() +  dottheme + 
  theme(axis.text.x = element_text(angle = 45, face="bold",color = 'black', hjust=1), axis.text.y = element_text(color = 'black',face="bold"))+ 
  scale_colour_gradientn(colours = pal)+ 
  theme(legend.position="right")

ggsave("./Fig/Celltype/Markers_DotPlot_2.pdf",plot = p10, width =6, height = 5)




###############

##===结果可视化===##
##Stackbar细胞丰度柱状图
library(dplyr)

table(scRNA$celltype)

tmp <- dplyr::select(scRNA@meta.data, c("orig.ident", "celltype"))

df <- data.frame()
for(i in unique(tmp$orig.ident)){
  df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(celltype) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("sample", "celltype", "value")
  df <- rbind(df, df_i)
}
#source("../Resource/Resource/sc_function.R")

#按样本统计细胞类型

ylevels<-  levels(df$celltype)

xlevels<-c('Control_1', 'Control_2','Control_3', 'pSS_1', 'pSS_2','pSS_3')
library(ggplot2)

p10 <- ggplot(df, aes(x=sample, y=value, fill=celltype)) +
  geom_bar(stat= "identity", position = "fill") +
  scale_fill_manual(values = color.cluster) +
  scale_x_discrete(limits=xlevels)+
  scale_y_discrete(limits=ylevels)+
  theme_bw()+
  theme(plot.title = element_blank())+
  theme(panel.grid =element_blank())+ # 删去网格线
  theme(axis.title.x=element_blank())+
  #theme(axis.title.y=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", 
                                    size=1.0, linetype="solid"))+
  scale_y_continuous(labels=scales::percent)+#变为百分比
  theme(legend.text = c())+#图例字体
  scale_y_continuous(labels=scales::percent)+
  labs(y = 'Relative Abundance') +
  #坐标轴字体
  theme(axis.text.x = element_text(size = 10, #family = "Arial",#不支持字体
                                   color = "black", 
                                   # face = "bold", 
                                   hjust = 1, 
                                   angle = 45))+ 
  theme(axis.text.y = element_text(size = 10, #family = "Arial",
                                   color = "black"))+
  labs(y = 'Relative Abundance') +
  #labs(x = 'Samples', y = 'Relative Abundance', title = 'celltypes composition') +
  guides(fill=guide_legend(title='celltype'))#图例标题

p10
ggsave('./Fig/Celltype/Stackbar_celltype.png', p10, width = 8, height = 6)
ggsave('./Fig/Celltype/Stackbar_celltype.pdf', p10, width = 8, height = 6)

p11<- ggplot(df, aes(x=celltype, y=value, fill=sample)) +
  geom_bar(stat= "identity", position = "fill") +
  scale_fill_manual(values = color.sample) +
  scale_y_discrete(limits=xlevels)+
  scale_x_discrete(limits=ylevels)+
  theme_bw()+
  theme(plot.title = element_blank())+
  theme(panel.grid =element_blank())+ # 删去网格线
  theme(axis.title.x=element_blank())+
  #theme(axis.title.y=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", 
                                    size=1.0, linetype="solid"))+
  scale_y_continuous(labels=scales::percent)+#变为百分比
  theme(legend.text = c())+#图例字体
  #坐标轴字体
  theme(axis.text.x = element_text(size = 10, #family = "Arial",#不支持字体
                                   color = "black", 
                                   # face = "bold", 
                                   hjust = 1, 
                                   angle = 45))+ 
  theme(axis.text.y = element_text(size = 10, #family = "Arial",
                                   color = "black"))+
  labs(y = 'Relative Abundance') +
  #labs(x = 'Samples', y = 'Relative Abundance', title = 'celltypes composition') +
  guides(fill=guide_legend(title='samples'))#图例标题
# theme_prism(base_size =10)+

p11

ggsave('./Fig/Celltype/Stackbar_sample.png', p11, width = 8, height = 6)
ggsave('./Fig/Celltype/Stackbar_sample.pdf', p11, width = 8, height = 6)


library(RColorBrewer)
mytheme <-theme_bw() +
  theme(panel.grid =element_blank(),
        panel.border = element_blank()) + ## 删去网格线
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

p12<- FeaturePlot(scRNA,
                   reduction="umap", 
                   features = checkgene,
                   ncol = 5)&theme(legend.position = "right")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))&
  mytheme
p12
ggsave("./Fig/Celltype/celltype_Markers_FeaturePlot.png",plot = p12, width =12, height = 6)
ggsave("./Fig/Celltype/celltype_Markers_FeaturePlot.pdf",plot = p12, width =12, height =6)


Idents(scRNA ) <- "celltype" 
table(scRNA $celltype)

Myeloid <- subset(scRNA, idents = "Myeloid")
saveRDS(Myeloid, file = './data/Myeloid.rds')

T_cells <- subset(scRNA, idents = "T_cells")
saveRDS(T_cells, file = './data/T_cells.rds')

B_cells <- subset(scRNA, idents = "B_cells")
saveRDS(B_cells, file = './data/B_cells.rds')

NK_cells <- subset(scRNA, idents = "NK_cells")
saveRDS(NK_cells, file = './data/NK_cells.rds')

Platelets <- subset(scRNA, idents = "Platelets")
saveRDS(Platelets, file = './data/Platelets.rds')