load("./sce2.Rdata")
library(Seurat)
library(dplyr)
library(ggplot2)
my9color <- c('#5470c6','#91cc75','#fac858','#ee6666','#73c0de','#3ba272','#fc8452','#9a60b4','#ea7ccc')
# total umap
sce@meta.data$immune_annotation
Idents(sce)<-"immune_annotation"
table(sce@meta.data$immune_annotation)
p1 <-DimPlot(sce, group.by = "immune_annotation", cols = my9color) +
  theme(
    plot.title = element_blank(),
    text = element_text(size = 20)
  )

Idents(sce)<-"immune_annotation"
genes<-"PTPRC"
p2 <- DotPlot(sce, features = genes) +
  coord_flip() +
  xlab("Genes") +
  ylab("Group") +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    text = element_text(size = 20),
    axis.text.x = element_text(size = 14)  # 将y轴上的字体大小设置为14
  ) 


Idents(sce)<-"singleR"
DimPlot(sce, group.by = "singleR", cols = my9color,label = T,label.size = 6) +
  theme(
    plot.title = element_text(size = 20),
    text = element_text(size = 20)
  ) +
  ggtitle("Total")

Idents(sce)<-"singleR"
genes<-c("C1QB","CD68","CD79A","MZB1","CD3E","CD3D","EPCAM","KRT18","COL1A1","COL1A2","CLDN5","CDH5")
#"C1QB","CD68"为巨噬细胞marker gene，以此类推CD79A CD19 B细胞，CD3E CD3D T细胞，EPCAM 表皮细胞，COL1A1,COL1A2成纤维细胞，CLDN5内皮细胞。
p3 <- DotPlot(sce, features = genes) +
  coord_flip() +
  xlab("Genes") +
  ylab("Cell_type") +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    text = element_text(size = 20),
    axis.text.x = element_text(size = 14)  # 将y轴上的字体大小设置为14
  ) 

#绘制细胞比例图~~堆积柱状图
library(plotrix)
cellnum <- table(sce$singleR,sce$Class)
cell.prop <- as.data.frame(prop.table(cellnum))
colnames(cell.prop)<-c("Celltype","Group","Proportion")
p4 <-ggplot(cell.prop,aes(Group,Proportion,fill=Celltype))+
  geom_bar(stat = "identity",position = "fill")+
  scale_fill_manual(values = my9color)+#"自定义fill的颜色"
  ggtitle("cell_portation")+
  theme_bw()+
  theme(axis.ticks.length = unit(0.5,"cm"),
        text = element_text(size=12),
        axis.text = element_text(size=14),
        plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend(title = NULL))

p5<-DimPlot(sce, group.by = 'singleR', split.by = 'Class', cols = my9color, label = T,label.size = 6) +#调节label大小
  theme(plot.title = element_blank(),
        text = element_text(size = 16),#调节图例字体大小
        strip.text = element_text(size = 18))#调节组别字体大小
#直接划分label字体太小，以下安照组别分别画图
Idents(sce)<-sce$Class
sce_Nor<-subset(sce,idents = "Normal")
sce_tumor<-subset(sce,idents = "Tumor")
plot1<-DimPlot(sce_Nor, group.by = "singleR", cols = my9color,label = T,label.size = 6) +
  theme(
    plot.title = element_text(size = 20),
    text = element_text(size = 20)
  ) +
  ggtitle("Normal")
plot2<-DimPlot(sce_tumor, group.by = "singleR", cols = my9color,label = T,label.size = 6) +
  theme(
    plot.title = element_text(size = 20),
    text = element_text(size = 20)
  ) +
  ggtitle("Tumor")

#寻找高变基因，做热图展示
Idents(sce) <- "singleR"
marker_big <- FindAllMarkers(sce,logfc.threshold = 0.5,only.pos = T)
#write.csv(marker_big,file = './GSE132465 human colorectal/SCI/Figure1/data/marker_big.csv')
marker_big <- read.csv("./SCI/Figure1/data/step 1/marker_big.csv")
big_top100 <- marker_big %>% group_by(cluster) %>% top_n(100,avg_log2FC)
sce_big100 <- ScaleData(sce,features = big_top100$gene)
big_top2 <- marker_big %>% group_by(cluster) %>% top_n(2,avg_log2FC)
sce_big2 <- ScaleData(sce,features = big_top2$gene)

p6 <- DoHeatmap(subset(sce_big100,downsample = 100),
                features = as.character(genes), 
                group.by = "singleR", assay = "RNA",
                group.colors = c("#C77CFF","#7CAE00","#00BFC4","#F8766D","#AB82FF","#90EE90")) +  
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))+
  theme(text = element_text(size = 16))
p7 <- DoHeatmap(subset(sce_big2,downsample = 100),
                features = as.character(unique(big_top2$gene)), 
                group.by = "singleR", assay = "RNA",
                group.colors = c("#C77CFF","#7CAE00","#00BFC4","#F8766D","#AB82FF","#90EE90")) +  
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))+
  theme(text = element_text(size = 16))
p7
#利用Feature包将基因可视化
color <-c("lightgrey","blue","seagreen2")
#features<-big_top2$gene
genes<-c("C1QB","CD68","CD79A","MZB1","CD3E","CD3D","EPCAM","KRT18","COL1A1","COL1A2","CLDN5","CDH5")
#"C1QB","CD68"为巨噬细胞marker gene，以此类推CD79A CD19 B细胞，CD3E CD3D T细胞，EPCAM 表皮细胞，COL1A1,COL1A2成纤维细胞，CLDN5内皮细胞。
p8<-FeaturePlot(sce,features =genes,cols =color,pt.size = 1)+
  theme(panel.border = element_rect(fill= NA,color = "black",size = 1,linetype = "solid"))
p8
ggsave("./GSE132465 human colorectal/SCI/Figure1/plot/Rplot10.png", plot2,width = 12,
       height = 6, dpi = 600)
