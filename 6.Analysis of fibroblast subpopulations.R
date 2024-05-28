library(Seurat)
library(tidyverse)
load("./SCI/cellchat/data/sce2.Rdata")
unique(sce$singleR)
Idents(sce)<-sce$singleR
sce_fibro<-subset(sce,idents = "Fibroblasts")

##降维
sce_fibro <- RunPCA(sce_fibro, features = VariableFeatures(object = sce_fibro))
# Examine and visualize PCA results a few different ways
print(sce_fibro[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sce_fibro, dims = 1:2, reduction = "pca")
DimPlot(sce_fibro, reduction = "pca") + NoLegend()
DimHeatmap(sce_fibro, dims = 1:5, cells = 500, balanced = TRUE)
ElbowPlot(sce_fibro)
#Cluster the cells(重新聚类)
sce_fibro <- FindNeighbors(sce_fibro, dims = 1:15)
sce_fibro <- FindClusters(sce_fibro, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(sce_fibro), 5)
sce_fibro <- RunUMAP(sce_fibro, dims = 1:10)
DimPlot(sce_fibro, reduction = "umap")
my9color<-c("#D52126","#88CCEE","#FEE52C","#03c383","#CC61B0")
DimPlot(sce_fibro,group.by = "seurat_clusters",split.by = 'Class',cols =my9color,label = TRUE)  + NoAxes() + theme(plot.title = element_blank() )
DimPlot(sce_fibro,group.by = "seurat_clusters",cols =my9color,label = TRUE,label.size = 5,reduction = "umap") + theme(plot.title = element_blank() )

## find all markers of cluster 2
markers <- FindAllMarkers(sce_fibro, only.pos = TRUE)
markers<-markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


color <-c("lightgrey","blue","seagreen2")
gene<-c("APOE","FBLN1",#Human iCAF(富含炎症反应和生长因子转录本的炎性CAF)marker gene
        "MCAM","CRIP1","TINAGL1",#vCAF 富含参与血管发育德转录本的血管CAF marker gene
        "POSTN","CTHRC1",#MyCAF 富含ECM转录本的肌肉成纤维细胞 marker gene
        "KRT8","KRT18","KRT19",#MesCAF 表达与EMT相关德转录本的EMT样CAF marker gene
        "HLA-DRB1","HLA-DRA")#apCAF 富含MHC II转录本的抗原提呈CAF marker gene
features<-gene
p1<-FeaturePlot(sce_fibro,features =features,cols =color,pt.size = 1)+
  theme(panel.border = element_rect(fill= NA,color = "black",size = 1,linetype = "solid"))
p1

library(dplyr)

cluster <- data.frame(colnames(sce_fibro), sce_fibro$seurat_clusters)
names(cluster) <- c("Sample", "seurat_clusters")

cluster <- cluster %>%
  mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  mutate(seurat_clusters = case_when(
    seurat_clusters %in% c("1", "3", "4", "5") ~ "iCAF",
    seurat_clusters %in% c("2", "7") ~ "vCAF",
    seurat_clusters %in% c("0", "6", "8") ~ "MyCAF",
    seurat_clusters %in% c("10") ~ "MesCAF",
    seurat_clusters %in% c("9") ~ "apCAF",
    TRUE ~ as.character(seurat_clusters)
  ))
class(sce_fibro$seurat_clusters)
class(cluster1$seurat_clusters)
# 用trimws()函数去除字符首尾空格和不可见字符
cluster$seurat_clusters <- trimws(cluster$seurat_clusters)
# 将字符型 vector 转换为 factor
cluster$seurat_clusters <- as.factor(cluster$seurat_clusters)
# 再次检查 class
class(cluster$seurat_clusters)
table(cluster$seurat_clusters)
sce_fibro$seurat_clusters<-cluster$seurat_clusters

## find all markers of cluster 2
Idents(sce_fibro)<-sce_fibro$seurat_clusters
my9color <- c('#5470c6','#91cc75','#fac858','#ee6666','#73c0de','#3ba272','#fc8452','#9a60b4','#ea7ccc')
DimPlot(sce_fibro,group.by = "seurat_clusters",cols =my9color,label = TRUE)  + NoAxes() + theme(plot.title = element_blank() )

markers <- FindAllMarkers(sce_fibro, only.pos = TRUE)
markers<-markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 2)

save(markers,file = "./fibro_sucelltype/markers.Rdata")


gene<-c("APOE","FBLN1",#Human iCAF(富含炎症反应和生长因子转录本的炎性CAF)marker gene
        "MCAM","CRIP1","TINAGL1",#vCAF 富含参与血管发育德转录本的血管CAF marker gene
        "POSTN","CTHRC1",#MyCAF 富含ECM转录本的肌肉成纤维细胞 marker gene
        "KRT8","KRT18","KRT19",#MesCAF 表达与EMT相关德转录本的EMT样CAF marker gene
        "CD74","HLA-DRB1","HLA-DRA")#apCAF 富含MHC II转录本的抗原提呈CAF marker gene

p3 <- DotPlot(sce_fibro , features = gene) +
  coord_flip() +
  xlab("Genes") +
  ylab("Cell_type") +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    text = element_text(size = 16),
    axis.text.x = element_text(size = 14)  # 将y轴上的字体大小设置为14
  ) +
  ggtitle("GSE132465 Fibroblasts subcelltype")
p3

save(sce_fibro,file="./fibro_sucelltype/fibro_subcelltype.Rdata")

####可视化
rm(list = ls());gc()
devtools::install_github("junjunlab/ClusterGVis")

library(ClusterGVis)
library(org.Hs.eg.db)
library(Seurat)
load("./fibro_sucelltype/data/fibro_subcelltype.Rdata")
fibro.markers <- Seurat::FindAllMarkers(sce_fibro,
                                        only.pos = TRUE,
                                        min.pct = 0.25,
                                        logfc.threshold = 0.25)

# get top 20 genes
fibro.markers_top20 <- fibro.markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 30, wt = avg_log2FC)

# check
head(fibro.markers_top20)
# prepare data from seurat object
st.data <- prepareDataFromscRNA(object = sce_fibro,
                                diffData = fibro.markers_top20,
                                showAverage = TRUE,
                                keep.uniqGene = FALSE,
                                sep = "_")

str(st.data)
# enrich for clusters
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 123456)
# check
head(enrich)

# add gene name
markGenes = unique(fibro.markers_top20$gene)[sample(1:length(unique(fibro.markers_top20$gene)),40,
                                                    replace = F)]

# line plot
visCluster(object = st.data,
           plot.type = "line")
#配置颜色
color<-jjAnno::useMyCol("stallion", n =5)
#[1] "#D51F26" "#272E6A" "#208A42" "#89288F" "#F47D2B"
# 设置需要使用的颜色名称
colors<-c("#D51F26","#D51F26","#D51F26","#D51F26","#D51F26",
          "#272E6A","#272E6A","#272E6A","#272E6A","#272E6A",
          "#208A42","#208A42","#208A42","#208A42","#208A42",
          "#89288F","#89288F","#89288F","#89288F","#89288F",
          "#F47D2B","#F47D2B","#F47D2B","#F47D2B","#F47D2B")
# heatmap plot
pdf('./fibro_sucelltype/plot/sc1.pdf',height = 10,width = 6,onefile = F)
visCluster(object = st.data,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           cluster.order = c(1:5))
dev.off()
# add bar plot
pdf('./fibro_sucelltype/plot/sc2.pdf',height = 10,width = 14,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 25,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = 1:5,
           go.col = colors,
           add.bar = TRUE)
dev.off()


