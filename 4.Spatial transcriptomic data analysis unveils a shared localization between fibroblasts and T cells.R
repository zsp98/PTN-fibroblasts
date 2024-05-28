library(devtools)
install_github("navinlabcode/CellTrek")
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readr)
library(CellTrek)
library(viridis)
library(ConsensusClusterPlus)
library(SeuratObject)
library(jsonlite)
library(png)
library(tidyverse)
library(ggpubr)


##根据链接下载空间转录组数据https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225857
data.dir = "./spatial analysis/GSM7058759_C4"
source("./spatial analysis/code/Load10X_Spatial_change.R")
T1 <- Load10X_Spatial_change(data.dir = data.dir,
                             filename = "filtered_feature_bc_matrix")
for (i in colnames((T1@images$slice1@coordinates))) {
  T1@images$slice1@coordinates[[i]] <- as.integer(T1@images$slice1@coordinates[[i]])
}
p1 <- SpatialDimPlot(T1,alpha = 0)
p2 <- SpatialFeaturePlot(T1, features = "nFeature_Spatial",
                         pt.size.factor = 3)
p1+p2

#数据预处理（过滤，NormalizeData）
#1）计算线粒体/核糖体/特定基因集的比例

mt.genes <- grep(pattern = "^MT-", x = rownames(T1), value = TRUE)
T1$percent.mito <- (Matrix::colSums(T1@assays$Spatial@counts[mt.genes, ])/Matrix::colSums(T1@assays$Spatial@counts))*100
#2）过滤/保留基因
#保留在所有spot中的表达量之和大于5的基因
genes_to_keep <- setdiff(names(which(Matrix::rowSums(T1@assays$Spatial@counts )>5)),mt.genes)
#subset过滤spot
{
  T1$percent.mito <- (Matrix::colSums(T1@assays[["Spatial"]]@counts[mt.genes, ])/Matrix::colSums(T1@assays[["Spatial"]]@counts))*100
  #2）过滤/保留基因
  #保留在所有spot中的表达量之和大于5的基因
  genes_to_keep <- setdiff(names(which(Matrix::rowSums(T1@assays[["Spatial"]]@counts)>5)),mt.genes)
}#Seurat则这样处理
T1 <- subset(T1,
             features =genes_to_keep,  #针对基因
             subset = nFeature_Spatial > 300 & percent.mito < 30 #针对spots
)
plot1 <- VlnPlot(T1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(T1, features = "nCount_Spatial",pt.size.factor = 3) + 
  theme(legend.position = "right")
wrap_plots(plot1, plot2)
load("./spatial analysis/data/2.01_sce_SCT.Rdata")
#load("./spatial analysis/data/2.1_sce_SCT.Rdata")
## Rename the cells/spots with syntactically valid names
T1 <- RenameCells(T1, new.names=make.names(Cells(T1)))

sce <- RenameCells(sce, new.names=make.names(Cells(sce)))
## Visualize the ST data
SpatialDimPlot(T1)
## Visualize the scRNA-seq data
#Idents(sce)<-sce$label
Idents(sce)<-sce$singleR
DimPlot(sce, label = T, label.size = 4.5)
#### 2. 使用 CellTrek 进行细胞映射
T1[["RNA"]]<-T1[["Spatial"]]
# We first co-embed ST and scRNA-seq datasets using traint
COAD_traint <- CellTrek::traint(st_data=T1, 
                                sc_data=sce, 
                                sc_assay='RNA', 
                                cell_names='singleR')

## 我们可以检查共嵌入的结果，以查看这两种数据模态之间是否存在重叠。
DimPlot(COAD_traint, group.by = "type") 
# 在共嵌入之后，我们可以将单个细胞映射到它们的空间位置。
# 在这里，我们使用非线性插值（intp = T，intp_lin=F）方法来增强ST的spots

COAD_celltrek <- CellTrek::celltrek(st_sc_int=COAD_traint, int_assay='traint', sc_data=sce, sc_assay = 'RNA', 
                                    reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000, 
                                    dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)$celltrek

# 在细胞映射完成后，我们可以使用 celltrek_vis 交互式可视化 CellTrek 的结果
COAD_celltrek$label <- factor(COAD_celltrek$label, 
                              levels=sort(unique(COAD_celltrek$label)))

CellTrek::celltrek_vis(COAD_celltrek@meta.data %>% 
                         dplyr::select(coord_x, coord_y, label:id_new),
                       COAD_celltrek@images$slice1@image, 
                       COAD_celltrek@images$slice1@scale.factors$lowres)
save(COAD_celltrek,COAD_traint,file = "./spatial analysis/data/2.15 T4_full_celltrek_data.Rdata")
#Step3. CellTrek细胞共定位分析
#### 3. 细胞共定位分析
# 基于 CellTrek 的结果，我们可以使用 SColoc 总结不同细胞类型之间的共定位模式。
# 在这里，我们以谷氨酸能神经元细胞类型作为示例（建议删除一些细胞类型，例如，n<20，细胞数量非常少的细胞类型）。
# 我们首先从我们的映射结果中子集化谷氨酸能神经元细胞类型
load("~/miniconda3/r4.3/spatial analysis/data/2.15 T4_full_celltrek_data.Rdata")
my9color <- c('#5470c6','#91cc75','#fac858','#ee6666','#73c0de','#3ba272','#fc8452','#9a60b4','#ea7ccc')
DimPlot(COAD_celltrek, group.by = 'singleR',cols = my9color,label = TRUE,label.size = 6)+
  theme(
    plot.title = element_blank(),
    text = element_text(size = 20)
  )

glut_cell <- unique(COAD_celltrek$singleR)
names(glut_cell) <- make.names(glut_cell)
COAD_celltrek_glut <- subset(COAD_celltrek, subset=singleR %in% glut_cell)

#然后使用scoloc进行共定位分析：
COAD_celltrek_glut$singleR <- factor(COAD_celltrek_glut$singleR,
                                     levels=glut_cell)

## 我们从图中提取最小生成树（MST）的结果。
COAD_sgraph_KL <- CellTrek::scoloc(COAD_celltrek_glut, col_cell='singleR', use_method='KL', eps=1e-50)
COAD_sgraph_KL_mst_cons <- COAD_sgraph_KL$mst_cons
rownames(COAD_sgraph_KL_mst_cons) <- colnames(COAD_sgraph_KL_mst_cons) <- glut_cell[colnames(COAD_sgraph_KL_mst_cons)]

## 然后，我们提取meta.data数据（包括细胞类型及其频率信息）。
COAD_cell_class <- COAD_celltrek@meta.data %>% dplyr::select(id=singleR) %>% unique
COAD_celltrek_count <- data.frame(freq = table(COAD_celltrek$singleR))
COAD_cell_class_new <- merge(COAD_cell_class, COAD_celltrek_count, by.x ="id", by.y = "freq.Var1")
# visualize the colocalization result
plot<-CellTrek::scoloc_vis(COAD_sgraph_KL_mst_cons,
                           meta_data=COAD_cell_class_new )
plotly::ggplotly(rssPlot$plot)
unique(COAD_celltrek$singleR)

#第二步，计算目标细胞与其他细胞的聚类：
Target_cell<-c("Endothelial cells","T cells","Epithelial cells","Monocytes","B cells")

output <- kdist(inp_df = inp_df, 
                ref = "Fibroblasts", #目标细胞类型
                ref_type = 'all', 
                que = Target_cell,  #其余的细胞
                k = 10, 
                new_name = "Fibroblasts vs Others",
                keep_nn = F)
head(output$kdist_df) #最后的结果
res <- output$kdist_df
res$barcode <- row.names(res)
inp_df$barcode <- row.names(inp_df)
res <- left_join(res, inp_df)
head(res)
#第三步，可视化：
library(ggpubr)
ggboxplot(data = res, 
          x = "cell_names",
          y = "Fibroblasts vs Others", 
          fill = "cell_names", 
          title = "K-distance to Fibroblasts cells")+ 
  stat_compare_means(method = "kruskal.test") +
  theme(plot.title = element_text(color="black",hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5), #,vjust = 0.5
        legend.position = "none") + labs(y = "Distance")


# 基于CellTrek的结果，我们可以使用SCoexp模块进一步研究感兴趣的细胞类型中的共表达模式。
# 在这里，我们将使用共识聚类(CC)方法以L5 IT单元为例。首先从图表结果中提取L5 IT单元格。
COAD_celltrek_F <- subset(COAD_celltrek, subset=singleR=='Fibroblasts')
DimPlot(COAD_celltrek_F, group.by = "label",cols = my9color,label = TRUE,label.size = 4)+
  theme(
    plot.title = element_blank(),
    text = element_text(size = 20)
  )


COAD_celltrek_F@assays$RNA@scale.data <- matrix(NA, 1, 1)
COAD_celltrek_F$cluster <- gsub('Fibroblasts VISp ', '', COAD_celltrek_F$singleR)
DimPlot(COAD_celltrek_F, group.by = 'cluster')
# 我们选择前2000个可变基因(不包括线粒体、核糖体和高零基因)
COAD_celltrek_F <- FindVariableFeatures(COAD_celltrek_F)
vst_df <- COAD_celltrek_F@assays$RNA@meta.features %>% data.frame %>% mutate(id=rownames(.))
nz_test <- apply(as.matrix(COAD_celltrek_F[['RNA']]@data), 1, function(x) mean(x!=0)*100)
hz_gene <- names(nz_test)[nz_test<20]
mt_gene <- grep('^Mt-', rownames(COAD_celltrek_F), value=T)
rp_gene <- grep('^Rpl|^Rps', rownames(COAD_celltrek_F), value=T)
vst_df <- vst_df %>% dplyr::filter(!(id %in% c(mt_gene, rp_gene, hz_gene))) %>% arrange(., -vst.variance.standardized)
feature_temp <- vst_df$id[1:2000]
# 我们使用 scoexp 进行空间加权基因共表达分析。
set.seed(1234)
COAD_celltrek_F_scoexp_res_cc <- CellTrek::scoexp(celltrek_inp=COAD_celltrek_F, 
                                                  assay='RNA', 
                                                  approach='cc', 
                                                  gene_select = feature_temp, 
                                                  sigm=140, 
                                                  avg_cor_min=.4, 
                                                  zero_cutoff=3, 
                                                  min_gen=50, max_gen=500)

# 我们可以使用热图可视化共表达式模块。
COAD_celltrek_F_k <- rbind(data.frame(gene=c(COAD_celltrek_F_scoexp_res_cc$gs[[1]]), G='K1'), 
                           data.frame(gene=c(COAD_celltrek_F_scoexp_res_cc$gs[[2]]), G='K2')) %>% 
  magrittr::set_rownames(.$gene) %>% dplyr::select(-1)
pheatmap::pheatmap(COAD_celltrek_F_scoexp_res_cc$wcor[rownames(COAD_celltrek_F_k), rownames(COAD_celltrek_F_k)], 
                   clustering_method='ward.D2', annotation_row=COAD_celltrek_F_k, show_rownames=F, show_colnames=F, 
                   treeheight_row=10, treeheight_col=10, annotation_legend = T, fontsize=8,
                   color=viridis(10), main='Fibroblasts cells spatial co-expression')                 

# 我们确定了2个不同的基因模块。基于我们识别的共表达模块，可以计算模块的分数。
COAD_celltrek_F <- AddModuleScore(COAD_celltrek_F, 
                                  features=COAD_celltrek_F_scoexp_res_cc$gs, 
                                  name='CC_', nbin=10, ctrl=50, seed=42)
## 可视化1
FeaturePlot(COAD_celltrek_F, 
            grep('CC_1', colnames(COAD_celltrek_F@meta.data), 
                 value=T), ncol = 1)
FeaturePlot(COAD_celltrek_F, 
            grep('CC_2', colnames(COAD_celltrek_F@meta.data), 
                 value=T), ncol = 1)
## 可视化3
SpatialFeaturePlot(COAD_celltrek_F, 
                   grep('CC_', colnames(COAD_celltrek_F@meta.data), value=T))


#将Fibroblast模块化，再进行分析
COAD_celltrek_F$Cluster<-COAD_celltrek_F$label
COAD_celltrek_F$Cluster <- ifelse(COAD_celltrek_F$Cluster == "iCAF", "Fibro_Cluster_1",
                                  ifelse(COAD_celltrek_F$Cluster %in% c("MesCAF", "MyCAF", "vCAF","apCAF"), "Fibro_Cluster_2",COAD_celltrek_F$Cluster))
DimPlot(COAD_celltrek_F, group.by = 'Cluster')

######将T细胞与成纤维细胞提取出来单独进行共定位分析
COAD_celltrek_TF <- subset(COAD_celltrek, subset = singleR %in% c('Fibroblasts',"T cells"))
unique(COAD_celltrek_TF$label)

COAD_celltrek_TF$Cluster<-COAD_celltrek_TF$label
unique(COAD_celltrek_TF$Cluster)
COAD_celltrek_TF$Cluster <- ifelse(COAD_celltrek_TF$Cluster == "iCAF", "Fibro_Cluster_1",
                                   ifelse(COAD_celltrek_TF$Cluster %in% c("MesCAF", "MyCAF", "vCAF","apCAF"), "Fibro_Cluster_2",
                                          ifelse(COAD_celltrek_TF$Cluster == "NK cells","T cells", COAD_celltrek_TF$Cluster)))
unique(COAD_celltrek_TF$Cluster)
glut_cell <- unique(COAD_celltrek_TF$Cluster)
names(glut_cell) <- make.names(glut_cell)
COAD_celltrek_glut <- subset(COAD_celltrek_TF, subset=Cluster %in% glut_cell)

#然后使用scoloc进行共定位分析：
COAD_celltrek_glut$Cluster <- factor(COAD_celltrek_glut$Cluster,
                                     levels=glut_cell)

## 我们从图中提取最小生成树（MST）的结果。
COAD_sgraph_KL <- CellTrek::scoloc(COAD_celltrek_glut, col_cell='Cluster', use_method='DT')
COAD_sgraph_KL_mst_cons <- COAD_sgraph_KL$mst_cons
rownames(COAD_sgraph_KL_mst_cons) <- colnames(COAD_sgraph_KL_mst_cons) <- glut_cell[colnames(COAD_sgraph_KL_mst_cons)]

## 然后，我们提取meta.data数据（包括细胞类型及其频率信息）。
COAD_cell_class <- COAD_celltrek_TF@meta.data %>% dplyr::select(id=Cluster) %>% unique
COAD_celltrek_count <- data.frame(freq = table(COAD_celltrek_TF$Cluster))
COAD_cell_class_new <- merge(COAD_cell_class, COAD_celltrek_count, by.x ="id", by.y = "freq.Var1")
# visualize the colocalization result
CellTrek::scoloc_vis(COAD_sgraph_KL_mst_cons,
                     meta_data=COAD_cell_class_new )
