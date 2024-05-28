library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(monocle)
library(readxl)
library(data.table)
library(readxl)
set.seed(20222202)

#######数据读取，聚类及注释

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465下载文件，并解压
# 读取注释文件（复制粘贴到名为cell_anno的表格中）GSE132465
cell_anno <- read_excel("cell_anno.xlsx")
head(cell_anno)
cell_anno <- as.data.frame(cell_anno)
rownames(cell_anno) <- cell_anno$Index
cell_anno$Index <- NULL
# 读取counts文件
counts <- fread("./GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz")
head(counts[1:4,1:4])
counts <- as.data.frame(counts)
rownames(counts) <- counts$Index
counts$Index <- NULL
head(counts[1:4,1:4])
#step2：创建seurat对象
sce <- CreateSeuratObject(counts = counts)
## 添加细胞注释文件
sce <- AddMetaData(sce,metadata = cell_anno)

head(sce@meta.data)
#--------------------------------- 5. 看看数据质量 -----------------------------
sce <- PercentageFeatureSet(sce,pattern = '^MT-',col.name = 'percent.mt')
fivenum(sce@meta.data$percent.mt)
head(sce@meta.data)
# 可视化基因和UMI(可以看到作者已经对基因进行了过滤: < 6000)
features <- c('nCount_RNA','nFeature_RNA')
VlnPlot(sce,features = features,pt.size = 0,group.by = 'orig.ident')&NoLegend()&labs(x = '')&theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
# 可视化线粒体基因(可以看到作者已经对线粒体基因进行了过滤: < 20%)
features <- c('percent.mt')
VlnPlot(sce,features = features,pt.size = 0,group.by = 'orig.ident')&NoLegend()&labs(x = '')&theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

# save(sce,file = 'sce.RData')

#降维聚类
# 因为作者已经进行了质控，所以接下来我们直接进行降维聚类
sce <- NormalizeData(sce) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA(npcs = 50)

ElbowPlot(sce,ndims = 50)
sce_ha <- sce

# 不整合看效果
sce <- RunUMAP(sce,dims = 1:30)
# save(sce,file = 'sce_no_har.Rdata')

DimPlot(sce,group.by = 'orig.ident',split.by = 'orig.ident',ncol = 5) + NoLegend() + NoAxes()
DimPlot(sce,group.by = 'Cell_type',split.by = 'Class') + NoAxes()

# harmony整合
sce <- RunHarmony(sce_ha,group.by.vars = 'orig.ident',plot_convergence = T,project.dim = F)
# 降维聚类
sce <- RunUMAP(sce,dims = 1:30,reduction = 'harmony')
save(sce,file = 'sce_har.Rdata')

#####细胞重新注释#####
rm(list = ls());gc()
load("./sce_har.Rdata")
Idents(sce)<-sce$Cell_subtype
# Specify genes  
genes = c("PTPRC")
# All on Dotplot 
p1 <- DotPlot(sce, features = genes) + coord_flip()
p1
# Annotate Immune vs Nonimmune clusters
# At this point we dont care for a more detailed annotation as we will annotate immune and non-immune separately later
dat=p1$data 
fivenum(cd45$avg.exp.scaled)
imm <- cd45[cd45$avg.exp.scaled > -0.6,]$id
imm
sce@meta.data$immune_annotation <- ifelse(sce@meta.data$Cell_subtype %in% imm, 'immune', 'non-immune')
# MAke a table 
table(sce@meta.data$immune_annotation)
#immune non-immune 
#40262      23427 
# Make and save relevant plots
p <- DimPlot(sce, group.by = 'immune_annotation')
p 
#将免疫细胞进行再分类
phe <- sce@meta.data
save(phe,file = './phe-of-immune-or-not.Rdata')

sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.5)
table(sce@meta.data$RNA_snn_res.0.5)  
table(phe$immune_annotation)
cells.use <- row.names(sce@meta.data)[which(phe$immune_annotation=='immune')]
length(cells.use)
sce_immu <-subset(sce, cells=cells.use)  
sce_immu
#SingleR自动注释
sce_for_SingleR <- GetAssayData(sce_immu, slot="data")
sce_for_SingleR
library(BiocParallel)
library(SingleR)
library(celldex)
#hpca.se <- HumanPrimaryCellAtlasData()
load("F:/biomatics/SC RNA/SC RNA/GSE132465 human colorectal/Single R annotation data/DatabaseImmuneCellExpressionData.Rda")#网太慢直接读取
DatabaseImmuneCellExpressionData
##开始注释
pred.DICE <- SingleR(test = sce_for_SingleR, ref = DatabaseImmuneCellExpressionData, labels = DatabaseImmuneCellExpressionData$label.main)
table(pred.DICE$labels)
celltype <- data.frame(ClusterID=rownames(pred.DICE), celltype=pred.DICE$labels, stringsAsFactors = F) 
clusters <- colnames(sce_immu)
sce_immu@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
DimPlot(sce_immu, reduction = "umap", group.by = "singleR",label = T)
DimPlot(sce_immu, reduction = "umap", group.by = "Cell_type")
phe <- sce_immu@meta.data
table(phe$singleR)
write.csv(phe,file = "./immune.txt")#wps合并小的cluster
phe2 <- read_excel("./2 immune.xls")
fix(phe2)
phe2<-as.data.frame(phe2)
row.names(phe2)<-phe2$...1
phe3<-phe2[,-1]
test.sce<-sce_immu
test.sce@meta.data<-phe3
DimPlot(test.sce,reduction = "umap", group.by = "singleR",label =T)
phe<-phe3
sce_immu <- test.sce
rm(phe2,phe3)
table(phe$singleR)
save(sce_immu,phe,file = './phe-of-subtypes-Immune-by-singleR.Rdata')

######################
#将非免疫细胞进行再分类
phe <- sce@meta.data
cells.use <- row.names(sce@meta.data)[which(phe$immune_annotation=='non-immune')]
length(cells.use)
sce_non_immu <-subset(sce, cells=cells.use)  
sce_non_immu
#SingleR自动注释
sce_for_SingleR <- GetAssayData(sce_non_immu, slot="data")
sce_for_SingleR
library(BiocParallel)
library(SingleR)
library(celldex)
##
load("F:/biomatics/SC RNA/SC RNA/GSE132465 human colorectal/Single R annotation data/SingleR_BlueprintEncode_bpe.se_human.RData")#网太慢直接读取
bpe.se
#
#pred.hesc <- SingleR(test = sce_for_SingleR, ref = hpca.se, labels = hpca.se$label.main)
#table(pred.hesc$labels)
pred.bpe <- SingleR(test = sce_for_SingleR, ref = bpe.se, labels = bpe.se$label.main)
table(pred.bpe$labels)
celltype <- data.frame(ClusterID=rownames(pred.bpe), celltype=pred.bpe$labels, stringsAsFactors = F) 
clusters <- colnames(sce_non_immu)
sce_non_immu@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
DimPlot(sce_non_immu, reduction = "umap", group.by = "singleR",label = T)
DimPlot(sce_non_immu, reduction = "umap", group.by = "Cell_type",label = T)
phe=sce_immu@meta.data
#write.csv(phe,file = "./GSE132465 前处理/data/explore immune-noimmue/non_immune.txt")#wps合并小的cluster
phe2 <- read_excel("./GSE132465 前处理/data/explore immune-noimmue/2 non_immune.xls")
phe2<-as.data.frame(phe2)
row.names(phe2)<-phe2$...1
phe3<-phe2[,-1]
test.sce<-sce_non_immu
test.sce@meta.data<-phe3
DimPlot(test.sce,reduction = "umap", group.by = "singleR",label =T)
phe<-phe3
sce_non_immu <- test.sce
rm(phe2,phe3,test.sce,tesr.sce);gc()
non_immu_phe<-phe
table(phe$singleR)
save(sce_non_immu,non_immu_phe,file = './phe-of-subtypes-noImmune-by-singleR.Rdata')
##########
#将免疫和非免疫子集重新合并,并将新的metadata赋值到原sce中
sce_im_nonimm <- merge(sce_immu, sce_non_immu,all = TRUE)
DimPlot(sce, reduction = "umap", group.by = "Cell_type",label = T)
DimPlot(sce_im_nonimm,reduction = "umap", group.by = "Cell_type",label = T)#可以看到合并之后是没有降维的
# 提取出两个数据集的metedata，进行匹配
phe<-sce@meta.data
phe1<-sce_im_nonimm@meta.data
clusters=colnames(sce)
celltype = data.frame(ClusterID=colnames(sce_im_nonimm), celltype=sce_im_nonimm$singleR, stringsAsFactors = F)
phe$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']#在celltype中搜索clusters中的名称，并返回对应的celltype名称。有点类似于WPS中的vlookup函数。
sce@meta.data<-phe
my9color <- c('#5470c6','#91cc75','#fac858','#ee6666','#73c0de','#3ba272','#fc8452','#9a60b4','#ea7ccc')
DimPlot(sce, reduction = "umap", group.by = "Cell_type",label = T)
DimPlot(sce, reduction = "umap", group.by = "singleR",label = T)
DimPlot(sce,group.by = 'singleR',split.by = 'Class',cols =my9color,label = TRUE)  + NoAxes() + theme(plot.title = element_blank())
save(sce,file = "././GSE132465 前处理/data/explore immune-noimmue/sce_immun_nonmun.Rdata")

##对单核细胞重新注释
load("./sce_immun_nonmun.Rdata")
Idents(sce)<-"singleR"
sce_mono<- subset(sce,idents ="Monocytes")

#hpca.se <- HumanPrimaryCellAtlasData()
load("F:/biomatics/SC RNA/SC RNA/GSE132465 human colorectal/Single R annotation data/SingleR_BlueprintEncode_bpe.se_human.RData")#网太慢直接读取
bpe.se
sce_for_SingleR <- GetAssayData(sce_mono, slot="data")
pred.bpe <- SingleR(test = sce_for_SingleR, ref = bpe.se, labels = bpe.se$label.main)
table(pred.bpe$labels)
celltype = data.frame(ClusterID=rownames(pred.bpe), celltype=pred.bpe$labels, stringsAsFactors = F) 
clusters=colnames(sce_mono)
sce_mono@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
DimPlot(sce_mono, reduction = "umap", group.by = "singleR",label = T)
phe=sce_mono@meta.data
write.csv(phe,file = "./GSE132465 前处理/data/explore immune-noimmue/mono.txt")#wps合并小的cluster
phe2 <- read_excel("./GSE132465 前处理/data/explore immune-noimmue/mono.xls")
phe2<-as.data.frame(phe2)
row.names(phe2)<-phe2$...1
phe2<-phe2[,-1]
test.sce<-sce_mono
test.sce@meta.data<-phe2

DimPlot(test.sce,reduction = "umap", group.by = "singleR",label =T)
sce_mono<-test.sce
rm(phe2,bpe.se,test.sce):gc()
#合并数据，方便下一步直接将重新注释的metadata导入到原sce文件中
Idents(sce)<-"singleR"
sce_1<- subset(sce,idents ="Epithelial cells")
sce_2<- subset(sce,idents ="Endothelial cells")
sce_3<- subset(sce,idents ="Fibroblasts")
sce_4<- subset(sce,idents ="T cells")
sce_5<- subset(sce,idents ="B cells")

sce_im_nonimm <- merge(sce_mono, sce_1,all = TRUE)
sce_im_nonimm2 <- merge(sce_2,sce_3,all = TRUE)
sce_im_nonimm3 <- merge(sce_4,sce_5,all = TRUE)
sce_12<-merge(sce_im_nonimm,sce_im_nonimm2,all = TRUE)
sce_123<-merge(sce_12,sce_im_nonimm3,all = TRUE)
rm(sce_1,sce_2,sce_3,sce_4,sce_5,sce_im_nonimm,sce_im_nonimm2,sce_im_nonimm3);gc()
dim(sce_123)
dim(sce)
phe<-sce@meta.data
phe1<-sce_123@meta.data
clusters<-colnames(sce)
celltype<-data.frame(ClusterID=colnames(sce_123), celltype=sce_123$singleR, stringsAsFactors = F)
phe$singleR<-celltype[match(clusters,celltype$ClusterID),'celltype']#在celltype中搜索clusters中的名称，并返回对应的celltype名称。有点类似于WPS中的vlookup函数。

immun_noim <-data.frame(ClusterID=colnames(sce_123), immun_noim=sce_123$immune_annotation, stringsAsFactors = F)
phe$immune_annotation=immun_noim[match(clusters,celltype$ClusterID),'immun_noim']#在celltype中搜索clusters中的名称，并返回对应的celltype名称。有点类似于WPS中的vlookup函数
sce@meta.data<-phe
my9color <- c('#5470c6','#91cc75','#fac858','#ee6666','#73c0de','#3ba272','#fc8452','#9a60b4','#ea7ccc')
DimPlot(sce, reduction = "umap", group.by = "singleR",label = T)
DimPlot(sce,group.by = 'singleR',split.by = 'Class',cols =my9color,label = TRUE)  + NoAxes() + theme(plot.title = element_blank())
DimPlot(sce,group.by = "immune_annotation",cols = my9color)  + NoAxes() + theme(plot.title = element_blank() )
save(sce,file = "./sce2.Rdata")

