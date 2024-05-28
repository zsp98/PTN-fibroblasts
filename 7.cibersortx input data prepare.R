library(Seurat)
library(SeuratObject)
library(dplyr)
#######################单细胞矩阵准备#############################
load("./fibro_sucelltype/data/fibro_subcelltype.Rdata")
sce_fibro$seurat_clusters <- ifelse(sce_fibro$seurat_clusters == "iCAF", "PTN+ Fibroblasts", "PTN- Fibroblasts")
dd<-unique(sce_fibro$seurat_clusters)
dd<-as.character(dd)
seuratObj <- sce_fibro#加载数据
#由于单细胞数据太大，而Cibersortx最多只支持1Gb的文件输入
#因此我们对每一个细胞进行抽样
set.seed(1234)
############ 1、从SeuratObject中抽样 #############
exp_matrix = seuratObj@assays$RNA@counts#表达矩阵
metadata = seuratObj@meta.data#metadata
identical(row.names(seuratObj@meta.data),colnames(exp_matrix))#TRUE
##开始抽样
cellnames = as.character(dd)
sample_name = c()
nsample = 150#请输入抽样数
####注意：样本量不够时请设置replace = T
#分层抽样
for (i in 1:length(cellnames)) {
  newmetadata = metadata[metadata$seurat_clusters == cellnames[i],]
  sample_name = c(sample_name,
                  row.names(newmetadata[sample(1:nrow(newmetadata),nsample,replace = F),]))
  sample_name = gsub("\\.[0-9]","",sample_name)
  rm(newmetadata)
}

############# 2. 获得矩阵#############
exp_matrix = seuratObj@assays$RNA@counts#表达矩阵
cell_names = metadata[sample_name,]$seurat_clusters#细胞名字
exp_matrix = exp_matrix[,sample_name]#提取抽样的矩阵
colnames(exp_matrix) = cell_names#细胞名命名矩阵列明
exp_matrix = as.data.frame(exp_matrix)

#colnames(exp_matrix) <- ifelse(colnames(exp_matrix) == "iCAF", "Fibro_Cluster_1",
#                               ifelse(colnames(exp_matrix) %in% c("MesCAF", "MyCAF", "vCAF"), "Fibro_Cluster_2", "Fibro_Cluster_3"))
#############
load("~/miniconda3/r4.3/SCI/cellchat/data/sce2.Rdata")
dd<-unique(Idents(sce))
dd<-as.character(dd)
dd<-dd[dd != "Fibroblasts"]
sce <- subset(sce, subset = singleR %in% dd)
seuratObj <-sce

set.seed(1234)
############ 1、从SeuratObject中抽样 #############
sce_exp_matrix = seuratObj@assays$RNA@counts#表达矩阵
sce_metadata = seuratObj@meta.data#metadata
identical(row.names(seuratObj@meta.data),colnames(sce_exp_matrix))#TRUE
##开始抽样
cellnames = as.character(dd)
sample_name = c()
nsample = 150#请输入抽样数
####注意：样本量不够时请设置replace = T
#分层抽样
for (i in 1:length(cellnames)) {
  newmetadata = sce_metadata[sce_metadata$singleR == cellnames[i],]
  sample_name = c(sample_name,
                  row.names(newmetadata[sample(1:nrow(newmetadata),nsample,replace = F),]))
  sample_name = gsub("\\.[0-9]","",sample_name)
  rm(newmetadata)
}

############# 2. 获得矩阵#############
sce_exp_matrix = seuratObj@assays$RNA@counts#表达矩阵
cell_names = sce_metadata[sample_name,]$singleR#细胞名字
sce_exp_matrix = sce_exp_matrix[,sample_name]#提取抽样的矩阵
colnames(sce_exp_matrix) = cell_names#细胞名命名矩阵列明
sce_exp_matrix = as.data.frame(sce_exp_matrix)

############## 3. 输出文件 #############
identical(rownames(exp_matrix),rownames(sce_exp_matrix))
final<-cbind(exp_matrix,sce_exp_matrix)
write.table(final,file = "./SCI/Cibersortx analysis/3.12 scRNA_assay_fibrocluster.txt",sep = "\t",
            row.names = TRUE, col.names = NA)
save(sce,sce_fibro,final,file = "./SCI/Cibersortx analysis/Sc reference data.Rdata")
