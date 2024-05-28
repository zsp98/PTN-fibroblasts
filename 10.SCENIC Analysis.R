library(Seurat)
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(SingleCellExperiment)
library(tidyverse)
library(dplyr)
load("./data/fibro_subcelltype.Rdata")
#提取出上调/下调的gene
results_list <- list()

for (cell_type in c("apCAF", "MesCAF", "MyCAF", "vCAF")) {
  markers <- FindMarkers(sce_fibro, ident.1 = cell_type, ident.2 = "iCAF", only.pos = FALSE)
  results_list[[cell_type]] <- markers
}

# 给results_list中的每个细胞类型的差异基因结果加一列基因名
for (cell_type in c("apCAF", "MesCAF", "MyCAF", "vCAF")) {
  # 获取差异基因结果
  markers <- results_list[[cell_type]]
  
  # 将行名（基因名）添加为新的列，并命名为“gene”
  markers$gene <- rownames(markers)
  
  # 更新results_list中的差异基因结果
  results_list[[cell_type]] <- markers
}



apCAF<-results_list[["apCAF"]]
MesCAF<-results_list[["MesCAF"]]
MyCAF<-results_list[["MyCAF"]]
vCAF<-results_list[["vCAF"]]
# 创建一个空的列表来存储每种细胞类型的上调和下调基因列表
upregulated_genes <- list()
downregulated_genes <- list()

# 遍历之前的结果列表
for (cell_type in c("apCAF", "MesCAF", "MyCAF", "vCAF")) {
  # 提取每种细胞类型相对于ident.2的差异基因列表
  markers <- results_list[[cell_type]]
  up_genes <- markers$gene[markers$avg_log2FC > 0]
  down_genes <- markers$gene[markers$avg_log2FC < 0]
  
  # 存储每种细胞类型的上调和下调基因列表
  upregulated_genes[[cell_type]] <- up_genes
  downregulated_genes[[cell_type]] <- down_genes
}

# 对每种细胞类型的上调和下调基因列表进行交集操作
install.packages("ggVennDiagram")
##加载r包
library(ggplot2)
library(ggVennDiagram)
# 使用 ggVennDiagram 函数创建基因集的 Venn 图
ggVennDiagram(upregulated_genes) + 
  
  # 使用 scale_fill_gradient 函数设置图的填充颜色为渐变色，从灰色到红色
  scale_fill_gradient(low = "grey90",high = "red")
ggVennDiagram(downregulated_genes) + 
  
  # 使用 scale_fill_gradient 函数设置图的填充颜色为渐变色，从灰色到红色
  scale_fill_gradient(low = "grey90",high = "blue")


intersection_up <- Reduce(intersect, upregulated_genes)
intersection_down <- Reduce(intersect, downregulated_genes)
write.csv(intersection_up,file = "./fibro_sucelltype/code/intersection_up.txt")
write.csv(intersection_down,file = "./fibro_sucelltype/code/intersection_down.txt")

# Create a data frame with two columns using the two vectors
gene<-union(intersection_up, intersection_down)
save(gene,file = "./fibro_sucelltype/data/Alter_gene.Rdata")

#提取出表达矩阵及注释信息
exprMat <- as.matrix(sce_fibro@assays$RNA@counts)
exprMat <-exprMat[gene,]
dim(exprMat)
cellInfo <- as.data.frame(sce_fibro@meta.data)
#构建scenicOptions 对象
myDatasetTitle <- "Fibroblasts SCENIC" # choose a name for your analysis
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
scenicOptions <- initializeScenic(org="hgnc",#mouse填'mgi', human填'hgnc',fly填'dmel') 
                                  dbDir="./SCI/SCENIC/data/hg19",
                                  datasetTitle=myDatasetTitle,
                                  nCores=8)#这里可以设置并行计算
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
###计算共表达网络

#加载细胞通讯得出的差异基因
dim(exprMat)
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)#筛选出267个基因
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
# 基因共表达网络计算

mymethod <- 'runGenie3' # 'grnboost2'
#mymethod <- 'runGenie3' # 'grnboost2'
library(reticulate)
if(mymethod=='runGenie3'){
  runGenie3(exprMat_filtered_log, scenicOptions)
}else{
  arb.algo = import('arboreto.algo')
  tf_names = getDbTfs(scenicOptions)
  tf_names = Seurat::CaseMatch(
    search = tf_names,
    match = rownames(exprMat_filtered))
  adj = arb.algo$grnboost2(
    as.data.frame(t(as.matrix(exprMat_filtered))),
    tf_names=tf_names, seed=2023L
  )
  colnames(adj) = c('TF','Target','weight')
  saveRDS(adj,file=getIntName(scenicOptions,
                              'genie3ll'))
}
#3.2. 构建并计算基因调控网络活性(GRN score)
exprMat_log <-log2(exprMat+1) 
#exprMat_log <- log2(exprMat+1)
#scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
#3.3. 一些下游的可视化与探索
tsneAUC(scenicOptions, aucType="AUC") #利用AUCell score运行tsne
saveRDS(scenicOptions, file="./fibro_sucelltype/SCENIC explore/scenicOptions.Rds")#保存结果 
#if(!file.exists('output/scenic.loom')){export2loom(scenicOptions, exprMat_up)}
# motif富集的相关信息在：
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
unique(motifEnrichment_selfMotifs_wGenes$highlightedTFs)
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="FOS"]#这样查看Sox8的信息
viewMotifs(tableSubset) 
#查看motif和对应基因
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
unique(regulonTargetsInfo$TF)
#[1] "FOS"  "FOSB" "JUN"  "JUNB"
tableSubset <- regulonTargetsInfo[TF=="FOS" & highConfAnnot==TRUE]
viewMotifs(regulonTargetsInfo) 
viewMotifs(tableSubset)
#计算并展示每种细胞类型特异性的(celltype Cell-type specific regulators (RSS))，
#用热图探究细胞与regulon的关系
#AUCell score
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
# 去除extened regulons
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seurat_clusters),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,as.character(cells)]))#按照每种celltype求regulon活性均值
# Scale函数是对列进行归一化，所以要把regulonActivity_byGroup转置成细胞为行，基因为列
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))#scale用于绘制热图
dim(regulonActivity_byCellType_Scaled)
regulonActivity_byCellType_Scaled<-na.omit(regulonActivity_byCellType_Scaled)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")#绘制热图

#输出调控基因活性
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)#转换为长格式
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]#查看motif
viewTable(topRegulators)#输出列表

#看激活比例
minPerc <- 0.7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$seurat_clusters), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]#这里的minPerc不知如何描述,把每种细胞中对应的转录因子激活比例求和，至少有一种细胞大于0.7regulon的被保留
ComplexHeatmap::Heatmap(binaryActPerc_subset, name="Regulon activity (%)", col = c("white","pink","red"))#取minPerc后的数据

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Binarized, name="Regulon activity (%)", col = c("white","pink","red"))#未取minPerc的数据
#计算Regulon Specificity Score (RSS)评估细胞特异性
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
# 去除extened regulons
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "seurat_clusters"])
rss <- na.omit(rss) 
rssPlot <- plotRSS(rss)
## RSS > 0.01被用于可视化
plotly::ggplotly(rssPlot$plot)
plotRSS_oneSet(rss, setName = "interneurons")
