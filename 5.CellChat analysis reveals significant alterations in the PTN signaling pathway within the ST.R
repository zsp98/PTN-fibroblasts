class(COAD_celltrek)
options(stringsAsFactors = FALSE)
library(CellChat)
library(Seurat)
library(tidyverse)
library(viridis)
library(RColorBrewer)
#载入数据
load("./spatial analysis/data/2.15 T4_full_celltrek_data.Rdata")
#准备输入文件
unique(COAD_celltrek$label)
COAD_celltrek$label <- ifelse(COAD_celltrek$label == "iCAF", "Fibro_Cluster_1",
                              ifelse(COAD_celltrek$label %in% c("MesCAF", "MyCAF", "vCAF","apCAF"), "Fibro_Cluster_2",
                                     ifelse(COAD_celltrek$label == "NK cells","T cells", COAD_celltrek$label)))
unique(COAD_celltrek$label)

#矩阵信息
data.input = Seurat::GetAssayData(COAD_celltrek, slot = "data", assay = "RNA") 
#meta信息
meta = COAD_celltrek@meta.data
# 空间图像信息
spatial.locs = Seurat::GetTissueCoordinates(COAD_celltrek, scale = NULL, 
                                            cols = c("imagerow", "imagecol")) 
# Scale factors and spot diameters 信息 
scale.factors = jsonlite::fromJSON(txt = 
                                     file.path("./spatial analysis/GSM7058759_C4/spatial", 'scalefactors_json.json'))
scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
)
#创建CellChat对象

cellchat <- createCellChat(object = data.input, 
                           meta = meta, 
                           group.by = "label", #前面的meta ，定义的名字是labels
                           datatype = "spatial", ###
                           coordinates = spatial.locs, 
                           scale.factors = scale.factors)
#设置参考数据库
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact")#取出相应分类用作分析数据库
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")#取出相应分类用作分析数据库
cellchat@DB <- CellChatDB.use#将数据库内容载入cellchat对象中
#cellchat@DB <-CellChatDB
#CellChat预处理
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat,features = NULL) # This step is necessary even if using the whole database
future::plan("multisession", workers = 1) #笔记本可以选1
##识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
#识别过表达配体受体对
cellchat <- identifyOverExpressedInteractions(cellchat)

#二~推断cell-cell network 
#1，推断细胞通讯网络
cellchat <- computeCommunProb(cellchat, 
                              type = "truncatedMean", trim = 0.1, 
                              distance.use = TRUE, 
                              scale.distance = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = 10)
#计算每个信号通路相关的所有配体-受体相互作用的通信结果
cellchat <- computeCommunProbPathway(cellchat)
#计算整合的细胞类型之间通信结果
cellchat <- aggregateNet(cellchat)

save(cellchat,file="~/miniconda3/r4.3/spatial analysis/data/3.22 COAD_celltrek_cellchat.Rdata")

#三 CellChat 可视化
#1，celltype之间通讯结果
load("./spatial analysis/data/3.22 COAD_celltrek_cellchat.Rdata")
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#根据使用netVisual_heatmap显示任意两个celltype之间的通讯次数（左）或总通讯强度(右)
p1 <- netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")

p2 <- netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")

p1 + p2

mat <- cellchat@net$count
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= T,edge.weight.max = max(mat), title.name =paste(rownames(mat)[i],"Number of interactions"))
}

#2，单个信号通路可视化
cellchat@netP$pathways
pathways.show <- c("PTN")
levels(cellchat@idents)   
#[1] "C9" "C3" "C5" "C4" "C0" "C7" "C6" "C8" "C1" "C2"     
vertex.receiver = c(1,2,3,6,7)  #选择的是levels(cellchat@idents) 中的
#左图中间的Target是vertex.receiver选定的细胞类型，右图是除vertex.receiver选中之外的另外的细胞类型。
netVisual_aggregate(cellchat, signaling = pathways.show,                      
                    vertex.receiver = vertex.receiver,
                    layout = "hierarchy")
#2）和弦图
#可以额外绘制空间转录组版本的和弦图，添加layout = "spatial" 。
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, 
                    layout = "circle",
                    sources.use = c(4,5),
                    targets.use = c(1,2,3,6,7))
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

##单个通路中不同配受体的贡献度
netAnalysis_contribution(cellchat, signaling = pathways.show)
# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", 
                    sources.use = c(4,5,6),
                    targets.use = c(7),
                    edge.width.max = 2, vertex.size.max = 1, 
                    alpha.image = 0.2, vertex.label.cex = 3.5)

#3） network centrality scores
#识别细胞组的信号角色（例如，占主导地位的发送器、接收器）以及主要贡献信号
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow=c(1,1))
centrality <- netAnalysis_computeCentrality(cellchat)
netAnalysis_signalingRole_network(centrality, signaling = pathways.show, 
                                  width = 8, height = 2.5, font.size = 10)
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("PTN"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("PTN"))
ht <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "all")


#3，绘制配体受体气泡图
#1）指定受体-配体细胞类型
#绘制指定受体-配体细胞类型中的全部配体受体结果的气泡图，通过sources.use 和 targets.use指定。

#指定受体-配体细胞类型
netVisual_bubble(cellchat, sources.use = c(4,5), 
                 targets.use = c(7), remove.isolate = FALSE)
#2）指定受体-配体细胞类型 且 指定通路
#同时通过signaling指定展示通路
netVisual_bubble(cellchat, sources.use = c(4,5), targets.use = c(7),                  
                 signaling = c("PTN"), remove.isolate = FALSE)
#3）ligand-receptor pair 表达
# Take an input of a ligand-receptor pair and show expression in binary
spatialFeaturePlot(cellchat, pairLR.use = "PTN_NCL", point.size = 1, 
                   do.binary = TRUE, cutoff = 0.05, enriched.only = F, 
                   color.heatmap = "Reds", direction = 1)
#通路贡献图
p1 <- netAnalysis_contribution(cellchat, signaling = "PTN",
                               title =  "PTN")#展现对特定通路的贡献程度
#批量写出通路贡献图
vertex.receiver = seq(7)
pathways.show.all <- cellchat@netP$pathways
dir.create('./spatial analysis/04.pathwat.show')
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  # netVisual(cellchat, signaling = pathways.show.all[i], 
  #           vertex.receiver = vertex.receiver, layout = "hierarchy")#直接出pdf与svg
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0('./spatial analysis/04.pathwat.show/',pathways.show.all[i], "_L-R_contribution.pdf"),
         plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}
###使用小提琴/点图绘制信号基因表达分布
library(Seurat) 
library(ggsignif)
library(ggplot2)
library(ggpubr)
# 数据读入和提取
DefaultAssay(COAD_celltrek)="RNA" 
Idents(COAD_celltrek) <- "label"
features <-c("PTN","NCL","SDC1","SDC2","SDC3","SDC4")  
expr = t(as.matrix(COAD_celltrek@assays$RNA@data[features,]))
posi= match(rownames(expr),rownames(COAD_celltrek@meta.data))
ts = cbind(expr,COAD_celltrek@meta.data[posi,c("Sample","label")])
# 批量绘制表达图
# Create a list of column names to loop through
cols <- colnames(ts)[1:6]
# Create a directory to store the plots
dir.create("./spatial analysis/plot/cellchat plot/expression/", showWarnings = FALSE)

# Loop through each column and create a plot
for (col in cols) {
  p <- ggplot(ts, aes(x = label, y = ts[[col]], fill = label))+
    geom_violin(scale = "width", alpha = 0.8, width = 0.6, size = 0.8)+
    xlab("")+
    ylab("Expression Level")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.2),
          axis.text.x = element_text(angle = 45, size = 10, vjust = 1, hjust = 1, color = "black"),
          axis.text.y = element_text(size = 10))+
    labs(fill = col) +
    ggtitle(col) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save the plot as a PDF file
  pdf(file = paste0("./spatial analysis/plot/cellchat plot/expression/", col, "_expression.pdf"))
  print(p)
  dev.off()
}
##只画成纤维细胞中的
COAD_celltrek_F<-subset(COAD_celltrek, subset=singleR=='Fibroblasts')
DefaultAssay(COAD_celltrek_F)="RNA"
Idents(COAD_celltrek_F) <- "label"
features <-c("PTN","NCL","SDC1","SDC2","SDC3","SDC4")  
features <-c("MDK","NCL","SDC2","SDC1","SDC4","LRP1","ITGS4","ITGA6","ITGB1")
expr = t(as.matrix(COAD_celltrek_F@assays$RNA@data[features,]))
posi= match(rownames(expr),rownames(COAD_celltrek_F@meta.data))
ts = cbind(expr,COAD_celltrek_F@meta.data[posi,c("Sample","label")])
# 批量绘制表达图
# Create a list of column names to loop through
cols <- colnames(ts)[1:6]
# Create a directory to store the plots
dir.create("./spatial analysis/plot/cellchat plot/expression/Only Fibro/", showWarnings = FALSE)

# Loop through each column and create a plot
for (col in cols) {
  # Perform t-test for the current column
  t_test_result <- t.test(ts[[col]] ~ ts$label)
  
  # Create the ggplot with t-test result for the current column
  p <- ggplot(ts, aes(x = label, y = .data[[col]], fill = label))+
    geom_violin(scale = "width", alpha = 0.8, width = 0.6, size = 0.8)+
    xlab("")+
    ylab("Expression Level")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.2),
          axis.text.x = element_text(angle = 45, size = 10, vjust = 1, hjust = 1, color = "black"),
          axis.text.y = element_text(size = 10))+
    labs(fill = col) +
    ggtitle(paste(col, "\n", "t-test p-value:", format(t_test_result$p.value, digits = 2))) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save the plot as a PDF file
  pdf(file = paste0("./spatial analysis/plot/cellchat plot/expression/Only Fibro/", col, "_expression.pdf"))
  print(p)
  dev.off()
}
#确定全局通信模式，探索多个细胞类型和信号通路如何协调在一起
library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
nPatterns = 2
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")
# 识别和可视化目标细胞的传入通信模式
selectK(cellchat, pattern = "incoming")
nPatterns = 7
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")
#信号网络的多重和分类学习分析
#此外，CellChat 能够量化所有重要信号通路之间的相似性，然后根据其CellChat 网络的相似性对其进行分组。分组可以基于功能或结构相似性进行。
#功能相似性：功能相似度高表示主要发送器和接收器相似，可解释为两个信号通路或两个配体受体对具有相似的作用。功能相似性分析要求两个数据集之间的细胞群组成相同。
#结构相似性：结构相似性用于比较其信号网络结构，而不考虑发送器和接收器的相似性。
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
library(future)
# 指定futures策略为sequential
plan(multiprocess)
# 运行netClustering函数
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)