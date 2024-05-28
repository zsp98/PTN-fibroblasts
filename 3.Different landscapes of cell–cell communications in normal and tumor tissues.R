##cellchat数据准备
load("./sce2.Rdata")
data.input <- sce[["RNA"]]@data # normalized data matrix
meta <- sce@meta.data # a dataframe with rownames containing cell mata data
unique(meta$Class)
meta <-select(meta,c(4:8))
cell.use <- rownames(meta)[meta$Class == "Normal"] # 按指定的变量提取细胞
data.input <- data.input[, cell.use]#取出对应细胞
meta = meta[cell.use, ]#取出对应细胞的meta信息
# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
unique(meta$Cell_subtype)
save(data.input,meta,cell.use,file = './SCI/cellchat/data/COAD_Nor_cellchat.Rdata')

rm(list = ls());gc()

load("./sce2.Rdata")
data.input <- sce[["RNA"]]@data # normalized data matrix
meta <- sce@meta.data # a dataframe with rownames containing cell mata data
unique(meta$Class)
meta <-select(meta,c(4:8))
cell.use <- rownames(meta)[meta$Class == "Tumor"] # 按指定的变量提取细胞
data.input <- data.input[, cell.use]#取出对应细胞
meta = meta[cell.use, ]#取出对应细胞的meta信息
# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
unique(meta$Cell_subtype)
save(data.input,meta,cell.use,file = './SCI/cellchat/data/COAD_Tumor_cellchat.Rdata')
#####
rm(list = ls());gc()
load("./SCI/data/COAD_Normal_cellchat.Rdata")
library(CellChat)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "singleR")
levels(cellchat@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
groupSize

#载入数据库并开始计算
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
#展示细胞组成的比例:
dplyr::glimpse(CellChatDB$interaction)#展示互作的记录
CellChatDB.use <- CellChatDB
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")#取出相应分类用作分析数据库
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use#将数据库内容载入cellchat对象中
#表达量预处理
cellchat <- subsetData(cellchat,features = NULL)#取出表达数据

cellchat <- identifyOverExpressedGenes(cellchat)#寻找高表达的基因#
#FindVariableFeatures()
cellchat <- identifyOverExpressedInteractions(cellchat)#寻找高表达的通路
cellchat <- projectData(cellchat, PPI.human)#投影到PPI
cellchat <- computeCommunProb(cellchat, raw.use = T)#默认计算方式为type = "truncatedMean",
#默认cutoff的值为20%，即表达比例在25%以下的基因会被认为是0， trim = 0.1可以调整比例阈值
cellchat <- filterCommunication(cellchat, min.cells = 10)
#去掉通讯数量很少的细胞
#df.net <- subsetCommunication(cellchat)#将细胞通讯预测结果以数据框的形式取出
#DT::datatable(df.net)
#df.net <- subsetCommunication(cellchat,slot.name = "netP")
##这种方式只取通路，数据结构更简单
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
#source可以以细胞类型的名称定义，也可以按照细胞名称中的顺序以数值向量直接取
#指定输入与输出的细胞集群
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))#指定通路提取
cellchat <- computeCommunProbPathway(cellchat)
#每对配受体的预测结果存在net中，每条通路的预测结果存在netp中
cellchat <- aggregateNet(cellchat)
#计算联路数与通讯概率，可用sources.use and targets.use指定来源与去向
save(cellchat,file = './SCI/data/2_COAD_Normal_cellchat.Rdata')

#同理整理Tumor数据

load("./SCI/data/COAD_Tumor_cellchat.Rdata")
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "singleR")
#Idents(pbmc) <- 'group'
levels(cellchat@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
groupSize

#载入数据库并开始计算
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
#展示细胞组成的比例:
dplyr::glimpse(CellChatDB$interaction)#展示互作的记录
CellChatDB.use <- CellChatDB
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")#取出相应分类用作分析数据库
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use#将数据库内容载入cellchat对象中
#表达量预处理
cellchat <- subsetData(cellchat,features = NULL)#取出表达数据

cellchat <- identifyOverExpressedGenes(cellchat)#寻找高表达的基因#
#FindVariableFeatures()
cellchat <- identifyOverExpressedInteractions(cellchat)#寻找高表达的通路
cellchat <- projectData(cellchat, PPI.human)#投影到PPI
cellchat <- computeCommunProb(cellchat, raw.use = T)#默认计算方式为type = "truncatedMean",
#默认cutoff的值为20%，即表达比例在25%以下的基因会被认为是0， trim = 0.1可以调整比例阈值
cellchat <- filterCommunication(cellchat, min.cells = 10)
#去掉通讯数量很少的细胞
#df.net <- subsetCommunication(cellchat)#将细胞通讯预测结果以数据框的形式取出
#DT::datatable(df.net)
#df.net <- subsetCommunication(cellchat,slot.name = "netP")
##这种方式只取通路，数据结构更简单
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
#source可以以细胞类型的名称定义，也可以按照细胞名称中的顺序以数值向量直接取
#指定输入与输出的细胞集群
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))#指定通路提取
cellchat <- computeCommunProbPathway(cellchat)
#每对配受体的预测结果存在net中，每条通路的预测结果存在netp中
cellchat <- aggregateNet(cellchat)
#计算联路数与通讯概率，可用sources.use and targets.use指定来源与去向
save(cellchat,file = './SCI/data/2_COAD_Tumor_cellchat.Rdata')

#######
#多组别分析
rm(list = ls());gc()
library(ggplot2)
library(CellChat)
library(patchwork)
library(cowplot)
load('./SCI/data/2_COAD_Normal_cellchat.Rdata')
cellchat.Normal <-cellchat
load('./SCI/data/2_COAD_Tumor_cellchat.Rdata')
cellchat.Tumor <- cellchat
rm(cellchat)
object.list <- list(Normal = cellchat.Normal, Tumor = cellchat.Tumor)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
save(object.list,cellchat,file = "./SCI/cellchat fibro Tcell/data/cellchat.Rdata")
#rm(cellchat.Normal,cellchat.Tumor);gc()
load("./SCI/cellchat/data/cellchat.Rdata")
#相当于Seurat对象的merge
cellchat
#最简单的展示，查看细胞互作的数量在不同条件下是否有差异
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
plot<- gg1 + gg2
ggsave("./SCI/cellchat fibro Tcell/plot/Rplot1.png", plot,width = 6,
       height = 3,dpi = 600)
#查看细胞通路在两组间的富集程度
gg3 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg4 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg3 + gg4
plot2<- gg3 + gg4
ggsave("./SCI/cellchat fibro Tcell/plot/Rplot2.png", plot2,width = 12,
       height = 6, dpi = 600)

#以circle plot的形式展示第二个组别中相较于第一个组别细胞通讯发生的变化，红色为上调蓝色为下调
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#以上是直接展示二组“相减”的结果，当然你也可以直接将两组分开展示
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


#同理，用heatmap也可以进行展示
gg5 <- netVisual_heatmap(cellchat)
gg6 <- netVisual_heatmap(cellchat, measure = "weight")
gg5+gg6

#### 展示特定通路， 老三样
#cicle
df.net <- subsetCommunication(cellchat)
df.net$Tumor$pathway_name#查看通路名称有哪些
unique(df.net$Tumor$pathway_name)
pathways.show <- c("PTN") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, 
                      layout = "circle", edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

#heatmap
pathways.show <- c("PTN") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, 
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

#弦图
pathways.show <- c("PTN") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, 
                      layout = "chord", 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}
#基因表达情况
cellchat@meta$datasets = factor(cellchat@meta$datasets, 
                                levels = c("Normal", "Tumor")#设定组别的level
)
plotGeneExpression(cellchat, signaling = "GDF", split.by = "datasets", colors.ggplot = T)


#研究细胞类型之间的通讯差异

levels(object.list[[1]]@idents) 
levels(object.list[[2]]@idents) 
group.cellType <- c("B cells","Endothelial cells","Epithelial cells","Fibroblasts","Monocytes","T cells" )
group.cellType <- factor(group.cellType, levels = c("B cells","Endothelial cells","Epithelial cells","Fibroblasts","Monocytes","T cells" ))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})#合并细胞类型
cellchat <- mergeCellChat(object.list, add.names = names(object.list))#合并两组数据

#直接展示相减的结果
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)


#两组分别展示
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#展示两组间各类细胞incoming与outcoming通讯的强度
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) +
    colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  # 计算网络中心性
  centrality <- netAnalysis_computeCentrality(object.list[[i]])
  
  # 绘制信号角色散点图
  gg[[i]] <- netAnalysis_signalingRole_scatter(centrality,
                                               title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax)
}
## Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
## Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways

patchwork::wrap_plots(plots = gg)

{
  #展示两组间的差异
  #source('netAnalysis_signalingChanges_scatter.Rd')
  
  centrality <- netAnalysis_computeCentrality(object.list[[i]])
  gg1 <- netAnalysis_signalingChanges_scatter(centrality,
                                              idents.use = "Epithelial cells",
                                              signaling.exclude = "CXCL")
  
  gg1 <- netAnalysis_signalingChanges_scatter(centrality,
                                              idents.use = "B cells",
                                              signaling.exclude = "MIF")
  
  
  
  
  gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "cDC1",
                                              signaling.exclude = c("MIF"))
  patchwork::wrap_plots(plots = list(gg1,gg2))
  
}
gg1 +gg2
#### 继续探索outgoing与incoming的组间模式差异
suppressMessages(library(ComplexHeatmap))

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
centrality <- netAnalysis_computeCentrality(object.list[[i]])
ht1 = netAnalysis_signalingRole_heatmap(centrality, pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
centrality2 <- netAnalysis_computeCentrality(object.list[[i+1]])
ht2 = netAnalysis_signalingRole_heatmap(centrality2, pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


#incoming
ht1 = netAnalysis_signalingRole_heatmap(centrality, pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(centrality2, pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#overall
ht1 = netAnalysis_signalingRole_heatmap(centrality, pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(centrality2, pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")

draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))



##### 从配受体的角度

netVisual_bubble(cellchat, sources.use = c(1:6), 
                 targets.use =c(1:6),comparison =c(1,2), angle.x = 45)
####局部进行放大
interaction_name<-c("PTN_NCL","PTN_SDC1","PTN_SDC2","PTN_SDC4")
interaction_name<-as.data.frame(interaction_name)
netVisual_bubble(cellchat, sources.use = c(4), pairLR.use = interaction_name,
                 targets.use =c(1:6),comparison =c(1,2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Tumor", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Normal", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2


##### 以上的计算依赖的都是细胞通讯的可能性(强度)，接下来我们将通过配受体对基因表达的角度进一步研究
#以通路为单位
pos.dataset = "Tumor"
features.name = pos.dataset
#差异计算：
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, thresh.pc = 0.1, 
                                       thresh.fc = 0.1, thresh.p = 1)

#提取细胞通讯预测数据
net <- netMappingDEG(cellchat, features.name = features.name)

#提取Tumor组上调的配体以及normal上调的受体，
net.up <- subsetCommunication(cellchat, net = net, datasets = "Tumor",
                              ligand.logFC = 0.2, receptor.logFC = NULL)

#反之亦然
net.down <- subsetCommunication(cellchat, net = net, datasets = "Tumor",
                                ligand.logFC = -0.2, receptor.logFC = -0.1)

#提取其中的差异基因
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

save(net.up,net.down,gene.up,gene.down,file = "./SCI/data/cellchat_pathway_gene.Rdata")

#可视化：
#气泡图
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, 
                        sources.use = 5, targets.use = c(1:6), 
                        comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", 
                                            names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, 
                        sources.use = 4, targets.use = c(1,5,6), 
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in ", 
                                            names(object.list)[2]))
gg1 + gg2

#弦图
par(mfrow = c(1,2), xpd=TRUE)

netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(1:4), 
                     slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Down-regulated signaling in ", 
                                         names(object.list)[1]))
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(1:4), 
                     slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Up-regulated signaling in ", 
                                         names(object.list)[2]))

