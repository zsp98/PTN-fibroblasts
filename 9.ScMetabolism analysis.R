install.packages(c("devtools", "data.table", "wesanderson", "Seurat", "devtools", "AUCell", "GSEABase", "GSVA", "ggplot2","rsvd"))
devtools::install_github("YosefLab/VISION@v2.1.0") #Please note that the version would be v2.1.
devtools::install_github("wu-yc/scMetabolism")
library(scMetabolism)
library(Seurat)
library(ggplot2)
library(rsvd)
load("F:/biomatics/SC RNA/SC RNA/GSE132465 human colorectal/scmetabolism/fibro_subcelltype.Rdata")
unique(sce_fibro$seurat_clusters)
sce_fibro$seurat_clusters <- ifelse(sce_fibro$seurat_clusters == "iCAF", "PTN+ Fibroblasts", "PTN- Fibroblasts")
unique(sce_fibro$seurat_clusters)
table(sce_fibro$seurat_clusters)
Idents(sce_fibro)<-sce_fibro$seurat_clusters
countexp.Seurat<-sc.metabolism.Seurat(obj = sce_fibro, 
                                      method = "AUCell", 
                                      imputation = F, ncores = 2, 
                                      metabolism.type = "KEGG")


#MedBioInfoCloud: rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[1:6]
#[1] "Glycolysis / Gluconeogenesis"            
#[2] "Citrate cycle (TCA cycle)"               
#[3] "Pentose phosphate pathway"               
#[4] "Pentose and glucuronate interconversions"
#[5] "Fructose and mannose metabolism"         
#[6] "Galactose metabolism"
DimPlot.metabolism(obj = countexp.Seurat, 
                   pathway = "Glycolysis / Gluconeogenesis", 
                   dimention.reduction.type = "umap", 
                   dimention.reduction.run = F, size = 1)

input.pathway<-rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[c(1,2,10,15)]
DotPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, 
                   phenotype = "seurat_clusters", norm = "y")
# #Box plot
input.pathway <- rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[1:3]
BoxPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "seurat_clusters", ncol = 5)

meta_score<-countexp.Seurat@assays[["METABOLISM"]]
meta_score<-as.data.frame(meta_score)
meta_score<-t(meta_score)
meta_score<-as.data.frame(meta_score)
cluster1_vals <- countexp.Seurat[["seurat_clusters"]] 
clustre_1<-subset(cluster1_vals,cluster1_vals$seurat_clusters == "PTN- Fibroblasts")
clustre_2<-subset(cluster1_vals,cluster1_vals$seurat_clusters == "PTN+ Fibroblasts")

rownames<-rownames(meta_score)
rownames<- substring(rownames, first = 7)
rownames <- gsub("\\.", "-", rownames)
rownames(meta_score)<-rownames

meta_score_1<-subset(meta_score,rownames(meta_score) %in% rownames(clustre_1))
meta_score_2<-subset(meta_score,rownames(meta_score) %in% rownames(clustre_2))

t_test_result <- t.test(meta_score_1$`Glycolysis / Gluconeogenesis`,meta_score_2$`Glycolysis / Gluconeogenesis` )

df<-cbind(cluster1_vals,meta_score)
colnames(df)[which(names(df) == "seurat_clusters")] <- "group"

# 使用ggplot绘制箱线图
df<-df[,c(1,2,3,11,16)]
plist <- list()

for (i in 2:5) {
  bar_tmp <- df[,c(colnames(df)[i], "group")]  # 从 'Exp_plot' 中提取当前基因的表达信息和样本组
  colnames(bar_tmp) <- c("Metabolism Score","group")
  theme <- theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5,size = 30),
          axis.text.x = element_text(hjust = 0.5,size = 15), 
          axis.text.y = element_text(hjust = 0.5,size = 15), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 16), 
          axis.line = element_line(size = 1), 
          legend.position = "none")
  p1 <- ggviolin(bar_tmp, x = 'group', y = 'Metabolism Score', fill = 'group',
                 palette = c("#4979b6","#d9352a"),
                 add = 'boxplot', add.params = list(fill = "white")) + theme +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 10,
                       bracket.size = 0.5, tip.length = 0.02,vjust = 0.4, method = 't.test') + theme(axis.text.y = element_text(size = 10)) + 
    ggtitle(colnames(df)[i]) + 
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) 
  plist[[i]] <- p1 # 将生成的图形存储在 'plist' 中
}
plot_grid(plist[[2]], plist[[3]],plist[[4]],plist[[5]],ncol = 2)
##MCT4/MCT1乳酸代谢流
genes<-c("SLC16A1","SLC16A3")
Idents(sce_fibro)<-sce_fibro$seurat_clusters
p2 <- DotPlot(sce_fibro, features = genes) +
  coord_flip() +
  xlab("Genes") +
  ylab("Group") +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    text = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_blank()# 将y轴上的字体大小设置为14
  ) 
DT<-p2$data
class(DT)
DT<-t(DT)
p2$data
ratio_1<-1.1868512/0.3183077
ratio_2<-0.1414773/0.2760017
ratio_3<-ratio_1/ratio_2
r_df<-data.frame(ratio_1,ratio_2)
colnames(r_df)<-c("PTN- Fibroblats","PTN+ Fibroblasts")
rownames(r_df)<-"SLC16A3/SLC16A1 Ratio"
r_df<-t(r_df)
r_df<-as.data.frame(r_df)
r_df$group<-rownames(r_df)

pb1 <- ggplot(r_df, aes(x = group, y = `SLC16A3/SLC16A1 Ratio`, color = `SLC16A3/SLC16A1 Ratio`, label = round(`SLC16A3/SLC16A1 Ratio`, 3))) +
  geom_point(size = 5) +
  scale_color_gradient(low = "blue", high = "red") +
  theme(axis.text.x = element_text(size = 15),  # 设置x轴标签字体大小
        axis.text.y = element_text(size = 15),  # 设置y轴标签字体大小
        axis.title.x = element_blank()) +  #隐藏x轴标题
  geom_text(position = position_dodge(0.9), vjust = -1, hjust = 0.5) +  # 将数值标记在点上
  theme(axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = 15))

pb1