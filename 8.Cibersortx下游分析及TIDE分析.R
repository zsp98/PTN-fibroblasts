## 在线下载COAD数据并整理成表达量矩阵
BiocManager::install("remotes")

BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")

BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")

BiocManager::install("SummarizedExperiment")

library(SummarizedExperiment)
library(TCGAbiolinks)

setwd("F:\\biomatics\\TCGA")

query <- GDCquery(project = "TCGA-COAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",   
                  workflow.type = "STAR - Counts", 
                  legacy = FALSE)

GDCdownload(query = query)
mydata=GDCprepare(query = query)
mydata2=as.data.frame(rowRanges(mydata))
geneexp <- assay(mydata,i = "unstrand")#tpm_unstrand fpkm_unstrand unstranded
geneexp=as.data.frame(geneexp)
geneexp1=cbind(gene_type=mydata2$gene_type,gene_name=mydata2$gene_name,geneexp)
COAD_Counts<-geneexp1
save(COAD_Counts,file = "./data/COAD_counts.rda")
##整理metadata
pkgs_in("jsonlite")

json_file <- "./COAD-bio/data/metadata.cart.2023-04-30.json"

metadata <- jsonlite::read_json(path = json_file, simplifyVector = T)
metadata <- tibble::tibble(
  file_name = metadata$file_name,
  md5sum = metadata$md5sum,
  TCGA_id_full = bind_rows(metadata$associated_entities)$entity_submitter_id,
  TCGA_id = stringr::str_sub(TCGA_id_full, 1, 16),
  patient_id = stringr::str_sub(TCGA_id, 1, 12),
  tissue_type_id = stringr::str_sub(TCGA_id, 14, 15),
  tissue_type = sapply(tissue_type_id, function(x) {
    switch(x,
           "01" = "Primary Solid Tumor",
           "02" = "Recurrent Solid Tumor",
           "03" = "Primary Blood Derived Cancer - Peripheral Blood",
           "05" = "Additional - New Primary",
           "06" = "Metastatic",
           "07" = "Additional Metastatic",
           "11" = "Solid Tissue Normal")}),   
  group = ifelse(tissue_type_id == "11", "Normal", "Tumor"))
save(metadata,file = "data/COAD_metadata.rda")

###整理临床信息
#2获得样本信息的表达矩阵
json <- jsonlite::fromJSON("./data/metadata.cart.2023-04-30.json")
entity_submitter_id <- sapply(json$associated_entities,function(x){x[,1]})
case_id <- sapply(json$associated_entities,function(x){x[,3]})
sample_case <- t(rbind(entity_submitter_id,case_id))

clinical <- read_tsv('./data/clinical/clinical.tsv')
clinical <- as.data.frame(clinical[duplicated(clinical$case_id),])
matrix <- merge(sample_case,clinical,by="case_id",all.x=T)
#提取关键的临床信息，如样本ID、年龄、性别、肿瘤、生存状态、病理分期等。
demo <- c("entity_submitter_id","age_at_index","ethnicity","gender","race",
          "vital_status","days_to_death","days_to_last_follow_up",
          "ajcc_pathologic_stage","ajcc_pathologic_t","ajcc_pathologic_m",
          "ajcc_pathologic_n","treatment_type")

matrix = matrix[,demo] #筛选需要的临床信息

colnames(matrix) <- c("TCGA_id_full","Age","Ethnicity","Gender","Race",
                      "Status","days_to_death","days_to_last_follow_up",
                      "Stage","T","M","N","Treatment") 
matrix = matrix[matrix$Status %in% c('Alive','Dead'),] #排除结局为"Not Reported"的Sample
#提取days_to_death、days_to_last_follow_up的信息用于死亡及存活患者生存时间的计算；同时，提取病理AJCC的分级Stage和TNM分期以及治疗信息，并重命名。
##读取metadata数据剔除非肿瘤样本
load("./data/COAD_metadata.rda")
###剔除非肿瘤样本
metadata<-subset(metadata,metadata$group=="Tumor")
table(metadata$group)
matrix<-subset(matrix,matrix$TCGA_id_full %in% metadata$TCGA_id_full)
#计算生存时间
####----03_simple_survival_analysis----####
matrix$days_to_last_follow_up[is.na(matrix$days_to_last_follow_up)] = 0 
matrix$days_to_death[is.na(matrix$days_to_death)] = 0   
matrix$days <- ifelse(matrix$Status=='Alive',matrix$days_to_last_follow_up,matrix$days_to_death)

matrix$days <- as.numeric(unlist(matrix$days))
matrix$month=round(matrix$days/30,0) #以month为单位，小数不保留
colnames(matrix)
rownames(matrix)<-matrix$TCGA_id_full

matrix<-select(matrix,c(1,2,4,6,9,10:12,14,15))
str(matrix)
matrix$Status<-ifelse(matrix$Status=="Alive",0,1)
matrix$Status <- factor(matrix$Status,levels = c(0,1),labels = c('Alive','Dead'))

matrix <-matrix[,-1]
save(matrix,file ="./data/2_matrix.Rdata" )
rm(list = ls());gc()

###TCGA-COAD数据TIDE分析
{
load("./data/COAD_counts.Rdata")
library(limma)
### 1 去除重复的基因名，且表达量取平均值，此处利用了limma包
exp <- COAD_counts

exp <- avereps(exp[,-1],      ##高频操作，     
               ID = exp$gene_name)#去除重复基因名，重复基因取平均值
exp <- as.data.frame(exp)
rownames(exp)<-exp$gene_name
exp[, 2:ncol(exp)] <- apply(exp[, 2:ncol(exp)], 2, as.numeric) # 将矩阵中的每列转换为数值型
exp<-exp[,-1]
#colnames(exp) = str_sub(colnames(exp),1,12)#缩短列名取前列名的前12位数字此针对TCGA数据
exp <- as.data.frame(exp)
exp<-log(exp+1)
###均值标准化处理
exp <- t(apply(exp, 1, function(x){x-(mean(x))})) #均值标准化

write.table(exp,file = "./codes/TIDE/TCGA/TCGA_COAD TIDE.txt", sep = "\t", row.names = TRUE,col.names = NA)

rm(list = ls());gc()
###下游分析可视化###

library(tidyverse)
library(survival)
library(survminer)
library(ggplot2)
library(readr)
load("./data/2_matrix.Rdata")
data <- read_csv("./data/3.12 CIBERSORTx_Fibrocluster_PTN_Adjusted_100.csv")
load("./data/COAD_metadata.rda")

#读取TIDE结果
TIDE<-read.csv(file = "./data/TCGA_TIDE_export_results.csv",header = T)
rownames(TIDE)<-TIDE$Patient
TIDE<-select(TIDE,c(4,6,11,12))

##合并CIbersortx数据和matrix（临床）数据
data<-as.data.frame(data)
data$Mixture<-gsub("\\.", "-", data$Mixture)

rownames(data)<-data$Mixture
data<-data[,-1]

data$dd<-rownames(data)
metadata$dd<-metadata$TCGA_id_full
TIDE$dd<-rownames(TIDE)

pdata<-merge(metadata,data,by="dd")
pdata<-merge(pdata,TIDE,by="dd")

rownames(pdata)<-pdata$dd
###
df<-pdata
table(df$group)
options(repr.plot.width=5, repr.plot.height=7)
ggplot(df, aes(x = df$group, y = df$`PTN- Fibroblasts`)) + 
  geom_boxplot(aes(color = df$group), outlier.shape = NA, width = 0.5) + 
  theme_classic() + guides(color = "none") + 
  scale_color_manual(values = c("#377EB8", "#E41A1C"),
                     labels = c("N (n = 41)","T (n = 483)")) + 
  geom_signif(comparisons = list(c("Normal", "Tumor")),
              map_signif_level = F,
              textsize = 6,
              test = "wilcox.test",
              step_increase = 0.1) + 
  geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.7) + 
  labs(x = "", y = "Percentage of infiltration") + 
  theme(
    legend.title = element_blank(),
    legend.position = c(0.2, 0.75),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) + 
  ggtitle("TCGA-COAD PTN- Fibroblasts")

####PTN+
ggplot(df, aes(x = df$group, y = df$`PTN+ Fibroblasts`)) + 
  geom_boxplot(aes(color = df$group), outlier.shape = NA, width = 0.5) + 
  theme_classic() + guides(color = "none") + 
  scale_color_manual(values = c("#377EB8", "#E41A1C"),
                     labels = c("N (n = 41)","T (n = 483)")) + 
  geom_signif(comparisons = list(c("Normal", "Tumor")),
              map_signif_level = F,
              textsize = 6,
              test = "wilcox.test",
              step_increase = 0.1) + 
  geom_point(position = position_jitter(width = 0.2), size =2, alpha = 0.7) + 
  labs(x = "", y = "Percentage of infiltration") + 
  theme(
    legend.title = element_blank(),
    legend.position = c(0.2, 0.75),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) + 
  ggtitle("TCGA-COAD PTN+ Fibroblasts")




# 提取肿瘤标本
names<-rownames(matrix)
data<-data[names,]
matrix$dd<-rownames(matrix)
TIDE<-TIDE[names,]

pdata<-merge(matrix,data,by="dd")
pdata<-merge(pdata,TIDE,by="dd")

rownames(pdata)<-pdata$dd
pdata$name<-substr(pdata$dd,1,12)
duplicate_rows <- duplicated(pdata$name)
# 提取不重复的行
pdata$Status<- ifelse(pdata$Status == 'Alive',0,1)
pdata <- pdata[!duplicate_rows, ]
df<-pdata
####另一种去重方式

#surv_cutpoint设置截取点
res.cut <- surv_cutpoint(df, time = "month", event = "Status",
                         variables = "PTN- Fibroblasts")
plot(res.cut, "PTN- Fibroblasts")

df$my_group <- ifelse(
  df$`PTN- Fibroblasts` > res.cut$cutpoint$cutpoint, "high", "low")

sfit <- survfit(Surv(month, Status)~my_group, data=df)
table(df$my_group)
p <- ggsurvplot(sfit, conf.int=F, pval=TRUE, 
                risk.table = T,
                pval.size = 6,
                palette = c('#FF2800', "#0000FF"),
                ggtheme = theme_classic(),
                legend.labs = c("high (n = 60)","low (n = 396)"))

options(repr.plot.width=5, repr.plot.height=5)
pg <- p$plot + labs(x = "Months", y = "Overall survival") + theme(
  axis.text = element_text(size = 15),
  axis.title = element_text(size = 22),
  legend.position = c(.2, .35),
  legend.text = element_text(size = 12),
  legend.title = element_blank(),
  legend.background = element_rect(fill = "transparent")
)
pg
table(df$my_group)
}

####GSE39582数据整理####
library(GEOquery)

GSE39582_rawdata <- getGEO('GSE39582',
                           
                           GSEMatrix= TRUE,
                           
                           getGPL = F)

###提取临床数据
GSE39582_pdata=pData(GSE39582_rawdata[[1]])
head(GSE39582_pdata) ##查看临床数据

##提取表达数据
GSE39582_exprSet=exprs(GSE39582_rawdata[[1]])
GSE39582_exprSet[1:5,1:5]

save(GSE39582_rawdata,GSE39582_exprSet,GSE39582_pdata,
     file = './data/Rawdata_GSE39582.RData')
#### Step II Annotation GSE ####
rm(list = ls())

library(tidyverse)

load("./data/anno_GPL570.RData")
load("./data/Rawdata_GSE39582.RData")

colnames(anno)
GSE39582_exprSet <- as.data.frame(GSE39582_exprSet)
GSE39582_exprSet$ID <- rownames(GSE39582_exprSet)
exprdf <- merge(anno,GSE39582_exprSet,by ='ID')
exprdf[1:6,1:6]
exprdf <- exprdf[,-1]

#删除没有gene symbol的探针组
exprdf<-exprdf[exprdf$Gene.Symbol!="",]
#有的探针组对应多个基因，用“///”分隔基因名，删掉这样的行
#当然也可以不删，通过分隔符'///'分割基因名，将第一个作为最后的基因名。
exprdf<-exprdf[!grepl("///", exprdf$Gene.Symbol),]

#有多个探针组对应同一个基因，取中值
#如果想取平均值，就把median改为mean
GSE39582_exprdf_uniq<-aggregate(.~Gene.Symbol,exprdf,median)

GSE39582_exprdf_uniq[1:6,1:6]

save(GSE39582_exprdf_uniq,file = './data/GSE39582_exprdf_uniq.RData')
rm(list = ls())
###整理临床数据##
rm(list = ls())

load("./data/Rawdata_GSE39582.RData")

colnames(GSE39582_pdata)

GSE39582_cli <- GSE39582_pdata[,c(2,12:18,21:24)]

colnames(GSE39582_cli) <- c('Sample_ID',
                            'Gender',
                            'Age',
                            'TNM',
                            'Tstage',
                            'Nstage',
                            'Mstage',
                            'tumor.location',
                            'RFS_status',
                            'RFS_time',
                            'OS_status',
                            'OS_time')

for (i in 2:ncol(GSE39582_cli)) {
  
  GSE39582_cli[,i] <- as.character(GSE39582_cli[,i])
  
  temp <- sapply(GSE39582_cli[,i], function(x){
    name <- unlist(strsplit(x,': '))
    x <- name[2]
  })
  
  temp[which(temp == 'N/A')] <- NA
  
  GSE39582_cli[,i] <- temp
  
  cat(i)
  
  rm(temp)
  
}
rm(i)
save(GSE39582_cli,file = './data/GSE39582_cli.RData')

####GSE39582数据Cibersortx及TIDE分析
library(tidyverse) ## 加载R包

options(stringsAsFactors = FALSE) #禁止chr转成factor
###1——整理临床数据
{
  load("./data/GSE39582_cli.RData")
  GSE39582_cli <- GSE39582_cli[GSE39582_cli$OS_time > 1, ]
  table(GSE39582_cli$Tstage)
  table(GSE39582_cli$Nstage)
  table(GSE39582_cli$Mstage)
  
  GSE39582_cli <- GSE39582_cli[GSE39582_cli$Tstage != 'Tis',]
  GSE39582_cli$Nstage <- ifelse(GSE39582_cli$Nstage == 'N0','N0','N+')
  GSE39582_cli <- GSE39582_cli[GSE39582_cli$Mstage != 'MX',]
  
  GSE39582_cli$Tstage <- ifelse(GSE39582_cli$Tstage == 'T0' |
                                  GSE39582_cli$Tstage == 'T1' | 
                                  GSE39582_cli$Tstage == "T2",'T0_T2','T3_T4')
  
  ## 查看空缺值（NA）
  source('./codes/Function/FindoutNA.R', encoding = "utf-8") ## 载入本课程的R function
  FindoutNA(GSE39582_cli)
  
  GSE39582_cli_omit <- na.omit(GSE39582_cli) ## 删除空缺值（注意：有空缺的数据整行删除）
  
  ## 比较数据集
  dim(GSE39582_cli)
  dim(GSE39582_cli_omit)
  
  nrow(GSE39582_cli) - nrow(GSE39582_cli_omit)
  
  rm(GSE39582_cli) ## 删除原始数据
  rm(FindoutNA) ## 删除Function
  
  ## 描述性分析 Table 1
  
  GSE39582_cli_omit$TNM %>% table()
  
  GSE39582_cli_omit$TNM <- factor(GSE39582_cli_omit$TNM,levels = c(0,1,2,3,4),labels = c('Stage0','StageI','StageII','StageIII','StageIV'))
  
  GSE39582_cli_omit <- GSE39582_cli_omit[,-1]
  
  str(GSE39582_cli_omit)
  
  GSE39582_cli_omit$Age <- as.numeric(GSE39582_cli_omit$Age)
  
  GSE39582_cli_omit$RFS_status <- factor(GSE39582_cli_omit$RFS_status,levels = c(0,1),labels = c('non-Relapse','Relapse'))
  
  GSE39582_cli_omit$RFS_time <- as.numeric(GSE39582_cli_omit$RFS_time)
  
  GSE39582_cli_omit$OS_status <- factor(GSE39582_cli_omit$OS_status,levels = c(0,1),labels = c('Alive','Death'))
  
  GSE39582_cli_omit$OS_time <- as.numeric(GSE39582_cli_omit$OS_time)
  
  save(GSE39582_cli_omit,file = ".//data/GSE39582_cli_omit.Rdata")
}
####2——整理表达数据
{
  load("./data/anno_GPL570.RData")
  load("./data/Rawdata_GSE39582.RData")
  
  colnames(anno)
  GSE39582_exprSet <- as.data.frame(GSE39582_exprSet)
  GSE39582_exprSet$ID <- rownames(GSE39582_exprSet)
  exprdf <- merge(anno,GSE39582_exprSet,by ='ID')
  exprdf[1:6,1:6]
  exprdf <- exprdf[,-1]
  
  #删除没有gene symbol的探针组
  exprdf<-exprdf[exprdf$Gene.Symbol!="",]
  #有的探针组对应多个基因，用“///”分隔基因名，删掉这样的行
  #当然也可以不删，通过分隔符'///'分割基因名，将第一个作为最后的基因名。
  exprdf<-exprdf[!grepl("///", exprdf$Gene.Symbol),]
  
  #有多个探针组对应同一个基因，取中值
  #如果想取平均值，就把median改为mean
  GSE39582_exprdf_uniq<-aggregate(.~Gene.Symbol,exprdf,median)
  
  GSE39582_exprdf_uniq[1:6,1:6]
  
  save(GSE39582_exprdf_uniq,file = './data/GSE39582_exprdf_uniq.RData')
  
}
rm(list = ls())

load("./data/GSE39582_exprdf_uniq.RData")

rownames(GSE39582_exprdf_uniq) <- GSE39582_exprdf_uniq$Gene.Symbol

GSE39582_exprdf_uniq <- GSE39582_exprdf_uniq[,-1]

##去处经过化疗治疗病人患者的数据
unique(GSE39582_pdata$characteristics_ch1.9)
GSE39582_pdata$characteristics_ch1.9<-as.character(GSE39582_pdata$characteristics_ch1.9)
pdata<-filter(GSE39582_pdata,GSE39582_pdata$characteristics_ch1.9 == "chemotherapy.adjuvant: N")
name<-rownames(pdata)
GSE39582<-t(GSE39582_exprdf_uniq)
GSE39582<-as.data.frame(GSE39582)
GSE39582<-filter(GSE39582,rownames(GSE39582) %in% name)
GSE39582<-t(GSE39582) 
GSE39582<-as.data.frame(GSE39582)
GSE39582 <- log2(GSE39582+1)#很多教程认为这一步不应该运行
###均值标准化处理，TIDE分析前必须
GSE39582 <- t(apply(GSE39582, 1, function(x){x-(mean(x))})) #均值标准化
write.table(GSE39582, "./codes/TIDE/GEO/GSE39582/input/GSE39582_matrix.txt",sep = "\t",
            row.names = TRUE, col.names = NA)

###Cibersortx及TIDE可视化
library(tidyverse)
library(survival)
library(survminer)
library(readr)
load("./data/GSE39582_cli.RData")
#data <- read.table("./codes/Cibersortx/GEO/GSE39582/data/CIBERSORTx_Job21_Adjusted.txt", header = TRUE, sep = "\t")
data <-read_csv("./output/CIBERSORTx_GSE39582_Fibrocluster_Adjusted.csv")

GSE39582_cli$Age <- as.numeric(GSE39582_cli$Age)

GSE39582_cli$RFS_status <- factor(GSE39582_cli$RFS_status,levels = c(0,1),labels = c('non-Relapse','Relapse'))

GSE39582_cli$RFS_time <- as.numeric(GSE39582_cli$RFS_time)

GSE39582_cli$OS_status <- factor(GSE39582_cli$OS_status,levels = c(0,1),labels = c('Alive','Death'))

GSE39582_cli$OS_time <- as.numeric(GSE39582_cli$OS_time)

load("./data/Rawdata_GSE39582.RData")
GSE39582_pdata<-select(GSE39582_pdata,c(2,11))
colnames(GSE39582_pdata)<-c("dd","group")
GSE39582_pdata$group <- ifelse(GSE39582_pdata$group %in% c("dataset: discovery", "dataset: validation"), "Tumor", "Normal")
table(GSE39582_pdata$group)
GSE39582_cli$dd<-rownames(GSE39582_cli)
GSE39582_cli<-merge(GSE39582_cli,GSE39582_pdata,by="dd")
##合并CIbersortx数据和matrix（临床）数据
data<-as.data.frame(data)
rownames(data)<-data$Mixture
data$dd<-rownames(data)
pdata<-merge(GSE39582_cli,data,by="dd")
rownames(pdata)<-pdata$dd
###
df<-pdata
table(df$group)
###PTN-
ggplot(df, aes(x = df$group, y = df$`PTN- Fibroblasts`)) + 
  geom_boxplot(aes(color = df$group), outlier.shape = NA, width = 0.5) + 
  theme_classic() + guides(color = "none") + 
  scale_color_manual(values = c("#377EB8", "#E41A1C"),
                     labels = c("N (n = 19)","T (n = 566)")) + 
  geom_signif(comparisons = list(c("Normal", "Tumor")),
              map_signif_level = F,
              textsize = 6,
              test = "wilcox.test",
              step_increase = 0.1) + 
  geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.7) + 
  labs(x = "", y = "Percentage of infiltration") + 
  theme(
    legend.title = element_blank(),
    legend.position = c(0.2, 0.75),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) + 
  ggtitle("GSE39582 PTN- Fibroblasts")

####PTN+
ggplot(df, aes(x = df$group, y = df$`PTN+ Fibroblasts`)) + 
  geom_boxplot(aes(color = df$group), outlier.shape = NA, width = 0.5) + 
  theme_classic() + guides(color = "none") + 
  scale_color_manual(values = c("#377EB8", "#E41A1C"),
                     labels = c("N (n = 19)","T (n = 566)")) + 
  geom_signif(comparisons = list(c("Normal", "Tumor")),
              map_signif_level = F,
              textsize = 6,
              test = "wilcox.test",
              step_increase = 0.1) + 
  geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.7) + 
  labs(x = "", y = "Percentage of infiltration") + 
  theme(
    legend.title = element_blank(),
    legend.position = c(0.2, 0.75),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) + 
  ggtitle("GSE39582 PTN+ Fibroblasts")

#剔除19个非肿瘤样本
pdata<-subset(pdata,pdata$group == "Tumor")
source('./codes/Function/FindoutNA.R', encoding = "utf-8") ## 载入本课程的R function
#OS_pdata<-pdata[,11:25]
OS_pdata<-pdata[,c(1,11:24)]
OS_pdata <- na.omit(OS_pdata)

OS_pdata$OS_status<- ifelse(OS_pdata$OS_status == 'Alive',0,1)

#RFS_pdata<-pdata[,c(9,10,13:25)]
RFS_pdata<-pdata[,c(9,10,13:20)]

RFS_pdata <- na.omit(RFS_pdata)
RFS_pdata$RFS_status<- ifelse(RFS_pdata$RFS_status == 'non-Relapse',0,1)

df<-OS_pdata
# 设置分界点另一种方法。

df$OS_status
df$OS_time
df$`PTN- Fibroblasts`
res.cut <- surv_cutpoint(df, time = "OS_time", event = "OS_status",
                         variables = "PTN- Fibroblasts")
plot(res.cut, "PTN- Fibroblasts")

df$my_group <- ifelse(
  df$`PTN- Fibroblasts` > res.cut$cutpoint$cutpoint, "high", "low")

sfit <- survfit(Surv(OS_time, OS_status)~my_group, data=df)
table(df$my_group)
p <- ggsurvplot(sfit, conf.int=F, pval=TRUE, 
                risk.table = T,
                pval.size = 6,
                palette = c('#FF2800', "#0000FF"),
                ggtheme = theme_classic(),
                legend.labs = c("high (n = 71)","low (n = 491)"))

options(repr.plot.width=5, repr.plot.height=5)
pg <- p$plot + labs(x = "Months", y = "Overall survival") + theme(
  axis.text = element_text(size = 15),
  axis.title = element_text(size = 22),
  legend.position = c(.2, .35),
  legend.text = element_text(size = 12),
  legend.title = element_blank(),
  legend.background = element_rect(fill = "transparent")
)
pg
table(df$my_group)  
