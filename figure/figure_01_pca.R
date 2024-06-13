library(Biobase)
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(svglite)
library(showtext) # 解决字体问题
library('DESeq2')
library("factoextra") # get_pca_var
library('ggrepel')
library('ggplot2')
library('RColorBrewer')
library('limma')
library('sva')
library('reshape2') # melt


sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")

source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# constants
useVersion   <- "version1"
colorScheme  <- getColorScheme()

# set output folder
outputFolder <- file.path("figures", "output", "figure_01")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

theme(text=element_text("serif"))
# 蛋白组
###########################################################################################
# load proteins (highly variable)
proteins <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "all","CLCs","all", "proteins.RDS"))
proteins <- proteins %>% .[fData(.)$hasHighVariance, ] %>% removeNAsFromESet()
# pca
pca <- prcomp(x = t(exprs(proteins)),
              scale. = TRUE)

variances <- pca %>%
  summary() %>%
  .$importance %>%
  .[1, ] %>%
  .^2

variances_explained <- round(100 * variances / sum(variances), 2)

df <- pca$x %>%
  as.data.frame() %>%
  rownames_to_column("proteomics_id") %>%
  left_join(dplyr::select(pData(proteins), proteomics_id = ID,everything()),by = "proteomics_id") %>%
  dplyr::select(PC1, PC2, name, group)

# plot
p1 <- df %>%
  ggplot(aes(x = PC1, y = PC2, color = group)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#AAAAAA") +
  geom_hline(yintercept = 0, linetype = 2, color = "#AAAAAA") +
  geom_point(size = 2.5, stroke = 0, alpha = 0.8) +
  scale_color_manual(values = colorScheme$type) +
  xlab(paste0("PC1 (", variances_explained[["PC1"]], "%)")) +
  ylab(paste0("PC2 (", variances_explained[["PC2"]], "%)")) +
  theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1))
p1

# export
ggsave(
  "pca_highly_variable_proteins.pdf",
  plot = p1,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 6,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
ggsave(
  "pca_highly_variable_proteins.jpg",
  plot = p1,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 6,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)

############################################################
# 这里不用高可变蛋白，用的是所有蛋白的表达
proteins <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "all","CLCs","all", "proteins.RDS"))

# pca
pca <- prcomp(x = t(exprs(proteins)),
              scale. = TRUE)

variances <- pca %>%
  summary() %>%
  .$importance %>%
  .[1, ] %>%
  .^2

variances_explained <- round(100 * variances / sum(variances), 2)

df <- pca$x %>%
  as.data.frame() %>%
  rownames_to_column("proteomics_id") %>%
  left_join(dplyr::select(pData(proteins), proteomics_id = ID,everything()),by = "proteomics_id") %>%
  dplyr::select(PC1, PC2, name, group)

# plot
p2 <- df %>%
  ggplot(aes(x = PC1, y = PC2, color = group)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#AAAAAA") +
  geom_hline(yintercept = 0, linetype = 2, color = "#AAAAAA") +
  geom_point(size = 4, stroke = 0, alpha = 0.8) +
  scale_color_manual(values = colorScheme$type) +
  xlab(paste0("PC1 (", variances_explained[["PC1"]], "%)")) +
  ylab(paste0("PC2 (", variances_explained[["PC2"]], "%)"))+
  theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1))
p2

ggsave(
  "pca_all_gene_proteins.pdf",
  plot = p2,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 6,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
ggsave(
  "pca_all_gene_proteins.jpg",
  plot = p2,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 6,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)

############################################################
# 这里用rna和蛋白质共有基因看
proteomics_common <- readRDS(file = file.path("/data/home/mali/ke_data/Protein/RDS/proteomics_common.RDS"))
#colnames(proteomics_common) <- pData(proteomics_common)$ID
#saveRDS(proteomics_common,file.path("/data/home/mali/ke_data/Protein/RDS/proteomics_common.RDS"))
proteins <- proteomics_common
# pca
pca <- prcomp(x = t(exprs(proteins)),
              scale. = TRUE)

variances <- pca %>%
  summary() %>%
  .$importance %>%
  .[1, ] %>%
  .^2

variances_explained <- round(100 * variances / sum(variances), 2)

df <- pca$x %>%
  as.data.frame() %>%
  rownames_to_column("proteomics_id") %>%
  left_join(dplyr::select(pData(proteins), proteomics_id = ID,everything()),by = "proteomics_id") %>%
  dplyr::select(PC1, PC2, name, group)

# plot
p3 <- df %>%
  ggplot(aes(x = PC1, y = PC2, color = group),label = TRUE) +
  geom_vline(xintercept = 0, linetype = 2, color = "#AAAAAA") +
  geom_hline(yintercept = 0, linetype = 2, color = "#AAAAAA") +
  geom_point(size = 4, stroke = 0, alpha = 0.8) +
  scale_color_manual(values = colorScheme$type) +
  xlab(paste0("PC1 (", variances_explained[["PC1"]], "%)")) +
  ylab(paste0("PC2 (", variances_explained[["PC2"]], "%)"))+
  theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1))
p3

ggsave(
  "pca_common_gene_proteins.pdf",
  plot = p3,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 6,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
ggsave(
  "pca_common_gene_proteins.jpg",
  plot = p3,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 6,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
######################################### 修改名字
type <- "RNA_protein"
name <- "common_gene"
prefix <- "RNA_protein_common"
#########################################
#> length( rownames(exprs(proteins)) %>% unlist %>% unique)
#[1] 3136

# 用keyType= 'SYMBOL'会没有结果，一定要转换
markers <- bitr(rownames(exprs(proteins)) %>% unlist %>% unique, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
goBP <- enrichGO(gene = markers$ENTREZID,
                     OrgDb = 'org.Mm.eg.db',
                     keyType= 'ENTREZID',
                     ont = 'BP',
                     universe = keys(org.Mm.eg.db),
                     pAdjustMethod = 'BH',
                     pvalueCutoff=0.01,
                     qvalueCutoff=0.01)
goBPSim <- clusterProfiler::simplify(goBP,
                                     cutoff = 0.5,
                                     by = 'p.adjust',
                                     select_fun = min)
## check and plot
outputFolder2 <- file.path(outputFolder, "clusterbc",type,name)
if (!is.null(outputFolder2) && !dir.exists(outputFolder2)) dir.create(outputFolder2, recursive = TRUE)
savepath <- outputFolder2
write.csv(as.data.frame(goBPSim), paste0(prefix, '_cluster_cp_BP_sim.csv') %>% file.path(savepath, .))
write.csv(as.data.frame(goBP), paste0(prefix, '_cluster_cp_BP.csv') %>% file.path(savepath, .))

#-------------------
# 绘制GO dotplot
dotplot(goBPSim, showCategory = 10, font.size =6)
ggsave(file.path(outputFolder,paste0(type,"_",name,"_",prefix,'_cp_BP_dotplot_Sim_10.jpg')), width = 5,height = 4)
ggsave(file.path(outputFolder,paste0(type,"_",name,"_",prefix,'_cp_BP_dotplot_Sim_10.pdf')), width = 5,height = 4)
dotplot(goBP, showCategory = 10, font.size =6)
ggsave(file.path(outputFolder,paste0(type,"_",name,"_",prefix,'_cp_BP_dotplot_10.jpg')), width = 5,height = 4)
ggsave(file.path(outputFolder,paste0(type,"_",name,"_",prefix,'_cp_BP_dotplot_10.pdf')), width = 5,height = 4)

###############################################################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('/data/home/mali/ke_data/transcription/RDS/degres_condi_Mock.RData')

rld <- rlog(degres)
pca <- prcomp(x = t(rldData),
              scale. = TRUE)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rld)[, 1], ID = rownames(colData(rld)))
colnames(pcaData) <- c("PC1","PC2","group","name")
# plot
p4 <- pcaData %>%
  ggplot(aes(x = PC1, y = PC2, color = group)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#AAAAAA") +
  geom_hline(yintercept = 0, linetype = 2, color = "#AAAAAA") +
  geom_point(size = 4, stroke = 0, alpha = 0.8) +
  scale_color_manual(values = colorScheme$type) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1))
p4

ggsave(
  "pca_all_gene_RNA.pdf",
  plot = p4,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 6,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
ggsave(
  "pca_all_gene_RNA.jpg",
  plot = p4,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 6,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)

#########---------------------------------------- PCA enrich

RNA_eSet_Scale <- readRDS("/data/home/mali/ke_data/transcription/RDS/RNA_rld_Scale.RDS")
colnames(RNA_eSet_Scale) <- RNA_eSet_Scale$name
clcString <- "median_expression_"
RNA_eSet_Scale <- addMedianExpressions(eSet = RNA_eSet_Scale,
                                 use = "group",
                                 removeReplicates = FALSE,
                                 recognitionString = clcString)


# 提取变量的分析结果，加载包出现在函数第一次出现前，知道函数与包的关系
var <- get_pca_var(pca)
head(var$contrib)
# 取pca1和pca2的，从大到小排序，看什么基因影响
pca1 <- sort(var$contrib[,1],decreasing = T)
pca2 <- sort(var$contrib[,2],decreasing = T)
head(pca1)
head(pca2)
length(pca1)
length(pca2)
pca1_data <- data.frame(gene = names(pca1), value = pca1)
pca2_data <- data.frame(gene = names(pca2), value = pca2)
write.table(pca1_data, file = file.path(outputFolder,"RNA_pca1.txt"), row.names =F ,col.name = T,quote = F,sep = "\t")
write.table(pca2_data, file = file.path(outputFolder,"RNA_pca2.txt"), row.names =F ,col.name = T,quote = F,sep = "\t")


##---------------------------------------
# 取前100个贡献基因
n = 2000
pca1_name <- names(pca1[1:n])
pca2_name <- names(pca2[1:n])
# 这些基因的表达量计算，样本取平均
gene_exp <-dplyr::select(fData(RNA_eSet_Scale),starts_with("median_expression_"))
colnames(gene_exp) <- c("YV","OV","O","Y")
gene_mat_pca1 <- gene_exp[rownames(gene_exp) %in% pca1_name,]
gene_mat_pca2 <- gene_exp[rownames(gene_exp) %in% pca2_name,]

gene_mat_pca1$gene <- rownames(gene_mat_pca1)
gene_mat_pca2$gene <- rownames(gene_mat_pca2)
# 数据变长型
gene_mat_pca1_mean <- melt(gene_mat_pca1)
colnames(gene_mat_pca1_mean) <- c("gene","sample","value")
head(gene_mat_pca1_mean)
#> head(gene_mat_pca1_mean)
# gene sample     value
# 1  19241     YV -1.386244
# 2 117158     YV -1.367452
# 3  11668     YV -1.370400
# 4  67971     YV -1.414893
# 5 109672     YV -1.376841
# 6  14862     YV -1.304813
gene_mat_pca2_mean <- melt(gene_mat_pca2)
colnames(gene_mat_pca2_mean) <- c("gene","sample","value")
# 每个分组的平均表达
# pc1
gene_mat_pca1_mean_sample <- NULL
gene_mat_pca1_mean_sample_tmp <- gene_mat_pca1[,1:4]
colnames(gene_mat_pca1[,1:4])
#[1] "median_expression_YV" "median_expression_OV" "median_expression_O"
#[4] "median_expression_Y"
gene_mat_pca1_mean_sample$sample <- c("YV","OV","O","Y")
gene_mat_pca1_mean_sample$value <- apply(gene_mat_pca1_mean_sample_tmp,2,mean)
gene_mat_pca1_mean_sample <- data.frame(gene_mat_pca1_mean_sample)
#gene_mat_pca1_mean_sample$sample <- factor( gene_mat_pca1_mean_sample$sample, levels = c("Y","O","OV","YV"))
# pc2
gene_mat_pca2_mean_sample <- NULL
gene_mat_pca2_mean_sample_tmp <- gene_mat_pca2[,1:4]
colnames(gene_mat_pca2[,1:4])
#[1] "median_expression_YV" "median_expression_OV" "median_expression_O"
#[4] "median_expression_Y"
gene_mat_pca2_mean_sample$sample <- c("YV","OV","O","Y")
gene_mat_pca2_mean_sample$value <- apply(gene_mat_pca2_mean_sample_tmp,2,mean)
gene_mat_pca2_mean_sample <- data.frame(gene_mat_pca2_mean_sample)
#gene_mat_pca2_mean_sample$sample <- factor( gene_mat_pca2_mean_sample$sample, levels = c("YV","Y","OV","O"))

pdf(file.path(outputFolder,"pca_line_RNA.pdf"),width=3,height=3)
p4_pc1 <- ggplot()+
  geom_line(data = gene_mat_pca1_mean,aes(x=sample,y=value,group = gene),color ="grey60",size=0.5)+
  geom_line(data = gene_mat_pca1_mean_sample, aes(x=sample,y=value,group = 1),color = "black",size=1)+
  xlab("PC1")+#横坐标名称
  ylab("normalised mean expression")+#纵坐标名称
  theme_bw() +#去掉背景灰色
   #ylim(0,15000) +
   #NoLegend() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA),strip.background = element_rect(fill = NA,colour = NA), axis.ticks.x = element_blank())
p4_pc2 <- ggplot()+
  geom_line(data = gene_mat_pca2_mean,aes(x=sample,y=value,group = gene),color ="grey60",size=0.5)+
  geom_line(data = gene_mat_pca2_mean_sample, aes(x=sample,y=value,group = 1),color = "black",size=1)+
  xlab("PC2")+#横坐标名称
  ylab("normalised mean expression")+#纵坐标名称
  theme_bw() +#去掉背景灰色
   #ylim(0,15000) +
   #NoLegend() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA),strip.background = element_rect(fill = NA,colour = NA), axis.ticks.x = element_blank())

p4_pc1
p4_pc2
dev.off()

#########----------------------------------------pca cluster
# 取前2000个贡献基因
n = 2000
pca1_name <- names(pca1[1:n])
pca2_name <- names(pca2[1:n])
# 这些基因的表达量计算，样本取平均
gene_exp <-dplyr::select(fData(RNA_eSet_Scale),starts_with("median_expression_"))
colnames(gene_exp) <- c("YV","OV","O","Y")
gene_mat_pca1 <- gene_exp[rownames(gene_exp) %in% pca1_name,]
gene_mat_pca2 <- gene_exp[rownames(gene_exp) %in% pca2_name,]

########################################### 修改
# 转录组
type <- "RNA"
scaleCount <- gene_mat_pca2
name <- "PC2"

#pca1_data <- read.table(file = file.path(outputFolder,"RNA_pca1.txt"),header = T)
pca1_data <- read.table(file = file.path(outputFolder,"RNA_pca2.txt"),header = T)

###########################################
##-------do cluster
clNum <- 2
kClust10 <- kmeans(scaleCount, centers = clNum, algorithm= 'MacQueen', nstart = 1000, iter.max = 10)

kClust10$ID <- names(kClust10$cluster)
kClust10_used <- cbind.data.frame(ID = kClust10$ID,cl = kClust10$cluster)

scaleCount <- cbind(scaleCount,kClust10_used)
##-------
write_csv(scaleCount, paste0("/data/home/mali/ke_data/transcription/count/type","_",name,"_",n,"_kmeans2.csv"))
##-------

#-------pca cluster plot
kmeansRes <- read_csv(paste0("/data/home/mali/ke_data/transcription/count/type","_",name,"_",n,"_kmeans2.csv"))
##-------
# pc_1
gene_mat_pca1_1 <- kmeansRes[which(kmeansRes$cl == "1"),1:5]
gene_mat_pca1_1 <- gene_mat_pca1_1 %>% column_to_rownames("ID")
gene_mat_pca1_1$gene <- rownames(gene_mat_pca1_1)
# 数据变长型
gene_mat_pca1_1_mean <- melt(gene_mat_pca1_1)
head(gene_mat_pca1_1_mean)
table(gene_mat_pca1_1_mean$variable)
colnames(gene_mat_pca1_1_mean) <- c("gene","sample","value")
# 每个分组的平均表达
gene_mat_pca1_1_mean_sample <- NULL
gene_mat_pca1_1_mean_sample_tmp <- gene_mat_pca1_1[,1:4]
colnames(gene_mat_pca1_1[,1:4])
#[1] "median_expression_YV" "median_expression_OV" "median_expression_O"
#[4] "median_expression_Y"
gene_mat_pca1_1_mean_sample$sample <- c("YV","OV","O","Y")
gene_mat_pca1_1_mean_sample$value <- apply(gene_mat_pca1_1_mean_sample_tmp,2,mean)
gene_mat_pca1_1_mean_sample <- data.frame(gene_mat_pca1_1_mean_sample)

# pc1_2
gene_mat_pca1_2 <- kmeansRes[which(kmeansRes$cl == "2"),1:5]
gene_mat_pca1_2 <- gene_mat_pca1_2 %>% column_to_rownames("ID")
gene_mat_pca1_2$gene <- rownames(gene_mat_pca1_2)
# 数据变长型
gene_mat_pca1_2_mean <- melt(gene_mat_pca1_2)
head(gene_mat_pca1_2_mean)
table(gene_mat_pca1_2_mean$variable)
colnames(gene_mat_pca1_2_mean) <- c("gene","sample","value")
# 每个分组的平均表达
gene_mat_pca1_2_mean_sample <- NULL
gene_mat_pca1_2_mean_sample_tmp <- gene_mat_pca1_2[,1:4]
colnames(gene_mat_pca1_2[,1:4])
#[1] "median_expression_YV" "median_expression_OV" "median_expression_O"
#[4] "median_expression_Y"
gene_mat_pca1_2_mean_sample$sample <- c("YV","OV","O","Y")
gene_mat_pca1_2_mean_sample$value <- apply(gene_mat_pca1_2_mean_sample_tmp,2,mean)
gene_mat_pca1_2_mean_sample <- data.frame(gene_mat_pca1_2_mean_sample)

pdf(file.path(outputFolder,paste0(name,"_",n,"_line_cluster_RNA.pdf")),width=3,height=3)
p4_1 <- ggplot()+
  geom_line(data = gene_mat_pca1_1_mean,aes(x=sample,y=value,group = gene),color ="grey60",size=0.5)+
  geom_line(data = gene_mat_pca1_1_mean_sample, aes(x=sample,y=value,group = 1),color = "black",size=1)+
  xlab("PC1")+#横坐标名称
  ylab("normalised mean expression")+#纵坐标名称
  theme_bw() +#去掉背景灰色
   #ylim(0,15000) +
   #NoLegend() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA),strip.background = element_rect(fill = NA,colour = NA), axis.ticks.x = element_blank())
p4_2 <- ggplot()+
  geom_line(data = gene_mat_pca1_2_mean,aes(x=sample,y=value,group = gene),color ="grey60",size=0.5)+
  geom_line(data = gene_mat_pca1_2_mean_sample, aes(x=sample,y=value,group = 1),color = "black",size=1)+
  xlab("PC2")+#横坐标名称
  ylab("normalised mean expression")+#纵坐标名称
  theme_bw() +#去掉背景灰色
   #ylim(0,15000) +
   #NoLegend() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA),strip.background = element_rect(fill = NA,colour = NA), axis.ticks.x = element_blank())

p4_1
p4_2
dev.off()



### 基因及富集分析
library(msigdbr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(enrichplot)
library(GSEABase)
library(tidyverse)

LOGFC_col <- n
# 导出影响pca的前n个基因
pca1_1_name_50 <- NULL
pca1_1_name_50$ENTREZID <- gene_mat_pca1_1$gene
pca1_1_name_50$exp <- pca1_data[pca1_data$gene %in% gene_mat_pca1_1$gene,"value"]
pca1_1_name_50 <- as.data.frame(pca1_1_name_50)
pca1_1_name_50 <- pca1_1_name_50[order(pca1_1_name_50$exp,decreasing = T),]

pca1_2_name_50 <- NULL
pca1_2_name_50$ENTREZID <- gene_mat_pca1_2$gene
pca1_2_name_50$exp <- pca1_data[pca1_data$gene %in% gene_mat_pca1_2$gene,"value"]
pca1_2_name_50 <- as.data.frame(pca1_2_name_50)
pca1_2_name_50 <- pca1_2_name_50[order(pca1_2_name_50$exp,decreasing = T),]

write.table(pca1_1_name_50, file = file.path(outputFolder,paste0(type,"_",name,"_",LOGFC_col,"_1.txt")), row.names = F,col.name = F,quote = F,sep = "\t")
write.table(pca1_2_name_50, file = file.path(outputFolder,paste0(type,"_",name,"_",LOGFC_col,"_2.txt")), row.names = F,col.name = F,quote = F,sep = "\t")


# 自定义基因集读取
geneset_immune <- read.gmt("/data/home/mali/biosoft/MSigDB/geneset_immune.gmt")

#开始ID转换， 制作能输入到GSEA分析函数中的 genelist（pca1，pca2原本是SYMBOL作为名字，制作 ENTREZID作为名字）
pca1_1_name_change <- bitr(pca1_1_name_50$ENTREZID, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Mm.eg.db")
pca1_1_name_50_tmp <- left_join(pca1_1_name_50,pca1_1_name_change, by="ENTREZID")
pca1_1_name_50_exp <- pca1_1_name_50_tmp$exp
#pca1_1_name_50_exp <- data.frame(pca1_1_name_50_exp)
#rownames(pca1_1_name_50_exp) <- toupper(pca1_1_name_50_tmp$SYMBOL)
#colnames(pca1_1_name_50_exp) <- "exp"
names(pca1_1_name_50_exp) <- toupper(pca1_1_name_50_tmp$SYMBOL)

pca1_2_name_change <- bitr(pca1_2_name_50$ENTREZID, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Mm.eg.db")
pca1_2_name_50_tmp <- left_join(pca1_2_name_50,pca1_2_name_change, by="ENTREZID")
pca1_2_name_50_exp <- pca1_2_name_50_tmp$exp
#pca1_2_name_50_exp <- data.frame(pca1_2_name_50_exp)
#rownames(pca1_2_name_50_exp) <- toupper(pca1_2_name_50_tmp$SYMBOL)
#colnames(pca1_2_name_50_exp) <- "exp"
names(pca1_2_name_50_exp) <- toupper(pca1_2_name_50_tmp$SYMBOL)
head(pca1_2_name_50_exp)


write.table(pca1_1_name_50_exp, file = file.path(outputFolder,paste0(type,"_",name,"_",LOGFC_col,"_1_symbol.txt")), row.names = T,col.name = F,quote = F,sep = "\t")
write.table(pca1_2_name_50_exp, file = file.path(outputFolder,paste0(type,"_",name,"_",LOGFC_col,"_2_symbol.txt")), row.names = T,col.name = F,quote = F,sep = "\t")


  pca1_1_gsea_h <- GSEA(pca1_1_name_50_exp, TERM2GENE=geneset_immune, verbose=FALSE,minGSSize = 2, maxGSSize = 200,pvalueCutoff = 0.5)
  pca1_1_gsea_results <- pca1_1_gsea_h@result
  pca1_1_gsea_results <- pca1_1_gsea_results[order(pca1_1_gsea_results$qvalue, decreasing = F),]
  write.csv(pca1_1_gsea_results,file = file.path(outputFolder,paste0(type,"_",name,"_",LOGFC_col,"_1_gsea.csv")))
  pca1_1_plot<- pca1_1_gsea_results[1:10,]
  pca1_1_plot$ID <- factor(pca1_1_plot$ID,levels =pca1_1_plot$ID )

  pca1_2_gsea_h <- GSEA(pca1_2_name_50_exp, TERM2GENE=geneset_immune, verbose=FALSE,minGSSize = 2, maxGSSize = 200,pvalueCutoff = 0.5)
  pca1_2_gsea_results <- pca1_2_gsea_h@result
  pca1_2_gsea_results<-pca1_2_gsea_results[order(pca1_2_gsea_results$qvalue, decreasing = F),]
  write.csv(pca1_2_gsea_results,file = file.path(outputFolder,paste0(type,"_",name,"_",LOGFC_col,"_2_gsea.csv")))
  pca1_2_plot<- pca1_2_gsea_results[1:10,]
  pca1_2_plot$ID <- factor(pca1_2_plot$ID,levels =pca1_2_plot$ID )


pdf(file.path(outputFolder,paste0(type,"_",name,"_",LOGFC_col,"_1_gsea.pdf")),width=5,height=4)
  p4_pca1_1 <- ggplot(na.omit(pca1_1_plot), aes(x = ID)) +
   geom_col(aes(y=-log10(p.adjust)), colour = "black", width = 0.6)+
   geom_line(aes(y=NES/1, group =1), colour = "grey60") +
   scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Normalized Enrichment Score")) +
   scale_colour_manual(values = c("black", "grey60")) +
   labs(y = "-log10(p.adjust)",x= "pca1_1",colour = "black")+
   #coord_flip()+
   #NoLegend() +
   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA),strip.background = element_rect(fill = NA,colour = NA),axis.text.x = element_text(angle = 90, hjust = 1),axis.line = element_line(colour = "black"))
p4_pca1_1
dev.off()


pdf(file.path(outputFolder,paste0(type,"_",name,"_",LOGFC_col,"_2_gsea.pdf")),width=5,height=4)
  p4_pca1_2 <- ggplot(pca1_2_plot, aes(x = ID)) +
   geom_col(aes(y=-log10(p.adjust)), colour = "black", width = 0.6)+
   geom_line(aes(y=NES/4, group =1), colour = "grey60") +
   scale_y_continuous(sec.axis = sec_axis(~.*4, name = "Normalized Enrichment Score")) +
   scale_colour_manual(values = c("black", "grey60")) +
   labs(y = "-log10(p.adjust)",x= "pca1_2",colour = "black")+
   #coord_flip()+
   #NoLegend() +
   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA),strip.background = element_rect(fill = NA,colour = NA),axis.text.x = element_text(angle = 90, hjust = 1),axis.line = element_line(colour = "black"))
p4_pca1_2
dev.off()

##################============================================================== 正常的GO富集分析
prefix <- "kmeans2"
for (i in kmeansRes$cl %>% unique) {
  ## BP
  goBP <- enrichGO(gene = kmeansRes %>% filter(cl == i) %>% .$ID %>% unlist %>% unique,
                   OrgDb = 'org.Mm.eg.db',
                   keyType= 'ENTREZID',
                   ont = 'BP',
                   universe = keys(org.Mm.eg.db),
                   pAdjustMethod = 'BH',
                   pvalueCutoff=0.01,
                   qvalueCutoff=0.01)


  goBPSim <- clusterProfiler::simplify(goBP,
                                       cutoff = 0.5,
                                       by = 'p.adjust',
                                       select_fun = min)
  ## check and plot
  outputFolder2 <- file.path(outputFolder, "clusterbc",type,name)
  if (!is.null(outputFolder2) && !dir.exists(outputFolder2)) dir.create(outputFolder2, recursive = TRUE)

  savepath <- outputFolder2
  write.csv(as.data.frame(goBPSim),
            paste0(prefix, '_cluster', i, '_cp_BP.csv') %>% file.path(savepath, .))

}

kall <- lapply(kmeansRes$cl %>% unique, function(x) {

  eachG <- kmeansRes %>% filter(cl == x) %>% .$ID %>% unlist %>% unique

  return(eachG)

}) %>%
  set_names(kmeansRes$cl %>% unique %>% paste0('cluster', .))

save(kall, file = file.path(savepath,paste0(type,"_",name,"_",prefix,'_kall.RData')))

#-------------------
# 绘制compareCluster的GO dotplot
kallGOBP <- compareCluster(geneCluster = kall,
                           fun = 'enrichGO',
                           OrgDb = 'org.Mm.eg.db',
                           keyType= 'ENTREZID',
                           ont = 'BP',
                           universe = keys(org.Mm.eg.db),
                           pAdjustMethod = 'BH',
                           pvalueCutoff=0.01,
                           qvalueCutoff=0.1)

# 通过去除富集GO项的冗余来简化富集GO和gseGO的输出
kallGOBPSim <- clusterProfiler::simplify(kallGOBP,
                                         cutoff = 0.9,
                                         by = 'p.adjust',
                                         select_fun = min)

dotplot(kallGOBPSim, showCategory = 10)
ggsave(file.path(outputFolder,paste0(type,"_",name,"_",prefix,'_cp_BP_dotplot_Sim_10.jpg')), width = 6,height = 5)
ggsave(file.path(outputFolder,paste0(type,"_",name,"_",prefix,'_cp_BP_dotplot_Sim_10.pdf')), width = 6,height = 5)

dotplot(kallGOBP, showCategory = 10)
ggsave(file.path(outputFolder,paste0(type,"_",name,"_",prefix,'_cp_BP_dotplot_10.jpg')), width = 6,height = 5)
ggsave(file.path(outputFolder,paste0(type,"_",name,"_",prefix,'_cp_BP_dotplot_10.pdf')), width = 6,height = 5)

# 保存到csv
kallGOBP %>%
  as.data.frame %>%
  write_csv(file.path(savepath,paste0(type,"_",name,"_",prefix,'_cp_BP.csv')))
kallGOBPSim %>%
  as.data.frame %>%
  write_csv(file.path(savepath,paste0(type,"_",name,"_",prefix,'_cp_BP_sim.csv')))


###########============================================== FIG1 All

library(ggplot2)
library(grid)
pdf(file = file.path(outputFolder,"Fig1_ALL_rna_2.pdf"), height =2.5, width = 8)
# 创建一个3x3的网格布局
gl <- grid.layout(nrow=2, ncol=3)
# 将小图放置在网格中
pushViewport(viewport(layout=gl))
vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col=y)
print(p4, vp=vplayout(1:2,1))
#print(p4_1, vp=vplayout(1,2))
print(p4_2, vp=vplayout(1:2,2))
#print(p4_pca1_1, vp=vplayout(1,3))
print(p4_pca1_2, vp=vplayout(1:2,3))

dev.off()

pdf(file = file.path(outputFolder,"Fig1_ALL_rna_1.pdf"), height =2.5, width = 8)
# 创建一个3x3的网格布局
gl <- grid.layout(nrow=2, ncol=3)
# 将小图放置在网格中
pushViewport(viewport(layout=gl))
vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col=y)
print(p4, vp=vplayout(1:2,1))
print(p4_1, vp=vplayout(1:2,2))
print(p4_pca1_1, vp=vplayout(1:2,3))

dev.off()


pdf(file = file.path(outputFolder,"Fig1_ALL_protein.pdf"), height =2.5, width = 8)
# 创建一个3x3的网格布局
gl <- grid.layout(nrow=2, ncol=3)
# 将小图放置在网格中
pushViewport(viewport(layout=gl))
vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col=y)
print(p2, vp=vplayout(1:2,1))
#print(p4_1, vp=vplayout(1,2))
print(p4_2, vp=vplayout(1:2,2))
#print(p4_pca1_1, vp=vplayout(1,3))
print(p4_pca1_2, vp=vplayout(1:2,3))

#grid.draw()
dev.off()
