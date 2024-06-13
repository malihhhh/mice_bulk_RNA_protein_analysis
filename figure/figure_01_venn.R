library("biomaRt")
library("org.Mm.eg.db")
library('clusterProfiler')

library('reshape')
library('ggplot2')
library('dplyr')
library('gridExtra')
library(xlsx)
library('readxl')

library('ggthemes')
library(VennDiagram)

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

#################################################### 转录组组数据
input_folder <- "/data/home/mali/ke_data/transcription/count/"
transcript_data <- read.csv(paste0(input_folder, "kmeans10.csv"))

log2FoldChange <- c("YV_vs_Y_log2FoldChange","OV_vs_O_log2FoldChange","O_vs_Y_log2FoldChange","OV_vs_YV_log2FoldChange")
logFC_cutoff <- with(transcript_data,mean(abs(YV_vs_Y_log2FoldChange)) + 2*sd(abs(YV_vs_Y_log2FoldChange)))
print(paste0(log2FoldChange[1]," logFC_cutoff: ",logFC_cutoff))
logFC_cutoff <- with(transcript_data,mean(abs(OV_vs_O_log2FoldChange)) + 2*sd(abs(OV_vs_O_log2FoldChange)))
print(paste0(log2FoldChange[2]," logFC_cutoff: ",logFC_cutoff))
logFC_cutoff <- with(transcript_data,mean(abs(O_vs_Y_log2FoldChange)) + 2*sd(abs(O_vs_Y_log2FoldChange)))
print(paste0(log2FoldChange[3]," logFC_cutoff: ",logFC_cutoff))
logFC_cutoff <- with(transcript_data,mean(abs(OV_vs_YV_log2FoldChange)) + 2*sd(abs(OV_vs_YV_log2FoldChange)))
print(paste0(log2FoldChange[4]," logFC_cutoff: ",logFC_cutoff))
###########################
# 转录组 计算差异基因数量； 标准padj < 0.05； |log2FoldChange| > 1
logFC_cutoff <- 1
###########################

# 用一列 *_change 来表示上调，下调或不显著
k1 <- (transcript_data$YV_vs_Y_padj < 0.05)&(transcript_data$YV_vs_Y_log2FoldChange < - logFC_cutoff)
k2 <- (transcript_data$YV_vs_Y_padj < 0.05)&(transcript_data$YV_vs_Y_log2FoldChange > logFC_cutoff)
transcript_data$YV_vs_Y_change <- ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
k1 <- (transcript_data$OV_vs_O_padj < 0.05)&(transcript_data$OV_vs_O_log2FoldChange < - logFC_cutoff)
k2 <- (transcript_data$OV_vs_O_padj < 0.05)&(transcript_data$OV_vs_O_log2FoldChange > logFC_cutoff)
transcript_data$OV_vs_O_change <- ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
k1 <- (transcript_data$O_vs_Y_padj < 0.05)&(transcript_data$O_vs_Y_log2FoldChange < - logFC_cutoff)
k2 <- (transcript_data$O_vs_Y_padj < 0.05)&(transcript_data$O_vs_Y_log2FoldChange > logFC_cutoff)
transcript_data$O_vs_Y_change <- ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
k1 <- (transcript_data$OV_vs_YV_padj < 0.05)&(transcript_data$OV_vs_YV_log2FoldChange < - logFC_cutoff)
k2 <- (transcript_data$OV_vs_YV_padj < 0.05)&(transcript_data$OV_vs_YV_log2FoldChange > logFC_cutoff)
transcript_data$OV_vs_YV_change <- ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))

tj <- data.frame(YV_vs_Y = as.integer(table(transcript_data$YV_vs_Y_change)),
                OV_vs_O = as.integer(table(transcript_data$OV_vs_O_change)),
                O_vs_Y = as.integer(table(transcript_data$O_vs_Y_change)),
                OV_vs_YV = as.integer(table(transcript_data$OV_vs_YV_change)),
            row.names = c("down","not","up"))


#--------上调
up <- list(YV_vs_Y = na.omit(transcript_data$Gene[transcript_data$YV_vs_Y_change=="UP"]),
            OV_vs_O = na.omit(transcript_data$Gene[transcript_data$OV_vs_O_change=="UP"]),
            O_vs_Y = na.omit(transcript_data$Gene[transcript_data$O_vs_Y_change=="UP"]),
            OV_vs_YV = na.omit(transcript_data$Gene[transcript_data$OV_vs_YV_change=="UP"]))
up_dataframe <- do.call(cbind, lapply(lapply(up, unlist), `length<-`, max(lengths(up))))
#组间交集元素获得
inter <- get.venn.partitions(up)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], file.path(outputFolder,'RNA_VENN_up.txt'), row.names = FALSE, sep = '\t', quote = FALSE)

#保存与再读入
write.csv(up_dataframe,file=file.path(outputFolder,"RNA_VENN_up.csv"),row.names = F)
dat <- read.csv(file.path(outputFolder,'RNA_VENN_up.csv'), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
up <- list(YV_vs_Y = na.omit(dat$YV_vs_Y),
                  OV_vs_O = na.omit(dat$OV_vs_O),
                  O_vs_Y = na.omit(dat$O_vs_Y),
                  OV_vs_YV = na.omit(dat$OV_vs_YV))
#绘图
R_venn_up <- venn.diagram(up, filename = NULL,main = "Up regulated genes", main.cex = 2, #imagetype = 'tiff',
    fill = colorScheme$type, alpha = 0.50,
    cat.col = colorScheme$type, cat.cex = 1.5, cat.fontfamily = 'serif',
    col = colorScheme$type, cex = 1.5, fontfamily = 'serif',
    output=TRUE)

pdf(file.path(outputFolder,'RNA_VENN_up.pdf'))
grid.draw(R_venn_up)
dev.off()
#--------

#--------下调
down <- list(YV_vs_Y = na.omit(transcript_data$Gene[transcript_data$YV_vs_Y_change=="DOWN"]),
            OV_vs_O = na.omit(transcript_data$Gene[transcript_data$OV_vs_O_change=="DOWN"]),
            O_vs_Y = na.omit(transcript_data$Gene[transcript_data$O_vs_Y_change=="DOWN"]),
            OV_vs_YV = na.omit(transcript_data$Gene[transcript_data$OV_vs_YV_change=="DOWN"]))
down_dataframe <- do.call(cbind, lapply(lapply(down, unlist), `length<-`, max(lengths(down))))
#组间交集元素获得
inter <- get.venn.partitions(down)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], file.path(outputFolder,'RNA_VENN_down.txt'), row.names = FALSE, sep = '\t', quote = FALSE)

#保存与再读入
write.csv(down_dataframe,file=file.path(outputFolder,"RNA_VENN_down.csv"),row.names = F)
dat <- read.csv(file.path(outputFolder,'RNA_VENN_down.csv'), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
down <- list(YV_vs_Y = na.omit(dat$YV_vs_Y),
                  OV_vs_O = na.omit(dat$OV_vs_O),
                  O_vs_Y = na.omit(dat$O_vs_Y),
                  OV_vs_YV = na.omit(dat$OV_vs_YV))
#绘图
R_venn_down <- venn.diagram(down, filename = NULL, main = "Down regulated genes", main.cex = 2,#imagetype = 'tiff',
    fill = colorScheme$type, alpha = 0.50,
    cat.col = colorScheme$type, cat.cex = 1.5, cat.fontfamily = 'serif',
    col = colorScheme$type, cex = 1.5, fontfamily = 'serif',
    output=TRUE)

pdf(file.path(outputFolder,'RNA_VENN_down.pdf'))
grid.draw(R_venn_up)
dev.off()
#--------
####################################################

#################################################### 蛋白组数据
O_vs_Y <- read_xlsx(paste0("/data/home/mali/ke_data/Protein/", "Differentially_expressed_statistics.xlsx"),sheet=2)
OV_vs_YV <- read_xlsx(paste0("/data/home/mali/ke_data/Protein/", "Differentially_expressed_statistics.xlsx"),sheet=3)
YV_vs_Y <- read_xlsx(paste0("/data/home/mali/ke_data/Protein/", "Differentially_expressed_statistics.xlsx"),sheet=4)
OV_vs_O <- read_xlsx(paste0("/data/home/mali/ke_data/Protein/", "Differentially_expressed_statistics.xlsx"),sheet=5)

tj_protein <- data.frame(O_vs_Y = as.integer(table(O_vs_Y$'Regulated Type')),
                        OV_vs_YV = as.integer(table(OV_vs_YV$'Regulated Type')),
                        YV_vs_Y = as.integer(table(YV_vs_Y$'Regulated Type')),
                        OV_vs_O = as.integer(table(OV_vs_O$'Regulated Type')),
                        row.names = c("Down","Up"))
#--------上调
up <- list(O_vs_Y = na.omit(O_vs_Y$'Protein accession'[O_vs_Y$'Regulated Type'=="Up"]),
            OV_vs_YV = na.omit(OV_vs_YV$'Protein accession'[OV_vs_YV$'Regulated Type'=="Up"]),
            YV_vs_Y = na.omit(YV_vs_Y$'Protein accession'[YV_vs_Y$'Regulated Type'=="Up"]),
            OV_vs_O = na.omit(OV_vs_O$'Protein accession'[OV_vs_O$'Regulated Type'=="Up"]))
up_dataframe <- do.call(cbind, lapply(lapply(up, unlist), `length<-`, max(lengths(up))))
#组间交集元素获得
inter <- get.venn.partitions(up)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], file.path(outputFolder,'Protein_VENN_up.txt'), row.names = FALSE, sep = '\t', quote = FALSE)

# 蛋白转换的geneID
up_gene <- list(O_vs_Y = na.omit(O_vs_Y$'Gene name'[O_vs_Y$'Regulated Type'=="Up"]),
            OV_vs_YV = na.omit(OV_vs_YV$'Gene name'[OV_vs_YV$'Regulated Type'=="Up"]),
            YV_vs_Y = na.omit(YV_vs_Y$'Gene name'[YV_vs_Y$'Regulated Type'=="Up"]),
            OV_vs_O = na.omit(OV_vs_O$'Gene name'[OV_vs_O$'Regulated Type'=="Up"]))
word_to_remove <- "--"
# 循环遍历列表中的每个元素
for (i in seq_along(up_gene)) {
  # 使用逻辑判断来删除包含指定单词的元素
  up_gene[[i]] <- up_gene[[i]][up_gene[[i]] != word_to_remove]
}
#组间交集元素获得
inter <- get.venn.partitions(up_gene)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], file.path(outputFolder,'Protein_VENN_up_gene.txt'), row.names = FALSE, sep = '\t', quote = FALSE)


#保存与再读入
write.csv(up_dataframe,file=file.path(outputFolder,"Protein_VENN_up.csv"),row.names = F)
dat <- read.csv(file.path(outputFolder,'Protein_VENN_up.csv'), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
up <- list(YV_vs_Y = na.omit(dat$YV_vs_Y),
                  OV_vs_O = na.omit(dat$OV_vs_O),
                  O_vs_Y = na.omit(dat$O_vs_Y),
                  OV_vs_YV = na.omit(dat$OV_vs_YV))
#绘图
P_venn_up <- venn.diagram(up, filename = NULL,main = "Up regulated genes", main.cex = 2, #imagetype = 'tiff',
    fill = colorScheme$type, alpha = 0.50,
    cat.col = colorScheme$type, cat.cex = 1.5, cat.fontfamily = 'serif',
    col = colorScheme$type, cex = 1.5, fontfamily = 'serif',
    output=TRUE)

pdf(file.path(outputFolder,'Protein_VENN_up.pdf'))
grid.draw(P_venn_up)
dev.off()
#--------

#--------下调
# 蛋白ID
down <- list(O_vs_Y = na.omit(O_vs_Y$'Protein accession'[O_vs_Y$'Regulated Type'=="Down"]),
            OV_vs_YV = na.omit(OV_vs_YV$'Protein accession'[OV_vs_YV$'Regulated Type'=="Down"]),
            YV_vs_Y = na.omit(YV_vs_Y$'Protein accession'[YV_vs_Y$'Regulated Type'=="Down"]),
            OV_vs_O = na.omit(OV_vs_O$'Protein accession'[OV_vs_O$'Regulated Type'=="Down"]))
down_dataframe <- do.call(cbind, lapply(lapply(down, unlist), `length<-`, max(lengths(down))))
#组间交集元素获得
inter <- get.venn.partitions(down)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], file.path(outputFolder,'Protein_VENN_down.txt'), row.names = FALSE, sep = '\t', quote = FALSE)

# 蛋白转换的geneID
down_gene <- list(O_vs_Y = na.omit(O_vs_Y$'Gene name'[O_vs_Y$'Regulated Type'=="Down"]),
            OV_vs_YV = na.omit(OV_vs_YV$'Gene name'[OV_vs_YV$'Regulated Type'=="Down"]),
            YV_vs_Y = na.omit(YV_vs_Y$'Gene name'[YV_vs_Y$'Regulated Type'=="Down"]),
            OV_vs_O = na.omit(OV_vs_O$'Gene name'[OV_vs_O$'Regulated Type'=="Down"]))
word_to_remove <- "--"
# 循环遍历列表中的每个元素
for (i in seq_along(down_gene)) {
  # 使用逻辑判断来删除包含指定单词的元素
  down_gene[[i]] <- down_gene[[i]][down_gene[[i]] != word_to_remove]
}
#组间交集元素获得
inter <- get.venn.partitions(down_gene)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], file.path(outputFolder,'Protein_VENN_down_gene.txt'), row.names = FALSE, sep = '\t', quote = FALSE)

#保存与再读入
write.csv(down_dataframe,file=file.path(outputFolder,"Protein_VENN_down.csv"),row.names = F)
dat <- read.csv(file.path(outputFolder,'Protein_VENN_down.csv'), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
down <- list(YV_vs_Y = na.omit(dat$YV_vs_Y),
                  OV_vs_O = na.omit(dat$OV_vs_O),
                  O_vs_Y = na.omit(dat$O_vs_Y),
                  OV_vs_YV = na.omit(dat$OV_vs_YV))
#绘图
P_venn_down <- venn.diagram(down, filename = NULL,main = "down regulated genes", main.cex = 2, #imagetype = 'tiff',
    fill = colorScheme$type, alpha = 0.50,
    cat.col = colorScheme$type, cat.cex = 1.5, cat.fontfamily = 'serif',
    col = colorScheme$type, cex = 1.5, fontfamily = 'serif',
    output=TRUE)

pdf(file.path(outputFolder,'Protein_VENN_down.pdf'))
grid.draw(P_venn_down)
dev.off()
#--------


###################------------ 差异表达基因数量柱状图
deg_data <- c()
name <- c("YV_vs_Y","OV_vs_O","O_vs_Y","OV_vs_YV")
for (i in 1:length(name)) {
  col_padj <- paste0(name[i],"_padj")
  col_FC <- paste0(name[i],"_log2FoldChange")
  result <- base::subset(transcript_data, as.numeric(transcript_data[,col_padj]) < 0.05 & abs(transcript_data[,col_FC]) > 1)
  uplen <- nrow(subset(result, result[,col_FC] > 0))
  downlen <- -nrow(subset(result, result[,col_FC] < 0))
  tem <- data.frame(group = name[i], Up_regulated = uplen, Down_regulated = downlen)
  deg_data <- rbind(deg_data,tem)
}
deg_data$group  <-  factor(deg_data$group, levels = c("O_vs_Y","OV_vs_YV","YV_vs_Y","OV_vs_O"))

# 蛋白组和代谢组有整理好的差异基因数量表格
# 蛋白组公司标准（FoldChange>1.2)
# 代谢组公司标准（Fold Change ≤ 0.5 和Fold Change ≤ 0.5；VIP ≥ 1  ）
protion_data <- read_xlsx(paste0("/data/home/mali/ke_data/Protein/", "Differentially_expressed_statistics.xlsx"),sheet=1)

colnames(protion_data) <- c("group","Up_regulated","Down_regulated")
protion_data$group <- c("O_vs_Y","OV_vs_YV","YV_vs_Y","OV_vs_O")
protion_data$group  <-  factor(protion_data$group, levels = c("O_vs_Y","OV_vs_YV","YV_vs_Y","OV_vs_O"))
protion_data$Down_regulated <- -protion_data$Down_regulated


# RNA table
df_a <- as.data.frame(c(441,129,106,24,149,49,225,110))
colnames(df_a) <- "DEG_of_RNA"
df_a$Type <- c("YV_vs_Y","YV_vs_Y","OV_vs_O","OV_vs_O","O_vs_Y","O_vs_Y","OV_vs_YV","OV_vs_YV")
df_a$Peptide <- c("up","down","up","down","up","down","up","down")

#protein table
df_b <- as.data.frame(c(269,154,331,204,280,283,367,352))
colnames(df_b) <- "DEG_of_protein"
df_b$Type <- c("YV_vs_Y","YV_vs_Y","OV_vs_O","OV_vs_O","O_vs_Y","O_vs_Y","OV_vs_YV","OV_vs_YV")
df_b$Peptide <- c("up","down","up","down","up","down","up","down")


#figures
theme_set(all_theme())
width <- 4
height <- 4.5

p_peps <- ggplot(df_a,aes(factor(Type),DEG_of_RNA,fill=factor(Peptide,c("up","down"))))+
  geom_bar(stat="identity")+
  geom_text(label= df_a$DEG_of_RNA, color = "black", vjust=-1)+
  scale_fill_manual(values = c("#742c55","#443158"),name="Quantification")+
  xlab(element_blank())+
  scale_y_continuous(labels = comma)+
 theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=0),legend.position ="none")

ggsave(
  filename = "bars_DEG_RNA.pdf",
  plot = p_peps,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = width,
  height = height+0.7,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)


p_proteins <- ggplot(df_b,aes(factor(Type),DEG_of_protein,fill=factor(Peptide,c("up","down"))))+
  geom_bar(stat="identity")+
  geom_text(label= df_b$DEG_of_protein, color = "black", vjust=0)+
  scale_fill_manual(values = c("#742c55","#443158"),name="Quantification")+
  xlab(element_blank())+
  scale_y_continuous(labels = comma)+
 theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=0),legend.position ="none")

ggsave(
  filename = "bars_DEG_protein.pdf",
  plot = p_proteins,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = width,
  height = height+0.7,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
