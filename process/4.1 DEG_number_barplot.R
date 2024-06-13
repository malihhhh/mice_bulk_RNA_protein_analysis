library("biomaRt")
library("org.Mm.eg.db")
library('clusterProfiler')
library('readxl')
library('reshape')
library('ggplot2')
library('dplyr')
library('gridExtra')
library(xlsx)
##########==================Code to generate bar plot with number of DEGs as in figure S3A----
compare_name = "DEG_number_barplot"
n_cores <- 4 # number of cores to use for parallelization (8 cores crashed , 4 cores took 6 hours)

useVersion   <- "version1"
outputFolder <- file.path("/data/home/mali/ke_data/output", useVersion, "DEG_number_barplot")
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)


setwd("/data/home/mali/ke_data")

# 转录组组数据
input_folder <- "/data/home/mali/ke_data/transcription/count/"
transcript_data <- read.csv(paste0(input_folder, "kmeans10.csv"))
# 转录组 计算差异基因数量； 标准padj < 0.05； |log2FoldChange| > 1
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
pst_protion_data <- read_xlsx(paste0("/data/home/mali/ke_data/PST_Protein/", "Differentially_expressed_statistics.xlsx"),sheet=1)
plasma_metabolome_data <- read_xlsx(paste0("/data/home/mali/ke_data/plasma_metabolome/", "sigMetabolitesCount.xlsx"),sheet=1)
lung_metabolome_data <- read_xlsx(paste0("/data/home/mali/ke_data/metabolome/", "sigMetabolitesCount.xlsx"),sheet=1)

colnames(protion_data) <- c("group","Up_regulated","Down_regulated")
protion_data$group <- c("O_vs_Y","OV_vs_YV","YV_vs_Y","OV_vs_O")
protion_data$group  <-  factor(protion_data$group, levels = c("O_vs_Y","OV_vs_YV","YV_vs_Y","OV_vs_O"))
protion_data$Down_regulated <- -protion_data$Down_regulated

colnames(pst_protion_data) <- c("group","Type","Up_regulated","Down_regulated")
pst_protion_data$group <- rep(c("O_vs_Y","OV_vs_YV","YV_vs_Y","OV_vs_O"),each=2)
pst_protion_data$group  <-  factor(pst_protion_data$group, levels = c("O_vs_Y","OV_vs_YV","YV_vs_Y","OV_vs_O"))
pst_protion_data$Down_regulated <- -pst_protion_data$Down_regulated

pst_protion_data_site <- pst_protion_data[c(1,3,5,7),c(1,3,4)]
pst_protion_data_protein <- pst_protion_data[c(2,4,6,8),c(1,3,4)]


plasma_metabolome_data <- plasma_metabolome_data[,c(1,4,3)] # 代谢组比对的是反的，把上调和下调的反着来
colnames(plasma_metabolome_data) <- c("group","Up_regulated","Down_regulated")
plasma_metabolome_data$group <- c("O_vs_Y","OV_vs_YV","YV_vs_Y","OV_vs_O")
plasma_metabolome_data$group  <-  factor(plasma_metabolome_data$group, levels = c("O_vs_Y","OV_vs_YV","YV_vs_Y","OV_vs_O"))
plasma_metabolome_data$Down_regulated <- -plasma_metabolome_data$Down_regulated

lung_metabolome_data <- lung_metabolome_data[,c(1,4,3)] # 代谢组比对的是反的，把上调和下调的反着来
colnames(lung_metabolome_data) <- c("group","Up_regulated","Down_regulated")
lung_metabolome_data$group <- c("O_vs_Y","OV_vs_YV","YV_vs_Y","OV_vs_O")
lung_metabolome_data$group  <-  factor(lung_metabolome_data$group, levels = c("O_vs_Y","OV_vs_YV","YV_vs_Y","OV_vs_O"))
lung_metabolome_data$Down_regulated <- -lung_metabolome_data$Down_regulated


##Reformat data
deg_data_plot <- melt(deg_data %>% as.data.frame)
protion_data_plot <- melt(protion_data %>% as.data.frame)
pst_protion_data_site_plot <- melt(pst_protion_data_site %>% as.data.frame)
pst_protion_data_protein_plot <- melt(pst_protion_data_protein %>% as.data.frame)
plasma_metabolome_data_plot <- melt(plasma_metabolome_data %>% as.data.frame)
lung_metabolome_data_plot <- melt(lung_metabolome_data %>% as.data.frame)
range(deg_data_plot$value)
range(protion_data_plot$value)
range(pst_protion_data_site_plot$value)
range(pst_protion_data_protein_plot$value)
range(plasma_metabolome_data_plot$value)
range(lung_metabolome_data_plot$value)

##--------------------------------------Generate bar plot



P1 <- ggplot(data = deg_data_plot) + geom_col(aes(x=group, y=value, fill=variable), width = 0.8, size = 0.1, colour = 'black')+ coord_flip()+ theme(legend.position='none') + ylim(-600,600)
P1 <- P1 + scale_fill_manual(values=c('#990000', '#004C99', '#994C99'))
P1 <- P1 + xlab("Group") + ylab("Number of significant DEG") + scale_x_discrete(labels = deg_data_plot$group)
P1 <- P1 +theme(axis.text=element_text(size=7, color = "black"), panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


P2 <- ggplot(data = lung_metabolome_data_plot) + geom_col(aes(x=group, y=value, fill=variable), width = 0.8, size = 0.1, colour = 'black')+ coord_flip()+ theme(legend.position='none') +ylim(-600,600)
P2 <- P2 + scale_fill_manual(values=c('#990000', '#004C99', '#994C99'))
P2 <- P2 + xlab(NULL)+ ylab(NULL) + scale_x_discrete(labels = lung_metabolome_data_plot$group)
P2 <- P2 +theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=7, color = "black"),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


P3 <- ggplot(data = plasma_metabolome_data_plot) + geom_col(aes(x=group, y=value, fill=variable), width = 0.8, size = 0.1, colour = 'black')+ coord_flip()+ theme(legend.position='none') + ylim(-600,600)
P3 <- P3 + scale_fill_manual(values=c('#990000', '#004C99', '#994C99'))
P3 <- P3 + xlab(NULL)+ ylab(NULL) + scale_x_discrete(labels = plasma_metabolome_data_plot$group)
P3 <- P3 +theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=7, color = "black"),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


P4 <- ggplot(data = protion_data_plot) + geom_col(aes(x=group, y=value, fill=variable), width = 0.8, size = 0.1, colour = 'black')+ coord_flip()+ theme(legend.position='none') + ylim(-600,600)
P4 <- P4 + scale_fill_manual(values=c('#990000', '#004C99', '#994C99'))
P4 <- P4 + xlab(NULL)+ ylab(NULL) + scale_x_discrete(labels = pst_protion_data_site$group)
P4 <- P4 +theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=7, color = "black"),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


P5 <- ggplot(data = pst_protion_data_site_plot) + geom_col(aes(x=group, y=value, fill=variable), width = 0.8, size = 0.1, colour = 'black')+ coord_flip()+ theme(legend.position='none') +  ylim(-1100,1100)
P5 <- P5 + scale_fill_manual(values=c('#990000', '#004C99', '#994C99'))
P5 <- P5 + xlab(NULL)+ ylab(NULL) + scale_x_discrete(labels = pst_protion_data_protein$group)
P5 <- P5 +theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=7, color = "black"),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


P6 <- ggplot(data = pst_protion_data_protein_plot) + geom_col(aes(x=group, y=value, fill=variable), width = 0.8, size = 0.1, colour = 'black')+ coord_flip()+ theme(legend.position='none') + ylim(-600,600)
P6 <- P6 + scale_fill_manual(values=c('#990000', '#004C99', '#994C99'))
P6 <- P6 + xlab(NULL)+ ylab(NULL) + scale_x_discrete(labels = pst_protion_data_protein$group)
P6 <- P6 +theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=7, color = "black"),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

pdf(file.path(outputFolder,"deg_impulse_number_comparison.pdf"), width = 6, height = 3)
par(mfrow = c(1,6))
P1 + P2 + P3 + P4 + P5 + P6
dev.off()

#g <- grid.arrange(P1,
#             P2,
#             P3,
#             P4,
#             P5,
#             P6,
#             nrow = 1,
#             ncol = 6,
#             widths = c(3, 3, 3, 3,3, 3) %>% {. / sum(.)})
#ggsave(file.path(outputFolder,"deg_impulse_number_comparison.pdf"), plot = g, width = 10, height = 3)


# 保存source data
data_all <- cbind(deg_data_plot,lung_metabolome_data_plot,plasma_metabolome_data_plot,protion_data_plot,protion_data_plot,pst_protion_data_site_plot,pst_protion_data_protein_plot)
colnames(data_all) <- paste0(rep(c("deg_data_plot","lung_metabolome_data_plot","plasma_metabolome_data_plot","protion_data_plot","protion_data_plot","pst_protion_data_site_plot","pst_protion_data_protein_plot"),each=3),"_",rep(colnames(deg_data_plot),6))
write.xlsx(x = data_all, file = file.path(outputFolder,"FIG1_source_data.xlsx"),
        sheetName = "figb", row.names = FALSE)









###------------------------------------韦恩图
library(VennDiagram)
setwd("/data/home/mali/bd_scRNAseq_PBMC_SARSCOV2_patient_Guan/child/DEG/Basic")
input_folder <- "/data/home/mali/bd_scRNAseq_PBMC_SARSCOV2_patient_Guan/child/DEG/Basic/"
file_names <- dir("./", pattern = glob2rx("*.csv"))

#------------------轻症比健康对照
comparison <- "_"
compartment <-  "markers"
files <- file_names[grep(comparison,file_names)]
path_save <- "/data/home/mali/bd_scRNAseq_PBMC_SARSCOV2_patient_Guan/child/Figure/"

deg_data <- data.frame(timepoint=as.character(), upgenes=as.integer(), downgenes=as.integer())
levels(covsc_combined)
celltype <- c("T_CD4","T_CD4_naive","T_CD8","T_CD8_naive","T_CD8_memory","T_DN","gdT","NK_CD16","NK_CD56","NKT","MO_CD14","MO_CD16","mDC","B","pDC","Plasmablasts")

##Get number of DEGs at each pseudotime compared to controls
for (i in 1:length(celltype)) {
  comparison <- "C_H"
  file_import <- grep(paste0(comparison,"_",celltype[[i]],"_",compartment,".csv"),files,value = T)
  result1 <- read.csv(paste0(input_folder, file_import))
  result1 <- subset(result1, result1$p_val_adj < 0.05 & abs(result1$avg_log2FC) > 0.25)
  rownames(result1) <- result1$gene
  uplen1 <- rownames(subset(result1, result1$avg_log2FC > 0))
  downlen1 <- rownames(subset(result1, result1$avg_log2FC < 0))
  comparison <- "C_R"
  file_import <- grep(paste0(comparison,"_",celltype[[i]],"_",compartment,".csv"),files,value = T)
  result2 <- read.csv(paste0(input_folder, file_import))
  result2 <- subset(result2, result2$p_val_adj < 0.05 & abs(result2$avg_log2FC) > 0.25)
  rownames(result2) <- result2$gene
  uplen2 <- rownames(subset(result2, result2$avg_log2FC > 0))
  downlen2 <- rownames(subset(result2, result2$avg_log2FC < 0))
  comparison <- "R_H"
  file_import <- grep(paste0(comparison,"_",celltype[[i]],"_",compartment,".csv"),files,value = T)
  result3 <- read.csv(paste0(input_folder, file_import))
  result3 <- subset(result3, result3$p_val_adj < 0.05 & abs(result3$avg_log2FC) > 0.25)
  rownames(result3) <- result3$gene
  uplen3 <- rownames(subset(result3, result3$avg_log2FC > 0))
  downlen3 <- rownames(subset(result3, result3$avg_log2FC < 0))
  venn.diagram(list( up_C_H=uplen1 ,
                   up_C_R=uplen2,
                   up_R_H=uplen3 ),
             fill=c("red","green","blue"),
             filename=paste0(path_save,celltype[i],"_up_venn.tiff"))
  venn.diagram(list( down_C_H=downlen1 ,
                   down_C_R=downlen2,
                   down_R_H=downlen3 ),
             fill=c("red","green","blue"),
             filename=paste0(path_save,celltype[i],"_down_venn.tiff"))
  #deg_data <- rbind(deg_data, data.frame(celltype[i], uplen, downlen))
}


#########3==============================Plot volcano plot for deseq result from non survivors and color genes by cohort 1 ##############==========================Code to generate volcano plot------------
#-------------数据准备
library(ggrepel)
library(ggplot2)

# 转录组组数据
input_folder <- "/data/home/mali/ke_data/transcription/count/"
transcript_data <- read.csv(paste0(input_folder, "kmeans10.csv"))
# 转录组 计算差异基因数量； 标准padj < 0.05； |log2FoldChange| > 1
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

### Read each deseq out file and filter for significant differentially expressed genes
plot_data <- data.frame()
for (i in 1:length(celltype)) {
  file_import <- grep(paste0("*_",celltype[[i]],"_",compartment,".csv"),files,value = T)
  result <- read.csv(paste0(input_folder, file_import))
  result <- subset(result, result$p_val_adj < 0.05 & abs(result$avg_log2FC) > 0.5)
  result$celltype <- celltype[i]
  result <- rownames_to_column(result, "Gene_ID")
  plot_data <- rbind(plot_data, result)

}

#Add direction of regulation to plot_data
plot_data$direction <- ifelse(plot_data$avg_log2FC < 0, "down", "up")

#Calculate the number of comparisons in which a gene is significant
gene_frequency <- table(plot_data$Gene_ID)
plot_data$frequency <- as.vector(gene_frequency[as.character(plot_data$Gene_ID)])

#Select interesting genes to label
INFA.list <- read.table("/data/home/mali/bd_scRNAseq_PBMC_SARSCOV2_patient_Guan/mild_cmild/geneset/HALLMARK_INTERFERON_ALPHA_RESPONSE.txt",header= F,skip=2)
INFG.list <- read.table("/data/home/mali/bd_scRNAseq_PBMC_SARSCOV2_patient_Guan/mild_cmild/geneset/HALLMARK_INTERFERON_GAMMA_RESPONSE.txt",header= F,skip=2)
TNF.list <- read.table("/data/home/mali/bd_scRNAseq_PBMC_SARSCOV2_patient_Guan/mild_cmild/geneset/geneset_TNF_NFKB.txt",header= F,skip=2)
inf.list <- read.table("/data/home/mali/bd_scRNAseq_PBMC_SARSCOV2_patient_Guan/mild_cmild/geneset/geneset_inflammation.txt",header= F,skip=2)

gene_to_label_INFA <- list(as.character(INFA.list$V1))
gene_to_label_INFG <- list(as.character(INFG.list$V1))
gene_to_label_TNF <- list(as.character(TNF.list$V1))
gene_to_label_inf <- list(as.character(inf.list$V1))

gene_to_label <- c("HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA−DOB", "HLA−DPA1", "HLA-DPB1", "HLA−DQA1", "HLA-DQA2", "HLA-DQB1", "HLA−DQB2", "HLA-DRA", "HLA-DRB1", "HLA-DRB5")
gene_to_label2 <- c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G")


#-------------数据准备完成

#----------Code to plot 可以在一个火山图上画多种细胞类群或者其他分类方法
#cairo_ps(file = "deg_celltype_lfc.ps", width = 9, height = 7)
p <- ggplot(data = plot_data, aes(x=celltype, y=avg_log2FC, fill = direction))
p <- p + geom_point(position=position_jitter(width=0.1), aes(size=-1*log10(p_val_adj), alpha=frequency), pch=21, color="black", stroke=0.1)
p <- p + geom_text_repel(data = subset(plot_data, plot_data$gene %in% gene_to_label & abs(plot_data$avg_log2FC) > 2), aes(x=celltype, y=avg_log2FC, label=gene),
                         color="black", size=3.5, nudge_x = -0.5)
p <- p + geom_text_repel(data = subset(plot_data, plot_data$gene %in% gene_to_label2& abs(plot_data$avg_log2FC) > 1), aes(x=celltype, y=avg_log2FC, label=gene),
                         color="black", size=3.5, nudge_x = -0.5)
p <- p + scale_color_manual(values=c('#004C99', '#990000')) + scale_fill_manual(values=c('#004C99', '#990000'))
p <- p + xlab("celltype") + ylab("Log fold change") + scale_x_discrete(labels = c("0 vs 1", "0 vs 2", "0 vs 3", "0 vs 4", "0 vs 5", "0 vs 6"))
p <- p + theme_bw() + theme(axis.text=element_text(size=14, color = "black"), axis.title=element_text(size=16), legend.text = element_text(size=12),
                            plot.title = element_text(size=20, face="bold", hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))
p
dev.off()
#----------

#-------------------------------------目前使用的普通火山图绘制
#划分上下调显著
plot_data$Significant <- ifelse(plot_data$p_val_adj < 0.05 & abs(plot_data$avg_log2FC) >= 1,
                            ifelse(plot_data$avg_log2FC > 1, "Up", "Down"), "Stable")

pdf(paste0(path_save,celltype,"_Volcano_plot.pdf"), height=7, width=7)
p <- ggplot(data = plot_data,aes(x=avg_log2FC, y=-log10(p_val_adj),color = Significant)) +
  geom_point(alpha=0.7, size = 1,shape=16)
p <- p + geom_vline(xintercept = 0, lty=3) + geom_vline(xintercept = 1, lty=4) + geom_vline(xintercept = -1, lty=4)
p <- p + geom_hline(yintercept = 3, lty=4)
p <- p + scale_color_manual(values = c("blue","grey", "red")) + scale_alpha_manual(values = c(0.3, 1))
p <- p + xlab("Expression (Log2 Fold Change)") + ylab("- Log10 p_val_adj")
p <- p + theme_bw() + theme(axis.text = element_text(size = 12, color="black"), axis.title = element_text(size = 12))

# 将需要标记的基因放置在label列(avg_log2FC >= 2，p_val_adj< 0.001)
# 展示avg_log2FC和p_val_adj筛选后的基因
plot_data$label <- ifelse(plot_data$p_val_adj< 0.001 & abs(plot_data$avg_log2FC) >= 1.5,
                      as.character(plot_data$gene), "")

#带框
#p + geom_label_repel(data = plot_data, aes(x = avg_log2FC,y = -log10(p_val_adj),  label = label),size = 3,  show.legend = FALSE)
#不带框
p + geom_text_repel(data = plot_data, aes(x = avg_log2FC,y = -log10(p_val_adj),  label = label),size = 3, color="black", show.legend = FALSE)

dev.off()


#MHC-II类分子是否存在
p <- p + geom_text_repel(data = subset(plot_data, plot_data$gene %in% gene_to_label & abs(plot_data$avg_log2FC) > 0.5), aes(x=plot_data$avg_log2FC, y = -log10(plot_data$p_val_adj), label=gene),
                         color="black", size=3.5, nudge_x = -0.5)
#MHC-I类分子是否存在
p <- p + geom_text_repel(data = subset(plot_data, plot_data$gene %in% gene_to_label2 & abs(plot_data$avg_log2FC) > 0.5), aes(x=plot_data$avg_log2FC, y = -log10(plot_data$p_val_adj), label=gene),
                         color="black", size=3.5, nudge_x = -0.5)
#INFA基因集基因是否存在
p <- p + geom_text_repel(data = subset(plot_data, plot_data$gene %in% gene_to_label_INFA & abs(plot_data$avg_log2FC) > 0.5), aes(x=plot_data$avg_log2FC, y = -log10(plot_data$p_val_adj), label=gene),
                         color="black", size=3.5, nudge_x = -0.5)
#INFG基因集基因是否存在
p <- p + geom_text_repel(data = subset(plot_data, plot_data$gene %in% gene_to_label_INFG & abs(plot_data$avg_log2FC) > 0.5), aes(x=plot_data$avg_log2FC, y = -log10(plot_data$p_val_adj), label=gene),
                         color="black", size=3.5, nudge_x = -0.5)
#TNF基因集基因是否存在
p <- p + geom_text_repel(data = subset(plot_data, plot_data$gene %in% gene_to_label_TNF & abs(plot_data$avg_log2FC) > 0.5), aes(x=plot_data$avg_log2FC, y = -log10(plot_data$p_val_adj), label=gene),
                         color="black", size=3.5, nudge_x = -0.5)
#INF基因集基因是否存在
p <- p + geom_text_repel(data = subset(plot_data, plot_data$gene %in% gene_to_label_inf & abs(plot_data$avg_log2FC) > 0.5), aes(x=plot_data$avg_log2FC, y = -log10(plot_data$p_val_adj), label=gene),
                         color="black", size=3.5, nudge_x = -0.5)



####### ======================== DO terms plot ====================================
# select top 10 GO terms
```{r}
library(tidyverse)
library(forcats) #forecats 库是 tidyverse 中的一个库，专门用于处理 R 中的因子
library(dplyr)
library(stringdist)
library(GOSemSim)

celltype <- "pDC"
#setwd("/data/home/mali/bd_scRNAseq_PBMC_SARSCOV2_patient_Guan/child/DEG/GO_DOWN")  #下调
#input_folder <- "/data/home/mali/bd_scRNAseq_PBMC_SARSCOV2_patient_Guan/child/DEG/GO_DOWN/" #下调
setwd("/data/home/mali/bd_scRNAseq_PBMC_SARSCOV2_patient_Guan/child/DEG/GO_UP") #上调
input_folder <- "/data/home/mali/bd_scRNAseq_PBMC_SARSCOV2_patient_Guan/child/DEG/GO_UP/" #上调
file_names <- dir("./", pattern = glob2rx("*.csv"))
comparison <- "C_H"
compartment <-  "markers_enrichGO"
files <- file_names[grep(comparison,file_names)]
file_import <- grep(paste0("*_",celltype,"_",compartment,".csv"),files,value = T)
result <- read.csv(paste0(input_folder, file_import))


#根据p adj值进行排序
topGOResults_top <- result[order(result$p.adjust),]
topGOResults_top<-head(topGOResults_top,50)

#筛选相似度，BP
topGOResults_top <- GO_SIM(topGOResults_top,hsGO,sim_value=0.8)
#去重，根据geneID
#topGOResults_top <- topGOResults_top[!duplicated(topGOResults_top$geneID),]
topGOResults_top <- GO_SIM_gene(topGOResults_top,sim_value=15)

#选取前多少个展示
topGOResults_top <- head(topGOResults_top,10)

plot_data <- separate(data = topGOResults_top, col = GeneRatio, into = c("gene1", "gene2"), sep = "/")
plot_data$GeneRatio <- as.integer(plot_data$gene1)/as.integer(plot_data$gene2)

####  绘点图
#pdf(paste0(path_save,celltype,"_GO_DOWN_plot.pdf"),width=5.2,height=3) #下调
pdf(paste0(path_save,celltype,"_GO_plot.pdf"),width=5.2,height=3) #上调
 p <- plot_data %>%
  mutate(Description = fct_reorder(Description, GeneRatio)) %>%
  ggplot(aes(x=GeneRatio, y=Description, size=Count, color=-log10(p.adjust))) + geom_point() + scale_colour_gradient(high="#990000", low="#FF9999", name='-log10(P-value)')  + xlim(min(range(plot_data$GeneRatio))-0.01,max(range(plot_data$GeneRatio))+0.01)
 p <- p + theme(legend.title=element_text(size=8, color = "black"),legend.key.size = unit(10, "pt"))
p <- p + theme(axis.text.x = element_text(size=7, color = "black", angle=45, hjust=1), axis.text.y = element_text(size=8, color = "black"),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p
dev.off()

##### 挑选感兴趣的 GO term 时使用
# GOs of interest
GO_interest_MO_CD14 <- c("regulation of tumor necrosis factor production", "neutrophil activation", "tumor necrosis factor superfamily cytokine production", "phagocytosis", "leukocyte cell-cell adhesion", "actin polymerization or depolymerization","superoxide anion generation","regulation of leukocyte degranulation", "positive regulation of cell activation", "neutrophil mediated immunity", "positive regulation of leukocyte cell-cell adhesion", "antigen processing and presentation", "antigen processing and presentation of peptide antigen","immune response-activating cell surface receptor signaling pathway")
plot_data <- result[which(result$Description %in% GO_interest_MO_CD14),]

plot_data <- separate(data = plot_data, col = GeneRatio, into = c("gene1", "gene2"), sep = "/")
plot_data$GeneRatio <- as.integer(plot_data$gene1)/as.integer(plot_data$gene2)
####### ============================================================


####### ======================== DO terms plot ====================================
# select top 10 GO terms
```{r}
library(tidyverse)
library(forcats) #forecats 库是 tidyverse 中的一个库，专门用于处理 R 中的因子
library(dplyr)
library(stringdist)
library(GOSemSim)
hsGO <- godata('org.Mm.eg.db', ont="BP")

gene_result <- read.csv("/data/home/mali/ke_data/output/version1/DEG_number_barplot/venn_OV_VS_YV_DOWN_GNEN.csv")
protein_result <- read.csv("/data/home/mali/ke_data/output/version1/DEG_number_barplot/venn_OV_VS_YV_DOWN_protein.csv")

result <- gene_result
#根据p adj值进行排序
topGOResults_top <- result[order(result$p.adjust),]
topGOResults_top<-head(topGOResults_top,50)
#筛选相似度，BP
topGOResults_top <- GO_SIM(topGOResults_top,hsGO,sim_value=0.8)
x = topGOResults_top
sim_value=0.8

                    Y <- data.frame()
                    Y <- x #[which(x$ONTOLOGY == "BP"),]
                    go1 <- Y$ID
                    go2 <- go1
                    sim_out <- mgoSim(go1, go2, semData=hsGO, measure="Wang", combine=NULL)
                    sim_out[upper.tri(sim_out)] <- 0
                    diag(sim_out) <- 0
                    col_max <-  apply(sim_out,2,max)
                    #quantile(col_max)
                    remain_ID <- names(which(col_max < sim_value)) #取相似度大于0.8的
                    Y <- subset(Y,subset= ID %in% remain_ID)
topGOResults_top <- Y


#去重，根据geneID
#topGOResults_top <- topGOResults_top[!duplicated(topGOResults_top$geneID),]
topGOResults_top <- GO_SIM_gene(topGOResults_top,sim_value=15)

#选取前多少个展示
topGOResults_top <- head(topGOResults_top,10)

plot_data <- separate(data = topGOResults_top, col = GeneRatio, into = c("gene1", "gene2"), sep = "/")
plot_data$GeneRatio <- as.integer(plot_data$gene1)/as.integer(plot_data$gene2)

####  绘点图
#pdf(paste0(path_save,celltype,"_GO_DOWN_plot.pdf"),width=5.2,height=3) #下调
pdf("/data/home/mali/ke_data/output/version1/DEG_number_barplot/venn_OV_VS_YV_DOWN_GNEN_GO_plot.pdf",width=5.2,height=3) #上调
 p <- plot_data %>%
  mutate(Description = fct_reorder(Description, GeneRatio)) %>%
  ggplot(aes(x=GeneRatio, y=Description, size=Count, color=-log10(p.adjust))) + geom_point() + scale_colour_gradient(high="#990000", low="#FF9999", name='-log10(P-value)')  + xlim(min(range(plot_data$GeneRatio))-0.01,max(range(plot_data$GeneRatio))+0.01)
 p <- p + theme(legend.title=element_text(size=8, color = "black"),legend.key.size = unit(10, "pt"))
p <- p + theme(axis.text.x = element_text(size=7, color = "black", angle=45, hjust=1), axis.text.y = element_text(size=8, color = "black"),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p
dev.off()

result <- protein_result
#根据p adj值进行排序
topGOResults_top <- result[order(result$p.adjust),]
topGOResults_top<-head(topGOResults_top,50)
#选取前多少个展示
topGOResults_top <- head(topGOResults_top,10)

plot_data <- separate(data = topGOResults_top, col = GeneRatio, into = c("gene1", "gene2"), sep = "/")
plot_data$GeneRatio <- as.integer(plot_data$gene1)/as.integer(plot_data$gene2)

####  绘点图
#pdf(paste0(path_save,celltype,"_GO_DOWN_plot.pdf"),width=5.2,height=3) #下调
pdf("/data/home/mali/ke_data/output/version1/DEG_number_barplot/venn_OV_VS_YV_DOWN_protein_GO_plot.pdf",width=5.2,height=3) #上调
 p <- plot_data %>%
  mutate(Description = fct_reorder(Description, GeneRatio)) %>%
  ggplot(aes(x=GeneRatio, y=Description, size=Count, color=-log10(p.adjust))) + geom_point() + scale_colour_gradient(high="#990000", low="#FF9999", name='-log10(P-value)')  + xlim(min(range(plot_data$GeneRatio))-0.01,max(range(plot_data$GeneRatio))+0.01)
 p <- p + theme(legend.title=element_text(size=8, color = "black"),legend.key.size = unit(10, "pt"))
p <- p + theme(axis.text.x = element_text(size=7, color = "black", angle=45, hjust=1), axis.text.y = element_text(size=8, color = "black"),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p
dev.off()

##### 挑选感兴趣的 GO term 时使用
# GOs of interest
GO_interest_MO_CD14 <- c("regulation of tumor necrosis factor production", "neutrophil activation", "tumor necrosis factor superfamily cytokine production", "phagocytosis", "leukocyte cell-cell adhesion", "actin polymerization or depolymerization","superoxide anion generation","regulation of leukocyte degranulation", "positive regulation of cell activation", "neutrophil mediated immunity", "positive regulation of leukocyte cell-cell adhesion", "antigen processing and presentation", "antigen processing and presentation of peptide antigen","immune response-activating cell surface receptor signaling pathway")
plot_data <- result[which(result$Description %in% GO_interest_MO_CD14),]

plot_data <- separate(data = plot_data, col = GeneRatio, into = c("gene1", "gene2"), sep = "/")
plot_data$GeneRatio <- as.integer(plot_data$gene1)/as.integer(plot_data$gene2)
####### ============================================================

