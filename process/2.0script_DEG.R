
## originally by Yulong Niu
## yulong.niu@hotmail.com
## 

###########################DEGs##################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~
#############################################################
checkZeros <- function(v, threshold) {
  res <- ifelse(sum(v == 0) > threshold, FALSE, TRUE)
  return(res)
}

checkFlg22 <- function(v, threshold) {

  require('magrittr')

  res <- v %>%
    split(c(rep(1,each=4),rep(2 : 4, each = 3))) %>%
    sapply(checkZeros, threshold) %>%
    all

  return(res)
}
##############################################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd('/data/home/mali/ke_data/transcription/')

#######################################################################
# 创建 DESeq数据
#######################################################################
library('tximport')
library('rhdf5')
library('magrittr')
library('DESeq2')
library('tidyverse')

useVersion   <- "version1"
outputFolder <- file.path("output", useVersion, "DESeq_RNA")
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)


#anno <- read_csv('/data/home/mali/genome/mouse/Ensembl_mmu_Anno_AddGeneName.csv',
#                 col_types = cols(Chromosome = col_character())) %>%
#  mutate_all(list(~replace(., is.na(.), '')))


condi <- unique(pData(RNA_eSet)$group)
kres <- exprs(RNA_eSet)
## sampleTable
sampleTable <- data.frame(condition = factor(c(rep(condi[1], each = 4),rep(condi[2:4], each = 3)))) %>%
  set_rownames(pData(RNA_eSet)$name)
sampleTable$condition %<>% relevel(ref = 'Y')
#> sampleTable
#    condition
#YV1        YV
#YV2        YV
#YV3        YV
#YV4        YV
#OV1        OV
#OV2        OV
#OV3        OV
#O1          O
#O2          O
#O3          O
#Y1          Y
#Y2          Y
#Y3          Y

# 创建DESeq数据 # 非整数会报错 some values in assay are not integers
degres <- DESeqDataSetFromMatrix(round(kres), sampleTable, ~condition)

#> degres
#class: DESeqDataSet 
#dim: 46670 13 
#metadata(1): version
#assays(1): counts
#rownames(46670): 20389 22287 ... novel.313 novel.322
#rowData names(0):
#colnames(13): YV1 YV2 ... Y2 Y3
#colData names(1): condition


## remove 0|0|x, 0|0|0
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkFlg22, 1) %>%
  degres[., ]
# 少了很多数据
#> degres
#class: DESeqDataSet 
#dim: 14258 13 
#metadata(1): version
#assays(1): counts
#rownames(14258): 20389 22287 ... 102632773 54722
#rowData names(0):
#colnames(13): YV1 YV2 ... Y2 Y3
#colData names(1): condition

degres <- DESeq(degres)
res <- results(degres)
> head(res)
#log2 fold change (MLE): condition YV vs Y 
#Wald test p-value: condition YV vs Y 
#DataFrame with 6 rows and 6 columns
#           baseMean log2FoldChange     lfcSE      stat      pvalue        padj
#          <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#20389      11453.66      -1.459546  0.168349  -8.66978 4.32945e-18 3.18045e-16
#22287      10099.72      -2.538680  0.234913 -10.80691 3.19259e-27 4.49517e-25
#100503605   6589.40      -0.909074  0.524278  -1.73395 8.29264e-02 3.34238e-01
#15122       5726.30      -1.036623  0.544419  -1.90409 5.68983e-02 2.60529e-01
#16071       4426.04      -1.550906  0.411652  -3.76752 1.64880e-04 2.15287e-03
#17105       3781.03      -1.043459  0.159305  -6.55009 5.75014e-11 2.27206e-09

## count transformation
rld <- rlog(degres) 
# 往往会利用到DESeq2的标准化函数：rlog和vst函数进行标准化，而这两种标准化方式是有一定区别的
# rlog更适用于样本量小于30的数据，相应的vst负责的是样本量大于30的数据
# 实质上我们通过rlog标准化的raw count数据，实质上是先将row count文件做一个log标准化处理，然后对每一个基因做一个负二项回归，其响应变量为基因的表达量，决策变量为不同的分组信息。
# 一般来说，包的作者建议使用VST，这取决于您如何运行代码，rlog可能会过度缩小组之间的更改。
ntd <- normTransform(degres)
##~~~~~~~~~~~~~~~~~~~~~~~~~~hidden batch effect~~~~~~~~~~~~~~~~~~~~~
# 代理变量法（Surrogate variable analysis，SVA）
# 未知的批次
# sva输入的数据需要进行转换，例如取对数（rlog或logCPM）。不过改变原始矩阵
library('sva')
library('ggplot2')

dat <- rld %>%
  assay %>%
  {.[rowMeans(.) > 1, ]}

#> head(dat,2)
#              Mock_1   Mock_2   Mock_3 Mock_Flg22_1 Mock_Flg22_2 Mock_Flg22_3
#AT1G18320.1 4.038552 3.919625 4.071504     3.828454      3.86371     3.894067
#AT5G11100.1 8.253179 8.001945 7.935380     8.443196      8.60240     8.436905
#            SynCom33_Flg22_1 SynCom33_Flg22_2 SynCom33_Flg22_3 SynCom35_Flg22_1
#AT1G18320.1         3.953679         3.875093         3.911458         4.016012
#AT5G11100.1         7.796617         7.881670         7.937094         7.728326
#            SynCom35_Flg22_2 SynCom35_Flg22_3
#AT1G18320.1         3.874353         3.860681
#AT5G11100.1         7.635167         7.612592

mod <- model.matrix(~ condition, colData(degres))
mod0 <- model.matrix(~ 1, colData(degres))

## ## manual detect surrogate variance （手动检测代理方差）
## svnum <- 4
## svseq <- svaseq(dat, mod, mod0, n.sv = svnum)

## auto detect sv （自动检测代理方差）
svobj <- sva(dat, mod, mod0)
svnum <- svobj$sv %>% ncol

svobj$sv %>%
  set_colnames(paste0('sv', seq_len(svnum))) %>%
  as_tibble %>%
  gather(key = 'sv', value = 'value') %>%
  mutate(condition = colData(degres) %>%
           .$condition %>%
           rep(svnum) %>%
           as.character,
         sample = rep(colnames(degres), svnum)) %>%
  mutate(group = paste(sv, condition, sep = '_')) %>%
  ggplot(aes(sample, value, colour = sv, group = group)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90))
ggsave(file.path(outputFolder,'auto_sv.jpg'))
ggsave(file.path(outputFolder,'auto_sv.pdf'))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~DEG Mock vs. 3 conditions~~~~~~~~~~~~~~~~~
degres$sv1 <- svobj$sv[, 1]
degres$sv2 <- svobj$sv[, 2]
degres$sv3 <- svobj$sv[, 3]
degres$sv4 <- svobj$sv[, 4]
design(degres) <- ~sv1 + sv2 + sv3 + sv4 + condition
# design() 这些泛型函数为计数数据集的操作和数据访问提供了基本接口。
# These generic functions provide basic interfaces to operations on and data access to count datasets.

# 很多人以为去除批次效应是要改变你的表达矩阵，新的表达矩阵然后去走差异分析流程，其实大部分的差异分析流程包里面，人家内置好了考虑你的批次效应这样的混杂因素的函数用法设计。例如：
# 构建DESeq2对象时的设计公式：design = ~ batch + conditions

degres <- DESeq(degres)
# 和之前结果不一样了
#> res <- results(degres)
#> head(res)
#log2 fold change (MLE): condition YV vs Y 
#Wald test p-value: condition YV vs Y 
#DataFrame with 6 rows and 6 columns
#           baseMean log2FoldChange     lfcSE      stat      pvalue        padj
#          <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#20389      11453.66       -1.51835  0.251302  -6.04194 1.52273e-09 7.53499e-08
#22287      10099.72       -2.45214  0.290387  -8.44437 3.05694e-17 3.31696e-15
#100503605   6589.40       -1.05057  0.567882  -1.84997 6.43176e-02 3.93241e-01
#15122       5726.30       -1.20691  0.566485  -2.13052 3.31289e-02 2.52419e-01
#16071       4426.04       -1.51110  0.495781  -3.04793 2.30427e-03 3.11129e-02
#17105       3781.03       -1.02250  0.219581  -4.65662 3.21449e-06 9.07515e-05

# 设置对比组
cond <- list(c('YV', 'Y'),
             c('OV', 'O'),
             c('O', 'Y'),
             c('OV', 'YV'))

# 获得差异表达
resRaw <- lapply(cond,
                 function(x) {
                   degres %>%
                     results(contrast = c('condition', x)) %T>%
                     summary %>%
                     as_tibble %>%
                     select(pvalue, padj, log2FoldChange) %>%
                     rename_all(.funs = list(~paste0(paste(x, collapse = '_vs_'), '_', .)))
                 }) %>% bind_cols

#out of 14258 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 609, 4.3%
#LFC < 0 (down)     : 251, 1.8%
#outliers [1]       : 0, 0%
#low counts [2]     : 5252, 37%
#(mean count < 5)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results#
#

#out of 14258 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 179, 1.3%
#LFC < 0 (down)     : 42, 0.29%
#outliers [1]       : 0, 0%
#low counts [2]     : 7464, 52%
#(mean count < 9)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results#
#

#out of 14258 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 283, 2%
#LFC < 0 (down)     : 88, 0.62%
#outliers [1]       : 0, 0%
#low counts [2]     : 9122, 64%
#(mean count < 13)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results#
#

#out of 14258 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 398, 2.8%
#LFC < 0 (down)     : 222, 1.6%
#outliers [1]       : 0, 0%
#low counts [2]     : 2488, 17%
#(mean count < 2)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

# > colnames(resRaw)
#  [1] "YV_vs_Y_pvalue"          "YV_vs_Y_padj"           
#  [3] "YV_vs_Y_log2FoldChange"  "OV_vs_O_pvalue"         
#  [5] "OV_vs_O_padj"            "OV_vs_O_log2FoldChange" 
#  [7] "O_vs_Y_pvalue"           "O_vs_Y_padj"            
#  [9] "O_vs_Y_log2FoldChange"   "OV_vs_YV_pvalue"        
# [11] "OV_vs_YV_padj"           "OV_vs_YV_log2FoldChange"


# 读取注释文件
geneanno <-  readxl::read_xlsx(path = "/data/home/mali/genome/mouse/gene.xlsx",sheet="gene")
geneanno <- geneanno %>% drop_na(gene_id)
# 直接从注释文件中获得基因名字
geneanno <- geneanno[,c(1:9)]
colnames(geneanno)[1] <- 'ID'
colnames(geneanno)[2] <- 'Gene'
colnames(geneanno)[9] <- 'Description'
geneanno$ID <- as.character(geneanno$ID)
head(geneanno)
## 获得总表，既有注释，原始表达矩阵，也有差异表达logFC，p等信息
res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), counts(degres, normalize = TRUE), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as_tibble %>%
  bind_cols(resRaw) %>%  # 相当于cbind
  inner_join(geneanno, by = 'ID') %>%  # 结果可以理解为a、b的交集，R语言中的merge函数也可以实现
  dplyr::select(ID, Gene : Description, YV1 : OV_vs_YV_log2FoldChange) # %>%
  # arrange(Mock_Flg22_vs_Mock_padj)
  # ‘arrange()’ orders the rows of a data frame by the values of selected columns.

#> head(res)
## A tibble: 6 × 34
#  ID        Gene  gene_…¹ gene_…² gene_…³ gene_…⁴ gene_…⁵ gene_…⁶ Descr…⁷    YV1
#  <chr>     <chr> <chr>     <dbl>   <dbl> <chr>     <dbl> <chr>   <chr>    <dbl>
#1 20389     Sftpc NC_000…  7.08e7  7.08e7 -           804 protei… surfac… 8762. 
#2 22287     Scgb… NC_000…  9.06e6  9.07e6 -           468 protei… secret… 2567. 
#3 100503605 Hbb-… NC_000…  1.03e8  1.03e8 -           639 protei… hemogl… 8422. 
#4 15122     Hba-… NC_000…  3.22e7  3.22e7 +           559 protei… hemogl… 6853. 
#5 16071     Igkc  NC_000…  7.07e7  7.07e7 +           320 C_regi… - && P…   91.8
#6 17105     Lyz2  NC_000…  1.17e8  1.17e8 -          1057 protei… lysozy… 2751. 
## … with 24 more variables: YV2 <dbl>, YV3 <dbl>, YV4 <dbl>, OV1 <dbl>,
##   OV2 <dbl>, OV3 <dbl>, O1 <dbl>, O2 <dbl>, O3 <dbl>, Y1 <dbl>, Y2 <dbl>,
##   Y3 <dbl>, YV_vs_Y_pvalue <dbl>, YV_vs_Y_padj <dbl>,
##   YV_vs_Y_log2FoldChange <dbl>, OV_vs_O_pvalue <dbl>, OV_vs_O_padj <dbl>,
##   OV_vs_O_log2FoldChange <dbl>, O_vs_Y_pvalue <dbl>, O_vs_Y_padj <dbl>,
##   O_vs_Y_log2FoldChange <dbl>, OV_vs_YV_pvalue <dbl>, OV_vs_YV_padj <dbl>,
##   OV_vs_YV_log2FoldChange <dbl>, and abbreviated variable names ¹​gene_chr, …
## ℹ Use `colnames()` to see all variable names
#> colnames(res)
# [1] "ID"                      "Gene"                   
# [3] "gene_chr"                "gene_start"             
# [5] "gene_end"                "gene_strand"            
# [7] "gene_length"             "gene_biotype"           
# [9] "Description"             "YV1"                    
#[11] "YV2"                     "YV3"                    
#[13] "YV4"                     "OV1"                    
#[15] "OV2"                     "OV3"                    
#[17] "O1"                      "O2"                     
#[19] "O3"                      "Y1"                     
#[21] "Y2"                      "Y3"                     
#[23] "YV_vs_Y_pvalue"          "YV_vs_Y_padj"           
#[25] "YV_vs_Y_log2FoldChange"  "OV_vs_O_pvalue"         
#[27] "OV_vs_O_padj"            "OV_vs_O_log2FoldChange" 
#[29] "O_vs_Y_pvalue"           "O_vs_Y_padj"            
#[31] "O_vs_Y_log2FoldChange"   "OV_vs_YV_pvalue"        
#[33] "OV_vs_YV_padj"           "OV_vs_YV_log2FoldChange"

write_csv(res, "/data/home/mali/ke_data/transcription/count/eachGroup_vs_Mock_k.csv")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('ggrepel')
library('ggplot2')
library('RColorBrewer')
library('limma')
library('sva')
library('reshape2') # melt
library("factoextra") # get_pca_var
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## remove low count
dat <- rld %>%
  assay %>%
  {.[rowMeans(.) > 1, ]}

group <- sampleTable$condition
design <- model.matrix(~ group)
rldData <- dat %>%
  removeBatchEffect(covariates = svobj$sv,
                    design = design)
# removeBatchEffect属于limma包
# removeBatchEffect或ComBat。但是处理之后会改变counts矩阵
# 这个函数用于进行聚类或无监督分析之前，移除与杂交时间或其他技术变异相关的批次效应。
# 它是针对芯片设计的，因此不要直接使用read counts，数据需要经过一定的标准化操作，如log转化。
# removeBatchEffect只用于衔接聚类、PCA等可视化展示，不要在线性建模之前使用。
# 因为用矫正后的数据进行差异表达分析，有两个缺陷：
#  一是批次因素和分组因素可能重叠，所以直接对原数据矫正批次可能会抵消一部分真实生物学因素；
#  二是低估了误差。
# 所以，如果想做差异表达分析，但数据中又有已知的批次问题，则最好把批次效应纳入线性模型中。


# 去除批次效应前后，两者由细微差别
# > head(dat,2)
#               Mock_1   Mock_2   Mock_3 Mock_Flg22_1 Mock_Flg22_2 Mock_Flg22_3
# AT1G18320.1 4.038552 3.919625 4.071504     3.828454      3.86371     3.894067
# AT5G11100.1 8.253179 8.001945 7.935380     8.443196      8.60240     8.436905
#             SynCom33_Flg22_1 SynCom33_Flg22_2 SynCom33_Flg22_3 SynCom35_Flg22_1
# AT1G18320.1         3.953679         3.875093         3.911458         4.016012
# AT5G11100.1         7.796617         7.881670         7.937094         7.728326
#             SynCom35_Flg22_2 SynCom35_Flg22_3
# AT1G18320.1         3.874353         3.860681
# AT5G11100.1         7.635167         7.612592
# > head(rldData,2)
#               Mock_1   Mock_2   Mock_3 Mock_Flg22_1 Mock_Flg22_2 Mock_Flg22_3
# AT1G18320.1 4.011298 3.973795 4.063460     3.839003     3.914161     3.828008
# AT5G11100.1 8.137680 8.068464 7.985115     8.365431     8.573095     8.522416
#             SynCom33_Flg22_1 SynCom33_Flg22_2 SynCom33_Flg22_3 SynCom35_Flg22_1
# AT1G18320.1         3.894719         3.891920         3.894698         3.999920
# AT5G11100.1         7.757249         7.841183         7.882388         7.797953
#             SynCom35_Flg22_2 SynCom35_Flg22_3
# AT1G18320.1         3.926128         3.870080
# AT5G11100.1         7.691977         7.641521
dim(rldData)
# [1] 12150    13

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# cols <- colData(rldData)[, 1] %>% factor(levels = c('#a6cee3', '#1f78b4', '#e31a1c', '#6a3d9a')) # rldData是一个matrix，不能用colData方法；colData(rld)可以
# 目的是获取列名，修改成下面这个
#cols <- colData(rld)[, 1]  %>% factor(levels = c('#a6cee3', '#1f78b4', '#e31a1c', '#6a3d9a'))
# 还是不对，但是下面也没用呀 (T _ T)
#> cols
# [1] <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA>
#Levels: #a6cee3 #1f78b4 #e31a1c #6a3d9a

#‘colData(x)’, ‘colData(x) <- value’: Get or set the column data.
#          ‘value’ is a DataFrame object. Row names of ‘value’ must be
#          NULL or consistent with the existing column names of ‘x’.
#SummarizedExperiment objects

#Description:
#     The SummarizedExperiment class is a matrix-like container where
#     rows represent features of interest (e.g. genes, transcripts,
#     exons, etc...)  and columns represent samples (with sample data
#     summarized as a DataFrame). A SummarizedExperiment object contains
#     one or more assays, each represented by a matrix-like object of
#     numeric or other mode.
#     
useVersion   <- "version1"
outputFolder <- file.path("/data/home/mali/ke_data/output", useVersion, "DESeq_RNA")
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)

RNA_eSet_Scale <- readRDS("/data/home/mali/ke_data/transcription/RDS/RNA_rld_Scale.RDS")
colnames(RNA_eSet_Scale) <- RNA_eSet_Scale$name 
clcString <- "median_expression_"
RNA_eSet_Scale <- addMedianExpressions(eSet = RNA_eSet_Scale,
                                 use = "group",
                                 removeReplicates = FALSE,
                                 recognitionString = clcString)

# 对应文章中的 Extended Data Fig. 6a
pca <- prcomp(t(rldData))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rld)[, 1], ID = rownames(colData(rld))) %>%
  mutate(age = rep(c('Young', 'Old', 'Young'), c(4, 6, 3))) %>%
  mutate(infection = c('with','without') %>% rep(c(7, 6)) %>% factor)

ggplot(pcaData, aes(x = PC1, y = PC2, colour = age)) +
  geom_point(aes(shape = infection), size = 4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_colour_manual(values = c('#000000', '#377eb8', '#e41a1c')) +
  scale_shape_manual(values = c(17, 1, 16),
                     name = 'Experimental\nConditions') +
  stat_ellipse(aes(x = PC1, y = PC2, group = Group, linetype = infection), type = 't', level = 0.7) +
  scale_linetype_manual(values = c(1, 2), guide = FALSE) +
  coord_fixed(1) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.text.align = 0,
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        legend.text=element_text(size= 13),
        legend.title = element_text(size = 14))

# 用的是rldData，用removeBatchEffect处理的批次效应，不应该叫PCA_sva
ggsave(file.path(outputFolder,'PCA_sva.pdf'), width = 13)
ggsave(file.path(outputFolder,'PCA_sva.jpg'), width = 13)




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
# 取前100个贡献基因
pca1_name <- names(pca1[1:50])
pca2_name <- names(pca2[1:50])
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
#> head(gene_mat_pca1_mean)
#                sample      value
#1 median_expression_YV -0.8844326
#2 median_expression_YV -1.3752713
#3 median_expression_YV  0.2433362
#4 median_expression_YV -0.9020273
#5 median_expression_YV -1.3393768
#6 median_expression_YV  0.8093122
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


pdf(file.path(outputFolder,"pca_line.pdf"),width=3,height=3)
p1 <- ggplot()+
  geom_line(data = gene_mat_pca1_mean,aes(x=sample,y=value,group = gene),color ="grey60",size=0.5)+
  geom_line(data = gene_mat_pca1_mean_sample, aes(x=sample,y=value,group = 1),color = "black",size=1)+
  xlab("PC1")+#横坐标名称
  ylab("normalised mean expression")+#纵坐标名称
  theme_bw() +#去掉背景灰色
   #ylim(0,15000) +
   #NoLegend() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA),strip.background = element_rect(fill = NA,colour = NA), axis.ticks.x = element_blank(), axis.text = element_text( size = 10), axis.title = element_text(size = 10))
p2 <- ggplot()+
  geom_line(data = gene_mat_pca2_mean,aes(x=sample,y=value,group = gene),color ="grey60",size=0.5)+
  geom_line(data = gene_mat_pca2_mean_sample, aes(x=sample,y=value,group = 1),color = "black",size=1)+
  xlab("PC2")+#横坐标名称
  ylab("normalised mean expression")+#纵坐标名称
  theme_bw() +#去掉背景灰色
   #ylim(0,15000) +
   #NoLegend() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA),strip.background = element_rect(fill = NA,colour = NA), axis.ticks.x = element_blank(), axis.text = element_text( size = 10), axis.title = element_text(size = 10))

p1
p2
dev.off()


### 基因及富集分析
library(msigdbr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(enrichplot)
library(GSEABase)

type <- "RNA"
LOGFC_col <- "TOP50"
# 导出影响pca的前50个基因
pca1_name_50 <- NULL
pca1_name_50$ENTREZID <- names(pca1[1:50] )
#pca1_name_50$ENTREZID <- names(pca1)
pca1_name_50$exp <- pca1[1:50] 
#pca1_name_50$exp <- pca1
pca1_name_50 <- as.data.frame(pca1_name_50)


pca2_name_50 <- NULL
pca2_name_50$ENTREZID <- names(pca2[1:50] )
#pca2_name_50$ENTREZID <- names(pca2)
pca2_name_50$exp <- pca2[1:50] 
#pca2_name_50$exp <- pca2
pca2_name_50 <- data.frame(pca2_name_50)
pca2_name_50 <- as.data.frame(pca2_name_50)


write.table(pca1_name_50, file = file.path(outputFolder,paste0(type,"_",LOGFC_col,"_pca1.txt")), row.names = F,col.name = F,quote = F,sep = "\t")
write.table(pca2_name_50, file = file.path(outputFolder,paste0(type,"_",LOGFC_col,"_pca2.txt")), row.names = F,col.name = F,quote = F,sep = "\t")


# 自定义基因集读取
geneset_immune <- read.gmt("/data/home/mali/biosoft/MSigDB/geneset_immune.gmt")

#开始ID转换， 制作能输入到GSEA分析函数中的 genelist（pca1，pca2原本是SYMBOL作为名字，制作 ENTREZID作为名字）
pca1_name_change <- bitr(pca1_name_50$ENTREZID, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Mm.eg.db")
pca1_name_50_tmp <- left_join(pca1_name_50,pca1_name_change, by="ENTREZID")
pca1_name_50_exp <- pca1_name_50_tmp$exp
names(pca1_name_50_exp) <- pca1_name_50_tmp$SYMBOL
pca1_name_50_exp <- data.frame(pca1_name_50_exp)
colnames(pca1_name_50_exp) <- "exp"
head(pca1_name_50_exp)

#开始ID转换， 制作能输入到GSEA分析函数中的 genelist（pca1，pca2原本是SYMBOL作为名字，制作 ENTREZID作为名字）
pca2_name_change <- bitr(pca2_name_50$ENTREZID, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Mm.eg.db")
pca2_name_50_tmp <- left_join(pca2_name_50,pca2_name_change, by="ENTREZID")
pca2_name_50_exp <- pca2_name_50_tmp$exp
names(pca2_name_50_exp) <- pca2_name_50_tmp$SYMBOL
pca2_name_50_exp <- data.frame(pca2_name_50_exp)
colnames(pca2_name_50_exp) <- "exp"
head(pca2_name_50_exp)

write.table(pca1_name_50_exp, file = file.path(outputFolder,paste0(type,"_",LOGFC_col,"_pca1_symbol.txt")), row.names = T,col.name = F,quote = F,sep = "\t")
write.table(pca2_name_50_exp, file = file.path(outputFolder,paste0(type,"_",LOGFC_col,"_pca2_symbol.txt")), row.names = T,col.name = F,quote = F,sep = "\t")


  pca1_gsea_h <- GSEA(pca1_name_50_exp[,"exp"], TERM2GENE=geneset_immune, verbose=FALSE,minGSSize = 2, maxGSSize = 200,pvalueCutoff = 0.5)
  pca1_gsea_results <- pca1_gsea_h@result
  pca1_gsea_results<-pca1_gsea_results[order(pca1_gsea_results$qvalues, decreasing = F),]
  write.csv(pca1_gsea_results,file = file.path(outputFolder,paste0(type,"_",LOGFC_col,"_pca_gsea.csv")))
  pca1_plot<- pca1_gsea_results[1:10,]

  pca2_gsea_h <- GSEA(pca2_name_50_exp[,"exp"], TERM2GENE=geneset_immune, verbose=FALSE,minGSSize = 2, maxGSSize = 200,pvalueCutoff = 0.5)
  pca2_gsea_results <- pca2_gsea_h@result
  pca2_gsea_results<-pca2_gsea_results[order(pca2_gsea_results$qvalues, decreasing = F),]
  write.csv(pca2_gsea_results,file = file.path(outputFolder,paste0(type,"_",LOGFC_col,"_pca_gsea.csv")))
  pca2_plot<- pca2_gsea_results[1:10,]
  pdf(file.path(outputFolder,paste0(type,"_",LOGFC_col,"_pca_gsea.pdf")),width=5,height=4)
  p1 <- ggplot(pca1_plot, aes(x = ID)) +
   geom_col(aes(y=-log10(qvalues)), colour = "black", width = 0.6)+
   geom_line(aes(y=NES/4, group =1), colour = "grey60") +
   scale_y_continuous(sec.axis = sec_axis(~.*4, name = "Normalized Enrichment Score")) +
   scale_colour_manual(values = c("black", "grey60")) +
   labs(y = "-log10(q-value)",x= "pca1",colour = "black")+
   #coord_flip()+ 
   #NoLegend() +
   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA),strip.background = element_rect(fill = NA,colour = NA), axis.text = element_text( size = 10), axis.title = element_text(size = 10),axis.text.x = element_text(angle = 90, hjust = 1),axis.line = element_line(colour = "black"))
  p2 <- ggplot(pca2_plot, aes(x = ID)) +
   geom_col(aes(y=-log10(qvalues)), colour = "black", width = 0.6)+
   geom_line(aes(y=NES/4, group =1), colour = "grey60") +
   scale_y_continuous(sec.axis = sec_axis(~.*4, name = "Normalized Enrichment Score")) +
   scale_colour_manual(values = c("black", "grey60")) +
   labs(y = "-log10(q-value)",x= "pca2",colour = "black")+
   #coord_flip()+ 
   #NoLegend() +
   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA),strip.background = element_rect(fill = NA,colour = NA), axis.text = element_text( size = 10), axis.title = element_text(size = 10),axis.text.x = element_text(angle = 90, hjust = 1),axis.line = element_line(colour = "black"))
  
  p1
  p2
  dev.off()


# 保存用sva去除批次效应的DESeq文件，以及用rlog变换后removeBatchEffect后的表达矩阵；
save(degres, rldData, file = '/data/home/mali/ke_data/transcription/RDS/degres_condi_Mock.RData')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################
