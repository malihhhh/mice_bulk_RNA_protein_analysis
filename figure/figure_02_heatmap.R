
## originally by Yulong Niu
## yulong.niu@hotmail.com

###########################replot heatmap##############################
library('tidyverse')
library('ComplexHeatmap')
library('limma')
library('DESeq2')
library('RColorBrewer')
library('magrittr')
library('xlsx')
library('biomaRt')
library('dplyr')
setwd('/data/home/mali/ke_data/')

useVersion   <- "version1"
outputFolder <- file.path("/data/home/mali/ke_data/output", useVersion, "replot_heatmap")
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)

`%!in%` <- compose(`!`, `%in%`)

load('/data/home/mali/ke_data/transcription/RDS/degres_condi_Mock.RData')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~select DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wholeDEG <- read_csv('./transcription/count/eachGroup_vs_Mock_k.csv')
dim(wholeDEG)
kmeansRes <- read_csv('./transcription/count/kmeans10.csv') %>%
  dplyr::select(ID, cl)
head(kmeansRes)
dim(kmeansRes)

# 加载gene注释文件
anno <- read_csv('/data/home/mali/genome/mouse/Ensembl_mmu_Anno_AddGeneName.csv',
                    col_types = cols(Chromosome = col_character())) %>%
  mutate_all(list(~replace(., is.na(.), ''))) %>%
  mutate(GeneID = entrezgene_id) %>%  # 因为要跟下面的数据对应上，下面用的是entrezgene_id
  dplyr::select(GeneID, Gene, Description) %>%
  dplyr::slice(which(!duplicated(.)))

fcsig <- wholeDEG %>%
  dplyr::select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > log2(2) ~ 1,
                                 . < -log2(2) ~ -1,
                                 TRUE ~ 0)))
dim(fcsig)


# > head(fcsig)
# # A tibble: 6 × 4
#   YV_vs_Y_log2FoldChange OV_vs_O_log2FoldChange O_vs_Y_log2FoldChange OV_vs_YV…¹
#                    <dbl>                  <dbl>                 <dbl>      <dbl>
# 1                     -1                     -1                     0          0
# 2                     -1                      0                     1          1
# 3                     -1                     -1                    -1         -1
# 4                     -1                     -1                    -1         -1
# 5                     -1                      0                     1          1
# 6                     -1                      0                     0          1
# # … with abbreviated variable name ¹​OV_vs_YV_log2FoldChange

padjsig <- wholeDEG %>%
  dplyr::select(ends_with('padj')) %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))
# > head(padjsig)
# # A tibble: 6 × 4
#   YV_vs_Y_padj OV_vs_O_padj O_vs_Y_padj OV_vs_YV_padj
#   <lgl>        <lgl>        <lgl>       <lgl>
# 1 TRUE         TRUE         FALSE       FALSE
# 2 TRUE         FALSE        FALSE       TRUE
# 3 FALSE        FALSE        FALSE       TRUE
# 4 FALSE        FALSE        FALSE       TRUE
# 5 TRUE         FALSE        TRUE        TRUE
# 6 TRUE         FALSE        FALSE       TRUE

heatsig <- (padjsig * fcsig) %>%
  as_tibble %>%
  abs %>%
  rowSums %>%
  {. >= 1} %>%
  which %>%
  dplyr::slice(wholeDEG, .) %>%
  inner_join(kmeansRes)
dim(heatsig)

# 保存显著的差异基因
heatsig %>%
  mutate_at(c('Gene', 'Description'), .funs = list(~if_else(is.na(.), '', .))) %>%
  write_csv('/data/home/mali/ke_data/transcription/count/kmeans10_sig_2.csv')


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~heatmap~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## rlog transformed
heatsig$ID <- as.character(heatsig$ID)
rawC <- rldData %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  inner_join(heatsig %>% dplyr::select(ID, cl))
  ## inner_join(kmeansRes) ## all transcripts

# > head(rawC)
# # A tibble: 6 × 15
#   ID       YV1   YV2   YV3   YV4   OV1   OV2   OV3    O1    O2    O3    Y1    Y2
#   <chr>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 20389  13.2  13.0   13.1  13.1  13.1  13.2  13.0  13.5  13.6  13.6  13.8  13.8
# 2 22287  12.4  12.3   12.3  12.4  13.5  13.3  13.4  13.5  13.5  13.5  12.9  13.3
# 3 10050… 12.8  12.5   12.6  12.5  11.7  12.3  11.7  12.4  12.1  12.2  13.1  13.0
# 4 15122  12.5  12.3   12.4  12.2  11.4  12.0  11.4  12.2  11.9  12.0  13.0  12.8
# 5 16071   9.84  9.62  10.1  10.1  11.5  12.4  12.1  11.9  12.4  12.1  10.2  10.1
# 6 17105  11.6  11.4   11.5  11.5  11.8  11.9  11.8  12.0  12.1  12.1  11.9  11.9
# # … with 2 more variables: Y3 <dbl>, cl <dbl>
# # ℹ Use `colnames()` to see all variable names

## scale counts
scaleC <- rawC %>%
  dplyr::select(contains(c('O','Y'))) %>% # match('O|Y') 或者 contains(c('O','Y'))
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% dplyr::select(ID, cl)) %>%
  mutate(GeneID = ID) %>% inner_join(anno)

# > head(scaleC)
# # A tibble: 6 × 15
#       OV1    OV2     OV3     O1     O2     O3    YV1    YV2    YV3     YV4
#     <dbl>  <dbl>   <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl>
# 1 -1.01   -0.466 -1.06    0.474  0.729  0.757 -0.408 -1.18  -0.753 -1.02
# 2  0.925   0.500  0.714   0.908  0.818  0.976 -1.25  -1.45  -1.47  -1.30
# 3 -1.67   -0.364 -1.64   -0.119 -0.756 -0.485  0.639  0.174  0.368 -0.0147
# 4 -1.72   -0.524 -1.62   -0.124 -0.564 -0.387  0.491  0.126  0.360 -0.0161
# 5  0.505   1.29   1.05    0.871  1.33   1.00  -1.02  -1.22  -0.789 -0.776
# 6 -0.0823  0.481 -0.0646  0.883  1.27   1.18  -0.774 -1.73  -1.30  -1.37
# # … with 5 more variables: Y1 <dbl>, Y2 <dbl>, Y3 <dbl>, ID <chr>, cl <dbl>
# # ℹ Use `colnames()` to see all variable names


## 1. DEG annotation
scaleC$ID <- as.double(scaleC$ID)
sigMat <- (padjsig * fcsig) %>%
  as_tibble %>%
  setNames(names(.) %>% substr(., start = 1, stop = nchar(.) - 5)) %>%
  dplyr::mutate(ID = wholeDEG$ID) %>%
  inner_join(scaleC %>% dplyr::select(ID), .) %T>%
  {(sum(.$ID == scaleC$ID) == nrow(.)) %>% print} %>%
  transmute_at(.var = vars(contains('vs')),
               list(~ case_when(. == -1 ~ 'down',
                                . == 0 ~'no',
                                . == 1 ~ 'up'))) %>%
  as.matrix


#> head(sigMat)
#     YV_vs_Y OV_vs_O O_vs_Y OV_vs_YV
#[1,] "down"  "down"  "no"   "no"
#[2,] "down"  "no"    "no"   "up"
#[3,] "no"    "no"    "no"   "down"
#[4,] "no"    "no"    "no"   "down"
#[5,] "down"  "no"    "up"   "up"
#[6,] "down"  "no"    "no"   "up"

colnames(sigMat) <- c('YV vs.Y',
                      'OV vs. O',
                      'O vs. Y',
                      'OV vs. YV')

# 这里是人的基因集
ISG.list <- read.csv("/data/home/mali/bd_scRNAseq_PBMC_SARSCOV2_patient_Guan/FuXian/scripts/isg.list.csv")
inf.list <- read.table("/data/home/mali/bd_scRNAseq_PBMC_SARSCOV2_patient_Guan/mild_cmild/geneset/geneset_inflammation.txt",header= F,skip=2)
TNF.list <- read.table("/data/home/mali/bd_scRNAseq_PBMC_SARSCOV2_patient_Guan/mild_cmild/geneset/geneset_TNF_NFKB.txt",header= F,skip=2)
INFA.list <- read.table("/data/home/mali/bd_scRNAseq_PBMC_SARSCOV2_patient_Guan/mild_cmild/geneset/HALLMARK_INTERFERON_ALPHA_RESPONSE.txt",header= F,skip=2)
INFG.list <- read.table("/data/home/mali/bd_scRNAseq_PBMC_SARSCOV2_patient_Guan/mild_cmild/geneset/HALLMARK_INTERFERON_GAMMA_RESPONSE.txt",header= F,skip=2)
MHC_GENE <- read.table("/data/home/mali/genome/mouse/geneset/MHC_GENE_MOUSE.txt",header = F,skip=1)
Complement_and_coagulation <- read.table("/data/home/mali/genome/mouse/geneset/Complement_and_coagulation.txt",header = F,skip=1)

# 变换成小鼠的
mart = useMart('ensembl')
human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl", verbose = TRUE, host = "dec2021.archive.ensembl.org")
mouse <- useMart('ensembl',dataset = "mmusculus_gene_ensembl", verbose = TRUE, host = "dec2021.archive.ensembl.org")


############################################################# 修改
# 在这里选定要看的基因集

sigMat_geneset <- c()
#############################################################
mark_gene <- c(ISG.list[[1]])
name <- "ISG"
hsa2mus_all <- getLDS(attributes = c("hgnc_symbol"),
                filters = "hgnc_symbol",
                values = mark_gene,
                mart = human,
                attributesL = c("mgi_symbol"),
                martL = mouse,uniqueRows = T)
sigMat_geneset_tmp <- scaleC %>%
                  dplyr::select(Gene) %>%
                  transmute_at(.var = vars(c('Gene')),
                               list(~ case_when(. %in% hsa2mus_all$MGI.symbol ~ 'Contain',
                                                . %!in% hsa2mus_all$MGI.symbol ~ 'no'))) %>%
                                                as.matrix
sigMat_geneset$ISG <- sigMat_geneset_tmp

mark_gene <- c(INFA.list[[1]])
name <- "INFA"
hsa2mus_all <- getLDS(attributes = c("hgnc_symbol"),
                filters = "hgnc_symbol",
                values = mark_gene,
                mart = human,
                attributesL = c("mgi_symbol"),
                martL = mouse,uniqueRows = T)
sigMat_geneset_tmp <- scaleC %>%
                  dplyr::select(Gene) %>%
                  transmute_at(.var = vars(c('Gene')),
                               list(~ case_when(. %in% hsa2mus_all$MGI.symbol ~ 'Contain',
                                                . %!in% hsa2mus_all$MGI.symbol ~ 'no'))) %>%
                                                as.matrix
sigMat_geneset$INFA <- sigMat_geneset_tmp


mark_gene <- c(INFG.list[[1]])
name <- "INFG"
hsa2mus_all <- getLDS(attributes = c("hgnc_symbol"),
                filters = "hgnc_symbol",
                values = mark_gene,
                mart = human,
                attributesL = c("mgi_symbol"),
                martL = mouse,uniqueRows = T)
sigMat_geneset_tmp <- scaleC %>%
                  dplyr::select(Gene) %>%
                  transmute_at(.var = vars(c('Gene')),
                               list(~ case_when(. %in% hsa2mus_all$MGI.symbol ~ 'Contain',
                                                . %!in% hsa2mus_all$MGI.symbol ~ 'no'))) %>%
                                                as.matrix
sigMat_geneset$INFG <- sigMat_geneset_tmp

mark_gene <- MHC_GENE$V1
name <- "MHC"
sigMat_geneset_tmp <- scaleC %>%
                  dplyr::select(Gene) %>%
                  transmute_at(.var = vars(c('Gene')),
                               list(~ case_when(. %in% mark_gene ~ 'Contain',
                                                . %!in% mark_gene ~ 'no'))) %>%
                                                as.matrix
sigMat_geneset$MHC <- sigMat_geneset_tmp

#mark_gene <- Complement_and_coagulation$V1
#name <- "Complement"
#sigMat_geneset_tmp <- scaleC %>%
#                  dplyr::select(Gene) %>%
#                  transmute_at(.var = vars(c('Gene')),
#                               list(~ case_when(. %in% hsa2mus_all$mark_gene ~ 'Contain',
#                                                . %!in% hsa2mus_all$mark_gene ~ 'no'))) %>%
#                                                as.matrix
#sigMat_geneset$Complement <- sigMat_geneset_tmp

sigMat_geneset <- as.data.frame(sigMat_geneset)
colnames(sigMat_geneset) <- c("ISG","INFA","INFG","MHC")
###-------------------------------绘图
ht_list <- Heatmap(matrix = scaleC %>% dplyr::select(contains(c('O','Y'),ignore.case=F)),
                   name = 'Scaled Counts',
                   ## row_order = order(scaleC$cl) %>% rev,
                   row_split = scaleC$cl,
                   row_gap = unit(2, "mm"),
                   column_order = 1 : 13, # Order of column. # 改
                   column_split = rep(c('OV', 'O', 'YV', 'Y'), c(3,3,4,3)), # 改
                   show_column_names = FALSE,
                   col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -6, -7)])(100),
                   column_title_gp = gpar(fontsize = 7),
                   use_raster = FALSE) +
  Heatmap(sigMat,
          col = c('down' = 'blue', 'no' = 'white', 'up' = 'red'),
          column_names_gp = gpar(fontsize = 5),
          heatmap_legend_param = list(title = 'DEGs'),
          cluster_columns = FALSE,
          use_raster = FALSE) +
  Heatmap(sigMat_geneset,
          col = c('Contain' = 'black', 'no' = 'white'),
          column_names_gp = gpar(fontsize = 5),
          heatmap_legend_param = list(title = 'Select geneset'),
          cluster_columns = FALSE,
          use_raster = FALSE)

filePrefix <- 'kmeans10_heatmap_sig_DEG2_geneset_all'

pdf(file.path(outputFolder,paste0(filePrefix, ".pdf")))
draw(ht_list)
dev.off()
# 和redult中的结果顺序不一样，大致结果是一样的；

# 保存source data
data_all <- cbind(scaleC,sigMat,sigMat_geneset %>% as_tibble)
write.xlsx(x = data_all, file = file.path(outputFolder,"FIG_heatmap_source_data_geneset_all.xlsx"),  # 确定好图片位置后修改
        sheetName = "gene_set", row.names = FALSE) # 确定好图片位置后修改

# Linux之convert命令 强大的convert命令 convert命令可以用来转换图像的格式，支持JPG, BMP, PCX, GIF, PNG, TIFF, XPM和XWD等类型，还可以改变文件大小
# 没装这个，有pdf文件就够了
#system(paste0('convert -density 1200 ', paste0(filePrefix, '.pdf'), ' ', paste0(filePrefix, '.jpg')))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~box plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## kmeansRes <- read_csv('kmeans10.csv') %>%
##   select(ID, cl)

kmeansResSig <- read_csv('/data/home/mali/ke_data/transcription/count/kmeans10_sig.csv') %>%
  dplyr::select(ID, cl)
kmeansResSig$ID <- as.character(kmeansResSig$ID)


colnames(rldData) # 注意根据列名调整meanFlg22函数的分组,以及sampleN
meanFlg22 <- function(v) {

  require('magrittr')

  res <- v %>%
    split(c(rep(1, each = 4),rep(2:4, each = 3))) %>% # 改
    sapply(mean, na.rm = TRUE)

  return(res)
}
sampleN <- c('YV', 'OV', 'O', 'Y')

boxplotData <- rldData %>%
  t %>%
  scale %>%
  t %>%
  apply(1, meanFlg22) %>%
  t %>%
  set_colnames(sampleN) %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  inner_join(kmeansResSig)

# > head(boxplotData)
# # A tibble: 6 × 6
#   ID            OV       O      YV      Y    cl
#   <chr>      <dbl>   <dbl>   <dbl>  <dbl> <dbl>
# 1 20389     -0.779 -0.832   0.0486  1.17      8
# 2 22287     -1.39   0.0428  0.813   0.403     5
# 3 100503605  0.394 -0.681  -0.839   0.845     6
# 4 15122      0.326 -0.752  -0.768   0.896     6
# 5 16071     -1.01   0.340   1.08   -0.311     1
# 6 17105     -1.27  -0.325   0.695   0.675     4

comparisons <- list(c('YV', 'Y'),
             c('OV', 'O'),
             c('O', 'Y'),
             c('OV', 'YV'))

for (i in 1:length(unique(boxplotData$cl))) {
  boxplotData %>%
    filter(cl == i) %>%
    select(-ID, -cl) %>%
    gather(key = 'Conditions', value = 'ScaleCounts') %>%
    mutate(Group = case_when(
             str_detect(Conditions, 'YV') ~ 'Y',
             str_detect(Conditions, 'OV') ~ 'O',
             str_detect(Conditions, 'O') ~ 'Y',
             str_detect(Conditions, 'OV') ~ 'YV', # 改
           )) %>%
    mutate(Conditions = Conditions %>% factor(levels = sampleN)) %>%
    mutate(Group = Group %>% factor(levels = c('Y', 'YV', 'O', 'OV'))) %>%
    ggplot(aes(x = Conditions, y = ScaleCounts, fill = Conditions)) +
    ## geom_boxplot(position = position_dodge2(preserve = 'single')) +
    geom_boxplot(width = 0.45) +
    geom_jitter(width=0.2, alpha=0.1, size = 1) +
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF" ,"#3C5488FF")) + # 改颜色
    ggpubr::stat_compare_means( comparisons = comparisons, label = "p.signif") +
    ylim(-2, 2) +
    ylab('Scaled counts') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0),
          plot.title = element_text(hjust = 0.5, size = 8, face = 'bold'),
          legend.text.align = 0,
          axis.text = element_text(size = 6),
          axis.title = element_text(size = 7),
          legend.text=element_text(size= 6),
          legend.title = element_text(size = 7))

  ggsave(file.path(outputFolder,paste0('kmeans10_boxplot', i, '.pdf')),width = 3,height=3)
  ggsave(file.path(outputFolder,paste0('kmeans10_boxplot', i, '.jpeg')),width = 3,height=3)
}
# 和result已有的结果一摸一样

#######################################################################
