library(data.table)
library(readxl)
library(tidyverse)
library(matrixStats)
library(Biobase)
library(biomaRt)
library(tidyr)
library('readr')
library('magrittr')
library('tibble')
library('gplots')
library('dendextend')
library('dynamicTreeCut')
library('ggplot2')
library('tidyr')
library('DESeq2')
library('dplyr')
library('RColorBrewer')
library('gridExtra')
library('cluster')
library('scales')


sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")






####################################################################################
message("Make protein ExpressionSets")
# 蛋白组
####################################################################################
dataDirectory <- "/data/home/mali/ke_data/Protein"
outputFolder <- "/data/home/mali/ke_data/Protein/RDS"


expressionsFileName <- file.path(dataDirectory, "Annotation_combine.xlsx")
expressionsSheetName <- "Annotation_Combine"

# Assay data
 RawData <- read_excel(path = expressionsFileName,
                        sheet = expressionsSheetName)
exprs <- RawData[,5:16]
rownames(exprs) <- RawData$'Protein accession'

# Sample annotation  (phenoData)
pDataFile <- file.path(dataDirectory, "metadata_P.csv")
pData <- read.csv(pDataFile)
rownames(pData) <- pData$ID

# 需要注意，assay data和sample data的关系，样本应对应一致，否则Expression Set构建会出错
# order meta rows according to protein columns
quantCols <- colnames(exprs)
pData.ord <- pData[quantCols, ]
identical(rownames(pData.ord),colnames(exprs))

# (featureData)
featureData <- RawData[,-(5:16)]
rownames(featureData) <- RawData$'Protein accession'
# build up eSet
expressions.raw <- ExpressionSet(assayData = as.matrix(exprs),
                                   phenoData = AnnotatedDataFrame(data = pData.ord,
                                                                  varMetadata = data.frame(labelDescription = colnames(pData.ord))),
                                   featureData = AnnotatedDataFrame(data = featureData,
                                                                    varMetadata = data.frame(labelDescription = colnames(featureData))))


#######---------------------------------------------
message("Remove NA in the ExpressionSets")
dim(expressions.raw)
# 去除含有NA的行，函数放在functions文件夹里
expressions.rmNA <- removeNAsFromESet(expressions.raw)
dim(expressions.rmNA)
#######---------------------------------------------
protien_raw <- expressions.raw
protien <- expressions.rmNA

####################################################################################
message(name, "Store the proteins and PST_Protein objects as RDS")
####################################################################################

saveRDS(object = protien,
        file = file.path(outputFolder, "proteins.RDS"))
saveRDS(object = protien_raw,
        file = file.path(outputFolder, "protien_raw.RDS"))






####################################################################################
message("Store the proteins as TXT")
####################################################################################
proteins <- readRDS(file = file.path(outputFolder, "proteins.RDS"))
protien_raw <- readRDS(file = file.path(outputFolder, "protien_raw.RDS"))
saveESetAsTxt(eSet = protien_raw,
              file = file.path(outputFolder, "proteins_all.txt"),
              colNameColumn = "ID",
              rowNameColumn = "Protein accession",
              additionalRows = c("ID", "name", "group"), # phenoData中的列名
              additionalColumns = c("Protein accession","Gene name"), # featureData中的列名
              extendedColnames = FALSE)

saveESetAsTxt(eSet = removeNAsFromESet(protien_raw, na_ratio = 0),
              file = file.path(outputFolder, "proteins_full_overlap.txt"),
              colNameColumn = "ID",
              rowNameColumn = "Protein accession",
              additionalRows = c("ID", "name", "group"),
              additionalColumns = c("Protein accession","Gene name"),
              extendedColnames = FALSE)

saveESetAsTxt(eSet = protien_raw[, grepl("Y", pData(protien_raw)$group)] %>% removeNAsFromESet(na_ratio = 0),
              file = file.path(outputFolder, "proteins_Y_full_overlap.txt"),
              colNameColumn = "ID",
              rowNameColumn = "Protein accession",
              additionalRows = c("ID", "name", "group"),
              additionalColumns = c("Protein accession","Gene name"),
              extendedColnames = FALSE)

saveESetAsTxt(eSet = protien_raw[, grepl("O", pData(protien_raw)$group)] %>% removeNAsFromESet(na_ratio = 0),
              file = file.path(outputFolder, "proteins_O_full_overlap.txt"),
              colNameColumn = "ID",
              rowNameColumn = "Protein accession",
              additionalRows = c("ID", "name", "group"),
              additionalColumns = c("Protein accession","Gene name"),
              extendedColnames = FALSE)






####################################################################################
message("Make protein ExpressionSets")
# 磷酸化蛋白组  和上面一样
####################################################################################
dataDirectory <- "/data/home/mali/ke_data/PST_Protein"
outputFolder <- "/data/home/mali/ke_data/PST_Protein/RDS"

expressionsFileName <- file.path(dataDirectory, "Annotation_combine.xlsx")
expressionsSheetName <- "Annotation_Combine"

# Assay data
 RawData <- read_excel(path = expressionsFileName,
                        sheet = expressionsSheetName)
rownames(RawData) <- paste0(RawData$'Protein accession',"_",RawData$Position)
exprs <- RawData[,6:17]
rownames(exprs) <- rownames(RawData)

# Sample annotation
pDataFile <- file.path(dataDirectory, "metadata_PST.csv")
pData <- read.csv(pDataFile)
rownames(pData) <- pData$ID

# 需要注意，assay data和sample data的关系，样本应对应一致，否则Expression Set构建会出错
# order meta rows according to protein columns
quantCols <- colnames(exprs)
pData.ord <- pData[quantCols, ]
identical(rownames(pData.ord),colnames(exprs))

# featureData
featureData <- RawData[,-(6:17)]
rownames(featureData) <- rownames(RawData)
# build up eSet
  expressions.raw <- ExpressionSet(assayData = as.matrix(exprs),
                                   phenoData = AnnotatedDataFrame(data = pData.ord,
                                                                  varMetadata = data.frame(labelDescription = colnames(pData.ord))),
                                   featureData = AnnotatedDataFrame(data = featureData,
                                                                    varMetadata = data.frame(labelDescription = colnames(featureData))))


#######---------------------------------------------
message("Remove NA in the ExpressionSets")
dim(expressions.raw)
# 去除含有NA的行，函数放在functions文件夹里
expressions.rmNA <- removeNAsFromESet(expressions.raw)
dim(expressions.rmNA)
#######---------------------------------------------
PST_protien_raw <- expressions.raw
PST_protien <- expressions.rmNA

# featureData中增加一列列名
PST_protien_raw@featureData@data$Protein_accession_Position <- paste0(PST_protien_raw@featureData@data$'Protein accession',"_",PST_protien_raw@featureData@data$Position)
PST_protien@featureData@data$Protein_accession_Position <- paste0(PST_protien@featureData@data$'Protein accession',"_",PST_protien@featureData@data$Position)

####################################################################################
message(name, "Store the proteins and PST_Protein objects as RDS")
####################################################################################

saveRDS(object = PST_protien,
        file = file.path(outputFolder, "PST_proteins.RDS"))

saveRDS(object = PST_protien_raw,
        file = file.path(outputFolder, "PST_protien_raw.RDS"))


####################################################################################


message(name, "Store the proteins as TXT")
####################################################################################
outputFolder <- "/data/home/mali/ke_data/PST_Protein/RDS"
PST_protien_raw <- readRDS(file = file.path(outputFolder, "PST_protien_raw.RDS"))
PST_protien <- readRDS(file = file.path(outputFolder, "PST_proteins.RDS"))


saveESetAsTxt(eSet = PST_protien_raw,
              file = file.path(outputFolder, "proteins_all.txt"),
              colNameColumn = "ID",
              rowNameColumn = "Protein_accession_Position",
              additionalRows = c("ID", "name", "group"), # phenoData中的列名
              additionalColumns = c("Protein accession","Position","Amino acid","Gene name"),# featureData中的列名
              extendedColnames = FALSE)

saveESetAsTxt(eSet = PST_protien,
              file = file.path(outputFolder, "proteins_full_overlap.txt"),
              colNameColumn = "ID",
              rowNameColumn = "Protein_accession_Position",
              additionalRows = c("ID", "name", "group"),
              additionalColumns = c("Protein accession","Position","Amino acid","Gene name"),
              extendedColnames = FALSE)

saveESetAsTxt(eSet = PST_protien[, grepl("Y", pData(PST_protien)$group)],
              file = file.path(outputFolder, "proteins_Y_full_overlap.txt"),
              colNameColumn = "ID",
              rowNameColumn = "Protein_accession_Position",
              additionalRows = c("ID", "name", "group"),
              additionalColumns = c("Protein accession","Position","Amino acid","Gene name"),
              extendedColnames = FALSE)

saveESetAsTxt(eSet = PST_protien[, grepl("O", pData(PST_protien)$group)],
              file = file.path(outputFolder, "proteins_O_full_overlap.txt"),
              colNameColumn = "ID",
              rowNameColumn = "Protein_accession_Position",
              additionalRows = c("ID", "name", "group"),
              additionalColumns = c("Protein accession","Position","Amino acid","Gene name"),
              extendedColnames = FALSE)



####################################################################################
message( "Clear environment")
####################################################################################

rm(list = ls())





####################################################################################
message("Make RNA ExpressionSets")
# 转录组  和上面一样
####################################################################################

dataDirectory <- "/data/home/mali/ke_data/transcription/count"
outputFolder <- "/data/home/mali/ke_data/transcription/RDS"
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)

# 读取注释文件
geneanno <-  read_xlsx(path = "/data/home/mali/genome/mouse/gene.xlsx",sheet="gene")
geneanno <- geneanno %>% drop_na(gene_id)

# > head(geneanno)
# # A tibble: 6 × 95
#     gene_id gene_name    gene_…¹ gene_…² gene_…³ gene_…⁴ gene_…⁵ gene_…⁶ gene_…⁷
#       <dbl> <chr>        <chr>     <dbl>   <dbl> <chr>     <dbl> <chr>   <chr>  
# 1 115487594 Gm26206      NC_000… 3172239 3172348 +           110 snRNA   U6 spl…
# 2    497097 Xkr4         NC_000… 3269956 3741733 -         14824 protei… X-link…
# 3 118567655 LOC118567655 NC_000… 3363896 3448035 +         24933 lncRNA  unchar…
# 4 108167595 LOC108167595 NC_000… 3741733 3742701 -           969 pseudo… - && Q…
# 5 115487633 Gm27396      NC_000… 3854099 3854156 -            58 pseudo… - && -…
# 6     19888 Rp1          NC_000… 4185896 4479489 -         12639 protei… retini…
# # … with 86 more variables: tf_family <chr>, ...11 <chr>, ...12 <chr>,
# #   ...13 <chr>, ...14 <chr>, ...15 <chr>, ...16 <chr>, ...17 <chr>,
# #   ...18 <chr>, ...19 <chr>, ...20 <chr>, ...21 <chr>, ...22 <chr>,
# #   ...23 <chr>, ...24 <chr>, ...25 <chr>, ...26 <chr>, ...27 <chr>,
# #   ...28 <chr>, ...29 <chr>, ...30 <chr>, ...31 <chr>, ...32 <chr>,
# #   ...33 <chr>, ...34 <chr>, ...35 <chr>, ...36 <chr>, ...37 <chr>,
# #   ...38 <chr>, ...39 <chr>, ...40 <lgl>, ...41 <lgl>, ...42 <lgl>, …
# # ℹ Use `colnames()` to see all variable names
# > dim(geneanno)
# [1] 46331    95

# 读取基因表达矩阵
expressionsFileName <- file.path(dataDirectory, "2.gene_matrix_anno.csv")


# Assay data
RawData <- read_csv(expressionsFileName)
#> dim(RawData)
#[1] 46670    18

#———————————————————————————————————————————————— biomaRt会漏掉很多基因
## 添加gene symbol
## start biomaRt
#ensembl <- useEnsembl(biomart = "ensembl",
#                      dataset = "mmusculus_gene_ensembl")#

## get translations to gene names
## 该函数是主要的biomaRt查询函数。给定一组过滤器和相应的值，它从连接到的BioMart数据库中检索用户指定的属性。
#lookup.raw <- getBM(attributes = c("entrezgene_id","ensembl_gene_id", "external_gene_name","ensembl_transcript_id"),
#                    filters = "entrezgene_id",
#                    values = RawData$gene_id,
#                    mart = ensembl)
#write.csv(lookup.raw, file.path(dataDirectory, "lookup.raw.csv"),row.names =FALSE)#

## 有很多基因没有注释上
#dim(lookup.raw)
##[1] 27640     3  
#sum(is.na(lookup.raw$ensembl_gene_id))
#sum(is.na(lookup.raw$external_gene_name))
## [1] 0#

## make a lookup
#lookup <- lookup.raw[,c("external_gene_name","ensembl_gene_id")]
#row.names(lookup) <- lookup.raw$entrezgene_id
# get gene names in tpm table
#RawData$external_gene_name <- geneanno[RawData$gene_id,external_gene_name]
#RawData$ensembl_gene_id <- geneanno[RawData$gene_id,ensembl_gene_id]
#———————————————————————————————————————————————— 

# 直接从注释文件中获得基因名字
RawData_anno <- plyr::join(RawData, geneanno[,c(1:10)],by = "gene_id")

# set the rownames
row.names(RawData_anno) <- RawData_anno$gene_id

# 表达矩阵
exprs <- RawData_anno[,2:14] # 13个样本
rownames(exprs) <- rownames(RawData_anno)

# Sample annotation
pDataFile <- file.path(dataDirectory, "metadata_R.csv")
pData <- read.csv(pDataFile)
rownames(pData) <- pData$ID

# 需要注意，assay data和sample data的关系，样本应对应一致，否则Expression Set构建会出错
# order meta rows according to protein columns
quantCols <- colnames(exprs)
pData.ord <- pData[quantCols, ]
identical(rownames(pData.ord),colnames(exprs))

# featureData
dim(RawData_anno)
# > dim(RawData_anno)
# [1] 46670    27
featureData <- RawData_anno[,-(2:14)]
rownames(featureData) <- rownames(RawData_anno)
head(featureData)
# build up eSet
  expressions.raw <- ExpressionSet(assayData = as.matrix(exprs),
                                   phenoData = AnnotatedDataFrame(data = pData.ord,
                                                                  varMetadata = data.frame(labelDescription = colnames(pData.ord))),
                                   featureData = AnnotatedDataFrame(data = featureData,
                                                                    varMetadata = data.frame(labelDescription = colnames(featureData))))


#######---------------------------------------------
message("Remove NA in the ExpressionSets")
dim(expressions.raw)
# 去除含有NA的行，函数放在functions文件夹里
expressions.rmNA <- removeNAsFromESet(expressions.raw)
dim(expressions.rmNA)
# 两者一样
RNA  <- expressions.raw
#######---------------------------------------------

####################################################################################
message(name, "Store the proteins and PST_Protein objects as RDS")
####################################################################################

saveRDS(object = RawData_anno,
        file = file.path(outputFolder, "RawData_anno.RDS"))
saveRDS(object = expressions.raw,
        file = file.path(outputFolder, "RNA.RDS"))
####################################################################################

message(name, "Store the proteins as TXT")
####################################################################################
outputFolder <- "/data/home/mali/ke_data/transcription/RDS"
dataDirectory <- "/data/home/mali/ke_data/transcription/count"
RNA <- readRDS(file = file.path(outputFolder, "RNA.RDS"))


saveESetAsTxt(eSet = RNA,
              file = file.path(outputFolder, "RNA_raw_anno.txt"),
              colNameColumn = "ID",
              rowNameColumn = "gene_id",
              additionalRows = c("ID", "name", "group"), # phenoData中的列名
              additionalColumns = c("gene_id","gene_name","gene_biotype"),# featureData中的列名
              extendedColnames = FALSE)



####################################################################################
message( "Clear environment")
####################################################################################

rm(list = ls())



####################################################################################
message("Make RNA rld Data raw and scale ExpressionSets")
# 转录组scaledata 生成 eSet
####################################################################################

dataDirectory <- "/data/home/mali/ke_data/transcription/count"
outputFolder <- "/data/home/mali/ke_data/transcription/RDS"
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)

# 读取script_cluster.R 生成的表达矩阵 kmeans10.csv
kmeans10 <-  read_csv("/data/home/mali/ke_data/transcription/count/kmeans10.csv")

####################################################################################
## scale
# 读取基因表达矩阵
kmeans10_scale <- kmeans10 %>%
  dplyr::select(ID : Description, starts_with('Scale'))

expressionsFileName <- kmeans10_scale %>% 
  dplyr::select(ID , starts_with('Scale'))  %>% 
  column_to_rownames(var = 'ID')
#> dim(expressionsFileName)
#[1] 12088    13
#> sum(is.na(expressionsFileName))
#[1] 0
# Assay data
exprs <- expressionsFileName
colnames(exprs)
# phenoData
# Sample annotation
pDataFile <- file.path(dataDirectory, "metadata_R_rld_Scale.csv")
pData <- read.csv(pDataFile)
rownames(pData) <- pData$ID
head(pData)
pData.ord <- pData
# featureData
featureData <-  kmeans10_scale %>% 
  dplyr::select(ID : Description) %>% 
  column_to_rownames(var = 'ID') %>%
  mutate(ID = rownames(.)) %>%
  dplyr::select(ID,everything())
head(featureData)
# build up eSet
  expressions.rld.Scale <- ExpressionSet(assayData = as.matrix(exprs),
                                   phenoData = AnnotatedDataFrame(data = pData.ord,
                                                                  varMetadata = data.frame(labelDescription = colnames(pData.ord))),
                                   featureData = AnnotatedDataFrame(data = featureData,
                                                                    varMetadata = data.frame(labelDescription = colnames(featureData))))

head(expressions.rld.Scale)
#######---------------------------------------------
message("Remove NA in the ExpressionSets")
dim(expressions.rld.Scale)
# 去除含有NA的行，函数放在functions文件夹里
expressions.rld.Scale.rmNA <- removeNAsFromESet(expressions.rld.Scale)
dim(expressions.rld.Scale.rmNA)
# 两者一样
####################################################################################
## Raw
# 读取基因表达矩阵
kmeans10_Raw <- kmeans10 %>%
  dplyr::select(ID : Description, starts_with('Raw'))

expressionsFileName <- kmeans10_Raw %>% 
  dplyr::select(ID , starts_with('Raw'))  %>% 
  column_to_rownames(var = 'ID')
#> dim(expressionsFileName)
#[1] 12088    13
#> sum(is.na(expressionsFileName))
#[1] 0
# Assay data
exprs <- expressionsFileName
colnames(exprs)
# phenoData
# Sample annotation
pDataFile <- file.path(dataDirectory, "metadata_R_rld_Raw.csv")
pData <- read.csv(pDataFile)
rownames(pData) <- pData$ID
head(pData)
pData.ord <- pData
# featureData
featureData <-  kmeans10_Raw %>% 
  dplyr::select(ID : Description) %>% 
  column_to_rownames(var = 'ID') %>%
  mutate(ID = rownames(.)) %>%
  dplyr::select(ID,everything())
head(featureData)
# build up eSet
  expressions.rld.Raw <- ExpressionSet(assayData = as.matrix(exprs),
                                   phenoData = AnnotatedDataFrame(data = pData.ord,
                                                                  varMetadata = data.frame(labelDescription = colnames(pData.ord))),
                                   featureData = AnnotatedDataFrame(data = featureData,
                                                                    varMetadata = data.frame(labelDescription = colnames(featureData))))

head(expressions.rld.Raw)
#######---------------------------------------------
message("Remove NA in the ExpressionSets")
dim(expressions.rld.Raw)
# 去除含有NA的行，函数放在functions文件夹里
expressions.rld.Raw.rmNA <- removeNAsFromESet(expressions.rld.Raw)
dim(expressions.rld.Raw.rmNA)
# 两者一样

#######---------------------------------------------

####################################################################################
message(name, "Store the objects as RDS")
####################################################################################

saveRDS(object = expressions.rld.Scale,
        file = file.path(outputFolder, "RNA_rld_Scale.RDS"))
saveRDS(object = expressions.rld.Raw,
        file = file.path(outputFolder, "RNA_rld_Raw.RDS"))
####################################################################################
####################################################################################
message( "Clear environment")
####################################################################################

rm(list = ls())









####################################################################################
message("Make lung tissue metabolome ExpressionSets")
# 肺组织代谢组 生成 eSet
####################################################################################

dataDirectory <- "/data/home/mali/ke_data/metabolome"
outputFolder <- "/data/home/mali/ke_data/metabolome/RDS"
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)

# 读取代谢物表达矩阵 kmeans10.csv
all_raw_data <-   read_xlsx(file.path(dataDirectory,"ALL_sample_data.xlsx"))
all_raw_data <- all_raw_data[,1:46]
dim(all_raw_data)
# 删除重复行
all_raw_data <- all_raw_data[!duplicated(all_raw_data$Index), ]
dim(all_raw_data)
#> colnames(all_raw_data)
# [1] "Index"                 "Compounds"             "物质"                 
# [4] "Class I"               "物质一级分类"          "Class II"             
# [7] "物质二级分类"          "Q1 (Da)"               "Molecular weight (Da)"
#[10] "Ionization model"      "Formula"               "3W-Mock-1"            
#[13] "3W-Mock-2"             "3W-Mock-3"             "3W-Mock-4"            
#[16] "3W-Mock-5"             "3W-Mock-6"             "3W-Mock-7"            
#[19] "3W-Mock-8"             "3W-Mock-9"             "3W-3dpi-1"            
#[22] "3W-3dpi-2"             "3W-3dpi-3"             "3W-3dpi-4"            
#[25] "3W-3dpi-5"             "3W-3dpi-6"             "3W-3dpi-7"            
#[28] "3W-3dpi-8"             "3W-3dpi-9"             "3W-3dpi-10"           
#[31] "3W-3dpi-11"            "20M-Mock-1"            "20M-Mock-2"           
#[34] "20M-Mock-3"            "20M-Mock-4"            "20M-Mock-5"           
#[37] "20M-Mock-6"            "20M-Mock-7"            "20M-3dpi-1"           
#[40] "20M-3dpi-2"            "20M-3dpi-3"            "20M-3dpi-4"           
#[43] "20M-3dpi-5"            "20M-3dpi-6"            "20M-3dpi-7"           
#[46] "20M-3dpi-8"  
file_names <- dir(dataDirectory, pattern = glob2rx("*_info.xlsx"))

for (i in 1:length(file_names)) {
  name <- strsplit(file_names[i],".",fixed = T)
  name <- substr(name[[1]][1],1,nchar(name[[1]][1])-5) # 删除最后5个字符_info
  data <- read_xlsx(file.path(dataDirectory,file_names[i]))
  data <- data[,c("Index","P-value", "VIP","Log2FC", "Type")]
  colnames(data) <- c("Index",paste0(name[[1]][1],"_",c("pvalue","Variable_Importance","Log2FoldChange","ifsig")))

  all_raw_data  <- left_join(all_raw_data,data,by="Index")
}


####################################################################################

## scale
## 原始表达矩阵计算z score
rawCount <- all_raw_data %>%
    dplyr::select(Index, "3W-Mock-1":"20M-3dpi-8")%>%
  as.data.frame #%>%
  #rownames_to_column('Index') # 这里总显示有重复，但是并没有重复
rownames(rawCount) <- rawCount$Index
rawCount$Index  <- NULL

# 改原始矩阵名字
rawC <- rawCount %>%
  rename_at(1:length(colnames(.)), .funs = list(~paste0('Raw_', .)))

# scale mean count
scaleC <- rawCount %>%
  t %>%
  scale %>%
  t %>%
  as.data.frame %>%
  dplyr::rename_at(1:length(colnames(.)), .funs = list(~paste0('Scale_', .)))

scaleC_index <- scaleC %>%
  mutate(Index = rownames(.)) %>%
  dplyr::select(Index,everything())

####################################################################################
## Raw
# Assay data
exprs <- rawCount
# phenoData
# Sample annotation
pDataFile <- file.path(dataDirectory, "metadata_l_metabolome.csv")
pData <- read.csv(pDataFile)
rownames(pData) <- pData$ID
head(pData)
pData.ord <- pData

# featureData
featureData <-  all_raw_data %>% 
  dplyr::select(Index : Formula, "20M-Mock_vs_20M-3dpi_pvalue":"3W-Mock_vs_3W-3dpi_ifsig") %>% 
  column_to_rownames(var = 'Index') %>%
  mutate(ID = rownames(.)) %>%
  dplyr::select(ID,everything())
head(featureData)


# build up eSet
  expressions.rld.Raw <- ExpressionSet(assayData = as.matrix(exprs),
                                   phenoData = AnnotatedDataFrame(data = pData.ord,
                                                                  varMetadata = data.frame(labelDescription = colnames(pData.ord))),
                                   featureData = AnnotatedDataFrame(data = featureData,
                                                                    varMetadata = data.frame(labelDescription = colnames(featureData))))

head(expressions.rld.Raw)


## scale
# Assay data
exprs <- scaleC
# phenoData
# Sample annotation
pDataFile <- file.path(dataDirectory, "metadata_l_metabolome.csv")
pData <- read.csv(pDataFile)
pData$ID <- paste0('Scale_', pData$ID)
rownames(pData) <- pData$ID
head(pData)
pData.ord <- pData
# featureData
featureData <-  all_raw_data %>% 
  dplyr::select(Index : Formula, "20M-Mock_vs_20M-3dpi_pvalue":"3W-Mock_vs_3W-3dpi_ifsig") %>% 
  column_to_rownames(var = 'Index') %>%
  mutate(ID = rownames(.)) %>%
  dplyr::select(ID,everything())
head(featureData)
# build up eSet
  expressions.rld.Scale <- ExpressionSet(assayData = as.matrix(exprs),
                                   phenoData = AnnotatedDataFrame(data = pData.ord,
                                                                  varMetadata = data.frame(labelDescription = colnames(pData.ord))),
                                   featureData = AnnotatedDataFrame(data = featureData,
                                                                    varMetadata = data.frame(labelDescription = colnames(featureData))))

head(expressions.rld.Scale)



# 总表，有注释，raw count，logFC， scale count
data_all <- inner_join(all_raw_data,
  scaleC_index,by="Index") %>%
  write_csv(file.path(dataDirectory, "eachGroup_vs_Mock_k.csv"))

#######---------------------------------------------
message("Remove NA in the ExpressionSets")
dim(expressions.rld.Scale)
# 去除含有NA的行，函数放在functions文件夹里
expressions.rld.Scale.rmNA <- removeNAsFromESet(expressions.rld.Scale)
dim(expressions.rld.Scale.rmNA)
# 两者一样
####################################################################################

#######---------------------------------------------
message("Remove NA in the ExpressionSets")
dim(expressions.rld.Raw)
# 去除含有NA的行，函数放在functions文件夹里
expressions.rld.Raw.rmNA <- removeNAsFromESet(expressions.rld.Raw)
dim(expressions.rld.Raw.rmNA)
# 两者一样

#######---------------------------------------------

####################################################################################
message(name, "Store the objects as RDS")
####################################################################################

saveRDS(object = expressions.rld.Scale,
        file = file.path(outputFolder, "lung_metabolome_Scale.RDS"))
saveRDS(object = expressions.rld.Raw,
        file = file.path(outputFolder, "lung_metabolome_Raw.RDS"))
####################################################################################
####################################################################################
message( "Clear environment")
####################################################################################

rm(list = ls())



####################################################################################
message("Make plasma tissue metabolome ExpressionSets")
# 血浆代谢组 生成 eSet
####################################################################################

dataDirectory <- "/data/home/mali/ke_data/plasma_metabolome"
outputFolder <- "/data/home/mali/ke_data/plasma_metabolome/RDS"
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)

# 读取代谢物表达矩阵 kmeans10.csv
all_raw_data <-  read_xlsx(file.path(dataDirectory,"ALL_sample_data.xlsx"))
all_raw_data <- all_raw_data[,1:46]
dim(all_raw_data)
# 删除重复行
all_raw_data <- all_raw_data[!duplicated(all_raw_data$Index), ]
dim(all_raw_data)
#> colnames(all_raw_data)
# [1] "Index"                 "Compounds"             "物质"                 
# [4] "Class I"               "物质一级分类"          "Class II"             
# [7] "物质二级分类"          "Q1 (Da)"               "Molecular weight (Da)"
#[10] "Ionization model"      "Formula"               "3W-Mock-1"            
#[13] "3W-Mock-2"             "3W-Mock-3"             "3W-Mock-4"            
#[16] "3W-Mock-5"             "3W-Mock-6"             "3W-Mock-7"            
#[19] "3W-Mock-8"             "3W-Mock-9"             "3W-3dpi-1"            
#[22] "3W-3dpi-2"             "3W-3dpi-3"             "3W-3dpi-4"            
#[25] "3W-3dpi-5"             "3W-3dpi-6"             "3W-3dpi-7"            
#[28] "3W-3dpi-8"             "3W-3dpi-9"             "3W-3dpi-10"           
#[31] "3W-3dpi-11"            "20M-Mock-1"            "20M-Mock-2"           
#[34] "20M-Mock-3"            "20M-Mock-4"            "20M-Mock-5"           
#[37] "20M-Mock-6"            "20M-Mock-7"            "20M-3dpi-1"           
#[40] "20M-3dpi-2"            "20M-3dpi-3"            "20M-3dpi-4"           
#[43] "20M-3dpi-5"            "20M-3dpi-6"            "20M-3dpi-7"           
#[46] "20M-3dpi-8"  
file_names <- dir(dataDirectory, pattern = glob2rx("*_info.xlsx"))

for (i in 1:length(file_names)) {
  name <- strsplit(file_names[i],".",fixed = T)
  name <- substr(name[[1]][1],1,nchar(name[[1]][1])-5) # 删除最后5个字符_info
  data <- data <- read_xlsx(file.path(dataDirectory,file_names[i]))
  data <- data[,c("Index","P-value", "VIP","Log2FC", "Type")]
  colnames(data) <- c("Index",paste0(name[[1]][1],"_",c("pvalue","Variable_Importance","Log2FoldChange","ifsig")))

  all_raw_data  <- left_join(all_raw_data,data,by="Index")
}


####################################################################################

## scale
## 原始表达矩阵计算z score
rawCount <- all_raw_data %>%
    dplyr::select(Index, "3W-Mock-1":"20M-3dpi-8")%>%
  as.data.frame #%>%
  #rownames_to_column('Index') # 这里总显示有重复，但是并没有重复
rownames(rawCount) <- rawCount$Index
rawCount$Index  <- NULL

# 改原始矩阵名字
rawC <- rawCount %>%
  rename_at(1:length(colnames(.)), .funs = list(~paste0('Raw_', .)))

# scale mean count
scaleC <- rawCount %>%
  t %>%
  scale %>%
  t %>%
  as.data.frame %>%
  dplyr::rename_at(1:length(colnames(.)), .funs = list(~paste0('Scale_', .)))

scaleC_index <- scaleC %>%
  mutate(Index = rownames(.)) %>%
  dplyr::select(Index,everything())

####################################################################################
## Raw
# Assay data
exprs <- rawCount
# phenoData
# Sample annotation
pDataFile <- file.path(dataDirectory, "metadata_p_metabolome.csv")
pData <- read.csv(pDataFile)
rownames(pData) <- pData$ID
head(pData)
pData.ord <- pData

# featureData
featureData <-  all_raw_data %>% 
  dplyr::select(Index : Formula, "20M-Mock_vs_20M-3dpi_pvalue":"3W-Mock_vs_3W-3dpi_ifsig") %>% 
  column_to_rownames(var = 'Index') %>%
  mutate(ID = rownames(.)) %>%
  dplyr::select(ID,everything())
head(featureData)


# build up eSet
  expressions.rld.Raw <- ExpressionSet(assayData = as.matrix(exprs),
                                   phenoData = AnnotatedDataFrame(data = pData.ord,
                                                                  varMetadata = data.frame(labelDescription = colnames(pData.ord))),
                                   featureData = AnnotatedDataFrame(data = featureData,
                                                                    varMetadata = data.frame(labelDescription = colnames(featureData))))

head(expressions.rld.Raw)


## scale
# Assay data
exprs <- scaleC
# phenoData
# Sample annotation
pDataFile <- file.path(dataDirectory, "metadata_p_metabolome.csv")
pData <- read.csv(pDataFile)
pData$ID <- paste0('Scale_', pData$ID)
rownames(pData) <- pData$ID
head(pData)
pData.ord <- pData
# featureData
featureData <-  all_raw_data %>% 
  dplyr::select(Index : Formula, "20M-Mock_vs_20M-3dpi_pvalue":"3W-Mock_vs_3W-3dpi_ifsig") %>% 
  column_to_rownames(var = 'Index') %>%
  mutate(ID = rownames(.)) %>%
  dplyr::select(ID,everything())
head(featureData)
# build up eSet
  expressions.rld.Scale <- ExpressionSet(assayData = as.matrix(exprs),
                                   phenoData = AnnotatedDataFrame(data = pData.ord,
                                                                  varMetadata = data.frame(labelDescription = colnames(pData.ord))),
                                   featureData = AnnotatedDataFrame(data = featureData,
                                                                    varMetadata = data.frame(labelDescription = colnames(featureData))))

head(expressions.rld.Scale)



# 总表，有注释，raw count，logFC， scale count
data_all <- inner_join(all_raw_data,
  scaleC_index,by="Index") %>%
  write_csv(file.path(dataDirectory, "eachGroup_vs_Mock_k.csv"))

#######---------------------------------------------
message("Remove NA in the ExpressionSets")
dim(expressions.rld.Scale)
# 去除含有NA的行，函数放在functions文件夹里
expressions.rld.Scale.rmNA <- removeNAsFromESet(expressions.rld.Scale)
dim(expressions.rld.Scale.rmNA)
# 两者一样
####################################################################################

#######---------------------------------------------
message("Remove NA in the ExpressionSets")
dim(expressions.rld.Raw)
# 去除含有NA的行，函数放在functions文件夹里
expressions.rld.Raw.rmNA <- removeNAsFromESet(expressions.rld.Raw)
dim(expressions.rld.Raw.rmNA)
# 两者一样

#######---------------------------------------------

####################################################################################
message(name, "Store the objects as RDS")
####################################################################################

saveRDS(object = expressions.rld.Scale,
        file = file.path(outputFolder, "plasma_metabolome_Scale.RDS"))
saveRDS(object = expressions.rld.Raw,
        file = file.path(outputFolder, "plasma_metabolome_Raw.RDS"))
####################################################################################
####################################################################################
message( "Clear environment")
####################################################################################

rm(list = ls())

