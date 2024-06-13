
#######################GO analysis############################

library('goseq') #BiocManager::install("goseq")
# 用于寻找GO terms，即基因富集分析。此方法基于 Wallenius non-central hyper-geometric distribution。相对于普通的超几何分布(Hyper-geometric distribution)，此分布的特点是从某个类别中抽取个体的概率与从某个类别之外抽取一个个体的概率是不同的，这种概率的不同是通过对基因长度的偏好性进行估计得到的，从而能更为准确地计算出 GO term 被差异基因富集的概率。
library('GO.db')  #所有这些后缀为.db的R包，其本质都为一个sqlite数据库，一种轻量级的关系型数据库，只不过是通过R来进行访问。
library('foreach')
library('doMC') # 使用并行包的多核功能为%dopar%函数提供一个并行后端。
library('KEGGAPI')
#依赖  devtools::install_github("YulongNiu/ParaMisc")
# devtools::install_github("YulongNiu/KEGGAPI")
library('BioCycAPI')
# 依赖KEGGAPI
# devtools::install_github("YulongNiu/BioCycAPI")
library('magrittr')
library('dplyr')
library('tibble')
library('readr')
library('stringr')
library('RCurl')
library('xml2')
library('urltools')

library('tidyverse')
library('DESeq2')
library('ComplexHeatmap')
library('RColorBrewer')
library('circlize')
library('org.Mm.eg.db') # library('org.At.tair.db')  # BiocManager::install("org.At.tair.db")
library('clusterProfiler')
library('magrittr')
library('tidyverse')
library('enrichplot')
library('ggnewscale')
registerDoMC(12)

load('/data/home/mali/genome/mouse/mmuGO.RData') # athGO
load('/data/home/mali/genome/mouse/mmuKEGG.RData') # athKEGG
load('/data/home/mali/genome/mouse/mouseBioCyc_entrez.RData') # athBioCyc_entrez

kmeansRes <- read_csv('./transcription/count/kmeans10_sig.csv',
                      col_types = cols(Chromosome = col_character()))
kmeansBkg <- read_csv('./transcription/count/kmeans10.csv',
                      col_types = cols(Chromosome = col_character()))
> head(kmeansRes)
# A tibble: 6 × 35
         ID Gene  gene_…¹ gene_…² gene_…³ gene_…⁴ gene_…⁵ gene_…⁶ Descr…⁷    YV1
      <dbl> <chr> <chr>     <dbl>   <dbl> <chr>     <dbl> <chr>   <chr>    <dbl>
1     20389 Sftpc NC_000…  7.08e7  7.08e7 -           804 protei… surfac… 8762.
2     22287 Scgb… NC_000…  9.06e6  9.07e6 -           468 protei… secret… 2567.
3 100503605 Hbb-… NC_000…  1.03e8  1.03e8 -           639 protei… hemogl… 8422.
4     15122 Hba-… NC_000…  3.22e7  3.22e7 +           559 protei… hemogl… 6853.
5     16071 Igkc  NC_000…  7.07e7  7.07e7 +           320 C_regi… - && P…   91.8
6     17105 Lyz2  NC_000…  1.17e8  1.17e8 -          1057 protei… lysozy… 2751.
# … with 25 more variables: YV2 <dbl>, YV3 <dbl>, YV4 <dbl>, OV1 <dbl>,
#   OV2 <dbl>, OV3 <dbl>, O1 <dbl>, O2 <dbl>, O3 <dbl>, Y1 <dbl>, Y2 <dbl>,
#   Y3 <dbl>, YV_vs_Y_pvalue <dbl>, YV_vs_Y_padj <dbl>,
#   YV_vs_Y_log2FoldChange <dbl>, OV_vs_O_pvalue <dbl>, OV_vs_O_padj <dbl>,
#   OV_vs_O_log2FoldChange <dbl>, O_vs_Y_pvalue <dbl>, O_vs_Y_padj <dbl>,
#   O_vs_Y_log2FoldChange <dbl>, OV_vs_YV_pvalue <dbl>, OV_vs_YV_padj <dbl>,
#   OV_vs_YV_log2FoldChange <dbl>, cl <dbl>, and abbreviated variable names …
# ℹ Use `colnames()` to see all variable names

# > dim(kmeansBkg)
# [1] 12088    61
# > dim(kmeansRes)
# [1] 1065   35
#> length(athGO)
#[1] 19306

useVersion   <- "version1"
outputFolder <- file.path("/data/home/mali/ke_data/output", useVersion, "geneset")
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)

savepath <- outputFolder

##~~~~~~~~~~~~~~~select genesets~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# athGO中kmeansBkg$ID存在的部分
athGO %<>%
  lapply(function(x) {x[x %in% kmeansBkg$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }
length(athGO)
#[1] 17002
name_length <- foreach(i = 1:length(athGO), .combine = rbind) %dopar% {
  eachMat <- names(athGO)[i]
  return(eachMat)
} %>% {.[!is.na(.)]}
length(name_length)
# [1] 17002
GOMat <- foreach(i = 1:length(athGO), .combine = rbind) %dopar% {
  eachMat <- cbind(athGO[[i]], names(athGO)[i])
  return(eachMat)
} %>% as.data.frame
# > head(GOMat)
#       V1         V2
# 1 170942 GO:0000001
# 2 225523 GO:0000001
# 3 269023 GO:0000001
# 4 107022 GO:0000001
# 5 110695 GO:0000001
# 6  16906 GO:0000001


athKEGG %<>%
  lapply(function(x) {x[x %in% kmeansBkg$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }
length(athKEGG)
# [1] 347
name_length <- foreach(i = 1:length(athKEGG), .combine = rbind) %dopar% {
  eachMat <- names(athKEGG)[i]
  return(eachMat)
} %>% {.[!is.na(.)]}
length(name_length)
# [1] 347

KEGGMat <- foreach(i = 1:length(athKEGG), .combine = rbind) %dopar% {
  eachMat <- cbind(athKEGG[[i]], names(athKEGG)[i])
  return(eachMat)
} %>% as.data.frame
# > head(KEGGMat)
#      V1            V2
# 1 15277 path:mmu00010
# 2 11522 path:mmu00010
# 3 16828 path:mmu00010
# 4 11532 path:mmu00010
# 5 18648 path:mmu00010
# 6 11670 path:mmu00010

athBioCyc_entrez %<>%
  lapply(function(x) {x[x %in% kmeansBkg$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }
length(athBioCyc_entrez)
# [1] 295
name_length <- foreach(i = 1:length(athBioCyc_entrez), .combine = rbind) %dopar% {
  eachMat <- names(athBioCyc_entrez)[i]
  return(eachMat)
} %>% {.[!is.na(.)]}
length(name_length)
# [1] 295

BioCycMat <- foreach(i = 1:length(athBioCyc_entrez), .combine = rbind) %dopar% {
  eachMat <- cbind(athBioCyc_entrez[[i]], names(athBioCyc_entrez)[i])
  return(eachMat)
} %>% as.data.frame
# > head(BioCycMat)
#       V1                           V2
# 1 110460   MOUSE:ACETOACETATE-DEG-PWY
# 2 110446   MOUSE:ACETOACETATE-DEG-PWY
# 3  67041   MOUSE:ACETOACETATE-DEG-PWY
# 4 238505 MOUSE:ADENOSYLHOMOCYSCAT-PWY
# 5  66925          MOUSE:AERESPDON-PWY
# 6  67680          MOUSE:AERESPDON-PWY
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~whole cluster gene-set~~~~~~~~~~~~~~~~~    需要重新加载一下 getKEGGPathAnno 函数，覆盖掉KEGGAPI中的函数才行

##' KEGG Database API - Get the whole pathway ID from KEGG database
##'
##' Get the pathway ID and annoation of a given KEGG species ID.
##' @title List pathway of a given species ID
##' @param specID KEGSS species org code or T number , for example "hsa" or "T01001".
##' @return A matrix of pathway ID and annotation.
##' @examples
##' hasPath <- getKEGGPathAnno('hsa')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @references \url{http://www.kegg.jp/kegg/rest/}
##' @export
getKEGGPathAnno <- function(specID){

  ## get KEGG pathway annotation list
  #url <- paste('http://rest.kegg.jp/list/pathway/', specID, sep = '') # 下面报错是因为这里网址错误，http + s 解决
  url <- paste('https://rest.kegg.jp/list/pathway/', specID, sep = '')
  pathAnno <- webTable(url, ncol = 2) # 是这一句报错

  colnames(pathAnno) <- c('pathID', 'Annotation')

  return(pathAnno)
}


##' BioCyc Database API - Pathway
##'
##' \itemize{
##'   \item \code{getCycPathway()}: Get whole pathways list.
##'   \item \code{getCycGenesfPathway()}: Get genes from a given pathway ID.
##' }
##'
##' It may take more than 10 minutes to retrieve the xml file.
##' @title Pathway
##' @inheritParams getCycGenes
##' @return
##' \itemize{
##'   \item \code{getCycPathway()}:  A \code{tbl_df} contains 1st column is pathway ID and 2nd column is pathway annotation.
##'   \item \code{getCycGenesfPathway()}: A \code{list} indicates genes, proteins, and common names.
##' }
##'
##' @examples
##' ## Candida albicans SC5314 pathways
##' calTU <- getCycPathway('CALBI')
##'
##' ## genes of pathway "CALBI:PWY3B3-8"
##' genes <- getCycGenesfPathway('CALBI:PWY3B3-8')
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom urltools url_encode
##' @importFrom xml2 read_xml xml_text xml_find_all
##' @importFrom magrittr %>%
##' @importFrom tibble tibble
##' @rdname pathway
##' @export
##'
getCycPathway <- function(speID,type = 'pathways') {
  url <- paste0('https://websvc.biocyc.org/xmlquery?%5bx:x%3C-', speID,'^^',type,'%5d')
  file_name <- paste0(speID,"_",type,".xml")
  shell_cmd <-paste0("./biocyc_curl_get.sh ","'", url,"' ","'", file_name,"'")
  system(shell_cmd, intern = FALSE, ignore.stderr=TRUE)
  pathxml <- read_xml(file_name)

  pathnodes <- pathxml %>%
    xml_find_all('/ptools-xml/Pathway')

  pathanno <- sapply(pathnodes, function(x){
    eachanno <- xml_find_all(x, 'common-name') %>%
      xml_text %>%
      ifelse(length(.) == 0, '', .)
    return(eachanno)
  })

  res <- tibble(pathID = xml_attr(pathnodes, 'ID'),
                pathAnno = pathanno)

  return(res)
}


####----------------------------------------------------
# 每个cluster生成GO,KEGG,biocyc的csv文件
####----------------------------------------------------

### kegg注释

pathAnno <- getKEGGPathAnno('mmu') %>%  # Get the whole pathway ID from KEGG database； 'ath'： species org code or T number , for example "hsa" or "T01001".
  as_tibble
# head(pathAnno)
pathAnno_KEGG <- pathAnno %>%
  dplyr::mutate(Annotation = Annotation %>% substr(., 1, nchar(.) - 29)) # 字符数量是根据这个调整的 >>> - Mus musculus (house mouse)
# > head(pathAnno_KEGG)
# # A tibble: 6 × 2
#   pathID        Annotation
#   <chr>         <chr>
# 1 path:mmu00010 Glycolysis / Gluconeogenesis
# 2 path:mmu00020 Citrate cycle (TCA cycle)
# 3 path:mmu00030 Pentose phosphate pathway
# 4 path:mmu00040 Pentose and glucuronate interconversions
# 5 path:mmu00051 Fructose and mannose metabolism
# 6 path:mmu00052 Galactose metabolism

### BioCyc注释

pathAnno_Cyc <- getCycPathway('MOUSE') %>%  # Get whole pathways list. speID: The BioCyc species ID, for example "ECOLI" is for"Escherichia coli K-12 substr. MG1655".
  dplyr::rename(Annotation = pathAnno) %>%
  mutate(Annotation = Annotation %>% str_replace_all('<.*?>', ''))
# > head(pathAnno_Cyc)
# # A tibble: 6 × 2
#   pathID          Annotation
#   <chr>           <chr>
# 1 MOUSE:URSIN-PWY ureide biosynthesis
# 2 MOUSE:PWY-5138  fatty acid &beta;-oxidation IV (unsaturated, even number)
# 3 MOUSE:PWY-4983  citrulline-nitric oxide cycle
# 4 MOUSE:PWY-5329  L-cysteine degradation III
# 5 MOUSE:PWY-6573  chondroitin sulfate degradation (metazoa)
# 6 MOUSE:PWY0-662  PRPP biosynthesis I


for (i in kmeansBkg$cl %>% unique) {

  prefix <- 'kmeans10'

  degVec <- (kmeansBkg$cl == i) %>%
    as.integer %>%
    set_names(kmeansBkg$ID)

  pwf <- nullp(degVec, bias.data = kmeansBkg$gene_length) # 这个函数要查一下含义
#  > head(pwf)
#            DEgenes bias.data       pwf
#  20389           0       804 0.2105356
#  22287           0       468 0.2494002
#  100503605       0       639 0.2293362
#  15122           0       559 0.2386593
#  16071           1       320 0.2671154
#  17105           0      1057 0.1831204
# > dim(pwf)
# [1] 12088     3
  ## GO
  GOTestWithCat <- goseq(pwf, gene2cat = GOMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    filter(!is.na(ontology))
# For 544 genes, we could not find any categories. These genes will be excluded.
  termCat <- c('BP', 'MF', 'CC')
  for (j in termCat) {
    write.csv(GOTestWithCat %>% filter(ontology == j),
              paste0(prefix, '_cluster', i, '_', j, '_unbkg.csv') %>% file.path(savepath, .))
  }

  ## KEGG

  KEGGTestWithCat <- goseq(pwf, gene2cat = KEGGMat, use_genes_without_cat = FALSE) %>% # pwf: An object containing gene names, DE calls, the probability weighting function. Usually generated by ‘nullp’.
    as_tibble %>%
    inner_join(., pathAnno_KEGG, by = c('category' = 'pathID')) %>%
    dplyr::mutate(ontology = 'KEGG')
# For 7128 genes, we could not find any categories. These genes will be excluded
  write.csv(KEGGTestWithCat,
            paste0(prefix, '_cluster', i, '_KEGG_unbkg.csv') %>% file.path(savepath, .))

  ## 代谢组
  BioCycTestWithCat <- goseq(pwf, gene2cat = BioCycMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    inner_join(., pathAnno_Cyc, by = c('category' = 'pathID')) %>%
    mutate(ontology = 'BioCyc')
# For 11316 genes, we could not find any categories. These genes will be excluded.
  write.csv(BioCycTestWithCat,
            paste0(prefix, '_cluster', i, '_BioCyc_unbkg.csv') %>% file.path(savepath, .))
}
#

# > head(GOTestWithCat)
# # A tibble: 6 × 7
#   category   over_represented_pvalue under_repre…¹ numDE…² numIn…³ term  ontol…⁴
#   <chr>                        <dbl>         <dbl>   <int>   <int> <chr> <chr>
# 1 GO:0010833               0.0000416          1.00       4       6 telo… BP
# 2 GO:0051224               0.000965           1.00       3       6 nega… BP
# 3 GO:0006486               0.00104            1.00      13     122 prot… BP
# 4 GO:0097502               0.00106            1.00       4      12 mann… BP
# 5 GO:0097035               0.00121            1.00       3       6 regu… BP
# 6 GO:2000642               0.00136            1          2       2 nega… BP
# # … with abbreviated variable names ¹​under_represented_pvalue, ²​numDEInCat,
# #   ³​numInCat, ⁴​ontology

#> colnames(GOTestWithCat)
#[1] "category"                 "over_represented_pvalue"
#[3] "under_represented_pvalue" "numDEInCat"
#[5] "numInCat"                 "term"
#[7] "ontology"    #

#head(KEGGTestWithCat) # 和上面的一样
#head(BioCycTestWithCat)
#colnames(BioCycTestWithCat)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



##~~~~~~~~~~~~~~~~~~~whole cluster gene-set with background~~~~~~~~~~
# 指定background，生成每个cluster 的go，kegg的csv文件
for (i in kmeansRes$cl %>% unique) {

  prefix <- 'kmeans10'

  eachRes <- kmeansRes %>%
    filter(cl == i) %>%
    {.$ID}
  eachBkg <- kmeansBkg %>%
    filter(cl == i) %>%
    {.$ID}
  eachLength <- kmeansBkg %>%
    filter(cl == i) %>%
    {.$gene_length}

  degVec <- rep(0, length(eachBkg)) %>%
    set_names(eachBkg)
  degVec[match(eachRes, eachBkg)] <- 1

  pwf <- nullp(degVec, bias.data = eachLength)

  ## GO
  GOTestWithCat <- goseq(pwf, gene2cat = GOMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    filter(!is.na(ontology))

  termCat <- c('BP', 'MF', 'CC')
  for (j in termCat) {
    write.csv(GOTestWithCat %>% filter(ontology == j),
              paste0(prefix, '_cluster', i, '_', j, '.csv') %>% file.path(savepath, .))
  }

  ## KEGG
  KEGGTestWithCat <- goseq(pwf, gene2cat = KEGGMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    inner_join(., pathAnno_KEGG, by = c('category' = 'pathID')) %>%
    mutate(ontology = 'KEGG')

  write.csv(KEGGTestWithCat,
            paste0(prefix, '_cluster', i, '_KEGG.csv') %>% file.path(savepath, .))

  ## BioCyc
  BioCycTestWithCat <- goseq(pwf, gene2cat = BioCycMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    inner_join(., pathAnno_Cyc, by = c('category' = 'pathID')) %>%
    mutate(ontology = 'BioCyc')
  write.csv(BioCycTestWithCat,
            paste0(prefix, '_cluster', i, '_BioCyc.csv') %>% file.path(savepath, .))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~choose sig~~~~~~~~~~~~~~~~~~~~~
## padj < 0.05 & |log2FC| > log2(1.5)
fcsig <- kmeansRes %>%
  dplyr::select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > 1 ~ 1,
                                 . < -1 ~ -1,
                                 TRUE ~ 0)))
# > colnames(fcsig)
# [1] "YV_vs_Y_log2FoldChange"  "OV_vs_O_log2FoldChange"
# [3] "O_vs_Y_log2FoldChange"   "OV_vs_YV_log2FoldChange"

padjsig <- kmeansRes %>%
  dplyr::select(ends_with('padj')) %>%
  abs %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))
#> colnames(padjsig)
# [1] "YV_vs_Y_padj"  "OV_vs_O_padj"  "O_vs_Y_padj"   "OV_vs_YV_padj"

sig <- (padjsig * fcsig) %>%
  as_tibble %>%
  mutate(ID = kmeansRes$ID, cl = kmeansRes$cl) %>%
  dplyr::select(ID, everything()) %>%
  dplyr::rename(YoungInfected_vs_YoungMock = YV_vs_Y_padj,
         OldInfected_vs_OldMock = OV_vs_O_padj,
         OldMock_vs_YoungMock = O_vs_Y_padj,
         OldInfected_vs_YoungInfected = OV_vs_YV_padj)
# 选择想要分析的组，改组名，新名字=旧名字

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 选择想要分析的组
vsGroup <- c('YoungInfected_vs_YoungMock', 'OldInfected_vs_OldMock', 'OldMock_vs_YoungMock','OldInfected_vs_YoungInfected')
cln <- 1:10 # 所有cluster结果都有
cln <- 4  # 可以具体看到哪个cluster的聚类

for (i in vsGroup) {
  for (j in cln) {
    degVec <- sig %>%
      transmute((!!as.name(i)) != 0 &
                cl == j) %>%
      unlist %>%
      as.integer
    names(degVec) <- sig$ID

    pwf <- nullp(degVec, bias.data = kmeansRes$gene_length)

    GOTestWithCat <- goseq(pwf, gene2cat = GOMat, use_genes_without_cat = FALSE) %>%
      as_tibble %>%
      filter(!is.na(ontology)) %>%
      dplyr::rename(Annotation = term)

    termCat <- c('BP', 'MF', 'CC')
    for (k in termCat) {
      write.csv(GOTestWithCat %>% filter(ontology == k),
                paste0(prefix,'_', i, '_cluster', j, '_', k, '.csv') %>% file.path(savepath, .))
    }
  }
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~KEGG~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (i in vsGroup) {
  for (j in cln) {
    degVec <- sig %>%
      transmute((!!as.name(i)) != 0 &
                cl == j) %>%
      unlist %>%
      as.integer
    names(degVec) <- sig$ID

    pwf <- nullp(degVec, bias.data = kmeansRes$gene_length)

    KEGGTestWithCat <- goseq(pwf, gene2cat = KEGGMat, use_genes_without_cat = FALSE) %>%
      as_tibble %>%
      inner_join(., pathAnno_KEGG, by = c('category' = 'pathID')) %>%
      mutate(ontology = 'KEGG')

    write.csv(KEGGTestWithCat,
              paste0(prefix,'_', i, '_cluster', j, '_KEGG', '.csv') %>% file.path(savepath, .))
  }
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

###~~~~~~~~~~~~~~~~~~~~~~~~~BioCyc~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (i in vsGroup) {
  for (j in cln) {
    degVec <- sig %>%
      transmute((!!as.name(i)) != 0 &
                cl == j) %>%
      unlist %>%
      as.integer
    names(degVec) <- sig$ID#
    pwf <- nullp(degVec, bias.data = kmeansRes$gene_length)#
    BioCycTestWithCat <- goseq(pwf, gene2cat = BioCycMat, use_genes_without_cat = FALSE) %>%
      as_tibble %>%
      inner_join(., pathAnno_Cyc, by = c('category' = 'pathID')) %>%
      mutate(ontology = 'BioCyc')#
    write.csv(BioCycTestWithCat,
              paste0(prefix,'_', i, '_cluster', j, '_BioCyc', '.csv') %>% file.path(savepath, .))
  }
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################################################


###############################cluster profiler#####################
library('org.Mm.eg.db') # library('org.At.tair.db')  # BiocManager::install("org.At.tair.db")
library('clusterProfiler')
library('magrittr')
library('tidyverse')
library('enrichplot')
library('ggnewscale')

useVersion   <- "version1"
outputFolder <- file.path("/data/home/mali/ke_data/output", useVersion, "clusterbc")
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)

savepath <- outputFolder

kmeansRes <- read_csv('./transcription/count/kmeans10_sig.csv',
                      col_types = cols(Chromosome = col_character()))
kmeansBkg <- read_csv('./transcription/count/kmeans10.csv',
                      col_types = cols(Chromosome = col_character()))



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
  write.csv(as.data.frame(goBPSim),
            paste0(prefix, '_cluster', i, '_cp_BP.csv') %>% file.path(savepath, .))

  ## KEGG
  kk2 <- enrichKEGG(gene = kmeansRes %>% filter(cl == i) %>% .$ID %>% unlist %>% unique,
                    organism = 'mmu',
                    pvalueCutoff = 0.05)

  write.csv(as.data.frame(kk2),
            paste0(prefix, '_cluster', i, '_cp_KEGG.csv') %>% file.path(savepath, .))
}

kall <- lapply(kmeansRes$cl %>% unique, function(x) {

  eachG <- kmeansRes %>% filter(cl == x) %>% .$ID %>% unlist %>% unique

  return(eachG)

}) %>%
  set_names(kmeansRes$cl %>% unique %>% paste0('cluster', .))

save(kall, file = paste0('./transcription/count/',prefix,'_kall.RData'))
# > head(kall)
# $cluster8
#  [1]  20389  20387 227753  76293 110257  12389  74318 211623  55990  67888
# [11]  21824  56492  20324  14199  54598  19285  18247  74761 114249  66109
# [21]  20250  11522 239739  12520 210992
# ......
# $cluster5
# ......

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
ggsave(file.path(savepath,paste0(prefix,'_cp_BP_dotplot_Sim_10.jpg')), width = 20)
ggsave(file.path(savepath,paste0(prefix,'_cp_BP_dotplot_Sim_10.pdf')), width = 20)

dotplot(kallGOBP, showCategory = 10)
ggsave(file.path(savepath,paste0(prefix,'_cp_BP_dotplot_10.jpg')), width = 20)
ggsave(file.path(savepath,paste0(prefix,'_cp_BP_dotplot_10.pdf')), width = 20)

# 保存到csv
kallGOBP %>%
  as.data.frame %>%
  write_csv(paste0('./transcription/count/',prefix,'_cp_BP.csv'))
kallGOBPSim %>%
  as.data.frame %>%
  write_csv(paste0('./transcription/count/',prefix,'_cp_BP_sim.csv'))
# 保存完整结果
save(kallGOBP, file = paste0('./transcription/count/',prefix,'_cp_BP.RData'))
save(kallGOBPSim, file = paste0('./transcription/count/',prefix,'_cp_BP_sim.RData'))

# 绘制网络图
kallGOBP2 <- pairwise_termsim(kallGOBP)
emapplot(kallGOBP2,
         showCategory = 10,
         pie='count',
         pie_scale=1.5,
         layout='nicely')
ggsave(file.path(savepath,paste0(prefix,'_cp_BP_network_10.jpg')), width = 20, height = 20)
ggsave(file.path(savepath,paste0(prefix,'_cp_BP_network_10.pdf')), width = 20, height = 20)
#-------------------
#-------------------
# 绘制compareCluster的kegg dotplot
kallKEGG <- compareCluster(geneCluster = kall,
                           fun = 'enrichKEGG',
                           organism = 'mmu',
                           pvalueCutoff = 0.05)
dotplot(kallKEGG, showCategory = 10)
ggsave(file.path(savepath,paste0(prefix,'_cp_KEGG_dotplot_10.jpg')), width = 20, height = 20)
ggsave(file.path(savepath,paste0(prefix,'_cp_KEGG_dotplot_10.pdf')), width = 20, height = 20)

# simplify only work for GO...
#kallKEGGSim 无

# 保存到csv
kallKEGG %>%
  as.data.frame %>%
  write_csv(paste0('./transcription/count/',prefix,'_cp_KEGG.csv'))

# 保存完整结果
save(kallKEGG, file = paste0('./transcription/count/',prefix,'_cp_KEGG.RData'))


# 绘制网络图
kallKEGG2 <- pairwise_termsim(kallKEGG)
emapplot(kallKEGG2,
         showCategory = 10,
         pie='count',
         pie_scale=1.5,
         layout='nicely')
ggsave(file.path(savepath,paste0(prefix,'_cp_KEGG_network_10.jpg')), width = 20, height = 20)
ggsave(file.path(savepath,paste0(prefix,'_cp_KEGG_network_10.pdf')), width = 20, height = 20)
#-------------------
#######################################################################

###################################plot###########################
# 读取csv文件，批量绘图，每个cluster 的go，kegg图
library('ggplot2')
library('readr')
library('dplyr')
library('magrittr')
library('foreach')
library('stringr')
cln <- 1:10
#cln <- 4
geneset <- c('BP', 'MF', 'CC', 'KEGG', 'BioCyc')

useVersion   <- "version1"
outputFolder <- file.path("/data/home/mali/ke_data/output", useVersion, "geneset")
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)
savepath <- outputFolder
prefix <- 'kmeans10'


for (i in cln) {
  for (j in geneset) {

    vsGroup <- c('YoungInfected_vs_YoungMock', 'OldInfected_vs_OldMock', 'OldMock_vs_YoungMock','OldInfected_vs_YoungInfected')

    pathPlot <- foreach(k = seq_along(vsGroup), .combine = bind_rows) %do% {
      vsGroup[k] %>%
        {paste0(prefix,'_', ., '_cluster', i, '_', j, '.csv') %>% file.path(savepath, .)} %>%
        read_csv %>%
        dplyr::select(Annotation, over_represented_pvalue, numDEInCat, numInCat) %>%
        dplyr::rename(pvalue = over_represented_pvalue) %>%
        mutate(group = vsGroup[k], ratio = numDEInCat / numInCat) %>%
        filter(pvalue < 0.05 &
               numDEInCat >= 1)
    }
    # 判断是否还有数据
    if(length(pathPlot$pvalue)==0){
      print("No data remain!")
    }else{
      # 去掉0值
      # 检查第一列中是否有0值
      if (identical(min(pathPlot[,"pvalue"]),0)) {
        pathPlot <- pathPlot[-which(pathPlot[,"pvalue"] == 0),]
        }else{
        pathPlot <- pathPlot
        }
      # 判断是否还有数据
      if(length(pathPlot$pvalue)==0){
        print("No data remain!")
      }else{

      colorPal <- colorRampPalette(rev(c('red', 'yellow', 'cyan', 'blue')), bias=1)(10)
      ggplot(pathPlot, aes(x = group, y = Annotation)) +
        geom_point(aes(size = ratio, colour = -log10(pvalue))) +
        scale_colour_gradientn(name = '-log10(P-value)', limits=c(0, max(-log10(pathPlot$pvalue))), colours = colorPal) +
        ## scale_x_discrete(labels = c('Flg22', 'Flg22+SynCom33', 'Flg22+SynCom35')) +
        ylab(j) +
        xlab('') +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
      ggsave(file.path(savepath,paste0(prefix,'_cluster', i, '_', j, '.pdf')),width = 20,height = 20)
      ggsave(file.path(savepath,paste0(prefix,'_cluster', i, '_', j, '.jpg')),width = 20,height = 20)
      }
    }
  }
}
##################################################################


##############################Plot GO heatmap##########################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~~~~
meanFlg22 <- function(v) {

  require('magrittr')

  res <- v %>%
    split(c(rep(c(1:2), each = 3),rep(3, each = 4),rep(4, each = 3))) %>%
    sapply(mean, na.rm = TRUE)

  return(res)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library('tidyverse')
library('DESeq2')
library('ComplexHeatmap')
library('RColorBrewer')
library('circlize')
library('org.Mm.eg.db') # library('org.At.tair.db')  # BiocManager::install("org.At.tair.db")
library('clusterProfiler')
library('magrittr')
library('enrichplot')
library('ggnewscale')
library(AnnotationDbi)


prefix <- "kmeans10"
useVersion   <- "version1"
outputFolder <- file.path("/data/home/mali/ke_data/output", useVersion, "clusterbc")
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)

savepath <- outputFolder

load(paste0('./transcription/count/',prefix,'_cp_BP.RData'))

topGONum <- 10
# 加载原始蛋白表达水平
protein_eSet <- readRDS("/data/home/mali/ke_data/output/version1/consensus_leukemic_clustering/all/CLCs/all/proteins.RDS")
colnames(protein_eSet) <- protein_eSet$name
# 加载 DESeq处理后的eSet数据
RNA_eSet_Scale <- readRDS("/data/home/mali/ke_data/transcription/RDS/RNA_rld_Scale.RDS")
colnames(RNA_eSet_Scale) <- RNA_eSet_Scale$name
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Sig terms~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 加载gene注释文件
anno <- read_csv('/data/home/mali/genome/mouse/Ensembl_mmu_Anno_AddGeneName.csv',
                    col_types = cols(Chromosome = col_character())) %>%
  mutate_all(list(~replace(., is.na(.), ''))) %>%
  mutate(GeneID = entrezgene_id) %>%  # 因为要跟下面的数据对应上，下面用的是entrezgene_id
  dplyr::select(GeneID, Gene, Description) %>%
  dplyr::slice(which(!duplicated(.)))

# > head(anno)
## A tibble: 6 × 3
#  GeneID Gene  Description
#  <chr>  <chr> <chr>
#1 14679  Gnai3 guanine nucleotide binding protein (G protein), alpha inhibiting…
#2 54192  Pbsn  probasin [Source:MGI Symbol;Acc:MGI:1860484]
#3 15417  Hoxb9 homeobox B9 [Source:MGI Symbol;Acc:MGI:96190]
#4 12544  Cdc45 cell division cycle 45 [Source:MGI Symbol;Acc:MGI:1338073]
#5 16002  Igf2  insulin-like growth factor 2 [Source:MGI Symbol;Acc:MGI:96434]
#6 11818  Apoh  apolipoprotein H [Source:MGI Symbol;Acc:MGI:88058]

cpBP <- "fortify"(kallGOBP ,showCategory = topGONum) %>%  # convert compareClusterResult to a data.frame that ready for plot
  as_tibble %>%
  dplyr::mutate(geneName = sapply(geneID, function(x) {
    strsplit(x, split = '/', fixed = TRUE) %>%
      unlist %>%
      tibble(GeneID = .) %>%
      inner_join(anno) %>%
      .$Gene
  }))

# > head(cpBP)
# # A tibble: 6 × 11
#   Cluster  ID    Descr…¹ GeneR…² BgRatio   pvalue p.adjust   qvalue geneID Count
#   <fct>    <chr> <fct>     <dbl> <chr>      <dbl>    <dbl>    <dbl> <chr>  <int>
# 1 "cluste… GO:0… extrac…  0.177  322/28… 2.18e-13 1.43e-10 1.03e-10 76293…    14
# 2 "cluste… GO:0… extrac…  0.177  323/28… 2.27e-13 1.43e-10 1.03e-10 76293…    14
# 3 "cluste… GO:0… extern…  0.177  324/28… 2.37e-13 1.43e-10 1.03e-10 76293…    14
# 4 "cluste… GO:0… extrac…  0.0633 44/289… 1.33e- 7 6.00e- 5 4.35e- 5 76293…     5
# 5 "cluste… GO:1… regula…  0.0380 10/289… 2.32e- 6 8.36e- 4 6.06e- 4 12389…     3
# 6 "cluste… GO:0… toxin …  0.0380 12/289… 4.23e- 6 1.09e- 3 7.90e- 4 55990…     3
# # … with 1 more variable: geneName <named list>, and abbreviated variable names
# #   ¹​Description, ²​GeneRatio
# # ℹ Use `colnames()` to see all variable names
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~heatmap matrix~~~~~~~~~~~~~~~~~~~~~~
load('/data/home/mali/ke_data/transcription/RDS/degres_condi_Mock.RData')

wholeDEG <- read_csv('/data/home/mali/ke_data/transcription/count/eachGroup_vs_Mock_k.csv')
kmeansRes <- read_csv('/data/home/mali/ke_data/transcription/count/kmeans10.csv') %>%
  select(ID, cl)

fcsig <- wholeDEG %>%
  select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > log2(1.5) ~ 1,
                                 . < -log2(1.5) ~ -1,
                                 TRUE ~ 0)))
padjsig <- wholeDEG %>%
  select(ends_with('padj')) %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))

heatsig <- (padjsig * fcsig) %>%
  as_tibble %>%
  abs %>%
  rowSums %>%
  {. >= 1} %>%
  which %>%
  dplyr::slice(wholeDEG, .) %>%
  inner_join(kmeansRes)

heatsig$ID <- as.character(heatsig$ID)
rawC <- rldData %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  inner_join(heatsig %>% select(ID, cl))
## inner_join(kmeansRes) ## all transcripts

## scale counts
scaleC <- rawC %>%
  select(contains(c('O','Y'),ignore.case = F)) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID, cl)) %>%
  mutate(GeneID = ID)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## go term
interesGO <- list(MHC_II = c('GO:0002399', 'GO:0002503', 'GO:0002495', 'GO:0002504','GO:0019886'),
                  MHC_I = c('GO:0002474','GO:0019885','GO:0042590'),
                  Humoral_immune = c("GO:0019882",'GO:0019724', 'GO:0002449', 'GO:0006959', 'GO:0016064', 'GO:0009404',"GO:0042113","GO:0050853"),
                  T_cell_immune = c('GO:0002456','GO:0002709','GO:0002711','GO:0050863','GO:0050870','GO:0002286','GO:0002292'),
                  Cytotoxicity = c('GO:0001909','GO:0001913','GO:0001910','GO:0001912','GO:0001916','GO:0001914','GO:0042267','GO:0042269','GO:0070942'),
                  defense_virus = c('GO:0031349','GO:0009615', 'GO:0051607', 'GO:0006826', 'GO:0046916', 'GO:0015706', 'GO:0010167', 'GO:0000041', 'GO:0031668', 'GO:0098771'),
                  response_stimulus = c('GO:0140546','GO:0002831','GO:0032103','GO:0002833','GO:0071216','GO:0031349'),
                  'Inf_β_response' = c('GO:0035456', 'GO:0035458', 'GO:0032728','GO:0032608','GO:0032648'),
                  'Inf_γ_response' = c('GO:0034341','GO:0071346','GO:0060333','GO:0032609','GO:0060330','GO:0060334','GO:0032649'),
                  'Inf_α_response' = c('GO:0035455','GO:0035457','GO:0032647','GO:0032607','GO:0032727'),
                  Inf_I_response= c('GO:0034340','GO:0060337','GO:0071357','GO:0032606','GO:0032479','GO:0032481','GO:0060338'),
                  TNF_response = c('GO:0032680','GO:0032640','GO:0032760','GO:0034612','GO:0071356','GO:0033209'),
                  cytokine_response = c('GO:0019221','GO:0060759','GO:0001959','GO:0002718','GO:0001961','GO:0002367','GO:0060759','GO:0060760'),
                  TNF_cytokine = c('GO:1903555','GO:0071706','GO:1903557'),
                  Inflammatory_response = c('GO:0050727','GO:0050729','GO:0002526','GO:0002532','GO:0002437','GO:0002438','GO:1900015','GO:0002861'),
                  Innate_immune = c("GO:0045088",'GO:0045089','GO:0002218'),
                  Adaptive_immune = c('GO:0002819','GO:0002821'),
                  complement_activation = c('GO:0006958','GO:0045917','GO:0045958','GO:0045960'),
                  PPAR = c('GO:0035360'),
                  Complement = c('GO:0006958'))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##----------------------- 看cpBP富集的通路有是否指定的interesGO
# 没有interesGene结果，报错，可以手动指定i绘图
for (i in seq_along(interesGO)) {

  interesGene <- cpBP %>%
    dplyr::filter(ID %in% interesGO[[i]]) %>%
    .$geneID %>%
    strsplit(split = '/', fixed = TRUE) %>%
    unlist %>%
    unique
  #> head(interesGene)
  #[1] "16149" "14960" "14961" "14969" "14999" "15000"
  interesMat <- scaleC %>%
    dplyr::filter(GeneID %in% interesGene) %>%
    inner_join(anno)
  #> head(interesMat)
  ## A tibble: 6 × 18
  #    OV1   OV2   OV3      O1     O2     O3    YV1    YV2    YV3    YV4     Y1
  #  <dbl> <dbl> <dbl>   <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
  #1 0.637 0.775 0.739 -0.0450 -0.227 -0.197  0.683  0.820  0.799  0.861 -1.70
  #2 0.223 0.913 0.464  1.52    1.37   1.33  -0.700 -1.10  -0.935 -0.736 -0.900
  #3 0.210 0.731 0.390  1.54    1.43   1.44  -0.799 -1.04  -0.858 -0.804 -0.885
  #4 0.314 1.01  0.467  1.56    1.21   1.31  -0.670 -0.912 -0.779 -0.758 -1.06
  #5 0.572 0.878 0.924 -0.290  -0.199 -0.342  0.481  0.800  0.925  0.910 -1.53
  #6 0.241 0.922 0.525  1.48    1.36   1.31  -0.941 -1.10  -0.920 -0.870 -0.719
  ## … with 7 more variables: Y2 <dbl>, Y3 <dbl>, ID <chr>, cl <dbl>,
  ##   GeneID <chr>, Gene <chr>, Description <chr>
  ## ℹ Use `colnames()` to see all variable names

  matcol <- colorRamp2(seq(min(scaleC %>% select(contains(c('O','Y'),ignore.case = F))), max(scaleC %>% select(contains(c('O','Y'),ignore.case = F))), length = 100), colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -6, -7)])(100))

  dim(interesMat) %>% print
  #> dim(interesMat) %>% print
  #[1] 16 17

  # 判断是否还有数据
  if(0 %in% dim(interesMat)){
    print("No data remain!")
  }else{
    interesMat <- as.data.frame(interesMat)
    duplicated_rows <- duplicated(interesMat$Gene)
    print(interesMat$Gene[duplicated_rows])
    # 删掉重复行
  if(any(duplicated_rows) == TRUE){interesMat <- interesMat[!duplicated_rows, ]}
    rownames(interesMat) <- interesMat$Gene
  ht_list <- Heatmap(matrix = interesMat %>%
                       select(contains(c('O','Y'),ignore.case = F)) %>%
                       apply(1, meanFlg22) %>%
                       t,
                     name = 'Scaled Counts',
                     row_order = order(interesMat$cl) %>% rev,
                     row_split = interesMat$cl,
                     row_gap = unit(2, "mm"),
                     column_order = 1 : 4,
                     column_split = rep(c('OV', 'O', 'YV', 'Y'), c(1, 1, 1, 1)),
                     show_column_names = FALSE,
                     col = matcol,
                     use_raster = FALSE)

  pdf(paste0(savepath, '/', 'GO_', names(interesGO)[i], '_1stadd_', topGONum, '.pdf'))
  draw(ht_list)
  dev.off()
  jpeg(paste0(savepath, '/', 'GO_', names(interesGO)[i], '_1stadd_', topGONum, '.jpg'))
  draw(ht_list)
  dev.off()
  }
}
##-----------------------
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##----------------------- 看interesGO中指定的go term 所有基因的表达情况
# 没有interesGene结果，报错，可以手动指定i绘图
## 在这里指定
interesGO <- list(PPAR = c('GO:0035360'),
                  Complement = c('GO:0006958'))

for (i in seq_along(interesGO)) {

  GOgeneID <- get(interesGO[[i]], org.Mm.egGO2ALLEGS) %>% mget(org.Mm.egSYMBOL) %>% unlist() %>% data.frame() %>% rownames_to_column(.)
  colnames(GOgeneID) <- c("ID","symbol")

  interesMat <- scaleC %>% inner_join(anno) %>%
    dplyr::filter(Gene %in% GOgeneID$symbol)
  #> head(interesMat)
  ## A tibble: 6 × 18
  #    OV1   OV2   OV3      O1     O2     O3    YV1    YV2    YV3    YV4     Y1
  #  <dbl> <dbl> <dbl>   <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
  #1 0.637 0.775 0.739 -0.0450 -0.227 -0.197  0.683  0.820  0.799  0.861 -1.70
  #2 0.223 0.913 0.464  1.52    1.37   1.33  -0.700 -1.10  -0.935 -0.736 -0.900
  #3 0.210 0.731 0.390  1.54    1.43   1.44  -0.799 -1.04  -0.858 -0.804 -0.885
  #4 0.314 1.01  0.467  1.56    1.21   1.31  -0.670 -0.912 -0.779 -0.758 -1.06
  #5 0.572 0.878 0.924 -0.290  -0.199 -0.342  0.481  0.800  0.925  0.910 -1.53
  #6 0.241 0.922 0.525  1.48    1.36   1.31  -0.941 -1.10  -0.920 -0.870 -0.719
  ## … with 7 more variables: Y2 <dbl>, Y3 <dbl>, ID <chr>, cl <dbl>,
  ##   GeneID <chr>, Gene <chr>, Description <chr>
  ## ℹ Use `colnames()` to see all variable names

  matcol <- colorRamp2(seq(min(scaleC %>% select(contains(c('O','Y'),ignore.case = F))), max(scaleC %>% select(contains(c('O','Y'),ignore.case = F))), length = 100), colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -6, -7)])(100))

  dim(interesMat) %>% print
  #> dim(interesMat) %>% print
  #[1] 16 17

  # 判断是否还有数据
  if(0 %in% dim(interesMat)){
    print("No data remain!")
  }else{
    interesMat <- as.data.frame(interesMat)
    duplicated_rows <- duplicated(interesMat$Gene)
    print(interesMat$Gene[duplicated_rows])
    # 删掉重复行
  if(any(duplicated_rows) == TRUE){interesMat <- interesMat[!duplicated_rows, ]}
    rownames(interesMat) <- interesMat$Gene
  ht_list <- Heatmap(matrix = interesMat %>%
                       select(contains(c('O','Y'),ignore.case = F)) %>%
                       apply(1, meanFlg22) %>%
                       t,
                     name = 'Scaled Counts',
                     row_order = order(interesMat$cl) %>% rev,
                     row_split = interesMat$cl,
                     row_gap = unit(2, "mm"),
                     column_order = 1 : 4,
                     column_split = rep(c('OV', 'O', 'YV', 'Y'), c(1, 1, 1, 1)),
                     show_column_names = FALSE,
                     col = matcol,
                     use_raster = FALSE)

  pdf(paste0(savepath, '/', 'GO_', names(interesGO)[i], '_1stadd_', 'all_gene_.pdf'))
  draw(ht_list)
  dev.off()
  jpeg(paste0(savepath, '/', 'GO_', names(interesGO)[i], '_1stadd_', 'all_gene_.jpg'))
  draw(ht_list)
  dev.off()
  }
}

## 指定基因列表，colnames = symbol
MHC_GENE <- read.table("/data/home/mali/genome/mouse/geneset/MHC_GENE_MOUSE.txt",header = T)
colnames(MHC_GENE) <- 'symbol'
Complement_and_coagulation <- read.table("/data/home/mali/genome/mouse/geneset/Complement_and_coagulation.txt",header = T)
colnames(Complement_and_coagulation) <- 'symbol'

gene_set <- Complement_and_coagulation
save_name <- "Complement_and_coagulation"

interesMat <- scaleC %>% inner_join(anno) %>%
    dplyr::filter(Gene %in% gene_set$symbol)

matcol <- colorRamp2(seq(min(scaleC %>% select(contains(c('O','Y'),ignore.case = F))), max(scaleC %>% select(contains(c('O','Y'),ignore.case = F))), length = 100), colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -6, -7)])(100))

dim(interesMat) %>% print
  #> dim(interesMat) %>% print
  #[1] 16 17

  # 判断是否还有数据
  if(0 %in% dim(interesMat)){
    print("No data remain!")
  }else{
    interesMat <- as.data.frame(interesMat)
    duplicated_rows <- duplicated(interesMat$Gene)
    print(interesMat$Gene[duplicated_rows])
    # 删掉重复行
  if(any(duplicated_rows) == TRUE){interesMat <- interesMat[!duplicated_rows, ]}
    rownames(interesMat) <- interesMat$Gene
  ht_list <- Heatmap(matrix = interesMat %>%
                       select(contains(c('O','Y'),ignore.case = F)) %>%
                       apply(1, meanFlg22) %>%
                       t,
                     name = 'Scaled Counts',
                     row_order = order(interesMat$cl) %>% rev,
                     row_split = interesMat$cl,
                     row_gap = unit(2, "mm"),
                     column_order = 1 : 4,
                     column_split = rep(c('OV', 'O', 'YV', 'Y'), c(1, 1, 1, 1)),
                     show_column_names = FALSE,
                     col = matcol,
                     use_raster = FALSE)

  pdf(paste0(savepath, '/', 'GO_',save_name,'_gene_.pdf'))
  draw(ht_list)
  dev.off()
  jpeg(paste0(savepath, '/', 'GO_',save_name,'gene_.jpg'))
  draw(ht_list)
  dev.off()
  }

#######################################################################
