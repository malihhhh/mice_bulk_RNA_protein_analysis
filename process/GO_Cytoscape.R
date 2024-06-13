
## originally by Yulong Niu
## yulong.niu@hotmail.com

# 产生了节点文件，应该还需要去cytoscape中绘图，这里没有产生图的代码
##########################cytoscape GO plot###########################
library('tidyverse')
library('reshape2')
library('magrittr')
library('foreach')
library('dplyr')
library('ggplot2')
library('enrichplot')
useVersion   <- "version1"
outputFolder <- file.path("/data/home/mali/ke_data/output", useVersion, "Cytoscape")
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)


##~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comb2_internal <- function(x) {

  l <- length(x)

  if (l < 2) {
    m <- matrix(0, nrow = 0, ncol = 2)
  } else {
    fV <- rep(x[1:(l - 1)], ((l - 1):1))
    tVIdx <- unlist(mapply(seq, 2:l, l))
    tV <- x[tVIdx]

    m <- cbind(fV, tV)
  }

  return(m)
}

combWhole_internal <- function(x, y, self = FALSE, bidirect = FALSE) {

  ## check x
  ## first unique x
  x <- unique(x)
  x <- x[x %in% y]

  l <- length(x)

  if (l >= length(y)) {
    m <- comb2_internal(x)
  } else {
    yLeft <- y[!(y %in% x)]
    m <- cbind(rep(x, each = length(yLeft)),
               rep(yLeft, l))
    m <- rbind(m,
               comb2_internal(x))
  }

  if (bidirect) {
    m <- rbind(m, m[, 2:1])
  } else {}

  if (self) {
    m <- rbind(m, combSelf_internal(x))
  } else {}

  colnames(m) <- NULL

  return(m)
}

JacSim <- function (x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x, y)))
}

Shrinkage <- function(x, minNew, maxNew) {
  ## INPUT: `x` is numeric vector. `minNew` and `maxNew` is the new min and new max.
  ## OUTPUT: A shrink vector.

  unit <- (maxNew - minNew) / (max(x) - min(x))

  res <- minNew + (x - min(x)) * unit

  return(res)
}

GOCytoEdge <- function(cpRes, JacSimThres = 0.2) {
  ## INTPUT: `cpRes` is the clusterProfiler table. `JacSimThres` is the threshold of Jaccard simialrity.
  ## OUTPUT: A tibble of edge matrix.
  ## USAGE: Generate the edge table for Cytoscape.

  require('tidyverse')

  ## step1: remove duplicated terms
  cpRes %<>%
    dplyr::select(ID, geneID, Description) %>%
    group_by(ID) %>%
    summarise(geneID = paste(geneID, collapse = '/'),
              Description = sample(Description, 1)) %>%
    ungroup

  termNum <- nrow(cpRes)

  ## step2: find intersection
  cpResList <- strsplit(cpRes$geneID, split = '/', fixed = TRUE) %>%
    sapply(unique) %>%
    setNames(cpRes$ID)

  interMat <- combWhole_internal(1 : termNum, 1:termNum) %>%
    set_colnames(c('from', 'to')) %>%
    as.data.frame %>%
    as_tibble

  interMat %<>%
    mutate(jacSim = apply(., 1, function(x) {
      eachJacSim <- JacSim(cpResList[[x[1]]], cpResList[[x[2]]])
      return(eachJacSim)
    })) %>%
    filter(jacSim >= JacSimThres) %>% ## filter by jaccard similarity
    mutate(edgeWidth = Shrinkage(jacSim, 2, 5)) %>% ## edge width
    mutate(SOURCE = cpRes$ID[from], TARGET = cpRes$ID[to]) %>%
    mutate(fromDesc = cpRes$Description[from], toDesc = cpRes$Description[to]) %>%
    dplyr::select(-from, -to)

  return(interMat)
}

GOCytoNode <- function(cpRes) {
  ## INTPUT: `cpRes` is the clusterProfiler table.
  ## OUTPUT: A tibble of node matrix.
  ## USAGE: Generate the node table for Cytoscape.

  require('reshape2')

  ## unique GO terms
  cpUniq <- cpRes %>%
    dplyr::select(ID, geneID, geneName, Description) %>%
    group_by(ID) %>%
    summarise(geneID = paste(geneID, collapse = '/'),
              Description = sample(Description, 1),
              geneID = sample(geneID, 1),
              geneName = sample(geneName, 1)) %>%
    ungroup

  nodeMat <- cpRes %>%
    mutate(Cluster = str_extract(Cluster, 'cluster\\d+')) %>%
    dcast(ID ~ Cluster, value.var = 'Count') %>%
    mutate_all(~ifelse(is.na(.), 0, .)) %>%
    mutate(nodeSize = dplyr::select(., -ID) %>% rowSums %>% Shrinkage(30, 50)) %>%
    inner_join(cpUniq) %>%
    dplyr::rename(SOURCE=ID)

  return(nodeMat)

}

CytoGeneEdge_ <- function(GOGene, geneAnno) {

  ## INPUT: `GOgene` is the 'GO-gene' matrix. `geneAnno` is the gene annotation.
  ## OUTPUT: A `tibble` of interaction matrix.

  res <- tibble(SOURCE = strsplit(GOGene$geneName, split = '/', fixed = TRUE) %>% unlist,
                TARGET = GOGene$ID,
                toDesc = GOGene$Description) %>%
    mutate(jacSim = 1,
           edgeWidth = 1) %>%
    mutate(fromDesc = anno$Description[match(SOURCE, anno$Gene)]) %>%
    dplyr::select(jacSim, edgeWidth, SOURCE, TARGET, fromDesc, toDesc)

  return(res)
}


GOCytoGeneEdge <- function(cpRes, ...) {
  ## INPUT: `cpRes` is the clusterProfiler table.
  ## OUTPUT: A tibble of gene-GO matrix.
  ## USAGE: Generate the node table for Cytoscape.

  require('foreach')

  res <- foreach(i = seq_len(nrow(cpRes)), .combine = bind_rows) %do% {
    CytoGeneEdge_(cpRes[i, ], ...)
  }

  return(res)
}

CytoGeneNode_ <- function(GOGene, anno) {
  ## INPUT: `GOgene` is the 'GO-gene' matrix. `geneAnno` is the gene annotation.
  ## OUTPUT: A `tibble` of node matrix.

  GOGene %<>%
    mutate(Cluster = str_extract(Cluster, 'cluster\\d+'))

  res <- tibble(SOURCE = strsplit(GOGene$geneName, split = '/', fixed = TRUE) %>% unlist,
                Cluster = GOGene$Cluster) %>%
    mutate(nodeSize = 10,
      geneID = strsplit(GOGene$geneID, split = '/', fixed = TRUE) %>% unlist,
      Description = anno$Description[match(SOURCE, anno$Gene)],
      geneName = SOURCE)

  return(res)
}


GOCytoGeneNode <- function(cpRes, ...) {
  ## INPUT: `cpRes` is the clusterProfiler table.
  ## OUTPUT: A tibble of gene-GO edge.
  ## USAGE: Generate the node table for Cytoscape.

  require('foreach')
  require('reshape2')

  nodeMat <- foreach(i = seq_len(nrow(cpRes)), .combine = bind_rows) %do% {
    print(i)
    #CytoGeneNode_(cpRes[i, ], ...)
    CytoGeneNode_(cpRes[i, ], anno)

  } %>%
    dplyr::slice(which(!duplicated(.)))

  clusterMat <- nodeMat %>%
    mutate(Val = 1) %>%
    dcast(SOURCE ~ Cluster, value.var = 'Val', fill = 0) %>%
    mutate_at(.vars = vars(contains('cluster')),
              .funs = ~ifelse(. > 0, 1, .))

  res <- inner_join(clusterMat, nodeMat %>% dplyr::select(-Cluster)) %>%
    as_tibble

  return(res)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# kallGOBP
load('/data/home/mali/ke_data/transcription/count/kmeans10_cp_BP.RData')
load('/data/home/mali/ke_data/transcription/count/kmeans10_cp_BP_sim.RData')

anno <- read_csv('/data/home/mali/genome/mouse/Ensembl_mmu_Anno_AddGeneName.csv',
                    col_types = cols(Chromosome = col_character())) %>%
  mutate_all(list(~replace(., is.na(.), ''))) %>%
  mutate(GeneID = entrezgene_id) %>%  # 因为要跟下面的数据对应上，下面用的是entrezgene_id
  dplyr::select(GeneID, Gene, Description) %>%
  dplyr::slice(which(!duplicated(.)))


#  > head(anno)
#  # A tibble: 6 × 3
#    GeneID Gene  Description
#    <chr>  <chr> <chr>
#  1 14679  Gnai3 guanine nucleotide binding protein (G protein), alpha inhibiting…
#  2 54192  Pbsn  probasin [Source:MGI Symbol;Acc:MGI:1860484]
#  3 15417  Hoxb9 homeobox B9 [Source:MGI Symbol;Acc:MGI:96190]
#  4 12544  Cdc45 cell division cycle 45 [Source:MGI Symbol;Acc:MGI:1338073]
#  5 16002  Igf2  insulin-like growth factor 2 [Source:MGI Symbol;Acc:MGI:96434]
#  6 11818  Apoh  apolipoprotein H [Source:MGI Symbol;Acc:MGI:88058]

cpBP <- "fortify"(kallGOBP,showCategory = 10000) %>%
  as_tibble %>%
  dplyr::mutate(.,geneName = sapply(geneID, function(x) {
    strsplit(x, split = '/', fixed = TRUE) %>%
      unlist %>%
      tibble(GeneID = .) %>%
      left_join(anno) %>%
      .$Gene %>%
        paste(collapse = '/')
}))

cpBPSim <- "fortify"(kallGOBPSim,showCategory = 10000) %>%
  as_tibble %>%
  dplyr::mutate(.,geneName = sapply(geneID, function(x) {
    strsplit(x, split = '/', fixed = TRUE) %>%
      unlist %>%
      tibble(GeneID = .) %>%
      left_join(anno) %>%
      .$Gene %>%
        paste(collapse = '/')
}))

# > head(cpBP)
# # A tibble: 6 × 11
#   Cluster  ID    Description GeneRatio BgRatio   pvalue p.adjust   qvalue geneID
#   <fct>    <chr> <fct>           <dbl> <chr>      <dbl>    <dbl>    <dbl> <chr>
# 1 "cluste… GO:0… extracellu…    0.177  322/28… 2.18e-13 1.43e-10 1.03e-10 76293…
# 2 "cluste… GO:0… extracellu…    0.177  323/28… 2.27e-13 1.43e-10 1.03e-10 76293…
# 3 "cluste… GO:0… external e…    0.177  324/28… 2.37e-13 1.43e-10 1.03e-10 76293…
# 4 "cluste… GO:0… extracellu…    0.0633 44/289… 1.33e- 7 6.00e- 5 4.35e- 5 76293…
# 5 "cluste… GO:1… regulation…    0.0380 10/289… 2.32e- 6 8.36e- 4 6.06e- 4 12389…
# 6 "cluste… GO:0… toxin meta…    0.0380 12/289… 4.23e- 6 1.09e- 3 7.90e- 4 55990…
# # ℹ 2 more variables: Count <int>, geneName <named list>
# > dim(cpBP)
# [1] 1615   11
# > dim(cpBPSim)
# [1] 1557   11

# 这里对cpBP进行一下筛选，只保留每个cluster最具代表性的前10个goterm
new_cpBP <- data.frame()
for (type in unique(cpBP$Cluster)) {
  type_rows <- cpBP[cpBP$Cluster == type, ]
  type_rows <- type_rows[1:min(nrow(type_rows), 10), ]
  new_cpBP <- rbind(new_cpBP, type_rows)
}
new_cpBP <- tibble(new_cpBP)

new_cpBPSim <- data.frame()
for (type in unique(cpBPSim$Cluster)) {
  type_rows <- cpBPSim[cpBPSim$Cluster == type, ]
  type_rows <- type_rows[1:min(nrow(type_rows), 10), ]
  new_cpBPSim <- rbind(new_cpBPSim, type_rows)
}
new_cpBPSim <- tibble(new_cpBPSim)

##########################################
name <- "cpBP"
##########################################
## GO network
GOCytoEdge(cpBP) %>% write_csv(file.path(outputFolder, paste0(name,'_GOCytoEdge.csv')))
GOCytoNode(cpBP) %>% write_csv(file.path(outputFolder,paste0(name,'_GOCytoNode.csv')))
cpBP %>% write_csv(file.path(outputFolder,paste0(name,'_original.csv')))

## GO-gene network
bind_rows(GOCytoEdge(cpBP),
          GOCytoGeneEdge(cpBP, anno)) %>% write_csv(file.path(outputFolder,'tmp4.csv'))

bind_rows(GOCytoNode(cpBP),
          GOCytoGeneNode(cpBP, anno)) %>% write_csv(file.path(outputFolder,'tmp5.csv'))

##########################################
name <- "cpBPSim"
##########################################
## GO network
GOCytoEdge(cpBPSim) %>% write_csv(file.path(outputFolder, paste0(name,'_GOCytoEdge.csv')))
GOCytoNode(cpBPSim) %>% write_csv(file.path(outputFolder,paste0(name,'_GOCytoNode.csv')))
cpBPSim %>% write_csv(file.path(outputFolder,paste0(name,'_original.csv')))

## GO-gene network
bind_rows(GOCytoEdge(cpBPSim),
          GOCytoGeneEdge(cpBPSim, anno)) %>% write_csv(file.path(outputFolder,paste0(name,'_GOCytoGeneEdge.csv')))

bind_rows(GOCytoNode(cpBPSim),
          GOCytoGeneNode(cpBPSim, anno)) %>% write_csv(file.path(outputFolder,paste0(name,'_GOCytoGeneNode.csv')))

##########################################
name <- "new_cpBP"
##########################################
## GO network
GOCytoEdge(new_cpBP) %>% write_csv(file.path(outputFolder, paste0(name,'_GOCytoEdge.csv')))
GOCytoNode(new_cpBP) %>% write_csv(file.path(outputFolder,paste0(name,'_GOCytoNode.csv')))
new_cpBP %>% write_csv(file.path(outputFolder,paste0(name,'_original.csv')))

## GO-gene network
bind_rows(GOCytoEdge(new_cpBP),
          GOCytoGeneEdge(new_cpBP, anno)) %>% write_csv(file.path(outputFolder,paste0(name,'_GOCytoGeneEdge.csv')))

bind_rows(GOCytoNode(new_cpBP),
          GOCytoGeneNode(new_cpBP, anno)) %>% write_csv(file.path(outputFolder,paste0(name,'_GOCytoGeneNode.csv')))

##########################################
name <- "new_cpBPSim"
##########################################
## GO network
GOCytoEdge(new_cpBPSim) %>% write_csv(file.path(outputFolder, paste0(name,'_GOCytoEdge.csv')))
GOCytoNode(new_cpBPSim) %>% write_csv(file.path(outputFolder,paste0(name,'_GOCytoNode.csv')))
new_cpBPSim %>% write_csv(file.path(outputFolder,paste0(name,'_original.csv')))

## GO-gene network
bind_rows(GOCytoEdge(new_cpBPSim),
          GOCytoGeneEdge(new_cpBPSim, anno)) %>% write_csv(file.path(outputFolder,paste0(name,'_GOCytoGeneEdge.csv')))

bind_rows(GOCytoNode(new_cpBPSim),
          GOCytoGeneNode(new_cpBPSim, anno)) %>% write_csv(file.path(outputFolder,paste0(name,'_GOCytoGeneNode.csv')))



##~~~~~~~~~~~~~~~~interesting genes cluster~~~~~~~~~~~~~~~~~~~~~~~~~~
## AT2G19190 -- SIRK/FRK1 -- 3
## AT1G14550 -- PER5 -- 3
## AT1G18570 -- MYB51 -- 5
## AT2G30750 -- CYP71A12 -- 5
## AT1G19250 -- FMO1 -- 5
## AT1G73805 -- SARD1 -- 1
## AT4G28460 -- PIP1 -- 5
## AT4G37290 -- PIP2 -- 5
## AT3G48090 -- EDS1 -- 1
## AT3G07040 -- RPM1 -- 1
## AT3G50950 -- ZAR1/RPP13L4 -- 1
## AT4G11170 -- -- 5
## AT5G41540 -- -- 1
interesGene <- c('Il6',
                 'AT1G14550',
                 'AT1G18570',
                 'AT2G30750',
                 'AT1G73805',
                 'AT4G28460',
                 'AT4G37290',
                 'AT3G48090',
                 'AT3G07040',
                 'AT3G50950',
                 'AT4G11170',
                 'AT5G41540')

gene_info <- str_detect(GOCytoNode(cpBP)$geneName, interesGene %>% paste(collapse = '|')) %>%
  GOCytoNode(cpBP)[., c('SOURCE', 'Description')]

#> head(gene_info,2)
#      SOURCE         Description
#2 GO:0000165        MAPK cascade
#3 GO:0001666 response to hypoxia

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####################################################################
