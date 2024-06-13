if (!requireNamespace("readr", quietly = TRUE)) install.packages('readr')
library(data.table)
library(readr)
# -------------------------------------------------------------------------
# Installation
# -------------------------------------------------------------------------
if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}
library(RCy3)

if(!"RColorBrewer" %in% installed.packages()){
  install.packages("RColorBrewer")
}
library(RColorBrewer)
# The whole point of RCy3 is to connect with Cytoscape. You will need to
# install and launch Cytoscape. Then "ping" Cytoscape to verify:
## ------------------------------------------------------------------------
cytoscapePing()


# -------------------------------------------------------------------------
# Getting Networks
# -------------------------------------------------------------------------
# From dataframes
name <- "new_cpBPSim"
name <- "cpBPSim"
name <- "cpBP"
name <- "new_cpBP"
# -------------------------------------------------------------------------


outputFolder <- "E:\\ke_data\\Cytoscape"
nodes_data <- read_csv(file.path(outputFolder,paste0(name,'_GOCytoGeneNode.csv')))
edges_data <- read_csv(file.path(outputFolder,paste0(name,'_GOCytoGeneEdge.csv')))
cluster_data <- read_csv(file.path(outputFolder,paste0(name,'_original.csv')))
nodes_data <- na.omit(nodes_data)
edges_data <- na.omit(edges_data)
cluster_data <- na.omit(cluster_data)
cluster_data$logP <- -log10(cluster_data$p.adjust)
range(cluster_data$logP)
dim(edges_data)
# [1] 1496    6
dim(nodes_data)
# [1] 461  14
dim(cluster_data)
# [1] 90 11
head(cluster_data)
head(nodes_data)

setnames(nodes_data, old=c ("SOURCE"), new=c ("id"))
setnames(edges_data, old=c ("SOURCE","TARGET"), new=c ("source","target"))
createNetworkFromDataFrames(nodes_data, edges_data,title =name)

getLayoutPropertyNames(layout.name='force-directed')
layoutNetwork(paste('force-directed',
                    'defaultSpringCoefficient=0.00001',
                    'defaultSpringLength=50',
                    sep=' '))
# 给节点添加信息
?loadTableData
loadTableData(cluster_data,  
              data.key.column = "ID",
              table = "node",
              table.key.column = "name")  
#default data.frame key is row.names

# 之前配置好的话直接用
style.name = "myStyle"
setVisualStyle(style.name)

# 配置样式
createVisualStyle(style.name)
setVisualStyle(style.name)

getTableColumnNames('node')
setNodeSizeMapping("nodeSize", style.name=style.name)
setNodeShapeDefault("ellipse", style.name=style.name) 
# 查看RColorBrewer的颜色
display.brewer.all(type = "div")
display.brewer.all(type = "qual")
setNodeCustomPieChart(c("cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8","cluster9"),colors = c(brewer.pal(8,"Dark2"), "#1F78B4" ), style.name=style.name)

setNodeFillOpacityDefault(255, style.name=style.name)

setNodeBorderColorMapping("logP",c(0,20),colors=c("#999999", '#E41A1C'),style.name=style.name)
setNodeBorderWidthMapping("logP",widths=c(2,20),style.name=style.name)
setNodeBorderWidthBypass(c("#999999", '#E41A1C'), 20)
setNodeBorderColorBypass(start, "#999999")
setNodeBorderColorBypass(finish, '#E41A1C')


setNodeLabelDefault('',style.name=style.name)
setNodeLabelMapping('name',style.name=style.name)

setEdgeLineWidthDefault(2, style.name=style.name)
#setNodeLabelMapping('edge','edgeWidth', style.name=style.name)
setEdgeLineWidthMapping("edgeWidth", widths=c(2,10),mapping.type = "c", default.width=2,style.name=style.name)

nodeSize <- getTableColumns('node','nodeSize')
nodeSize

'''
# Filter for neighbors of high degree nodes
createDegreeFilter(filter.name = "logP",
                   criterion = c(0,15),
                   predicate = "IS_NOT_BETWEEN")
selectFirstNeighbors() # expand selection to first neighbors
createSubnetwork(subnetwork.name = "first neighbors of high degree nodes")

# Filter for high edge betweenness
createColumnFilter(filter.name = "edge betweenness",
                   type = "edges",
                   column = "EdgeBetweenness",
                   4000,
                   "GREATER_THAN",
                   network = net.suid)
createSubnetwork(subnetwork.name = "high edge betweenness")
'''

########################################### Saving and Exporting Networks
# Saving sessions
saveSession(file.path(outputFolder,name)) #.cys
## Leave filename blank to update previously saved session file

# Exporting images and networks
exportNetwork(filename =file.path(outputFolder,name)) #.sif
## Optionally specify filename, default is network name
## Optionally specify type: SIF(default), CX, cyjs, graphML, NNF, SIF, xGMML
exportImage(filename = file.path(outputFolder,name),type='pdf',height=20,width=20) #.png
## Optionally specify filename, default is network name
## Optionally specify type: PNG (default), JPEG, PDF, PostScript, SVG 


################################################################################
# 绘制cluster相关的热图
################################################################################


library(data.table)
library(readr)
library(RCy3)
library(RColorBrewer)
library(dplyr)

# 加载cytoscape绘图的数据
openSession("E:\\ke_data\\Cytoscape\\new_cpBP.cys")
outputFolder <- "E:\\ke_data\\Cytoscape"
name <- "new_cpBP"
nodes_data <- read_csv(file.path(outputFolder,paste0(name,'_GOCytoGeneNode.csv')))
nodes_data <- na.omit(nodes_data)
head(nodes_data)
node_cluster_all <- list()
node_cluster_gene <- list()
node_cluster_go <- list()
for (i in 1:9){
  n <- paste0("cluster",i)
  node_cluster_all[[i]] <- nodes_data[which(nodes_data[[n]] != 0),]
  node_cluster_go[[i]] <- node_cluster_all[[i]] %>% filter(grepl('GO:', SOURCE))
  node_cluster_gene[[i]] <- node_cluster_all[[i]] %>% filter(!grepl('GO:', SOURCE))
}
node_cluster_all
node_cluster_gene
node_cluster_go
saveRDS(node_cluster_all, file.path(outputFolder,paste0(name,'_node_cluster_all.rds')))
saveRDS(node_cluster_gene, file.path(outputFolder,paste0(name,'_node_cluster_gene.rds')))
saveRDS(node_cluster_go, file.path(outputFolder,paste0(name,'_node_cluster_go.rds')))
# 将这3个文件转移到服务器/data/home/mali/ke_data/output/version1/Cytoscape/
node_cluster_go[[9]][,c(1,11,13)]
node_cluster_go[[2]][,c(1,11,13)]
node_cluster_go[[1]][,c(1,11,13)]
node_cluster_go[[5]][,c(1,11,13)]
node_cluster_go[[4]][,c(1,11,13)]
node_cluster_go[[6]][,c(1,11,13)]
node_cluster_go[[7]][,c(1,11,13)]
node_cluster_go[[3]][,c(1,11,13)]
node_cluster_go[[8]][,c(1,11,13)]
################################################################################
# 在服务器上运行
################################################################################
library(dplyr)
library(RColorBrewer)
library(circlize)

useVersion   <- "version1"
outputFolder <- file.path("/data/home/mali/ke_data/output", useVersion, "Cytoscape")
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)

name <- "new_cpBP"
node_cluster_all <- readRDS(file.path(outputFolder,paste0(name,'_node_cluster_all.rds')))
node_cluster_gene <- readRDS(file.path(outputFolder,paste0(name,'_node_cluster_gene.rds')))
node_cluster_go <- readRDS(file.path(outputFolder,paste0(name,'_node_cluster_go.rds')))

# 加载在基因表达矩阵
load('/data/home/mali/ke_data/transcription/RDS/degres_condi_Mock.RData')
wholeDEG <- read_csv('/data/home/mali/ke_data/transcription/count/eachGroup_vs_Mock_k.csv')
kmeansRes <- read_csv('/data/home/mali/ke_data/transcription/count/kmeans10.csv') %>%
  dplyr::select(ID, cl)

fcsig <- wholeDEG %>%
  dplyr::select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > log2(1.5) ~ 1,
                                 . < -log2(1.5) ~ -1,
                                 TRUE ~ 0)))
padjsig <- wholeDEG %>%
  dplyr::select(ends_with('padj')) %>%
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
  inner_join(heatsig %>% dplyr::select(ID, cl))
## inner_join(kmeansRes) ## all transcripts

## scale counts
scaleC <- rawC %>%
  dplyr::select(contains(c('O','Y'),ignore.case = F)) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% dplyr::select(ID, cl)) %>%
  mutate(GeneID = ID)


meanFlg22 <- function(v) {
  
  require('magrittr')
  
  res <- v %>%
    split(c(rep(c(1:2), each = 3),rep(3, each = 4),rep(4, each = 3))) %>%
    sapply(mean, na.rm = TRUE)
  
  return(res)
}


## 指定基因列表，colnames = symbol
for (i in 1:9){
  gene <- node_cluster_gene[[i]][1]
  colnames(gene) <- 'symbol'
  save_name <- paste0("cluster",i)
  interesMat <- scaleC %>% inner_join(anno) %>% dplyr::filter(Gene %in% gene$symbol)
  matcol <- colorRamp2(seq(min(scaleC %>% dplyr::select(contains(c('O','Y'),ignore.case = F))), max(scaleC %>% dplyr::select(contains(c('O','Y'),ignore.case = F))), length = 100), colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -6, -7)])(100))
  dim(interesMat) %>% print
  
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
                         dplyr::select(contains(c('O','Y'),ignore.case = F)) %>%
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
    
    pdf(paste0(outputFolder, '/', 'GO_',save_name,'_gene_.pdf'))
    draw(ht_list)
    dev.off()
    jpeg(paste0(outputFolder, '/', 'GO_',save_name,'_gene_.jpg'))
    draw(ht_list)
    dev.off()
  }
}

MHC_GENE <- read.table("/data/home/mali/genome/mouse/geneset/MHC_GENE_MOUSE.txt",header = T)
colnames(MHC_GENE) <- 'symbol'
Complement_and_coagulation <- read.table("/data/home/mali/genome/mouse/geneset/Complement_and_coagulation.txt",header = T)
colnames(Complement_and_coagulation) <- 'symbol'

gene_set <- Complement_and_coagulation
save_name <- "Complement_and_coagulation"

interesMat <- scaleC %>% inner_join(anno) %>%
  dplyr::filter(Gene %in% gene_set$symbol)

matcol <- colorRamp2(seq(min(scaleC %>% dplyr::select(contains(c('O','Y'),ignore.case = F))), max(scaleC %>% dplyr::select(contains(c('O','Y'),ignore.case = F))), length = 100), colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -6, -7)])(100))

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
                       dplyr::select(contains(c('O','Y'),ignore.case = F)) %>%
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
