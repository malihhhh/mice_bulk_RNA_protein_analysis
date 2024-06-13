
## originally by Yulong Niu
## yulong.niu@hotmail.com

######################hierarchical clustering####################
setwd('/data/home/mali/ke_data/')
# 我自己演示得到的结果保存在了~/plant/MPIPZ_microbe-host_homeostasis-master/results/removeZero/end
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

useVersion   <- "version1"
outputFolder <- file.path("/data/home/mali/ke_data/output", useVersion, "cluster_RNA")
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)

load('/data/home/mali/ke_data/transcription/RDS/degres_condi_Mock.RData')
deganno <- read_csv('/data/home/mali/ke_data/transcription/count/eachGroup_vs_Mock_k.csv',
                   col_types = cols(Chromosome = col_character()))
#> dim(deganno)
#[1] 14139    34
##~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meanFlg22 <- function(v) {

  require('magrittr')

  res <- v %>%
    split(c(rep(1, each = 4),rep(2 : 4, each = 3))) %>%
    sapply(mean, na.rm = TRUE)

  return(res)
}

##p value calculation from WGCNA
corPvalueStudent <- function(cor, nSamples) {

  ## ref: https://courses.lumenlearning.com/introstats1/chapter/testing-the-significance-of-the-correlation-coefficient/
  T <- sqrt(nSamples - 2) * cor / sqrt(1 - cor^2)

  p <- apply(T, 1:2, function(x) {
    if (x < 0) {
      eachp <- 1 -  pt(x, nSamples - 2, lower.tail = FALSE)
    } else {
      eachp <- pt(x, nSamples - 2, lower.tail = FALSE)
    }
  })

  return(p)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~prepare counts~~~~~~~~~~~~~~~~~~~~~~~~~~
rawCount <- rldData
dim(rldData)
# [1] 12150    13
## mean value of normalized count
sampleN <- c('YV', 'OV', 'O', 'Y')
meanCount <- rawCount %>%
  apply(1, meanFlg22) %>%
  t
colnames(meanCount) <- sampleN
head(meanCount)
#> head(meanCount)
#                 YV       OV        O        Y
#20389     13.112135 13.10988 13.57573 13.77995
#22287     12.343966 13.41149 13.50753 13.15494
#100503605 12.602466 11.88846 12.25138 13.07212
#15122     12.349530 11.58976 12.05158 12.88871
#16071      9.915731 11.99895 12.12904 10.13740
#17105     11.518689 11.84820 12.08211 11.94058
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~cluster~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## scale
scaleCount <- meanCount %>%
  t %>%
  scale %>%
  t
scaleCount %<>% .[complete.cases(.), ]
head(scaleCount)
#head(scaleCount)
#                  YV           OV          O           Y
#20389     -0.8358732 -0.842563055  0.5368694  1.14156681
#22287     -1.4392678  0.581008299  0.7627709  0.09548855
#100503605  0.2947868 -1.119151757 -0.4004710  1.22483595
#15122      0.2380712 -1.157224789 -0.3091010  1.22825463
#16071     -0.9564479  0.807521998  0.9176789 -0.76875292
#17105     -1.3734883  0.003354552  0.9807741  0.38935964
dim(scaleCount)
#[1] 12147     4

## Cluster rows by Pearson correlation
hr <- scaleCount %>%
  t %>%
  cor(method = 'pearson') %>%
  {1 - .} %>%
  as.dist %>%
  hclust(method = 'complete')

## Clusters columns by Spearman correlation
hc <- scaleCount %>%
  cor(method = 'spearman') %>%
  {1 - .} %>%
  as.dist %>%
  hclust(method = 'complete')

cairo_pdf(file.path(outputFolder,'heatmap_script_cluster.pdf'))
heatmap.2(meanCount,
          Rowv = as.dendrogram(hr),
          Colv = as.dendrogram(hc),
          col = redgreen(100),
          scale = 'row',
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          main = 'Heatmap.2',
          trace = 'none')
dev.off()

hc %>%
  as.dendrogram(method = 'average') %>%
  plot(main = 'Sample Clustering',
       ylab = 'Height')

hr %>%
  as.dendrogram(method = 'average') %>%
  plot(leaflab = 'none',
       main = 'Gene Clustering',
       ylab = 'Height')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~dynamic tree cut~~~~~~~~~~~~~~~~~~~~
clusDyn <- scaleCount %>%
  t %>%
  cor %>%
  {1 - .} %>%
  as.dist %T>%
  gc %>%
  as.matrix %T>%
  gc %>%
  cutreeDynamic(hr, distM = ., method = 'hybrid')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~K-means cluster~~~~~~~~~~~~~~~~~~~~~~
z_var <- apply(meanCount, 1, var)
z_mean <- apply(meanCount, 1, mean)
plot(log2(z_mean), log2(z_var), pch = '.')
abline(h = 1, col='red')
abline(v = 1, col='red')
text(x = 13,
     y = 23,
     labels = 'variance > 1 &\n mean > 1',
     col = 'red')

## filter
## meanCount %<>% .[which(z_var > 0 & z_mean > 0), ]

## choose groups
## 1. sum of squared error
wss <- (nrow(scaleCount) - 1) * sum(apply(scaleCount, 2, var))

for (i in 2:20) {
  wss[i] <- sum(kmeans(scaleCount,
                       centers=i,
                       algorithm = 'MacQueen')$withinss)
}

ggplot(tibble(k = 1:20, wss = wss), aes(k, wss)) +
  geom_point(colour = '#D55E00', size = 3) +
  geom_line(linetype = 'dashed') +
  xlab('Number of clusters') +
  ylab('Sum of squared error')
ggsave(file.path(outputFolder,'kmeans_sse.pdf'))
ggsave(file.path(outputFolder,'kmeans_sse.jpg'))

## 2. Akaike information criterion
kmeansAIC = function(fit){
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + 2*m*k)
}

aic <- numeric(20)
for (i in 1:20) {
  fit <- kmeans(x = scaleCount, centers = i, algorithm = 'MacQueen')
  aic[i] <- kmeansAIC(fit)
}

ggplot(tibble(k = 1:20, aic = aic), aes(k, wss)) +
  geom_point(colour = '#009E73', size = 3) +
  geom_line(linetype = 'dashed') +
  xlab('Number of clusters') +
  ylab('Akaike information criterion')
ggsave(file.path(outputFolder,'kmeans_AIC.pdf'))
ggsave(file.path(outputFolder,'kmeans_AIC.jpg'))

## execute
clNum <- 10
kClust10 <- kmeans(scaleCount, centers = clNum, algorithm= 'MacQueen', nstart = 1000, iter.max = 20)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~cut trees by height ~~~~~~~~~~~~~~~~~~~~~~~~~
hclusth2.0 <- cutree(hr, h = 2.0)
hclusth1.5 <- cutree(hr, h = 1.5)
hclusth1.0 <- cutree(hr, h = 1.0)
hclusth0.5 <- cutree(hr, h = 0.5)

cairo_pdf(file.path(outputFolder,'genetree.pdf'))
treeR <- hr %>%
  as.dendrogram(method = 'average')

plot(treeR,
     leaflab = 'none',
     main = 'Gene Clustering',
     ylab = 'Height')

cbind(hclusth0.5, hclusth1.0, hclusth1.5,hclusth2.0, clusDyn, kClust10$cluster) %>%
  colored_bars(treeR,
               sort_by_labels_order = TRUE,
               y_shift = -0.1,
               rowLabels = c('h=0.5','h=1.0','h=1.5','h=2.0', 'Dynamic', 'k-means(k=10)'),
               cex.rowLabels=0.7)
abline(h=2.0, lty = 2, col='grey')
abline(h=1.5, lty = 2, col='grey')
abline(h=1.0, lty = 2, col='grey')
abline(h=0.5, lty = 2, col='grey')
dev.off()

#cgenes <- c('AT1G14550.1', 'AT2G30750.1', 'AT2G19190.1', 'AT5G24110.1')
#hclusth1.5[cgenes]
#hclusth1.0[cgenes]
#hclusth0.5[cgenes]
#kClust10$cluster[cgenes]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~plot patterns~~~~~~~~~~~~~~~~~~~~~~~~
## join cluster and scaled normalized counts
kClust10$ID <- names(kClust10$cluster)
kClust10_used <- cbind.data.frame(ID = kClust10$ID,cl = kClust10$cluster)
#> head(kClust10_used)
#                 ID cluster
#20389         20389       8
#22287         22287       5
#100503605 100503605       6
#15122         15122       6
#16071         16071       1
#17105         17105       4
#> table(kClust10_used$cluster)
#   1    2    3    4    5    6    7    8    9   10
#1071 1224  493 1320 1201 1483 1424 1612 1677  642

deganno$ID <- as.character(deganno$ID)
deganno_add <- left_join(deganno,kClust10_used, by = 'ID')
write_csv(deganno_add, "/data/home/mali/ke_data/transcription/count/kmeans10.csv")



kmeansRes <- read_csv('/data/home/mali/ke_data/transcription/count/kmeans10.csv',
                      col_types = cols(Chromosome = col_character())) %>%
  dplyr::select(ID, cl) %>%
  dplyr::rename(clreal = cl)

#cl <- kmeansRes$clreal[match(names(kClust10$cluster), kmeansRes$ID)] %>%
#  set_names(names(kClust10$cluster))

cl <- kClust10$cluster
prefix <- paste0('kmeans', clNum)

# scale mean count + cluster
clusterGene <- scaleCount %>%
  as.data.frame %>%
  rownames_to_column(var = 'ID') %>%
  as_tibble %>%
  {
    cl <- as.data.frame(cl) %>%
      rownames_to_column(var = 'ID')
    inner_join(., cl)
  }

## plot core cluster
clusterCore <- clusterGene %>%
  group_by(cl) %>%
  summarise_at(-1, mean, na.rm = TRUE) %>% ## mean of each cluster
  mutate(cl = paste0('cluster_', cl) %>%
           factor(levels = paste0('cluster_', cl))) %>%
  gather(Sample, NorExpress, -1) %>%
  mutate(Sample = Sample %>% factor(levels = sampleN, ordered = TRUE))

ggplot(clusterCore, aes(Sample, NorExpress, col = cl, group = cl)) +
  geom_point() +
  geom_line() +
  facet_wrap(. ~ cl, ncol = 2) +
  ylab('Scaled counts') +
  scale_color_manual(values = hue_pal()(clNum),
                     breaks = kClust10$cluster %>%
                       table %>%
                       names %>%
                       paste0('cluster_', .),
                     labels = kClust10$cluster %>%
                       table %>%
                       {paste0('cluster_', names(.), ' ', .)},
                     guide = guide_legend(title = paste0('kmeans (k = ',clNum, ')'))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(outputFolder,paste0(prefix, '.pdf')))
ggsave(file.path(outputFolder,paste0(prefix, '.jpg')))

## plot all genes
clusterGenePlot <- clusterGene %>%
  gather(Sample, NorExpress, -ID, -cl) %>%
  mutate(cl = paste0('cluster_', cl) %>%
           factor(levels = paste0('cluster_', sort(unique(cl))))) %>%
  mutate(Sample = Sample %>% factor(levels = sampleN, ordered = TRUE))

clusterCorePlot <- clusterCore %>%
  dplyr::mutate(ID = 1 : nrow(clusterCore))
ggplot(clusterGenePlot, aes(Sample, NorExpress, group = ID)) +
  geom_line(color = 'grey30', alpha = 0.01) +
  facet_wrap(. ~ cl, ncol = 2) +
  geom_point(data = clusterCorePlot, aes(Sample, NorExpress, col = cl, group = ID)) +
  geom_line(data = clusterCorePlot, aes(Sample, NorExpress, group = cl, col = cl)) +
  ylab('Scaled counts') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(colour = guide_legend(title = paste0('kmeans (k = ',clNum, ')')))
ggsave(file.path(outputFolder,paste0(prefix, '_genes.pdf')), width = 10, dpi = 320)
ggsave(file.path(outputFolder,paste0(prefix, '_genes.jpg')), width = 10, dpi = 320)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~cluster cor phenotype~~~~~~~~~~~~~~~~~
traits <- data.frame(age = c(0, 1, 1, 0),
                     infection = c(1, 1, 0, 0))

cores <- clusterGene %>%
  group_by(cl) %>%
  summarise_at(2:5, mean, na.rm = TRUE) %>%
  mutate(cl = cl %>% paste0('cluster_', .)) %>%
  column_to_rownames(var = 'cl') %>%
  t
#> cores
#    cluster_1   cluster_2  cluster_3  cluster_4  cluster_5  cluster_6
#YV -0.7456057  0.06783868 -0.1623165 -1.1890152 -1.0410784  0.5980541
#OV  0.2666544  1.03513371  0.8035259 -0.2265719  1.0713386 -0.9477306
#O   1.1924870  0.11229097 -1.1081243  0.7844676  0.4149209 -0.6294178
#Y  -0.7135357 -1.21526336  0.4669149  0.6311195 -0.4451810  0.9790944
#     cluster_7   cluster_8  cluster_9 cluster_10
#YV  1.31306243 -0.46021674  1.0280185  0.5389170
#OV -0.38205417 -0.86383051  0.5115264 -1.1927424
#O  -0.83315047  0.04944028 -0.5433020  0.7022405
#Y  -0.09664146  1.27460697 -0.9962430 -0.0484151

moduleTraitCor <- cor(cores, traits, use = 'p')
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(traits))

traitPPlot <- moduleTraitPvalue %>%
  as.data.frame %>%
  rownames_to_column('cluster') %>%
  gather(trait, pvalue, age : infection) %>%
  as_tibble

traitCorPlot <- moduleTraitCor %>%
  as.data.frame %>%
  rownames_to_column('cluster') %>%
  gather(trait, correlation, age : infection) %>%
  as_tibble %>%
  mutate(x = rep(0 : (ncol(traits) - 1), each = ncol(cores))) %>%
  mutate(y = rep((ncol(cores) - 1) : 0, ncol(traits))) %>%
  inner_join(traitPPlot) %>%
  mutate(addtext = paste0(round(correlation, digit = 2),
                          '\n',
                          '(',
                          round(pvalue, digit = 2),
                          ')'))

ggplot(traitCorPlot, aes(x = x, y = y, fill = correlation)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name = 'RdYlBu')))(100),
                       breaks = seq(-1, 1, 0.5),
                       labels = format(seq(-1, 1, 0.5)),
                       limits = c(-1, 1)) +
  geom_text(aes(label = addtext)) +
  scale_x_continuous(breaks = 0 : (ncol(traits) - 1), labels = c("age","infection")) +
  scale_y_continuous(breaks = 0 : (ncol(cores) - 1), labels = paste0('cluster_', (ncol(cores)):1)) +
  xlab('Trait') +
  ylab('Cluster')
ggsave(file.path(outputFolder,paste0(prefix, '_trait.jpg')))
ggsave(file.path(outputFolder,paste0(prefix, '_trait.pdf')))
# 这个可以看到与表型相关性高的cluster，很有用；
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~heat map~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 不是mean
scaleC <- rawCount %>%
  t %>%
  scale %>%
  t %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  rename_at(-1, .funs = list(~paste0('Scale_', .)))

#> head(scaleC)
## A tibble: 6 × 14
#  ID     Scale…¹ Scale…² Scale…³ Scale…⁴ Scale…⁵ Scale…⁶ Scale…⁷ Scale…⁸ Scale…⁹
#  <chr>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#1 20389   -0.408  -1.18   -0.753 -1.02   -1.01    -0.466 -1.06     0.474   0.729
#2 22287   -1.25   -1.45   -1.47  -1.30    0.925    0.500  0.714    0.908   0.818
#3 10050…   0.639   0.174   0.368 -0.0147 -1.67    -0.364 -1.64    -0.119  -0.756
#4 15122    0.491   0.126   0.360 -0.0161 -1.72    -0.524 -1.62    -0.124  -0.564
#5 16071   -1.02   -1.22   -0.789 -0.776   0.505    1.29   1.05     0.871   1.33
#6 17105   -0.774  -1.73   -1.30  -1.37   -0.0823   0.481 -0.0646   0.883   1.27
#> colnames(scaleC)
# [1] "ID"        "Scale_YV1" "Scale_YV2" "Scale_YV3" "Scale_YV4" "Scale_OV1"
# [7] "Scale_OV2" "Scale_OV3" "Scale_O1"  "Scale_O2"  "Scale_O3"  "Scale_Y1"
#[13] "Scale_Y2"  "Scale_Y3"

rawC <- rawCount %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  rename_at(-1, .funs = list(~paste0('Raw_', .)))
# tibble是data.frame的一种形式，语法和data.frame较为相似，其实大部分功能都与data.frame相似。tibble来自于tidyverse生态系统中的tibble包。
#> head(rawC)
## A tibble: 6 × 14
#  ID       Raw_YV1 Raw_YV2 Raw_YV3 Raw_YV4 Raw_OV1 Raw_OV2 Raw_OV3 Raw_O1 Raw_O2
#  <chr>      <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>  <dbl>  <dbl>
#1 20389      13.2    13.0     13.1    13.1    13.1    13.2    13.0   13.5   13.6
#2 22287      12.4    12.3     12.3    12.4    13.5    13.3    13.4   13.5   13.5
#3 1005036…   12.8    12.5     12.6    12.5    11.7    12.3    11.7   12.4   12.1
#4 15122      12.5    12.3     12.4    12.2    11.4    12.0    11.4   12.2   11.9
#5 16071       9.84    9.62    10.1    10.1    11.5    12.4    12.1   11.9   12.4
#6 17105      11.6    11.4     11.5    11.5    11.8    11.9    11.8   12.0   12.1
## … with 4 more variables: Raw_O3 <dbl>, Raw_Y1 <dbl>, Raw_Y2 <dbl>,
##   Raw_Y3 <dbl>
## ℹ Use `colnames()` to see all variable names
#> colnames(rawC)
# [1] "ID"      "Raw_YV1" "Raw_YV2" "Raw_YV3" "Raw_YV4" "Raw_OV1" "Raw_OV2"
# [8] "Raw_OV3" "Raw_O1"  "Raw_O2"  "Raw_O3"  "Raw_Y1"  "Raw_Y2"  "Raw_Y3"
#

degresC <- deganno %>%
  dplyr::select(ID, YV_vs_Y_pvalue : OV_vs_YV_log2FoldChange)
#> colnames(degresC)
# [1] "ID"                      "YV_vs_Y_pvalue"
# [3] "YV_vs_Y_padj"            "YV_vs_Y_log2FoldChange"
# [5] "OV_vs_O_pvalue"          "OV_vs_O_padj"
# [7] "OV_vs_O_log2FoldChange"  "O_vs_Y_pvalue"
# [9] "O_vs_Y_padj"             "O_vs_Y_log2FoldChange"
#[11] "OV_vs_YV_pvalue"         "OV_vs_YV_padj"
#[13] "OV_vs_YV_log2FoldChange"

#> dim(rawC)
#[1] 12150    14
#> dim(scaleC)
#[1] 12150    14
#> dim(degresC)
#[1] 14139    13

heatPlot <- rawC %>%
  inner_join(scaleC) %>%
  inner_join(degresC) %>%
  {
    cl <- as.data.frame(cl) %>%
      rownames_to_column(var = 'ID')
    inner_join(., cl)
  } %T>%
  {(sum(names(cl) == .$ID) == nrow(.)) %>% print} %>% ## check cl names and degresC row names
   dplyr::slice(cl %>% order)
#> heatPlot
## A tibble: 12,088 × 40
#   ID    Raw_YV1 Raw_YV2 Raw_YV3 Raw_YV4 Raw_OV1 Raw_OV2 Raw_OV3 Raw_O1 Raw_O2
#   <chr>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>  <dbl>  <dbl>
# 1 16071    9.84    9.62   10.1    10.1    11.5    12.4    12.1   11.9   12.4
# 2 13627   10.9    10.9    10.9    10.9    11.1    11.1    11.1   11.1   11.1
# 3 56040   10.7    10.7    10.7    10.7    10.9    11.0    11.0   11.1   11.0
# 4 16019    8.37    8.23    8.56    8.62   10.0    10.4    10.4   10.4   10.8
# 5 22070   10.0    10.0    10.0    10.0    10.2    10.2    10.2   10.3   10.3
# 6 20044    9.97    9.90    9.93    9.93   10.1    10.1    10.1   10.2   10.2
# 7 16149    9.60    9.46    9.52    9.58    9.90   10.1     9.98  10.3   10.3
# 8 22121    9.68    9.71    9.71    9.72    9.90    9.98    9.95   9.98   9.97
# 9 20091    9.59    9.58    9.61    9.61    9.78    9.78    9.80   9.99  10.0
#10 13010    9.69    9.68    9.77    9.71    9.72    9.85    9.78   9.85   9.87
## … with 12,078 more rows, and 30 more variables: Raw_O3 <dbl>, Raw_Y1 <dbl>,
##   Raw_Y2 <dbl>, Raw_Y3 <dbl>, Scale_YV1 <dbl>, Scale_YV2 <dbl>,
##   Scale_YV3 <dbl>, Scale_YV4 <dbl>, Scale_OV1 <dbl>, Scale_OV2 <dbl>,
##   Scale_OV3 <dbl>, Scale_O1 <dbl>, Scale_O2 <dbl>, Scale_O3 <dbl>,
##   Scale_Y1 <dbl>, Scale_Y2 <dbl>, Scale_Y3 <dbl>, YV_vs_Y_pvalue <dbl>,
##   YV_vs_Y_padj <dbl>, YV_vs_Y_log2FoldChange <dbl>, OV_vs_O_pvalue <dbl>,
##   OV_vs_O_padj <dbl>, OV_vs_O_log2FoldChange <dbl>, O_vs_Y_pvalue <dbl>, …
## ℹ Use `print(n = ...)` to see more rows, and `colnames()` to see all variable names

heatRawPlot <- heatPlot %>%
  select(ID, starts_with('Raw')) %>%
  gather(sample, raw, -1) %>%
  mutate(x = rep(0 : 12, each = nrow(heatPlot))) %>% # 这里的数值和下面的数值要根据starts_with('Raw')的样本是变化
  mutate(y = rep(0 : (nrow(heatPlot) - 1), 13))

heatScalePlot <- heatPlot %>%
  select(ID, starts_with('Scale')) %>%
  gather(sample, scale, -1) %>%
  mutate(x = rep(0 : 12, each = nrow(heatPlot))) %>%# 这里的数值和下面的数值要根据starts_with('Raw')的样本是变化
  mutate(y = rep(0 : (nrow(heatPlot) - 1), 13))

heatlog2FCPlot <- heatPlot %>%
  select(ID, ends_with('FoldChange')) %>%
  gather(sample, log2FC, -1) %>%
  mutate(x = rep(0 : 3, each = nrow(heatPlot))) %>%# 这里的数值和下面的数值要根据starts_with('Raw')的样本是变化
  mutate(y = rep(0 : (nrow(heatPlot) - 1), 4))

## sig |FC| > 1 and padj < 0.05
fcsig <- heatPlot %>%
  select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > 1 ~ 1,
                                . < -1 ~ -1,
                                TRUE ~ 0)))

padjsig <- heatPlot %>%
  select(ends_with('padj')) %>%
  abs %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))

heatsigPlot <- (padjsig * fcsig) %>%
  as_tibble %>%
  gather(sample, sig, 1:4) %>% # 这里的数值和下面的数值要根据starts_with('Raw')的样本是变化
  mutate(x = rep(0 : 3, each = nrow(heatPlot))) %>%
  mutate(y = rep(0 : (nrow(heatPlot) - 1), 4))

heatGroupPlot <- heatPlot %>%
  select(ID, cluster = cl) %>%
  mutate(x = 0) %>%
  mutate(y = 0 : (nrow(heatPlot) - 1))

theme_flg22 <- function(...) {
  theme_bw() %+replace%
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks.length = unit(0, 'mm'),
          axis.line = element_blank(),
          panel.spacing = unit(0, 'mm'),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), 'line'),
          legend.spacing = unit(0, 'mm'),
          ...)
}

ggplot(heatRawPlot, aes(x = x, y = y, fill = log2(raw))) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 7, name = 'GnBu'))(100), name = 'log2(count)') +
  scale_x_continuous(breaks = 0 : 12, # 这里的数值和下面的数值要根据starts_with('Raw')的样本是变化
                     labels = rep(sampleN, c(4,3,3,3)) %>% # 加了each后：each = c(4,3,3,3)，只会用第一个4
                       paste(c(rep(1 :4, 1),rep(1:3, 3)), sep = '_')) +
  theme_flg22(legend.position = 'left',
              axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(outputFolder,paste0(prefix, '_heatmap_raw.jpg')))
ggsave(file.path(outputFolder,paste0(prefix, '_heatmap_raw.pdf')))

ggplot(heatScalePlot, aes(x = x, y = y, fill = scale)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name = 'RdYlBu')))(100), name = 'scale(count)') +
  scale_x_continuous(breaks = 0 : 12,# 这里的数值和下面的数值要根据starts_with('Raw')的样本是变化
                     labels = rep(sampleN, c(4,3,3,3)) %>%
                       paste(c(rep(1 :4, 1),rep(1:3, 3)), sep = '_')) +
  theme_flg22(legend.position = 'left',
              axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(outputFolder,paste0(prefix, '_heatmap_scale.jpg')))
ggsave(file.path(outputFolder,paste0(prefix, '_heatmap_scale.pdf')))

ggplot(heatlog2FCPlot, aes(x = x, y = y, fill = log2FC)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 8, name = 'PiYG')))(100), name = 'log2(FoldChange)') +
  scale_x_continuous(breaks = 0 : 3,# 这里的数值和下面的数值要根据starts_with('Raw')的样本是变化
                     labels = c("YV_vs_Y","OV_vs_O","O_vs_Y","OV_vs_YV")) +
  theme_flg22(legend.position = 'left',
              axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(outputFolder,paste0(prefix, '_heatmap_logFC.jpg')))
ggsave(file.path(outputFolder,paste0(prefix, '_heatmap_logFC.pdf')))

ggplot(heatsigPlot, aes(x = x, y = y)) +
  geom_tile(aes(fill = factor(sig))) +
  scale_fill_manual(name = 'Significant', labels = c('Down', 'No', 'Up'), values = c('green', 'grey90', 'red')) +
  scale_x_continuous(breaks = 0 : 3,
                     labels = c("YV_vs_Y","OV_vs_O","O_vs_Y","OV_vs_YV")) +
  theme_flg22(legend.position = 'left',
              axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(outputFolder,paste0(prefix, '_heatmap_sig.jpg')))
ggsave(file.path(outputFolder,paste0(prefix, '_heatmap_sig.pdf')))

ggplot(heatGroupPlot, aes(x = x, y = y)) +
  geom_tile(aes(fill = factor(cluster))) +
  scale_x_continuous(breaks = 0,
                     labels = 'group') +
  theme_flg22(legend.position = 'left',
              axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(outputFolder,paste0(prefix, '_heatmap_group.jpg')))
ggsave(file.path(outputFolder,paste0(prefix, '_heatmap_group.pdf')))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~merge all plots~~~~~~~~~~~~~~~~~~~~
groupe <- ggplot(heatGroupPlot, aes(x = x, y = y)) +
  geom_tile(aes(fill = factor(cluster))) +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_flg22(title = element_blank(),
              #legend.position = 'none'
              )

groupne <- heatGroupPlot %>%
  group_by(cluster) %>%
  summarise(y = median(y)) %>%
  mutate(x = 0, cluster = cluster %>% paste0('cluster', .)) %>%
  ggplot(aes(x = x, y = y, label = cluster)) +
  geom_text() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0),  limits = c(0, nrow(heatGroupPlot)), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_flg22(title = element_blank(),
              #legend.position = 'none'
              )

rawe <- ggplot(heatRawPlot, aes(x = x, y = y, fill = log2(raw))) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 7, name = 'GnBu'))(100), name = 'log2(count)') +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_flg22(title = element_blank(),
              #legend.position = 'none'
              )

scalee <- ggplot(heatScalePlot, aes(x = x, y = y, fill = scale)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral')))(50), name = 'scale(count)') +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_flg22(title = element_blank(),
              #legend.position = 'none'
              )

fce <- ggplot(heatlog2FCPlot, aes(x = x, y = y, fill = log2FC)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 8, name = 'PiYG')))(100), name = 'log2(FoldChange)') +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_flg22(title = element_blank(),
              #legend.position = 'none'
              )

sige <- ggplot(heatsigPlot, aes(x = x, y = y)) +
  geom_tile(aes(fill = factor(sig))) +
  scale_fill_manual(name = 'Significant', labels = c('Down', 'No', 'Up'), values = c('green', 'grey90', 'red')) +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_flg22(title = element_blank(),
              #legend.position = 'none'
              )

sigte <- heatGroupPlot %>%
  select(cluster, y) %>%
  inner_join(heatsigPlot) %>%
  select(sample, sig, x, y, cluster) %>%
  group_by(sample, cluster) %>%
  count(sig) %>%
  spread(sig, n) %>%
  {
    loc <- heatGroupPlot %>%
      group_by(cluster) %>%
      summarise(y = median(y))
    inner_join(., loc)
  } %>%
  rename('down' = `-1`, 'no' = `0`, 'up' = `1`) %>%
  ungroup %>%
  mutate_at(c('down', 'no', 'up'), .funs = list(~if_else(is.na(.), 0L, .))) %>%
  mutate(x = rep(c(0.2, 0.4, 0.6, 0.8), each = max(cl))) %>%
  mutate(signum = paste0(down, '/', no, '/', up)) %>%
  select(signum, x, y) %>%
  ggplot(aes(x = x, y = y, label = signum)) +
  geom_text() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0),  limits = c(0, nrow(heatGroupPlot)), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.1, 0.5), breaks = NULL) +
  theme_flg22(title = element_blank(),
              #legend.position = 'none'
              )


blanke <- ggplot(tibble(x = 0, y = 0 : (nrow(heatPlot) - 1)),
                 aes(x = x, y = y)) +
  geom_tile(colour = 'white') +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_flg22(title = element_blank(),
              #legend.position = 'none'
              )

g <- grid.arrange(groupne,
             groupe,
             blanke,
             rawe,
             blanke,
             scalee,
             blanke,
             fce,
             blanke,
             sige,
             sigte,
             nrow = 1,
             ncol = 11,
             widths = c(3.5, 1, 0.5, 5, 0.5, 13, 0.5, 3, 0.5, 3, 10) %>% {. / sum(.)})

ggsave(file.path(outputFolder,file = paste0(prefix, '_heatmap_all_element.pdf')), plot = g, width=15,height=15)
ggsave(file.path(outputFolder,file = paste0(prefix, '_heatmap_all_element.jpg')), plot = g,width=15,height=15)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


g <- grid.arrange(groupne,
                  #groupe,
                  blanke,
                  scalee,
                  blanke,
                  sige,
                  nrow = 1,
                  ncol = 5,
                  widths = c(3.5, 0.5, 13, 0.5, 5) %>% {. / sum(.)})
ggsave(file.path(outputFolder,file = paste0(prefix, '_heatmap_all.pdf')), plot = g,width=15,height=20)
ggsave(file.path(outputFolder,file = paste0(prefix, '_heatmap_all.jpg')), plot = g,width=15,height=20)

## write the cluster file
# kmeans10.csv文件生成
inner_join(deganno, heatPlot) %>%
  mutate_at(c('Gene', 'Description'), .funs = list(~if_else(is.na(.), '', .))) %>%
  write_csv(paste0("/data/home/mali/ke_data/transcription/count/",prefix, '.csv'))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################################################################
# 选择感兴趣的基因，看表达量，还会用红色标注出平均值，但是没有显著检验；
##########################plot genes############################
cgenes <- c('21926', '16193', '16176', '20296','12981') # 分别代表TNF基因,IL6，IL-1B,CCL2, GM-CSF

gene_exp <- deganno %>%
  filter(ID %in% cgenes) %>%
  select(Gene, YV1 : Y3) %>%
  gather(ID, NormCount, -1) %>%
  mutate(ID = ID %>% substring(1, nchar(.) - 1)) %>% # 去掉ID最后一个字符，这样就能够方便分组
  mutate(ID = factor(ID, levels = sampleN)) %>%
  mutate(Gene = factor(Gene))

genePlot %>%
  ggplot(aes(x = ID, y = NormCount)) +
  geom_dotplot(binaxis='y', stackdir='center') +
  stat_summary(fun.y = mean, geom = 'point', color='red') +
  facet_wrap(. ~ Gene, ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(file.path(outputFolder,file = paste0(prefix, '_selectgene.pdf')))
ggsave(file.path(outputFolder,file = paste0(prefix, '_selectgene.jpg')))
##################################################################


# 这里没有用到
########################separate DEGs############################
library('tibble')
library('readr')
library('dplyr')

setwd('~/plant/MPIPZ_microbe-host_homeostasis-master/results/removeZero/')
library('icesTAF')
mkdir("kmeanssig")
savepath <- '~/plant/MPIPZ_microbe-host_homeostasis-master/results/removeZero/kmeanssig/'

## cgenes <- c('AT1G14550.1', 'AT2G30750.1', 'AT2G19190.1')
## kres %>%
##   slice(which(kres$ID %in% cgenes))

kres <- read_csv('kmeans10.csv',
                 col_types = cols(Chromosome = col_character()))

## base columns
g <- c('Mock_Flg22', 'SynCom33_Flg22', 'SynCom35_Flg22')

for (i in g) {

  eachcols <- c(paste(i, 1:3, sep = '_'),
                paste0(i, c('_vs_Mock_pvalue', '_vs_Mock_padj', '_vs_Mock_log2FoldChange')))

  eachres <- kres %>%
    select(ID : Mock_3, eachcols, cl) %>%
    filter(eachcols[6] %>% get %>% abs %>% {. > 1}) %>%
    arrange(eachcols[6] %>% get %>% desc)

  cls <- eachres$cl %>% unique

  for (j in cls) {
    fname <- i %>% paste0('_vs_Mock_cluster', j, '.csv') %>% file.path(savepath, .)
    eachres %>%
      filter(cl == j) %>%
      mutate(Gene = Gene %>% coalesce('')) %>%
      mutate(Description = Description %>% coalesce('')) %>%
      write_csv(fname)
  }
}
#################################################################
