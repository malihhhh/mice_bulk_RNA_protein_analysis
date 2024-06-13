library('tidyverse')
library('ComplexHeatmap')
library('limma')
library('DESeq2')
library('RColorBrewer')
library('magrittr')
library('xlsx')
library('biomaRt')

sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")

source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

colorScheme  <- getColorScheme()

# set output folder
outputFolder <- file.path("figures", "output", "figure_02")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

theme(text=element_text("serif"))


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
  test <- boxplotData %>%
    filter(cl == i) %>%
    dplyr::select(-ID, -cl) %>%
    gather(key = 'Conditions', value = 'ScaleCounts') %>%
    mutate(Conditions = Conditions %>% factor(levels = sampleN)) %>%
    as.data.frame()

    ggplot(test,aes_string(x = "Conditions", y = "ScaleCounts", fill = "Conditions")) +
    geom_boxplot(alpha = 0.25, outlier.color = NA) +
    geom_point(size = 0.7,alpha = 0.4, position = position_jitter(width = 0.25),aes_string(x = "Conditions", y = "ScaleCounts", color = "Conditions"))+
    #geom_jitter(width=0.2, alpha=0.1, size = 1) +
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF" ,"#3C5488FF")) + # 改颜色
    scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF" ,"#3C5488FF")) + # 改颜色
    ggpubr::stat_compare_means( comparisons = comparisons, label = "p.signif",hide.ns = TRUE,vjust = 0.5,size = 2,na.rm = TRUE) +
    #ggpubr::stat_compare_means( comparisons = comparisons, label = "p.signif") +
    ylab('Scaled counts') +
    #theme_classic() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0),
          plot.title = element_text(hjust = 0.5, size = 8, face = 'bold'),
          legend.text.align = 0,
          axis.text = element_text(size = 6),
          axis.title = element_text(size = 7),
          legend.text=element_text(size= 6),
          legend.title = element_text(size = 7),
          legend.position ="none")
  ggsave(file.path(outputFolder,paste0('kmeans10_boxplot', i, '.pdf')),width = 2,height=2.5)
  ggsave(file.path(outputFolder,paste0('kmeans10_boxplot', i, '.jpeg')),width = 2,height=2.5)
}

# 和result已有的结果一摸一样

#######################################################################
