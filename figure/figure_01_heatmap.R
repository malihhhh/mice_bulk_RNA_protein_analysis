####################################################################################
####################################################################################
# heatmap selected genes
####################################################################################
####################################################################################

# dependencies
library(Biobase)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(randomcoloR)
library(svglite)

sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")

# constants
name         <- "[Produce figure 1D] "
useVersion   <- "version1"
colorScheme  <- getColorScheme()

# set output folder
outputFolder <- file.path("figures", "output", "figure_01")
# create output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)


####################################################################################
message(name, "Heatmap with all samples for selected proteins with replicates")
####################################################################################

# load data
proteins <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "all","CLCs","all", "proteins.RDS")) 

# 可以突出想要展示的基因，目前没有
proteins_filt <- proteins
# top annotation
ha_top <- HeatmapAnnotation(df = proteins_filt %>%
                              pData() %>%
                              dplyr::select(group #,name
                                ),
                            col = list(group = colorScheme$type #,
                                       #name = colorScheme$subtype.alt
                                       ),
                            show_annotation_name = TRUE,
                            annotation_name_gp = gpar(fontsize = 7),
                            show_legend = TRUE,
                            annotation_label = c("group"#, "Subtype"
                              ))

# bottom annotation
ha_bottom <- HeatmapAnnotation(
  text = anno_text(x = pData(proteins_filt)$name,
                   rot = 90,
                   just = "right",
                   gp = gpar(fontsize = 7),
                   location = unit(2.2, "cm")))

# generate data matrix
mat <- exprs(proteins_filt)
range(mat)
mat_log <- log2(mat)
range(mat_log)
mat <- remove_problematic_combs(mat, 1)

# the heatmap
hm <- Heatmap(name = "log2 protein level",
              matrix = mat_log,
              #name = "protein level",
              #matrix = mat,
              col = colorRamp2(c(-1.5, 0, 1.5),
                               c(colorScheme$blue,
                                 colorScheme$white,
                                 colorScheme$red)),
              na_col = "gray90",
              show_row_names = TRUE,
              row_names_gp = gpar(fontsize = 7),
              show_column_names = FALSE,
              cluster_rows = TRUE,
              clustering_distance_rows = "pearson",
              clustering_method_rows = "ward.D2",
              cluster_columns = TRUE,
              clustering_distance_columns = "pearson",
              clustering_method_columns = "ward.D2",
              show_heatmap_legend = TRUE,
              # column_split = factor(pData(proteins_filt)$hierarCluster_pearson_ward.D2, levels = c(1, 2, 4, 6, 3, 7, 5)),
              cluster_column_slices = FALSE,
              top_annotation = ha_top,
              bottom_annotation = ha_bottom)

#view plot
hm

pdf(file = file.path(outputFolder, "heatmap_selected_proteins.pdf"), width = 10)
hm
dev.off()


####################################################################################
message(name, "Heatmap with all samples for highly variable proteins with replicates")
####################################################################################

# load proteins (highly variable)
proteins <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "all","CLCs","all", "proteins.RDS")) 
proteins <- proteins %>% .[fData(.)$hasHighVariance, ] %>% removeNAsFromESet()

# 可以突出想要展示的基因，目前没有
proteins_filt <- proteins
# top annotation
ha_top <- HeatmapAnnotation(df = proteins_filt %>%
                              pData() %>%
                              dplyr::select(group #,name
                                ),
                            col = list(group = colorScheme$type #,
                                       #name = colorScheme$subtype.alt
                                       ),
                            show_annotation_name = TRUE,
                            annotation_name_gp = gpar(fontsize = 7),
                            show_legend = TRUE,
                            annotation_label = c("group"#, "Subtype"
                              ))

# bottom annotation
ha_bottom <- HeatmapAnnotation(
  text = anno_text(x = pData(proteins_filt)$name,
                   rot = 90,
                   just = "right",
                   gp = gpar(fontsize = 7),
                   location = unit(2.2, "cm")))

# generate data matrix
mat <- exprs(proteins_filt)
range(mat)
mat_log <- log2(mat)
range(mat_log)
mat <- remove_problematic_combs(mat, 1)
mat_log <- remove_problematic_combs(mat_log, 1)
# the heatmap
hm <- Heatmap(name = "log2 protein level",
              matrix = mat_log,
              #name = "protein level",
              #matrix = mat,
              col = colorRamp2(c(-1.5, 0, 1.5),
                               c(colorScheme$blue,
                                 colorScheme$white,
                                 colorScheme$red)),
              na_col = "gray90",
              show_row_names = TRUE,
              row_names_gp = gpar(fontsize = 7),
              show_column_names = FALSE,
              cluster_rows = TRUE,
              clustering_distance_rows = "pearson",
              clustering_method_rows = "ward.D2",
              cluster_columns = TRUE,
              clustering_distance_columns = "pearson",
              clustering_method_columns = "ward.D2",
              show_heatmap_legend = TRUE,
              # column_split = factor(pData(proteins_filt)$hierarCluster_pearson_ward.D2, levels = c(1, 2, 4, 6, 3, 7, 5)),
              cluster_column_slices = FALSE,
              top_annotation = ha_top,
              bottom_annotation = ha_bottom)

#view plot
hm

pdf(file = file.path(outputFolder, "heatmap_highly_variable_proteins.pdf"), width = 10)
hm
dev.off()
