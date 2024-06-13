library(ggplot2)
library(scales)
library(ggthemes)
library(svglite)

source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

width <- 5
height <- 4.5

# set output folder
outputFolder <- file.path("figures", "output", "figure_01")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# peptide table
df_a <- as.data.frame(c(54094))
colnames(df_a) <- "Peptides"
df_a$Type <- c("Peptides")
df_a$Peptide <- c(">1 sets")

#protein table
df_b <- as.data.frame(c(4840, 4098, 3136))
colnames(df_b) <- "Gene Symbols"
df_b$Type <- c("Total Proteins","Proteins \n in all sets","Overlap \n Proteins-mRNA")

#mrna table
df_c <- as.data.frame(c(12715, 623, 472, 130))
colnames(df_c)<- "mRNA"
df_c$Type <- c("Protein Coding","lncRNA","Pseudogenes","snoRNA")
df_c$label <- c("White","White","White","Black")

#gene table
df_d <- as.data.frame(c(32770-14139, 14139))
colnames(df_d) <- "Gene"
df_d$Type <- c("Genes", "Genes")
df_d$Peptide <- c(">1 sets", "missing \n <1 sets")


#figures
p_peps <- ggplot(df_a,aes(factor(Type),Peptides,fill=factor(Peptide,c(">1 sets"))))+
  geom_bar(stat="identity")+
  geom_text(label= df_a$Peptides, color = "black", vjust=2)+
  scale_fill_manual(values = c("#742c55"),name="Quantification")+
  xlab(element_blank())+
  scale_y_continuous(labels = comma)+
 theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=0))

ggsave(
  filename = "bars_peptides.pdf",
  plot = p_peps,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = width + 0.7,
  height = height+0.7,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)


p_prots <- ggplot(df_b,aes(factor(Type,levels = c("Total Proteins","Proteins \n in all sets","Overlap \n Proteins-mRNA")),`Gene Symbols`,fill=factor(Type,levels = c("Overlap \n Proteins-mRNA","Proteins \n in all sets","Total Proteins"))))+
  geom_bar(stat="identity",show.legend = F)+
  geom_text(label= df_b$`Gene Symbols`, color = "black", vjust=1.5)+
  scale_fill_manual(values = c("#474c64","#443158","#742c55"))+
  xlab(element_blank())+
 theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=0))


ggsave(
  filename = "bars_proteins.pdf",
  plot = p_prots,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = width - 1.3,
  height = height+0.7,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)


p_mrna <- ggplot(df_c,aes(factor(Type,levels = c("Protein Coding","lncRNA","Pseudogenes","snoRNA")),mRNA,fill=factor(Type,levels =  c("snoRNA","Pseudogenes","lncRNA","Protein Coding"))))+
  geom_bar(stat="identity",show.legend = F)+
  geom_text(data = df_c[df_c$mRNA != "130",], label = df_c$mRNA[1:3],color = "black", vjust=1.5)+
  geom_text(data= df_c[df_c$mRNA == "130",], label = df_c$mRNA[4],color = "black", vjust=-0.5)+
  scale_fill_manual(values = c("#f3e5ab","#cad7ae","#a2cab2","#7abcb5"))+
  xlab(element_blank())+
 theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=0))

ggsave(
  filename = "bars_mrna.pdf",
  plot = p_mrna,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = width,
  height = height+0.7,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)

  p_genes <- ggplot(df_d,aes(factor(Type),Gene,fill=factor(Peptide,c(">1 sets", "missing \n <1 sets"))))+
    geom_bar(stat="identity")+
    geom_text(label= df_d$Gene, color = "black", vjust=2)+
    scale_fill_manual(values = c("#742c55","#443158"),name="Quantification")+
    xlab(element_blank())+
    scale_y_continuous(labels = comma)+
   theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=0))

  ggsave(
    filename = "bars_genes.pdf",
    plot = p_genes,
    device = NULL,
    path = outputFolder,
    scale = 1,
    width = width + 0.7,
    height = height+ 0.7,
    units = c("cm"),
    dpi = 600,
    limitsize = FALSE)
