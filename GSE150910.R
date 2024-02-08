##########################################################################################################################################
Pipeline to GEO GSE150910 (RNA-seq)
##########################################################################################################################################
# Install packages
##########################################################################################################################################

BiocManager::install("ComplexHeatmap")
BiocManager::install("circlize")
BiocManager::install("EnhancedVolcano")
BiocManager::install("DESeq2")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")

#Plotting
install.packages("RColorBrewer")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyr")

##########################################################################################################################################
# Load libraries
##########################################################################################################################################

library(devtools)
library(remotes)

#General packages of Bioconductor
library(tidyverse)
library(GEOquery)
library(ComplexHeatmap)
library(circlize)
library(EnhancedVolcano)
library(DESeq2) # otra opción es limma, edgeR
library(AnnotationDbi)
library(org.Hs.eg.db)
library(RColorBrewer)

##########################################################################################################################################
#01.- Load phenodata
##########################################################################################################################################

setwd("GSE150910")
#Phenodata <- read.table("PhenodatawithHP.txt", stringsAsFactors = F, sep = "", header =  T)
Phenodata <- read.table("PhenodataonlyHP.txt", stringsAsFactors = F, sep = "", header =  T)

raw_counts <- read.delim("GSE150910_counts.csv", stringsAsFactors = F, sep = ",", row.names = 1)

#Filter with dplyr
<<<<<<< HEAD
raw_counts <- raw_counts %>%
dplyr::select(contains("HP"))
=======
raw_counts <- raw_counts %>% dplyr::select(contains("HP"))
>>>>>>> bd52abf7e8740bd39ec1d91a03b8267f0a69a305

raw_counts <- as.matrix(raw_counts)
class(raw_counts)

# Order the column by name 
raw_counts <- raw_counts[, order(colnames(raw_counts))]
colnames(raw_counts)

##########################################################################################################################################

rownames(Phenodata) <- Phenodata$SampleName

Phenodata <- Phenodata2

all(rownames(Phenodata) %in% colnames(raw_counts))
#(sort(rownames(Phenodata)))==(sort(colnames(raw_counts)))

all(colnames(raw_counts) %in% rownames(Phenodata))

Phenodata$Group <- as.factor(Phenodata$Group)
Phenodata$cluster <- as.factor(Phenodata$cluster)

Phenodata <- Phenodata[ order(rownames(Phenodata)),]

dds = DESeqDataSetFromMatrix(countData = raw_counts,
                             colData = Phenodata,
                             design = ~Group)

dds = DESeqDataSetFromMatrix(countData = raw_counts,
                             colData = Phenodata,
                             design = ~cluster)

##########################################################################################################################################
#02.- PCA
##########################################################################################################################################

rld <- vst(dds, blind=FALSE)
print(rld)

pcaData <- plotPCA(rld, intgroup = "Group", returnData = TRUE)
pcaData <- plotPCA(rld, intgroup = "cluster", returnData = TRUE)
print(pcaData)

percentVar <- round(100 * attr(pcaData, "percentVar"))

#Generate the PCA plot
pdf("PCA_GSE150910_HP.pdf")
ggplot(pcaData, aes(x = PC1, y = PC2, label=pcaData$name, color = Group)) +
#ggplot(pcaData, aes(x = PC1, y = PC2, label=pcaData$name, color = cluster)) +
  geom_point(size = 3) +
  geom_point(size = 3) +
  geom_text(size=3)+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("GSE150910 | PCA plot of the log-transformed norm expression data") +
  theme(plot.title = element_text(hjust = 0.2))+
  scale_shape_manual(values = c(4,15)) +
<<<<<<< HEAD
  scale_color_manual(values = c("#0000ff", "green", "#fb0007")) + theme(plot.title = element_text(face = "bold"))
=======
  scale_color_manual(values = c("#0000ff", "#fb0007", "green")) + 
  theme(plot.title = element_text(face = "bold"))
>>>>>>> bd52abf7e8740bd39ec1d91a03b8267f0a69a305
dev.off()

clust <- kmeans(pcaData[,1:2], centers=2)$cluster
table(clust)
group1 <- clust[clust == 1]
group2 <- clust[clust == 2]

group1 <- data.frame(SampleName=names(group1), cluster=group1, row.names=NULL)
group2 <- data.frame(SampleName=names(group2), cluster=group2, row.names=NULL)
cluster <- rbind(group1,group2)

Phenodata2 <- merge(Phenodata, cluster, by.y="SampleName")
rownames(Phenodata2) <- Phenodata2$SampleName

library(plotly)

fig <- plot_ly(data = pcaData, x = ~PC1, y = ~PC2,
               type = 'scatter', mode = 'markers',
               text = ~paste('Sample: ', name))
fig


##########################################################################################################################################
#03.- Clustering
##########################################################################################################################################
library("FactoMineR")
library("factoextra")

# assay(rld) # To obtain normalized counts
count_matrix <- assay(rld)

### Compute pairwise correlation values
rld_cor <- cor(count_matrix)

head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames

### Load pheatmap package
library(pheatmap)

### Plot heatmap
pheatmap(rld_cor)

heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor, color = heat.colors, border_color=NA, fontsize = 10, 
         fontsize_row = 10, height=20)


##########################################################################################################################################
#04.- Differential expression with DEseq2
##########################################################################################################################################

dds <- DESeq(dds)

model.matrix(~Phenodata$Group)
dds$condition <- relevel(dds$Group, ref = "Control")

#IPF vs Control
res <- results(dds, contrast=c("Group", "IPF", "Control"))
#HP vs Control
res <- results(dds, contrast=c("Group", "HP", "Control"))
#HP vs IPF
res <- results(dds, contrast=c("Group", "IPF", "Control"))
res <- results(dds, contrast=c("Group", "HP", "Control"))
res <- results(dds, contrast=c("Group", "HP", "IPF"))

table_final <- data.frame(res)

columns(org.Hs.eg.db)

Annot <- AnnotationDbi::select(
  org.Hs.eg.db, keys=rownames(table_final),
  columns=c("ENSEMBL","SYMBOL","GENENAME","GENETYPE" ), keytype="SYMBOL")

table_final <- cbind(rownames(table_final),
                     data.frame(table_final, row.names=NULL))

#Create the final table
AnnotFinal <- merge(x=Annot,y=table_final,by.x='SYMBOL',by.y='rownames(table_final)')

#Eliminate the NA
AnnotFinal <- AnnotFinal %>% drop_na

#Export table
write.table(table_final, file = "DEG_GSE150910.txt", quote = F, row.names = F, sep="\t", col.names = T)

data_filtered <- AnnotFinal %>% filter(padj < 0.05 & (log2FoldChange > 1.5 | log2FoldChange < -1.5 ))

nrow(data_filtered)

write.table(data_filtered, file="./06_RESULTADOS/DEG_FILTER_GSE150910.txt", sep="\t", row.names=T, col.names=T, quote=F)

##########################################################################################################################################
#05.- Volcano plot
##########################################################################################################################################

keyvals <- ifelse(AnnotFinal$log2FoldChange <= -1.5 & AnnotFinal$padj <0.05,
                  '#0000ff',  ifelse(AnnotFinal$log2FoldChange >=1.5 & AnnotFinal$padj <0.05,
                                     '#fb0007', '#e3deeb'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#fb0007'] <- 'Up'
names(keyvals)[keyvals == '#e3deeb'] <- 'NS'
names(keyvals)[keyvals == '#0000ff'] <- 'Down'

#Volcano Script
pdf("04_GRAFICOS/Volcano_GSE150910.pdf")
EnhancedVolcano(AnnotFinal,lab = AnnotFinal$SYMBOL,
                x = 'log2FoldChange',
                y = 'padj',
                ylab = bquote(~padj),
                xlab = bquote(~log2FoldChange),
                ylim = c(0, 75),
                xlim = c(-5, 5),
                pCutoff = 0.05,
                FCcutoff = 1.5,
                colCustom = keyvals,
                cutoffLineType = 'solid',
                cutoffLineCol = 'coral4',
                hlineWidth = 0.2,
                cutoffLineWidth = 0.2,
                axisLabSize = 8,
                pointSize = 1.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabels = c("NS", expression(log2FoldChange), "adj.P.Val", expression(padj ~ and
                                                                                           ~ log2FoldChange)),
                legendLabSize = 8,
                title = "     GSE150910 | Volcano plot of differential expression",
                titleLabSize = 13,
                subtitleLabSize = 9,
                captionLabSize = 9,
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
                caption = NULL,
                subtitle = "Comparation between IPF vs Control                          pCutoff = 0.05,     FCcutoff = 1.5",
                labSize = 1.0)
dev.off()


##########################################################################################################################################
#06.- Heatmap
##########################################################################################################################################

cnorm <- counts(dds, normalized = T)

IPF <- as.data.frame(cnorm) %>% dplyr::select(-contains("control"))
Control <- as.data.frame(cnorm) %>% dplyr::select(-contains("ipf"))

cnorm_sort <- data.frame(Control,IPF)
colnames(cnorm_sort)

table_final_filtered <- data.frame(table_final, cnorm_sort)

data_filtered <- table_final_filtered %>% filter(padj < 0.05 & (log2FoldChange > 1.5 | log2FoldChange < -1.5 ))

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
dataZ <- t(apply(data_filtered[,8:213], 1, cal_z_score))
data_filtered2 <- cbind(data_filtered[,2:7], dataZ[,1:103], dataZ[,104:206])
data_filtered3 <- data_filtered2 %>% drop_na

col_fun <- colorRamp2(seq(min(data_filtered3[,7:212]), max(data_filtered3[,7:212]), length = 3), c("#0000ff", "white", "#fb0007"))
col_fun

rwb <- colorRampPalette(colors = c("#0000ff", "white", "#fb0007"))(30)

lgd <- Legend(col_fun = col_fun, title = "Row Z-Score")

pdf("./04_GRAFICOS/Heatmap_GSE150910.pdf")
Heatmap(as.matrix(data_filtered3[,7:212]),
        name = "Z-Score", column_title = "GSE150910 | Differential gene expression heatmap",  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
        col = rwb,
        column_order = order(as.numeric(gsub("column", "", colnames(data_filtered3[,7:212])))),
        clustering_distance_rows = "euclidean",
        row_names_gp = gpar(fontsize = 0),
        column_names_gp = gpar(fontsize = 3, fontface = "bold"),
        clustering_method_rows = "ward.D2",
        width = unit(9, "cm"),
        height = unit(15, "cm"))
dev.off()
