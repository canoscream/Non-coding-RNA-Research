#####################################################################################
# Pipeline to Agilent-019118 Human miRNA Microarray 2.0 G4470B (Probe Name version)
# Date: Agust/2021
#####################################################################################

#####################################################################
# Packages
#####################################################################

library(GEOquery)
library(limma)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)

#####################################################################
# Download dataset from GEO 
#####################################################################
setwd("GSE21394")


#getGEOSuppFiles("GSE21394")
#list.files("GSE21394")


untar("GSE21394/GSE21394_RAW.tar", exdir ="GSE21394")
list.files("GSE21394/data")

targetinfo <- readTargets("Phenodata_HP_VS_IPF_VS_CONTROL.txt",
                          row.names="Sample",
                          sep="")

targetinfo <- readTargets("Phenodata_HP_VS_IPF.txt",
                          row.names="Sample",
                          sep="")
project <- read.maimages(targetinfo$Sample,
                         path="GSE21394/TXT", 
                         source="agilent",
                         green.only=TRUE)

# Quality Analysis
#####################################################################
# Boxplot raw data
#####################################################################
par(mar=c(8,8,5,5), cex=0.3, cex.axis=1.0, cex.lab=1.3)

boxplot(project$E,
        main="Boxplot of log2-intensitites for the Raw data",
        xlab="", ylab=bquote(~Log[2]~expression),
        names= targetinfo$Filename,
        col = c(rep("#b3e2cd", times=6),rep("#fdcdac",times=6),"#cbd5e8",
                "#fdcdac","#cbd5e8",rep("#fdcdac",times=2)),
        las=2,
        boxwex = 0.6,
        staplewex = 0.6,
        outline=F)

boxplot(project$E,
        main="Boxplot of log2-intensitites for the Raw data",
        xlab="", ylab=bquote(~Log[2]~expression),
        names= targetinfo$Filename,
        col = c(rep("#fdcdac",times=9), rep("#cbd5e8",times=2)),
        las=2,
        boxwex = 0.6,
        staplewex = 0.6,
        outline=F)

#####################################################################
# Background correction
#####################################################################

project.bgc <- backgroundCorrect(project, method="normexp", offset=16)

#####################################################################
# Boxplot (Background correction)
#####################################################################

par(mar=c(8,8,5,5), cex=0.3, cex.axis=1.0, cex.lab=1.3)

boxplot(project.bgc$E,
        main="Boxplot of log2-intensitites for the BGC data",
        xlab="", ylab=bquote(~Log[2]~expression),
        names= targetinfo$Filename,
        col= c(rep("#b3e2cd", times=6),rep("#fdcdac",times=6),"#cbd5e8",
          "#fdcdac","#cbd5e8",rep("#fdcdac",times=2)),
        las=2,
        boxwex = 0.6,
        staplewex = 0.6,
        outline=F)

boxplot(project.bgc$E,
        main="Boxplot of log2-intensitites for the BGC data",
        xlab="", ylab=bquote(~Log[2]~expression),
        names= targetinfo$Filename,
        col = c(rep("#fdcdac",times=9), rep("#cbd5e8",times=2)),
        las=2,
        boxwex = 0.6,
        staplewex = 0.6,
        outline=F)

#####################################################################
# Normalization
#####################################################################

project.NormData <-normalizeBetweenArrays(project.bgc,method=
                                            "cyclicloess")

#####################################################################
# Boxplot (After normalization)
#####################################################################

par(mar=c(8,8,5,5), cex=0.3, cex.axis=1.0, cex.lab=1.3)

boxplot(project.NormData$E,
        main="Boxplot of log2-intensitites for the BGC data",
        xlab="", ylab=bquote(~Log[2]~expression),
        names= targetinfo$Filename,
        col= c(rep("#b3e2cd", times=6),rep("#fdcdac",times=6),"#cbd5e8",
               "#fdcdac","#cbd5e8",rep("#fdcdac",times=2)),
        las=2,
        boxwex = 0.6,
        staplewex = 0.6,
        outline=F)

boxplot(project.NormData$E,
        main="Boxplot of log2-intensitites for the BGC data",
        xlab="", ylab=bquote(~Log[2]~expression),
        names= targetinfo$Filename,
        col = c(rep("#fdcdac",times=9), rep("#cbd5e8",times=2)),
        las=2,
        boxwex = 0.6,
        staplewex = 0.6,
        outline=F)

#####################################################################
# PCA (NormData)
#####################################################################

PCA <- prcomp(t(project.NormData$E), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Disease =
                       Biobase::pData(targets)$Group)

ggplot(dataGG, aes(PC1, PC2, label=row.names(dataGG))) +
  geom_text(size=1)+geom_point(aes(colour = targetinfo$Group)) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("#b3e2cd", "#fdcdac", "#cbd5e8" ))

#####################################################################
# Annotation
#####################################################################

Annot1 <- project.NormData$genes
Annot2 <- Annot1[1:13737, c(4, 5, 6)]
geneMatrix <- cbind(Annot2, project.NormData$E)

Annot3 <- (geneMatrix[geneMatrix$'ControlType' == "0",])

write.table(project.NormData$genes, file="E.txt", 
            sep="\t", row.names=F, col.names=TRUE)


write.table(project.NormData$genes, file="H.txt", 
            sep="\t", row.names=F, col.names=TRUE)
#####################################################################
#Diferential Expression
#####################################################################

casos <- as.factor(targetinfo$Group)
design = model.matrix(~ 0+casos)
colnames(design) = levels(casos)
design

cont.matrix<- makeContrasts(IPF-Control, HP-Control, HP-IPF, levels=design)

cont.matrix<- makeContrasts(HP-IPF, levels=design)
cont.matrix

fit<-lmFit(Annot3[4:20],design)
fit2= contrasts.fit(fit,cont.matrix)
fit3= eBayes(fit2)
fit3
head(fit3)
table_CD <- topTable(fit3, number = Inf)

fit<-lmFit(Annot3[4:14],design)
fit2= contrasts.fit(fit,cont.matrix)
fit3= eBayes(fit2)
fit3
head(fit3)
table_CD <- topTable(fit3, number = Inf)

head(table_CD)
#Save tables
table_CD = topTable(fit3, coef=3, number=nrow(fit), 
                    sort.by= "none",adjust="fdr")

#write.table(table_CD, file="DEG_GSE21394_IPF_vs_Ctrl.txt", sep="\t", 
#            row.names=TRUE, col.names=TRUE, quote=FALSE)

#Final Tables
FinalTable <- cbind(Annot3, table_CD)

write.table(FinalTable, file="FinalTable_GSE21394_IPF_vs_Ctrl.txt", sep="\t", 
            row.names = FALSE)

write.table(FinalTable, file="FinalTable_GSE21394_HP_vs_Ctrl.txt", sep="\t", 
            row.names = FALSE)

write.table(FinalTable, file="FinalTable_GSE21394_HP_vs_IPF.txt", sep="\t", 
            row.names = FALSE)
write.table(FinalTable, file="FinalTable_GSE21394_HP_vs_IPFsinCTRLS.txt", sep="\t", 
            row.names = FALSE)

#####################################################################
#Volcano plot
#####################################################################
EnhancedVolcano(FinalTable,lab = FinalTable$"GeneName",
                x = 'logFC',
                y = 'adj.P.Val',
                ylab = bquote(~adj.P.Val),
                xlab = bquote(~logFC),
                ylim = c(0, 3),
                xlim = c(-5, 5),
                pCutoff = 0.05,
                FCcutoff = 1,
                col=c('#c2bcbc', '#a6a2a2', '#f59289', '#eb4334'),
                cutoffLineType = 'solid',
                cutoffLineCol = 'coral4',
                hlineWidth = 0.2,
                cutoffLineWidth = 0.2,
                axisLabSize = 8,
                pointSize = 1.0,
                colAlpha = 1,
                legendPosition = 'bottom',
                legendLabels = c("NS", expression(LogFC), "adj.P.Val", expression(adj.P.Val ~ and
                                                                                  ~ logFC)),
                legendLabSize = 8,
                title = "          Volcano plot of differential expresion of genes",
                titleLabSize = 13,
                subtitleLabSize = 9,
                captionLabSize = 9,
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
                caption = paste0("Total = ", nrow(FinalTable), " genes"),
                subtitle = "Comparation between IPF vs Control                                        pCutoff = 0.05,     FCcutoff = 1 ",
                labSize = 1.0)


topgenes = FinalTable[FinalTable[, "adj.P.Val"] < 0.05, ]

dim(topgenes)
head(topgenes)

#To IPFvsCTRL
write.table(topgenes, file="TopGenes_GSE21394IPFvsctrl.txt", sep="\t", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)
#To HPvsCTRL
write.table(topgenes, file="TopGenes_GSE21394HPvsctrl.txt", sep="\t", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)
#TO HPvsIPF
write.table(topgenes, file="TopGenes_GSE21394HPvsIPF.txt", sep="\t", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)
#Analysis without ctrls
write.table(topgenes, file="TopGenes_GSE21394sinctrls.txt", sep="\t", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)

#Genes up and down regulated
topups = topgenes[topgenes[, "logFC"] > 1, ]
dim(topups)

#To IPF vs ctrl
write.table(topgenes, file="TopUps_GSE21394IPFvsctrls.txt", sep="\t", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)
#To HP vs ctrl
write.table(topgenes, file="TopUps_GSE21394HPvsctrls.txt", sep="\t", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)
#To HP vs IPF
write.table(topgenes, file="TopUps_GSE21394HPvsIPF.txt", sep="\t", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)
#To analysis without ctrls
write.table(topgenes, file="TopUps_GSE21394sinctrls.txt", sep="\t", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)

#Genes up and down regulated
topdowns = topgenes[topgenes[, "logFC"] < -1, ]
dim(topdowns)
#To IPF vs ctrls
write.table(topgenes, file="TopDowns_GSE21394IPFvsctrls.txt", sep="\t", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)
#To HP vs ctrls
write.table(topgenes, file="TopDowns_GSE21394HPvsctrls.txt", sep="\t", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)
#To HP vs ctrls
write.table(topgenes, file="TopDowns_GSE21394HPvsIPF.txt", sep="\t", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)
#To analysis without ctrls
write.table(topgenes, file="TopDowns_GSE21394sin ctrls.txt", sep="\t", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)

#####################################################################
#Heatmap
#####################################################################

#Select data 
glimpse(FinalTable)
data_filtered <- FinalTable %>% filter(adj.P.Val < 0.05 & (logFC > 1 | logFC < -1 ))
nrow(data_filtered)

#Calculate Z value
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
dataZ <- t(apply(data_filtered[,4:18], 1, cal_z_score))
data_filtered2 <- cbind(data_filtered[,2:3], dataZ[,1:6], dataZ[,7:15])
data_filtered3 <- data_filtered2 %>% drop_na

#Color selection
col_fun <- colorRamp2(seq(min(data_filtered3[,3:8]), max(data_filtered3[,3:17]), length = 3), c("#f1a340", "#f7f7f7", "#998ec3"))
col_fun

#Saving
calc_ht_size = function(ht, unit = "inch") {
  pdf()
  ht = draw(ht)
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.on()
  c(w, h)
}

#Legends
lgd <- Legend(col_fun = col_fun, title = "Row Z-Score")

#Heatmap without genes
ht <- Heatmap(as.matrix(data_filtered3[,3:17]),
              name = "Z-Score", column_title = "Differential gene expression heatmap (GSE21394)",  column_title_gp = gpar(fontsize = 17, fontface = "bold"),
              col = col_fun,
              column_order = order(as.numeric(gsub("column", "", colnames(data_filtered3[,3:17])))),
              clustering_distance_rows = "euclidean",
              row_names_gp = gpar(fontsize = 1),
              column_names_gp = gpar(fontsize = 5, fontface = "bold"),
              clustering_method_rows = "ward.D2",
              width = unit(13, "cm"),
              height = unit(15, "cm"))

#Hetmap with genes

ht2 <- Heatmap(as.matrix(data_filtered3[,3:17]),
               name = "Z-Score", column_title = "Differential gene expression heatmap (GSE21394",  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
               col = col_fun,
               column_order = order(as.numeric(gsub("column", "", colnames(data_filtered3[,3:17])))),
               clustering_distance_rows = "euclidean",
               clustering_method_rows = "ward.D2",
               row_names_max_width = max_text_width(
                 rownames(data_filtered3$GeneName),
                 gp = gpar(fontsize = 12)),
               row_labels = data_filtered3$GeneName,
               width = unit(13, "cm"),
               height = unit(16, "cm"),
               row_names_gp = gpar(fontsize = 1, fontface = "bold"),
               column_names_gp = gpar(fontsize = 5, fontface = "bold"),
               show_heatmap_legend = TRUE)

#Plot PDF
size <- calc_ht_size(ht2)
size
pdf("test.pdf", width = size[1], height = size[1])
ht
dev.off()
