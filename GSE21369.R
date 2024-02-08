 Pipeline for Human Genome U133 Plus 2.0 Array Affimetrix
#Author: Marco Antonio Espina Ordoñez
#Date: June/2021

#####################################################################
# To install packages
#####################################################################
#install.packages("devtools")
#devtools::install_github("r-lib/remotes")
#install.packages("BiocManager")
#BiocManager::install("oligoClasses")
#BiocManager::install("oligo")
#BiocManager::install("limma")
#BiocManager::install("GEOquery")
#BiocManager::install("EnhancedVolcano")
#BiocManager::install("pheatmap")
#BiocManager::install("arrayQualityMetrics")
#BiocManager::install("ArrayExpress")
#BiocManager::install("pd.hg.u133.plus.2")
#BiocManager::install("hgu133plus2.db")
#install.packages("AnnotationDbi")
#install.packages("org.Hs.eg.db")
#Plotting and color options packages
#install.packages("RColorBrewer")
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("tidyr")

#####################################################################
# To load libaries 
#####################################################################
library(devtools)
library(remotes)
#General Bioconductor packages
library(oligoClasses)
library(oligo)
library(limma)
library(pheatmap)
library(EnhancedVolcano)
library(GEOquery)
library(arrayQualityMetrics)
library(ArrayExpress)
library(hgu133plus2.db)
#Plotting and color options packages)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)

#####################################################################
baseDir <- "GSE21369"
setwd (baseDir)

#####################################################################
# Daownload dataset from GEO 
#####################################################################
# getGEOSuppFiles("GSE21369")
# list.files("GSE21369")
# 
# untar("GSE21369/GSE21369_RAW.tar", exdir = "GSE21369/CEL")
# list.files("GSE21369/CEL")

pd = read.csv("GSE21369/FenodatasinIPF1,6,3,2ycontrol5.csv", header=TRUE, as.is=TRUE)
pd[,"Filename"] = paste(pd[,"Group"], pd[,"Replicate"], sep=".")
pd

#####################################################################
# To read CEL archives
# The sample GSM595417.CEL was excluded because its quality was low 
#####################################################################
celList <- list.celfiles("GSE21369/CEL", pattern = '*.CEL', full.names=TRUE, listGzipped=TRUE)
raw_data <- read.celfiles(celList)
pData(raw_data) = pd
sampleNames(raw_data) <- pd$Filename

#####################################################################
#Quality analysis (Image from each array)
#####################################################################
# for (i in 1:19)
# {
#   name = paste("image",i,".jpg",sep="")
#   jpeg(name)
#   image(raw_data[,i],main=pd$Filename[i])
#   dev.off()
# }

#####################################################################
#Quality analysis (Pseudoimages)
#####################################################################
Pset = fitProbeLevelModel(raw_data)

# for (i in 1:19)
# {
#   name = paste("pseudoimage",i,".jpg",sep="")
#   jpeg(name)
#   image(Pset,which=i,type="residuals",main=pd$Filename[i])
#   dev.off()
# }

#####################################################################
#Quality Analysis with arrayQualityMetrics (raw_data)
#####################################################################
arrayQualityMetrics(expressionset = raw_data[, 1:14],
                    outdir = "raw_data_report",
                    force = TRUE,
                    do.logtransform = TRUE)

#####################################################################
#Normalization
#####################################################################

norm_data <- oligo::rma(raw_data)
norm_matrix <- Biobase::exprs(norm_data)

#####################################################################
#Quality analysis with arrayQualityMetrics (norm_data)
#####################################################################
arrayQualityMetrics(expressionset = norm_data,
                    outdir = "report_for_normData",
                    force = TRUE)

#####################################################################
#PCA
#####################################################################

PCA <- prcomp(t(norm_matrix), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Group = pData(norm_data)$Group,
                     Sample = pData(norm_data)$Sample,
                     Replicate = pData(norm_data)$Replicate 
)


ggplot(dataGG, aes(PC1, PC2, label=row.names(dataGG))) +
  geom_point(aes( colour = Group))+
  geom_text(size=2)+
  ggtitle("PCA plot of the log-transformed Norm expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("#fa1d19", "#1937fa"))

#####################################################################
#Annotation
#####################################################################
man_threshold <- 4

no_of_samples <- table(paste0(pData(norm_data)$Group))

no_of_samples

samples_cutoff <- min(no_of_samples)

idx_man_threshold <- apply(Biobase::exprs(norm_data), 1,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)

annot_manfiltered <- subset(norm_data, idx_man_threshold)

columns(hgu133plus2.db)

annot_data <- AnnotationDbi::select(hgu133plus2.db,
                                    keys = (featureNames(annot_manfiltered)),
                                    columns = c("SYMBOL", "GENENAME", "ENTREZID"),
                                    keytype = "PROBEID")


annot_data <- subset(annot_data, !is.na(SYMBOL))

annot_grouped <- group_by(annot_data, PROBEID)

annot_summarized <- dplyr::summarize(annot_grouped, no_of_matches = n_distinct(SYMBOL))

head(annot_summarized)

annot_filtered <- filter(annot_summarized, no_of_matches > 1)

head(annot_filtered)

probe_stats <- annot_filtered

nrow(probe_stats)

ids_to_exlude <- (featureNames(annot_manfiltered) %in% probe_stats$PROBEID)

table(ids_to_exlude)

annot_final <- subset(annot_manfiltered, !ids_to_exlude)

validObject(annot_final)

head(annot_data)

fData(annot_final)$PROBEID <- rownames(fData(annot_final))

fData(annot_final) <- left_join(fData(annot_final), annot_data)


#####################################################################
# Differential expression
#####################################################################
casos <- as.factor(pd$Group)
design = model.matrix(~ 0+casos)
colnames(design) = levels(casos)
design

cont.matrix<- makeContrasts(IPF-Control, HP-Control, HP-IPF, levels=design)
cont.matrix

fit<-lmFit(annot_final,design)
fit2= contrasts.fit(fit,cont.matrix)
fit3= eBayes(fit2)
fit3
head(fit3)
table_CD <- topTable(fit3, number = Inf)

head(table_CD)

#Save tables
table_CD = topTable(fit3, coef=1, number=nrow(fit), sort.by= "p",adjust="fdr")
write.table(table_CD, file="DEG1_IPFvsCtrl.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

table_CD = topTable(fit3, coef=2, number=nrow(fit), sort.by= "p",adjust="fdr")
write.table(table_CD, file="DEG2_HPvsCtrl.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

table_CD = topTable(fit3, coef=3, number=nrow(fit), sort.by= "p",adjust="fdr")
write.table(table_CD, file="DEG3_HPvsIPF.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
#####################################################################
#Volcano
#####################################################################
EnhancedVolcano(table_CD,lab = table_CD$SYMBOL,
                x = 'logFC',
                y = 'adj.P.Val',
                ylab = bquote(~adj.P.Val),
                xlab = bquote(~logFC),
                ylim = c(0, 2.7),
                xlim = c(-6, 4.5),
                pCutoff = 0.01,
                FCcutoff = 1,
                col=c('#c2bcbc', '#c2bcbc', '#c2bcbc', '#eb4334'),
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
                title = "             Volcano plot of differential expresion of genes",
                titleLabSize = 13,
                subtitleLabSize = 9,
                captionLabSize = 9,
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
                caption = paste0("Total = ", nrow(table_CD), " variables"),
                subtitle = "Comparation between IPF vs Control                                         pCutoff = 0.01,     FCcutoff = 1 ",
                labSize = 1.0)

#####################################################################
#Differentially expressed genes
#####################################################################

topgenes = table_CD[table_CD[, "adj.P.Val"] < 0.01, ]
dim(topgenes)
head(topgenes)
write.table(topgenes, file="TopGenes_GSE21369.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#Genes genes up and down regulated
topups = topgenes[topgenes[, "logFC"] > 1, ]
dim(topups)
topdowns = topgenes[topgenes[, "logFC"] < -1, ]
dim(topdowns)

DEresults = decideTests(fit3,method='global',adjust.method="BH",p.value=0.01,lfc=1)
DEresults[1:10,]

ups = subset(DEresults, DEresults[,1]==1)
downs = subset(DEresults, DEresults[,1]==-1)
IDs.up = rownames(ups)
IDs.down = rownames(downs)

write.table(IDs.up,row.names=TRUE,col.names=FALSE,quote=TRUE,file="upIDs.txt")
write.table(IDs.down,row.names=TRUE,col.names=FALSE,quote=TRUE,file="downIDs.txt")

#####################################################################
#Venn Diagram
#####################################################################
vennDiagram(DEresults, mar=rep(1,4), names = NULL, cex=c(0.8,1,0.7),
            lwd=2, include=c("up", "down"), counts.col=c("red", "blue"),
            circle.col=c("red", "blue"))

#####################################################################

