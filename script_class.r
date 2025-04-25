## load edgeR
library('edgeR')

### Input raw count data into edgeR ###
rawdata <- read.table("counts.txt", header=TRUE, row.names=1)

dim(rawdata)
#
## Read a design file with sample information ###
design<- read.table("design.txt", sep ="\t", row.names=1, header =TRUE)

### Make class label ###
group <- factor(design$Group)

# Make DGEList object ###
y <- DGEList(counts=rawdata, group=group)

# Transformation form raw scale ###
### Before and after pre-processing ###
cpm <- cpm(y)
lcpm <- cpm(y, log=TRUE)

### Filter out lowly expressed genes ###
keep <- rowSums(cpm(y)>2) >= 3
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y)
#[1] 18207    11

# Extracting samplenames ###
samplenames <- substring(colnames(y), 1, nchar(colnames(y)))

#### Plot before and after filtering ###
#library(RColorBrewer)
#nsamples <- ncol(y)
#col <- brewer.pal(nsamples, "Paired")
#par(mfrow=c(1,2))
#plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, main="", xlab="")
#title(main="A. Raw data", xlab="Log-cpm")
##abline(v=0, lty=3)
#for (i in 2:nsamples){
# den <- density(lcpm[,i])
# lines(den$x, den$y, col=col[i], lwd=2)
#}
#legend("topright", samplenames, text.col=col, bty="n")
#lcpm <- cpm(y, log=TRUE)
#plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, main="", xlab="")
#title(main="B. Filtered data", xlab="Log-cpm")
#abline(v=0, lty=3)
#for (i in 2:nsamples){
#   den <- density(lcpm[,i])
#   lines(den$x, den$y, col=col[i], lwd=2)
#}
#legend("topright", samplenames, text.col=col, bty="n")
########################################

### Before normalization ###
preNorm <- y

### Normalization for composition bias
y <- calcNormFactors(y)
y$samples

#### Box-plot between normalized and unnormalized ###
### https://www.bioconductor.org/help/workflows/RNAseq123/ ####
#x2 <- y
#x2$samples$norm.factors <- 1
#x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
#x2$counts[,2] <- x2$counts[,2]*5

par(mfrow=c(1,2))
lcpm <- cpm(preNorm, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="Unnormalised data",ylab="Log-cpm")
#x2 <- calcNormFactors(x2)  
#x2$samples$norm.factors

## [1] 0.0547 6.1306 1.2293 1.1705 1.2149 1.0562 1.1459 1.2613 1.1170

lcpm <- cpm(y, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="Normalised data",ylab="Log-cpm")


### groups: control D7COC D7EV
############### Estimate dispersion ###
###y <- estimateDisp(y,deisgn)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

### Result for control - D7COC ####
et1 <- exactTest(y, pair=c("control","D7COC"))
res <- topTags(et1)

con.D7COC <- topTags(et1, n=nrow(et1$table))
write.table(con.D7COC, file="controlVSD7COC.txt", sep = "\t")

#### Heatmap clustering ####
logCPM <- cpm(y, prior.count=2, log=TRUE)
#rownames(logCPM) <- y$genes$Symbol
#colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")

o1 <- order(et1$table$PValue)
logCPM.con.D7COC <- logCPM[o1[1:30],]

logCPM.con.D7COC <- t(scale(t(logCPM.con.D7COC)))
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")
pdf("con_D7COC_heatmap.pdf")
heatmap.2(logCPM.con.D7COC, col=col.pan, Rowv=TRUE, scale="none", trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none", margin=c(10,9), lhei=c(2,10), lwid=c(2,6))
dev.off()

####### Comapring control vs D7EV

et2 <- exactTest(y, pair=c("control","D7EV"))
res <- topTags(et2)

con.D7EV <- topTags(et2, n=nrow(et2$table))
write.table(con.D7EV, file="controlVSD7EV.txt", sep = "\t")

## Heatmap clustering ####

o2 <- order(et2$table$PValue)
logCPM.con.D7EV <- logCPM[o2[1:30],]

logCPM.con.D7EV <- t(scale(t(logCPM.con.D7EV)))
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")
pdf("con_D7EV_heatmap.pdf")
heatmap.2(logCPM.con.D7EV, col=col.pan, Rowv=TRUE, scale="none", trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none", margin=c(10,9), lhei=c(2,10), lwid=c(2,6))
dev.off()

#### Result for D7COC vs D7EV

et3 <- exactTest(y, pair=c("D7COC","D7EV"))
res <- topTags(et3)

D7COC.D7EV <- topTags(et3, n=nrow(et3$table))
write.table(D7COC.D7EV, file="D7COCVSD7EV.txt", sep = "\t")

### Heatmap clustering ####
o3 <- order(et3$table$PValue)
logCPM.D7COC.D7EV <- logCPM[o3[1:30],]

logCPM.D7COC.D7EV <- t(scale(t(logCPM.D7COC.D7EV)))
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")
pdf("D7COC_D7EV_heatmap.pdf")
heatmap.2(logCPM.D7COC.D7EV, col=col.pan, Rowv=TRUE, scale="none", trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none", margin=c(10,9), lhei=c(2,10), lwid=c(2,6))
dev.off()


###### Addition of Day 14 samples
#### Result for Control vs D14COC


et4 <- exactTest(y, pair=c("control","D14COC"))
res <- topTags(et4)

control.D14COC <- topTags(et4, n=nrow(et4$table))
write.table(control.D14COC, file="controlVSD14COC.txt", sep = "\t")

### Heatmap clustering ####
o4 <- order(et4$table$PValue)
logCPM.control.D14COC <- logCPM[o4[1:30],]

logCPM.control.D14COC <- t(scale(t(logCPM.control.D14COC)))
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")
pdf("control_D14COC_heatmap.pdf")
heatmap.2(logCPM.control.D14COC, col=col.pan, Rowv=TRUE, scale="none", trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none", margin=c(10,9), lhei=c(2,10), lwid=c(2,6))
dev.off()

#### Result for D7COC vs D14COC

et5 <- exactTest(y, pair=c("D7COC","D14COC"))
res <- topTags(et5)

D7COC.D14COC <- topTags(et5, n=nrow(et5$table))
write.table(D7COC.D14COC, file="D7COCVSD14COC.txt", sep = "\t")

### Heatmap clustering ####
o5 <- order(et5$table$PValue)
logCPM.D7COC.D14COC <- logCPM[o5[1:30],]

logCPM.D7COC.D14COC <- t(scale(t(logCPM.D7COC.D14COC)))
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")
pdf("D7COC_D14COC_heatmap.pdf")
heatmap.2(logCPM.D7COC.D14COC, col=col.pan, Rowv=TRUE, scale="none", trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none", margin=c(10,9), lhei=c(2,10), lwid=c(2,6))
dev.off()

###### Heatmap of samples
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

###### PCA

plotPCA(vsd, intgroup = c("dex", "cell"))

###### plot MDS

lcpm <- cpm(y, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)

plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")

