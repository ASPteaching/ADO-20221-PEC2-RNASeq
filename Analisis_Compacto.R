## ----class.source = 'fold-hide', setup, include=FALSE-------------------------
library(knitr)
library(rmdformats)

## Global options
options(max.print="75")
opts_chunk$set(echo=FALSE,
	             cache=TRUE,
               prompt=FALSE,
               tidy=TRUE,
               comment=NA,
               message=FALSE,
               warning=FALSE)
opts_knit$set(width=75)


## ----instalaPaquetes, librerias, echo=TRUE, eval=FALSE------------------------
## if(!require(BiocManager)){
##       install.packages("BiocManager", dep=TRUE)
## }
## 
## installifnot <- function (pckgName, BioC=TRUE){
##   if(BioC){
##     if(!require(pckgName, character.only=TRUE)){
##       BiocManager::install(pckgName)
##     }
##   }else{
##     if(!require(pckgName, character.only=TRUE)){
##       install.packages(pckgName, dep=TRUE)
##     }
##   }
## }
## installifnot("limma")
## installifnot("edgeR")
## installifnot("org.Hs.eg.db")
## installifnot("clusterProfiler")
## installifnot("dplyr", BioC=FALSE)
## installifnot("gplots", BioC=FALSE)
## installifnot("ggvenn", BioC=FALSE)
## installifnot("pheatmap")
## installifnot("prettydoc", BioC=FALSE)
## installifnot("ggnewscale", BioC=FALSE)


## ----directorios--------------------------------------------------------------
if(!dir.exists("datos")) dir.create("datos")
if(!dir.exists("results")) dir.create("results") 


## -----------------------------------------------------------------------------
counts <- read.csv("datos/RawCounts.csv", row.names = 1)
selectedCounts <- as.matrix(counts)


## ----creaTargets--------------------------------------------------------------
sampleNames <- colnames(selectedCounts)
grupos <- c(rep("COVID", 17), rep("SANO", 17))
colores=c(rep("red", 17), rep("blue", 17))
selectedTargets <- data.frame(samples=sampleNames, group=grupos, cols=colores)
rownames(selectedTargets) <- selectedTargets[,1]


## ----getCPM-------------------------------------------------------------------
library(edgeR)
selectedCounts[1:5,1:6]
counts.CPM <- cpm(selectedCounts)
counts.CPM[1:5,1:6]


## ----subsetThreshold----------------------------------------------------------
thresh <- counts.CPM > 0
keep <- (rowSums(thresh[,1:10]) >= 3) &
        (rowSums(thresh[,11:20]) >= 3)
counts.keep <- counts.CPM[keep,]
dim(counts.CPM)
dim(counts.keep)


## ----comparaMatrices----------------------------------------------------------
head(counts.CPM)
head(counts.keep)


## ----makeDGEObj0--------------------------------------------------------------
dgeObj <- DGEList(counts = counts.keep, 
                  lib.size = colSums(counts.keep),
                  norm.factors = rep(1,ncol(counts.keep)), 
                  samples = selectedTargets,
                  group = selectedTargets$group, 
                  genes = rownames(counts.keep), 
                  remove.zeros = FALSE)
dgeObj


## ----makeDGEObj---------------------------------------------------------------
dim(dgeObj)
dgeObjShort<-dgeObj[,c(1:5, 11:15)]
# Library size information is stored in the samples slot
dgeObjShort$samples
colnames(dgeObjShort$counts)


## ----calcNormFactors----------------------------------------------------------
library(edgeR)
dgeObj_norm <- calcNormFactors(dgeObj)


## ----log2count_norm-----------------------------------------------------------
log2count_norm <- cpm(dgeObj_norm, log=TRUE)


## ----distriCounts1------------------------------------------------------------
par(mfrow=c(2,1))
rawCounts <- dgeObj_norm$counts
boxplot(rawCounts, ylab="CPM",las=2, xlab="", col = dgeObj$samples$cols, cex.axis=0.6, main="Distribución de contajes")
boxplot(log2count_norm, ylab="Log2-CPM",las=2, xlab="", col=dgeObj$samples$cols, cex.axis=0.6, main="Distribución de log(contajes)")
abline(h=median(log2count_norm), col="blue")
par(mfrow=c(1,1))


## -----------------------------------------------------------------------------
sampleDists <- dist(t(log2count_norm))
round(sampleDists,1)


## -----------------------------------------------------------------------------
library(factoextra)
fviz_dist(sampleDists)


## -----------------------------------------------------------------------------
hc <-hclust(sampleDists)
plot(hc,labels = colnames(log2count_norm),main = "Agrpamiento jerárquico de las muestras", cex=0.8)


## -----------------------------------------------------------------------------
col.status <- dgeObj_norm$samples$cols
plotMDS(log2count_norm,col=col.status, main="Status", cex=0.7)


## ----matrizDisenyo------------------------------------------------------------
group = as.factor(dgeObj_norm$samples$group)
design = model.matrix(~ 0 + group)
colnames(design) = gsub("group", "", colnames(design))
row.names(design) = sampleNames
design


## ----matrizContrastes---------------------------------------------------------
cont.matrix = makeContrasts(CONTROLvsCOVID = COVID - SANO,
levels=colnames(design))
cont.matrix


## ----voom---------------------------------------------------------------------
voomObj <- voom(dgeObj_norm, design)
voomObj


## ----ajusteLM-----------------------------------------------------------------
fit <- lmFit(voomObj)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)


## ----topTable-----------------------------------------------------------------
toptab <- topTable(fit.cont,coef=1,sort.by="p", number=nrow(fit.cont))
head(toptab)


## ----volcano------------------------------------------------------------------
volcanoplot(fit.cont,coef=1,highlight=100, main="COVID vs SANO")


## ----deg1VOOM-----------------------------------------------------------------
topGenesBas <- rownames(subset(toptab, (abs(logFC)> 2) & (adj.P.Val < 0.01)))
length(topGenesBas)


## ----mapaDeColores------------------------------------------------------------
library(pheatmap)
mat  <- log2count_norm[topGenesBas, ]
mat  <- mat - rowMeans(mat)
pheatmap(mat)


## -----------------------------------------------------------------------------
y = estimateDisp(dgeObj_norm, design, robust=TRUE)
plotBCV(y)


## -----------------------------------------------------------------------------
fit <- glmQLFit(y, design, robust = TRUE)


## -----------------------------------------------------------------------------
res <- glmQLFTest(fit, contrast = cont.matrix)
head(res)


## -----------------------------------------------------------------------------
topTags_edge <- topTags(res, n=dim(log2count_norm)[1]) # todos los genes
head(topTags_edge)


## -----------------------------------------------------------------------------
topGenes_edge <- rownames(subset(topTags_edge$table, (abs(logFC)> 2) & (FDR < 0.01)))
length(topGenes_edge)


## ----comparaVoomEdgeR---------------------------------------------------------
library(ggvenn)
x = list(LimmaVoom = topGenesBas, edgeR = topGenes_edge)
ggvenn(x, fill_color = c("#0073C2FF", "#EFC000FF"), stroke_size = 0.5, set_name_size = 3)


## -----------------------------------------------------------------------------
topGenes <- union(topGenesBas, topGenes_edge)
length(topGenes)
universe <- rownames(toptab)
length(universe)


## ----anotaTop-----------------------------------------------------------------
library(org.Hs.eg.db)
AnnotationDbi::keytypes(org.Hs.eg.db)
topAnots = AnnotationDbi::select(org.Hs.eg.db, topGenes, c("SYMBOL", "ENTREZID", "GENENAME"),
keytype = "ENSEMBL")
head(topAnots)
dim(topAnots)


## ----anots--------------------------------------------------------------------
univAnots = AnnotationDbi::select(org.Hs.eg.db, universe, c("SYMBOL", "ENTREZID", "GENENAME"), keytype = "ENSEMBL")
head(univAnots)
dim(univAnots)


## ----enrichment---------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
ego = enrichGO(gene = topGenes, 
               universe=universe,
               keyType = "ENSEMBL", 
               OrgDb = org.Hs.eg.db,
               ont="BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE)
head(ego[,-8], n=5)


## ----viewEnrichment1----------------------------------------------------------
dotplot(ego, showCategory=7)


## ----viewEnrichment2----------------------------------------------------------
library(ggplot2)
ego2 = simplify(ego)
cnetplot(ego2, showCategory = 3, cex_category =0.3, cex_label_category =0.7, 
         cex_gene=0.2, cex_label_gene=0.4, circular=TRUE, colorEdge=TRUE)


## ----viewEnrichment3----------------------------------------------------------
library(enrichplot)
goplot(ego2, showCategory=6, cex=0.1)


## ----viewEnrichment5----------------------------------------------------------
heatplot(ego2)


## ----viewEnrichment6----------------------------------------------------------
term_similarity_matrix = pairwise_termsim(ego)
emapplot(term_similarity_matrix, showCategory = 15, group_category=TRUE, group_legend=TRUE)

