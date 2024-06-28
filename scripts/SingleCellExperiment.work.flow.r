
library(SingleCellExperiment)

library(SingleR)
library(scry)
library(scran)
library(scater)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(multcompView)

sc_exon   = read.csv('../data/GSE123454_cells_exon_counts.csv.gz', row.names = 1)
sc_intron = read.csv('../data/GSE123454_cells_intron_counts.csv.gz', row.names = 1)

sc_exon <- CreateSeuratObject(counts = sc_exon, project = "s_exon")
sc_intron <- CreateSeuratObject(counts = sc_intron, project = "s_introl")


sc_exon@meta.data$orig.ident <- "s_exon"
Idents(sc_exon) <- "orig.ident"

sc_intron@meta.data$orig.ident <- "s_exon"
Idents(sc_intron) <- "orig.ident"

common_features <- intersect(row.names(sc_exon), row.names(sc_intron))

sc_exon <- subset(sc_exon, features = common_features)
sc_intron <- subset(sc_intron, features = common_features)

sc_exon <- GetAssayData(sc_exon, slot = "counts") + GetAssayData(sc_intron, slot = "counts") 


sc_exon <- CreateSeuratObject(counts = sc_exon, project = "both")


sc_exon@meta.data$orig.ident <- "both"
Idents(sc_exon) <- "orig.ident"

sc_exon <- as.SingleCellExperiment(sc_exon)

sc_exon <- logNormCounts(sc_exon)
sc_exon <- runPCA(sc_exon)
sc_exon
cl <- quickCluster(sc_exon)
sc_exon <- computeSumFactors(sc_exon, clusters = cl)
sc_exon <- logNormCounts(sc_exon)
sc_exon <- runUMAP(sc_exon, dimred = "PCA")


tiff("librarySizeFactors.tiff", units="in", width=7, height=7, res=300)

plot(librarySizeFactors(sc_exon), sizeFactors(sc_exon), xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16,
     col=as.integer(factor(cl)))

dev.off()


allzero <- rowMeans(counts(sc_exon)==0)==1

sc_exon <- sc_exon[which(!allzero),]

set.seed(1001)
dec.pbmc <- modelGeneVar(sc_exon)
head(dec.pbmc)



tiff("Mean Variance log expression.tiff", units="in", width=7, height=7, res=300)

plot(dec.pbmc$mean, dec.pbmc$total, pch=16, cex=0.5,
     xlab="Mean of log-expression", ylab="Variance of log-expression")
curfit <- metadata(dec.pbmc)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)

dev.off()

top.pbmc <- getTopHVGs(dec.pbmc, n=3000)


sc_exon <- runPCA(sc_exon, subset_row=top.pbmc)

plotPCA(sc_exon)

sc_exon <- runUMAP(sc_exon, dimred="PCA")
plotUMAP(sc_exon)

g <- buildSNNGraph(sc_exon, k=10, use.dimred = 'PCA')


library(GGally)
ggnet2(g, size=1)

library(igraph)
clust <- igraph::cluster_louvain(g)




pppp <- ggnet2(g, color=membership(clust), size=1)
ggsave("PCA cluste.tiff", dpi = 300, width = 7, height = 7)



sc_exon$Louvain <- factor(membership(clust))


tiff("Umap1.tiff", units="in", width=7, height=7, res=300)


plotUMAP(sc_exon, colour_by="Louvain")

dev.off()

markers <- findMarkers(sc_exon, sc_exon$Louvain)


mm <- unique(unlist(lapply(markers, function(x) rownames(x)[1:10])))



library(SingleR)
library(celldex)

library(scRNAseq)

sce <- fetchDataset("zeisel-brain-2015", "2023-12-14")

sce <- sce[,!is.na(sce$level2class)]

library(scuttle)
sce <- logNormCounts(sce)

pred.grun <- SingleR(test=sc_exon, ref=sce, labels=sce$level2class, de.method="wilcox")


tiff("Cell HeatMap.tiff", units="in", width=7, height=7, res=300)

plotScoreHeatmap(pred.grun)

dev.off()

sc_exon$singler <- factor(pred.grun$pruned.labels)


tiff("Umap cluster.tiff", units="in", width=7, height=7, res=300)

plotUMAP(sc_exon, colour_by = "singler")

dev.off()
