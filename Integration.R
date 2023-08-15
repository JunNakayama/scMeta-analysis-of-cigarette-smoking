library(Seurat)
library(harmony)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(MAST)
library(DoubletFinder)

### Import dataset as seurat objects
GSE173896 <- readRDS("~/Analysis/Tobacco/meta/COPD.rds")
GSE122960 <- readRDS("~/Analysis/Tobacco/meta/GSE122960.rds")
GSE123405 <- readRDS("~/Analysis/Tobacco/meta/GSE123405.rds")
GSE130148 <- readRDS("~/Analysis/Tobacco/meta/GSE130148.rds")
GSE131907 <- readRDS("~/Analysis/Tobacco/meta/GSE131907.rds")
GSE135893 <- readRDS("~/Analysis/Tobacco/meta/GSE135893.rds")
GSE136831 <- readRDS("~/Analysis/Tobacco/meta/GSE136831.rds")
EGA00001004082 <- readRDS("~/Analysis/Tobacco/meta/EGA00001004082.rds")


D <- merge(GSE173896, y = c(GSE122960, GSE123405, GSE130148, GSE131907, GSE135893, GSE136831, EGA00001004082), merge.data = TRUE)
rm(GSE173896, GSE122960, GSE123405, GSE130148, GSE131907, GSE135893, GSE136831, EGA00001004082)


######################
######################
######################

## Integration of meta-data


metadata <- read.csv("metadata2.csv")
rownames(metadata) <- metadata[,1]

D@meta.data <- metadata
cl = c("Non-smoker", "Smoker", "COPD")
D <- subset(D, subset = Class %in% cl)

# > D
# An object of class Seurat 
# 78147 features across 374658 samples within 1 assay 
# Active assay: RNA (78147 features, 0 variable features)

COPD <- readRDS("~/Analysis/Tobacco/meta/COPD.rds")
NAME <- rownames(COPD)
rm(COPD)
#> length(NAME)
#[1] 33538
D <- D[rownames(D) %in% NAME,]

D[["percent.mt2"]] <- PercentageFeatureSet(D, pattern = c("MT-"))

VlnPlot(D, features = "percent.mt2", group.by = "Dataset")

D <- subset(D, subset = nFeature_RNA > 1000 & percent.mt2 < 20)

D <- NormalizeData(D)
D <- FindVariableFeatures(D)
D <- ScaleData(D, vars.to.regress = "percent.mt2")
D <- RunPCA(D, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)



######################
######################
######################

### Removal of Doublet by DoubletFinder
DCY  <- SplitObject(D, split.by = "Dataset")

CY1 <- SplitObject(DCY[[1]], split.by = "orig.ident")
CY1 <- lapply(X = CY1, FUN = function(x) {

	sweep.res.list_CY <- paramSweep_v3(x, PCs = 1:15, sct = FALSE)
	sweep.stats_CY <- summarizeSweep(sweep.res.list_CY, GT = FALSE)
	bcmvn_CY <- find.pK(sweep.stats_CY)

	homotypic.prop <- modelHomotypic(x@active.ident) 
	nExp_poi <- round(0.075*nrow(x@meta.data))
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	x <- doubletFinder_v3(x, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	list <- x@meta.data[,length(x@meta.data)]
	x@meta.data <- x@meta.data[,-c(16,17)]
	x@meta.data$DoubletFinder <- list
	x
})
CY2 <- SplitObject(DCY[[2]], split.by = "orig.ident")
CY2 <- lapply(X = CY2, FUN = function(x) {

	sweep.res.list_CY <- paramSweep_v3(x, PCs = 1:15, sct = FALSE)
	sweep.stats_CY <- summarizeSweep(sweep.res.list_CY, GT = FALSE)
	bcmvn_CY <- find.pK(sweep.stats_CY)

	homotypic.prop <- modelHomotypic(x@active.ident) 
	nExp_poi <- round(0.075*nrow(x@meta.data))
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	x <- doubletFinder_v3(x, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	list <- x@meta.data[,length(x@meta.data)]
	x@meta.data <- x@meta.data[,-c(16,17)]
	x@meta.data$DoubletFinder <- list
	x
})
CY3 <- SplitObject(DCY[[3]], split.by = "orig.ident")
CY3 <- lapply(X = CY3, FUN = function(x) {

	sweep.res.list_CY <- paramSweep_v3(x, PCs = 1:15, sct = FALSE)
	sweep.stats_CY <- summarizeSweep(sweep.res.list_CY, GT = FALSE)
	bcmvn_CY <- find.pK(sweep.stats_CY)

	homotypic.prop <- modelHomotypic(x@active.ident) 
	nExp_poi <- round(0.075*nrow(x@meta.data))
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	x <- doubletFinder_v3(x, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	list <- x@meta.data[,length(x@meta.data)]
	x@meta.data <- x@meta.data[,-c(16,17)]
	x@meta.data$DoubletFinder <- list
	x
})
CY4 <- SplitObject(DCY[[4]], split.by = "orig.ident")
CY4 <- lapply(X = CY4, FUN = function(x) {

	sweep.res.list_CY <- paramSweep_v3(x, PCs = 1:15, sct = FALSE)
	sweep.stats_CY <- summarizeSweep(sweep.res.list_CY, GT = FALSE)
	bcmvn_CY <- find.pK(sweep.stats_CY)

	homotypic.prop <- modelHomotypic(x@active.ident) 
	nExp_poi <- round(0.075*nrow(x@meta.data))
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	x <- doubletFinder_v3(x, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	list <- x@meta.data[,length(x@meta.data)]
	x@meta.data <- x@meta.data[,-c(16,17)]
	x@meta.data$DoubletFinder <- list
	x
})
CY5 <- SplitObject(DCY[[5]], split.by = "orig.ident")
CY5 <- lapply(X = CY5, FUN = function(x) {

	sweep.res.list_CY <- paramSweep_v3(x, PCs = 1:15, sct = FALSE)
	sweep.stats_CY <- summarizeSweep(sweep.res.list_CY, GT = FALSE)
	bcmvn_CY <- find.pK(sweep.stats_CY)

	homotypic.prop <- modelHomotypic(x@active.ident) 
	nExp_poi <- round(0.075*nrow(x@meta.data))
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	x <- doubletFinder_v3(x, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	list <- x@meta.data[,length(x@meta.data)]
	x@meta.data <- x@meta.data[,-c(16,17)]
	x@meta.data$DoubletFinder <- list
	x
})
CY6 <- SplitObject(DCY[[6]], split.by = "orig.ident")
CY6 <- lapply(X = CY6, FUN = function(x) {

	sweep.res.list_CY <- paramSweep_v3(x, PCs = 1:15, sct = FALSE)
	sweep.stats_CY <- summarizeSweep(sweep.res.list_CY, GT = FALSE)
	bcmvn_CY <- find.pK(sweep.stats_CY)

	homotypic.prop <- modelHomotypic(x@active.ident) 
	nExp_poi <- round(0.075*nrow(x@meta.data))
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	x <- doubletFinder_v3(x, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	list <- x@meta.data[,length(x@meta.data)]
	x@meta.data <- x@meta.data[,-c(16,17)]
	x@meta.data$DoubletFinder <- list
	x
})

CY7 <- SplitObject(DCY[[7]], split.by = "orig.ident")

CY8 <- SplitObject(DCY[[8]], split.by = "orig.ident")
CY8 <- lapply(X = CY8, FUN = function(x) {

	sweep.res.list_CY <- paramSweep_v3(x, PCs = 1:15, sct = FALSE)
	sweep.stats_CY <- summarizeSweep(sweep.res.list_CY, GT = FALSE)
	bcmvn_CY <- find.pK(sweep.stats_CY)

	homotypic.prop <- modelHomotypic(x@active.ident) 
	nExp_poi <- round(0.075*nrow(x@meta.data))
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	x <- doubletFinder_v3(x, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	list <- x@meta.data[,length(x@meta.data)]
	x@meta.data <- x@meta.data[,-c(16,17)]
	x@meta.data$DoubletFinder <- list
	x
})



for(i in 1:15){
	x <- CY7[[i]]
	sweep.res.list_CY <- paramSweep_v3(x, PCs = 1:15, sct = FALSE)
	sweep.stats_CY <- summarizeSweep(sweep.res.list_CY, GT = FALSE)
	bcmvn_CY <- find.pK(sweep.stats_CY)

	homotypic.prop <- modelHomotypic(x@active.ident) 
	nExp_poi <- round(0.075*nrow(x@meta.data))
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	x <- doubletFinder_v3(x, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	list <- x@meta.data[,length(x@meta.data)]
	x@meta.data <- x@meta.data[,-c(16,17)]
	x@meta.data$DoubletFinder <- list
	CY7[[i]] <- x

}

for(i in 17:length(CY7)){
	x <- CY7[[i]]
	sweep.res.list_CY <- paramSweep_v3(x, PCs = 1:15, sct = FALSE)
	sweep.stats_CY <- summarizeSweep(sweep.res.list_CY, GT = FALSE)
	bcmvn_CY <- find.pK(sweep.stats_CY)

	homotypic.prop <- modelHomotypic(x@active.ident) 
	nExp_poi <- round(0.075*nrow(x@meta.data))
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	x <- doubletFinder_v3(x, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	list <- x@meta.data[,length(x@meta.data)]
	x@meta.data <- x@meta.data[,-c(16,17)]
	x@meta.data$DoubletFinder <- list
	CY7[[i]] <- x

}
CY7[[16]]@meta.data$DoubletFinder <- NA


rm(DCY,D)



Drev <- CY1[[1]]
for(i in 2:length(CY1)){
	Drev <- merge(Drev, y = CY1[[i]])
}
for(i in 1:length(CY2)){
	Drev <- merge(Drev, y = CY2[[i]])
}
for(i in 1:length(CY3)){
	Drev <- merge(Drev, y = CY3[[i]])
}
for(i in 1:length(CY4)){
	Drev <- merge(Drev, y = CY4[[i]])
}
for(i in 1:length(CY5)){
	Drev <- merge(Drev, y = CY5[[i]])
}
for(i in 1:length(CY6)){
	Drev <- merge(Drev, y = CY6[[i]])
}
for(i in 1:length(CY7)){
	Drev <- merge(Drev, y = CY7[[i]])
}
for(i in 1:length(CY8)){
	Drev <- merge(Drev, y = CY8[[i]])
}


rm(CY1,CY2,CY3,CY4,CY5,CY6,CY7,CY8)



######################
######################
######################
######################
######################
######################



s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Drev <- CellCycleScoring(Drev, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)



######################
######################
######################
## SCT normalization

D.list <- SplitObject(Drev, split.by = "orig.ident")


D.list <- lapply(X = D.list , FUN = function(x) {
    x <- SCTransform(x, vars.to.regress = "percent.mt2")
    x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = D.list)
D.list <- PrepSCTIntegration(object.list = D.list, anchor.features = features)

D.list <- lapply(X = D.list, FUN = function(x) {
    x <- RunPCA(x, verbose = FALSE, approx = FALSE, features = features, npcs = 10)
})



Drev <- D.list[[1]]
for(i in 2:length(D.list)){
	Drev <- merge(Drev, y = D.list[[i]])
}




Drev <- FindVariableFeatures(Drev, verbose = FALSE, assay = "RNA")
Drev <- ScaleData(Drev, assay = "RNA", verbose = FALSE)
Drev <- RunPCA(Drev, assay = "RNA", npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
Drev <- RunHarmony(Drev, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, epsilon.harmony = -Inf)
Drev <- RunUMAP(Drev, reduction = "harmony", dims = 1:100, min.dist = 0.5)
Drev <- FindNeighbors(Drev, dims = 1:15, reduction = "harmony")
Drev <- FindClusters(Drev, resolution = 1)



saveRDS(Drev, file = "Drev.rds")


# > Drev
# An object of class Seurat 
# 33538 features across 257663 samples within 1 assay 
# Active assay: RNA (33538 features, 2000 variable features)
#  3 dimensional reductions calculated: pca, harmony, umap


sel = c("Singlet", NA)
sDrev <- subset(Drev, subset = DoubletFinder %in% sel)

sDrev <- RunHarmony(sDrev, "orig.ident", assay.use = "SCT", plot_convergence = TRUE, epsilon.harmony = -Inf)
sDrev <- RunUMAP(sDrev, reduction = "harmony", dims = 1:100, min.dist = 0.5)
sDrev <- FindNeighbors(sDrev, dims = 1:15, reduction = "harmony")
sDrev <- FindClusters(sDrev, resolution = 1)


DefaultAssay(sDrev) <- "SCT"


# > sDrev
# An object of class Seurat 
# 33538 features across 238338 samples within 1 assay 
# Active assay: RNA (33538 features, 2000 variable features)
#  3 dimensional reductions calculated: pca, harmony, umap


saveRDS(sDrev, file = "sDrev.rds")


sDrev <- PrepSCTFindMarkers(sDrev, assay = "SCT", verbose = TRUE)
sDrev.markers <- FindAllMarkers(sDrev, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST", assay = "SCT")
write.table(sDrev.markers, file = 'sDrevMarker-genes.tsv', sep='	')

######################
######################
######################
######################
######################

library(lisi)

hoge <- DimPlot(sDrev, reduction = "umap", raster = FALSE) + NoLegend()
x <- matrix(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2)
meta_data <- data.frame(sDrev@meta.data$Class6, sDrev@meta.data$Dataset)

res <- compute_lisi(x, meta_data, c("sDrev.meta.data.Class6", "sDrev.meta.data.Dataset"))




######################
######################
######################
######################
######################
######################
######################
######################
# Visiualization of the integrated atlas
######################
######################
######################
######################
######################
######################
######################
hoge <- DimPlot(sDrev, reduction = "umap", raster = FALSE) + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$ident)
data <- data[order(data$ident, decreasing = FALSE),]


png("UMAP.png", width = 1000, height = 1000)

ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_brewer(palette = "Dark2") + ## Dark2だとたりない
	NoLegend()

dev.off()


## Density plot

png("density.png", width = 1000, height = 1000)

ggplot(data, mapping = aes(UMAP1, UMAP2)) +
  geom_hex(bins = 80) +
  scale_fill_continuous(type = "viridis") +
  theme_classic() +
  geom_density2d(size = 0.4, bins = 10, colour = "red2")

dev.off()


hoge <- DimPlot(sDrev, reduction = "umap", group.by = "Class6", raster = FALSE) + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Class6)
data <- data[order(data$ident, decreasing = FALSE),]


png("umap.png", width = 1000, height = 1000)

ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_manual(values = c("firebrick", "mediumblue",  "limegreen", "darkorchid", "skyblue2", "magenta2")) + 
	NoLegend()

dev.off()


hoge1 <- RunUMAP(sDrev, reduction = "pca", dims = 1:100, min.dist = 0.5)
hoge <- DimPlot(hoge1, reduction = "umap", group.by = "Class6", raster = FALSE) + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Class6)
data <- data[order(data$ident, decreasing = FALSE),]


png("umap-pca.png", width = 1000, height = 1000)

ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_manual(values = c("firebrick", "mediumblue",  "limegreen", "darkorchid", "skyblue2", "magenta2")) + 
	NoLegend()

dev.off()


hoge <- DimPlot(hoge1, reduction = "umap", group.by = "Dataset", raster = FALSE) + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Dataset)
data <- data[order(data$ident, decreasing = FALSE),]

png("umap-pca-dataset.png", width = 1000, height = 1000)

ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.2, size = 0.25) +
	theme_classic() +
	scale_color_brewer(palette = "Dark2") + 
	NoLegend()

dev.off()


rm(hoge, hoge1)



## Smoking
hoge <- DimPlot(sDrev, reduction = "umap", group.by = "Smoking", raster = FALSE) + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Smoking)
data <- data[order(data$ident, decreasing = FALSE),]


png("smoking.png", width = 1000, height = 1000)

ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_manual(values = c("skyblue2", "red3")) + 
	NoLegend()

dev.off()


hoge <- DimPlot(subset(sDrev, subset = Smoking == "Never"), reduction = "umap", group.by = "Smoking", raster = FALSE) + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Smoking)
data <- data[order(data$ident, decreasing = FALSE),]


png("smokingnev.png", width = 1000, height = 1000)

ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_manual(values = c("skyblue2")) + 
	NoLegend()

dev.off()


hoge <- DimPlot(subset(sDrev, subset = Smoking == "Smoker"), reduction = "umap", group.by = "Smoking", raster = FALSE) + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Smoking)
data <- data[order(data$ident, decreasing = FALSE),]


png("smokingsmo.png", width = 1000, height = 1000)

ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_manual(values = c("red3")) + 
	NoLegend()

dev.off()



hoge <- DimPlot(sDrev, reduction = "umap", group.by = "Dataset", raster = FALSE) + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Dataset)
data <- data[order(data$ident, decreasing = FALSE),]

png("dataset.png", width = 1000, height = 1000)

ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.2, size = 0.25) +
	theme_classic() +
	scale_color_brewer(palette = "Dark2") + 
	NoLegend()

dev.off()







##########################
# EGA00001004082 
# GSE122960
# GSE123405
# GSE130148
# GSE131907
# GSE135893
# GSE136831
# GSE173896

hoge <- DimPlot(subset(sDrev, subset = Dataset == "EGA00001004082"), reduction = "umap", group.by = "Class6", raster = FALSE) + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Class6)
data <- data[order(data$ident, decreasing = FALSE),]

png("EGAS00001004082.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_manual(values = c("firebrick", "mediumblue",  "limegreen", "darkorchid", "skyblue2", "magenta2")) + 
	NoLegend()
dev.off()

hoge <- DimPlot(subset(sDrev, subset = Dataset == "GSE122960"), reduction = "umap", group.by = "Class6", raster = FALSE) + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Class6)
data <- data[order(data$ident, decreasing = FALSE),]

png("GSE122960.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_manual(values = c("firebrick", "mediumblue",  "limegreen", "darkorchid", "skyblue2", "magenta2")) + 
	NoLegend()
dev.off()

hoge <- DimPlot(subset(sDrev, subset = Dataset == "GSE123405"), reduction = "umap", group.by = "Class6", raster = FALSE) + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Class6)
data <- data[order(data$ident, decreasing = FALSE),]

png("GSE123405.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_manual(values = c("firebrick", "mediumblue",  "limegreen", "darkorchid", "skyblue2", "magenta2")) + 
	NoLegend()
dev.off()

hoge <- DimPlot(subset(sDrev, subset = Dataset == "GSE130148"), reduction = "umap", group.by = "Class6", raster = FALSE) + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Class6)
data <- data[order(data$ident, decreasing = FALSE),]

png("GSE130148.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_manual(values = c("firebrick", "mediumblue",  "limegreen", "darkorchid", "skyblue2", "magenta2")) + 
	NoLegend()
dev.off()

hoge <- DimPlot(subset(sDrev, subset = Dataset == "GSE131907"), reduction = "umap", group.by = "Class6", raster = FALSE) + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Class6)
data <- data[order(data$ident, decreasing = FALSE),]

png("GSE131907.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_manual(values = c("firebrick", "mediumblue",  "limegreen", "darkorchid", "skyblue2", "magenta2")) + 
	NoLegend()
dev.off()

hoge <- DimPlot(subset(sDrev, subset = Dataset == "GSE135893"), reduction = "umap", group.by = "Class6", raster = FALSE) + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Class6)
data <- data[order(data$ident, decreasing = FALSE),]

png("GSE135893.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_manual(values = c("firebrick", "mediumblue",  "limegreen", "darkorchid", "skyblue2", "magenta2")) + 
	NoLegend()
dev.off()

hoge <- DimPlot(subset(sDrev, subset = Dataset == "GSE136831"), reduction = "umap", group.by = "Class6", raster = FALSE) + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Class6)
data <- data[order(data$ident, decreasing = FALSE),]

png("GSE136831.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_manual(values = c("firebrick", "mediumblue",  "limegreen", "darkorchid", "skyblue2", "magenta2")) + 
	NoLegend()
dev.off()

hoge <- DimPlot(subset(sDrev, subset = Dataset == "GSE173896"), reduction = "umap", group.by = "Class6", raster = FALSE) + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Class6)
data <- data[order(data$ident, decreasing = FALSE),]

png("GSE173896.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_manual(values = c("firebrick", "mediumblue",  "limegreen", "darkorchid", "skyblue2", "magenta2")) + 
	NoLegend()
dev.off()




######################
## Drawing plots with the selected genes
hoge <- FeaturePlot(sDrev, features = "COL1A2", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]

png("COL1A2all.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey", high = "limegreen") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "PTPRC", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]

png("PTPRCall.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey", high = "darkorchid") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "EPCAM", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]

png("EPCAMall.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey", high = "mediumblue") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "CLDN5", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]

png("CLDN5all.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey", high = "firebrick") + 
	NoLegend()
dev.off()

######################
######################
######################

## FeaturePlot
hoge <- FeaturePlot(sDrev, features = "ACTA2", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allACTA2.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "COL1A1", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allCOL1A1.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "MKI67", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allMKI67.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "CD79A", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allCD79A.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "AGER", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allAGER.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "SFTPC", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allSFTPC.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "KRT5", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allKRT5.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "TP63", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allTP63.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "CD3E", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allCD3E.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "CD8A", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allCD8A.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "MSR1", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allMSR1.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "MARCO", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allMARCO.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "CD14", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allCD14.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "IL3RA", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allIL3RA.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "S100A8", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allS100A8.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "S100A12", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allS100A12.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "CD1C", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allCD1C.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "FOXJ1", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allFOXJ1.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "FOXI1", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allFOXI1.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "MUC5B", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allMUC5B.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "MUC5AC", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allMUC5AC.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "PROX1", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allPROX1.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "ACKR1", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allACKR1.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "MS4A2", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allMS4A2.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "CD8A", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allCD8A.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "SERPINB3", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allSERPINB3.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "PDGFRB", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allPDGFRB.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "FGFR4", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allFGFR4.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

hoge <- FeaturePlot(sDrev, features = "CLEC9A", raster = FALSE)
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]
png("allCLEC9A.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey90", high = "magenta") + 
	NoLegend()
dev.off()

######################
######################
######################
######################
######################

## Cell number count matrix

TT <- data.frame(table(EPI@meta.data$orig.ident, EPI@meta.data$EPIClass))
hoge <- data.frame(table(STR@meta.data$orig.ident, STR@meta.data$STRClass))
TT <- rbind(TT, hoge)
hoge <- data.frame(table(END@meta.data$orig.ident, END@meta.data$ENDClass))
TT <- rbind(TT, hoge)
hoge <- data.frame(table(LYM@meta.data$orig.ident, LYM@meta.data$LYMClass))
TT <- rbind(TT, hoge)
hoge <- data.frame(table(MYE@meta.data$orig.ident, MYE@meta.data$MYEClass))
TT <- rbind(TT, hoge)


TT = dcast(TT, Var1 ~ Var2, value.var = "Freq")

write.table(TT, file = "Table-origident-cellclus-matrix.csv", sep = ",")





#####################################
#####################################
#####################################


sDrev <- RunUMAP(sDrev, reduction = "harmony", dims = 1:100, min.dist = 0.5)
sDrev <- FindNeighbors(sDrev, dims = 1:15, reduction = "harmony")
sDrev <- FindClusters(sDrev, resolution = 1)



#####################################
#####################################
#####################################
STR <- subset(sDrev, subset = Class6 == c("Fibroblastic"))

DefaultAssay(STR) <- "SCT"
STR <- RunUMAP(STR, reduction = "harmony", dims = 1:100, min.dist = 0.5)
STR <- FindNeighbors(STR, dims = 1:15, reduction = "harmony")
STR <- FindClusters(STR, resolution = 1)


DimPlot(STR, label = TRUE, label.size = 10) + NoLegend()
DimPlot(STR, group.by = "Smoking") + NoLegend()




## Removal of Epi/Imm contaminated 
STR <- subset(STR, ident = c(0,1,2,3,4,5,6,7,9,10,11,12,15,16))
DefaultAssay(STR) <- "SCT"
STR <- RunUMAP(STR, reduction = "harmony", dims = 1:100, min.dist = 0.5)
STR <- FindNeighbors(STR, dims = 1:15, reduction = "harmony")
STR <- FindClusters(STR, resolution = 1)

DimPlot(STR, label = TRUE, label.size = 10) + NoLegend()
DimPlot(STR, group.by = "Smoking") + NoLegend()




STR <- PrepSCTFindMarkers(STR, assay = "SCT", verbose = TRUE)
STR.markers <- FindAllMarkers(STR, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST")
write.table(STR.markers, file = 'STRMarker-genes.tsv', sep='	')





## Change the cluster name

modmeta <- STR@meta.data
IDENT <- as.character(STR@active.ident)

clus = c(0)
IDENT[IDENT %in% clus] <- "Pericyte"
clus = c(1,2,4,6,9,10)
IDENT[IDENT %in% clus] <- "AdvFibroblast"
clus = c(5,7,8,14)
IDENT[IDENT %in% clus] <- "AlvFibroblast"
clus = c(11)
IDENT[IDENT %in% clus] <- "Lipofibroblast"
clus = c(3)
IDENT[IDENT %in% clus] <- "Myofibroblast"
clus = c(12)
IDENT[IDENT %in% clus] <- "SmoothMuscle"
clus = c(13)
IDENT[IDENT %in% clus] <- "Methothelial"


modmeta <- cbind(modmeta, STRClass = IDENT)
STR@meta.data <- modmeta


DimPlot(STR, reduction = "umap", group.by = "STRClass", label = TRUE) + NoLegend()
DimPlot(STR, reduction = "umap", group.by = "STRClass") + NoLegend()



saveRDS(STR, file = "STR.rds")



colclass = c("deepskyblue1", "red3")
VlnPlot(STR, features = "ACTA2", split.by = "Smoking", group.by = "STRClass", pt.size = 0, cols = colclass)



hoge <- data.frame(table(STR@meta.data$orig.ident, STR@meta.data$STRClass, STR@meta.data$Smoking, STR@meta.data$Phase))
write.table(hoge, file = "TableSTR-origident-Smo-STRClass-Phase.csv", sep = ",")



##########################











#####################################
EPI <- subset(sDrev, subset = Class6 %in% c("Epithilia","ProliferatingEpithelia"))

DefaultAssay(EPI) <- "SCT"
EPI <- RunUMAP(EPI, reduction = "harmony", dims = 1:100, min.dist = 0.5)
EPI <- FindNeighbors(EPI, dims = 1:15, reduction = "harmony")
EPI <- FindClusters(EPI, resolution = 1)


DimPlot(EPI, label = TRUE, label.size = 10) + NoLegend()
DimPlot(EPI, group.by = "Smoking") + NoLegend()



EPI <- subset(EPI, ident = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24))
DefaultAssay(EPI) <- "SCT"
EPI <- RunUMAP(EPI, reduction = "harmony", dims = 1:100, min.dist = 0.5)
EPI <- FindNeighbors(EPI, dims = 1:15, reduction = "harmony")
EPI <- FindClusters(EPI, resolution = 1)


DimPlot(EPI, label = TRUE, label.size = 10) + NoLegend()
DimPlot(EPI, group.by = "Smoking") + NoLegend()


EPIc <- EPI
Idents(EPIc) <- EPIc@meta.data$EPIClass
EPI.markers <- FindAllMarkers(EPIc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST", assay = "SCT", recorrect_umi = FALSE)
write.table(EPI.markers, file = 'EPIMarker-genes.tsv', sep='	')




## Change the cluster name

modmeta <- EPI@meta.data
IDENT <- as.character(EPI@active.ident)

clus = c(0,10,3,19,22,24)
IDENT[IDENT %in% clus] <- "AT2"
clus = c(6)
IDENT[IDENT %in% clus] <- "Mucous"
clus = c(5,14,18,23)
IDENT[IDENT %in% clus] <- "Club"
clus = c(11)
IDENT[IDENT %in% clus] <- "AT1"
clus = c(1,2)
IDENT[IDENT %in% clus] <- "Basal"
clus = c(21)
IDENT[IDENT %in% clus] <- "Basal-px"
clus = c(4,15)
IDENT[IDENT %in% clus] <- "Basal-d"
clus = c(16)
IDENT[IDENT %in% clus] <- "ProliferatingEpithelia"
clus = c(17,25)
IDENT[IDENT %in% clus] <- "Ionocyte"
clus = c(9)
IDENT[IDENT %in% clus] <- "Goblet"
clus = c(13)
IDENT[IDENT %in% clus] <- "Serous"
clus = c(7,8,12)
IDENT[IDENT %in% clus] <- "Cil-px"
clus = c(20)
IDENT[IDENT %in% clus] <- "Cilia"

modmeta <- cbind(modmeta, EPIClass = IDENT)
EPI@meta.data <- modmeta


DimPlot(EPI, reduction = "umap", group.by = "EPIClass", label = TRUE) + NoLegend()
DimPlot(EPI, reduction = "umap", group.by = "EPIClass") + NoLegend()



saveRDS(EPI, file = "EPI.rds")




#####################################
END <- subset(sDrev, subset = Class6 == c("Endotheila"))

DefaultAssay(END) <- "SCT"
END <- RunUMAP(END, reduction = "harmony", dims = 1:100, min.dist = 0.5)
END <- FindNeighbors(END, dims = 1:15, reduction = "harmony")
END <- FindClusters(END, resolution = 1)


DimPlot(END, label = TRUE, label.size = 10) + NoLegend()
DimPlot(END, group.by = "Smoking") + NoLegend()



## Removal of Epi/Imm contaminated
END <- subset(END, ident = c(0,1,2,3,4,5,6,7,8,11))
DefaultAssay(END) <- "SCT"
END <- RunUMAP(END, reduction = "harmony", dims = 1:100, min.dist = 0.5)
END <- FindNeighbors(END, dims = 1:15, reduction = "harmony")
END <- FindClusters(END, resolution = 1)

DimPlot(END, label = TRUE, label.size = 10) + NoLegend()
DimPlot(END, group.by = "Smoking") + NoLegend()



## Chnge the cluster name

modmeta <- END@meta.data
IDENT <- as.character(END@active.ident)

clus = c(0,1,2,5,10)
IDENT[IDENT %in% clus] <- "Capillary"
clus = c(8)
IDENT[IDENT %in% clus] <- "Artery"
clus = c(6,9)
IDENT[IDENT %in% clus] <- "Capillary-a"
clus = c(4)
IDENT[IDENT %in% clus] <- "Lymphatic"
clus = c(3,7)
IDENT[IDENT %in% clus] <- "Vein"


modmeta <- cbind(modmeta, ENDClass = IDENT)
END@meta.data <- modmeta


DimPlot(END, reduction = "umap", group.by = "ENDClass", label = TRUE) + NoLegend()



saveRDS(END, file = "END.rds")



END <- PrepSCTFindMarkers(END, assay = "SCT", verbose = TRUE)
END.markers <- FindAllMarkers(END, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST")
write.table(END.markers, file = 'ENDMarker-genes.tsv', sep='	')


##########################
##########################
##########################
##########################
##########################
##########################



IMM <- subset(sDrev, subset = Class6 %in% c("Immune", "ProliferatingImmunecell"))

IMM <- RunUMAP(IMM, reduction = "harmony", dims = 1:100, min.dist = 0.5)
IMM <- FindNeighbors(IMM, dims = 1:15, reduction = "harmony")
IMM <- FindClusters(IMM, resolution = 1)

DimPlot(IMM, label = TRUE, label.size = 10) + NoLegend()
DimPlot(IMM, group.by = "Smoking") + NoLegend()




##########################
##########################
##########################
##########################
##########################
##########################

## Separation of Lymphoid and Myeloid
LYM <- subset(IMM, ident = c(0,2,5,17,24,25))

LYM <- RunUMAP(LYM, reduction = "harmony", dims = 1:100, min.dist = 0.5)
LYM <- FindNeighbors(LYM, dims = 1:15, reduction = "harmony")
LYM <- FindClusters(LYM, resolution = 1)

DimPlot(LYM, label = TRUE, label.size = 10) + NoLegend()
DimPlot(LYM, group.by = "Smoking") + NoLegend()



LYM.markers <- FindAllMarkers(LYM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST")
write.table(LYM.markers, file = 'LYMMarker-genes.tsv', sep='	')



## Change the cluster name

modmeta <- LYM@meta.data
IDENT <- as.character(LYM@active.ident)

clus = c(2,7)
IDENT[IDENT %in% clus] <- "NK"
clus = c(4,8,11,14)
IDENT[IDENT %in% clus] <- "NKT"
clus = c(12,13)
IDENT[IDENT %in% clus] <- "Bplasma"
clus = c(10)
IDENT[IDENT %in% clus] <- "Bcell"
clus = c(0,5)
IDENT[IDENT %in% clus] <- "CD8+Tcell"
clus = c(1,3,6,9,15)
IDENT[IDENT %in% clus] <- "Tcell"


modmeta <- cbind(modmeta, LYMClass = IDENT)
LYM@meta.data <- modmeta


DimPlot(LYM, reduction = "umap", group.by = "LYMClass", label = TRUE) + NoLegend()



saveRDS(LYM, file = "LYM.rds")




##########################
MYE <- subset(IMM, ident = c(1,3,4,6,7,8,9,10,11,12,13,14,15,18,19,20,21,22,23))

MYE <- RunUMAP(MYE, reduction = "harmony", dims = 1:100, min.dist = 0.5)
MYE <- FindNeighbors(MYE, dims = 1:15, reduction = "harmony")
MYE <- FindClusters(MYE, resolution = 1)

DimPlot(MYE, label = TRUE, label.size = 10) + NoLegend()
DimPlot(MYE, group.by = "Smoking") + NoLegend()




MYE <- subset(MYE, ident = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))
MYE <- RunUMAP(MYE, reduction = "harmony", dims = 1:100, min.dist = 0.5)
MYE <- FindNeighbors(MYE, dims = 1:15, reduction = "harmony")
MYE <- FindClusters(MYE, resolution = 1)

DimPlot(MYE, label = TRUE, label.size = 10) + NoLegend()
DimPlot(MYE, group.by = "Smoking") + NoLegend()






MYE.markers <- FindAllMarkers(MYE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST")
write.table(MYE.markers, file = 'MYEMarker-genes.tsv', sep='	')



## Change the cluster name

modmeta <- MYE@meta.data
IDENT <- as.character(MYE@active.ident)

clus = c(6,9)
IDENT[IDENT %in% clus] <- "cMonocyte"
clus = c(5,12)
IDENT[IDENT %in% clus] <- "IntermediateMonocyte"
clus = c(2,7,10,13,16,20)
IDENT[IDENT %in% clus] <- "Macrophage"
clus = c(0,1,3,4,15,17,18)
IDENT[IDENT %in% clus] <- "CD68Macrophage"
clus = c(8)
IDENT[IDENT %in% clus] <- "Neutrophil"
clus = c(11)
IDENT[IDENT %in% clus] <- "DC"
clus = c(19)
IDENT[IDENT %in% clus] <- "pDC"
clus = c(14)
IDENT[IDENT %in% clus] <- "MAST"
clus = c(21)
IDENT[IDENT %in% clus] <- "Basophil"

modmeta <- cbind(modmeta, MYEClass = IDENT)
MYE@meta.data <- modmeta


DimPlot(MYE, reduction = "umap", group.by = "MYEClass", label = TRUE) + NoLegend()



saveRDS(MYE, file = "MYE.rds")





## Density plot
hoge <- DimPlot(EPI, reduction = "umap") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$ident)
data <- data[order(data$ident, decreasing = FALSE),]

png("EPIdensity.png", width = 1000, height = 1000)
ggplot(data, mapping = aes(UMAP1, UMAP2)) +
  geom_hex(bins = 80) +
  scale_fill_continuous(type = "viridis") +
  theme_classic() +
  geom_density2d(size = 0.4, bins = 10, colour = "red2")
dev.off()

hoge <- DimPlot(STR, reduction = "umap") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$ident)
data <- data[order(data$ident, decreasing = FALSE),]

png("STRdensity.png", width = 1000, height = 1000)
ggplot(data, mapping = aes(UMAP1, UMAP2)) +
  geom_hex(bins = 80) +
  scale_fill_continuous(type = "viridis") +
  theme_classic() +
  geom_density2d(size = 0.4, bins = 10, colour = "red2")
dev.off()

hoge <- DimPlot(END, reduction = "umap") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$ident)
data <- data[order(data$ident, decreasing = FALSE),]

png("ENDdensity.png", width = 1000, height = 1000)
ggplot(data, mapping = aes(UMAP1, UMAP2)) +
  geom_hex(bins = 80) +
  scale_fill_continuous(type = "viridis") +
  theme_classic() +
  geom_density2d(size = 0.4, bins = 10, colour = "red2")
dev.off()

hoge <- DimPlot(LYM, reduction = "umap") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$ident)
data <- data[order(data$ident, decreasing = FALSE),]

png("LYMdensity.png", width = 1000, height = 1000)
ggplot(data, mapping = aes(UMAP1, UMAP2)) +
  geom_hex(bins = 80) +
  scale_fill_continuous(type = "viridis") +
  theme_classic() +
  geom_density2d(size = 0.4, bins = 10, colour = "red2")
dev.off()

hoge <- DimPlot(MYE, reduction = "umap") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$ident)
data <- data[order(data$ident, decreasing = FALSE),]

png("MYEdensity.png", width = 1000, height = 1000)
ggplot(data, mapping = aes(UMAP1, UMAP2)) +
  geom_hex(bins = 80) +
  scale_fill_continuous(type = "viridis") +
  theme_classic() +
  geom_density2d(size = 0.4, bins = 10, colour = "red2")
dev.off()




## UMAP
hoge <- DimPlot(EPI, reduction = "umap", raster = FALSE, group.by = "EPIClass") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$EPIClass)
data <- data[order(data$ident, decreasing = FALSE),]

png("EPIUMAP.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",  "#FB9A99",  "#CAB2D6", "#B2DF8A","#FF7F00", "#A6CEE3", 
			 "#6A3D9A","#B15928", "#FDBF6F", "#B15928","#A6CEE3", "#6A3D9A",  "#B15928", "#66A61E", "#E31A1C","#FDBF6F", "#666666", "#1F78B4", 
			  "#E6AB02", "#FFFF99", "#FFFF99","#A6761D", "#33A02C", "#FF7F00")) +
	NoLegend()
dev.off()
png("EPIUMAP2.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",  "#FB9A99",  "#CAB2D6", "#B2DF8A","#FF7F00", "#A6CEE3", 
			 "#6A3D9A","#B15928", "#FDBF6F", "#B15928","#A6CEE3", "#6A3D9A",  "#B15928", "#66A61E", "#E31A1C","#FDBF6F", "#666666", "#1F78B4", 
			  "#E6AB02", "#FFFF99", "#FFFF99","#A6761D", "#33A02C", "#FF7F00")) 
dev.off()

hoge <- DimPlot(STR, reduction = "umap", raster = FALSE, group.by = "STRClass") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$STRClass)
data <- data[order(data$ident, decreasing = FALSE),]

png("STRUMAP.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",  "#FB9A99",  "#CAB2D6", "#B2DF8A","#FF7F00", "#A6CEE3", 
			 "#6A3D9A","#B15928", "#FDBF6F", "#B15928","#A6CEE3", "#6A3D9A",  "#B15928", "#66A61E", "#E31A1C","#FDBF6F", "#666666", "#1F78B4", 
			  "#E6AB02", "#FFFF99", "#FFFF99","#A6761D", "#33A02C", "#FF7F00")) +
	NoLegend()
dev.off()
png("STRUMAP2.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",  "#FB9A99",  "#CAB2D6", "#B2DF8A","#FF7F00", "#A6CEE3", 
			 "#6A3D9A","#B15928", "#FDBF6F", "#B15928","#A6CEE3", "#6A3D9A",  "#B15928", "#66A61E", "#E31A1C","#FDBF6F", "#666666", "#1F78B4", 
			  "#E6AB02", "#FFFF99", "#FFFF99","#A6761D", "#33A02C", "#FF7F00")) 
dev.off()

hoge <- DimPlot(END, reduction = "umap", raster = FALSE, group.by = "ENDClass") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$ENDClass)
data <- data[order(data$ident, decreasing = FALSE),]

png("ENDUMAP.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",  "#FB9A99",  "#CAB2D6", "#B2DF8A","#FF7F00", "#A6CEE3", 
			 "#6A3D9A","#B15928", "#FDBF6F", "#B15928","#A6CEE3", "#6A3D9A",  "#B15928", "#66A61E", "#E31A1C","#FDBF6F", "#666666", "#1F78B4", 
			  "#E6AB02", "#FFFF99", "#FFFF99","#A6761D", "#33A02C", "#FF7F00"))  +
	NoLegend()
dev.off()
png("ENDUMAP2.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",  "#FB9A99",  "#CAB2D6", "#B2DF8A","#FF7F00", "#A6CEE3", 
			 "#6A3D9A","#B15928", "#FDBF6F", "#B15928","#A6CEE3", "#6A3D9A",  "#B15928", "#66A61E", "#E31A1C","#FDBF6F", "#666666", "#1F78B4", 
			  "#E6AB02", "#FFFF99", "#FFFF99","#A6761D", "#33A02C", "#FF7F00")) 
dev.off()

hoge <- DimPlot(LYM, reduction = "umap", raster = FALSE, group.by = "LYMClass") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$LYMClass)
data <- data[order(data$ident, decreasing = FALSE),]

png("LYMUMAP.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",  "#FB9A99",  "#CAB2D6", "#B2DF8A","#FF7F00", "#A6CEE3", 
			 "#6A3D9A","#B15928", "#FDBF6F", "#B15928","#A6CEE3", "#6A3D9A",  "#B15928", "#66A61E", "#E31A1C","#FDBF6F", "#666666", "#1F78B4", 
			  "#E6AB02", "#FFFF99", "#FFFF99","#A6761D", "#33A02C", "#FF7F00")) +
	NoLegend()
dev.off()
png("LYMUMAP2.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",  "#FB9A99",  "#CAB2D6", "#B2DF8A","#FF7F00", "#A6CEE3", 
			 "#6A3D9A","#B15928", "#FDBF6F", "#B15928","#A6CEE3", "#6A3D9A",  "#B15928", "#66A61E", "#E31A1C","#FDBF6F", "#666666", "#1F78B4", 
			  "#E6AB02", "#FFFF99", "#FFFF99","#A6761D", "#33A02C", "#FF7F00")) 
dev.off()

hoge <- DimPlot(MYE, reduction = "umap", raster = FALSE, group.by = "MYEClass") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$MYEClass)
data <- data[order(data$ident, decreasing = FALSE),]

png("MYEUMAP.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",  "#FB9A99",  "#CAB2D6", "#B2DF8A","#FF7F00", "#A6CEE3", 
			 "#6A3D9A","#B15928", "#FDBF6F", "#B15928","#A6CEE3", "#6A3D9A",  "#B15928", "#66A61E", "#E31A1C","#FDBF6F", "#666666", "#1F78B4", 
			  "#E6AB02", "#FFFF99", "#FFFF99","#A6761D", "#33A02C", "#FF7F00"))  +
	NoLegend()
dev.off()
png("MYEUMAP2.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",  "#FB9A99",  "#CAB2D6", "#B2DF8A","#FF7F00", "#A6CEE3", 
			 "#6A3D9A","#B15928", "#FDBF6F", "#B15928","#A6CEE3", "#6A3D9A",  "#B15928", "#66A61E", "#E31A1C","#FDBF6F", "#666666", "#1F78B4", 
			  "#E6AB02", "#FFFF99", "#FFFF99","#A6761D", "#33A02C", "#FF7F00")) 
dev.off()





## Smoking
hoge <- DimPlot(EPI, reduction = "umap", group.by = "Smoking") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Smoking)
data <- data[order(data$ident, decreasing = FALSE),]
png("EPISmoking.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_manual(values = c("skyblue2", "red3")) + 
	NoLegend()
dev.off()

hoge <- DimPlot(STR, reduction = "umap", group.by = "Smoking") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Smoking)
data <- data[order(data$ident, decreasing = FALSE),]
png("STRSmoking.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_manual(values = c("skyblue2", "red3")) + 
	NoLegend()
dev.off()

hoge <- DimPlot(END, reduction = "umap", group.by = "Smoking") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Smoking)
data <- data[order(data$ident, decreasing = FALSE),]
png("ENDSmoking.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_manual(values = c("skyblue2", "red3")) + 
	NoLegend()
dev.off()

hoge <- DimPlot(LYM, reduction = "umap", group.by = "Smoking") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Smoking)
data <- data[order(data$ident, decreasing = FALSE),]
png("LYMSmoking.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_manual(values = c("skyblue2", "red3")) + 
	NoLegend()
dev.off()

hoge <- DimPlot(MYE, reduction = "umap", group.by = "Smoking") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Smoking)
data <- data[order(data$ident, decreasing = FALSE),]
png("MYESmoking.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_manual(values = c("skyblue2", "red3")) + 
	NoLegend()
dev.off()



## Dataset
hoge <- DimPlot(EPI, reduction = "umap", group.by = "Dataset") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Dataset)
data <- data[order(data$ident, decreasing = FALSE),]
png("EPIdataset.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_brewer(palette = "Dark2") + 
	NoLegend()
dev.off()

hoge <- DimPlot(STR, reduction = "umap", group.by = "Dataset") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Dataset)
data <- data[order(data$ident, decreasing = FALSE),]
png("STRdataset.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_brewer(palette = "Dark2") + 
	NoLegend()
dev.off()

hoge <- DimPlot(END, reduction = "umap", group.by = "Dataset") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Dataset)
data <- data[order(data$ident, decreasing = FALSE),]
png("ENDdataset.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_brewer(palette = "Dark2") + 
	NoLegend()
dev.off()

hoge <- DimPlot(LYM, reduction = "umap", group.by = "Dataset") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Dataset)
data <- data[order(data$ident, decreasing = FALSE),]
png("LYMdataset.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_brewer(palette = "Dark2") + 
	NoLegend()
dev.off()

hoge <- DimPlot(MYE, reduction = "umap", group.by = "Dataset") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Dataset)
data <- data[order(data$ident, decreasing = FALSE),]
png("MYEdataset.png", width = 1000, height = 1000)
ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.75) +
	theme_classic() +
	scale_color_brewer(palette = "Dark2") + 
	NoLegend()
dev.off()





#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
# Basal-d

Bd <- subset(EPI, subset = EPIClass == c("Basal-d"))

Bd <- RunUMAP(Bd, reduction = "harmony", dims = 1:100, min.dist = 0.5)
Bd <- FindNeighbors(Bd, dims = 1:15, reduction = "harmony")
Bd <- FindClusters(Bd, resolution = 1)


colclass = c("deepskyblue3", "red3")
DimPlot(Bd, label = TRUE, label.size = 10) + NoLegend()
DimPlot(Bd, group.by = "Smoking", pt.size = 1.25, cols = colclass) + NoLegend()



Idents(Bd) <- Bd@meta.data$Smoking
Bd.markers <- FindAllMarkers(Bd, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST", assay = "SCT", recorrect_umi = FALSE)
write.table(Bd.markers, file = 'BdMarker-genes.csv', sep=',', row.names = FALSE)

colclass = c("red3","deepskyblue3")
png("Basal-d_Vln_FOS.png", width = 1000, height = 1000)
	VlnPlot(Bd, "FOS", cols = colclass) + NoLegend()
dev.off()

png("Basal-d_Vln_JUN.png", width = 1000, height = 1000)
	VlnPlot(Bd, "JUN", cols = colclass) + NoLegend()
dev.off()

png("Basal-d_Vln_ATF3.png", width = 1000, height = 1000)
	VlnPlot(Bd, "ATF3", cols = colclass) + NoLegend()
dev.off()

png("Basal-d_Vln_ATP5PO.png", width = 1000, height = 1000)
	VlnPlot(Bd, "AATP5PO", cols = colclass) + NoLegend()
dev.off()

png("Basal-d_Vln_ELOB.png", width = 1000, height = 1000)
	VlnPlot(Bd, "ELOB", cols = colclass) + NoLegend()
dev.off()

png("Basal-d_Vln_AH3F3A.png", width = 1000, height = 1000)
	VlnPlot(Bd, "H3F3A", cols = colclass) + NoLegend()
dev.off()

png("Basal-d_Vln_KRT5.png", width = 1000, height = 1000)
	VlnPlot(Bd, "KRT5", cols = colclass) + NoLegend()
dev.off()

#####################################
#####################################
#####################################
#####################################


#####################################
# Basophil

Bp <- subset(MYE, subset = MYEClass == c("Basophil"))

Bp <- RunUMAP(Bp, reduction = "harmony", dims = 1:5)
Bp <- FindNeighbors(Bp, dims = 1:15, reduction = "harmony")
Bp <- FindClusters(Bp, resolution = 0.3)


colclass = c("deepskyblue3", "red3")
DimPlot(Bp, label = TRUE, label.size = 10) + NoLegend()
DimPlot(Bp, group.by = "Smoking", pt.size = 1.25, cols = colclass) + NoLegend()



Idents(Bp) <- Bp@meta.data$Smoking
Bp.markers <- FindAllMarkers(Bp, only.pos = TRUE, assay = "SCT", recorrect_umi = FALSE)
write.table(Bp.markers, file = 'BpMarker-genes-notMAST.csv', sep=',', row.names = FALSE)




colclass = c("red3","deepskyblue3")


hoge <- VlnPlot(Bp, "JUND", cols = colclass) + NoLegend()
write.table(hoge$data, file = 'Bp-expresssion-JUND.csv', sep=',', row.names = FALSE)

hoge <- VlnPlot(Bp, "FOSB", cols = colclass) + NoLegend()
write.table(hoge$data, file = 'Bp-expresssion-FOSB.csv', sep=',', row.names = FALSE)

hoge <- VlnPlot(Bp, "IGHG1", cols = colclass) + NoLegend()
write.table(hoge$data, file = 'Bp-expresssion-IGHG1.csv', sep=',', row.names = FALSE)

hoge <- VlnPlot(Bp, "IGHG3", cols = colclass) + NoLegend()
write.table(hoge$data, file = 'Bp-expresssion-IGHG3.csv', sep=',', row.names = FALSE)

hoge <- VlnPlot(Bp, "S100A8", cols = colclass) + NoLegend()
write.table(hoge$data, file = 'Bp-expresssion-S100A8.csv', sep=',', row.names = FALSE)

hoge <- VlnPlot(Bp, "FCGR3A", cols = colclass) + NoLegend()
write.table(hoge$data, file = 'Bp-expresssion-FCGR3A.csv', sep=',', row.names = FALSE)

hoge <- VlnPlot(Bp, "TPSB2", cols = colclass) + NoLegend()
write.table(hoge$data, file = 'Bp-expresssion-TPSB2.csv', sep=',', row.names = FALSE)

hoge <- VlnPlot(Bp, "IGHA1", cols = colclass) + NoLegend()
write.table(hoge$data, file = 'Bp-expresssion-IGHA1.csv', sep=',', row.names = FALSE)



png("Basop_Vln_JUND.png", width = 1000, height = 1000)
	VlnPlot(Bp, "JUND", cols = colclass) + NoLegend()
dev.off()

png("Basop_Vln_FOSB.png", width = 1000, height = 1000)
	VlnPlot(Bp, "FOSB", cols = colclass) + NoLegend()
dev.off()

png("Basop_Vln_.png", width = 1000, height = 1000)
	VlnPlot(Bp, "IGHG1", cols = colclass) + NoLegend()
dev.off()

png("Basop_Vln_IGHG3.png", width = 1000, height = 1000)
	VlnPlot(Bp, "IGHG3", cols = colclass) + NoLegend()
dev.off()

png("Basop_Vln_S100A8.png", width = 1000, height = 1000)
	VlnPlot(Bp, "S100A8", cols = colclass) + NoLegend()
dev.off()

png("Basop_Vln_FCGR3A.png", width = 1000, height = 1000)
	VlnPlot(Bp, "FCGR3A", cols = colclass) + NoLegend()
dev.off()

png("Basop_Vln_TPSB2.png", width = 1000, height = 1000)
	VlnPlot(Bp, "TPSB2", cols = colclass) + NoLegend()
dev.off()

png("Basop_Vln_IGHA1.png", width = 1000, height = 1000)
VlnPlot(Bp, "IGHA1", cols = colclass) + NoLegend()
dev.off()




#####################################
#####################################
#####################################
#####################################


write.table(EPI@meta.data, file = "METADATA-EPI.csv", sep = ",")
write.table(STR@meta.data, file = "METADATA-STR.csv", sep = ",")
write.table(END@meta.data, file = "METADATA-END.csv", sep = ",")
write.table(LYM@meta.data, file = "METADATA-LYM.csv", sep = ",")
write.table(MYE@meta.data, file = "METADATA-MYE.csv", sep = ",")



# Heatmap of numeric matrix of selected marker genes

gene = c("ACTA2", "CNN1",  "CSPG4", "PDGFRB",  "MSLN", "PLIN2", "APOE",  "CC3", "COL1A1", "FGFR4",  "SFRP2")


data <- STR[gene,]
Idents(data) <- STR@meta.data$STRClass
mat <- AverageExpression(data)

pheatmap(log(mat$RNA+1))
pheatmap(scale(log(mat$RNA+1)), cluster_rows = FALSE)


write.table(scale(log(mat$RNA+1)), file = "STRmarkergene-average.csv", sep = ",")





gene = c("ACKR1", "CPE", "GJA5", "DKK2", "EDNRB", "S100A3",  "PROX1", "IL7R", "SLC6A4", "CA4")


data <- END[gene,]
Idents(data) <- END@meta.data$ENDClass
mat <- AverageExpression(data)

pheatmap(log(mat$RNA+1))
pheatmap(scale(log(mat$RNA+1)), cluster_rows = FALSE)


write.table(scale(log(mat$RNA+1)), file = "ENDmarkergene-average.csv", sep = ",")






gene = c("CEACAM6", "AGER", "TOP2A","MKI67", "SFTPC", "SCGB3A2","HHIP", "SERPINB3", "KRT5", "KRT17", "TP63", "LTF", "PRR4", 
	  "SCGB3A1",  "SCGB1A1", "MUC5B", "MUC5AC","FOXI1","FOXJ1", "EPPIN")


data <- EPI[gene,]
Idents(data) <- EPI@meta.data$EPIClass
mat <- AverageExpression(data)

pheatmap(log(mat$RNA+1))
pheatmap(log(mat$RNA+1), cluster_rows = FALSE)
pheatmap(scale(log(mat$RNA+1)), cluster_rows = FALSE)


write.table(scale(log(mat$RNA+1)), file = "EPImarkergene-average.csv", sep = ",")





gene = c("CEACAM6", "AGER", "SFTPC", "SFTPD")


data <- EPI[gene,]
Idents(data) <- EPI@meta.data$EPIClass
data <- subset(data, idents = c("AT1", "AT2"))
mat <- AverageExpression(data)

pheatmap(log(mat$RNA+1))
pheatmap(log(mat$RNA+1), cluster_rows = FALSE)
pheatmap(scale(log(mat$RNA+1)), cluster_rows = FALSE)


write.table(scale(log(mat$RNA+1)), file = "EPIALVmarkergene-average.csv", sep = ",")





gene = c( "KRT5", "KRT17", "TP63","TOP2A", "SERPINB3", "SCGB3A2", "HHIP", "FOXJ1", "EPPIN", "LTF", "PRR4","FOXI1", "MUC5B", "MUC5AC")


data <- EPI[gene,]
Idents(data) <- EPI@meta.data$EPIClass
data <- subset(data, idents = c("Serous", "Club", "ProliferatingEpithelia", "Cil-px", "Basal-px", "Basal", "Basal-d", "Mucous", "Goblet", "Ionocyte", "Cilia"))
mat <- AverageExpression(data)

pheatmap(log(mat$RNA+1))
pheatmap(log(mat$RNA+1), cluster_rows = FALSE)
pheatmap(scale(log(mat$RNA+1)), cluster_rows = FALSE)


write.table(scale(log(mat$RNA+1)), file = "EPIBROmarkergene-average.csv", sep = ",")









gene = c("CD79A", "MS4A1", "CD27",  "NKG7", "FCER1G", "TYROBP","CD8A", "CD3G", "IL7R")


data <- LYM[gene,]
Idents(data) <- LYM@meta.data$LYMClass
mat <- AverageExpression(data)

pheatmap(log(mat$RNA+1))
pheatmap(log(mat$RNA+1), cluster_rows = FALSE)
pheatmap(scale(log(mat$RNA+1)), cluster_rows = FALSE)


write.table(scale(log(mat$RNA+1)), file = "LYMmarkergene-average.csv", sep = ",")








gene = c("TPSB2", "MS4A2", "KIT",  "CD68", "MARCO", "IL3RA", "CD1C", "CD14", "THBS1", "VCAN", "IFITM2")


data <- MYE[gene,]
Idents(data) <- MYE@meta.data$MYEClass
mat <- AverageExpression(data)

pheatmap(log(mat$RNA+1))
pheatmap(log(mat$RNA+1), cluster_rows = FALSE)
pheatmap(scale(log(mat$RNA+1)), cluster_rows = FALSE)


write.table(scale(log(mat$RNA+1)), file = "MYEmarkergene-average.csv", sep = ",")















##################################
##################################










