library(Seurat)
library(harmony)
library(ComplexHeatmap)
library(scater)
library(Matrix)
library(tidyverse)
library(dplyr)
library(tidyr)
library(reshape2)
library(pheatmap)
library(Rtsne)
library(cowplot)
library(ggplot2)
library(ggsci)
library(scales)
library(MAST)
library(DOSE)
library(patchwork)
library(plotly)
library(monocle)
library(MASS)
library(loomR)
library(RColorBrewer)
library(grDevices)
library(colorRamps)
library(data.table)
library(hexbin)


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

## Integration of meta-data
metadata <- read.csv("metadata2.csv")
rownames(metadata) <- metadata[,1]
D@meta.data <- metadata
cl = c("Non-smoker", "Smoker", "COPD")
D <- subset(D, subset = Class %in% cl)

## Normalization and selection
D <- NormalizeData(D)
D[["percent.mt2"]] <- PercentageFeatureSet(D, pattern = c("MT-"))
D <- subset(D, subset = nFeature_RNA > 1000 & percent.mt2 < 20)
D <- FindVariableFeatures(D)
D <- ScaleData(D, vars.to.regress = "percent.mt2")
D <- RunPCA(D, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
D <- RunHarmony(D, "orig.ident")
D <- RunUMAP(D, reduction = "harmony", dims = 1:100, min.dist = 0.5)
D <- FindNeighbors(D, dims = 1:15, reduction = "harmony")
D <- FindClusters(D, resolution = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
D <- CellCycleScoring(D, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Idents(D) <- D@meta.data$seurat_clusters

saveRDS(D, file = "D.rds")




DimPlot(D, reduction = "umap", group.by = "Class4", label = TRUE) + NoLegend()
DimPlot(D, reduction = "umap", group.by = "Class38", label = TRUE) + NoLegend()



write.table(table(D@active.ident, D@meta.data$Class38), file = "Table-ident-Class38.csv", sep = ",")
write.table(table(D@meta.data$Smoking, D@meta.data$Class38), file = "Table-Smo-Class38.csv", sep = ",")





#### Visualization
hoge <- DimPlot(D, reduction = "umap") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$ident)
data <- data[order(data$ident, decreasing = FALSE),]

ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	NoLegend()

hoge <- DimPlot(D, reduction = "umap", group.by = "Class4") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Class4)
data <- data[order(data$ident, decreasing = FALSE),]

ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_manual(values = c("firebrick", "mediumblue", "darkorchid", "skyblue2", "magenta2", "limegreen")) + 
	NoLegend()

hoge <- DimPlot(D, reduction = "umap", group.by = "Class38") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Class38)
data <- data[order(data$ident, decreasing = FALSE),]

ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic()  + 
	NoLegend()

## Density plot
ggplot(data, mapping = aes(UMAP1, UMAP2)) +
  geom_hex(bins = 80) +
  scale_fill_continuous(type = "viridis") +
  theme_classic() +
  geom_density2d(size = 0.4, bins = 10, colour = "red2")



## Smoking umap
hoge <- DimPlot(D, reduction = "umap", group.by = "Smoking") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Smoking)
data <- data[order(data$ident, decreasing = FALSE),]

ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_manual(values = c("skyblue2", "red3")) + 
	NoLegend()


hoge <- DimPlot(subset(D, subset = Smoking == "Never"), reduction = "umap", group.by = "Smoking") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Smoking)
data <- data[order(data$ident, decreasing = FALSE),]

ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_manual(values = c("skyblue2")) + 
	NoLegend()


hoge <- DimPlot(subset(D, subset = Smoking == "Smoker"), reduction = "umap", group.by = "Smoking") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Smoking)
data <- data[order(data$ident, decreasing = FALSE),]

ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_manual(values = c("red3")) + 
	NoLegend()


## Dataset umap

hoge <- DimPlot(D, reduction = "umap", group.by = "Dataset") + NoLegend()
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data$Dataset)
data <- data[order(data$ident, decreasing = FALSE),]

ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 0.5, size = 0.25) +
	theme_classic() +
	scale_color_brewer(palette = "Dark2") + 
	NoLegend()



## Feature visulaization
hoge <- FeaturePlot(D, features = "EPCAM")
data <- data.frame(UMAP1 = hoge$data$UMAP_1, UMAP2 = hoge$data$UMAP_2, ident = hoge$data[,4])
data <- data[order(data$ident, decreasing = FALSE),]

ggplot(data, aes(UMAP1, UMAP2, color = ident)) +
	geom_point(position = "jitter", alpha = 1, size = 0.25) +
	theme_classic() +
	scale_color_gradient(low = "grey", high = "mediumblue") + 
	NoLegend()






#### Separation of dataset in each cell type

