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




#### Separation of dataset into each cell type
## Epithelia
EPI <- subset(D, subset = Class4 %in% c("Epithelia", "ProliferatingEpithelia"))
EPI <- RunUMAP(EPI, reduction = "harmony", dims = 1:100, min.dist = 0.5)
EPI <- FindNeighbors(EPI, dims = 1:15, reduction = "harmony")
EPI <- FindClusters(EPI, resolution = 1)

DimPlot(EPI, label = TRUE, label.size = 10) + NoLegend()
DimPlot(EPI, group.by = "Smoking") + NoLegend()

EPI.markers <- FindAllMarkers(EPI, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST")
write.table(EPI.markers, file = 'EPIMarker-genes.tsv', sep='	')

DimPlot(EPI, reduction = "umap", group.by = "EPIClass", label = TRUE) + NoLegend()

saveRDS(EPI, file = "EPI.rds")



# Basal-px
Bpx <- subset(EPI, subset = EPIClass == c("Basal-px"))
Bpx <- RunUMAP(Bpx, reduction = "harmony", dims = 1:100, min.dist = 0.5)
Bpx <- FindNeighbors(Bpx, dims = 1:15, reduction = "harmony")
Bpx <- FindClusters(Bpx, resolution = 1)

colclass = c("deepskyblue3", "red3")
DimPlot(Bpx, label = TRUE, label.size = 10) + NoLegend()
DimPlot(Bpx, group.by = "Smoking", pt.size = 1.25, cols = colclass) + NoLegend()

Idents(Bpx) <- Bpx@meta.data$Smoking
Bpx.markers <- FindAllMarkers(Bpx, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST")
write.table(Bpx.markers, file = 'BpxMarker-genes.csv', sep=',')



## Immune-cell
IMM <- subset(D, subset = Class4 %in% c("ImmuneCell", "ProliferatingImmunecell"))
IMM <- RunUMAP(IMM, reduction = "harmony", dims = 1:100, min.dist = 0.5)
IMM <- FindNeighbors(IMM, dims = 1:15, reduction = "harmony")
IMM <- FindClusters(IMM, resolution = 1)

DimPlot(IMM, label = TRUE, label.size = 10) + NoLegend()
DimPlot(IMM, group.by = "Smoking") + NoLegend()

IMM.markers <- FindAllMarkers(IMM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST")
write.table(IMM.markers, file = 'IMMMarker-genes.tsv', sep='	')

saveRDS(IMM, file = "IMM.rds")

## Dividing into Lym and Mye
LYM <- subset(IMM, ident = c(1,7,14,5,18))
MYE <- subset(IMM, ident = c(2,9,13,19,3,6,4,15,0,20,10,12,11,16,22))
LYM <- RunUMAP(LYM, reduction = "harmony", dims = 1:100, min.dist = 0.5)
LYM <- FindNeighbors(LYM, dims = 1:15, reduction = "harmony")
LYM <- FindClusters(LYM, resolution = 1)

DimPlot(LYM, label = TRUE, label.size = 10) + NoLegend()
DimPlot(LYM, group.by = "Smoking") + NoLegend()
DimPlot(LYM, reduction = "umap", group.by = "LYMClass", label = TRUE) + NoLegend()


LYM.markers <- FindAllMarkers(LYM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST")
write.table(LYM.markers, file = 'LYMMarker-genes.tsv', sep='	')

saveRDS(LYM, file = "LYM.rds")



MYE <- RunUMAP(MYE, reduction = "harmony", dims = 1:100, min.dist = 0.5)
MYE <- FindNeighbors(MYE, dims = 1:15, reduction = "harmony")
MYE <- FindClusters(MYE, resolution = 1)

DimPlot(MYE, label = TRUE, label.size = 10) + NoLegend()
DimPlot(MYE, group.by = "Smoking") + NoLegend()
DimPlot(MYE, reduction = "umap", group.by = "MYEClass", label = TRUE) + NoLegend()

MYE.markers <- FindAllMarkers(MYE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST")
write.table(MYE.markers, file = 'MYEMarker-genes.tsv', sep='	')

saveRDS(MYE, file = "MYE.rds")




## Fibroblastic
STR <- subset(D, subset = Class4 == c("Stroma"))
STR <- RunUMAP(STR, reduction = "harmony", dims = 1:100, min.dist = 0.5)
STR <- FindNeighbors(STR, dims = 1:15, reduction = "harmony")
STR <- FindClusters(STR, resolution = 1)

DimPlot(STR, label = TRUE, label.size = 10) + NoLegend()
DimPlot(STR, group.by = "Smoking") + NoLegend()

STR.markers <- FindAllMarkers(STR, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST")
write.table(STR.markers, file = 'STRMarker-genes.tsv', sep='	')

DimPlot(STR, reduction = "umap", group.by = "STRClass", label = TRUE) + NoLegend()

saveRDS(STR, file = "STR.rds")




# Adventital Fibroblast ## ACTA2
AdvF <- subset(STR, subset = STRClass == c("AdvFibroblast"))
AdvF <- RunUMAP(AdvF, reduction = "harmony", dims = 1:100, min.dist = 0.5)
colclass = c("deepskyblue3", "red3")
DimPlot(AdvF, label = TRUE, label.size = 10) + NoLegend()
DimPlot(AdvF, group.by = "Smoking", pt.size = 1.25, cols = colclass) + NoLegend()
VlnPlot(AdvF, features = "ACTA2", cols = colclass, group.by = "Smoking")

AdvF <- subset(AdvF, subset = Smoking == c("Smoker"))
IDENT <- AdvF@assays$RNA@scale.data["ACTA2",]
IDENT[IDENT > 0] <- "high"
IDENT[IDENT <= 0] <- "low"
Idents(AdvF) <- IDENT
AdvF.markers <- FindAllMarkers(AdvF, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST")
write.table(AdvF.markers, file = 'AdvF-ACTA2-HandL-Marker-genes.csv', sep=',')

gene <- AdvF.markers$gene
DotPlot(AdvF, features = gene, cols = colclass, dot.scale = 12) + RotatedAxis() + coord_flip()

data <- VlnPlot(AdvF, features = "ACTA2")
write.table(data$data, file = "Vln-STRAdvFibSmo-ACTA2.csv", sep = ",")


## Endothelia
END <- subset(D, subset = Class4 == c("Endothelium"))
END <- RunUMAP(END, reduction = "harmony", dims = 1:100, min.dist = 0.5)
END <- FindNeighbors(END, dims = 1:15, reduction = "harmony")
END <- FindClusters(END, resolution = 1)

DimPlot(END, label = TRUE, label.size = 10) + NoLegend()
DimPlot(END, group.by = "Smoking") + NoLegend()

END.markers <- FindAllMarkers(END, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST")
write.table(END.markers, file = 'ENDMarker-genes.tsv', sep='	')

DimPlot(END, reduction = "umap", group.by = "ENDClass", label = TRUE) + NoLegend()

saveRDS(END, file = "END.rds")



## Visualization by violin plot
colclass = c("deepskyblue1", "red3")
VlnPlot(END, features = "COL18A1", split.by = "Smoking", group.by = "ENDClass", pt.size = 0, cols = colclass)

data <- VlnPlot(END, features = "COL18A1", split.by = "Smoking", group.by = "ENDClass", pt.size = 0, cols = colclass)
write.table(data$data, file = "Vln-END-COL18A1.csv", sep = ",")




# Lymphatic ## ANTGPT2
LymV <- subset(END, subset = ENDClass == c("Lymphatic"))
LymV <- RunUMAP(LymV, reduction = "harmony", dims = 1:100, min.dist = 0.5)

colclass = c("deepskyblue3", "red3")
DimPlot(LymV, label = TRUE, label.size = 10) + NoLegend()
DimPlot(LymV, group.by = "Smoking", pt.size = 1.25, cols = colclass) + NoLegend()
VlnPlot(LymV, features = "ANGPT2", cols = colclass, group.by = "Smoking")

LymV <- subset(LymV, subset = Smoking == c("Smoker"))
IDENT <- LymV@assays$RNA@scale.data["ANGPT2",]
IDENT[IDENT > 0] <- "high"
IDENT[IDENT <= 0] <- "low"
Idents(LymV) <- IDENT
LymV.markers <- FindAllMarkers(LymV, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST")
write.table(LymV.markers, file = 'Lymphatic-ANGPT2-HandL-Marker-genes.csv', sep=',')

gene <- LymV.markers$gene
DotPlot(LymV, features = gene, cols = colclass, dot.scale = 12) + RotatedAxis() + coord_flip()

data <- VlnPlot(LymV, features = "ANGPT2")
write.table(data$data, file = "Vln-ENDLymSmo-ANGPT2.csv", sep = ",")



## Marker visualization
# Marker Matrix for Heatmap
gene = c("VEGFD", "FGFR4", "PDFDRA", "DES", "PDGFRB",  "MSLN", "UPK3B", "MKI67", "SFRP2", "PDGFRL", "ASPN")
data <- STR[gene,]
Idents(data) <- STR@meta.data$STRClass
mat <- AverageExpression(data)
pheatmap(log(mat$RNA+1))
pheatmap(scale(log(mat$RNA+1)), cluster_rows = FALSE)
pheatmap(scale(log(mat$RNA+1)), cluster_rows = FALSE, color = c("cyan","red", "red2", "red4"))

write.table(scale(log(mat$RNA+1)), file = "STRmarkergene-average.csv", sep = ",")

gene = c("ACKR1", "CPE", "IL7R", "SLC6A4", "CA4", "GJA5", "S100A3", "EDNRB", "PROX1")
data <- END[gene,]
Idents(data) <- END@meta.data$ENDClass
mat <- AverageExpression(data)
pheatmap(log(mat$RNA+1))
pheatmap(scale(log(mat$RNA+1)), cluster_rows = FALSE)
pheatmap(scale(log(mat$RNA+1)), cluster_rows = FALSE, color = c("cyan", "red2", "red4"))

write.table(scale(log(mat$RNA+1)), file = "ENDmarkergene-average.csv", sep = ",")

gene = c("SFTPC", "CEACAM6", "AGER", "KRT5", "KRT17","SERPINB3", "MKI67", "PCNA", "SCGB3A2", "SCGB3A1", "SCGB1A1", "LTF", "PRR4", "MUC5B", "MUC5AC",
		 "AQP5", "FOXJ1", "MUC12", "EPPIN","FOXI1")
data <- EPI[gene,]
Idents(data) <- EPI@meta.data$EPIClass
mat <- AverageExpression(data)
pheatmap(log(mat$RNA+1))
pheatmap(scale(log(mat$RNA+1)), cluster_rows = FALSE)
pheatmap(scale(log(mat$RNA+1)), cluster_rows = FALSE, color = c("cyan", "red2", "red4"))

write.table(scale(log(mat$RNA+1)), file = "EPImarkergene-average.csv", sep = ",")

gene = c("COTL1", "LDHB", "CD3G", "LEF1", "IL7R", "CD8A", "FCER1G", "MS4A1", "CD27", "NKG7", "KLRD1", "GZMB", "GZMK", "GZMH")
data <- LYM[gene,]
Idents(data) <- LYM@meta.data$LYMClass
mat <- AverageExpression(data)
pheatmap(log(mat$RNA+1))
pheatmap(scale(log(mat$RNA+1)), cluster_rows = FALSE)
pheatmap(scale(log(mat$RNA+1)), cluster_rows = FALSE, color = c("cyan", "red2", "red4"))

write.table(scale(log(mat$RNA+1)), file = "LYMmarkergene-average.csv", sep = ",")

gene = c("MS4A2", "IL3RA", "KIT", "CD1C", "CD14", "FCGR3A", "S100A8", "S100A9", "IFITM2", "CD68", "MARCO", "LPL", "KAZN")
data <- MYE[gene,]
Idents(data) <- MYE@meta.data$MYEClass
mat <- AverageExpression(data)
pheatmap(log(mat$RNA+1))
pheatmap(scale(log(mat$RNA+1)), cluster_rows = FALSE)
pheatmap(scale(log(mat$RNA+1)), cluster_rows = FALSE, color = c("cyan", "red2", "red4"))

write.table(scale(log(mat$RNA+1)), file = "MYEmarkergene-average.csv", sep = ",")














