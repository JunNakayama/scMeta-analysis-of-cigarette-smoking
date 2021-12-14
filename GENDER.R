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

#### GENDER analysis
#### Only Smoker


GEN = EPI
Res = list()
	GEN = subset(GEN, subset = Gender %in% c("Male", "Female"))
	GEN = subset(GEN, subset = Smoking == c("Smoker"))

n = unique(GEN@meta.data$EPIClass)

for(i in 1:length(n)){
	c <- n[i]
	CELL <- subset(GEN, subset = EPIClass == c)

	Idents(CELL) <- CELL@meta.data$Gender

	markers <- FindAllMarkers(CELL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST")
	markers <- data.frame(markers, Cellident = c)
	Res <- rbind(Res, markers)
}
Res <- Res[Res$p_val_adj < 0.05,]

write.table(Res, "GENDERmarkers-EPI.csv", sep = ",")


