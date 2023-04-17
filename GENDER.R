
## This analysis was removed from the revised manuscript






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




## Enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ReactomePA)

h.all.genes = rownames(EPI)
h.all.genes.entrez = bitr(h.all.genes, fromType="SYMBOL", 
                         toType="ENTREZID", OrgDb="org.Hs.eg.db")
h.all.genes.entrez = h.all.genes.entrez[,2]

pmarker.entrez =  bitr(Res$gene, fromType="SYMBOL", 
                         toType="ENTREZID", OrgDb="org.Hs.eg.db")
clus = data.frame(Res$gene, Res$cluster, Res$Cellident)
pmarkers.list = merge(pmarker.entrez, clus, by = 1)
sortlist = order(pmarkers.list[,3])
pmarkers.list = pmarkers.list[sortlist,]




# Pathway
nd = unique(pmarkers.list$Res.cluster)
ent.list = NULL
PathRes = list()
for(i in 1:length(nd)){
	cl <- nd[i]
	ent.list = pmarkers.list[pmarkers.list$Res.cluster == cl, ]
	ent.list = ent.list[ent.list$Res.Cellident == "Basal", ]
	ent.list = ent.list[,2]

	GOenrich <- enrichPathway(gene = ent.list, universe = h.all.genes.entrez, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = T)
	PathRes[i] = GOenrich

}
names(PathRes) = paste(unique(nd), "Basal")
Path <- list()
cl <- nd[1]
ho <- as.data.frame(PathRes[cl])
cname <- colnames(ho)
for(i in 1:length(nd)){
	cl <- nd[i]
	ho <- as.data.frame(PathRes[cl])
	colnames(ho) <- cname
	if(length(ho) == 0){
		print(i)
		}else{
	ho <- mutate(ho, i, cl)
	Path <- rbind(Path, ho)
	}
}
write.table(Path, file = "Reactome-GEN-EPI-Basal-pathway.csv", sep = ",")

nd = unique(pmarkers.list$Res.cluster)
ent.list = NULL
PathRes = list()
for(i in 1:length(nd)){
	cl <- nd[i]
	ent.list = pmarkers.list[pmarkers.list$Res.cluster == cl, ]
	ent.list = ent.list[ent.list$Res.Cellident == "Cilia", ]
	ent.list = ent.list[,2]

	GOenrich <- enrichPathway(gene = ent.list, universe = h.all.genes.entrez, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = T)
	PathRes[i] = GOenrich

}
names(PathRes) = paste(unique(nd), "Cilia")
Path <- list()
cl <- nd[1]
ho <- as.data.frame(PathRes[cl])
cname <- colnames(ho)
for(i in 1:length(nd)){
	cl <- nd[i]
	ho <- as.data.frame(PathRes[cl])
	colnames(ho) <- cname
	if(length(ho) == 0){
		print(i)
		}else{
	ho <- mutate(ho, i, cl)
	Path <- rbind(Path, ho)
	}
}
write.table(Path, file = "Reactome-GEN-EPI-Cilia-pathway.csv", sep = ",")

nd = unique(pmarkers.list$Res.cluster)
ent.list = NULL
PathRes = list()
for(i in 1:length(nd)){
	cl <- nd[i]
	ent.list = pmarkers.list[pmarkers.list$Res.cluster == cl, ]
	ent.list = ent.list[ent.list$Res.Cellident == "TracheaBasal", ]
	ent.list = ent.list[,2]

	GOenrich <- enrichPathway(gene = ent.list, universe = h.all.genes.entrez, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = T)
	PathRes[i] = GOenrich

}
names(PathRes) = paste(unique(nd), "TracheaBasal")
Path <- list()
cl <- nd[1]
ho <- as.data.frame(PathRes[cl])
cname <- colnames(ho)
for(i in 1:length(nd)){
	cl <- nd[i]
	ho <- as.data.frame(PathRes[cl])
	colnames(ho) <- cname
	if(length(ho) == 0){
		print(i)
		}else{
	ho <- mutate(ho, i, cl)
	Path <- rbind(Path, ho)
	}
}
write.table(Path, file = "Reactome-GEN-EPI-TracheaBasal-pathway.csv", sep = ",")

nd = unique(pmarkers.list$Res.cluster)
ent.list = NULL
PathRes = list()
for(i in 1:length(nd)){
	cl <- nd[i]
	ent.list = pmarkers.list[pmarkers.list$Res.cluster == cl, ]
	ent.list = ent.list[ent.list$Res.Cellident == "TracheaBasal-px", ]
	ent.list = ent.list[,2]

	GOenrich <- enrichPathway(gene = ent.list, universe = h.all.genes.entrez, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = T)
	PathRes[i] = GOenrich

}
names(PathRes) = paste(unique(nd), "TracheaBasal-px")
Path <- list()
cl <- nd[1]
ho <- as.data.frame(PathRes[cl])
cname <- colnames(ho)
for(i in 1:length(nd)){
	cl <- nd[i]
	ho <- as.data.frame(PathRes[cl])
	colnames(ho) <- cname
	if(length(ho) == 0){
		print(i)
		}else{
	ho <- mutate(ho, i, cl)
	Path <- rbind(Path, ho)
	}
}
write.table(Path, file = "Reactome-GEN-EPI-TracheaBasal-px-pathway.csv", sep = ",")

nd = unique(pmarkers.list$Res.cluster)
ent.list = NULL
PathRes = list()

for(i in 1:length(nd)){
	cl <- nd[i]
	ent.list = pmarkers.list[pmarkers.list$Res.cluster == cl, ]
	ent.list = ent.list[ent.list$Res.Cellident == "AT2", ]
	ent.list = ent.list[,2]

	GOenrich <- enrichPathway(gene = ent.list, universe = h.all.genes.entrez, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = T)
	PathRes[i] = GOenrich

}
names(PathRes) = paste(unique(nd), "AT2")
Path <- list()
cl <- nd[1]
ho <- as.data.frame(PathRes[cl])
cname <- colnames(ho)
for(i in 1:length(nd)){
	cl <- nd[i]
	ho <- as.data.frame(PathRes[cl])
	colnames(ho) <- cname
	if(length(ho) == 0){
		print(i)
		}else{
	ho <- mutate(ho, i, cl)
	Path <- rbind(Path, ho)
	}
}
write.table(Path, file = "Reactome-GEN-EPI-AT2-pathway.csv", sep = ",")

nd = unique(pmarkers.list$Res.cluster)
ent.list = NULL
PathRes = list()

for(i in 1:length(nd)){
	cl <- nd[i]
	ent.list = pmarkers.list[pmarkers.list$Res.cluster == cl, ]
	ent.list = ent.list[ent.list$Res.Cellident == "AT1", ]
	ent.list = ent.list[,2]

	GOenrich <- enrichPathway(gene = ent.list, universe = h.all.genes.entrez, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = T)
	PathRes[i] = GOenrich

}
names(PathRes) = paste(unique(nd), "AT1")
Path <- list()
cl <- nd[1]
ho <- as.data.frame(PathRes[cl])
cname <- colnames(ho)
for(i in 1:length(nd)){
	cl <- nd[i]
	ho <- as.data.frame(PathRes[cl])
	colnames(ho) <- cname
	if(length(ho) == 0){
		print(i)
		}else{
	ho <- mutate(ho, i, cl)
	Path <- rbind(Path, ho)
	}
}
write.table(Path, file = "Reactome-GEN-EPI-AT1-pathway.csv", sep = ",")

nd = unique(pmarkers.list$Res.cluster)
ent.list = NULL
PathRes = list()
for(i in 1:length(nd)){
	cl <- nd[i]
	ent.list = pmarkers.list[pmarkers.list$Res.cluster == cl, ]
	ent.list = ent.list[ent.list$Res.Cellident == "Goblet", ]
	ent.list = ent.list[,2]

	GOenrich <- enrichPathway(gene = ent.list, universe = h.all.genes.entrez, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = T)
	PathRes[i] = GOenrich

}
names(PathRes) = paste(unique(nd), "Goblet")
Path <- list()
cl <- nd[1]
ho <- as.data.frame(PathRes[cl])
cname <- colnames(ho)
for(i in 1:length(nd)){
	cl <- nd[i]
	ho <- as.data.frame(PathRes[cl])
	colnames(ho) <- cname
	if(length(ho) == 0){
		print(i)
		}else{
	ho <- mutate(ho, i, cl)
	Path <- rbind(Path, ho)
	}
}
write.table(Path, file = "Reactome-GEN-EPI-Goblet-pathway.csv", sep = ",")

nd = unique(pmarkers.list$Res.cluster)
ent.list = NULL
PathRes = list()
for(i in 1:length(nd)){
	cl <- nd[i]
	ent.list = pmarkers.list[pmarkers.list$Res.cluster == cl, ]
	ent.list = ent.list[ent.list$Res.Cellident == "ProliferatingEpithelia", ]
	ent.list = ent.list[,2]

	GOenrich <- enrichPathway(gene = ent.list, universe = h.all.genes.entrez, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = T)
	PathRes[i] = GOenrich

}
names(PathRes) = paste(unique(nd), "ProliferatingEpithelia")
Path <- list()
cl <- nd[1]
ho <- as.data.frame(PathRes[cl])
cname <- colnames(ho)
for(i in 1:length(nd)){
	cl <- nd[i]
	ho <- as.data.frame(PathRes[cl])
	colnames(ho) <- cname
	if(length(ho) == 0){
		print(i)
		}else{
	ho <- mutate(ho, i, cl)
	Path <- rbind(Path, ho)
	}
}
write.table(Path, file = "Reactome-GEN-EPI-ProliferatingEpithelia-pathway.csv", sep = ",")

nd = unique(pmarkers.list$Res.cluster)
ent.list = NULL
PathRes = list()
for(i in 1:length(nd)){
	cl <- nd[i]
	ent.list = pmarkers.list[pmarkers.list$Res.cluster == cl, ]
	ent.list = ent.list[ent.list$Res.Cellident == "Basal-px", ]
	ent.list = ent.list[,2]

	GOenrich <- enrichPathway(gene = ent.list, universe = h.all.genes.entrez, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = T)
	PathRes[i] = GOenrich

}
names(PathRes) = paste(unique(nd), "Basal-px")
Path <- list()
cl <- nd[1]
ho <- as.data.frame(PathRes[cl])
cname <- colnames(ho)
for(i in 1:length(nd)){
	cl <- nd[i]
	ho <- as.data.frame(PathRes[cl])
	colnames(ho) <- cname
	if(length(ho) == 0){
		print(i)
		}else{
	ho <- mutate(ho, i, cl)
	Path <- rbind(Path, ho)
	}
}
write.table(Path, file = "Reactome-GEN-EPI-Basal-px-pathway.csv", sep = ",")

nd = unique(pmarkers.list$Res.cluster)
ent.list = NULL
PathRes = list()
for(i in 1:length(nd)){
	cl <- nd[i]
	ent.list = pmarkers.list[pmarkers.list$Res.cluster == cl, ]
	ent.list = ent.list[ent.list$Res.Cellident == "Cil-px", ]
	ent.list = ent.list[,2]

	GOenrich <- enrichPathway(gene = ent.list, universe = h.all.genes.entrez, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = T)
	PathRes[i] = GOenrich

}
names(PathRes) = paste(unique(nd), "Cil-px")
Path <- list()
cl <- nd[1]
ho <- as.data.frame(PathRes[cl])
cname <- colnames(ho)
for(i in 1:length(nd)){
	cl <- nd[i]
	ho <- as.data.frame(PathRes[cl])
	colnames(ho) <- cname
	if(length(ho) == 0){
		print(i)
		}else{
	ho <- mutate(ho, i, cl)
	Path <- rbind(Path, ho)
	}
}
write.table(Path, file = "Reactome-GEN-EPI-Cil-px-pathway.csv", sep = ",")

nd = unique(pmarkers.list$Res.cluster)
ent.list = NULL
PathRes = list()
for(i in 1:length(nd)){
	cl <- nd[i]
	ent.list = pmarkers.list[pmarkers.list$Res.cluster == cl, ]
	ent.list = ent.list[ent.list$Res.Cellident == "Club", ]
	ent.list = ent.list[,2]

	GOenrich <- enrichPathway(gene = ent.list, universe = h.all.genes.entrez, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = T)
	PathRes[i] = GOenrich

}
names(PathRes) = paste(unique(nd), "Club")
Path <- list()
cl <- nd[1]
ho <- as.data.frame(PathRes[cl])
cname <- colnames(ho)
for(i in 1:length(nd)){
	cl <- nd[i]
	ho <- as.data.frame(PathRes[cl])
	colnames(ho) <- cname
	if(length(ho) == 0){
		print(i)
		}else{
	ho <- mutate(ho, i, cl)
	Path <- rbind(Path, ho)
	}
}
write.table(Path, file = "Reactome-GEN-EPI-Club-pathway.csv", sep = ",")

nd = unique(pmarkers.list$Res.cluster)
ent.list = NULL
PathRes = list()
for(i in 1:length(nd)){
	cl <- nd[i]
	ent.list = pmarkers.list[pmarkers.list$Res.cluster == cl, ]
	ent.list = ent.list[ent.list$Res.Cellident == "Serous", ]
	ent.list = ent.list[,2]

	GOenrich <- enrichPathway(gene = ent.list, universe = h.all.genes.entrez, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = T)
	PathRes[i] = GOenrich

}
names(PathRes) = paste(unique(nd), "Serous")
Path <- list()
cl <- nd[1]
ho <- as.data.frame(PathRes[cl])
cname <- colnames(ho)
for(i in 1:length(nd)){
	cl <- nd[i]
	ho <- as.data.frame(PathRes[cl])
	colnames(ho) <- cname
	if(length(ho) == 0){
		print(i)
		}else{
	ho <- mutate(ho, i, cl)
	Path <- rbind(Path, ho)
	}
}
write.table(Path, file = "Reactome-GEN-EPI-Serous-pathway.csv", sep = ",")

