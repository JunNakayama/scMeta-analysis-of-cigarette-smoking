library(Seurat)
library(harmony)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(MAST)
library(DoubletFinder)


## Inport gene list from Module candidate
list2 <- read_csv("genelist2.csv")
list2 = unlist(list2)
allgene = rownames(D)
list2 = list(allgene[allgene %in% unlist(list2)])

Res = list()
EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
Vhoge <- VlnPlot(EPImod, features = "Module.1", group.by = "EPIClass", split.by = "Smoking", pt.size = 0)
Res <- bind_rows(Res, Vhoge$data)

STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
Vhoge <- VlnPlot(STRmod, features = "Module.1", group.by = "STRClass", split.by = "Smoking", pt.size = 0)
Res <- bind_rows(Res, Vhoge$data)

ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
Vhoge <- VlnPlot(ENDmod, features = "Module.1", group.by = "ENDClass", split.by = "Smoking", pt.size = 0)
Res <- bind_rows(Res, Vhoge$data)

LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
Vhoge <- VlnPlot(LYMmod, features = "Module.1", group.by = "LYMClass", split.by = "Smoking", pt.size = 0)
Res <- bind_rows(Res, Vhoge$data)

MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")
Vhoge <- VlnPlot(MYEmod, features = "Module.1", group.by = "MYEClass", split.by = "Smoking", pt.size = 0)
Res <- bind_rows(Res, Vhoge$data)

write.table(Res, file = "Mod-melt-MODULE.csv", sep = ",")



## Module Mean extraction
list2 <- read_csv("genelist2.csv")
list2 = unlist(list2)
allgene = rownames(D)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")
ResDF <- list()

Cell <- unique(EPImod@meta.data$EPIClass)
Cell <- Cell[-14]
for(i in 1:length(Cell)){
	C1 <- Cell[i]
	mod1 <- subset(EPImod, subset = EPIClass == C1)

	a <- mean(subset(mod1, subset = Smoking == "Never")@meta.data$Module.1)
	b <- mean(subset(mod1, subset = Smoking == "Smoker")@meta.data$Module.1)

	Res <- data.frame(Cell = C1, Never = a, Smoker = b)
	ResDF <- rbind(ResDF, Res)
}
Cell <- unique(STRmod@meta.data$STRClass)
for(i in 1:length(Cell)){
	C1 <- Cell[i]
	mod1 <- subset(STRmod, subset = STRClass == C1)

	a <- mean(subset(mod1, subset = Smoking == "Never")@meta.data$Module.1)
	b <- mean(subset(mod1, subset = Smoking == "Smoker")@meta.data$Module.1)

	Res <- data.frame(Cell = C1, Never = a, Smoker = b)
	ResDF <- rbind(ResDF, Res)
}
Cell <- unique(ENDmod@meta.data$ENDClass)
for(i in 1:length(Cell)){
	C1 <- Cell[i]
	mod1 <- subset(ENDmod, subset = ENDClass == C1)

	a <- mean(subset(mod1, subset = Smoking == "Never")@meta.data$Module.1)
	b <- mean(subset(mod1, subset = Smoking == "Smoker")@meta.data$Module.1)

	Res <- data.frame(Cell = C1, Never = a, Smoker = b)
	ResDF <- rbind(ResDF, Res)
}
Cell <- unique(LYMmod@meta.data$LYMClass)
for(i in 1:length(Cell)){
	C1 <- Cell[i]
	mod1 <- subset(LYMmod, subset = LYMClass == C1)

	a <- mean(subset(mod1, subset = Smoking == "Never")@meta.data$Module.1)
	b <- mean(subset(mod1, subset = Smoking == "Smoker")@meta.data$Module.1)

	Res <- data.frame(Cell = C1, Never = a, Smoker = b)
	ResDF <- rbind(ResDF, Res)
}
Cell <- unique(MYEmod@meta.data$MYEClass)
for(i in 1:length(Cell)){
	C1 <- Cell[i]
	mod1 <- subset(MYEmod, subset = MYEClass == C1)

	a <- mean(subset(mod1, subset = Smoking == "Never")@meta.data$Module.1)
	b <- mean(subset(mod1, subset = Smoking == "Smoker")@meta.data$Module.1)

	Res <- data.frame(Cell = C1, Never = a, Smoker = b)
	ResDF <- rbind(ResDF, Res)
}

write.table(ResDF, file = "Mod-MeanMatrix-MODULE.csv", sep = ",")
