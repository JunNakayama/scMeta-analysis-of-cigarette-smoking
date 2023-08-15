library(Seurat)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(MAST)



## Inport gene list from Module candidate

list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")


ResDF <- list()
Res <- list()
Res <- data.frame(Cell = EPImod@meta.data$EPIClass, Smoking = EPImod@meta.data$Smoking, EMT = EPImod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = STRmod@meta.data$STRClass, Smoking = STRmod@meta.data$Smoking, EMT = STRmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = ENDmod@meta.data$ENDClass, Smoking = ENDmod@meta.data$Smoking, EMT = ENDmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = LYMmod@meta.data$LYMClass, Smoking = LYMmod@meta.data$Smoking, EMT = LYMmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = MYEmod@meta.data$MYEClass, Smoking = MYEmod@meta.data$Smoking, EMT = MYEmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)

write.table(ResDF, file = "Mod-allrawdata-HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.csv", sep = ",")


 ggplot(ResDF, aes(fill=Smoking, y=EMT, x=Cell)) + 
    geom_violin(position="dodge", alpha=0.8) +
    theme_classic()




list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/HALLMARK_HEME_METABOLISM.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")


ResDF <- list()
Res <- list()
Res <- data.frame(Cell = EPImod@meta.data$EPIClass, Smoking = EPImod@meta.data$Smoking, EMT = EPImod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = STRmod@meta.data$STRClass, Smoking = STRmod@meta.data$Smoking, EMT = STRmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = ENDmod@meta.data$ENDClass, Smoking = ENDmod@meta.data$Smoking, EMT = ENDmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = LYMmod@meta.data$LYMClass, Smoking = LYMmod@meta.data$Smoking, EMT = LYMmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = MYEmod@meta.data$MYEClass, Smoking = MYEmod@meta.data$Smoking, EMT = MYEmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)

write.table(ResDF, file = "Mod-allrawdata-HALLMARK_HEME_METABOLISM.csv", sep = ",")



list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/HOUSTIS_ROS.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")


ResDF <- list()
Res <- list()
Res <- data.frame(Cell = EPImod@meta.data$EPIClass, Smoking = EPImod@meta.data$Smoking, EMT = EPImod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = STRmod@meta.data$STRClass, Smoking = STRmod@meta.data$Smoking, EMT = STRmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = ENDmod@meta.data$ENDClass, Smoking = ENDmod@meta.data$Smoking, EMT = ENDmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = LYMmod@meta.data$LYMClass, Smoking = LYMmod@meta.data$Smoking, EMT = LYMmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = MYEmod@meta.data$MYEClass, Smoking = MYEmod@meta.data$Smoking, EMT = MYEmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)

write.table(ResDF, file = "Mod-allrawdata-HOUSTIS_ROS.csv", sep = ",")



list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/REACTOME_AUTOPHAGY.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")


ResDF <- list()
Res <- list()
Res <- data.frame(Cell = EPImod@meta.data$EPIClass, Smoking = EPImod@meta.data$Smoking, EMT = EPImod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = STRmod@meta.data$STRClass, Smoking = STRmod@meta.data$Smoking, EMT = STRmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = ENDmod@meta.data$ENDClass, Smoking = ENDmod@meta.data$Smoking, EMT = ENDmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = LYMmod@meta.data$LYMClass, Smoking = LYMmod@meta.data$Smoking, EMT = LYMmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = MYEmod@meta.data$MYEClass, Smoking = MYEmod@meta.data$Smoking, EMT = MYEmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)

write.table(ResDF, file = "Mod-allrawdata-REACTOME_AUTOPHAGY.csv", sep = ",")



list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/REACTOME_CELLULAR_SENESCENCE.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")


ResDF <- list()
Res <- list()
Res <- data.frame(Cell = EPImod@meta.data$EPIClass, Smoking = EPImod@meta.data$Smoking, EMT = EPImod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = STRmod@meta.data$STRClass, Smoking = STRmod@meta.data$Smoking, EMT = STRmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = ENDmod@meta.data$ENDClass, Smoking = ENDmod@meta.data$Smoking, EMT = ENDmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = LYMmod@meta.data$LYMClass, Smoking = LYMmod@meta.data$Smoking, EMT = LYMmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = MYEmod@meta.data$MYEClass, Smoking = MYEmod@meta.data$Smoking, EMT = MYEmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)

write.table(ResDF, file = "Mod-allrawdata-REACTOME_CELLULAR_SENESCENCE.csv", sep = ",")



list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/REACTOME_CIRCADIAN_CLOCK.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")


ResDF <- list()
Res <- list()
Res <- data.frame(Cell = EPImod@meta.data$EPIClass, Smoking = EPImod@meta.data$Smoking, EMT = EPImod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = STRmod@meta.data$STRClass, Smoking = STRmod@meta.data$Smoking, EMT = STRmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = ENDmod@meta.data$ENDClass, Smoking = ENDmod@meta.data$Smoking, EMT = ENDmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = LYMmod@meta.data$LYMClass, Smoking = LYMmod@meta.data$Smoking, EMT = LYMmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = MYEmod@meta.data$MYEClass, Smoking = MYEmod@meta.data$Smoking, EMT = MYEmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)

write.table(ResDF, file = "Mod-allrawdata-REACTOME_CIRCADIAN_CLOCK.csv", sep = ",")



list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/REACTOME_INTERFERON_SIGNALING.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")


ResDF <- list()
Res <- list()
Res <- data.frame(Cell = EPImod@meta.data$EPIClass, Smoking = EPImod@meta.data$Smoking, EMT = EPImod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = STRmod@meta.data$STRClass, Smoking = STRmod@meta.data$Smoking, EMT = STRmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = ENDmod@meta.data$ENDClass, Smoking = ENDmod@meta.data$Smoking, EMT = ENDmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = LYMmod@meta.data$LYMClass, Smoking = LYMmod@meta.data$Smoking, EMT = LYMmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = MYEmod@meta.data$MYEClass, Smoking = MYEmod@meta.data$Smoking, EMT = MYEmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)

write.table(ResDF, file = "Mod-allrawdata-REACTOME_INTERFERON_SIGNALING.csv", sep = ",")



list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/REACTOME_MITOPHAGY.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")


ResDF <- list()
Res <- list()
Res <- data.frame(Cell = EPImod@meta.data$EPIClass, Smoking = EPImod@meta.data$Smoking, EMT = EPImod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = STRmod@meta.data$STRClass, Smoking = STRmod@meta.data$Smoking, EMT = STRmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = ENDmod@meta.data$ENDClass, Smoking = ENDmod@meta.data$Smoking, EMT = ENDmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = LYMmod@meta.data$LYMClass, Smoking = LYMmod@meta.data$Smoking, EMT = LYMmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = MYEmod@meta.data$MYEClass, Smoking = MYEmod@meta.data$Smoking, EMT = MYEmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)

write.table(ResDF, file = "Mod-allrawdata-REACTOME_MITOPHAGY.csv", sep = ",")



list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/REACTOME_PYROPTOSIS.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")


ResDF <- list()
Res <- list()
Res <- data.frame(Cell = EPImod@meta.data$EPIClass, Smoking = EPImod@meta.data$Smoking, EMT = EPImod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = STRmod@meta.data$STRClass, Smoking = STRmod@meta.data$Smoking, EMT = STRmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = ENDmod@meta.data$ENDClass, Smoking = ENDmod@meta.data$Smoking, EMT = ENDmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = LYMmod@meta.data$LYMClass, Smoking = LYMmod@meta.data$Smoking, EMT = LYMmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = MYEmod@meta.data$MYEClass, Smoking = MYEmod@meta.data$Smoking, EMT = MYEmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)

write.table(ResDF, file = "Mod-allrawdata-REACTOME_PYROPTOSIS.csv", sep = ",")



list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/WP_FERROPTOSIS.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")


ResDF <- list()
Res <- list()
Res <- data.frame(Cell = EPImod@meta.data$EPIClass, Smoking = EPImod@meta.data$Smoking, EMT = EPImod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = STRmod@meta.data$STRClass, Smoking = STRmod@meta.data$Smoking, EMT = STRmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = ENDmod@meta.data$ENDClass, Smoking = ENDmod@meta.data$Smoking, EMT = ENDmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = LYMmod@meta.data$LYMClass, Smoking = LYMmod@meta.data$Smoking, EMT = LYMmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)
Res <- data.frame(Cell = MYEmod@meta.data$MYEClass, Smoking = MYEmod@meta.data$Smoking, EMT = MYEmod@meta.data$Module.1)
ResDF <- rbind(ResDF, Res)

write.table(ResDF, file = "Mod-allrawdata-WP_FERROPTOSIS.csv", sep = ",")










################################
## Calculation of Mean
Cell <- unique(EPImod@meta.data$EPIClass)
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

write.table(ResDF, file = "Mod-MeanMatrix-HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.csv", sep = ",")



list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/HALLMARK_HEME_METABOLISM.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")
ResDF <- list()

Cell <- unique(EPImod@meta.data$EPIClass)
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

write.table(ResDF, file = "Mod-MeanMatrix-HHALLMARK_HEME_METABOLISM.csv", sep = ",")



list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/HOUSTIS_ROS.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")
ResDF <- list()

Cell <- unique(EPImod@meta.data$EPIClass)
for(i in 1:length(Cell)){
	C1 <- Cell[i]
	mod1 <- subset(EPImod, subset = EPIClass == C1)
	mod1 <- subset(mod1, subset = Module.1 > 0)
	a <- mean(subset(mod1, subset = Smoking == "Never")@meta.data$Module.1)
	b <- mean(subset(mod1, subset = Smoking == "Smoker")@meta.data$Module.1)

	Res <- data.frame(Cell = C1, Never = a, Smoker = b)
	ResDF <- rbind(ResDF, Res)
}
Cell <- unique(STRmod@meta.data$STRClass)
for(i in 1:length(Cell)){
	C1 <- Cell[i]
	mod1 <- subset(STRmod, subset = STRClass == C1)
	mod1 <- subset(mod1, subset = Module.1 > 0)
	a <- mean(subset(mod1, subset = Smoking == "Never")@meta.data$Module.1)
	b <- mean(subset(mod1, subset = Smoking == "Smoker")@meta.data$Module.1)

	Res <- data.frame(Cell = C1, Never = a, Smoker = b)
	ResDF <- rbind(ResDF, Res)
}
Cell <- unique(ENDmod@meta.data$ENDClass)
for(i in 1:length(Cell)){
	C1 <- Cell[i]
	mod1 <- subset(ENDmod, subset = ENDClass == C1)
	mod1 <- subset(mod1, subset = Module.1 > 0)
	a <- mean(subset(mod1, subset = Smoking == "Never")@meta.data$Module.1)
	b <- mean(subset(mod1, subset = Smoking == "Smoker")@meta.data$Module.1)

	Res <- data.frame(Cell = C1, Never = a, Smoker = b)
	ResDF <- rbind(ResDF, Res)
}
Cell <- unique(LYMmod@meta.data$LYMClass)
for(i in 1:length(Cell)){
	C1 <- Cell[i]
	mod1 <- subset(LYMmod, subset = LYMClass == C1)
	mod1 <- subset(mod1, subset = Module.1 > 0)
	a <- mean(subset(mod1, subset = Smoking == "Never")@meta.data$Module.1)
	b <- mean(subset(mod1, subset = Smoking == "Smoker")@meta.data$Module.1)

	Res <- data.frame(Cell = C1, Never = a, Smoker = b)
	ResDF <- rbind(ResDF, Res)
}
Cell <- unique(MYEmod@meta.data$MYEClass)
for(i in 1:length(Cell)){
	C1 <- Cell[i]
	mod1 <- subset(MYEmod, subset = MYEClass == C1)
	mod1 <- subset(mod1, subset = Module.1 > 0)
	a <- mean(subset(mod1, subset = Smoking == "Never")@meta.data$Module.1)
	b <- mean(subset(mod1, subset = Smoking == "Smoker")@meta.data$Module.1)

	Res <- data.frame(Cell = C1, Never = a, Smoker = b)
	ResDF <- rbind(ResDF, Res)
}

write.table(ResDF, file = "Mod-MeanMatrix-HOUSTIS_ROS2.csv", sep = ",")


list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/REACTOME_AUTOPHAGY.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")
ResDF <- list()

Cell <- unique(EPImod@meta.data$EPIClass)
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

write.table(ResDF, file = "Mod-MeanMatrix-REACTOME_AUTOPHAGY.csv", sep = ",")


list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/REACTOME_CELLULAR_SENESCENCE.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")
ResDF <- list()

Cell <- unique(EPImod@meta.data$EPIClass)
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

write.table(ResDF, file = "Mod-MeanMatrix-REACTOME_CELLULAR_SENESCENCE.csv", sep = ",")


list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/REACTOME_CIRCADIAN_CLOCK.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")
ResDF <- list()

Cell <- unique(EPImod@meta.data$EPIClass)
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

write.table(ResDF, file = "Mod-MeanMatrix-REACTOME_CIRCADIAN_CLOCK.csv", sep = ",")


list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/REACTOME_INTERFERON_SIGNALING.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")
ResDF <- list()

Cell <- unique(EPImod@meta.data$EPIClass)
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

write.table(ResDF, file = "Mod-MeanMatrix-REACTOME_INTERFERON_SIGNALING.csv", sep = ",")


list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/REACTOME_MITOPHAGY.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")
ResDF <- list()

Cell <- unique(EPImod@meta.data$EPIClass)
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

write.table(ResDF, file = "Mod-MeanMatrix-REACTOME_MITOPHAGY.csv", sep = ",")


list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/REACTOME_PYROPTOSIS.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")
ResDF <- list()

Cell <- unique(EPImod@meta.data$EPIClass)
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

write.table(ResDF, file = "Mod-MeanMatrix-REACTOME_PYROPTOSIS.csv", sep = ",")


list2 <- read.delim2("~/Analysis/Tobacco/meta/Modulelist/WP_FERROPTOSIS.txt")
list2 <- list2[-1,]
list2 = unlist(list2)
allgene = rownames(EPI)
list2 = list(allgene[allgene %in% unlist(list2)])

EPImod <- AddModuleScore(EPI, features = list2, pool = allgene, name = "Module.")
STRmod <- AddModuleScore(STR, features = list2, pool = allgene, name = "Module.")
ENDmod <- AddModuleScore(END, features = list2, pool = allgene, name = "Module.")
LYMmod <- AddModuleScore(LYM, features = list2, pool = allgene, name = "Module.")
MYEmod <- AddModuleScore(MYE, features = list2, pool = allgene, name = "Module.")
ResDF <- list()

Cell <- unique(EPImod@meta.data$EPIClass)
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

write.table(ResDF, file = "Mod-MeanMatrix-WP_FERROPTOSIS.csv", sep = ",")
