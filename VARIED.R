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

#######################
## VARIED
#######################
#######################
## Normalized Closseness Centrarities
ResDF <- list()
nomcentres <- list()
snomcentres <- list()
cl = unique(EPI@meta.data$EPIClass)
for(i in 1:length(cl)){
	Cell <- cl[i]
	netcell = subset(EPI, subset = EPIClass == Cell)
	Nev <- subset(netcell, subset = Smoking == "Never")
	Smo <- subset(netcell, subset = Smoking == "Smoker")

	cNet <- cor(as.matrix(GetAssayData(Nev)))
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)

	cNet <- cor(as.matrix(GetAssayData(Smo)))
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	snomcentres[i] = list((n-1)/sta)

}
names(nomcentres) = cl
names(snomcentres) = cl
nomcentres <- data.frame(reshape2::melt(nomcentres), Smoking = "Never")
snomcentres <- data.frame(reshape2::melt(snomcentres), Smoking = "Smoker")
ResDF = rbind(ResDF, nomcentres)
ResDF = rbind(ResDF, snomcentres)

cl = unique(STR@meta.data$STRClass)
nomcentres = list()
snomcentres <- list()
for(i in 1:length(cl)){
	Cell <- cl[i]
	netcell = subset(STR, subset = STRClass == Cell)
	Nev <- subset(netcell, subset = Smoking == "Never")
	Smo <- subset(netcell, subset = Smoking == "Smoker")

	cNet <- cor(as.matrix(GetAssayData(Nev)))
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)

	cNet <- cor(as.matrix(GetAssayData(Smo)))
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	snomcentres[i] = list((n-1)/sta)

}
names(nomcentres) = cl
names(snomcentres) = cl
nomcentres <- data.frame(reshape2::melt(nomcentres), Smoking = "Never")
snomcentres <- data.frame(reshape2::melt(snomcentres), Smoking = "Smoker")
ResDF = rbind(ResDF, nomcentres)
ResDF = rbind(ResDF, snomcentres)

cl = unique(END@meta.data$ENDClass)
nomcentres = list()
snomcentres <- list()
for(i in 1:length(cl)){
	Cell <- cl[i]
	netcell = subset(END, subset = ENDClass == Cell)
	Nev <- subset(netcell, subset = Smoking == "Never")
	Smo <- subset(netcell, subset = Smoking == "Smoker")

	cNet <- cor(as.matrix(GetAssayData(Nev)))
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)

	cNet <- cor(as.matrix(GetAssayData(Smo)))
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	snomcentres[i] = list((n-1)/sta)

}
names(nomcentres) = cl
names(snomcentres) = cl
nomcentres <- data.frame(reshape2::melt(nomcentres), Smoking = "Never")
snomcentres <- data.frame(reshape2::melt(snomcentres), Smoking = "Smoker")
ResDF = rbind(ResDF, nomcentres)
ResDF = rbind(ResDF, snomcentres)

cl = unique(LYM@meta.data$LYMClass)
nomcentres = list()
snomcentres <- list()
for(i in 1:length(cl)){
	Cell <- cl[i]
	netcell = subset(LYM, subset = LYMClass == Cell)
	Nev <- subset(netcell, subset = Smoking == "Never")
	Smo <- subset(netcell, subset = Smoking == "Smoker")

	cNet <- cor(as.matrix(GetAssayData(Nev)))
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)

	cNet <- cor(as.matrix(GetAssayData(Smo)))
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	snomcentres[i] = list((n-1)/sta)

}
names(nomcentres) = cl
names(snomcentres) = cl
nomcentres <- data.frame(reshape2::melt(nomcentres), Smoking = "Never")
snomcentres <- data.frame(reshape2::melt(snomcentres), Smoking = "Smoker")
ResDF = rbind(ResDF, nomcentres)
ResDF = rbind(ResDF, snomcentres)

cl = unique(MYE@meta.data$MYEClass)
nomcentres = list()
snomcentres <- list()
for(i in 1:length(cl)){
	Cell <- cl[i]
	netcell = subset(MYE, subset = MYEClass == Cell)
	Nev <- subset(netcell, subset = Smoking == "Never")
	Smo <- subset(netcell, subset = Smoking == "Smoker")

	cNet <- cor(as.matrix(GetAssayData(Nev)))
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)

	cNet <- cor(as.matrix(GetAssayData(Smo)))
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	snomcentres[i] = list((n-1)/sta)

}
names(nomcentres) = cl
names(snomcentres) = cl
nomcentres <- data.frame(reshape2::melt(nomcentres), Smoking = "Never")
snomcentres <- data.frame(reshape2::melt(snomcentres), Smoking = "Smoker")
ResDF = rbind(ResDF, nomcentres)
ResDF = rbind(ResDF, snomcentres)

write.table(ResDF, file = "Cent-allcell.csv", sep = ',')


#######################
#######################
## Visualization
g <- ggplot(ResDF)
g + geom_boxplot(aes(x = L1, y = value, fill = Smoking), alpha = 0.75) +
	scale_fill_manual(values = c("skyblue2","red3")) +
	theme_classic() + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1))

g <- ggplot(subset(ResDF, subset = L1 %in% c("Basal", "Basal-px", "Basophil", "Serous", "TracheaBasal")))
g + geom_boxplot(aes(x = L1, y = value, fill = Smoking), alpha = 0.75) +
	scale_fill_manual(values = c("skyblue2","red3")) +
	theme_classic() + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1))


library(ggridges)
ridge_cent <- read_csv("ridge_cent.csv") ## Extract the specific cell type
data <- reshape2::melt(ridge_cent)

ggplot(data, aes(x = value, y = variable, fill = variable)) +
	geom_density_ridges(jittered_points = TRUE, position = position_points_jitter(width = 0.05, height = 0), 
								alpha = 0.75, point_size = 3,  point_shape = '|', point_alpha = 0.75) +
	scale_fill_cyclical(values = c("red3", "skyblue2")) +
	xlim(1,4) +
	theme_classic()
	
