library(Seurat)
library(harmony)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(MAST)
library(DoubletFinder)



#######################
#######################
####### VARIED ########
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

	cNet <- cor(as.matrix(GetAssayData(Nev, assay = "SCT")))
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)

	cNet <- cor(as.matrix(GetAssayData(Smo, assay = "SCT")))
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

	cNet <- cor(as.matrix(GetAssayData(Nev, assay = "SCT")))
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)

	cNet <- cor(as.matrix(GetAssayData(Smo, assay = "SCT")))
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

	cNet <- cor(as.matrix(GetAssayData(Nev, assay = "SCT")))
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)

	cNet <- cor(as.matrix(GetAssayData(Smo, assay = "SCT")))
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

	cNet <- cor(as.matrix(GetAssayData(Nev, assay = "SCT")))
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)

	cNet <- cor(as.matrix(GetAssayData(Smo, assay = "SCT")))
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

	cNet <- cor(as.matrix(GetAssayData(Nev, assay = "SCT")))
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)

	cNet <- cor(as.matrix(GetAssayData(Smo, assay = "SCT")))
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
#######################
#######################

library(ggridges)

ridge_cent <- read.csv("ridge_cent.csv")

data <- reshape2::melt(ridge_cent)


png("VARIED-Basald.png", width = 1000, height = 1000)
ggplot(data, aes(x = value, y = variable, fill = variable)) +
	geom_density_ridges(jittered_points = TRUE, position = position_points_jitter(width = 0.05, height = 0), 
								alpha = 0.75, point_size = 3,  point_shape = '|', point_alpha = 0.75) +
	scale_fill_cyclical(values = c("red3", "skyblue2")) +
	xlim(1,4) +
	theme_classic()
dev.off()


#######################
#######################
#######################
#######################


## VARIED Random sampling test

Basald <- subset(EPI, subset = EPIClass == c("Basal-d"))
mat <- as.matrix(GetAssayData(Basald, assay = "SCT"))

medmerge <- list()

# 10
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 10)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- med[,1]
write.table(ResDF, file = "Cent-random-Basald10.csv", sep = ',')

# 50
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 50)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- cbind(medmerge, med[,1])
write.table(ResDF, file = "Cent-random-Basald50.csv", sep = ',')

# 100
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 100)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- cbind(medmerge, med[,1])
write.table(ResDF, file = "Cent-random-Basald100.csv", sep = ',')

# 500
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 100)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- cbind(medmerge, med[,1])
write.table(ResDF, file = "Cent-random-Basald500.csv", sep = ',')

# 1000
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 100)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- cbind(medmerge, med[,1])
write.table(ResDF, file = "Cent-random-Basald1000.csv", sep = ',')


colnames(medmerge) <- paste(c(10, 50, 100, 500, 1000))
write.table(medmerge, file = "Cent-random-Basaldmed-merge.csv", sep = ',')




#######################
#######################
Basal <- subset(EPI, subset = EPIClass == c("Basal"))
mat <- as.matrix(GetAssayData(Basal, assay = "SCT"))

medmerge <- list()

# 10
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 10)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- med[,1]
write.table(ResDF, file = "Cent-random-Basal10.csv", sep = ',')

# 50
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 50)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- cbind(medmerge, med[,1])
write.table(ResDF, file = "Cent-random-Basal50.csv", sep = ',')

# 100
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 100)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- cbind(medmerge, med[,1])
write.table(ResDF, file = "Cent-random-Basal100.csv", sep = ',')

# 500
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 100)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- cbind(medmerge, med[,1])
write.table(ResDF, file = "Cent-random-Basal500.csv", sep = ',')

# 1000
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 100)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- cbind(medmerge, med[,1])
write.table(ResDF, file = "Cent-random-Basal1000.csv", sep = ',')


colnames(medmerge) <- paste(c(10, 50, 100, 500, 1000))
write.table(medmerge, file = "Cent-random-Basalmed-merge.csv", sep = ',')


#######################
#######################
AT1 <- subset(EPI, subset = EPIClass == c("AT1"))
mat <- as.matrix(GetAssayData(AT1, assay = "SCT"))

medmerge <- list()

# 10
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 10)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- med[,1]
write.table(ResDF, file = "Cent-random-AT1-10.csv", sep = ',')

# 50
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 50)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- cbind(medmerge, med[,1])
write.table(ResDF, file = "Cent-random-AT1-50.csv", sep = ',')

# 100
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 100)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- cbind(medmerge, med[,1])
write.table(ResDF, file = "Cent-random-AT1-100.csv", sep = ',')

# 500
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 100)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- cbind(medmerge, med[,1])
write.table(ResDF, file = "Cent-random-AT1-500.csv", sep = ',')

# 1000
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 100)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- cbind(medmerge, med[,1])
write.table(ResDF, file = "Cent-random-AT1-1000.csv", sep = ',')


colnames(medmerge) <- paste(c(10, 50, 100, 500, 1000))
write.table(medmerge, file = "Cent-random-AT1med-merge.csv", sep = ',')


#######################
#######################
AT2 <- subset(EPI, subset = EPIClass == c("AT2"))
mat <- as.matrix(GetAssayData(AT2, assay = "SCT"))

medmerge <- list()

# 10
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 10)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- med[,1]
write.table(ResDF, file = "Cent-random-AT2-10.csv", sep = ',')

# 50
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 50)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- cbind(medmerge, med[,1])
write.table(ResDF, file = "Cent-random-AT2-50.csv", sep = ',')

# 100
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 100)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- cbind(medmerge, med[,1])
write.table(ResDF, file = "Cent-random-AT2-100.csv", sep = ',')

# 500
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 100)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- cbind(medmerge, med[,1])
write.table(ResDF, file = "Cent-random-AT2-500.csv", sep = ',')

# 1000
ResDF <- list()
nomcentres <- list()
med <- list()
for(i in 1:100){
	test <- mat[,sample(1:ncol(mat), 100)]
	cNet <- cor(test)
	cNet <- 1 - cNet
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)
	nomcentres[i] = list((n-1)/sta)
	med[i] <- median(unlist(nomcentres[i]))
}
ResDF <- reshape2::melt(nomcentres)
med <- reshape2::melt(med)
medmerge <- cbind(medmerge, med[,1])
write.table(ResDF, file = "Cent-random-AT2-1000.csv", sep = ',')


colnames(medmerge) <- paste(c(10, 50, 100, 500, 1000))
write.table(medmerge, file = "Cent-random-AT2med-merge.csv", sep = ',')



	
