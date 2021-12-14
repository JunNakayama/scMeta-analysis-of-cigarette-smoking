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

##############
# Extraction with Age
metaage <- data.frame(ident = D@meta.data$orig.ident, Age = D@meta.data$Age)
metaage <- unique(metaage)
metaage <- na.omit(metaage)
rownames(metaage) <- metaage[,1]

ENDNev <- subset(END, subset = Smoking == "Never")
ENDSmo <- subset(END, subset = Smoking == "Smoker")
ENDNev <- subset(ENDNev, subset = Age > 0)
ENDSmo <- subset(ENDSmo, subset = Age > 0)

n <- unique(END@meta.data$ENDClass)
G <- END[rowSums(END) > 0]
G <- rownames(G)
ENDNev <- ENDNev[unlist(G),]
ENDSmo <- ENDSmo[unlist(G),]
ResDF <- list()

for(i in 1:length(n)){

	Cell <- n[i]

	# Never
	Nev <- subset(ENDNev, subset = ENDClass == Cell)
	Smo <- subset(ENDSmo, subset = ENDClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev)
	Nevgene <- data.frame(t(Nevgene[["RNA"]]))

	hoge <- apply(Nevgene, MARGIN = 2, function(y)
			{
			res <- sum(y > 0)
			return(res)
			}
		)

	hoge <- hoge[hoge > 10]
	hoge <- names(hoge)
	Nevgene <- Nevgene[,colnames(Nevgene) %in% unlist(hoge)]


	# Smoker
	Idents(Smo) <- Smo@meta.data$orig.ident
	Smogene <- AverageExpression(Smo)
	Smogene <- data.frame(t(Smogene[["RNA"]]))

	hoge <- apply(Smogene, MARGIN = 2, function(y)
			{
			res <- sum(y > 0)
			return(res)
			}
		)

	hoge <- hoge[hoge > 10]
	hoge <- names(hoge)
	Smogene <- Smogene[,colnames(Smogene) %in% unlist(hoge)]


	# Common genes in Never and Smoker
	NGen <- colnames(Nevgene)
	SGen <- colnames(Smogene)
	G2 <- intersect(NGen, SGen)


	# Merge with Age
	Nevgene <- cbind(Nevgene, metaage[rownames(Nevgene),])
	Smogene <- cbind(Smogene, metaage[rownames(Smogene),])


	## Regresstion analysis
	for(k in 1:length(G2)){
		gene <- G2[k]
	
		LM <- data.frame(Exp = Nevgene[gene], Age = Nevgene$Age)

			m <-lm(LM[,1] ~ LM$Age)
			a <- coef(m)[2]
			p <- summary(m)$coefficients[2,4] 

		LMs <- data.frame(Exp = Smogene[gene], Age = Smogene$Age)

			mm <-lm(LMs[,1] ~ LMs$Age)
			b <- coef(mm)[2]
			p2 <- summary(mm)$coefficients[2,4] 

		Res <- data.frame(gene = gene, Cell = Cell, Never = paste(a), NeverP = paste(p), Smoking = paste(b), SmokingP = paste(p2))
		ResDF <- rbind(ResDF, Res)


	}

}

write.table(ResDF, file = "Age-END-regression.csv", sep = ",")




### Regression analysis with genes
LM <- data.frame(Exp = Nevgene["MALAT1"], Age = Nevgene$Age)
m <-lm(LM[,1] ~ LM$Age)
plot(LM[,1] ~ LM$Age)

LMs <- data.frame(Exp = Smogene["MALAT1"], Age = Smogene$Age)
mm <-lm(LMs[,1] ~ LMs$Age)
plot(LMs[,1] ~ LMs$Age)


#### row data extraction
ENDNev <- subset(END, subset = Smoking == "Never")
ENDSmo <- subset(END, subset = Smoking == "Smoker")
ENDNev <- subset(ENDNev, subset = Age > 0)
ENDSmo <- subset(ENDSmo, subset = Age > 0)

gene = "MALAT1"
n <- unique(END@meta.data$ENDClass)
Res <- list()
for(i in 1:length(n)){
	Cell <- n[i]
	Nev <- subset(ENDNev, subset = ENDClass == Cell)
	Smo <- subset(ENDSmo, subset = ENDClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev, features = gene)
	Idents(Smo) <- Smo@meta.data$orig.ident
	Smogene <- AverageExpression(Smo, features = gene)

	NRes <- data.frame(Smoking = "Never", Cell = Cell, Nevgene)
	SRes <- data.frame(Smoking = "Smoker", Cell = Cell, Smogene)

	Res <- bind_rows(Res, NRes, SRes)
}
write.table(Res, file = "Age-END-MALAT1-regression.csv", sep = ",")


gene = "FTL"
n <- unique(END@meta.data$ENDClass)
Res <- list()
for(i in 1:length(n)){
	Cell <- n[i]
	Nev <- subset(ENDNev, subset = ENDClass == Cell)
	Smo <- subset(ENDSmo, subset = ENDClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev, features = gene)
	Idents(Smo) <- Smo@meta.data$orig.ident
	Smogene <- AverageExpression(Smo, features = gene)

	NRes <- data.frame(Smoking = "Never", Cell = Cell, Nevgene)
	SRes <- data.frame(Smoking = "Smoker", Cell = Cell, Smogene)

	Res <- bind_rows(Res, NRes, SRes)
}
write.table(Res, file = "Age-END-FTL-regression.csv", sep = ",")




##############
##############
### Hetamap Visualaiztion
metaage <- data.frame(ident = D@meta.data$orig.ident, Age = D@meta.data$Age)
metaage <- unique(metaage)
metaage <- na.omit(metaage)
rownames(metaage) <- metaage[,1]

G <- read_csv("genelist2.csv")
ENDNev <- subset(END, subset = Smoking == "Never")
ENDSmo <- subset(END, subset = Smoking == "Smoker")
ENDNev <- subset(ENDNev, subset = Age > 0)
ENDSmo <- subset(ENDSmo, subset = Age > 0)

n <- unique(END@meta.data$ENDClass)
ENDNev <- ENDNev[unlist(G),]
ENDSmo <- ENDSmo[unlist(G),]
ResDF <- list()

for(i in 1:length(n)){

	Cell <- n[i]

	# Nev
	Nev <- subset(ENDNev, subset = ENDClass == Cell)
	Smo <- subset(ENDSmo, subset = ENDClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev)
	Nevgene <- data.frame(t(Nevgene[["RNA"]]))

	hoge <- apply(Nevgene, MARGIN = 2, function(y)
			{
			res <- sum(y > 0)
			return(res)
			}
		)

	hoge <- hoge[hoge > 10]
	hoge <- names(hoge)
	Nevgene <- Nevgene[,colnames(Nevgene) %in% unlist(hoge)]


	# Smo
	Idents(Smo) <- Smo@meta.data$orig.ident
	Smogene <- AverageExpression(Smo)
	Smogene <- data.frame(t(Smogene[["RNA"]]))

	hoge <- apply(Smogene, MARGIN = 2, function(y)
			{
			res <- sum(y > 0)
			return(res)
			}
		)

	hoge <- hoge[hoge > 10]
	hoge <- names(hoge)
	Smogene <- Smogene[,colnames(Smogene) %in% unlist(hoge)]

	NGen <- colnames(Nevgene)
	SGen <- colnames(Smogene)
	G2 <- intersect(NGen, SGen)

	Nevgene <- cbind(Nevgene, metaage[rownames(Nevgene),])
	Smogene <- cbind(Smogene, metaage[rownames(Smogene),])


	if(length(Nevgene$ident) < 10){
		paste(Cell)
	}else{
		if(length(Smogene$ident) < 10){
			paste(Cell)
		}else{
			for(k in 1:length(G2)){
				gene <- G2[k]
			
				LM <- data.frame(Exp = Nevgene[gene], Age = Nevgene$Age)

				m <-lm(LM[,1] ~ LM$Age)
				a <- coef(m)[2]

				LMs <- data.frame(Exp = Smogene[gene], Age = Smogene$Age)

				mm <-lm(LMs[,1] ~ LMs$Age)
				b <- coef(mm)[2]

				Res <- data.frame(gene = gene, Cell = Cell, Never = paste(a), Smoking = paste(b))
				ResDF <- rbind(ResDF, Res)

			}
		}

	}

}
STRNev <- subset(STR, subset = Smoking == "Never")
STRSmo <- subset(STR, subset = Smoking == "Smoker")
STRNev <- subset(STRNev, subset = Age > 0)
STRSmo <- subset(STRSmo, subset = Age > 0)
n <- unique(STR@meta.data$STRClass)
STRNev <- STRNev[unlist(G),]
STRSmo <- STRSmo[unlist(G),]
for(i in 1:length(n)){

	Cell <- n[i]

	# Nev
	Nev <- subset(STRNev, subset = STRClass == Cell)
	Smo <- subset(STRSmo, subset = STRClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev)
	Nevgene <- data.frame(t(Nevgene[["RNA"]]))

	hoge <- apply(Nevgene, MARGIN = 2, function(y)
			{
			res <- sum(y > 0)
			return(res)
			}
		)

	hoge <- hoge[hoge > 10]
	hoge <- names(hoge)
	Nevgene <- Nevgene[,colnames(Nevgene) %in% unlist(hoge)]


	# Smo
	Idents(Smo) <- Smo@meta.data$orig.ident
	Smogene <- AverageExpression(Smo)
	Smogene <- data.frame(t(Smogene[["RNA"]]))

	hoge <- apply(Smogene, MARGIN = 2, function(y)
			{
			res <- sum(y > 0)
			return(res)
			}
		)

	hoge <- hoge[hoge > 10]
	hoge <- names(hoge)
	Smogene <- Smogene[,colnames(Smogene) %in% unlist(hoge)]

	NGen <- colnames(Nevgene)
	SGen <- colnames(Smogene)
	G2 <- intersect(NGen, SGen)

	Nevgene <- cbind(Nevgene, metaage[rownames(Nevgene),])
	Smogene <- cbind(Smogene, metaage[rownames(Smogene),])


	if(length(Nevgene$ident) < 10){
		paste(Cell)
	}else{
		if(length(Smogene$ident) < 10){
			paste(Cell)
		}else{
			for(k in 1:length(G2)){
				gene <- G2[k]
			
				LM <- data.frame(Exp = Nevgene[gene], Age = Nevgene$Age)

				m <-lm(LM[,1] ~ LM$Age)
				a <- coef(m)[2]

				LMs <- data.frame(Exp = Smogene[gene], Age = Smogene$Age)

				mm <-lm(LMs[,1] ~ LMs$Age)
				b <- coef(mm)[2]

				Res <- data.frame(gene = gene, Cell = Cell, Never = paste(a), Smoking = paste(b))
				ResDF <- rbind(ResDF, Res)

			}
		}

	}

}
EPINev <- subset(EPI, subset = Smoking == "Never")
EPISmo <- subset(EPI, subset = Smoking == "Smoker")
EPINev <- subset(EPINev, subset = Age > 0)
EPISmo <- subset(EPISmo, subset = Age > 0)
n <- unique(EPI@meta.data$EPIClass)
EPINev <- EPINev[unlist(G),]
EPISmo <- EPISmo[unlist(G),]
for(i in 1:length(n)){

	Cell <- n[i]

	# Nev
	Nev <- subset(EPINev, subset = EPIClass == Cell)
	Smo <- subset(EPISmo, subset = EPIClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev)
	Nevgene <- data.frame(t(Nevgene[["RNA"]]))

	hoge <- apply(Nevgene, MARGIN = 2, function(y)
			{
			res <- sum(y > 0)
			return(res)
			}
		)

	hoge <- hoge[hoge > 10]
	hoge <- names(hoge)
	Nevgene <- Nevgene[,colnames(Nevgene) %in% unlist(hoge)]


	# Smo
	Idents(Smo) <- Smo@meta.data$orig.ident
	Smogene <- AverageExpression(Smo)
	Smogene <- data.frame(t(Smogene[["RNA"]]))

	hoge <- apply(Smogene, MARGIN = 2, function(y)
			{
			res <- sum(y > 0)
			return(res)
			}
		)

	hoge <- hoge[hoge > 10]
	hoge <- names(hoge)
	Smogene <- Smogene[,colnames(Smogene) %in% unlist(hoge)]

	NGen <- colnames(Nevgene)
	SGen <- colnames(Smogene)
	G2 <- intersect(NGen, SGen)

	Nevgene <- cbind(Nevgene, metaage[rownames(Nevgene),])
	Smogene <- cbind(Smogene, metaage[rownames(Smogene),])


	if(length(Nevgene$ident) < 10){
		paste(Cell)
	}else{
		if(length(Smogene$ident) < 10){
			paste(Cell)
		}else{
			for(k in 1:length(G2)){
				gene <- G2[k]
			
				LM <- data.frame(Exp = Nevgene[gene], Age = Nevgene$Age)

				m <-lm(LM[,1] ~ LM$Age)
				a <- coef(m)[2]

				LMs <- data.frame(Exp = Smogene[gene], Age = Smogene$Age)

				mm <-lm(LMs[,1] ~ LMs$Age)
				b <- coef(mm)[2]

				Res <- data.frame(gene = gene, Cell = Cell, Never = paste(a), Smoking = paste(b))
				ResDF <- rbind(ResDF, Res)

			}
		}

	}

}
LYMNev <- subset(LYM, subset = Smoking == "Never")
LYMSmo <- subset(LYM, subset = Smoking == "Smoker")
LYMNev <- subset(LYMNev, subset = Age > 0)
LYMSmo <- subset(LYMSmo, subset = Age > 0)
n <- unique(LYM@meta.data$LYMClass)
LYMNev <- LYMNev[unlist(G),]
LYMSmo <- LYMSmo[unlist(G),]
for(i in 1:length(n)){

	Cell <- n[i]

	# Nev
	Nev <- subset(LYMNev, subset = LYMClass == Cell)
	Smo <- subset(LYMSmo, subset = LYMClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev)
	Nevgene <- data.frame(t(Nevgene[["RNA"]]))

	hoge <- apply(Nevgene, MARGIN = 2, function(y)
			{
			res <- sum(y > 0)
			return(res)
			}
		)

	hoge <- hoge[hoge > 10]
	hoge <- names(hoge)
	Nevgene <- Nevgene[,colnames(Nevgene) %in% unlist(hoge)]


	# Smo
	Idents(Smo) <- Smo@meta.data$orig.ident
	Smogene <- AverageExpression(Smo)
	Smogene <- data.frame(t(Smogene[["RNA"]]))

	hoge <- apply(Smogene, MARGIN = 2, function(y)
			{
			res <- sum(y > 0)
			return(res)
			}
		)

	hoge <- hoge[hoge > 10]
	hoge <- names(hoge)
	Smogene <- Smogene[,colnames(Smogene) %in% unlist(hoge)]

	NGen <- colnames(Nevgene)
	SGen <- colnames(Smogene)
	G2 <- intersect(NGen, SGen)

	Nevgene <- cbind(Nevgene, metaage[rownames(Nevgene),])
	Smogene <- cbind(Smogene, metaage[rownames(Smogene),])


	if(length(Nevgene$ident) < 10){
		paste(Cell)
	}else{
		if(length(Smogene$ident) < 10){
			paste(Cell)
		}else{
			for(k in 1:length(G2)){
				gene <- G2[k]
			
				LM <- data.frame(Exp = Nevgene[gene], Age = Nevgene$Age)

				m <-lm(LM[,1] ~ LM$Age)
				a <- coef(m)[2]

				LMs <- data.frame(Exp = Smogene[gene], Age = Smogene$Age)

				mm <-lm(LMs[,1] ~ LMs$Age)
				b <- coef(mm)[2]

				Res <- data.frame(gene = gene, Cell = Cell, Never = paste(a), Smoking = paste(b))
				ResDF <- rbind(ResDF, Res)

			}
		}

	}

}
MYENev <- subset(MYE, subset = Smoking == "Never")
MYESmo <- subset(MYE, subset = Smoking == "Smoker")
MYENev <- subset(MYENev, subset = Age > 0)
MYESmo <- subset(MYESmo, subset = Age > 0)
n <- unique(MYE@meta.data$MYEClass)
MYENev <- MYENev[unlist(G),]
MYESmo <- MYESmo[unlist(G),]
for(i in 1:length(n)){

	Cell <- n[i]

	# Nev
	Nev <- subset(MYENev, subset = MYEClass == Cell)
	Smo <- subset(MYESmo, subset = MYEClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev)
	Nevgene <- data.frame(t(Nevgene[["RNA"]]))

	hoge <- apply(Nevgene, MARGIN = 2, function(y)
			{
			res <- sum(y > 0)
			return(res)
			}
		)

	hoge <- hoge[hoge > 10]
	hoge <- names(hoge)
	Nevgene <- Nevgene[,colnames(Nevgene) %in% unlist(hoge)]


	# Smo
	Idents(Smo) <- Smo@meta.data$orig.ident
	Smogene <- AverageExpression(Smo)
	Smogene <- data.frame(t(Smogene[["RNA"]]))

	hoge <- apply(Smogene, MARGIN = 2, function(y)
			{
			res <- sum(y > 0)
			return(res)
			}
		)

	hoge <- hoge[hoge > 10]
	hoge <- names(hoge)
	Smogene <- Smogene[,colnames(Smogene) %in% unlist(hoge)]

	NGen <- colnames(Nevgene)
	SGen <- colnames(Smogene)
	G2 <- intersect(NGen, SGen)

	Nevgene <- cbind(Nevgene, metaage[rownames(Nevgene),])
	Smogene <- cbind(Smogene, metaage[rownames(Smogene),])



	if(length(Nevgene$ident) < 10){
		paste(Cell)
	}else{
		if(length(Smogene$ident) < 10){
			paste(Cell)
		}else{
			for(k in 1:length(G2)){
				gene <- G2[k]
			
				LM <- data.frame(Exp = Nevgene[gene], Age = Nevgene$Age)

				m <-lm(LM[,1] ~ LM$Age)
				a <- coef(m)[2]

				LMs <- data.frame(Exp = Smogene[gene], Age = Smogene$Age)

				mm <-lm(LMs[,1] ~ LMs$Age)
				b <- coef(mm)[2]

				Res <- data.frame(gene = gene, Cell = Cell, Never = paste(a), Smoking = paste(b))
				ResDF <- rbind(ResDF, Res)

			}
		}

	}

}

write.table(ResDF, file = "Age-AllCell-genelist-regression.csv", sep = ",")








