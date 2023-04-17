library(Seurat)
library(harmony)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(MAST)
library(DoubletFinder)


#### AGED
# Extraction of agiging data
# END
metaage <- data.frame(ident = sDrev@meta.data$orig.ident, Age = sDrev@meta.data$Age)
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
	Nevgene <- data.frame(t(Nevgene[["SCT"]]))

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
	Smogene <- data.frame(t(Smogene[["SCT"]]))

	hoge <- apply(Smogene, MARGIN = 2, function(y)
			{
			res <- sum(y > 0)
			return(res)
			}
		)

	hoge <- hoge[hoge > 10]
	hoge <- names(hoge)
	Smogene <- Smogene[,colnames(Smogene) %in% unlist(hoge)]


	# Commonn genes in Never-smoker and Smoker
	NGen <- colnames(Nevgene)
	SGen <- colnames(Smogene)
	G2 <- intersect(NGen, SGen)


	# Integration with Age
	Nevgene <- cbind(Nevgene, metaage[rownames(Nevgene),])
	Smogene <- cbind(Smogene, metaage[rownames(Smogene),])


	## Extraction and Cox analysis of common genes which has more 10 specimens
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


##################################
##################################
# STR

STRNev <- subset(STR, subset = Smoking == "Never")
STRSmo <- subset(STR, subset = Smoking == "Smoker")

STRNev <- subset(STRNev, subset = Age > 0)
STRSmo <- subset(STRSmo, subset = Age > 0)


n <- unique(STR@meta.data$STRClass)
G <- STR[rowSums(STR) > 0]
G <- rownames(G)
STRNev <- STRNev[unlist(G),]
STRSmo <- STRSmo[unlist(G),]

ResDF <- list()

for(i in 1:length(n)){

	Cell <- n[i]

	# Never
	Nev <- subset(STRNev, subset = STRClass == Cell)
	Smo <- subset(STRSmo, subset = STRClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev)
	Nevgene <- data.frame(t(Nevgene[["SCT"]]))

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
	Smogene <- data.frame(t(Smogene[["SCT"]]))

	hoge <- apply(Smogene, MARGIN = 2, function(y)
			{
			res <- sum(y > 0)
			return(res)
			}
		)

	hoge <- hoge[hoge > 10]
	hoge <- names(hoge)
	Smogene <- Smogene[,colnames(Smogene) %in% unlist(hoge)]


	# Commonn genes in Never-smoker and Smoker
	NGen <- colnames(Nevgene)
	SGen <- colnames(Smogene)
	G2 <- intersect(NGen, SGen)


	# Integration with Age
	Nevgene <- cbind(Nevgene, metaage[rownames(Nevgene),])
	Smogene <- cbind(Smogene, metaage[rownames(Smogene),])


	## Extraction and Cox analysis of common genes which has more 10 specimens
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

write.table(ResDF, file = "Age-STR-regression.csv", sep = ",")


##################################
##################################
# EPI

EPINev <- subset(EPI, subset = Smoking == "Never")
EPISmo <- subset(EPI, subset = Smoking == "Smoker")

EPINev <- subset(EPINev, subset = Age > 0)
EPISmo <- subset(EPISmo, subset = Age > 0)


n <- unique(EPI@meta.data$EPIClass)
G <- EPI[rowSums(EPI) > 0]
G <- rownames(G)
EPINev <- EPINev[unlist(G),]
EPISmo <- EPISmo[unlist(G),]

ResDF <- list()

for(i in 1:length(n)){

	Cell <- n[i]

	# Never
	Nev <- subset(EPINev, subset = EPIClass == Cell)
	Smo <- subset(EPISmo, subset = EPIClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev)
	Nevgene <- data.frame(t(Nevgene[["SCT"]]))

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
	Smogene <- data.frame(t(Smogene[["SCT"]]))

	hoge <- apply(Smogene, MARGIN = 2, function(y)
			{
			res <- sum(y > 0)
			return(res)
			}
		)

	hoge <- hoge[hoge > 10]
	hoge <- names(hoge)
	Smogene <- Smogene[,colnames(Smogene) %in% unlist(hoge)]


	# Commonn genes in Never-smoker and Smoker
	NGen <- colnames(Nevgene)
	SGen <- colnames(Smogene)
	G2 <- intersect(NGen, SGen)


	# Integration with Age
	Nevgene <- cbind(Nevgene, metaage[rownames(Nevgene),])
	Smogene <- cbind(Smogene, metaage[rownames(Smogene),])


	## Extraction and Cox analysis of common genes which has more 10 specimens
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

write.table(ResDF, file = "Age-EPI-regression.csv", sep = ",")



##################################
##################################
# LYM

LYMNev <- subset(LYM, subset = Smoking == "Never")
LYMSmo <- subset(LYM, subset = Smoking == "Smoker")

LYMNev <- subset(LYMNev, subset = Age > 0)
LYMSmo <- subset(LYMSmo, subset = Age > 0)


n <- unique(LYM@meta.data$LYMClass)
G <- LYM[rowSums(LYM) > 0]
G <- rownames(G)
LYMNev <- LYMNev[unlist(G),]
LYMSmo <- LYMSmo[unlist(G),]

ResDF <- list()

for(i in 1:length(n)){

	Cell <- n[i]

	# Never
	Nev <- subset(LYMNev, subset = LYMClass == Cell)
	Smo <- subset(LYMSmo, subset = LYMClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev)
	Nevgene <- data.frame(t(Nevgene[["SCT"]]))

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
	Smogene <- data.frame(t(Smogene[["SCT"]]))

	hoge <- apply(Smogene, MARGIN = 2, function(y)
			{
			res <- sum(y > 0)
			return(res)
			}
		)

	hoge <- hoge[hoge > 10]
	hoge <- names(hoge)
	Smogene <- Smogene[,colnames(Smogene) %in% unlist(hoge)]


	# Commonn genes in Never-smoker and Smoker
	NGen <- colnames(Nevgene)
	SGen <- colnames(Smogene)
	G2 <- intersect(NGen, SGen)


	# Integration with Age
	Nevgene <- cbind(Nevgene, metaage[rownames(Nevgene),])
	Smogene <- cbind(Smogene, metaage[rownames(Smogene),])


	## Extraction and Cox analysis of common genes which has more 10 specimens
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

write.table(ResDF, file = "Age-LYM-regression.csv", sep = ",")





##################################
##################################
# MYE

metaage <- data.frame(ident = D@meta.data$orig.ident, Age = D@meta.data$Age)
metaage <- unique(metaage)
metaage <- na.omit(metaage)
rownames(metaage) <- metaage[,1]



MYENev <- subset(MYE, subset = Smoking == "Never")
MYESmo <- subset(MYE, subset = Smoking == "Smoker")

MYENev <- subset(MYENev, subset = Age > 0)
MYESmo <- subset(MYESmo, subset = Age > 0)


n <- unique(MYE@meta.data$MYEClass)
G <- MYE[rowSums(MYE) > 0]
G <- rownames(G)
MYENev <- MYENev[unlist(G),]
MYESmo <- MYESmo[unlist(G),]

ResDF <- list()

for(i in 1:length(n)){

	Cell <- n[i]

	# Never
	Nev <- subset(MYENev, subset = MYEClass == Cell)
	Smo <- subset(MYESmo, subset = MYEClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev)
	Nevgene <- data.frame(t(Nevgene[["SCT"]]))

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
	Smogene <- data.frame(t(Smogene[["SCT"]]))

	hoge <- apply(Smogene, MARGIN = 2, function(y)
			{
			res <- sum(y > 0)
			return(res)
			}
		)

	hoge <- hoge[hoge > 10]
	hoge <- names(hoge)
	Smogene <- Smogene[,colnames(Smogene) %in% unlist(hoge)]


	# Commonn genes in Never-smoker and Smoker
	NGen <- colnames(Nevgene)
	SGen <- colnames(Smogene)
	G2 <- intersect(NGen, SGen)


	# Integration with Age
	Nevgene <- cbind(Nevgene, metaage[rownames(Nevgene),])
	Smogene <- cbind(Smogene, metaage[rownames(Smogene),])


	## Extraction and Cox analysis of common genes which has more 10 specimens
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

write.table(ResDF, file = "Age-MYE-regression.csv", sep = ",")






##################################
##################################
##################################
##################################
##################################
##################################
##################################



EPINev <- subset(EPI, subset = Smoking == "Never")
EPISmo <- subset(EPI, subset = Smoking == "Smoker")

gene = "MALAT1"
n <- unique(EPI@meta.data$EPIClass)
Res <- list()
for(i in 1:length(n)){
	Cell <- n[i]
	Nev <- subset(EPINev, subset = EPIClass == Cell)
	Smo <- subset(EPISmo, subset = EPIClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev, features = gene)
	Idents(Smo) <- Smo@meta.data$orig.ident
	Smogene <- AverageExpression(Smo, features = gene)

	NRes <- data.frame(Smoking = "Never", Cell = Cell, Nevgene)
	SRes <- data.frame(Smoking = "Smoker", Cell = Cell, Smogene)

	Res <- bind_rows(Res, NRes, SRes)
}

write.table(Res, file = "Age-EPI-MALAT1-regression.csv", sep = ",")


gene = "SCGB3A1"
n <- unique(EPI@meta.data$EPIClass)
Res <- list()
for(i in 1:length(n)){
	Cell <- n[i]
	Nev <- subset(EPINev, subset = EPIClass == Cell)
	Smo <- subset(EPISmo, subset = EPIClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev, features = gene)
	Idents(Smo) <- Smo@meta.data$orig.ident
	Smogene <- AverageExpression(Smo, features = gene)

	NRes <- data.frame(Smoking = "Never", Cell = Cell, Nevgene)
	SRes <- data.frame(Smoking = "Smoker", Cell = Cell, Smogene)

	Res <- bind_rows(Res, NRes, SRes)
}

write.table(Res, file = "Age-EPI-SCGB3A1-regression.csv", sep = ",")


gene = "SFTPB"
n <- unique(EPI@meta.data$EPIClass)
Res <- list()
for(i in 1:length(n)){
	Cell <- n[i]
	Nev <- subset(EPINev, subset = EPIClass == Cell)
	Smo <- subset(EPISmo, subset = EPIClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev, features = gene)
	Idents(Smo) <- Smo@meta.data$orig.ident
	Smogene <- AverageExpression(Smo, features = gene)

	NRes <- data.frame(Smoking = "Never", Cell = Cell, Nevgene)
	SRes <- data.frame(Smoking = "Smoker", Cell = Cell, Smogene)

	Res <- bind_rows(Res, NRes, SRes)
}

write.table(Res, file = "Age-EPI-SFTPB-regression.csv", sep = ",")


gene = "SFTPC"
n <- unique(EPI@meta.data$EPIClass)
Res <- list()
for(i in 1:length(n)){
	Cell <- n[i]
	Nev <- subset(EPINev, subset = EPIClass == Cell)
	Smo <- subset(EPISmo, subset = EPIClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev, features = gene)
	Idents(Smo) <- Smo@meta.data$orig.ident
	Smogene <- AverageExpression(Smo, features = gene)

	NRes <- data.frame(Smoking = "Never", Cell = Cell, Nevgene)
	SRes <- data.frame(Smoking = "Smoker", Cell = Cell, Smogene)

	Res <- bind_rows(Res, NRes, SRes)
}

write.table(Res, file = "Age-EPI-SFTPC-regression.csv", sep = ",")





LYMNev <- subset(LYM, subset = Smoking == "Never")
LYMSmo <- subset(LYM, subset = Smoking == "Smoker")

gene = "MALAT1"
n <- unique(LYM@meta.data$LYMClass)
Res <- list()
for(i in 1:length(n)){
	Cell <- n[i]
	Nev <- subset(LYMNev, subset = LYMClass == Cell)
	Smo <- subset(LYMSmo, subset = LYMClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev, features = gene)
	Idents(Smo) <- Smo@meta.data$orig.ident
	Smogene <- AverageExpression(Smo, features = gene)

	NRes <- data.frame(Smoking = "Never", Cell = Cell, Nevgene)
	SRes <- data.frame(Smoking = "Smoker", Cell = Cell, Smogene)

	Res <- bind_rows(Res, NRes, SRes)
}

write.table(Res, file = "Age-LYM-MALAT1-regression.csv", sep = ",")







MYENev <- subset(MYE, subset = Smoking == "Never")
MYESmo <- subset(MYE, subset = Smoking == "Smoker")

gene = "FTL"
n <- unique(MYE@meta.data$MYEClass)
Res <- list()
for(i in 1:length(n)){
	Cell <- n[i]
	Nev <- subset(MYENev, subset = MYEClass == Cell)
	Smo <- subset(MYESmo, subset = MYEClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev, features = gene)
	Idents(Smo) <- Smo@meta.data$orig.ident
	Smogene <- AverageExpression(Smo, features = gene)

	NRes <- data.frame(Smoking = "Never", Cell = Cell, Nevgene)
	SRes <- data.frame(Smoking = "Smoker", Cell = Cell, Smogene)

	Res <- bind_rows(Res, NRes, SRes)
}

write.table(Res, file = "Age-MYE-FTL-regression.csv", sep = ",")


gene = "TMSB4X"
n <- unique(MYE@meta.data$MYEClass)
Res <- list()
for(i in 1:length(n)){
	Cell <- n[i]
	Nev <- subset(MYENev, subset = MYEClass == Cell)
	Smo <- subset(MYESmo, subset = MYEClass == Cell)

	Idents(Nev) <- Nev@meta.data$orig.ident
	Nevgene <- AverageExpression(Nev, features = gene)
	Idents(Smo) <- Smo@meta.data$orig.ident
	Smogene <- AverageExpression(Smo, features = gene)

	NRes <- data.frame(Smoking = "Never", Cell = Cell, Nevgene)
	SRes <- data.frame(Smoking = "Smoker", Cell = Cell, Smogene)

	Res <- bind_rows(Res, NRes, SRes)
}

write.table(Res, file = "Age-MYE-TMSB4X-regression.csv", sep = ",")











