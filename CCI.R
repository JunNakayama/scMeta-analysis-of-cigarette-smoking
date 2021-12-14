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



### CCI
library(circlize)


## EPI-MYE
## Select the Immune-interactions
CCI_list <- read_csv("CCI_list.csv")
CCIgene <- read_csv("CCIgene.csv")

Chemo = CCI_list[CCI_list$AliasA %in% unlist(CCIgene),]
uni.chemo = list(Chemo[,1], Chemo[,2])
uni.chemo = unique(unlist(uni.chemo))
EPICCI = EPI[uni.chemo,]
MYECCI = MYE[uni.chemo,]
Idents(EPICCI) <- EPICCI@meta.data$EPIClass
Idents(MYECCI) <- MYECCI@meta.data$MYEClass
CCImerge = merge(MYECCI, c(EPICCI))
CCImerge.store = CCImerge


#######################
CCImerge = CCImerge.store
smoker = subset(CCImerge, subset = Smoking == "Smoker")
Never = subset(CCImerge, subset = Smoking == "Never")

#######################
CCImerge =  Never
unigenes = rownames(CCImerge)
mel = list()
for(i in 1: length(unigenes)){

	expr <- FetchData(object = CCImerge, vars = unigenes[i])

	if(any(expr > 2) >= 1){

		ge = CCImerge[, which(x = expr > 2)]
		ge = table(ge@active.ident)
		ho = melt(ge)
		ho = mutate(ho, Gene = unigenes[i])
		mel = rbind(mel, ho)

	}else{
		print(unigenes[i])
	}
	
}
tocell = table(CCImerge@active.ident)
Ratiomel = list()
for(i in 1:length(tocell)){

	cell = names(tocell[i])
	LS = mel[mel[,1] == cell,]
	LS[,2] = LS[,2]/tocell[i]
	Ratiomel = rbind(Ratiomel, LS)

}
Ratiomel = na.omit(Ratiomel)
Ratiomel = reshape2::dcast(Ratiomel, Var1 ~ Gene, value.var = "value")
write.table(Ratiomel, file = "CCI-RatioMAT-NeverMYE-IMMgene.txt", sep = "	")

ratiomelt = na.omit(melt(Ratiomel))
ratiomelt = ratiomelt[ratiomelt$value > 0.1,]
CCI_perlist = Chemo[Chemo$AliasA %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = CCI_perlist[CCI_perlist$AliasB %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = mutate(CCI_perlist, merge = paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "&"))
a = as.character(unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[CCI_perlist$merge %in% unique(CCI_perlist$merge),]
CCI_perlist = subset(CCI_perlist, CCI_perlist$merge %in% a)
a = unique(CCI_perlist$merge)
aa = str_split(a, "&", n = 2, simplify = TRUE)
CCI_perlist = data.frame(L = aa[,1], R = aa[,2], merge = unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[!(as.character(CCI_perlist$L) == as.character(CCI_perlist$R)),]

CCItotal = list()
for(i in 1:nrow(CCI_perlist)){

	geneset = CCI_perlist[i,]
	Ligand = as.character(unlist(geneset[1]))
	Receptor = as.character(unlist(geneset[2]))

	LigCCI = ratiomelt[ratiomelt$variable == Ligand, ]
	ResCCI = ratiomelt[ratiomelt$variable == Receptor, ]

	for(k in 1:nrow(LigCCI)){

		Lig = LigCCI[k,3]
		score = lapply(ResCCI$value, function(x){x * Lig})
		Res = data.frame( CCIscore = unlist(score), Recepter.Cell = ResCCI$Var1)
		Res = mutate(Res, Ligand.Cell = LigCCI[k,1], CCI = geneset)

		CCItotal = bind_rows(CCItotal, Res)
	}

}


CCItotal = CCItotal[CCItotal$Recepter.Cell != c("MT-1") & CCItotal$Ligand.Cell != c("MT-1"),]
CCItotal = CCItotal[CCItotal$Recepter.Cell != c("MT-2") & CCItotal$Ligand.Cell != c("MT-2"),]
GroupA = c(unlist(as.character(unique(EPICCI@active.ident))))
GroupB = c(unlist(as.character(unique(MYECCI@active.ident))))
AtoB = CCItotal[CCItotal$Ligand.Cell %in% GroupA & CCItotal$Recepter.Cell %in% GroupB,]
BtoA = CCItotal[CCItotal$Ligand.Cell %in% GroupB & CCItotal$Recepter.Cell %in% GroupA,]

LigA = unlist(as.character(unique(AtoB$Ligand.Cell)))
AtoB.CCI = list()
for(i in 1:length(LigA)){

	G <- LigA[i]
	set = AtoB[AtoB$Ligand.Cell == G, ]

	CELLset = unique(set$Recepter.Cell)

	for(k in 1:length(CELLset)){

		R <- CELLset[k]
		REset <- set[set$Recepter.Cell == R, ]
		CCIsum <- sum(REset$CCIscore)

		DF = data.frame(LigandCell = G, ReceptorCell = R, CCIscore.sum = CCIsum)

		AtoB.CCI = rbind(AtoB.CCI, DF)

	}

}

LigB = unlist(as.character(unique(BtoA$Ligand.Cell)))
BtoA.CCI = list()
for(i in 1:length(LigB)){

	G <- LigB[i]
	set = BtoA[BtoA$Ligand.Cell == G, ]

	CELLset = unique(set$Recepter.Cell)

	for(k in 1:length(CELLset)){

		R <- CELLset[k]
		REset <- set[set$Recepter.Cell == R, ]
		CCIsum <- sum(REset$CCIscore)

		DF = data.frame(LigandCell = G, ReceptorCell = R, CCIscore.sum = CCIsum)

		BtoA.CCI = rbind(BtoA.CCI, DF)

	}

}

write.table(AtoB, file = "CCI-interaction-score-AtoB-NeverMYE.txt", sep = "	")
write.table(BtoA, file = "CCI-interaction-score-BtoA-NeverMYE.txt", sep = "	")
write.table(AtoB.CCI, file = "CCI-interaction-score-AtoB-CELL-NeverMYE.txt", sep = "	")
write.table(BtoA.CCI, file = "CCI-interaction-score-BtoA-CELL-NeverMYE.txt", sep = "	")


#######################
#######################
CCImerge = smoker
unigenes = rownames(CCImerge)
mel = list()
for(i in 1: length(unigenes)){

	expr <- FetchData(object = CCImerge, vars = unigenes[i])

	if(any(expr > 2) >= 1){

		ge = CCImerge[, which(x = expr > 2)]
		ge = table(ge@active.ident)
		ho = melt(ge)
		ho = mutate(ho, Gene = unigenes[i])
		mel = rbind(mel, ho)

	}else{
		print(unigenes[i])
	}
	
}
tocell = table(CCImerge@active.ident)
Ratiomel = list()
for(i in 1:length(tocell)){

	cell = names(tocell[i])
	LS = mel[mel[,1] == cell,]
	LS[,2] = LS[,2]/tocell[i]
	Ratiomel = rbind(Ratiomel, LS)

}
Ratiomel = na.omit(Ratiomel)
Ratiomel = reshape2::dcast(Ratiomel, Var1 ~ Gene, value.var = "value")
write.table(Ratiomel, file = "CCI-RatioMAT-smokerMYE-IMMgener.txt", sep = "	")

A = list()
B = list()
CCItotal = list()
CCI_perlist = Chemo[Chemo$AliasA %in% unlist(unique(mel[,3])),]
CCI_perlist = CCI_perlist[CCI_perlist$AliasB %in% unlist(unique(mel[,3])),]
CCI_perlist = mutate(CCI_perlist, merge = paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "|"))
unique(list(paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "|")))
ratiomelt = na.omit(melt(Ratiomel))
ratiomelt = ratiomelt[ratiomelt$value > 0.1,]
CCI_perlist = Chemo[Chemo$AliasA %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = CCI_perlist[CCI_perlist$AliasB %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = mutate(CCI_perlist, merge = paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "&"))
a = as.character(unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[CCI_perlist$merge %in% unique(CCI_perlist$merge),]
CCI_perlist = subset(CCI_perlist, CCI_perlist$merge %in% a)
a = unique(CCI_perlist$merge)
aa = str_split(a, "&", n = 2, simplify = TRUE)
CCI_perlist = data.frame(L = aa[,1], R = aa[,2], merge = unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[!(as.character(CCI_perlist$L) == as.character(CCI_perlist$R)),]

CCItotal = list()
for(i in 1:nrow(CCI_perlist)){

	geneset = CCI_perlist[i,]
	Ligand = as.character(unlist(geneset[1]))
	Receptor = as.character(unlist(geneset[2]))

	LigCCI = ratiomelt[ratiomelt$variable == Ligand, ]
	ResCCI = ratiomelt[ratiomelt$variable == Receptor, ]

	for(k in 1:nrow(LigCCI)){

		Lig = LigCCI[k,3]
		score = lapply(ResCCI$value, function(x){x * Lig})
		Res = data.frame( CCIscore = unlist(score), Recepter.Cell = ResCCI$Var1)
		Res = mutate(Res, Ligand.Cell = LigCCI[k,1], CCI = geneset)

		CCItotal = bind_rows(CCItotal, Res)
	}

}
CCItotal = CCItotal[CCItotal$Recepter.Cell != c("MT-1") & CCItotal$Ligand.Cell != c("MT-1"),]
CCItotal = CCItotal[CCItotal$Recepter.Cell != c("MT-2") & CCItotal$Ligand.Cell != c("MT-2"),]
GroupA = c(unlist(as.character(unique(EPICCI@active.ident))))
GroupB = c(unlist(as.character(unique(MYECCI@active.ident))))
AtoB = CCItotal[CCItotal$Ligand.Cell %in% GroupA & CCItotal$Recepter.Cell %in% GroupB,]
BtoA = CCItotal[CCItotal$Ligand.Cell %in% GroupB & CCItotal$Recepter.Cell %in% GroupA,]

LigA = unlist(as.character(unique(AtoB$Ligand.Cell)))
AtoB.CCI = list()
for(i in 1:length(LigA)){

	G <- LigA[i]
	set = AtoB[AtoB$Ligand.Cell == G, ]

	CELLset = unique(set$Recepter.Cell)

	for(k in 1:length(CELLset)){

		R <- CELLset[k]
		REset <- set[set$Recepter.Cell == R, ]
		CCIsum <- sum(REset$CCIscore)

		DF = data.frame(LigandCell = G, ReceptorCell = R, CCIscore.sum = CCIsum)

		AtoB.CCI = rbind(AtoB.CCI, DF)

	}

}

LigB = unlist(as.character(unique(BtoA$Ligand.Cell)))
BtoA.CCI = list()
for(i in 1:length(LigB)){

	G <- LigB[i]
	set = BtoA[BtoA$Ligand.Cell == G, ]

	CELLset = unique(set$Recepter.Cell)

	for(k in 1:length(CELLset)){

		R <- CELLset[k]
		REset <- set[set$Recepter.Cell == R, ]
		CCIsum <- sum(REset$CCIscore)

		DF = data.frame(LigandCell = G, ReceptorCell = R, CCIscore.sum = CCIsum)

		BtoA.CCI = rbind(BtoA.CCI, DF)

	}

}

write.table(AtoB, file = "CCI-interaction-score-AtoBsmokerMYE.txt", sep = "	")
write.table(BtoA, file = "CCI-interaction-score-BtoAsmokerMYE.txt", sep = "	")
write.table(AtoB.CCI, file = "CCI-interaction-score-AtoB-CELLsmokerMYE.txt", sep = "	")
write.table(BtoA.CCI, file = "CCI-interaction-score-BtoA-CELLsmokerMYE.txt", sep = "	")






#######################
#######################
#######################
#######################
### EPI-LYM
Chemo = CCI_list[CCI_list$AliasA %in% unlist(CCIgene),]
uni.chemo = list(Chemo[,1], Chemo[,2])
uni.chemo = unique(unlist(uni.chemo))
EPICCI = EPI[uni.chemo,]
LYMCCI = LYM[uni.chemo,]
Idents(EPICCI) <- EPICCI@meta.data$EPIClass
Idents(LYMCCI) <- LYMCCI@meta.data$LYMClass
CCImerge = merge(LYMCCI, c(EPICCI))
CCImerge.store = CCImerge


#######################
CCImerge = CCImerge.store
smoker = subset(CCImerge, subset = Smoking == "Smoker")
Never = subset(CCImerge, subset = Smoking == "Never")
#######################
CCImerge =  Never
unigenes = rownames(CCImerge)
mel = list()
for(i in 1: length(unigenes)){

	expr <- FetchData(object = CCImerge, vars = unigenes[i])

	if(any(expr > 2) >= 1){

		ge = CCImerge[, which(x = expr > 2)]
		ge = table(ge@active.ident)
		ho = melt(ge)
		ho = mutate(ho, Gene = unigenes[i])
		mel = rbind(mel, ho)

	}else{
		print(unigenes[i])
	}
	
}
tocell = table(CCImerge@active.ident)
Ratiomel = list()
for(i in 1:length(tocell)){

	cell = names(tocell[i])
	LS = mel[mel[,1] == cell,]
	LS[,2] = LS[,2]/tocell[i]
	Ratiomel = rbind(Ratiomel, LS)

}
Ratiomel = na.omit(Ratiomel)
Ratiomel = reshape2::dcast(Ratiomel, Var1 ~ Gene, value.var = "value")
write.table(Ratiomel, file = "CCI-RatioMAT-NeverLYM-IMMgene.txt", sep = "	")

ratiomelt = na.omit(melt(Ratiomel))
ratiomelt = ratiomelt[ratiomelt$value > 0.1,]
CCI_perlist = Chemo[Chemo$AliasA %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = CCI_perlist[CCI_perlist$AliasB %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = mutate(CCI_perlist, merge = paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "&"))
a = as.character(unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[CCI_perlist$merge %in% unique(CCI_perlist$merge),]
CCI_perlist = subset(CCI_perlist, CCI_perlist$merge %in% a)
a = unique(CCI_perlist$merge)
aa = str_split(a, "&", n = 2, simplify = TRUE)
CCI_perlist = data.frame(L = aa[,1], R = aa[,2], merge = unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[!(as.character(CCI_perlist$L) == as.character(CCI_perlist$R)),]

CCItotal = list()
for(i in 1:nrow(CCI_perlist)){

	geneset = CCI_perlist[i,]
	Ligand = as.character(unlist(geneset[1]))
	Receptor = as.character(unlist(geneset[2]))

	LigCCI = ratiomelt[ratiomelt$variable == Ligand, ]
	ResCCI = ratiomelt[ratiomelt$variable == Receptor, ]

	for(k in 1:nrow(LigCCI)){

		Lig = LigCCI[k,3]
		score = lapply(ResCCI$value, function(x){x * Lig})
		Res = data.frame( CCIscore = unlist(score), Recepter.Cell = ResCCI$Var1)
		Res = mutate(Res, Ligand.Cell = LigCCI[k,1], CCI = geneset)

		CCItotal = bind_rows(CCItotal, Res)
	}

}
CCItotal = CCItotal[CCItotal$Recepter.Cell != c("MT-1") & CCItotal$Ligand.Cell != c("MT-1"),]
CCItotal = CCItotal[CCItotal$Recepter.Cell != c("MT-2") & CCItotal$Ligand.Cell != c("MT-2"),]
GroupA = c(unlist(as.character(unique(EPICCI@active.ident))))
GroupB = c(unlist(as.character(unique(LYMCCI@active.ident))))
AtoB = CCItotal[CCItotal$Ligand.Cell %in% GroupA & CCItotal$Recepter.Cell %in% GroupB,]
BtoA = CCItotal[CCItotal$Ligand.Cell %in% GroupB & CCItotal$Recepter.Cell %in% GroupA,]

LigA = unlist(as.character(unique(AtoB$Ligand.Cell)))
AtoB.CCI = list()
for(i in 1:length(LigA)){

	G <- LigA[i]
	set = AtoB[AtoB$Ligand.Cell == G, ]

	CELLset = unique(set$Recepter.Cell)

	for(k in 1:length(CELLset)){

		R <- CELLset[k]
		REset <- set[set$Recepter.Cell == R, ]
		CCIsum <- sum(REset$CCIscore)

		DF = data.frame(LigandCell = G, ReceptorCell = R, CCIscore.sum = CCIsum)

		AtoB.CCI = rbind(AtoB.CCI, DF)

	}

}

LigB = unlist(as.character(unique(BtoA$Ligand.Cell)))
BtoA.CCI = list()
for(i in 1:length(LigB)){

	G <- LigB[i]
	set = BtoA[BtoA$Ligand.Cell == G, ]

	CELLset = unique(set$Recepter.Cell)

	for(k in 1:length(CELLset)){

		R <- CELLset[k]
		REset <- set[set$Recepter.Cell == R, ]
		CCIsum <- sum(REset$CCIscore)

		DF = data.frame(LigandCell = G, ReceptorCell = R, CCIscore.sum = CCIsum)

		BtoA.CCI = rbind(BtoA.CCI, DF)

	}

}

write.table(AtoB, file = "CCI-interaction-score-AtoB-NeverLYM.txt", sep = "	")
write.table(BtoA, file = "CCI-interaction-score-BtoA-NeverLYM.txt", sep = "	")
write.table(AtoB.CCI, file = "CCI-interaction-score-AtoB-CELL-NeverLYM.txt", sep = "	")
write.table(BtoA.CCI, file = "CCI-interaction-score-BtoA-CELL-NeverLYM.txt", sep = "	")


#######################
#######################
CCImerge = smoker
unigenes = rownames(CCImerge)
mel = list()
for(i in 1: length(unigenes)){

	expr <- FetchData(object = CCImerge, vars = unigenes[i])

	if(any(expr > 2) >= 1){

		ge = CCImerge[, which(x = expr > 2)]
		ge = table(ge@active.ident)
		ho = melt(ge)
		ho = mutate(ho, Gene = unigenes[i])
		mel = rbind(mel, ho)

	}else{
		print(unigenes[i])
	}
	
}
tocell = table(CCImerge@active.ident)
Ratiomel = list()
for(i in 1:length(tocell)){

	cell = names(tocell[i])
	LS = mel[mel[,1] == cell,]
	LS[,2] = LS[,2]/tocell[i]
	Ratiomel = rbind(Ratiomel, LS)

}
Ratiomel = na.omit(Ratiomel)
Ratiomel = reshape2::dcast(Ratiomel, Var1 ~ Gene, value.var = "value")
write.table(Ratiomel, file = "CCI-RatioMAT-smokerLYM-IMMgener.txt", sep = "	")

A = list()
B = list()
CCItotal = list()
CCI_perlist = Chemo[Chemo$AliasA %in% unlist(unique(mel[,3])),]
CCI_perlist = CCI_perlist[CCI_perlist$AliasB %in% unlist(unique(mel[,3])),]
CCI_perlist = mutate(CCI_perlist, merge = paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "|"))
unique(list(paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "|")))
ratiomelt = na.omit(melt(Ratiomel))
ratiomelt = ratiomelt[ratiomelt$value > 0.1,]
CCI_perlist = Chemo[Chemo$AliasA %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = CCI_perlist[CCI_perlist$AliasB %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = mutate(CCI_perlist, merge = paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "&"))
a = as.character(unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[CCI_perlist$merge %in% unique(CCI_perlist$merge),]
CCI_perlist = subset(CCI_perlist, CCI_perlist$merge %in% a)
a = unique(CCI_perlist$merge)
aa = str_split(a, "&", n = 2, simplify = TRUE)
CCI_perlist = data.frame(L = aa[,1], R = aa[,2], merge = unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[!(as.character(CCI_perlist$L) == as.character(CCI_perlist$R)),]

CCItotal = list()
for(i in 1:nrow(CCI_perlist)){

	geneset = CCI_perlist[i,]
	Ligand = as.character(unlist(geneset[1]))
	Receptor = as.character(unlist(geneset[2]))

	LigCCI = ratiomelt[ratiomelt$variable == Ligand, ]
	ResCCI = ratiomelt[ratiomelt$variable == Receptor, ]

	for(k in 1:nrow(LigCCI)){

		Lig = LigCCI[k,3]
		score = lapply(ResCCI$value, function(x){x * Lig})
		Res = data.frame( CCIscore = unlist(score), Recepter.Cell = ResCCI$Var1)
		Res = mutate(Res, Ligand.Cell = LigCCI[k,1], CCI = geneset)

		CCItotal = bind_rows(CCItotal, Res)
	}

}
CCItotal = CCItotal[CCItotal$Recepter.Cell != c("MT-1") & CCItotal$Ligand.Cell != c("MT-1"),]
CCItotal = CCItotal[CCItotal$Recepter.Cell != c("MT-2") & CCItotal$Ligand.Cell != c("MT-2"),]
GroupA = c(unlist(as.character(unique(EPICCI@active.ident))))
GroupB = c(unlist(as.character(unique(LYMCCI@active.ident))))
AtoB = CCItotal[CCItotal$Ligand.Cell %in% GroupA & CCItotal$Recepter.Cell %in% GroupB,]
BtoA = CCItotal[CCItotal$Ligand.Cell %in% GroupB & CCItotal$Recepter.Cell %in% GroupA,]

LigA = unlist(as.character(unique(AtoB$Ligand.Cell)))
AtoB.CCI = list()
for(i in 1:length(LigA)){

	G <- LigA[i]
	set = AtoB[AtoB$Ligand.Cell == G, ]

	CELLset = unique(set$Recepter.Cell)

	for(k in 1:length(CELLset)){

		R <- CELLset[k]
		REset <- set[set$Recepter.Cell == R, ]
		CCIsum <- sum(REset$CCIscore)

		DF = data.frame(LigandCell = G, ReceptorCell = R, CCIscore.sum = CCIsum)

		AtoB.CCI = rbind(AtoB.CCI, DF)

	}

}

LigB = unlist(as.character(unique(BtoA$Ligand.Cell)))
BtoA.CCI = list()
for(i in 1:length(LigB)){

	G <- LigB[i]
	set = BtoA[BtoA$Ligand.Cell == G, ]

	CELLset = unique(set$Recepter.Cell)

	for(k in 1:length(CELLset)){

		R <- CELLset[k]
		REset <- set[set$Recepter.Cell == R, ]
		CCIsum <- sum(REset$CCIscore)

		DF = data.frame(LigandCell = G, ReceptorCell = R, CCIscore.sum = CCIsum)

		BtoA.CCI = rbind(BtoA.CCI, DF)

	}

}

write.table(AtoB, file = "CCI-interaction-score-AtoBsmokerLYM.txt", sep = "	")
write.table(BtoA, file = "CCI-interaction-score-BtoAsmokerLYM.txt", sep = "	")
write.table(AtoB.CCI, file = "CCI-interaction-score-AtoB-CELLsmokerLYM.txt", sep = "	")
write.table(BtoA.CCI, file = "CCI-interaction-score-BtoA-CELLsmokerLYM.txt", sep = "	")




