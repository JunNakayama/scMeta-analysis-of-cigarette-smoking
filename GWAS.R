library(Seurat)
library(harmony)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(MAST)
library(DoubletFinder)


gene = c("AJAP1", "NPHP4", "FUBP", "DNAJB4", "MUC1", "ADAM15", "THBS3", "NRXN1", "NUP35", "HIBCH", "INPP1", "PMS1", "STAT1", "TP63", "C3ORF21", "KCNIP4", "TERT",
		"CLPTM1L", "XRCC4", "PAHA2", "CSF2", "IL3", "SCL22A5", "ASCL6", "STK32A", "PPP2R2B", "DPYSL3", "HIST1H1E", "BAG6", "APOM", "TNXB", "MSH5", "BTNL2", "PRRC2A",
		"FKBPL", "HSPA1B", "FOXP4", "FOXP4-AS1", "GTF2H4", "LRFN2", "HLA-A", "HLA-DQB1", "DCBLD1", "ROS1", "RNASET2", "SP4", "DNAH11", "EPHX2", "CHRNA2", "NRG1",
		"CDKN2A", "CDKN2B", "CDKN2B-AS1", "MTAP", "GATA3", "FFAR4", "OBFC1", "VT11A", "MPZL3", "AMICA1", "RAD52", "ACVR1B", "NR4A1", "NR1H4", "SLC17A8", "SH2B3",
		"MIPEP", "TNFRSF19", "BRCA2", "GPC5", "SEMA6D", "SECISBP2L", "CHRNA5", "CHRNA3", "CHRNB4", "IREB2", "PSMA4", "HYKK", "BPTF", "FAM38B", "APCDD1", "NAPG",
		"GAREM", "TGFB1", "CYP2A6", "BPIFB1", "CYP24A1", "RTEL1", "CHEK2", "LIF", "HORMAD2", "MTMR3")

GWAS = EPI[gene,]
Gender <- unique(GWAS@meta.data$Gender)
matlist = list()

for(i in 1:length(Gender)){
	GEN <- Gender[i]
	cl <- subset(GWAS, subset = Gender == GEN)
	Idents(object = cl) <- cl@meta.data$EPIClass
	avelist = AverageExpression(cl, add.ident = "Smoking")
	avelist = data.frame(avelist$RNA)
	matlist = bind_cols(matlist, avelist, Gender = GEN)

}

write.table(matlist, "GWAS-EPI.csv", sep = ",")



GWAS = EPI[gene,]
Gender <- unique(GWAS@meta.data$Gender)
Res = list()
for(i in 1:length(Gender)){
	GEN <- Gender[i]
	cl <- subset(GWAS, subset = Gender == GEN)
	Idents(object = cl) <- cl@meta.data$EPIClass
	avelist = AverageExpression(cl, add.ident = "Smoking")
	avelist = data.frame(avelist$RNA, Genes = rownames(avelist$RNA))
	matlist <- reshape2::melt(avelist)
	matlist <- data.frame(matlist, Gender = GEN)
	Res <- bind_rows(Res, matlist)
}
GWAS = STR[gene,]
Gender <- unique(GWAS@meta.data$Gender)
for(i in 1:length(Gender)){
	GEN <- Gender[i]
	cl <- subset(GWAS, subset = Gender == GEN)
	Idents(object = cl) <- cl@meta.data$STRClass
	avelist = AverageExpression(cl, add.ident = "Smoking")
	avelist = data.frame(avelist$RNA, Genes = rownames(avelist$RNA))
	matlist <- reshape2::melt(avelist)
	matlist <- data.frame(matlist, Gender = GEN)
	Res <- bind_rows(Res, matlist)
}
GWAS = END[gene,]
Gender <- unique(GWAS@meta.data$Gender)
for(i in 1:length(Gender)){
	GEN <- Gender[i]
	cl <- subset(GWAS, subset = Gender == GEN)
	Idents(object = cl) <- cl@meta.data$ENDClass
	avelist = AverageExpression(cl, add.ident = "Smoking")
	avelist = data.frame(avelist$RNA, Genes = rownames(avelist$RNA))
	matlist <- reshape2::melt(avelist)
	matlist <- data.frame(matlist, Gender = GEN)
	Res <- bind_rows(Res, matlist)
}
GWAS = LYM[gene,]
Gender <- unique(GWAS@meta.data$Gender)
for(i in 1:length(Gender)){
	GEN <- Gender[i]
	cl <- subset(GWAS, subset = Gender == GEN)
	Idents(object = cl) <- cl@meta.data$LYMClass
	avelist = AverageExpression(cl, add.ident = "Smoking")
	avelist = data.frame(avelist$RNA, Genes = rownames(avelist$RNA))
	matlist <- reshape2::melt(avelist)
	matlist <- data.frame(matlist, Gender = GEN)
	Res <- bind_rows(Res, matlist)
}
GWAS = MYE[gene,]
Gender <- unique(GWAS@meta.data$Gender)
for(i in 1:length(Gender)){
	GEN <- Gender[i]
	cl <- subset(GWAS, subset = Gender == GEN)
	Idents(object = cl) <- cl@meta.data$MYEClass
	avelist = AverageExpression(cl, add.ident = "Smoking")
	avelist = data.frame(avelist$RNA, Genes = rownames(avelist$RNA))
	matlist <- reshape2::melt(avelist)
	matlist <- data.frame(matlist, Gender = GEN)
	Res <- bind_rows(Res, matlist)
}


write.table(Res, "GWAS-all.csv", sep = ",")
