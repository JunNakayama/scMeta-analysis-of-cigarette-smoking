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
Res = list()
	Idents(object = GWAS) <- GWAS@meta.data$EPIClass
	avelist = AverageExpression(GWAS, add.ident = "Smoking")
	avelist = data.frame(avelist$SCT, Genes = rownames(avelist$SCT))
	matlist <- reshape2::melt(avelist)
	matlist <- data.frame(matlist)
	Res <- bind_rows(Res, matlist)
GWAS = STR[gene,]
	Idents(object = GWAS) <- GWAS@meta.data$STRClass
	avelist = AverageExpression(GWAS, add.ident = "Smoking")
	avelist = data.frame(avelist$SCT, Genes = rownames(avelist$SCT))
	matlist <- reshape2::melt(avelist)
	matlist <- data.frame(matlist)
	Res <- bind_rows(Res, matlist)
GWAS = END[gene,]
	Idents(object = GWAS) <- GWAS@meta.data$ENDClass
	avelist = AverageExpression(GWAS, add.ident = "Smoking")
	avelist = data.frame(avelist$SCT, Genes = rownames(avelist$SCT))
	matlist <- reshape2::melt(avelist)
	matlist <- data.frame(matlist)
	Res <- bind_rows(Res, matlist)
GWAS = LYM[gene,]
	Idents(object = GWAS) <- GWAS@meta.data$LYMClass
	avelist = AverageExpression(GWAS, add.ident = "Smoking")
	avelist = data.frame(avelist$SCT, Genes = rownames(avelist$SCT))
	matlist <- reshape2::melt(avelist)
	matlist <- data.frame(matlist)
	Res <- bind_rows(Res, matlist)
GWAS = MYE[gene,]
	Idents(object = GWAS) <- GWAS@meta.data$MYEClass
	avelist = AverageExpression(GWAS, add.ident = "Smoking")
	avelist = data.frame(avelist$SCT, Genes = rownames(avelist$SCT))
	matlist <- reshape2::melt(avelist)
	matlist <- data.frame(matlist)
	Res <- bind_rows(Res, matlist)


write.table(Res, "GWAS-all-withoutGender.csv", sep = ",")


mat = reshape2::dcast(Res, variable ~ Genes, value.var = "value")

write.table(mat, "GWAS-all-withoutGenderMatrix.csv", sep = ",")



## Boxplot
hoge <- VlnPlot(EPI, "MUC1", group.by = "EPIClass", split.by = "Smoking") + NoLegend()
data <- data.frame(Gene = hoge$data$MUC1,ident = hoge$data$ident, Smoking = hoge$data$split)
data <- data[order(data$ident, decreasing = FALSE),]

write.table(data, "GWAS-EPI-MUC1.csv", sep = ",")


## Boxplot
hoge <- VlnPlot(MYE, "HLA-A", group.by = "MYEClass", split.by = "Smoking") + NoLegend()
data <- data.frame(Gene = hoge$data$`HLA-A`,ident = hoge$data$ident, Smoking = hoge$data$split)
data <- data[order(data$ident, decreasing = FALSE),]

write.table(data, "GWAS-MYE-HLA-A.csv", sep = ",")
