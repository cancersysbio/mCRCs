mCRC_dir <- Sys.getenv("mCRC_DIR")
setwd(file.path(mCRC_dir,"results"))

genes1 <- c("APC","KRAS","TP53","SMAD4","PIK3CA","PTEN","ATM","TCF7L2",
	 "PIK3R1","ERBB3","AXIN2","FN1","DNAH8","APOB","PTPRT","SYNE1")

genes2 <- c("APC","TP53","KRAS","PIK3CA","FBXW7","SMAD4","NRAS","TCF7L2",
	  "FAM123B","SMAD2","CTNNB1","KIAA1804","SOX9","ACVR1B","GPC6","EDNRB")

short_list <- unique(c(genes1,genes2))

tcga1 <- read.delim("../Genelists/coadread_630_tcga_Mutated Genes.txt",as.is=TRUE)
tcga2 <- read.delim("../Genelists/coadread_630_tcga_Copy Number Altered Genes.txt",as.is=TRUE)
intogen1 <- read.delim("../Genelists/intogen-COREAD-229-drivers-data.tsv",as.is=TRUE)

tcga1_100 <- tcga1[!is.na(tcga1$MutSig) & tcga1$X. >= 2,1:5]
tcga1_100 <- tcga1_100[ order(tcga1_100$MutSig),]
tcga1_100 <- tcga1_100[ tcga1_100$MutSig <= 1e-3,]

new_genes <- unique(c(tcga1_100$Gene,tcga2$Gene[!is.na(tcga2$Gistic)],
		      intogen1$SYMBOL))

short_list <- unique(c(short_list,new_genes))

cat(short_list,file="CRCGenes.txt")
