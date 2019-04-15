mCRC_dir <- Sys.getenv("mCRC_DIR")
setwd(file.path(mCRC_dir,"results"))

samples <- read.delim("alltumorlabels.txt",as.is=TRUE)

load("mCRC_All_TitanCNA.RData")

purity <- data.frame(Sample=samples$Label,Purity=rep(0,nrow(samples)),
		     stringsAsFactors = FALSE)

for (i in 1:nrow(samples)) {
  sn1 <- samples$Bam[i]
  cna1 <- alltitancnaresults[[sn1]]
  purity$Purity[i] <- round((1 - tail(cna1$convergeParams$n,1)) *
			      (1 - tail(cna1$convergeParams$s[1,],1)),3)
}

write.table(purity,file="mCRC_All_Purity.txt",
	    sep="\t",quote=FALSE,row.names=FALSE)

patients <- sapply(strsplit(purity$Sample,"_",fixed=TRUE),head,1)
tissues <- sapply(strsplit(purity$Sample,"_",fixed=TRUE),"[",2)

pts1 <- c("V46","V402","V514","V559","V750","V824","V855","V930","V953","V974","TB1","TB7")
pts2 <- c("V402","V750","V824","V930","V953","V974","TB1","TB7")

spaces <- rep(0.25,nrow(samples))
spaces[which(patients[-1] != patients[-length(patients)])+1] <- 1.25

pdf("mCRC_All_Purity.pdf",height=5,width=10)

par(mar=c(5,4,2,2))
idx1 <- which(patients %in% pts1)
barplot(purity$Purity[idx1],width=0.8,space=spaces[idx1],##col=cols,
	names.arg=purity$Sample[idx1],ylim=c(0,1),
	las=2,ylab="Purity",cex.names=0.6)

idx1 <- which(patients %in% pts2)
barplot(purity$Purity[idx1],width=0.8,space=spaces[idx1],##col=cols,
	names.arg=purity$Sample[idx1],ylim=c(0,1),
	las=2,ylab="Purity",cex.names=0.7)

dev.off()
