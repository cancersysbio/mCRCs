mCRC_dir <- Sys.getenv("mCRC_DIR")
setwd(file.path(mCRC_dir,"results"))

dir1 <- file.path(mCRC_dir,"bams/all/")
samples <- read.delim("alltumorlabels.txt",as.is=TRUE)

allmediancov <- matrix(0,nrow=nrow(samples),ncol=4)
colnames(allmediancov) <- c("All","RmDup","DupRatio","MisMatch")
rownames(allmediancov) <- samples$Label

for (i in 1:nrow(samples)) {
  sn1 <- samples$Bam[i]

  fn1 <- paste0(dir1,sn1,".bam.cov")
  fn2 <- paste0(dir1,sn1,".bam.rmdup.cov")
  fn3 <- paste0(dir1,sn1,".bam.stats")

  if (file.exists(fn1)) {
    coverage <- read.table(fn1)
    allmediancov[i,1] <- coverage$V2[ which(cumsum(coverage$V5) >= 0.5)[1] ]
    coverage <- read.table(fn2)
    allmediancov[i,2] <- coverage$V2[ which(cumsum(coverage$V5) >= 0.5)[1] ]
    
    a <- read.table(fn3,sep="\t",fill=TRUE,skip=6,nrow=22)
    allmediancov[i,3] <- round(a[12,3]/a[7,3],3)
    allmediancov[i,4] <- signif(a[22,3],3)
  }
}

write.table(allmediancov,file="mCRC_All_Coverage.txt",
	    sep="\t",quote=FALSE)

patients <- sapply(strsplit(rownames(allmediancov),"_",fixed=TRUE),"[",1)
pts1 <- c("V46","V402","V514","V559","V750","V824","V855","V930","V953","V974","TB1","TB7")
pts2 <- c("V402","V750","V824","V930","V953","V974","TB1","TB7")

spaces <- rep(0.25,nrow(allmediancov))
spaces[which(patients[-1] != patients[-length(patients)])+1] <- 1.25

allmediancov[,1] <- allmediancov[,1] - allmediancov[,2]

pdf("mCRC_All_Coverage.pdf",height=6,width=9)
par(mfrow=c(2,1))

idx1 <- patients %in% pts1
par(mar=c(3,4,1,0))
barplot(t(allmediancov[idx1,2:1]),width=0.8,space=spaces[idx1],##col=cols,
	legend.text=c("Removing duplicates","All"),main="",
	las=2,ylab="Median depth",cex.names=0.6)
par(mar=c(1,4,3,0))
barplot(allmediancov[idx1,3],names.arg=FALSE,ylim=c(1,0),##col=cols,
	ylab="Duplication rate",axes=TRUE,axisnames=FALSE,
	width=0.8,space=spaces[idx1])

idx1 <- patients %in% pts2

par(mar=c(3,4,1,0))
barplot(t(allmediancov[idx1,2:1]),width=0.8,space=spaces[idx1],##col=cols,
	legend.text=c("Removing duplicates","All"),main="",
	las=2,ylab="Median depth",cex.names=0.7)
par(mar=c(1,4,3,0))
barplot(allmediancov[idx1,3],names.arg=FALSE,ylim=c(1,0),##col=cols,
	ylab="Duplication rate",axes=TRUE,axisnames=FALSE,
	width=0.8,space=spaces[idx1])

dev.off()



allmeancov <- matrix(0,nrow=nrow(samples),ncol=4)
colnames(allmeancov) <- c("All","RmDup","DupRatio","MisMatch")
rownames(allmeancov) <- samples$Label

for (i in 1:nrow(samples)) {
  sn1 <- samples$Bam[i]

  fn1 <- paste0(dir1,sn1,".bam.cov")
  fn2 <- paste0(dir1,sn1,".bam.rmdup.cov")
  fn3 <- paste0(dir1,sn1,".bam.stats")

  if (file.exists(fn1)) {
    coverage <- read.table(fn1)
    allmeancov[i,1] <- weighted.mean(coverage$V2,coverage$V3)
    coverage <- read.table(fn2)
    allmeancov[i,2] <- weighted.mean(coverage$V2,coverage$V3) 
    
    a <- read.table(fn3,sep="\t",fill=TRUE,skip=6,nrow=22)
    allmeancov[i,3] <- round(a[12,3]/a[7,3],3)
    allmeancov[i,4] <- signif(a[22,3],3)
  }
}

write.table(allmeancov,file="mCRC_All_Mean_Coverage.txt",
	    sep="\t",quote=FALSE)


tissues <- apply(do.call("rbind",strsplit(rownames(allmeancov),"_",fixed=TRUE))[,1:2],1,paste,collapse="_")
pooled_cov <- xtabs(allmeancov[,2] ~ tissues)[unique(tissues)]

patients <- sapply(strsplit(names(pooled_cov),"_",fixed=TRUE),"[",1)
pts1 <- c("V46","V402","V514","V559","V750","V824","V855","V930","V953","V974","TB1","TB7")
pts2 <- c("V402","V750","V824","V930","V953","V974","TB1","TB7")

spaces <- rep(0.25,length(pooled_cov))
spaces[which(patients[-1] != patients[-length(patients)])+1] <- 1.25

pdf("mCRC_All_Mean_Coverage.pdf",height=4,width=8)

idx1 <- patients %in% pts1
par(mar=c(5,4,2,0))
barplot(pooled_cov[idx1],width=0.8,space=spaces[idx1],##col=cols,
	main="",las=2,ylab="Mean depth",cex.names=0.8)

idx1 <- patients %in% pts2
par(mar=c(5,4,2,0))
barplot(pooled_cov[idx1],width=0.8,space=spaces[idx1],##col=cols,
	main="",las=2,ylab="Mean depth",cex.names=0.9)

dev.off()
