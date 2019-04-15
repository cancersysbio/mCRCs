#!/srv/gsfs0/software/R/3.2.2/bin/Rscript

##Generate plots of PyClone results

inputpar <- commandArgs(TRUE)
if (length(inputpar) < 1 ) stop("Wrong number of input arguments: <result folder>")

dir2 <- inputpar[1]

samples <- read.delim(file.path(dir2,"alltumorlabels.txt"),as.is=TRUE)$Bam
samplelabels <- read.delim(file.path(dir2,"alltumorlabels.txt"),as.is=TRUE)$Label

##Compare to CHAT adjusted CCF

tissues <- unique(apply(sapply(strsplit(samples,"_"),"[",1:2),2,paste,collapse="_"))
tissuelabels <- unique(apply(sapply(strsplit(samplelabels,"_"),"[",1:2),2,paste,collapse="_"))

patients <- c("V402","V750","V824","V930","V953","V974","mCRCTB1","mCRCTB7","UchiCase2","CHET9","CHET40")
patients2 <- c("V402","V750","V824","V930","V953","V974","TB1","TB7","Uchi2","Kim1","Kim2")

pdf("mCRC_All_Filtered_SNV_CHAT_vs_PyClone.pdf",width=10,height=5)

for (pt1 in patients) {
  print(pt1)
  pt2 <- patients2[match(pt1,patients)]
  fn1 <- paste0(dir2,pt1,"_MuTectSNV_Coding_CHAT_Filtered_Clean.txt")
  fn2 <- paste0(pt2,".pyclone.out")
  count1 <- read.delim(fn1,as.is=TRUE)
  if (! file.exists(fn2)) next
  count2 <- read.delim(fn2,as.is=TRUE)
  
  ##Only keep SNVs after all filters
  count1 <- count1[apply(count1[,grep("Primary_vs",names(count1)),drop=FALSE]==2,1,any),]
  count1$Pos2 <- paste(count1$Chr,count1$Pos,sep=":")
  count2$Pos2 <- sapply(strsplit(count2$X,"_"),"[",3)
  if (! ( all(count1$Pos2 %in% count2$Pos2) & nrow(count1) == nrow(count2)))
    stop("Not match")

  count2 <- count2[match(count1$Pos2,count2$Pos2),]

  sns <- sub("_freq","",grep("freq$",names(count1),value=TRUE))[-1]
  sns <- sns[sns %in% samples]
  sns2 <- samplelabels[match(sns,samples)]

  tissues1 <- grep("average_freq",names(count1),value=TRUE)
  tissues2 <- paste(pt1,sub("_average_freq","",tissues1),sep="_")
  tissuelabels1 <- tissuelabels[match(tissues2,tissues)]

  primary_idx <- which(tissues1 == "Primary_average_freq")
  other_idx <- (1:length(tissues1))[-primary_idx]

  for (j in other_idx) {
    idx1 <- c(primary_idx,j)
    clean <- count1[,paste(sapply(strsplit(tissues1[idx1],"_",fixed=0),"[",1),collapse="_vs_")] == 2
    ccf1 <- count1[clean,grep("average_CCF_adj",names(count1))[idx1]]

    idx1a <- grep(tissuelabels1[idx1[1]],sns2)
    idx1b <- grep(tissuelabels1[idx1[2]],sns2)
    cover1a <- count1[clean,paste0(sns[idx1a],"_cover")]
    cover1b <- count1[clean,paste0(sns[idx1b],"_cover")]

    ccf2 <- cbind(rowSums(count2[clean,sns2[idx1a]]*cover1a)/rowSums(cover1a),
		  rowSums(count2[clean,sns2[idx1b]]*cover1b)/rowSums(cover1b))

    par(mfrow=c(1,2),pty="s",mar=c(4,4,2,2))
    plot(ccf1,xlim=c(0,1),ylim=c(0,1),main=paste(pt2,"CHAT"),
	 xlab=paste(tissuelabels1[idx1[1]],"CCF"),
         ylab=paste(tissuelabels1[idx1[2]],"CCF"),
	 col=count2$cluster_id[clean],pch=count2$cluster_id[clean])
    plot(ccf2,xlim=c(0,1),ylim=c(0,1),main=paste(pt2,"PyClone"),
	 xlab=paste(tissuelabels1[idx1[1]],"CCF"),
         ylab=paste(tissuelabels1[idx1[2]],"CCF"),
	 col=count2$cluster_id[clean],pch=count2$cluster_id[clean])
  }
}

dev.off()
