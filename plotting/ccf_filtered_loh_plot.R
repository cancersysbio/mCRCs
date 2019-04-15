##Remove SNVs in potential unique LOH regions but not shrink tumor site CCFs to 0 as done in the "Clean" steps

mCRC_dir <- Sys.getenv("mCRC_DIR")
setwd(file.path(mCRC_dir,"plot"))

dir1 <- "../results/"
samples <- read.delim("../results/alltumorlabels.txt",as.is=TRUE)$Bam
samplelabels <- read.delim("../results/alltumorlabels.txt",as.is=TRUE)$Label

tissues <- unique(apply(sapply(strsplit(samples,"_"),"[",1:2),2,paste,collapse="_"))
tissuelabels <- unique(apply(sapply(strsplit(samplelabels,"_"),"[",1:2),2,paste,collapse="_"))

##Cases with primary samples and multi-region sampling
patients <- c("V402","V750","V824","V930","V953","V974","mCRCTB1","mCRCTB7","UchiCase2","CHET9","CHET40")
patients2 <- c("V402","V750","V824","V930","V953","V974","TB1","TB7","Uchi2","Kim1","Kim2")

short_list <- scan("../results/CRCGenes.txt","c")

label_genes <- short_list

my.cols <- colorRampPalette(c("white", "yellow", "darkred"), space="rgb")

pdf("mCRC_All_Filtered_LOH_SNV_CCF.pdf",width=7.5,height=12.5)

layout(matrix(c(1:11,0,12:14),ncol=3,byrow=TRUE))
par(mar=c(3,3,3,3),mgp=c(2,1,0),pty="s")

for (i in 1:length(patients)) 
{
  pt1 <- patients[i]
  print(pt1)
  pt2 <- patients2[i]
  count1 <- read.delim(paste(dir1,pt1,"_MuTectSNV_Coding_CHAT_Filtered.txt",sep=""),as.is=TRUE)
  count2 <- read.delim(paste(dir1,pt1,"_MuTectSNV_Coding_CHAT_Filtered_Clean.txt",sep=""),as.is=TRUE)
  count1 <- count1[count1$CtoT_filter == "Pass",]
  count2 <- count2[count2$CtoT_filter == "Pass",]
  if (! all(count1$Pos == count2$Pos)) stop("Not match")

  tissues1 <- grep("average_freq",names(count1),value=TRUE)
  tissues2 <- paste(pt1,sub("_average_freq","",tissues1),sep="_")
  tissuelabels1 <- tissuelabels[match(tissues2,tissues)]

  primary_idx <- which(tissues1 == "Primary_average_freq")
  other_idx <- (1:length(tissues1))[-primary_idx]

  for (j in other_idx) {
    idx2 <- c(primary_idx,j)
    clean <- count2[,paste(sapply(strsplit(tissues1[idx2],"_",fixed=0),"[",1),collapse="_vs_")]
    count3 <- count1[clean==2,]
    
    freq1 <- count3[,grep("average_freq",names(count3))[idx2]]
    ccf1 <- count3[,grep("average_CCF_adj",names(count3))[idx2]]

    nonzero <- apply(freq1 > 0,1,any)
    gene_idx <- which(count3$Type != "synonymous" & 
			count3$Gene %in% label_genes & 
			  apply(freq1 >= 0.05,1,any))
    cols <- rep(1,nrow(count3))
    cols[gene_idx] <- 2
    pchs <- rep(1,nrow(count3))
    pchs[gene_idx] <- 2
 
    smoothScatter(freq1[nonzero,],xlim=c(-0.05,1),ylim=c(-0.05,1),
		  pch=1,cex=0.3,colramp=my.cols,nrpoints=0,
		  xlab=paste(tissuelabels1[idx2[1]],"VAF"),
		  ylab=paste(tissuelabels1[idx2[2]],"VAF"),
		  main=pt2)
    points(freq1[gene_idx,],cex=0.6,col=2,pch=2)
    text(freq1[gene_idx,],cex=0.7,pos=1,count3$Gene[gene_idx])
  }
}

for (i in 1:length(patients)) 
{
  pt1 <- patients[i]
  print(pt1)
  pt2 <- patients2[i]
  count1 <- read.delim(paste(dir1,pt1,"_MuTectSNV_Coding_CHAT_Filtered.txt",sep=""),as.is=TRUE)
  count2 <- read.delim(paste(dir1,pt1,"_MuTectSNV_Coding_CHAT_Filtered_Clean.txt",sep=""),as.is=TRUE)
  count1 <- count1[count1$CtoT_filter == "Pass",]
  count2 <- count2[count2$CtoT_filter == "Pass",]
  if (! all(count1$Pos == count2$Pos)) stop("Not match")

  tissues1 <- grep("average_freq",names(count1),value=TRUE)
  tissues2 <- paste(pt1,sub("_average_freq","",tissues1),sep="_")
  tissuelabels1 <- tissuelabels[match(tissues2,tissues)]

  primary_idx <- which(tissues1 == "Primary_average_freq")
  other_idx <- (1:length(tissues1))[-primary_idx]

  for (j in other_idx) {
    idx2 <- c(primary_idx,j)
    clean <- count2[,paste(sapply(strsplit(tissues1[idx2],"_",fixed=0),"[",1),collapse="_vs_")]
    count3 <- count1[clean==2,]

    freq1 <- count3[,grep("average_freq",names(count3))[idx2]]
    ccf1 <- count3[,grep("average_CCF_adj",names(count3))[idx2]]

    nonzero <- apply(freq1 > 0,1,any)
    gene_idx <- which(count3$Type != "synonymous" & 
			count3$Gene %in% label_genes & 
			  apply(freq1 >= 0.05,1,any))
    cols <- rep(1,nrow(count3))
    cols[gene_idx] <- 2
    pchs <- rep(1,nrow(count3))
    pchs[gene_idx] <- 2
    
    plot(freq1[nonzero,],xlim=c(-0.05,1),ylim=c(-0.05,1),cex=0.6,
	 col=cols[nonzero],pch=pchs[nonzero],
	 xlab=paste(tissuelabels1[idx2[1]],"VAF"),
	 ylab=paste(tissuelabels1[idx2[2]],"VAF"),
	 main=pt2)
    text(freq1[gene_idx,],cex=0.7,pos=1,count3$Gene[gene_idx])
  }
}

for (i in 1:length(patients)) 
{
  pt1 <- patients[i]
  print(pt1)
  pt2 <- patients2[i]
  count1 <- read.delim(paste(dir1,pt1,"_MuTectSNV_Coding_CHAT_Filtered.txt",sep=""),as.is=TRUE)
  count2 <- read.delim(paste(dir1,pt1,"_MuTectSNV_Coding_CHAT_Filtered_Clean.txt",sep=""),as.is=TRUE)
  count1 <- count1[count1$CtoT_filter == "Pass",]
  count2 <- count2[count2$CtoT_filter == "Pass",]
  if (! all(count1$Pos == count2$Pos)) stop("Not match")

  tissues1 <- grep("average_freq",names(count1),value=TRUE)
  tissues2 <- paste(pt1,sub("_average_freq","",tissues1),sep="_")
  tissuelabels1 <- tissuelabels[match(tissues2,tissues)]

  primary_idx <- which(tissues1 == "Primary_average_freq")
  other_idx <- (1:length(tissues1))[-primary_idx]

  for (j in other_idx) {
    idx2 <- c(primary_idx,j)
    clean <- count2[,paste(sapply(strsplit(tissues1[idx2],"_",fixed=0),"[",1),collapse="_vs_")]
    count3 <- count1[clean==2,]

    freq1 <- count3[,grep("average_freq",names(count3))[idx2]]
    ccf1 <- count3[,grep("average_CCF_adj",names(count3))[idx2]]

    nonzero <- apply(freq1 > 0,1,any)
    gene_idx <- which(count3$Type != "synonymous" & 
			count3$Gene %in% label_genes & 
			  apply(freq1 >= 0.05,1,any))
    cols <- rep(1,nrow(count3))
    cols[gene_idx] <- 2
    pchs <- rep(1,nrow(count3))
    pchs[gene_idx] <- 2
    
    top_cluster <- apply(ccf1 >= 0.8,1,all)
    gene_idx2 <-  which(count3$Type != "synonymous" & 
			  count3$Gene %in% label_genes & top_cluster)
    
    smoothScatter(ccf1[nonzero,],xlim=c(-0.05,1),ylim=c(-0.05,1),
		  pch=1,cex=0.3,colramp=my.cols,nrpoints=0,
		  xlab=paste(tissuelabels1[idx2[1]],"CCF"),
		  ylab=paste(tissuelabels1[idx2[2]],"CCF"),
		  main=pt2)
    points(ccf1[gene_idx,],cex=0.6,col=2,pch=2)
    text(ccf1[gene_idx,],cex=0.7,pos=1,count3$Gene[gene_idx])
    gene_idx3 <- gene_idx2[order(ccf1[gene_idx2,2],decreasing=TRUE)]
    mtext(paste(count3$Gene[gene_idx3],collapse="\n"),4,line=0.2,
	  las=2,at=1,padj=1,cex=0.5)

    gene_idx4 <-  which(count3$Type != "synonymous" & 
			  ccf1[,1] < 0.2 & ccf1[,2] >= 0.6)
    gene_idx4 <- gene_idx4[order(ccf1[gene_idx4,2],decreasing=TRUE)]
    ##mtext(paste(count3$Gene[gene_idx4],collapse="\n"),2,line=0.2,
    ##	  las=2,at=1,padj=1,cex=0.5)
  }
}
  
for (i in 1:length(patients)) 
{
  pt1 <- patients[i]
  print(pt1)
  pt2 <- patients2[i]
  count1 <- read.delim(paste(dir1,pt1,"_MuTectSNV_Coding_CHAT_Filtered.txt",sep=""),as.is=TRUE)
  count2 <- read.delim(paste(dir1,pt1,"_MuTectSNV_Coding_CHAT_Filtered_Clean.txt",sep=""),as.is=TRUE)
  count1 <- count1[count1$CtoT_filter == "Pass",]
  count2 <- count2[count2$CtoT_filter == "Pass",]
  if (! all(count1$Pos == count2$Pos)) stop("Not match")

  tissues1 <- grep("average_freq",names(count1),value=TRUE)
  tissues2 <- paste(pt1,sub("_average_freq","",tissues1),sep="_")
  tissuelabels1 <- tissuelabels[match(tissues2,tissues)]

  primary_idx <- which(tissues1 == "Primary_average_freq")
  other_idx <- (1:length(tissues1))[-primary_idx]

  for (j in other_idx) {
    idx2 <- c(primary_idx,j)
    clean <- count2[,paste(sapply(strsplit(tissues1[idx2],"_",fixed=0),"[",1),collapse="_vs_")]
    count3 <- count1[clean==2,]

    freq1 <- count3[,grep("average_freq",names(count3))[idx2]]
    ccf1 <- count3[,grep("average_CCF_adj",names(count3))[idx2]]

    nonzero <- apply(freq1 > 0,1,any)
    gene_idx <- which(count3$Type != "synonymous" & 
			count3$Gene %in% label_genes & 
			  apply(freq1 >= 0.05,1,any))
    cols <- rep(1,nrow(count3))
    cols[gene_idx] <- 2
    pchs <- rep(1,nrow(count3))
    pchs[gene_idx] <- 2
    
    top_cluster <- apply(ccf1 >= 0.8,1,all)
    gene_idx2 <-  which(count3$Type != "synonymous" & 
			  count3$Gene %in% label_genes & top_cluster)
    plot(ccf1[nonzero,],xlim=c(-0.05,1),ylim=c(-0.05,1),cex=0.6,
	 col=cols[nonzero],pch=pchs[nonzero],
	 xlab=paste(tissuelabels1[idx2[1]],"CCF"),
	 ylab=paste(tissuelabels1[idx2[2]],"CCF"),
	 main=pt2)
    text(ccf1[gene_idx,],cex=0.7,pos=1,count3$Gene[gene_idx])
    gene_idx3 <- gene_idx2[order(ccf1[gene_idx2,2],decreasing=TRUE)]
    mtext(paste(count3$Gene[gene_idx3],collapse="\n"),4,line=0.2,
	  las=2,at=1,padj=1,cex=0.5)

    gene_idx4 <-  which(count3$Type != "synonymous" & 
			  ccf1[,1] < 0.2 & ccf1[,2] >= 0.6)
    gene_idx4 <- gene_idx4[order(ccf1[gene_idx4,2],decreasing=TRUE)]
    ##mtext(paste(count3$Gene[gene_idx4],collapse="\n"),2,line=0.2,
    ##	  las=2,at=1,padj=1,cex=0.5)
   }
}

dev.off()
