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

common_ccf <- numeric(0)
private_ccf <- numeric(0)

for (i in 1:length(patients)) 
{
  pt1 <- patients[i]
  print(pt1)
  pt2 <- patients2[i]
  count1 <- read.delim(paste(dir1,pt1,"_MuTectSNV_Coding_CHAT_Filtered.txt",sep=""),as.is=TRUE)
  count1 <- count1[count1$CtoT_filter == "Pass",]

  ccf1 <- count1[,grep("average_CCF_adj",names(count1))]
  common_idx <- apply(ccf1>=0.2,1,all)

  common_ccf <- c(common_ccf,unlist(ccf1[common_idx,]))
  private_ccf <- c(private_ccf,unlist(ccf1[!common_idx,]))
}

private_ccf <- private_ccf[private_ccf>=0.1]

pdf("mCRC_CCF_Histograms.pdf")
par(mfrow=c(2,1),mar=c(4,4,3,1))
hist(common_ccf,20,xlim=c(0.1,1),main="Common SNVs",xlab="CCF")
hist(private_ccf,20,xlim=c(0.1,1),main="Other SNVs",xlab="CCF")
dev.off()
