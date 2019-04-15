#!/srv/gsfs0/software/R/3.2.2/bin/Rscript

##Generate input files for LicheE

inputpar <- commandArgs(TRUE)
if (length(inputpar) < 1 ) stop("Wrong number of input arguments: <result folder>")

dir1 <- inputpar[1]

samples <- read.delim(file.path(dir1,"alltumorlabels.txt"),as.is=TRUE)$Bam
samplelabels <- read.delim(file.path(dir1,"alltumorlabels.txt"),as.is=TRUE)$Label

##Cases with primary samples
patients <- unique(sapply(strsplit(samples,"_",fixed=TRUE),"[",1))

for (pt1 in patients) {
  fn1 <- paste(dir1,"/",pt1,"_MuTectSNV_Coding_CHAT_Filtered_Clean.txt",sep="")
  count1 <- read.delim(fn1,as.is=TRUE)
  ##SNVs passing all filters for all pairwise primary-met comparions
  count1 <- count1[apply(count1[,grep("Primary_vs",names(count1)),drop=FALSE]==2,1,all),]
  sns <- sub("_mutect","",grep("mutect",names(count1),value=TRUE))
  sns <- sns[sns %in% samples]
  sns2 <- samplelabels[match(sns,samples)]
  
  count1_allele <- paste(count1$Gene,count1$Chr,count1$Pos,count1$Ref,count1$Var,sep="_")
  mutect <- count1[,paste(sns,"_mutect",sep="")]
  freqs <- count1[,paste(sns,"_freq",sep="")]
  cover <- count1[,paste(sns,"_cover",sep="")]
  ccf <- count1[,paste(sns,"_CCF_adj",sep="")]
  loh <- count1[,paste(sns,"_LOH",sep="")]
  reads <- round(cover*freqs)
  colnames(mutect) <- colnames(freqs) <- 
    colnames(cover) <- colnames(reads) <- colnames(ccf) <- sns2
  
  freqs2 <- ccf/2
  present_code <- apply(0+cbind(FALSE,present),1,paste,collapse="")
  snv_name <- count1_allele[idx1]
  
  count2 <- cbind(count1[idx1,1:2],Name=snv_name,##Porfile=present_code,
		  Normal=0.0,freqs2)
  write.table(count2,file=paste(pt1,"_SNV_Lichee_Input.txt",sep=""),
	      sep="\t",quote=FALSE,row.names=FALSE)
}
