#!/srv/gsfs0/software/R/R-3.1.1/bin/Rscript

##Different ways of handling SNVs not adjusted by CHAT and SNVs in A1 or A2 lineage

library(CHAT)

mCRC_dir <- Sys.getenv("mCRC_DIR")
dir1 <- paste0(mCRC_dir,"/scripts/data_processing/")
source(paste0(dir1,"getSampleCCF2.R"))
source(paste0(dir1,"TitanCNA_04_TitanCNA2seg.R"))

inputpar <- commandArgs(TRUE)
if (length(inputpar) < 2 ) stop("Wrong number of input arguments: <mutation file> <TitanCNA result file>")

patient1 <- strsplit(basename(inputpar[1]),"_",fixed=TRUE)[[1]][1]
load(inputpar[[2]])

chrs <- 1:22

count1 <- read.delim(inputpar[1],as.is=TRUE)
normal <-  sub("_freq","",grep("freq$",names(count1),value=TRUE))[1]
samples <- sub("_freq","",grep("freq$",names(count1),value=TRUE))[-1]
count1 <- count1[count1$Chr %in% chrs,]

for (j in 1:length(samples)) {
  sn1 <- samples[j]
  print(sn1)
  cover1 <- count1[[paste(sn1,"_cover",sep="")]]
  freq1 <- count1[[paste(sn1,"_freq",sep="")]]
  cover0 <- count1[[paste(normal,"_cover",sep="")]]
  freq0 <- count1[[paste(normal,"_freq",sep="")]]
  
  ref1 <- round(cover1 * (1-freq1))
  var1 <- round(cover1 * freq1)
  snv1 <- data.frame(Chr=count1$Chr,Pos=count1$Pos,
		     Tumor_cover=cover1,Tumor_freq=freq1,
		     Normal_cover=cover0,Normal_freq=freq0,
		     stringsAsFactors=FALSE)
  ##snv1 <- snv1[snv1$Tumor_cover > 0,]
    
  cna <- alltitancnaresults[[sn1]]
  cnaseg <- titancna2seg(cna$results,cna$convergeParams)
  
  cna1 <- cnaseg[,c(1,2,3,5,4,9,4,11,7,6)]
  cna1[,8] <- cna1[,8]*(1-cnaseg$normalproportion[1])
  cna1[,1] <- as.character(cna1[,1])
  cna1$cellularprevalence[is.na(cna1$cellularprevalence)] <- 
    max(cna1$cellularprevalence,na.rm=TRUE)
  
  ccf1 <- getSampleCCF2(snv1,cna1)
  ccf1$CCF[is.na(ccf1$CCF)] <- 2*ccf1$Tumor_freq[is.na(ccf1$CCF)]
  ccf1$CCF[ccf1$Tumor_freq==0] <- 0

  ccf1$LOH <- rep(FALSE,nrow(ccf1))
  ccf1$LOH[ccf1$n_minor == 0 & ccf1$sAGP == max(ccf1$sAGP)] <- TRUE

  purity1 <- (1 - tail(cna$convergeParams$n,1)) *
    (1 - tail(cna$convergeParams$s[1,],1))
  ccf1a <- ccf1$CCF/purity1

  count1 <- cbind(count1,round(2*freq1,3),round(2*freq1/purity1,3),FALSE)
  colnames(count1)[ncol(count1)-(2:0)] <- paste(sn1,c("_CCF","_CCF_adj","_LOH"),sep="")
  idx1 <- match(paste(ccf1$Chr,ccf1$Pos,sep="_"),paste(count1$Chr,count1$Pos,sep="_"))
  count1[idx1,ncol(count1)-2] <- round(ccf1$CCF,3)
  count1[idx1,ncol(count1)-1] <- round(ccf1a,3)
  count1[idx1,ncol(count1)] <- ccf1$LOH
  
}

##Averaging within each tissue site
tissues <- sapply(strsplit(samples,"_",fixed=TRUE),"[",2)
tissue_types <- unique(tissues)
freq <- count1[,grep("freq",names(count1))[-1]]
ccfa <- count1[,grep("CCF_adj",names(count1))]
cover <- count1[,grep("cover",names(count1))[-1]]
loh <- count1[,grep("LOH",names(count1))]
for (tissue1 in tissue_types) {
  idx1 <- which(tissues %in% tissue1)
  freq1 <- freq[,idx1,drop=FALSE]
  ccfa1 <- ccfa[,idx1,drop=FALSE]
  cover1 <- cover[,idx1,drop=FALSE]
  loh1 <- loh[,idx1,drop=FALSE]
  count1[[paste(tissue1,"_average_freq",sep="")]] <- 
    round(rowSums(freq1 * cover1) / rowSums(cover1),3)
  count1[[paste(tissue1,"_average_CCF_adj",sep="")]] <- 
    round(rowSums(ccfa1 * cover1) / rowSums(cover1),3)
  count1[[paste(tissue1,"_average_LOH",sep="")]] <- apply(loh1,1,all)
}

##truncate at 1
idx <- grep("CCF_adj",names(count1))
for (idx1 in idx)
  count1[,idx1] <- pmin(1,count1[,idx1])

##re-order columns
count1 <- count1[,c(1:8,match(paste(rep(samples,each=6),c("cover","freq","mutect","CCF","CCF_adj","LOH"),sep="_"),names(count1)),
		    match(paste(rep(tissue_types,each=3),c("average_freq","average_CCF_adj","average_LOH"),sep="_"),names(count1)),
		    match(c("Snp138","Freq1000G","COSMIC"),names(count1)))]

write.table(count1,file=sub(".txt$","_CHAT.txt",inputpar[1]),
	    quote=FALSE,sep="\t",row.names=TRUE)
