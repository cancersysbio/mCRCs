#!/srv/gsfs0/software/R/3.2.2/bin/Rscript

##Read and tabulate filter results

##CtoT individual FFPE samples

filter <- read.delim("mCRC_Vienna_CtoT_1.pass.all.maf.annotated",
		     as.is=TRUE,comment.char="#")
filter[,2] <- as.character(filter[,2])		      
snvs <- apply(filter[,c(1,2,4,5,8)],1,paste,collapse="_")

filter2 <- unique(filter[,c(1,2,4,5,8)])
snvs2 <- apply(filter2,1,paste,collapse="_")
filter2$Patient <- sapply(strsplit(filter2$Matched_Norm_Sample_Barcode,"_",fixed=TRUE),"[",1)

filter2$CtoT_1 <- 1
filter2$CtoT_1[snvs2 %in% snvs[filter$oxoGCut == 0]] <- 0

fail_snvs <- apply(filter2[filter2$CtoT_1==1,c(6,1,2,3,4)],1,paste,collapse="_")

##Add to original mutation tables

fns1 <- c(list.files(".","Coding.txt"),
	  list.files(".","CHAT.txt"))

for (fn1 in fns1) {
  print(fn1)
  pt1 <- strsplit(basename(fn1),"_")[[1]][1]
  snv1 <- read.delim(fn1,as.is=TRUE)
  snv1$CtoT_filter <- rep("Pass",nrow(snv1))
  snvs1 <- paste(pt1,snv1[,1],snv1[,2],snv1[,3],snv1[,4],sep="_")
  snv1$CtoT_filter[snvs1 %in% fail_snvs] <- "Fail"
  write.table(snv1,file=sub(".txt","_Filtered.txt",fn1),quote=FALSE,sep="\t",row.names=FALSE)
}

pdf("mCRC_Vienna_FFPE_Artifact_Filters.pdf",width=8,height=4)

fns1 <- list.files(".","MuTectSNV_Coding_Filtered.txt",full.names=TRUE)
fns2 <- list.files(".","MuTectSNV_Coding_CHAT_Filtered.txt",full.names=TRUE)

for (i in 1:length(fns1)) {
    snv1 <- read.delim(fns1[i],as.is=TRUE)
    snv2 <- read.delim(fns2[i],as.is=TRUE)
    pt1 <- strsplit(basename(fns2[i]),"_")[[1]][1]

    idx2 <- snv2$CtoT == "Fail"
    idx1 <- !idx2

    tissues <- sub("_average_freq","",grep("_average_freq",names(snv2),value=TRUE))
    if (! "Primary" %in% tissues) next
    tissues <- tissues[tissues != "Primary"]

    freq1 <- snv2[,"Primary_average_freq"]
    ccf1<- snv2[,"Primary_average_CCF_adj"]
    for (tissue1 in tissues) {
    	freq2 <- snv2[,paste0(tissue1,"_average_freq")]
	ccf2 <- snv2[,paste0(tissue1,"_average_CCF_adj")]

	par(mfrow=c(1,2),mar=c(4,4,2,2),pty="s")

	plot(freq1[idx1],freq2[idx1],xlim=c(0,1),ylim=c(0,1),cex=0.5,
	     xlab=paste(pt1,"Primary"),ylab=paste(pt1,tissue1),main="VAF")
	points(freq1[idx2],freq2[idx2],cex=0.6,col=2,pch=2)
	legend("topright",
	       paste(c("Not filtered","Filtered C>T"),c(sum(idx1),sum(idx2))),
	       col=c(1,2,4),pch=c(1,2,3),
	       pt.cex=c(0.5,0.6,0.6),cex=0.8)
	
	plot(ccf1[idx1],ccf2[idx1],xlim=c(0,1),ylim=c(0,1),cex=0.5,
	     xlab=paste(pt1,"Primary"),ylab=paste(pt1,tissue1),main="CCF")
	points(ccf1[idx2],ccf2[idx2],cex=0.6,col=2,pch=2)
    }
}

dev.off()
