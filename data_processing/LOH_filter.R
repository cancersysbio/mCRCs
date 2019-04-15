#!/srv/gsfs0/software/R/3.2.2/bin/Rscript

##Generate "clean" lists of mutations used for computational analyses

samples <- read.delim("alltumorlabels.txt",as.is=TRUE)$Bam
samplelabels <- read.delim("alltumorlabels.txt",as.is=TRUE)$Label

patients <- unique(sapply(strsplit(samples,"_",fixed=TRUE),"[",1))
patients2 <- unique(sapply(strsplit(samplelabels,"_",fixed=TRUE),"[",1))

cutoff <- 0.6

pdf("mCRC_SNV_LOH_Clean.pdf")

for (pt1 in patients) {
  print(pt1)
  fn1 <- paste(pt1,"_MuTectSNV_Coding_CHAT_Filtered.txt",sep="")
  count1 <- read.delim(fn1,as.is=TRUE)
  ##count1 <- count1[count1$CtoT_filter == "Pass",]

  sns <- sub("_mutect","",grep("mutect",names(count1),value=TRUE))
  sites <- sapply(strsplit(sns,"_",fixed=TRUE),"[",2)
  site_types <- unique(sites)
  
  mutect <- count1[,paste(sns,"_mutect",sep="")]
  freqs <- count1[,paste(sns,"_freq",sep="")]
  cover <- count1[,paste(sns,"_cover",sep="")]
  ccf <- count1[,paste(sns,"_CCF_adj",sep="")]
  loh <- count1[,paste(sns,"_LOH",sep="")]
  reads <- round(cover*freqs)
  colnames(mutect) <- colnames(freqs) <- colnames(loh) <- 
    colnames(cover) <- colnames(reads) <- colnames(ccf) <- sns
  
  present <- (reads >= 3 & freqs >= 0.05) | freqs >= 0.1
  absent <- reads <= 2 & freqs < 0.05 & cover >= 10
  ##"Present" in at least one sample 
  ##idx1 <- !apply(present,1,any)

  ##Set average CCF to 0 for mutations with low variant read counts and VAFs in all samples from one tumor site
  present1 <- (reads >= 3 & freqs >= 0.03) | freqs >= 0.1 ##Lower cutoff for VAF
  for (site1 in site_types) {
    idx_1 <- sites == site1
    idx_np <- apply(!present1[,idx_1,drop=FALSE],1,all)
    if (any(idx_np))
      count1[idx_np,paste0(site1,"_average_CCF_adj")] <- 0
  }

  site1 <- "Primary"
  for (site2 in site_types[site_types != "Primary"]) {
    idx <- sites %in% c(site1,site2)
    idx_1 <- sites == site1
    idx_2 <- sites == site2

    ##"Present" in at least one sample from this pair of tumor sites
    idx1 <- !apply(present[,idx],1,any)
    
    ##Mutations absent in subset of samples potentially due to LOH
    idx2 <- (! apply(present[,idx] & loh[,idx],1,any)) & 
      ((apply(absent[,idx_1,drop=FALSE] & loh[,idx_1,drop=FALSE],1,any) & 
	  ( ! (apply(absent[,idx_1,drop=FALSE] & ! loh[,idx_1,drop=FALSE],1,any))) & 
	    ! (apply(loh[,idx_1,drop=FALSE],1,all) & count1[,paste0(site2,"_average_CCF_adj")] < cutoff)) |
	      (apply(absent[,idx_2,drop=FALSE] & loh[,idx_2,drop=FALSE],1,any) &
		 ( ! (apply(absent[,idx_2,drop=FALSE] & ! loh[,idx_2,drop=FALSE],1,any))) &
		   ! (apply(loh[,idx_2,drop=FALSE],1,all) & count1[,paste0(site1,"_average_CCF_adj")] < cutoff))) 

    ##Remove all SNVs with clonal LOH in one site and not clonal LOH in the other site
    idx2b <- rep(FALSE,length(idx2))
    idx2b[(apply(loh[,idx_1,drop=FALSE],1,all) & apply(!loh[,idx_2,drop=FALSE],1,any)) | 
	  (apply(loh[,idx_2,drop=FALSE],1,all) & apply(!loh[,idx_1,drop=FALSE],1,any))] <- TRUE


    ##Mutations with less than 20 coverage in one site
    cover2 <- apply(cover[,idx],1,function(x) xtabs(x~sites[idx]))
    idx3 <- apply(cover2<20,2,any) ##| apply(cover[,idx]<5,1,any)
    
    keep <- rep(TRUE,nrow(count1))
    keep[count1$CtoT_filter != "Pass"] <- FALSE
    keep[idx1 | idx2 | idx3] <- FALSE
    keep2 <- keep
    keep2[idx2b] <- FALSE

    count1[[paste0(site1,"_vs_",site2)]] <- 0
    count1[[paste0(site1,"_vs_",site2)]][keep] <- 1
    count1[[paste0(site1,"_vs_",site2)]][keep2] <- 2

  }
  
  write.table(count1,file=sub(".txt","_Clean.txt",fn1),
  	      sep="\t",quote=FALSE,row.names=FALSE)
}

dev.off()
