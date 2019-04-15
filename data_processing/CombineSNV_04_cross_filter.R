#!/srv/gsfs0/software/R/3.2.2/bin/Rscript

##Cross sample filter for SNVs observed with low AF in multiple cases

fns <- list.files(".","MuTectSNV_All.txt$")
mutect_fns <- list.files(".","mutect$")

pts <- sapply(strsplit(fns,"_",fixed=TRUE),"[",1)

cross_filter <- NULL
for (i in 1:length(pts)) {
  pt1 <- pts[i]
  fn1 <- fns[i]
  snv1 <- read.delim(fn1,as.is=TRUE)
  cross_filter <- rbind(cross_filter,cbind(Patient=pt1,snv1[,1:4]))
}

all_alleles <- paste(cross_filter[,2],cross_filter[,3],cross_filter[,4],cross_filter[,5],sep="_")

for(j in 1:length(mutect_fns)) {
  fn2 <- mutect_fns[j]
  pt2 <- strsplit(fn2,"_",fixed=TRUE)[[1]][1]
  print(fn2)
  mutect2 <- read.delim(fn2,comment.char="#",as.is=TRUE)
  mutect2_alleles <- paste(mutect2[,1],mutect2[,2],mutect2[,4],mutect2[,5],sep="_")
  mutect2 <- mutect2[mutect2_alleles %in% all_alleles,]
  mutect2_alleles <- paste(mutect2[,1],mutect2[,2],mutect2[,4],mutect2[,5],sep="_")
  if (is.null(cross_filter[[pt2]])) 
    cross_filter[[pt2]] <- 0
  idx1 <- which(all_alleles %in% mutect2_alleles)
  idx2 <- match(all_alleles[idx1], mutect2_alleles)
  cross_filter[[pt2]][idx1] <- pmax(cross_filter[[pt2]][idx1],mutect2$tumor_f[idx2])
  write.table(cross_filter,file="mCRC_Vienna_Cross_Filter.txt",
	      quote=FALSE,sep="\t",row.names=FALSE)
}

cross_filter$Count <- 0
for (i in 1:length(pts)) {
  idx1 <- cross_filter$Patient == pts[i]
  cross_filter$Count[idx1] <- rowSums(cross_filter[idx1,pts[-i]] >= 0.01 & cross_filter[idx1,pts[-i]] <= 0.05)
}

write.table(cross_filter,file="mCRC_Vienna_Cross_Filter.txt",
	    quote=FALSE,sep="\t",row.names=FALSE)

for (i in 1:length(pts)) {
  pt1 <- pts[i]
  fns1 <- list.files(".",paste0(pt1,"_MuTectSNV_Coding.txt"))
  fns1 <- grep("All",fns1,invert=TRUE,value=TRUE)
  filter1 <- cross_filter[cross_filter$Patient == pt1,]
  for (fn1 in fns1) {
    snv1 <- read.delim(fn1,as.is=TRUE)
    idx1 <- 
      match(paste(snv1$Chr,snv1$Pos,snv1$Ref,snv1$Var,sep="_"),
	    paste(filter1$Chr,filter1$Pos,filter1$Ref,filter1$Var,sep="_"))
    idx2 <- which(filter1$Count[idx1] >= 2)
    snv1 <- snv1[-idx2,]
    write.table(snv1,fn1,quote=FALSE,sep="\t",row.names=FALSE)
  }
}
