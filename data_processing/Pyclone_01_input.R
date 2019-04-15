#!/srv/gsfs0/software/R/3.2.2/bin/Rscript

##Generate input files for PyClone

inputpar <- commandArgs(TRUE)
if (length(inputpar) < 1 ) stop("Wrong number of input parameters: <result folder>")

mCRC_dir <- Sys.getenv("mCRC_DIR")
dir1 <- paste0(mCRC_dir,"/scripts/data_processing/")
source(paste0(dir1,"getSampleCCF2.R"))
source(paste0(dir1,"TitanCNA_04_TitanCNA2seg.R"))

dir2 <- inputpar[1]
load(file.path(dir2,"mCRC_All_TitanCNA.RData"))

samples <- read.delim(file.path(dir2,"alltumorlabels.txt"),as.is=TRUE)$Bam
samplelabels <- read.delim(file.patch(dir2,"alltumorlabels.txt",as.is=TRUE)$Label

patients <- unique(sapply(strsplit(samples,"_",fixed=TRUE),"[",1))
patients2 <- unique(sapply(strsplit(samplelabels,"_",fixed=TRUE),"[",1))

all_purity <- data.frame(Sample=samplelabels,Purity=rep(0,length(samplelabels)))

for (pt1 in patients) {
  print(pt1)
  fn1 <- paste(dir2,pt1,"_MuTectSNV_Coding_CHAT_Filtered_Clean.txt",sep="")
  count1 <- read.delim(fn1,as.is=TRUE)
  
  ##Only keep SNVs after all filters
  count1 <- count1[apply(count1[,grep("Primary_vs",names(count1)),drop=FALSE]==2,1,any),]

  normal <- sub("_freq","",grep("freq$",names(count1),value=TRUE))[1]
  sns <- sub("_freq","",grep("freq$",names(count1),value=TRUE))[-1]
  sns <- sns[sns %in% samples]
  sns2 <- samplelabels[match(sns,samples)]

  for (j in 1:length(sns)) {
    sn1 <- sns[j]
    sn2 <- sns2[j]
    print(sn1)

    cna1 <- alltitancnaresults[[sn1]]
    purity1 <- round((1 - tail(cna1$convergeParams$n,1)) *
		       (1 - tail(cna1$convergeParams$s[1,],1)),3)
    all_purity$Purity[match(sn2,all_purity$Sample)] <- purity1

    cnaseg <- titancna2seg(cna1$results,cna1$convergeParams)
    cna2 <- cnaseg[,c(1,2,3,5,4,9,4,11,7,6)]
    cna2[,8] <- cna2[,8]*(1-cnaseg$normalproportion[1])
    cna2[,1] <- as.character(cna2[,1])
    cna2$cellularprevalence[is.na(cna2$cellularprevalence)] <- 
      max(cna2$cellularprevalence,na.rm=TRUE)

    cover1 <- count1[[paste(sn1,"_cover",sep="")]]
    cover1[cover1 == 0] <- 1 ##dummy count
    freq1 <- count1[[paste(sn1,"_freq",sep="")]]
    cover0 <- count1[[paste(normal,"_cover",sep="")]]
    freq0 <- count1[[paste(normal,"_freq",sep="")]]
    
    ref1 <- round(cover1 * (1-freq1))
    var1 <- round(cover1 * freq1)
    snv1 <- data.frame(Chr=count1$Chr,Pos=count1$Pos,
		       Tumor_cover=cover1,Tumor_freq=freq1,
		       Normal_cover=cover0,Normal_freq=freq0,
		       stringsAsFactors=FALSE)
 
    ccf1 <- getSampleCCF2(snv1,cna2)
   
    cover1 <- count1[[paste(sn1,"_cover",sep="")]]
    freq1 <- count1[[paste(sn1,"_freq",sep="")]]
    ref1 <- round(cover1 * (1-freq1))
    var1 <- round(cover1 * freq1)
    snv1 <- data.frame(mutation_id=paste(count1$Gene,"_",count1$Ref,">",count1$Var,"_",
		       count1$Chr,":",count1$Pos,sep=""),
		       ref_counts=ref1,var_counts=var1,
		       normal_cn=2,
		       minor_cn=ccf1$n_minor,major_cn=ccf1$n_total-ccf1$n_minor)
    
    ##if (cna1 == 0) cna1 <- 1 ##temp fix
    
    write.table(snv1,file=paste(sn2,"_SNV.mut",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
  }
}

write.table(all_purity,file="mCRC_All_Purity.txt",
	    quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
