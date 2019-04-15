#!/srv/gsfs0/software/R/3.2.2/bin/Rscript

inputpar <- commandArgs(TRUE)
if (length(inputpar) < 2 ) stop("Wrong number of input arguments: <patient_name> <sample_list>")
if (!file.exists(inputpar[2]))
  stop("File ",inputpar[2]," doesn't exist")

patient <- inputpar[1]
bamfiles <- scan(inputpar[2],"c")

samplenames <- basename(bamfiles)
pileupfiles <- paste(samplenames,".pileup",sep="") 
mutectfiles <- paste(samplenames,".mutect",sep="")
mutectfilteredfiles <- paste(samplenames,".mutect.filtered.nopara.txt.varscan.filtered",sep="")
annovarfiles <- paste(samplenames,".mutect.filtered.nopara.txt.varscan.filtered.annovar",sep="") 
samplenames <- sapply(strsplit(samplenames,".",fixed=TRUE),"[",1)

n_samples <- length(samplenames)

snvdata <- vector("list",n_samples+1)
names(snvdata) <- c(paste(patient,"_Normal",sep=""),samplenames)

varscan <- read.table(paste(patient,".combined.varscanfilter",sep=""),as.is=TRUE)
varscan <- varscan[order(match(varscan$V1,c(1:22,"X","Y")),varscan$V2),]
npos <- nrow(varscan)
allsnvs <- varscan[,1:4]
names(allsnvs) <- c("Chr","Pos","Ref","Var")
allsnvs$Gene <- ""
allsnvs$Type <- ""
allsnvs_pos <- paste(allsnvs$Chr,allsnvs$Pos,sep="_")
allsnvs_alleles <- paste(allsnvs$Chr,allsnvs$Pos,allsnvs$Ref,allsnvs$Var,sep="_")
  
snvdata[[1]] <- data.frame(pcover=rep(0,npos),pfreq=rep(0,npos),
				mcover=rep(0,npos),mfreq=rep(0,npos))
for (j in 1:n_samples) {
  
  sample1 <- samplenames[j]
  bam1 <- bamfiles[j]
  mutect <- mutectfilteredfiles[j]
  mutect0 <- mutectfiles[j]
  pileup <- pileupfiles[j]

  snvdata[[j+1]] <- data.frame(pcover=rep(0,npos),pfreq=rep(0,npos),
				    mcover=rep(0,npos),mfreq=rep(0,npos),
				    mutect=rep("no",npos),
				    stringsAsFactors=FALSE)

  ##From bam files

  snv1 <- read.table(pileup,header=FALSE,as.is=TRUE,quote="",comment.char="",fill=TRUE,col.names=paste0("V",1:6))
  count1 <- data.frame(cover=rep(0,npos),read=rep(0,npos),freq=rep(0,npos))

  indels <- grep("[+-][0-9]+",snv1$V5) 
  if (length(indels) > 0) {
    for (i in indels) {
      bases <- strsplit(snv1$V5[i],"")[[1]]
      bases.indel <- gregexpr("[+-][0-9]+",snv1$V5[i])
      bases.indel.length <- abs(as.integer(regmatches(snv1$V5[i],bases.indel)[[1]]))
      for (k in length(bases.indel[[1]]):1) {
	k1 <- bases.indel[[1]][k]
	k2 <- bases.indel.length[k]
	bases <- bases[-c(k1:(k1+nchar(k2)+k2))]
      }
      snv1$V5[i] <- paste(bases,collapse="")
    }
  }

  ##Remove mapping quality scores in pileup file
  snv1$V5 <- gsub("\\^.","",snv1$V5)

  ##count alleles
  snv1$V5 <- toupper(snv1$V5)
  snv1.alleles <- matrix(0,nrow=nrow(snv1),ncol=4)
  colnames(snv1.alleles) <- c("A","C","G","T")
  for (k in 1:4) {
    snv1.alleles[,k] <- sapply(sapply(gregexpr(colnames(snv1.alleles)[k],snv1$V5),">",0),sum)
  }
  
  idx <- match(paste(allsnvs$Chr,allsnvs$Pos,allsnvs$Ref,sep="_"),
	       paste(snv1$V1,snv1$V2,snv1$V3,sep="_"))

  count1$cover[!is.na(idx)] <- snv1$V4[idx[!is.na(idx)]]
  count1$read[!is.na(idx)] <- 
    snv1.alleles[cbind(idx[!is.na(idx)],match(allsnvs$Var[!is.na(idx)],colnames(snv1.alleles)))]
  count1$freq[count1$cover > 0] <- round(count1$read[count1$cover > 0] / count1$cover[count1$cover > 0],3)

  snvdata[[j+1]]$pcover <- count1$cover
  snvdata[[j+1]]$pfreq <- count1$freq


  ##From unfiltered mutect files

  snv1 <- read.delim(mutect0,header=TRUE,comment.char="#",as.is=TRUE)
  snv1_alleles <- paste(snv1$contig,snv1$position,snv1$ref_allele,snv1$alt_allele,sep="_")
  idx <- match(snv1_alleles,allsnvs_alleles)
  idx1 <- idx[!is.na(idx)]
  idx2 <- which(!is.na(idx))
  snvdata[[j+1]]$mcover[idx1] <- snv1$t_ref_count[idx2] + snv1$t_alt_count[idx2]
  snvdata[[j+1]]$mfreq[idx1] <- round(snv1$tumor_f[idx2],3)
  snvdata[[j+1]]$mutect[idx1] <- "filtered"
  snvdata[[1]]$mcover[idx1] <- snv1$n_ref_count[idx2] + snv1$n_alt_count[idx2]
  snvdata[[1]]$mfreq[idx1] <- round(snv1$n_alt_count[idx2]/(snv1$n_ref_count[idx2] + snv1$n_alt_count[idx2]),3)

  
  ##From filtered mutect files

  snv1 <- read.delim(mutect,header=TRUE,as.is=TRUE)
  snv1_alleles <- paste(snv1$chr,snv1$pos,snv1$ref,snv1$var,sep="_")
  idx <- match(snv1_alleles,allsnvs_alleles)
  snvdata[[j+1]]$mcover[idx] <- snv1$tumor_reads
  snvdata[[j+1]]$mfreq[idx] <- round(snv1$tumor_varfrac,3)
  snvdata[[j+1]]$mutect[idx] <- "yes"
}


##Annotations

annovar <- NULL
for (i in 1:n_samples) {
  annovar1 <- read.delim(annovarfiles[i],as.is=TRUE)
  annovar <- rbind(annovar,annovar1)
}
annovar <- unique(annovar)
annovar_alleles <- paste(annovar$Chr,annovar$Start,annovar$Ref,annovar$Alt,sep="_")

idx <- match(allsnvs_alleles,annovar_alleles)
allsnvs$Gene <- annovar$Gene.refGene[idx]
allsnvs$Type <- annovar$Func.refGene[idx]
allsnvs$Type[allsnvs$Type == "exonic"] <- annovar$ExonicFunc.refGene[idx[allsnvs$Type == "exonic"]]
allsnvs$Type <- sub(" SNV","",allsnvs$Type)
allsnvs$Type[allsnvs$Type == "unknown"] <- "exonic_unknown"
allsnvs$Snp138 <- annovar$snp138[idx]
allsnvs$Freq1000G <- round(annovar$X1000g2015aug_all[idx],4)
allsnvs$Freq1000G[is.na(allsnvs$Freq1000G)] <- 0
allsnvs$COSMIC <- annovar$cosmic70[idx]

all_snv <- allsnvs[,1:6]
for (i in 1:length(snvdata)) {
  sample_snv <- snvdata[[i]]
  names(sample_snv) <- paste(names(snvdata)[i],"_",names(sample_snv),sep="")
  all_snv <- cbind(all_snv,sample_snv)
}
all_snv <- cbind(all_snv,allsnvs[,7:ncol(allsnvs)])

write.table(all_snv,paste(patient,"_MuTectSNV.txt",sep=""),
	    row.names=FALSE,quote=FALSE,sep="\t")

PON <- read.delim("PON/PON_mutect.txt",as.is=TRUE)
PON <- PON[PON$Count >= 2,]
all_snv_alleles <- paste(all_snv$Chr,all_snv$Pos,all_snv$Ref,all_snv$Var,sep="_")

##Keep non-coding SNVs
all_snv2 <- all_snv[all_snv$Freq1000G < 0.02 & varscan$V22 == "PASS" &
			 (! all_snv_alleles %in% PON$Allele),]

if (nrow(all_snv2) > 0) {
  freq <- all_snv2[,grep("mfreq",names(all_snv2))]
  cover <- all_snv2[,grep("mcover",names(all_snv2))]
  freq2 <- all_snv2[,grep("pfreq",names(all_snv2))]
  cover2 <- all_snv2[,grep("pcover",names(all_snv2))]
  freq[cover == 0] <- freq2[cover == 0]
  cover[cover == 0] <- cover2[cover == 0]
  all_snv2[,grep("mfreq",names(all_snv2))] <- freq
  all_snv2[,grep("mcover",names(all_snv2))] <- cover
}
all_snv2 <- all_snv2[,-c(grep("pfreq",names(all_snv2)),grep("pcover",names(all_snv2)))]
names(all_snv2) <- sub("mcover","cover",names(all_snv2))
names(all_snv2) <- sub("mfreq","freq",names(all_snv2))

write.table(all_snv2,paste(patient,"_MuTectSNV_All.txt",sep=""),
	    row.names=FALSE,quote=FALSE,sep="\t")


##Remove non-coding and exonic_unknown SNVs
all_snv3 <- all_snv[ annovar$Func.refGene[idx] %in% c("exonic","splicing") &
		     all_snv$Freq1000G < 0.02 & varscan$V22 == "PASS" & 
		       allsnvs$Type != "exonic_unknown" &
			 (! all_snv_alleles %in% PON$Allele),]

if (nrow(all_snv3) > 0) {
  freq <- all_snv3[,grep("mfreq",names(all_snv3))]
  cover <- all_snv3[,grep("mcover",names(all_snv3))]
  freq2 <- all_snv3[,grep("pfreq",names(all_snv3))]
  cover2 <- all_snv3[,grep("pcover",names(all_snv3))]
  freq[cover == 0] <- freq2[cover == 0]
  cover[cover == 0] <- cover2[cover == 0]
  all_snv3[,grep("mfreq",names(all_snv3))] <- freq
  all_snv3[,grep("mcover",names(all_snv3))] <- cover
}
all_snv3 <- all_snv3[,-c(grep("pfreq",names(all_snv3)),grep("pcover",names(all_snv3)))]
names(all_snv3) <- sub("mcover","cover",names(all_snv3))
names(all_snv3) <- sub("mfreq","freq",names(all_snv3))

write.table(all_snv3,paste(patient,"_MuTectSNV_Coding.txt",sep=""),
	    row.names=FALSE,quote=FALSE,sep="\t")
