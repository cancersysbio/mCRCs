#!/srv/gsfs0/software/R/3.2.2/bin/Rscript

##Salvage read counts from all samples in each case 

library(Rsamtools)
options(stringsAsFactors=FALSE)

inputpar <- commandArgs(TRUE)
if (length(inputpar) < 2 ) stop("Wrong number of input arguments: <patient_name> <sample_list>")
if (!file.exists(inputpar[2]))
  stop("File ",inputpar[2]," doesn't exist")

patient <- inputpar[1]
bamfiles <- scan(inputpar[2],"c")
samples <- sapply(strsplit(basename(bamfiles),".",fixed=TRUE),"[",1)
print(patient)

indels <- NULL

for (sn1 in samples) {
  fn1 <- paste0(sn1,".strelka.indels.annovar")
  indel1 <- read.delim(fn1,as.is=TRUE)
  indel1 <- indel1[indel1$Func.refGene %in% c("exonic","splicing"),]
  indel1$Gene <- indel1$Gene.refGene
  indel1$Type <- indel1$Func.refGene
  indel1$Type[indel1$Type=="exonic"] <- indel1$ExonicFunc.refGene[indel1$Type=="exonic"]
  indel1$Type <- sub(" deletion","",indel1$Type)
  indel1$Type <- sub(" insertion","",indel1$Type)
  indel1$Type <- paste0(indel1$Type,"_indel")
  indel1[[paste0(patient,"_Normal_cover")]] <- indel1$Normal_cover
  indel1[[paste0(patient,"_Normal_freq")]] <- round(indel1$Normal_alt_freq,3)
  indel1[[paste0(sn1,"_cover")]] <- indel1$Cover
  indel1[[paste0(sn1,"_freq")]] <- round(indel1$Alt_freq,3)
  indel1[[paste0(sn1,"_strelka")]] <- "yes"

  indel1 <- indel1[,c(1:5,20:26,17:19)]
  names(indel1)[5] <- "Var"
  names(indel1)[13:15] <- c("Snp138","Freq1000G","COSMIC")

  if (is.null(indels)) 
    indels <- indel1
  else
    indels <- merge(indels,indel1,all=TRUE)
}

indels <- indels[,c(1:9,13:ncol(indels),10:12)]
indels <- indels[order(match(indels$Chr,c(1:22,"X","Y")),indels$Start,indels$End),]

covers <- indels[,grep("cover",names(indels))[-1]]
freqs <- indels[,grep("freq",names(indels))[-1]]
reads <- round(covers * freqs)
indels <- indels[apply(covers >= 10 & reads >= 3,1,any,na.rm=TRUE),]
rownames(indels) <- NULL

##Salvage reads

indel_types <- rep("",nrow(indels))
indel_types[indels$Ref == "-"] <- "+"
indel_types[indels$Var == "-"] <- "-"
pos <- GRanges(seqnames=indels$Chr,
	       ranges=IRanges(indels$Start,indels$Start))
pos2 <- paste(indels$Chr,indels$Start,sep="_")

for (i in 1:length(samples)) {
  sn1 <- samples[i]
  bam1 <- bamfiles[i]
  print(sn1)

  indels[is.na(indels[,paste0(sn1,"_strelka")]),paste0(sn1,"_strelka")] <- "no"

  if (! file.exists(bam1))
    stop(paste(bam1,"Bam file doesn't exists"))
  
  pileup1 <- pileup(bam1,scanBamParam=ScanBamParam(which=pos),
		    pileupParam=PileupParam(max_depth=2000,
		      min_base_quality=20, min_mapq=40,
		      distinguish_strands=FALSE,
		      include_deletions=TRUE,include_insertions=TRUE))
  pileup1$pos2 <- paste(pileup1$seqnames,pileup1$pos,sep="_")
  
  idx1 <- which(is.na(indels[,paste0(sn1,"_freq")]))
  for (k in idx1) {
    indels[k,paste0(sn1,"_cover")] <- 
      sum(pileup1$count[pileup1$pos2 == pos2[k]])
    indels[k,paste0(sn1,"_freq")] <-
      round(sum(pileup1$count[pileup1$pos2 == pos2[k] & pileup1$nucleotide == indel_types[k]]) / 
	      sum(pileup1$count[pileup1$pos2 == pos2[k]]),3)
    if (indels[k,paste0(sn1,"_cover")] == 0)
      indels[k,paste0(sn1,"_freq")] <- 0
  }
}

if (any(indels$Freq1000G >= 0.02, na.rm=TRUE))
  indels <- indels[-which(indels$Freq1000G >= 0.02),] ##1000 Genomes VAF

indels <- indels[ ! indels$Gene == "SIRPB1",] ##artifacts due to different sequencing platforms for normal and tumor samples
indels <- indels[!grepl("unknown",indels$Type),]

write.table(indels,file=paste0(patient,"_Strelka_Indel_Coding.txt"),
	    quote=FALSE,sep="\t",row.names=FALSE)
