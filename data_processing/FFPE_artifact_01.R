#!/srv/gsfs0/software/R/3.2.2/bin/Rscript

##Generate input files for D-ToxoG filter
##_1.txt: individual samples

inputpar <- commandArgs(TRUE)
if (length(inputpar) < 1 ) stop("Wrong number of input arguments: <bam folder>")

library(Rsamtools)

bamdir <- inputpar[1]

fns <- list.files(".","_MuTectSNV_Coding.txt")
pts <- sub("_MuTectSNV_Coding.txt","",fns)

all_snvs1 <- all_snvs2 <- NULL

for (i in 1:length(pts)) {
  print(pts[i])
  snv_file <- fns[i]
  snv <- read.delim(snv_file,as.is=TRUE)
  normal <- sub("_cover","",grep("cover",names(snv),value=TRUE)[1])
  sns <- sub("_mutect","",grep("mutect",names(snv),value=TRUE))
  mutect_call <- snv[,grep("_mutect",names(snv))]
  
  bams <- paste0(sns,".bam")
  pos <- GRanges(seqnames=snv$Chr,
		 ranges=IRanges(snv$Pos,snv$Pos))
  pos2 <- paste(snv$Chr,snv$Pos,sep="_")
  allele1 <- paste(snv[,1],snv[,2],snv[,3],snv[,4],sep="_")
  context1 <- rep(NA,nrow(snv))
  
  for (j in 1:length(bams)) {
    print(bams[j])
    mutect <- read.delim(paste0(bams[j],".mutect"),
			 as.is=TRUE,comment.char="#")
    allele2 <- paste(mutect[,1],mutect[,2],mutect[,4],mutect[,5],sep="_")
    context2 <- mutect$context[match(allele1, allele2)]
    context1[!is.na(context2)] <- context2[!is.na(context2)]
  }
  substr(context1,4,4) <- snv$Ref
  
  for (j in 1:length(bams)) {
    print(bams[j])
    bamfile <- paste0(bamdir,bams[j])
    if (! file.exists(bamfile))
      stop("Bam file doesn't exists")

    ##1F2R reads
    flag1 <- scanBamFlag(isProperPair=TRUE,
			 isMinusStrand=FALSE,isFirstMateRead=TRUE)
    flag2 <- scanBamFlag(isProperPair=TRUE,
			 isMinusStrand=TRUE,isSecondMateRead=TRUE)
    ##2F1R reads
    flag3 <- scanBamFlag(isProperPair=TRUE,
			 isMinusStrand=TRUE,isFirstMateRead=TRUE)
    flag4 <- scanBamFlag(isProperPair=TRUE,
			 isMinusStrand=FALSE,isSecondMateRead=TRUE)
    
    flags <- list(flag1,flag2,flag3,flag4)
    pileup_tables <- vector("list",4)

    for (k in 1:4) {
      pileup1 <- pileup(bamfile,
			scanBamParam=ScanBamParam(flag=flags[[k]],which=pos),
			pileupParam=PileupParam(max_depth=2000,
			  min_base_quality=20, min_mapq=40,
			  min_nucleotide_depth=0, min_minor_allele_depth=0,
			  distinguish_strands=FALSE,
			  include_deletions=FALSE,include_insertions=FALSE))
      pileup1$pos2 <- factor(paste(pileup1$seqnames,pileup1$pos,sep="_"),pos2)
      pileup1_table <- xtabs(count~pos2+nucleotide,data=pileup1)[,c("A","C","G","T")]
      pileup_tables[[k]] <- pileup1_table
    }

    pileup_table12 <- pileup_tables[[1]] + pileup_tables[[2]]
    pileup_table21 <- pileup_tables[[3]] + pileup_tables[[4]]
 
    idx_ref <- cbind(pos2,snv$Ref)
    idx_alt <- cbind(pos2,snv$Var)
    
    snv1 <- snv[,c(1,2,2,3,4)]
    names(snv1) <- c("Chromosome","Start_position","End_position",
		     "Reference_Allele","Tumor_Seq_Allele2")
    snv1 <- cbind(snv1,i_picard_oxoQ=0)
    snv1 <- cbind(snv1,
		  Tumor_Sample_Barcode=sns[j],
		  Matched_Norm_Sample_Barcode=normal,
		  ref_context=context1,
		  i_t_ALT_F1R2=pileup_table12[idx_alt],
		  i_t_ALT_F2R1=pileup_table21[idx_alt],
		  i_t_REF_F1R2=pileup_table12[idx_ref],
		  i_t_REF_F2R1=pileup_table21[idx_ref],
		  i_t_Foxog=rep(0,nrow(snv)),
		  Variant_Type="SNP")

    snv1$i_t_Foxog[snv$Ref %in% c("C","A")] <- 
      (snv1$i_t_ALT_F2R1 / (snv1$i_t_ALT_F1R2 + snv1$i_t_ALT_F2R1))[snv$Ref %in% c("C","A")]
    snv1$i_t_Foxog[snv$Ref %in% c("G","T")] <- 
      (snv1$i_t_ALT_F1R2 / (snv1$i_t_ALT_F1R2 + snv1$i_t_ALT_F2R1))[snv$Ref %in% c("G","T")]
    
    snv1 <- snv1[mutect_call[,j] =="yes",]
    all_snvs1 <- rbind(all_snvs1,snv1)

    if (j == 1) {
      pileup_table12_pool <- pileup_table12
      pileup_table21_pool <- pileup_table21
    }
    else {
      pileup_table12_pool <- pileup_table12_pool + pileup_table12
      pileup_table21_pool <- pileup_table21_pool + pileup_table21
    }
  }
  
  idx_ref <- cbind(pos2,snv$Ref)
  idx_alt <- cbind(pos2,snv$Var)
  
  snv2 <- snv[,c(1,2,2,3,4)]
  names(snv2) <- c("Chromosome","Start_position","End_position",
		   "Reference_Allele","Tumor_Seq_Allele2")
  snv2 <- cbind(snv2,i_picard_oxoQ=0)
  snv2 <- cbind(snv2,
		Tumor_Sample_Barcode=paste0(pts[i],"_Tumor"),
		Matched_Norm_Sample_Barcode=normal,
		ref_context=context1,
		i_t_ALT_F1R2=pileup_table12_pool[idx_alt],
		i_t_ALT_F2R1=pileup_table21_pool[idx_alt],
		i_t_REF_F1R2=pileup_table12_pool[idx_ref],
		i_t_REF_F2R1=pileup_table21_pool[idx_ref],
		i_t_Foxog=rep(0,nrow(snv)),
		Variant_Type="SNP")
  
  snv2$i_t_Foxog[snv$Ref %in% c("C","A")] <- 
    (snv2$i_t_ALT_F2R1 / (snv2$i_t_ALT_F1R2 + snv2$i_t_ALT_F2R1))[snv$Ref %in% c("C","A")]
  snv2$i_t_Foxog[snv$Ref %in% c("G","T")] <- 
    (snv2$i_t_ALT_F1R2 / (snv2$i_t_ALT_F1R2 + snv2$i_t_ALT_F2R1))[snv$Ref %in% c("G","T")]
  ##snv2 <- snv2[ ! is.na(snv2$i_t_Foxog),]

  all_snvs2 <- rbind(all_snvs2,snv2)
}

all_snvs1$Tumor_Seq_Allele1 <- all_snvs1$Tumor_Seq_Allele2
all_snvs2$Tumor_Seq_Allele1 <- all_snvs2$Tumor_Seq_Allele2

##Individual samples
all_snvs1$i_t_Foxog <- 1 - all_snvs1$i_t_Foxog
write.table(all_snvs1,file="mCRC_Vienna_CtoT_1.txt",
	    quote=FALSE,sep="\t",row.names=FALSE)
