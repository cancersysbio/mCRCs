##Plot tumor site average CNA from MRS
##Remove small segments

library(GenomicRanges)

seg_cutoff <- 1e6

mCRC_dir <- Sys.getenv("mCRC_DIR")
setwd(file.path(mCRC_dir,"plot"))

mycolors <- c("#D55E00","#0072B2","#009E73","#E69F00","#CC79A7","#F0E442",
              "#999999","#FFFFFF")

short_list <- scan("../results/CRCGenes.txt","c")

samples <- read.delim("../results/alltumorlabels.txt",as.is=TRUE)

load("../results/mCRC_All_TitanCNA.RData")

source(file.path(mCRC_dir,"/scripts/data_processing/TitanCNA_05_cna2genes.R"))
source(file.path(mCRC_dir,"/scripts/data_processing/TitanCNA_04_TitanCNA2seg.R"))
chr_sizes <- read.table(file.path(mCRC_dir,"/ref_genomes/broad.Homo_sapiens_assembly19.fa.sizes",as.is=TRUE))[1:22,2]
chr_cumsize <- cumsum(c(0,as.numeric(chr_sizes)))

##adjust segmean for overall ploidy and purity

all_cna <- vector("list",length=nrow(samples))
names(all_cna) <- samples$Label

for (i in 1:nrow(samples)) {
  sn1 <- samples$Bam[i]
  cna <- titancna2seg(alltitancnaresults[[sn1]]$results,
		      alltitancnaresults[[sn1]]$convergeParams)
  cna$chrom <- as.integer(as.character(cna$chrom))
  cna$loc.end <- pmin(cna$loc.end,chr_sizes[cna$chrom])
  cna$allelicratio <- round(cna$allelicratio,3)
  cna$logcopynumberratio <- round(cna$logcopynumberratio,3)

  ploidy1 <- cna$ploidy[1]
  purity1 <- cna$normalproportion[1]
  ploidy2 <- ploidy1 * (1 - purity1) + 2 * purity1

  prev1 <- cna$cellularprevalence * (1 - purity1)
  purity2 <- max(prev1,na.rm=TRUE)
  prev1[is.na(prev1)] <- purity2
  seg_mean <- cna$seg.mean + log2(ploidy2/2)
  seg_mean_adj <- log2(pmax(2^-2,1+(2^seg_mean-1)/purity2))

  cna$seg.mean.adj <- round(seg_mean_adj,3)
  all_cna[[i]] <- cna
}

sns <- samples$Label
patients <- sapply(strsplit(sns,"_",fixed=TRUE),"[",1)
sites <- sapply(strsplit(sns,"_",fixed=TRUE),"[",2)

allpatients <- unique(patients)

all_diff_genes <- NULL

for (pt1 in allpatients) {
  print(pt1)
  all_gr <- NULL
  sns1 <- sns[patients == pt1]
  for (sn1 in sns1) {
    cna <- all_cna[[sn1]]
    gr <- GRanges(seqnames=cna$chrom,
		  ranges=IRanges(cna$loc.start,cna$loc.end),
		  strand="*")
    if (is.null(all_gr)) {
      all_gr <- gr
    }
    else {
      all_gr <- disjoin(c(gr,all_gr))
    }
  }
  all_gr <- all_gr[width(all_gr)>=seg_cutoff]
  
  cna_seg <- data.frame(Chr=seqnames(all_gr),Start=start(all_gr),
			End=end(all_gr))
  sites1 <- sites[patients == pt1]
  for (site1 in unique(sites1)) {
    n_site1 <- table(sites1)[site1]
    site_seg <- rep(0,nrow(cna_seg))
    for (sn1 in sns1[sites1 == site1]) {
      cna <- all_cna[[sn1]]
      gr <- GRanges(seqnames=cna$chrom,
		    ranges=IRanges(cna$loc.start,cna$loc.end),
		    strand="*")
      idx <- findOverlaps(all_gr,gr)
      site_seg[queryHits(idx)] <- site_seg[queryHits(idx)] + cna$seg.mean.adj[subjectHits(idx)]/n_site1
    }
    cna_seg <- cbind(cna_seg,site_seg,2^site_seg*2,
		     2^site_seg*2-median(2^site_seg*2))
    names(cna_seg)[ncol(cna_seg)-(2:0)] <- paste0(site1,c("_lr","_cn","_rcn"))
  }

  site0 <- "P"
  for (site1 in unique(sites1)[-1]) {
    idx1 <- which(cna_seg[[paste0(site1,"_cn")]] - cna_seg[[paste0(site0,"_cn")]] <= -0.8 & cna_seg[[paste0(site1,"_rcn")]] - cna_seg[[paste0(site0,"_rcn")]] <= -0.8 & cna_seg[[paste0(site1,"_cn")]] <= 1.2 & cna_seg[[paste0(site1,"_rcn")]] <= -0.8)
    idx2 <- which(cna_seg[[paste0(site1,"_cn")]] - cna_seg[[paste0(site0,"_cn")]] >= 0.8 & cna_seg[[paste0(site1,"_rcn")]] - cna_seg[[paste0(site0,"_rcn")]] >= 0.8 & cna_seg[[paste0(site1,"_cn")]] >= 2.8 & cna_seg[[paste0(site1,"_rcn")]] >= 0.8)
    idx <- rep(0,nrow(cna_seg))
    idx[idx1] <- -1
    idx[idx2] <- 1
    cna_seg <- cbind(cna_seg,idx)
    names(cna_seg)[ncol(cna_seg)] <- paste0(site1,"_diff")
  }

  diff_genes <- NULL
  for (site1 in unique(sites1)) {
    if (site1 != "P") {
      cndiff <- cna_seg[[paste0(site1,"_diff")]] != 0
      if (!any(cndiff)) next
      gns <- cna2genes(cna_seg[cndiff,1:4],short_list)
      gns <- gns[!is.na(gns[,5]),1:4]
      gns <- gns[!duplicated(gns[,1]),]
      diff_genes <- c(diff_genes,as.character(gns[,1]))
    }
  }
  diff_genes <- unique(diff_genes)
  all_diff_genes <- c(all_diff_genes,diff_genes)
}
all_diff_genes <- names(which(table(all_diff_genes)>=2))


pdf("mCRC_All_CNA_Diff_2.pdf",width=10,height=6)

for (pt1 in allpatients) {
  print(pt1)
  all_gr <- NULL
  sns1 <- sns[patients == pt1]
  for (sn1 in sns1) {
    cna <- all_cna[[sn1]]
    gr <- GRanges(seqnames=cna$chrom,
		  ranges=IRanges(cna$loc.start,cna$loc.end),
		  strand="*")
    if (is.null(all_gr)) {
      all_gr <- gr
    }
    else {
      all_gr <- disjoin(c(gr,all_gr))
    }
  }
  all_gr <- all_gr[width(all_gr)>=seg_cutoff]

  cna_seg <- data.frame(Chr=seqnames(all_gr),Start=start(all_gr),
			End=end(all_gr))
  sites1 <- sites[patients == pt1]
  for (site1 in unique(sites1)) {
    n_site1 <- table(sites1)[site1]
    site_seg <- rep(0,nrow(cna_seg))
    for (sn1 in sns1[sites1 == site1]) {
      cna <- all_cna[[sn1]]
      gr <- GRanges(seqnames=cna$chrom,
		    ranges=IRanges(cna$loc.start,cna$loc.end),
		    strand="*")
      idx <- findOverlaps(all_gr,gr)
      site_seg[queryHits(idx)] <- site_seg[queryHits(idx)] + cna$seg.mean.adj[subjectHits(idx)]/n_site1
    }
    cna_seg <- cbind(cna_seg,site_seg,2^site_seg*2,
		     2^site_seg*2-median(2^site_seg*2))
    names(cna_seg)[ncol(cna_seg)-(2:0)] <- paste0(site1,c("_lr","_cn","_rcn"))
  }

  cna_seg$Pos1 <- cna_seg$Start + chr_cumsize[cna_seg$Chr]
  cna_seg$Pos2 <- cna_seg$End + chr_cumsize[cna_seg$Chr]

  site0 <- "P"
  for (site1 in unique(sites1)[-1]) {
    idx1 <- which(cna_seg[[paste0(site1,"_cn")]] - cna_seg[[paste0(site0,"_cn")]] <= -0.8 & cna_seg[[paste0(site1,"_rcn")]] - cna_seg[[paste0(site0,"_rcn")]] <= -0.8 & cna_seg[[paste0(site1,"_cn")]] <= 1.2 & cna_seg[[paste0(site1,"_rcn")]] <= -0.8)
    idx2 <- which(cna_seg[[paste0(site1,"_cn")]] - cna_seg[[paste0(site0,"_cn")]] >= 0.8 & cna_seg[[paste0(site1,"_rcn")]] - cna_seg[[paste0(site0,"_rcn")]] >= 0.8 & cna_seg[[paste0(site1,"_cn")]] >= 2.8 & cna_seg[[paste0(site1,"_rcn")]] >= 0.8)
    idx <- rep(0,nrow(cna_seg))
    idx[idx1] <- -1
    idx[idx2] <- 1
    cna_seg <- cbind(cna_seg,idx)
    names(cna_seg)[ncol(cna_seg)] <- paste0(site1,"_diff")
  }

  if (length(unique(sites1)) == 3)
    layout(matrix(1:5,ncol=1),
	   height=c(1.5,4,4,4,1.5))
  else 
    layout(matrix(1:5,ncol=1),
	   height=c(1.5,4,4,1.5,4))    
  par(mar=c(0,4,3,4))
  plot(0,0,xlim=c(0,sum(as.numeric(chr_sizes))),
       type="n",xaxt="n",yaxt="n",bty="n",ylab="",xlab="",
       main=pt1,cex.main=1.5)
  seg_length <- cna_seg$End - cna_seg$Start
  for (site1 in unique(sites1)) {
    par(mar=c(0,4,0,4))
    plot(0,0,ylim=c(-2,2),xlim=c(0,sum(as.numeric(chr_sizes))),
	 type="n",xaxt="n",bty="n",ylab="Log ratio",xlab="",las=2)
    axis(4,0,site1,las=2,tick=FALSE,cex.axis=1.5)
    lr <- cna_seg[[paste0(site1,"_lr")]]
    idx1 <- lr >= 0
    idx2 <- lr < 0
    rect(cna_seg$Pos1[idx1],0,cna_seg$Pos2[idx1],lr[idx1],
	 col=mycolors[1],border=NA)
    if (any(idx2))
      rect(cna_seg$Pos1[idx2],0,cna_seg$Pos2[idx2],lr[idx2],
	   col=mycolors[2],border=NA)
    abline(h=weighted.mean(lr,seg_length),lty=2,lwd=2,col="grey")
  }
  par(mar=c(3,4,0,4))
  plot(0,0,xlim=c(0,sum(as.numeric(chr_sizes))),
       type="n",xaxt="n",yaxt="n",bty="n",ylab="",xlab="",las=2)
  axis(1,at=chr_cumsize,labels=FALSE)
  axis(1,at=chr_cumsize[1:22] + chr_sizes/2, labels=1:22, tick=FALSE)
}

dev.off()
