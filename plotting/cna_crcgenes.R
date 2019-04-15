mCRC_dir <- Sys.getenv("mCRC_DIR")
setwd(file.path(mCRC_dir,"plot"))

genes <- scan("../results/CRCGenes.txt","c")

chr_sizes <- read.table(file.path(mCRC_dir,"ref_genomes/ref_genomes/broad.Homo_sapiens_assembly19.fa.sizes",as.is=TRUE)[1:22,2]
chr_cumsize <- cumsum(c(0,as.numeric(chr_sizes)))

source(file.path(mCRC_dir,"scripts/data_processing/TitanCNA_04_TitanCNA2seg.R"))
source(file.path(mCRC_dir,"scripts/data_processing/TitanCNA_05_cna2genes.R"))
a <- read.delim("../results/alltumorlabels.txt",as.is=TRUE)
a <- a[grep("^V",a$Bam),]

load("../results/mCRC_All_TitanCNA.RData")

alltitancnaresults <- alltitancnaresults[match(a$Bam,names(alltitancnaresults))]
names(alltitancnaresults) <- a$Label

allsamples <- names(alltitancnaresults)

alltitancnagenes <- matrix(NA,ncol=length(allsamples)+3,nrow=length(genes))
colnames(alltitancnagenes) <- c("Chr","Start","End",allsamples)
rownames(alltitancnagenes) <- genes

for (i in 1:length(allsamples)) {
  print(allsamples[i])
  cna <- titancna2seg(alltitancnaresults[[i]]$results,alltitancnaresults[[i]]$convergeParams)
  cna$chrom <- as.integer(as.character(cna$chrom))
  cna$loc.end <- pmin(cna$loc.end,chr_sizes[cna$chrom])
  cna$allelicratio <- round(cna$allelicratio,3)
  cna$logcopynumberratio <- round(cna$logcopynumberratio,3)

  ##adjust seg.mean for overall ploidy and purity
  ploidy1 <- cna$ploidy[1]
  purity1 <- cna$normalproportion[1]
  ploidy2 <- ploidy1 * (1 - purity1) + 2 * purity1

  prev1 <- cna$cellularprevalence * (1 - purity1)
  prev1[is.na(prev1)] <- max(prev1,na.rm=TRUE)
  seg_mean <- cna$seg.mean + log2(ploidy2/2)
  seg_mean_adj <- log2(1+(2^seg_mean-1)/prev1)
  seg_mean_adj[1+(2^seg_mean-1)/prev1 <= 0] <- -2 ##Homozygous deletions
  cna$seg.mean.adj <- round(seg_mean_adj,3)

  cnagenes <- cna2genes(cna,genes)
  dup_genes <- names(which(table(cnagenes$gene)>1))
  if (length(dup_genes)>0) {
    for (gene1 in dup_genes) {
      idx1 <- cnagenes$gene == gene1
      cnagenes$seg.mean.adj[idx1] <- cnagenes$seg.mean.adj[idx1][which.max(abs(cnagenes$seg.mean.adj[idx1]))]
      cnagenes$loc.start[idx1] <- min(cnagenes$loc.start[idx1])
      cnagenes$loc.end[idx1] <- max(cnagenes$loc.end[idx1])
    }
  }
  cnagenes <- cnagenes[!duplicated(cnagenes$gene),]
  cnagenes$chrom <- sub("X","23",cnagenes$chrom)
  cnagenes$chrom <- as.integer(cnagenes$chrom)
  alltitancnagenes[match(cnagenes$gene,genes),1:3] <- as.matrix(cnagenes[,2:4])
  alltitancnagenes[match(cnagenes$gene,genes),i+3] <- cnagenes$seg.mean.adj
}

allcnagenes <- alltitancnagenes[,c(1:3,grep("^V",colnames(alltitancnagenes)))]
##allcnagenes[,4:ncol(allcnagenes)] <- 2*2^allcnagenes[,4:ncol(allcnagenes)]
allcnagenes <- allcnagenes[order(allcnagenes[,1],allcnagenes[,2]),]
allcnagenes <- allcnagenes[!is.na(allcnagenes[,4]),]

##write.table(allcnagenes,file="allcnagenes.txt",quote=FALSE,sep="\t",row.names=FALSE)
##allcnagenes <- as.matrix(read.delim("allcnagenes.txt"))

allcnagenes_p <- as.data.frame(allcnagenes[,c(1:3,grep("_P_",colnames(allcnagenes)))])
allcnagenes_bm <- as.data.frame(allcnagenes[,c(1:3,grep("_BM_",colnames(allcnagenes)))])

sns <- colnames(allcnagenes)[-(1:3)]
patients <- sapply(strsplit(sns,"_",fixed=TRUE),"[",1)
sites <- sapply(strsplit(sns,"_",fixed=TRUE),"[",2)
allpatients <- unique(patients)

site_counts <- xtabs(~patients+sites)
npatients <- length(unique(patients))
sample_weight <- rep(1,length(sns))
names(sample_weight) <- sns
for (i in 1:length(sns)) {
  sample_weight[i] <- 1/site_counts[match(patients[i],rownames(site_counts)),
                               match(sites[i],colnames(site_counts))]
}

cutoff <- c(log2(1.2/2),log2(3.6/2))

sns_p <- grep("_P_",names(allcnagenes_p),value=TRUE)
allcnagenes_p$n_amp <- rowSums((allcnagenes_p[,sns_p]>=cutoff[2])*
				 rep(sample_weight[sns_p],each=nrow(allcnagenes_p)))
allcnagenes_p$n_del <- rowSums((allcnagenes_p[,sns_p]<=cutoff[1])*
				 rep(sample_weight[sns_p],each=nrow(allcnagenes_p)))

sns_bm <- grep("_BM_",names(allcnagenes_bm),value=TRUE)
allcnagenes_bm$n_amp <- rowSums((allcnagenes_bm[,sns_bm]>=cutoff[2])*
				 rep(sample_weight[sns_bm],each=nrow(allcnagenes_bm)))
allcnagenes_bm$n_del <- rowSums((allcnagenes_bm[,sns_bm]<=cutoff[1])*
				 rep(sample_weight[sns_bm],each=nrow(allcnagenes_bm)))


mycolors <- c("#D55E00","#0072B2","#009E73","#E69F00","#CC79A7","#F0E442",
              "#999999","#FFFFFF")

samples <- read.delim("alltumorlabels.txt",as.is=TRUE)
samples <- samples[grep("^V",samples$Bam),]

##adjust segmean for overall ploidy and purity

all_cna <- vector("list",nrow(samples))
names(all_cna) <- samples$Label

for (i in 1:nrow(samples)) {
  sn1 <- samples$Label[i]
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

  cna2 <- cbind(sample=samples$Label[i],cna[,1:4],
		seg.mean=round(seg_mean_adj,3))
  all_cna[[i]] <- cna2
}

sns <- samples$Label
patients <- sapply(strsplit(sns,"_",fixed=TRUE),"[",1)
sites <- sapply(strsplit(sns,"_",fixed=TRUE),"[",2)

##weights for avareging (averaging within patients then across patients)

counts <- xtabs(~patients+sites)
npatients <- length(unique(patients))

sample_weight <- rep(1,length(sns))
names(sample_weight) <- sns
for (i in 1:length(sns)) {
  sample_weight[i] <- 1/counts[match(patients[i],rownames(counts)),
			       match(sites[i],colnames(counts))]
}
sample_weight <- sample_weight/npatients

##Log ratios

primary_lr <- NULL
brain_lr <- NULL

for (i in 1:nrow(samples)) {
  sn1 <- samples$Label[i]
  print(sn1)
  cna <- all_cna[[sn1]][,-1]

  lr <- GRanges(seqnames=cna$chrom,
		ranges=IRanges(cna$loc.start,cna$loc.end),
		strand="*",seg.mean=cna$seg.mean*sample_weight[sn1])
  
  if (sites[i] == "P") {
    if (is.null(primary_lr)) {
      primary_lr <- lr
    }
    else {
      lr2 <- disjoin(c(lr,primary_lr))
      lr2$seg.mean <- 0
      idx <- findOverlaps(lr2,lr)
      lr2$seg.mean[queryHits(idx)] <- lr2$seg.mean[queryHits(idx)] + lr$seg.mean[subjectHits(idx)]
      idx <- findOverlaps(lr2,primary_lr)
      lr2$seg.mean[queryHits(idx)] <- lr2$seg.mean[queryHits(idx)] + primary_lr$seg.mean[subjectHits(idx)]
      primary_lr <- lr2
    }
  }
      
  if (sites[i] == "BM") {
    if (is.null(brain_lr)) {
      brain_lr <- lr
    }
    else {
      lr2 <- disjoin(c(lr,brain_lr))
      lr2$seg.mean <- 0
      idx <- findOverlaps(lr2,lr)
      lr2$seg.mean[queryHits(idx)] <- lr2$seg.mean[queryHits(idx)] + lr$seg.mean[subjectHits(idx)]
      idx <- findOverlaps(lr2,brain_lr)
      lr2$seg.mean[queryHits(idx)] <- lr2$seg.mean[queryHits(idx)] + brain_lr$seg.mean[subjectHits(idx)]
      brain_lr <- lr2
    }
  }
}

cutoff <- 1e5
primary_lr <- primary_lr[width(primary_lr)>=cutoff]
brain_lr <- brain_lr[width(brain_lr)>=cutoff]

sample_weight <- sample_weight*npatients

pdf("mCRC_Vienna_CNA_Genes.pdf",width=10,height=8)
layout(matrix(1:2,nrow=2),height=c(0.4,0.6))
par(mar=c(0,4,4,1))

a1 <- as.data.frame(primary_lr)
a1$seqnames <- as.integer(as.character(a1$seqnames))
a1$start <- a1$start + chr_cumsize[a1$seqnames]
a1$end <- a1$end + chr_cumsize[a1$seqnames]
idx1 <- a1$seg.mean > 0
idx2 <- a1$seg.mean < 0

plot(0,0,ylim=c(-0.5,1.5),xlim=c(0,sum(as.numeric(chr_sizes))),
     type="n",xaxt="n",bty="n",ylab="Average log ratio",xlab="",
     main="Primary")
rect(a1$start[idx1],0,a1$end[idx1],a1$seg.mean[idx1],col=mycolors[1],border=NA)
rect(a1$start[idx2],0,a1$end[idx2],a1$seg.mean[idx2],col=mycolors[2],border=NA)
abline(v=chr_cumsize,lty=3)
axis(3,at=chr_cumsize[-23]+chr_sizes/2,labels=1:22,tick=FALSE,mgp=c(0,0,0))

par(mar=c(5,4,0,1))

idx1 <- which(allcnagenes_p$n_amp >= 6 | (allcnagenes_p$n_del >= 2 & allcnagenes_p$n_amp < 1) | allcnagenes_bm$n_amp >= 6 | (allcnagenes_bm$n_del >= 2 & allcnagenes_bm$n_amp < 1))

x_ratio <- sum(as.numeric(chr_sizes))/(length(idx1)+1)

plot(0,0,xlim=c(0,sum(as.numeric(chr_sizes))),ylim=c(-2,4),type="n",
     xaxt="n",yaxt="n",ylab="",xlab="",bty="n")
axis(2,at=c(-2,-1,0,1,2))
title(ylab="Log ratio",adj=1/3)
cols1 <- mycolors[1] ##colorRamp(c(mycolors[c(8,1)]))((1:npatients+5)/(npatients+5))
##cols1 <- rgb(red=cols1[,1],green=cols1[,2],blue=cols1[,3],maxColorValue=255)
cols2 <- mycolors[2] ##colorRamp(c(mycolors[c(8,2)]))((1:npatients+5)/(npatients+5))
##cols2 <- rgb(red=cols2[,1],green=cols2[,2],blue=cols2[,3],maxColorValue=255)

for (i in 1:length(idx1)) {
  lines(rep(x_ratio*(i+0.5),2),c(-2,2),lty=3)
  if (i == 1)
    lines(rep(x_ratio*(i-0.5),2),c(-2,2),lty=3)
  lines(c(x_ratio*i,allcnagenes_p[idx1[i],2]+chr_cumsize[allcnagenes_p[idx1[i],1]]),c(2.2,4),lty=2)
  order1 <- order(allcnagenes_p[idx1[i],sns_p])
  a1 <- cbind(as.numeric(allcnagenes_p[idx1[i],sns_p[order1]]),
	      sample_weight[sns_p[order1]])
  if (sum(a1[a1[,1]>=0.5,2])>=1) {
    a1 <- a1[(sum(a1[,1]<0.5)+1):nrow(a1),,drop=FALSE]
    a1[,2] <- round(rev(cumsum(rev(a1[,2]))))
    a1 <- rbind(c(0.5,a1[1,2]),a1)
    a1[,1] <- pmin(2,a1[,1])
    w1 <- a1[2:nrow(a1),2]/npatients*0.4
    rect(x_ratio*(rep(i,nrow(a1)-1)-w1),a1[1:(nrow(a1)-1),1],
	 x_ratio*(rep(i,nrow(a1)-1)+w1),a1[2:nrow(a1),1],border=NA,
	 col=cols1)
  }
  order1 <- order(allcnagenes_p[idx1[i],sns_p],decreasing=TRUE)
  a1 <- cbind(as.numeric(allcnagenes_p[idx1[i],sns_p[order1]]),
	      sample_weight[sns_p[order1]])
  if (sum(a1[a1[,1]<=-0.5,2])>=1) {
    a1 <- a1[(sum(a1[,1]>-0.5)+1):nrow(a1),,drop=FALSE]
    a1[,2] <- round(rev(cumsum(rev(a1[,2]))))
    a1 <- rbind(c(-0.5,a1[1,2]),a1)
    w1 <- a1[2:nrow(a1),2]/npatients*0.4
    rect(x_ratio*(rep(i,nrow(a1)-1)-w1),a1[1:(nrow(a1)-1),1],
	 x_ratio*(rep(i,nrow(a1)-1)+w1),a1[2:nrow(a1),1],,border=NA,
	 col=cols2)
  }
}
axis(1,at=x_ratio*(1:length(idx1)),rownames(allcnagenes_p)[idx1],
     las=2,tick=FALSE,cex.axis=1,mgp=c(0,0,0))
abline(h=0,lty=2)


par(mar=c(0,4,4,1))

a1 <- as.data.frame(brain_lr)
a1$seqnames <- as.integer(as.character(a1$seqnames))
a1$start <- a1$start + chr_cumsize[a1$seqnames]
a1$end <- a1$end + chr_cumsize[a1$seqnames]
idx1 <- a1$seg.mean > 0
idx2 <- a1$seg.mean < 0

plot(0,0,ylim=c(-0.5,1.5),xlim=c(0,sum(as.numeric(chr_sizes))),
     type="n",xaxt="n",bty="n",ylab="Average log ratio",xlab="Chromosome",
     main="Brain metastasis")
rect(a1$start[idx1],0,a1$end[idx1],a1$seg.mean[idx1],col=mycolors[1],border=NA)
rect(a1$start[idx2],0,a1$end[idx2],a1$seg.mean[idx2],col=mycolors[2],border=NA)
abline(v=chr_cumsize,lty=3)
axis(3,at=chr_cumsize[-23]+chr_sizes/2,labels=1:22,tick=FALSE,mgp=c(0,0,0))


par(mar=c(5,4,0,1))

idx1 <- which(allcnagenes_p$n_amp >= 6 | (allcnagenes_p$n_del >= 2 & allcnagenes_p$n_amp < 1) | allcnagenes_bm$n_amp >= 6 | (allcnagenes_bm$n_del >= 2 & allcnagenes_bm$n_amp < 1))

x_ratio <- sum(as.numeric(chr_sizes))/(length(idx1)+1)

plot(0,0,xlim=c(0,sum(as.numeric(chr_sizes))),ylim=c(-2,4),type="n",
     xaxt="n",yaxt="n",ylab="",xlab="",bty="n")
axis(2,at=c(-2,-1,0,1,2))
title(ylab="Log ratio",adj=1/3)
##cols1 <- colorRamp(c(mycolors[c(8,1)]))((1:npatients+5)/(npatients+5))
##cols1 <- rgb(red=cols1[,1],green=cols1[,2],blue=cols1[,3],maxColorValue=255)
##cols2 <- colorRamp(c(mycolors[c(8,2)]))((1:npatients+5)/(npatients+5))
##cols2 <- rgb(red=cols2[,1],green=cols2[,2],blue=cols2[,3],maxColorValue=255)

for (i in 1:length(idx1)) {
  lines(rep(x_ratio*(i+0.5),2),c(-2,2),lty=3)
  if (i == 1)
    lines(rep(x_ratio*(i-0.5),2),c(-2,2),lty=3)
  lines(c(x_ratio*i,allcnagenes_bm[idx1[i],2]+chr_cumsize[allcnagenes_bm[idx1[i],1]]),c(2.2,4),lty=2)
  order1 <- order(allcnagenes_bm[idx1[i],sns_bm])
  a1 <- cbind(as.numeric(allcnagenes_bm[idx1[i],sns_bm[order1]]),
	      sample_weight[sns_bm[order1]])
  if (sum(a1[a1[,1]>=0.5,2])>=1) {
    a1 <- a1[(sum(a1[,1]<0.5)+1):nrow(a1),,drop=FALSE]
    a1[,2] <- round(rev(cumsum(rev(a1[,2]))))
    a1 <- rbind(c(0.5,a1[1,2]),a1)
    a1[,1] <- pmin(2,a1[,1])
    w1 <- a1[2:nrow(a1),2]/npatients*0.4
    rect(x_ratio*(rep(i,nrow(a1)-1)-w1),a1[1:(nrow(a1)-1),1],
	 x_ratio*(rep(i,nrow(a1)-1)+w1),a1[2:nrow(a1),1],border=NA,
	 col=cols1)
  }
  order1 <- order(allcnagenes_bm[idx1[i],sns_bm],decreasing=TRUE)
  a1 <- cbind(as.numeric(allcnagenes_bm[idx1[i],sns_bm[order1]]),
	      sample_weight[sns_bm[order1]])
  if (sum(a1[a1[,1]<=-0.5,2])>=1) {
    a1 <- a1[(sum(a1[,1]>-0.5)+1):nrow(a1),,drop=FALSE]
    a1[,2] <- round(rev(cumsum(rev(a1[,2]))))
    a1 <- rbind(c(-0.5,a1[1,2]),a1)
    w1 <- a1[2:nrow(a1),2]/npatients*0.4
    rect(x_ratio*(rep(i,nrow(a1)-1)-w1),a1[1:(nrow(a1)-1),1],
	 x_ratio*(rep(i,nrow(a1)-1)+w1),a1[2:nrow(a1),1],border=NA,
	 col=cols2)
  }
}
axis(1,at=x_ratio*(1:length(idx1)),rownames(allcnagenes_bm)[idx1],
     las=2,tick=FALSE,cex.axis=1,mgp=c(0,0,0))
abline(h=0,lty=2)

dev.off()
