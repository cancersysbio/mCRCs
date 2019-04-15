##Combine all primary samples and all brain met samples

library(GenomicRanges)

mCRC_dir <- Sys.getenv("mCRC_DIR")
setwd(file.path(mCRC_dir,"plot"))

mycolors <- c("#D55E00","#0072B2","#009E73","#E69F00","#CC79A7","#F0E442",
              "#999999","#FFFFFF")

short_list <- scan("../results/CRCGenes.txt","c")

samples <- read.delim("../results/alltumorlabels.txt",as.is=TRUE)
samples <- samples[grep("^V",samples$Bam),]

load("../results/mCRC_All_TitanCNA.RData")

source(file.path(mCRC_dir,"/scripts/data_processing/TitanCNA_05_cna2genes.R"))
source(file.path(mCRC_dir,"/scripts/data_processing/TitanCNA_04_TitanCNA2seg.R"))
chr_sizes <- read.table(file.path(mCRC_dir,"/ref_genomes/broad.Homo_sapiens_assembly19.fa.sizes",as.is=TRUE))[1:22,2]
chr_cumsize <- cumsum(c(0,as.numeric(chr_sizes)))

##adjust segmean for overall ploidy and purity

all_cna <- vector("list",nrow(samples))
names(all_cna) <- samples$Label

for (i in 1:nrow(samples)) {
  sn1 <- samples$Bam[i]
  cna <- titancna2seg(alltitancnaresults[[sn1]]$results,
		      alltitancnaresults[[sn1]]$convergeParams)
  cna$chrom <- as.integer(as.character(cna$chrom))
  cna$loc.end <- pmin(cna$loc.end,chr_sizes[cna$chrom,2])
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

##Thresholds of -0.8 and 0.8

counts <- xtabs(~patients+sites)
npatients <- length(unique(patients))

sample_weight <- rep(1,length(sns))
names(sample_weight) <- sns
for (i in 1:length(sns)) {
  sample_weight[i] <- 1/counts[match(patients[i],rownames(counts)),
			       match(sites[i],colnames(counts))]
}

primary_lr2 <- as.data.frame(primary_lr)[,1:4]
primary_lr2$seqnames <- as.integer(as.character(primary_lr2$seqnames))
brain_lr2 <- as.data.frame(brain_lr)[,1:4]
brain_lr2$seqnames <- as.integer(as.character(brain_lr2$seqnames))

for (i in 1:nrow(samples)) {
  sn1 <- samples$Label[i]
  print(sn1)
  cna <- all_cna[[sn1]][,-1]

  lr <- GRanges(seqnames=cna$chrom,
		ranges=IRanges(cna$loc.start,cna$loc.end),
		strand="*",seg.mean=cna$seg.mean)
  
  if (sites[i] == "P") {
    primary_lr2$seg.mean <- 0
    idx <- findOverlaps(primary_lr,lr)
    primary_lr2$seg.mean[queryHits(idx)] <- cna$seg.mean[subjectHits(idx)]
    colnames(primary_lr2)[ncol(primary_lr2)] <- sn1
  }
      
  if (sites[i] == "BM") {
    brain_lr2$seg.mean <- 0
    idx <- findOverlaps(brain_lr,lr)
    brain_lr2$seg.mean[queryHits(idx)] <- cna$seg.mean[subjectHits(idx)]
    colnames(brain_lr2)[ncol(brain_lr2)] <- sn1
  }
}

idx <- grep("^V",names(primary_lr2))
primary_lr2$N_amp <- rowSums((primary_lr2[,idx]>=0.8)*rep(sample_weight[names(primary_lr2)[idx]],each=nrow(primary_lr2)))
primary_lr2$N_del <- rowSums((primary_lr2[,idx]<=-0.8)*rep(sample_weight[names(primary_lr2)[idx]],each=nrow(primary_lr2)))
primary_lr2 <- primary_lr2[primary_lr2$N_amp >= 2 | primary_lr2$N_del >= 2,]
primary_lr2[,idx] <- apply(primary_lr2[,idx],2,pmin,2)
primary_lr2$cn.amp <- 
  rowSums((2^primary_lr2[,idx]*2-2)*(primary_lr2[,idx]>=0.8)*rep(sample_weight[names(primary_lr2)[idx]],each=nrow(primary_lr2)))
primary_lr2$cn.del <- 
  rowSums((2^primary_lr2[,idx]*2-2)*(primary_lr2[,idx]<=-0.8)*rep(sample_weight[names(primary_lr2)[idx]],each=nrow(primary_lr2)))


idx <- grep("^V",names(brain_lr2))
brain_lr2$N_amp <- rowSums((brain_lr2[,idx]>=0.8)*rep(sample_weight[names(brain_lr2)[idx]],each=nrow(brain_lr2)))
brain_lr2$N_del <- rowSums((brain_lr2[,idx]<=-0.8)*rep(sample_weight[names(brain_lr2)[idx]],each=nrow(brain_lr2)))
brain_lr2 <- brain_lr2[brain_lr2$N_amp >= 2 | brain_lr2$N_del >= 2,]
brain_lr2[,idx] <- apply(brain_lr2[,idx],2,pmin,2)
brain_lr2$cn.amp <- 
  rowSums((2^brain_lr2[,idx]*2-2)*(brain_lr2[,idx]>=0.8)*rep(sample_weight[names(brain_lr2)[idx]],each=nrow(brain_lr2)))
brain_lr2$cn.del <- 
  rowSums((2^brain_lr2[,idx]*2-2)*(brain_lr2[,idx]<=-0.8)*rep(sample_weight[names(brain_lr2)[idx]],each=nrow(brain_lr2)))


pdf("mCRC_Vienna_CNA_Comparison1.pdf",width=10,height=8)
par(mfrow=c(2,1),mar=c(4,4,2,1))

a1 <- as.data.frame(primary_lr)
a1$seqnames <- as.integer(as.character(a1$seqnames))
a1$start <- a1$start + chr_cumsize[a1$seqnames]
a1$end <- a1$end + chr_cumsize[a1$seqnames]
idx1 <- a1$seg.mean > 0
idx2 <- a1$seg.mean < 0

plot(0,0,ylim=c(-0.5,1.5),xlim=c(0,sum(as.numeric(chr_sizes[,2]))),
     type="n",xaxt="n",bty="n",ylab="Average log ratio",xlab="Chromosome",
     main="Primary")
rect(a1$start[idx1],0,a1$end[idx1],a1$seg.mean[idx1],col=mycolors[1],border=NA)
rect(a1$start[idx2],0,a1$end[idx2],a1$seg.mean[idx2],col=mycolors[2],border=NA)
abline(v=chr_cumsize[-1],lty=3)
axis(1,at=chr_cumsize+chr_sizes[,2]/2,labels=1:22,tick=FALSE)

a1 <- as.data.frame(brain_lr)
a1$seqnames <- as.integer(as.character(a1$seqnames))
a1$start <- a1$start + chr_cumsize[a1$seqnames]
a1$end <- a1$end + chr_cumsize[a1$seqnames]
idx1 <- a1$seg.mean > 0
idx2 <- a1$seg.mean < 0

plot(0,0,ylim=c(-0.5,1.5),xlim=c(0,sum(as.numeric(chr_sizes[,2]))),
     type="n",xaxt="n",bty="n",ylab="Average log ratio",xlab="Chromosome",
     main="Brain metastasis")
rect(a1$start[idx1],0,a1$end[idx1],a1$seg.mean[idx1],col=mycolors[1],border=NA)
rect(a1$start[idx2],0,a1$end[idx2],a1$seg.mean[idx2],col=mycolors[2],border=NA)
abline(v=chr_cumsize[-1],lty=3)
axis(1,at=chr_cumsize+chr_sizes[,2]/2,labels=1:22,tick=FALSE)

dev.off()

pdf("mCRC_Vienna_CNA_Comparison2.pdf",width=10,height=8)
par(mfrow=c(2,1),mar=c(4,4,2,1))

a1 <- primary_lr2
a1$start <- a1$start + chr_cumsize[a1$seqnames]
a1$end <- a1$end + chr_cumsize[a1$seqnames]
idx1 <- a1$N_amp >= 2
idx2 <- a1$N_del >= 2

plot(0,0,ylim=c(-10,40),
     xlim=c(0,sum(as.numeric(chr_sizes[,2]))),
     type="n",xaxt="n",yaxt="n",bty="n",
     ylab="Relative copy number",xlab="Chromosome",main="Primary")
axis(2,at=c(-10,0,10,20,30,40),labels=c(-5,0,10,20,30,40))
cols <- colorRamp(c(mycolors[c(8,1)]))((a1$N_amp[idx1]+5)/(npatients+5))
cols <- rgb(red=cols[,1],green=cols[,2],blue=cols[,3],maxColorValue=255)
rect(a1$start[idx1],0,a1$end[idx1],a1$cn.amp[idx1],
     col=cols,border=NA)
cols <- colorRamp(c(mycolors[c(8,2)]))((a1$N_del[idx2]+5)/(npatients+5))
cols <- rgb(red=cols[,1],green=cols[,2],blue=cols[,3],maxColorValue=255)
rect(a1$start[idx2],0,a1$end[idx2],a1$cn.del[idx2]*2,
     col=cols,border=NA)
abline(v=chr_cumsize[-1],lty=3)
axis(1,at=chr_cumsize+chr_sizes[,2]/2,labels=1:22,tick=FALSE)

y <- seq(20,40,length.out=11)
cols1 <- colorRamp(c(mycolors[c(8,1)]))((1:10+5)/(npatients+5))
cols1 <- rgb(red=cols1[,1],green=cols1[,2],blue=cols1[,3],maxColorValue=255)
cols2 <- colorRamp(c(mycolors[c(8,2)]))((1:10+5)/(npatients+5))
cols2 <- rgb(red=cols2[,1],green=cols2[,2],blue=cols2[,3],maxColorValue=255)
for (i in 2:10) {
  lines(c(0,0),y[c(i,i+1)],col=cols1[i],lwd=5)
  lines(c(5e7,5e7),y[c(i,i+1)],col=cols2[i],lwd=5)
}
for (i in seq(2,10,2)) {
  text(1e8,mean(y[c(i,i+1)]),i,cex=0.8)
}


a1 <- brain_lr2
a1$start <- a1$start + chr_cumsize[a1$seqnames]
a1$end <- a1$end + chr_cumsize[a1$seqnames]
idx1 <- a1$N_amp >= 2
idx2 <- a1$N_del >= 2

plot(0,0,ylim=c(-10,40),
     xlim=c(0,sum(as.numeric(chr_sizes[,2]))),
     type="n",xaxt="n",yaxt="n",bty="n",
     ylab="Relative copy number",xlab="Chromosome",main="Brain metastasis")
axis(2,at=c(-10,0,10,20,30,40),labels=c(-5,0,10,20,30,40))
cols <- colorRamp(c(mycolors[c(8,1)]))((a1$N_amp[idx1]+5)/(npatients+5))
cols <- rgb(red=cols[,1],green=cols[,2],blue=cols[,3],maxColorValue=255)
rect(a1$start[idx1],0,a1$end[idx1],a1$cn.amp[idx1],
     col=cols,border=NA)
cols <- colorRamp(c(mycolors[c(8,2)]))((a1$N_del[idx2]+5)/(npatients+5))
cols <- rgb(red=cols[,1],green=cols[,2],blue=cols[,3],maxColorValue=255)
rect(a1$start[idx2],0,a1$end[idx2],a1$cn.del[idx2]*2,
     col=cols,border=NA)
abline(v=chr_cumsize[-1],lty=3)
axis(1,at=chr_cumsize+chr_sizes[,2]/2,labels=1:22,tick=FALSE)

y <- seq(20,40,length.out=11)
cols1 <- colorRamp(c(mycolors[c(8,1)]))((1:10+5)/(npatients+5))
cols1 <- rgb(red=cols1[,1],green=cols1[,2],blue=cols1[,3],maxColorValue=255)
cols2 <- colorRamp(c(mycolors[c(8,2)]))((1:10+5)/(npatients+5))
cols2 <- rgb(red=cols2[,1],green=cols2[,2],blue=cols2[,3],maxColorValue=255)
for (i in 2:10) {
  lines(c(0,0),y[c(i,i+1)],col=cols1[i],lwd=5)
  lines(c(5e7,5e7),y[c(i,i+1)],col=cols2[i],lwd=5)
}
for (i in seq(2,10,2)) {
  text(1e8,mean(y[c(i,i+1)]),i,cex=0.8)
}

dev.off()
