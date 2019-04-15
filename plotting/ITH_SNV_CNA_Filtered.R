##Horizontal orientation
##Add LOH tracks to CNA ITH plot

mycolors <- c("#D55E00","#0072B2","#009E73","#E69F00","#CC79A7","#F0E442",
	      "#999999","#FFFFFF")

mutation_order <- function(freqs,present,ambiv,tissues,tissue_unique)
{
  present_tissue <- matrix(FALSE,nrow=nrow(present),ncol=length(tissue_unique))
  for (i in 1:length(tissue_unique)) {
    present_tissue[,i] <- apply(present[,tissues %in% tissue_unique[i],drop=FALSE],1,any)
  }
  snv_type <- rep(10,nrow(present))
  snv_type[apply(present_tissue,1,all)] <- 1
  if (length(tissue_unique) == 3) {
    snv_type[apply(present_tissue[,c(1,2)],1,all) & ! present_tissue[,3]] <- 2
    snv_type[apply(present_tissue[,c(1,3)],1,all) & ! present_tissue[,2]] <- 3
    snv_type[apply(present_tissue[,c(2,3)],1,all) & ! present_tissue[,1]] <- 4
  }
  for (i in 1:length(tissue_unique)) {
    snv_type[present_tissue[,i] & apply(!present_tissue[,-i,drop=FALSE],1,all)] <- 4+i
  }

  n_sample <- -apply(present,1,sum)
  n_ambiv <- -apply(ambiv,1,sum)
  freqs2 <- freqs
  freqs2[!present] <- 0
  mean_freqs <- -apply(freqs2,1,mean,na.rm=TRUE)
  snv_order <- do.call(order,as.data.frame(cbind(snv_type,n_sample,!present,mean_freqs)))
  return(list(snv_order,snv_type))
}

mutation_plot <- function(freqs,present,ambiv,snv_type,genes,types,main,n_max=13,newplot=TRUE,addlegend=TRUE) 
{
  samples <- colnames(present)
  patient <- strsplit(samples[1],"_",fixed=TRUE)[[1]][1]
  ##samples <- sub(paste(patient,"_",sep=""),"",samples)
  n_samples <- length(samples)
  n_snv <- nrow(present)

  if (newplot) {
    par(mar=c(1,5,3,3),mgp=c(0,0,0))
    plot(1,1,xlim=c(0.5,n_snv+0.5),ylim=c(n_max+2,-0.5),type="n",
	 xlab="",ylab="",axes=FALSE)
  }
  ##axis(4,(n_samples+1)/2,main,tick=FALSE,cex.axis=1.2,font=2,line=1)
  tmp <- 0.15
  for (i in unique(snv_type)) {
    lines(c(sum(snv_type<i)+1,sum(snv_type<=i)),rep(n_samples+0.7,2),
	  col=i+1,lwd=4)
    text(sum(snv_type<i)+1+sum(snv_type==i)/2,n_samples+0.7+0.15-tmp,
	 sum(snv_type==i),pos=1,cex=0.9)
    tmp <- -tmp
  }
  abline(h=1:(n_samples-1)+0.5,lty=2)
  axis(2,at=1:n_samples,labels=samples,tick=FALSE,las=2,mgp=c(0,0,0))
  axis(4,at=1:n_samples,labels=colSums(present),tick=FALSE,las=2,mgp=c(0,0,0),cex.axis=0.9)
  if (addlegend)
    legend(0,n_samples+1.5,c(expression(VAF>=0.1),expression(VAF<0.1)),
	   fill=c(mycolors[3],"#7FCEB9"),
	   border=NA,bty="n",cex=1)
  crc_genes <- NULL
  crc_genes_at <- NULL
  all_cols <- matrix("",nrow=n_snv,ncol=n_samples)
  for (i in 1:n_snv) {
    for (j in 1:n_samples) {
      if (present[i,j] & freqs[i,j] >= 0.1) 
	col1 <- mycolors[3]
      else if (present[i,j] & freqs[i,j] < 0.1) 
	col1 <- "#7FCEB9"
      else if (ambiv[i,j] && snv_type[i]==1) 
	col1 <- "grey"
      else col1 <- "white"
      all_cols[i,j] <- col1
    }
    if (genes[i] %in% short_list && max(freqs[i,]) >= 0.1 & 
	! types[i] %in% c("exonic_unknown","synonymous")) {
      crc_genes <- c(crc_genes,genes[i])
      crc_genes_at <- c(crc_genes_at,i)
    }
  }
  for (j in 1:n_samples) {
    rect((1:n_snv)-0.5,rep(j+0.45,n_snv),(1:n_snv)+0.5,rep(j-0.45,n_snv),
	 col=all_cols[,j],border=NA)
  }
  if (length(crc_genes) > 0) {
    str_height <- strwidth("W",cex=0.8)*1.2
    x_pos <- -str_height
    all_x_pos <- NULL
    for (i in 1:length(crc_genes)){
      x_pos1 <- crc_genes_at[i]
      if (x_pos1 - x_pos < str_height) 
	x_pos1 <- x_pos + str_height
      x_pos <- x_pos1
      all_x_pos <- c(all_x_pos,x_pos)
      lines(c(crc_genes_at[i],x_pos),c(0.5,-0.5),lty=2)
      ##text(x_pos1,n_samples+1.5,labels=crc_genes[i],pos=1,cex=0.8,
      ##	   adj=c(0.5,0.5),srt=90)
    }
    axis(3,at=all_x_pos,crc_genes,cex.axis=0.8,tick=FALSE,
	 pos=-0.5,mgp=c(0,0,0),las=2)
  }
}


mCRC_dir <- Sys.getenv("mCRC_DIR")
setwd(file.path(mCRC_dir,"plot"))

dir1 <- "../results/"
samples <- read.delim("../results/alltumorlabels.txt",as.is=TRUE)$Bam
samplelabels <- read.delim("../results/alltumorlabels.txt",as.is=TRUE)$Label


##Cases with primary samples
patients <- c("V46","V402","V514","V559","V750","V824","V855","V930","V953","V974","mCRCTB1","mCRCTB7","UchiCase2","CHET9","CHET40")

idx <- sapply(strsplit(samples,"_",fixed=TRUE),"[",1) %in% patients
samples <- samples[idx]
samplelabels <- samplelabels[idx]

tissue_types <- c("P","LN","LI","LU","BM")

##ITH diagram

pdf("mCRC_All_Filtered_ITH_SNV_Indel.pdf",height=8,width=8)

for (pt1 in patients) {
  
  fn1 <- paste(dir1,pt1,"_MuTectSNV_Indel_Coding_Filtered.txt",sep="")
  count1 <- read.delim(fn1,as.is=TRUE)
  count1 <- count1[count1$CtoT_filter == "Pass",]
  normal <- sub("_freq","",grep("freq",names(count1),value=TRUE)[1])
  sns <- sub("_mutect","",grep("mutect",names(count1),value=TRUE))

  sns <- samplelabels[match(sns,samples)]
  sns2 <- sns[!is.na(sns)]
  pt2 <- strsplit(sns2[1],"_",fixed=TRUE)[[1]][1]

  mutect <- count1[,grep("mutect",names(count1))]
  freqs <- count1[,grep("freq",names(count1))[-1]]
  cover <- count1[,grep("cover",names(count1))[-1]]  
  reads <- round(freqs*cover)
  present <- (freqs >= 0.05 & reads >= 3) | freqs >= 0.1
  ##ambiv <- array(FALSE,dim=dim(present))
  ##present <- (freqs >= 0.02 & reads >= 2) | freqs >= 0.1 | mutect == "yes"
  ambiv <- (reads > 0 & ! present) | (reads == 0 & cover < 10)
  colnames(present) <- colnames(ambiv) <- colnames(freqs) <- sns
  present <- present[,sns2]
  ambiv <- ambiv[,sns2]
  freqs <- freqs[,sns2]

  tissues <- sapply(strsplit(sns2,"_",fixed=TRUE),"[",2)
  tissue_unique <- tissue_types[tissue_types %in% tissues]
  
  temp <- mutation_order(freqs,present,ambiv,tissues,tissue_unique)
  snv_order <- temp[[1]]
  snv_type <- temp[[2]]
  snv_order <- snv_order[snv_type[snv_order] != 10] ##SNV with freq < 0.05 in all samples

  mutation_plot(freqs[snv_order,],present[snv_order,],ambiv[snv_order,],
		snv_type[snv_order],
                count1$Gene[snv_order],count1$Type[snv_order],pt2)


  ##Only keep SNV with >= 10X in all samples
  count1 <- count1[apply(cover>=10,1,all),]

  mutect <- count1[,grep("mutect",names(count1))]
  freqs <- count1[,grep("freq",names(count1))[-1]]
  cover <- count1[,grep("cover",names(count1))[-1]]  
  reads <- round(freqs*cover)
  present <- (freqs >= 0.05 & reads >= 3) | freqs >= 0.1
  ##ambiv <- array(FALSE,dim=dim(present))
  ambiv <- (reads > 0 & ! present) | (reads == 0 & cover < 10)
  ##present <- (freqs >= 0.02 & reads >= 2) | freqs >= 0.1 | mutect == "yes"
  ##ambiv <- (reads == 0 & cover < 10) | (reads > 0 & ! present)
  colnames(present) <- colnames(ambiv) <- colnames(freqs) <- sns
  present <- present[,sns2]
  ambiv <- ambiv[,sns2]
  freqs <- freqs[,sns2]
  
  tissues <- sapply(strsplit(sns2,"_",fixed=TRUE),"[",2)
  tissue_unique <- tissue_types[tissue_types %in% tissues]
  
  temp <- mutation_order(freqs,present,ambiv,tissues,tissue_unique)
  snv_order <- temp[[1]]
  snv_type <- temp[[2]]
  snv_order <- snv_order[snv_type[snv_order] != 10] ##SNV with freq < 0.05 in all samples

  mutation_plot(freqs[snv_order,],present[snv_order,],ambiv[snv_order,],
		snv_type[snv_order],
                count1$Gene[snv_order],count1$Type[snv_order],
                paste(pt2,"( >= 10X )")) 

}

dev.off()



##CNA

load("../results/mCRC_All_TitanCNA.RData")
source("../scripts/data_processing/TitanCNA_04_TitanCNA2seg.R")
chr_sizes <- read.table("../ref_genomes/broad.Homo_sapiens_assembly19.fa.sizes",as.is=TRUE)[1:22,2]
chr_cumsize <- cumsum(c(0,as.numeric(chr_sizes)))

cna_plot <- function(cna,j) {
  cna$loc.start <- cna$loc.start + chr_cumsize[cna$chrom]
  cna$loc.end <- pmin(cna$loc.end,chr_sizes[cna$chrom])
  cna$loc.end <- cna$loc.end + chr_cumsize[cna$chrom]
  
  cols <- rep(0,nrow(cna))
  idx1 <- cna$seg.mean < 0
  cols2 <- colorRamp(c(mycolors[c(2,8)]))(pmax(1+cna$seg.mean[idx1]/2,0))
  cols[idx1] <- rgb(red=cols2[,1],green=cols2[,2],blue=cols2[,3],
		    maxColorValue=255)
  idx2 <- cna$seg.mean >= 0
  cols2 <- colorRamp(c(mycolors[c(1,8)]))(pmax(1-cna$seg.mean[idx2]/2,0))
  cols[idx2] <- rgb(red=cols2[,1],green=cols2[,2],blue=cols2[,3],
		    maxColorValue=255)

  rect(cna$loc.start,j-0.45,cna$loc.end,j+0.15,col=cols,border=NA)

  cols <- rep("lightgray",nrow(cna))
  cols[cna$loh] <- "black"
  rect(cna$loc.start,j+0.2,cna$loc.end,j+0.4,col=cols,border=NA)
}
  

pdf("mCRC_All_Filtered_ITH_CNA.pdf",height=8,width=8)

for (i in 1:length(patients)) 
{
  pt1 <- patients[i]
  print(pt1)
  sns <- grep(pt1,samples,value=TRUE)
  sns2 <- samplelabels[match(sns,samples)]
  n_samples <- length(sns)
  tissues <- sapply(strsplit(sns,"_",fixed=TRUE),"[",2)
  pt2 <- strsplit(sns2[1],"_",fixed=TRUE)[[1]][1]

  n_max <- 13
  par(mar=c(1,5,3,3),mgp=c(0,0,0))
  plot(1,1,xlim=c(0,max(chr_cumsize)),ylim=c(n_max+2,-0.5),type="n",
       xlab="",ylab="",axes=FALSE)
  ##axis(3,(n_samples+1)/2,pt2,tick=FALSE,cex.axis=1.2,font=2,line=1)  
  axis(1,at=chr_cumsize,labels=rep("",length(chr_cumsize)),pos=n_samples+0.7)
  chrs2 <- c(1:10,seq(12,22,2))
  axis(1,at=(chr_cumsize[-1]-chr_sizes/2)[chrs2],chrs2,tick=FALSE,pos=n_samples+1)
  axis(2,at=1:length(sns),labels=sns2,las=2,tick=FALSE)
  abline(h=1:(n_samples-1)+0.5,lty=2)
  legend(0,n_samples+1.5,c("Gain","Loss","LOH"),fill=c(mycolors[1:2],"black"),
	 border=NA,bty="n",cex=1)
  for (j in 1:length(sns)) {
    sn1 <- sns[j]
    cna1 <- titancna2seg(alltitancnaresults[[sn1]]$results,
			 alltitancnaresults[[sn1]]$convergeParams)
    cna1$chrom <- as.integer(as.character(cna1$chrom))
    ploidy1 <- cna1$ploidy[1]
    purity1 <- cna1$normalproportion[1]
    ploidy2 <- ploidy1 * (1 - purity1) + 2 * purity1
    prev1 <- cna1$cellularprevalence * (1 - purity1)
    purity2 <- max(prev1,na.rm=TRUE)
    prev1[is.na(prev1)] <- purity2
    seg_mean <- cna1$seg.mean + log2(ploidy2/2)
    seg_mean_adj <- log2(pmax(2^-2,1+(2^seg_mean-1)/purity2))
    ##seg_mean_adj[1+(2^seg_mean-1)/prev1 <= 0] <- -2 ##Homozygous deletions

    loh1 <- cna1$minor_cn == 0 & prev1 == max(prev1)

    cna2 <- cbind(cna1[,1:3],seg.mean=seg_mean_adj,loh=loh1)

    cna_plot(cna2,j)
  }
  
}

dev.off()

pdf("mCRC_All_Filtered_ITH_SNV_Indel_CNA.pdf",height=8,width=12)

for (i in 1:length(patients)) 
{
  pt1 <- patients[i]
  print(pt1)
  fn1 <- paste(dir1,pt1,"_MuTectSNV_Indel_Coding_Filtered.txt",sep="")
  count1 <- read.delim(fn1,as.is=TRUE)
  count1 <- count1[count1$CtoT_filter == "Pass",]
  normal <- sub("_freq","",grep("freq",names(count1),value=TRUE)[1])
  sns <- sub("_mutect","",grep("mutect",names(count1),value=TRUE))

  sns1 <- samplelabels[match(sns,samples)]
  sns2 <- sns1[!is.na(sns1)]
  sns2 <- samplelabels[samplelabels %in% sns2]
  pt2 <- strsplit(sns2[1],"_",fixed=TRUE)[[1]][1]

  mutect <- count1[,grep("mutect",names(count1))]
  freqs <- count1[,grep("freq",names(count1))[-1]]
  cover <- count1[,grep("cover",names(count1))[-1]]  
  reads <- round(freqs*cover)
  present <- (freqs >= 0.05 & reads >= 3) | freqs >= 0.1
  ##ambiv <- array(FALSE,dim=dim(present))
  ##present <- (freqs >= 0.02 & reads >= 2) | freqs >= 0.1 | mutect == "yes"
  ambiv <- (reads > 0 & ! present) | (reads == 0 & cover < 10)
  colnames(present) <- colnames(ambiv) <- colnames(freqs) <- sns1
  present <- present[,sns2]
  ambiv <- ambiv[,sns2]
  freqs <- freqs[,sns2]

  tissues <- sapply(strsplit(sns2,"_",fixed=TRUE),"[",2)
  tissue_unique <- tissue_types[tissue_types %in% tissues]
  
  temp <- mutation_order(freqs,present,ambiv,tissues,tissue_unique)
  snv_order <- temp[[1]]
  snv_type <- temp[[2]]
  snv_order <- snv_order[snv_type[snv_order] != 10] ##SNV with freq < 0.05 in all samples

  n_max <- 25
  layout(matrix(1:2,nrow=1),width=c(0.55,0.45))

  mutation_plot(freqs[snv_order,],present[snv_order,],ambiv[snv_order,],
		snv_type[snv_order],
                count1$Gene[snv_order],count1$Type[snv_order],pt2,
		n_max=n_max)

  sns <- grep(pt1,samples,value=TRUE)
  sns2 <- samplelabels[match(sns,samples)]
  n_samples <- length(sns)
  tissues <- sapply(strsplit(sns,"_",fixed=TRUE),"[",2)
  pt2 <- strsplit(sns2[1],"_",fixed=TRUE)[[1]][1]

  par(mar=c(1,0,3,1),mgp=c(0,0,0))
  plot(1,1,xlim=c(0,max(chr_cumsize)),ylim=c(n_max+2,-0.5),type="n",
       xlab="",ylab="",axes=FALSE)
  axis(1,at=chr_cumsize,labels=rep("",length(chr_cumsize)),pos=n_samples+0.7)
  chrs2 <- c(1:10,seq(12,22,2))
  axis(1,at=(chr_cumsize[-1]-chr_sizes/2)[chrs2],chrs2,tick=FALSE,pos=n_samples+1)
  abline(h=1:(n_samples-1)+0.5,lty=2)
  legend(0,n_samples+1.5,c("Gain","Loss","LOH"),fill=c(mycolors[1:2],"black"),
	 border=NA,bty="n",cex=1)

  for (j in 1:length(sns)) {
    sn1 <- sns[j]
    cna1 <- titancna2seg(alltitancnaresults[[sn1]]$results,
			 alltitancnaresults[[sn1]]$convergeParams)
    cna1$chrom <- as.integer(as.character(cna1$chrom))
    ploidy1 <- cna1$ploidy[1]
    purity1 <- cna1$normalproportion[1]
    ploidy2 <- ploidy1 * (1 - purity1) + 2 * purity1
    prev1 <- cna1$cellularprevalence * (1 - purity1)
    purity2 <- max(prev1,na.rm=TRUE)
    prev1[is.na(prev1)] <- purity2
    seg_mean <- cna1$seg.mean + log2(ploidy2/2)
    seg_mean_adj <- log2(pmax(2^-2,1+(2^seg_mean-1)/purity2))
    ##seg_mean_adj[1+(2^seg_mean-1)/prev1 <= 0] <- -2 ##Homozygous deletions

    loh1 <- cna1$minor_cn == 0 & prev1 == max(prev1)

    cna2 <- cbind(cna1[,1:3],seg.mean=seg_mean_adj,loh=loh1)

    cna_plot(cna2,j)
  }
}

dev.off()
