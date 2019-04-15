##Combine multi-sample data 
##One point for each tumor site
##Include indel calls from Strelka
##Use more stringent criteria to include genes
##Require a SNV to be in all samples for one tumor site and have AF >= 0.1 in at least one such sample

##Shorter list of genes
##Include all FFPE and FF samples. No Uchi samples
##Add ploidy barplots underneath the CNA panel

##Add LOH annotations for all mutations by using different symbols for homozygous mutations (due to LOH) and non-homozygous mutations 

mCRC_dir <- Sys.getenv("mCRC_DIR")
setwd(file.path(mCRC_dir,"plot"))

load("../results/mCRC_All_LOH_Mutations.RData")

mycolors <- c("#D55E00","#0072B2","#009E73","#E69F00","#CC79A7","#F0E442",
	      "#999999","#FFFFFF")

dir1 <- "../results/"
samples <- read.delim("../results/alltumorlabels.txt",as.is=TRUE)$Bam
samplelabels <- read.delim("../results/alltumorlabels.txt",as.is=TRUE)$Label

sitetypes <- c("P","LN","LI","LU","BM")

##Cases with primary samples
patients <- c("V46","V402","V514","V559","V750","V824","V855","V930","V953","V974","mCRCTB1","mCRCTB7")
patients2 <- c("V46","V402","V514","V559","V750","V824","V855","V930","V953","V974","TB1","TB7")

idx1 <- which(sapply(strsplit(samples,"_",fixed=TRUE),"[",1) %in% patients)
samples <- samples[idx1]
samplelabels <- samplelabels[idx1]

all_snvs <- vector("list",length(patients))
names(all_snvs) <- patients
all_snvs2 <- all_snvs
for (i in 1:length(patients)) {
  patient1 <- patients[i]
  loh1 <- loh_mutations[[patient1]]

  fn1 <- paste(dir1,patient1,"_MuTectSNV_Indel_Coding_Filtered.txt",sep="")
  patient1 <- sub("mCRCTB","TB",patient1)
  patient1 <- sub("Case","",patient1)  
  count1 <- read.delim(fn1,as.is=TRUE)
  count1 <- count1[count1$CtoT_filter == "Pass",]
  freqs <- count1[,grep("freq",names(count1))[-1]]
  count1 <- count1[apply(freqs,1,max)>=0.1,]
  rownames(count1) <- NULL

  idx_loh <- match(paste(loh1$Chr,loh1$Pos,sep="_"),
		   paste(count1$Chr,count1$Pos,sep="_"))

  normal <- sub("_cover","",grep("cover",names(count1),value=TRUE)[1])
  sns <- sub("_mutect","",grep("mutect",names(count1),value=TRUE))
  sns2 <- samplelabels[match(sns,samples)]

  freqs <- count1[,grep("freq",names(count1))[-1]]
  cover <- count1[,grep("cover",names(count1))[-1]]  
  reads <- round(freqs*cover)
  present <- (freqs >= 0.05 & reads >= 3) | freqs >= 0.1
  
  sites <- sapply(strsplit(sns2,"_",fixed=TRUE),"[",2)
  sites[! sns %in% samples] <- NA
  site_types <- sitetypes[sitetypes %in% sites]
  count2 <- count1[,1:6]

  for (site1 in site_types) {
    idx1 <- which(sites == site1)
    idx2 <- match(sns2[idx1],names(loh1))
    freq1 <-freqs[,idx1,drop=FALSE]
    site_present <- rep("no",nrow(count1))
    site_present[apply(present[,idx1,drop=FALSE],1,any)] <- "yes"
    site_present[apply(present[,idx1,drop=FALSE],1,all) & apply(freq1,1,max) >= 0.1] <- "all"
    count2[[paste(patient1,"_",site1,"_present",sep="")]] <- site_present

    site_loh <- rep("no",nrow(count1))
    for (j in idx_loh) {
      if (is.na(j)) next
      if (sum(present[j,idx1] & 
		grepl("A1",loh1[match(j,idx_loh),idx2]))>=length(idx1)/2)
	site_loh[j] <- "yes"
    }
    count2[[paste(patient1,"_",site1,"_LOH",sep="")]] <- site_loh
  }
  all_snvs[[i]] <- count2[count2$Type != "synonymous",]
  all_snvs2[[i]] <- count2
}

allcounts <- NULL

for (i in 1:length(patients)) {
  
  patient1 <- patients2[i]
  count1 <- all_snvs2[[i]]
  present <- count1[,grep("present",names(count1))] != "no"

  sns <- sub("_present","",grep("present",names(count1),value=TRUE))
  sites <- sapply(strsplit(sns,"_",fixed=TRUE),"[",2)
  
  snv_counts <- matrix(0,nrow=2,ncol=length(sites))
  colnames(snv_counts) <- paste(patient1,sites,sep="_")
  snv_counts[1,] <- sum(apply(present,1,all))
  snv_counts[2,] <- colSums(present) - sum(apply(present,1,all))
  allcounts <- cbind(allcounts,snv_counts)
}

##mCRC genes

short_list <- scan("../results/CRCGenes.txt","c")

mutated_genes <- NULL

for (i in 1:length(patients)) 
{
  count1 <- all_snvs[[i]]
  count1 <- count1[apply(count1[,grep("present",names(count1))]=="all",1,any),]
  mutated_genes <- c(mutated_genes,unique(count1$Gene[count1$Gene %in% short_list]))
}

genes <- short_list[short_list %in% mutated_genes]

nmutation <- matrix(0,nrow=length(genes),ncol=length(patients))
rownames(nmutation) <- genes
colnames(nmutation) <- patients2
nmutation2 <- nmutation

for (i in 1:length(patients)) {
  count1 <- all_snvs[[i]]
  count1 <- count1[apply(count1[,grep("present",names(count1))]=="all",1,any),]
  
  table1 <- table(count1$Gene[count1$Gene %in% genes])
  nmutation[names(table1),i] <- pmax(nmutation[names(table1),i],table1)
  count2 <- data.frame(Gene=count1$Gene,
		       N_sample=apply(count1[,grep("present",names(count1))] != "no",1,sum))
  table2 <- xtabs(count2$N_sample ~ count2$Gene)
  table2 <- table2[names(table2) %in% genes]
  nmutation2[names(table2),i] <- table2
}

genes <- genes[order(rowSums(nmutation>0),rowSums(nmutation),
		     rowSums(nmutation2),decreasing=TRUE)]
genes <- genes[order(match(genes, c("APC","KRAS","TP53","SMAD4","PIK3CA")))]
nmutation <- nmutation[genes,]
genes <- c(rownames(nmutation)[apply(nmutation>0,1,sum)>=2],"BRAF","PTEN","PIK3R1")

nmutation <- matrix(0,nrow=length(genes),ncol=length(patients))
rownames(nmutation) <- genes
colnames(nmutation) <- patients2
nmutation2 <- nmutation

for (i in 1:length(patients)) {
  count1 <- all_snvs[[i]]
  ##count1 <- count1[apply(count1[,grep("present",names(count1))]=="all",1,any),]
  
  table1 <- table(count1$Gene[count1$Gene %in% genes])
  nmutation[names(table1),i] <- pmax(nmutation[names(table1),i],table1)
  count2 <- data.frame(Gene=count1$Gene,
		       N_sample=apply(count1[,grep("present",names(count1))] != "no",1,sum))
  table2 <- xtabs(count2$N_sample ~ count2$Gene)
  table2 <- table2[names(table2) %in% genes]
  nmutation2[names(table2),i] <- table2
}

genes <- genes[order(rowSums(nmutation>0),rowSums(nmutation),
		     rowSums(nmutation2),decreasing=TRUE)]
genes <- genes[order(match(genes, c("APC","KRAS","TP53","SMAD4","PIK3CA")))]
nmutation <- nmutation[genes,]
max_mutation <- apply(nmutation,1,max)

present_code <- function(present) {
  code <- rep(0,nrow(present))
  ##present in all samples
  code[apply(present,1,all)] <- 1 
  return(code)
}

mutation_matrix <- matrix(0,nrow=sum(max_mutation),ncol=ncol(allcounts))
rownames(mutation_matrix) <- rep(genes,max_mutation)
colnames(mutation_matrix) <- colnames(allcounts)

for (i in 1:length(patients)) {
  patient1 <- patients2[i]
  count1 <- all_snvs[[i]]
  count1 <- count1[count1$Gene %in% genes,]
  ##count1 <- count1[apply(count1[,grep("present",names(count1))]=="all",1,any),]

  if (nrow(count1) == 0) next
  sns <- sub("_present","",grep("present",names(count1),value=TRUE))

  present <- count1[,grep("present",names(count1))] != "no"
  loh <- count1[,grep("LOH",names(count1))] == "yes"

  for (j in 1:length(genes)) {
    gene1 <- genes[j]
    if (gene1 %in% count1$Gene) {
      idx <- which(count1$Gene == gene1)
      if (length(idx) > 1)
	idx <- idx[order(rowSums(present[idx,]),decreasing=TRUE)]
      for (k in 1:length(idx)) {
	idx1 <- idx[k]
	idx2 <- which(rownames(mutation_matrix) == gene1)[k]
	for (kk in 1:length(sns)) {
	  mutation_matrix[idx2,sns[kk]] <- present[idx1,kk]+loh[idx1,kk]*20
	}
	if (count1$Type[idx1] %in% c("stopgain",
				     "stopgain_indel","frameshift_indel")) {
	  mutation_matrix[idx2,grep(patient1,colnames(mutation_matrix))] <- mutation_matrix[idx2,grep(patient1,colnames(mutation_matrix))] + 10
	}
	##if (grepl("indel",count1$Type[idx1])) {
	##  mutation_matrix[idx2,grep(patient1,colnames(mutation_matrix))] <- mutation_matrix[idx2,grep(patient1,colnames(mutation_matrix))] + 20
	##}
      }
    }
  }
}


gene_table <- rep(0,length(genes))
names(gene_table) <- genes

for (i in 1:length(patients)) {
  count1 <- all_snvs[[i]]
  count1 <- count1[count1$Gene %in% genes,]
  ##count1 <- count1[apply(count1[,grep("present",names(count1))]=="all",1,any),]
  gene_table[unique(count1$Gene)] <- gene_table[unique(count1$Gene)] + 1
}

mutation_matrix2 <- mutation_matrix[nrow(mutation_matrix):1,]

sns2 <- colnames(mutation_matrix2)
pts2 <- sapply(strsplit(sns2,"_",fixed=TRUE),head,1)

idx2 <- 1
for (i in 2:length(pts2)) {
  if (pts2[i] != pts2[i-1])
    idx2 <- c(idx2,NA)
  idx2 <- c(idx2,i)
}

mutation_matrix2 <- mutation_matrix2[,idx2]
sns2 <- colnames(mutation_matrix2)

mutation_point <- function(x,y,code) {
  if (code < 10 | (code >20 & code < 30)) ##missense
    xys <- cbind(c(x-0.4,x,x+0.4,x-0.4),
		 c(y-0.4,y+0.4,y-0.4,y-0.4))
  if ((code > 10 & code < 20) | code > 30) ##nonsense
    xys <- cbind(c(x-0.4,x-0.4,x+0.4,x+0.4,x-0.4),
		 c(y-0.4,y+0.4,y+0.4,y-0.4,y-0.4))

  if (code > 20) ##LOH
    cols <- mycolors[5] 
  else ##No LOH
    cols <- mycolors[3] 

  code <- code %% 10
  if (code != 0) {
    polygon(xys[1:4,],border=NA,col=cols)
  }
}


library(GenomicRanges)

load("../results/mCRC_All_TitanCNA.RData")

source("../scripts/data_processing/TitanCNA_04_TitanCNA2seg.R")
chr_sizes <- read.table("../ref_genomes/broad.Homo_sapiens_assembly19.fa.sizes",as.is=TRUE)[1:22,]
chr_cumsize <- cumsum(c(0,as.numeric(chr_sizes[-22,2])))

##adjust segmean for overall ploidy and purity

all_cna <- vector("list",length=length(samples))
names(all_cna) <- samplelabels
ploidys <- rep(0,length(samples))
names(ploidys) <- samplelabels

for (i in 1:length(samples)) {
  sn1 <- samples[i]
  cna <- titancna2seg(alltitancnaresults[[sn1]]$results,
		      alltitancnaresults[[sn1]]$convergeParams)
  cna$chrom <- as.integer(as.character(cna$chrom))
  cna$allelicratio <- round(cna$allelicratio,3)
  cna$logcopynumberratio <- round(cna$logcopynumberratio,3)

  ploidy1 <- cna$ploidy[1]

  max_prev <- max(cna$cellularprevalence,na.rm=TRUE)
  ploidys[i] <- (ploidy1-2*(1-max_prev))/max_prev

  purity1 <- cna$normalproportion[1]
  ploidy2 <- ploidy1 * (1 - purity1) + 2 * purity1

  prev1 <- cna$cellularprevalence * (1 - purity1)
  purity2 <- max(prev1,na.rm=TRUE)
  prev1[is.na(prev1)] <- purity2
  seg_mean <- cna$seg.mean + log2(ploidy2/2)
  seg_mean_adj <- log2(pmax(2^-2,1+(2^seg_mean-1)/purity2))
  ##seg_mean_adj[1+(2^seg_mean-1)/prev1 <= 0] <- -2 ##Homozygous deletions

  cna$seg.mean.adj <- round(seg_mean_adj,3)
  all_cna[[i]] <- cna
}
##save(all_cna,file="mCRC_All_TitanCNA_Adj.RData")

sns <- samplelabels
patients <- sapply(strsplit(sns,"_",fixed=TRUE),"[",1)
sites <- sapply(strsplit(sns,"_",fixed=TRUE),"[",2)

all_gr <- NULL
for (i in 1:length(samples)) {
  sn1 <- samplelabels[i]
  print(sn1)
  cna <- all_cna[[i]]

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

all_cna_seg <- data.frame(Chr=seqnames(all_gr),
			  Start=start(all_gr),
			  End=end(all_gr),
			  stringsAsFactors=FALSE)

for (i in 1:length(samples)) {
  sn1 <- samplelabels[i]
  print(sn1)
  cna <- all_cna[[i]]

  cna$cellularprevalence[cna$LOHcall == "HET"] <- 
    max(cna$cellularprevalence,na.rm=TRUE)

  seg.mean <- cna$seg.mean.adj

  cn1 <- rep(NA,length(all_gr))
  gr <- GRanges(seqnames=cna$chrom,
		ranges=IRanges(cna$loc.start,cna$loc.end),
		strand="*")
  idx <- findOverlaps(all_gr,gr)
  cn1[queryHits(idx)] <- seg.mean[subjectHits(idx)]
  
  all_cna_seg <- cbind(all_cna_seg,cn1)
  names(all_cna_seg)[length(all_cna_seg)] <- sn1
}


seg_list <- read.table("CNA_Cytoband.txt")
seg_list_gr <- GRanges(seqnames=seg_list$V1,
		       ranges=IRanges(seg_list$V2,seg_list$V3),
		       strand="*")

weight_matrix <- matrix(0,nrow=nrow(seg_list),ncol=nrow(all_cna_seg))
disjoin_gr <- disjoin(c(seg_list_gr,all_gr))
for (i in 1:nrow(seg_list)) {
  idx1 <- which(overlapsAny(disjoin_gr,seg_list_gr[i]))
  idx2 <- findOverlaps(all_gr,disjoin_gr[idx1])
  weight_matrix[i,queryHits(idx2)] <- width(disjoin_gr[idx1[subjectHits(idx2)]])
}
weight_matrix <- weight_matrix/rowSums(weight_matrix)

seg_list_matrix <- data.frame(Chr=seqnames(seg_list_gr),
			      Start=start(seg_list_gr),
			      End=end(seg_list_gr),
			      Label=seg_list$V4,
			      stringsAsFactors=FALSE)

for (i in 1:length(samples)) {
  sn1 <- samplelabels[i]
  seg_cn <- all_cna_seg[[sn1]]
  seg_cn[is.na(seg_cn)] <- 0
  seg_cn <- weight_matrix %*% seg_cn
  seg_list_matrix <- cbind(seg_list_matrix,cn=seg_cn)
  names(seg_list_matrix)[ncol(seg_list_matrix)] <- sn1
}

##One value for each tumor site
seg_list_table <- NULL
patients <- unique(patients)
ploidys2 <- numeric(0)

for (i in 1:length(patients)) {

  patient1 <- patients[i]
  sns <- grep(patient1,samplelabels,value=TRUE)

  cna_list1 <- seg_list_matrix[sns]
  
  sites <- sapply(strsplit(sns,"_",fixed=TRUE),"[",2)
  site_types <- unique(sites[!is.na(sites)])

  cna_table1 <- matrix(0,nrow=nrow(seg_list_matrix),ncol=length(sns))
  rownames(cna_table1) <- seg_list_matrix$Label
  colnames(cna_table1) <- sns
  cna_table1[,] <- as.matrix(seg_list_matrix[,sns])
  
  cna_table2 <-  matrix(0,nrow=nrow(seg_list_matrix),ncol=length(site_types))
  rownames(cna_table2) <- seg_list_matrix$Label
  colnames(cna_table2) <- paste(patient1,site_types,sep="_")

  for (j in 1:length(site_types)) {
    idx1 <- sites %in% site_types[j]
    cna_table2[,j] <- rowMeans(cna_table1[,idx1,drop=FALSE])
    ploidys2 <- c(ploidys2,mean(ploidys[grep(paste(patient1,site_types[j],sep="_"),names(ploidys))]))
    names(ploidys2)[length(ploidys2)] <- paste(patient1,site_types[j],sep="_")
  }
  seg_list_table <- cbind(seg_list_table,cna_table2)
}

seg_list_table <- seg_list_table[nrow(seg_list_table):1,]

sns2 <- colnames(seg_list_table)
pts2 <- sapply(strsplit(sns2,"_",fixed=TRUE),head,1)

idx2 <- 1
for (i in 2:length(pts2)) {
  if (pts2[i] != pts2[i-1])
    idx2 <- c(idx2,NA)
  idx2 <- c(idx2,i)
}

seg_list_table <- seg_list_table[,idx2]
ploidys2 <- ploidys2[idx2]

sns2 <- colnames(seg_list_table)

cna_point <- function(x,y,code) {
  ##code <- log2(round(2^code*2)/2)
  if (code < 0) {
    cols <- colorRamp(c(mycolors[c(2,8)]))(pmax(1+code/2,0))
    cols <- rgb(red=cols[1,1],green=cols[1,2],blue=cols[1,3],
		maxColorValue=255)
      ##rgb(red=pmax(1+code/2,0),green=pmax(1+code/2,0),blue=1)
  }
  else if (code >= 0) {
    cols <- colorRamp(c(mycolors[c(1,8)]))(pmax(1-code/2,0))
    cols <- rgb(red=cols[1,1],green=cols[1,2],blue=cols[1,3],
		maxColorValue=255)
      ##rgb(green=pmax(1-code/2,0),blue=pmax(1-code/2,0),red=1)
  }

  xys <- cbind(c(x-0.4,x-0.4,x+0.4,x+0.4,x-0.4),
	       c(y-0.4,y+0.4,y+0.4,y-0.4,y-0.4))

  if (code != 0) {
    polygon(xys[1:4,],border=NA,col=cols)
  }
}

col_labels <- rbind(colorRamp(c(mycolors[c(2,8)]))(c(0.5)),
		    colorRamp(c(mycolors[c(1,8)]))(1-log2((3:8)/2)/2))
col_labels <- rgb(red=col_labels[,1],green=col_labels[,2],
		  blue=col_labels[,3],
		  maxColorValue=255)
names(col_labels) <- c(1,3:8)


pdf("mCRC_All_Filtered_SNV_Indel_CNA_LOH_2_Matrix.pdf",width=8,height=7.5)

layout(matrix(c(3,1,4,5,2,6),nrow=3),widths=c(8.5,1.5),heights=c(1.3,7,1.7))

par(mar=c(0.5,5,0.5,0.5))
plot(1,1,ylim=c(1,nrow(mutation_matrix2)+nrow(seg_list_table)),
     xlim=c(1,ncol(mutation_matrix2)),
     xlab="",ylab="",type="n",axes=FALSE)
genes2 <- rownames(mutation_matrix)
genes2[-match(genes,genes2)] <- NA
genes2 <- rev(genes2)
axis(2,at=1:nrow(seg_list_table),labels=rownames(seg_list_table),
     las=2,cex.axis=0.9,tick=FALSE,mgp=c(0,0,0))
axis(2,at=1:nrow(mutation_matrix2)+nrow(seg_list_table),labels=genes2,
     las=2,cex.axis=0.9,tick=FALSE,mgp=c(0,0,0))
##axis(1,at=1:ncol(mutation_matrix2),labels=sns2,
##     las=2,cex.axis=0.7,tick=FALSE)

abline(v=which(is.na(mutation_matrix2[1,])),lty=3)
abline(h=which(rownames(mutation_matrix2)[-1] != rownames(mutation_matrix2)[-nrow(mutation_matrix2)])+0.5+nrow(seg_list_table),lty=3)
abline(h=1:nrow(seg_list_table)+0.5,lty=3)
  
for (i in 1:ncol(seg_list_table)) {
  if (is.na(seg_list_table[1,i])) next
  for (j in 1:nrow(seg_list_table)) {
    if (seg_list_table[j,i] != 0)
      cna_point(i,j,seg_list_table[j,i])
  }
}

for (i in 1:ncol(mutation_matrix2)) {
  if (is.na(mutation_matrix2[1,i])) next
  for (j in 1:nrow(mutation_matrix2)) {
    if (mutation_matrix2[j,i] > 0)
      mutation_point(i,j+nrow(seg_list_table),mutation_matrix2[j,i])
  }
}

par(mar=c(0.5,0,0.5,0.5))
par(mgp=c(-14,-16,-17))
idx3 <- rep(NA,nrow(mutation_matrix))
idx3[match(names(gene_table),rownames(mutation_matrix))] <- 1:length(gene_table)
gene_table2 <- gene_table[idx3]
gene_table2 <- gene_table2[length(gene_table2):1]
gene_table2 <- c(rep(NA,nrow(seg_list_table)),gene_table2)

barplot(gene_table2,width=0.8,space=0.25,axisnames=FALSE,
	ylim=c(0.6,length(gene_table2)-0.4),cex.axis=1,cex.names=1,
	xlab="No of affected\npatients",horiz=TRUE)

legend("bottomleft",
       c("Mutation","Missense","Nonsense","Missense w/ LOH","Nonsense w/ LOH",
	 "Copy number",names(col_labels)),
       pch=c(0,17,15,17,15,0,rep(15,7)),
       col=c("white",mycolors[c(3,3,5,5)],"white",col_labels),
       text.font=c(2,1,1,1,1,2,1,1,1,1,1,1,1),
       border=NA,bty="n",cex=0.9,pt.cex=1.2)


par(mar=c(0,5,0.5,0.5))
par(mgp=c(3,1,0))

allcounts2 <- allcounts[,idx2]

xlim <- c(0.6,ncol(allcounts2)-0.4)
ylim <- c(0,max(colSums(allcounts2)+100,na.rm=TRUE))

barplot(allcounts2,width=0.8,space=0.25,
	legend.text=c("Shared","Not Shared"), col=c("black","lightgrey"),
	names.arg=sns2,xlim=xlim,ylim=ylim,ylab="No of mutations",
	cex.axis=1,cex.names=1,las=2,
	args.legend=list(x="topleft",cex=1,bty="n"),
	axes=TRUE,axisnames=FALSE)


par(mar=c(5,5,0,0.5))
barplot(ploidys2,width=0.8,space=0.25,xlim=xlim,ylim=c(4,0),
	names.arg=sns2,ylab="Ploidy",
	cex.axis=1,cex.names=1,las=2,
	axes=TRUE,axisnames=TRUE)
lines(c(xlim[1]-0.5,xlim[2]+0.5),c(2,2),lty=3,lwd=2)


par(mar=c(0.5,0,0.5,0.5))

plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))

dev.off()


pdf("mCRC_All_Filtered_SNV_Indel_LOH_2_Matrix.pdf",width=8,height=5)

layout(matrix(c(3,1,4,2),nrow=2),widths=c(8.5,1.5),heights=c(2.5,7.5))

par(mar=c(5.0,5,0.5,0.5))
plot(1,1,ylim=c(1,nrow(mutation_matrix2)),
     xlim=c(1,ncol(mutation_matrix2)),
     xlab="",ylab="",type="n",axes=FALSE)
genes2 <- rownames(mutation_matrix)
genes2[-match(genes,genes2)] <- NA
genes2 <- rev(genes2)
axis(2,at=1:nrow(mutation_matrix2),labels=genes2,
     las=2,cex.axis=0.9,tick=FALSE,mgp=c(0.5,0.5,0))
axis(1,at=1:ncol(mutation_matrix2),labels=colnames(mutation_matrix2),
     las=2,cex.axis=0.9,tick=FALSE,mgp=c(0.5,0.5,0))

abline(v=which(is.na(mutation_matrix2[1,])),lty=3)
abline(h=which(rownames(mutation_matrix2)[-1] != rownames(mutation_matrix2)[-nrow(mutation_matrix2)])+0.5,lty=3)
  
for (i in 1:ncol(mutation_matrix2)) {
  if (is.na(mutation_matrix2[1,i])) next
  for (j in 1:nrow(mutation_matrix2)) {
    if (mutation_matrix2[j,i] > 0)
      mutation_point(i,j,mutation_matrix2[j,i])
  }
}

par(mar=c(5.0,0,0.5,0.5))
par(mgp=c(3,1,0))
idx3 <- rep(NA,nrow(mutation_matrix))
idx3[match(names(gene_table),rownames(mutation_matrix))] <- 1:length(gene_table)
gene_table2 <- gene_table[idx3]
gene_table2 <- gene_table2[length(gene_table2):1]

barplot(gene_table2,width=0.8,space=0.25,axisnames=FALSE,xlim=c(0,12),
	ylim=c(0.6,length(gene_table2)-0.4),cex.axis=1,cex.names=1,
	xlab="No of affected\npatients",horiz=TRUE)

par(mar=c(0,5,0.5,0.5))
par(mgp=c(3,1,0))

allcounts2 <- allcounts[,idx2]

xlim <- c(0.6,ncol(allcounts2)-0.4)
ylim <- c(0,max(colSums(allcounts2)+100,na.rm=TRUE))

barplot(allcounts2,width=0.8,space=0.25,
	legend.text=c("Shared","Not Shared"), col=c("black","lightgrey"),
	names.arg=sns2,xlim=xlim,ylim=ylim,ylab="No of mutations",
	cex.axis=1,cex.names=1,las=2,
	args.legend=list(x="topleft",cex=1,bty="n"),
	axes=TRUE,axisnames=FALSE)

par(mar=c(0.5,0,0.5,0.5))

plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
legend("bottomleft",
       c("Missense","Nonsense","Missense w/ LOH","Nonsense w/ LOH"),
       pch=c(17,15,17,15),col=c(mycolors[c(3,3,5,5)]),
       text.font=c(1,1,1,1),border=NA,bty="n",cex=0.8,pt.cex=1.2)

dev.off()



pdf("mCRC_All_Filtered_CNA_Matrix.pdf",width=8,height=4)

layout(matrix(c(3,1,4,2),nrow=2),widths=c(8.5,1.5),heights=c(2.5,7.5))

par(mar=c(5,5.5,0.5,0.5))
plot(1,1,ylim=c(1,nrow(seg_list_table)),
     xlim=c(1,ncol(seg_list_table)),
     xlab="",ylab="",type="n",axes=FALSE)
axis(2,at=1:nrow(seg_list_table),labels=rownames(seg_list_table),
     las=2,cex.axis=0.9,tick=FALSE,mgp=c(0.5,0.5,0))
axis(1,at=1:ncol(seg_list_table),labels=colnames(seg_list_table),
     las=2,cex.axis=0.9,tick=FALSE,mgp=c(0.5,0.5,0))

abline(v=which(is.na(seg_list_table[1,])),lty=3)
abline(h=1:nrow(seg_list_table)+0.5,lty=3)
  
for (i in 1:ncol(seg_list_table)) {
  if (is.na(seg_list_table[1,i])) next
  for (j in 1:nrow(seg_list_table)) {
    if (seg_list_table[j,i] != 0)
      cna_point(i,j,seg_list_table[j,i])
  }
}

par(mar=c(5,0,0.5,0.5))
plot(0,0,type="n",axes=FALSE,xlab="",ylab="")

legend("left",c("Copy number",names(col_labels)),
       pch=c(0,rep(15,7)),col=c("white",col_labels),
       text.font=c(2,1,1,1,1,1,1,1),
       border=NA,bty="n",cex=0.9,pt.cex=1.2)

par(mar=c(0,5.5,0.5,0.5))
par(mgp=c(3,1,0))

barplot(ploidys2,width=0.8,space=0.25,xlim=xlim,ylim=c(0,4),
	axisnames=FALSE,ylab="Ploidy",cex.axis=1,cex.names=1,las=2,axes=TRUE)
lines(c(xlim[1]-0.5,xlim[2]+0.5),c(2,2),lty=3,lwd=2)

dev.off()
