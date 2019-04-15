###Mutations grouped based on treatment
###Compare primary/brain clonal, shared (observed in >= 2) and private (observed in 1) mutations

mCRC_dir <- Sys.getenv("mCRC_DIR")
setwd(file.path(mCRC_dir,"plot"))

dir1 <- "../results/"
samples <- read.delim("../results/alltumorlabels.txt",as.is=TRUE)$Bam
samplelabels <- read.delim("../results/alltumorlabels.txt",as.is=TRUE)$Label

mycolors <- c("#D55E00","#0072B2","#009E73","#E69F00","#CC79A7","#F0E442",
	      "#999999","#FFFFFF")

##Cases with primary samples
patients <- c("V46","V402","V514","V559","V750","V824","V855","V930","V953","V974")
FFPE <- c("V46","V402","V514","V559","V750","V824","V855","V930","V953","V974")
FF <- c("mCRCTB1","mCRCTB7","UchiCase2")

chemo <- c("V46","V514","V559","V750","V855","V930")
nochemo <- c("V402","V824","V953","V974")

get_spec <- function(snv) {
  change <- c("CA","CT","CG",
	      "GT","GA","GC",
	      "AC","AG","AT",
	      "TG","TC","TA")
  code <- c(1,2,3,1,2,3,4,5,6,4,5,6)
  recode <- code[match(paste(snv[,1],snv[,2],sep=""),change)]
  return(tabulate(recode,6))
}

spec <- matrix(0,nrow=6,ncol=7)
colnames(spec) <- c("C|G->A|T","C|G->T|A","C|G->G|C",
		    "A|T->C|G","A|T->G|C","A|T->T|A",
		    "Indel")
rownames(spec) <- paste(rep(c("Chemo","No chemo"),each=3),c("S","P","BM"))

count <- matrix(0,nrow=6,ncol=2)
rownames(count) <- rownames(spec)
colnames(count) <- c("SNV","Indel")

for (i in 1:length(patients)) {
  patient1 <- patients[i]
  print(patient1)
  fn1 <- paste(dir1,patient1,"_MuTectSNV_Indel_Coding_CHAT_Filtered.txt",sep="")
  count1 <- read.delim(fn1,as.is=TRUE)
  count1 <- count1[count1$CtoT_filter == "Pass",]
  sns <- sub("_mutect","",grep("mutect",names(count1),value=TRUE))
  ##if (!all(sns %in% samples)) print(sns[! sns %in% samples])
  sns1 <- sns[sns %in% samples]
  sns2 <- samplelabels[match(sns1,samples)]
  tissues <- sapply(strsplit(sns2,"_"),"[",2)

  indel <- grepl("indel",count1$Type)
  snv <- !indel
  
  mutect <- count1[,grep("mutect",names(count1))]
  freqs <- count1[,grep("freq",names(count1))[-1]][,1:length(sns)]
  cover <- count1[,grep("cover",names(count1))[-1]]  
  reads <- round(cover*freqs)
  present <- (freqs >= 0.05 & reads >= 3) | freqs >= 0.1
  colnames(present) <- sub("_freq","",colnames(present))
  present <- present[,sns1]

  idx_p <- grep("P",tissues)
  idx_b<- grep("BM",tissues)

  idx1 <- apply(present[,idx_p,drop=FALSE],1,any) & 
    apply(present[,idx_b,drop=FALSE],1,any)
  idx2 <- apply(present[,idx_p,drop=FALSE],1,any) & 
    (! apply(present[,idx_b,drop=FALSE],1,any))
  idx3 <- (! apply(present[,idx_p,drop=FALSE],1,any)) & 
    apply(present[,idx_b,drop=FALSE],1,any)

  if (patient1 %in% chemo) {
    spec[1,1:6] <- spec[1,1:6] + get_spec(count1[idx1,3:4])
    spec[1,7] <- spec[1,7] + sum(idx1 & indel)
    spec[2,1:6] <- spec[2,1:6] + get_spec(count1[idx2,3:4])
    spec[2,7] <- spec[2,7] + sum(idx2 & indel)
    spec[3,1:6] <- spec[3,1:6] + get_spec(count1[idx3,3:4])
    spec[3,7] <- spec[3,7] + sum(idx3 & indel)

    count[1,1] <- count[1,1] + sum(idx1 & snv)
    count[1,2] <- count[1,2] + sum(idx1 & indel)
    count[2,1] <- count[2,1] + mean(colSums(present[(!idx1) & snv,idx_p,drop=FALSE]))
    count[2,2] <- count[2,2] + mean(colSums(present[(!idx1) & indel,idx_p,drop=FALSE]))
    count[3,1] <- count[3,1] + mean(colSums(present[(!idx1) & snv,idx_b,drop=FALSE]))
    count[3,2] <- count[3,2] + mean(colSums(present[(!idx1) & indel,idx_b,drop=FALSE]))    
  }
  
  if (patient1 %in% nochemo) {
    spec[4,1:6] <- spec[4,1:6] + get_spec(count1[idx1,3:4])
    spec[4,7] <- spec[4,7] + sum(idx1 & indel)
    spec[5,1:6] <- spec[5,1:6] + get_spec(count1[idx2,3:4])
    spec[5,7] <- spec[5,7] + sum(idx2 & indel)
    spec[6,1:6] <- spec[6,1:6] + get_spec(count1[idx3,3:4])
    spec[6,7] <- spec[6,7] + sum(idx3 & indel)

    count[4,1] <- count[4,1] + sum(idx1 & snv)
    count[4,2] <- count[4,2] + sum(idx1 & indel)
    count[5,1] <- count[5,1] + mean(colSums(present[(!idx1) & snv,idx_p,drop=FALSE]))
    count[5,2] <- count[5,2] + mean(colSums(present[(!idx1) & indel,idx_p,drop=FALSE]))
    count[6,1] <- count[6,1] + mean(colSums(present[(!idx1) & snv,idx_b,drop=FALSE]))
    count[6,2] <- count[6,2] + mean(colSums(present[(!idx1) & indel,idx_b,drop=FALSE]))    
  }
}

count[1:3,] <- count[1:3,]/length(chemo)
count[4:6,] <- count[4:6,]/length(nochemo)

pdf("mCRC_All_Filtered_SNV_Indel_Group_Spectrum.pdf",width=6,height=6)
par(mfrow=c(2,1))

par(mar=c(0.5,4,3,1))

spec2 <- spec/rowSums(spec)
idx <- c("C|G->T|A","C|G->A|T","C|G->G|C","A|T->G|C","A|T->C|G","A|T->T|A","Indel")
spec2 <- spec2[,idx]
spec3 <- spec2[,-7] / rowSums(spec2[,-7])

spaces <- c(0.25,0.25,0.25,1.5,0.25,0.25)

barplot(t(spec2),width=0.8,space=spaces,col=mycolors,
	xlim=c(0,(nrow(spec2)+1)*1.4),axisnames=FALSE,
	legend.text=TRUE,las=1,cex.names=1,
	ylab="Proportion",main="Coding SNVs and indels",
	args.legend=list("toprigh",cex=1,bty="n"))

par(mar=c(4,4,0.5,1))

a <- barplot(t(count),width=0.8,space=spaces,
	     xlim=c(0,(nrow(spec2)+1)*1.4),ylim=c(max(count),0),
	     names.arg=rep(c("S","P","BM"),2),
	     legend.text=TRUE,las=1,cex.names=1,
	     ylab="Average count",main="",
	     args.legend=list("toprigh",cex=1,bty="n"))

axis(1,at=a[c(2,5)],c("Chemo","No chemo"),tick=FALSE,line=2)


par(mar=c(0.5,4,3,1))

a <- barplot(t(spec3),width=0.8,space=spaces,col=mycolors,
	     xlim=c(0,(nrow(spec2)+1)*1.4),axisnames=FALSE,
	     legend.text=TRUE,las=1,cex.names=1,
	     ylab="Proportion",main="Coding SNVs",
	     args.legend=list("toprigh",cex=1,bty="n"))

par(mar=c(4,4,0.5,1))

a <- barplot(count[,1],width=0.8,space=spaces,
	     xlim=c(0,(nrow(spec2)+1)*1.4),ylim=c(max(count),0),
	     names.arg=rep(c("S","P","BM"),2),
	     legend.text=FALSE,las=1,cex.names=1,
	     ylab="Average count",main="",
	     args.legend=list("toprigh",cex=1,bty="n"))

axis(1,at=a[c(2,5)],c("Chemo","No chemo"),tick=FALSE,line=2)

dev.off()
