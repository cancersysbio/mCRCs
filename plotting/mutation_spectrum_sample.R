###Sample level mutational spectrum

mCRC_dir <- Sys.getenv("mCRC_DIR")
setwd(file.path(mCRC_dir,"plot"))

dir1 <- "../results/"
samples <- read.delim("../results/alltumorlabels.txt",as.is=TRUE)$Bam
samplelabels <- read.delim("../results/alltumorlabels.txt",as.is=TRUE)$Label

mycolors <- c("#D55E00","#0072B2","#009E73","#E69F00","#CC79A7","#F0E442",
	      "#999999","#FFFFFF")

##Cases with primary samples
patients <- c("V46","V402","V514","V559","V750","V824","V855","V930","V953","V974","mCRCTB1","mCRCTB7")

idx <- sapply(strsplit(samples,"_",fixed=TRUE),"[",1) %in% patients
samples <- samples[idx]
samplelabels <- samplelabels[idx]

get_spec <- function(snv) {
  change <- c("CA","CT","CG",
	      "GT","GA","GC",
	      "AC","AG","AT",
	      "TG","TC","TA")
  code <- c(1,2,3,1,2,3,4,5,6,4,5,6)
  recode <- code[match(paste(snv[,1],snv[,2],sep=""),change)]
  return(tabulate(recode,6))
}

spec <- matrix(0,nrow=length(samples),ncol=7)
colnames(spec) <- c("C|G->A|T","C|G->T|A","C|G->G|C",
		    "A|T->C|G","A|T->G|C","A|T->T|A","Indel")
rownames(spec) <- samples

for (i in 1:length(patients)) {
  patient1 <- patients[i]
  fn1 <- paste(dir1,patient1,"_MuTectSNV_Indel_Coding_Filtered.txt",sep="")
  count1 <- read.delim(fn1,as.is=TRUE)
  count1 <- count1[count1$CtoT_filter == "Pass",]
  sns <- sub("_mutect","",grep("mutect",names(count1),value=TRUE))
  sns1 <- sns[sns %in% samples]
  sns2 <- samplelabels[match(sns1,samples)]

  indel <- grepl("indel",count1$Type)
  snv <- !indel

  mutect <- count1[,grep("mutect",names(count1))]
  freqs <- count1[,grep("freq",names(count1))[-1]]
  cover <- count1[,grep("cover",names(count1))[-1]]  
  reads <- round(cover*freqs)
  present <- (freqs >= 0.05 & reads >= 3) | freqs >= 0.1
  colnames(present) <- sub("_freq","",colnames(present))
  present <- present[,sns1]

  ##Sample level
  for (j in 1:length(sns1)) {
    sn1 <- sns1[j]
    spec[sn1,1:6] <- get_spec(count1[present[,sn1],3:4])
    spec[sn1,7] <- sum(present[indel,sn1])
  }
}

spec2 <- spec/rowSums(spec)
idx <- c("C|G->T|A","C|G->A|T","C|G->G|C","A|T->G|C","A|T->C|G","A|T->T|A","Indel")
spec2 <- spec2[,idx]

patients <- sapply(strsplit(samples,"_",fixed=TRUE),head,1)
spaces <- rep(0.25,length(samples))
spaces[which(patients[-1] != patients[-length(patients)])+1] <- 1.5
spaces[1] <- 0.25
#spaces <- c(1/9,11/9,spaces)

pdf("mCRC_All_Filtered_SNV_Indel_Spectrum.pdf",width=14,height=5)
par(mar=c(5,4,1,1))
sns2 <- samplelabels
a <- barplot(t(spec2),width=0.8,space=spaces,col=mycolors,
	     legend.text=TRUE,names.arg=sns2,
	     las=2,cex.names=0.7,xlim=c(1,nrow(spec2)+sum(spaces>1)+6),
	     ylab="Proportion",##main="Coding SNV",
	     args.legend=list("toprigh",cex=0.8,bty="n"),ylim=c(0,1))
##text(a+0.3,1.01,rowSums(spec),pos=3,srt=90,cex=0.7,offset=0.7)
dev.off()
