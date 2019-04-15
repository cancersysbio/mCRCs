###Mutations grouped based on sharing patterns
###Compare primary/brain clonal, shared (observed in >= 2) and private (observed in 1) mutations

mCRC_dir <- Sys.getenv("mCRC_DIR")
setwd(file.path(mCRC_dir,"plot"))

dir1 <- "../results/"
samples <- read.delim("../results/alltumorlabels.txt",as.is=TRUE)$Bam
samplelabels <- read.delim("../results/alltumorlabels.txt",as.is=TRUE)$Label

##Cases with primary samples
patients <- c("V46","V402","V514","V559","V750","V824","V855","V930","V953","V974","mCRCTB1","mCRCTB7","UchiCase2")
FFPE <- c("V46","V402","V514","V559","V750","V824","V855","V930","V953","V974")
FF <- c("mCRCTB1","mCRCTB7","UchiCase2")

get_spec <- function(snv) {
  change <- c("CA","CT","CG",
	      "GT","GA","GC",
	      "AC","AG","AT",
	      "TG","TC","TA")
  code <- c(1,2,3,1,2,3,4,5,6,4,5,6)
  recode <- code[match(paste(snv[,1],snv[,2],sep=""),change)]
  return(tabulate(recode,6))
}

spec <- matrix(0,nrow=12,ncol=7)
colnames(spec) <- c("C|G->A|T","C|G->T|A","C|G->G|C",
		    "A|T->C|G","A|T->G|C","A|T->T|A",
		    "Indel")
rownames(spec) <- paste(rep(c("FFPE","FF"),each=6),
			c("Primary C","Primary S","Primary P",
			  "Brain C","Brain S","Brain P"))
rownames(spec) <- sub("FF Brain","FF Liver",rownames(spec))

spec_f <- spec ##FFPE artifacts filtered

for (i in 1:length(patients)) {
  patient1 <- patients[i]
  print(patient1)
  fn1 <- paste(dir1,patient1,"_MuTectSNV_Indel_Coding_CHAT_Filtered.txt",sep="")
  count1 <- read.delim(fn1,as.is=TRUE)
  sns <- sub("_mutect","",grep("mutect",names(count1),value=TRUE))
  if (!all(sns %in% samples)) print(sns[! sns %in% samples])
  sns2 <- samplelabels[match(sns,samples)]
  tissues <- sapply(strsplit(sns2,"_"),"[",2)
  
  mutect <- count1[,grep("mutect",names(count1))]
  freqs <- count1[,grep("freq",names(count1))[-1]][,1:length(sns)]
  cover <- count1[,grep("cover",names(count1))[-1]]  
  reads <- round(cover*freqs)
  present <- (freqs >= 0.05 & reads >= 3) | freqs >= 0.1

  if (length(grep("P",tissues))>=2 && patient1 %in% FFPE) {
    idx <- grep("P",tissues)
    idx1 <- count1$Primary_average_CCF_adj >= 0.6
    idx2 <- count1$Primary_average_CCF_adj < 0.6 & rowSums(present[,idx]) > 1
    idx3 <- count1$Primary_average_CCF_adj < 0.6 & 
      rowSums(present[,idx]) == 1 &
	rowSums(present[,-idx]) == 0
    
    spec[1,1:6] <- spec[1,1:6] + get_spec(count1[idx1,3:4])
    spec[1,7] <- spec[1,7] + length(grep("indel",count1$Type[idx1]))
    spec[2,1:6] <- spec[2,1:6] + get_spec(count1[idx2,3:4])
    spec[2,7] <- spec[2,7] + length(grep("indel",count1$Type[idx2]))
    spec[3,1:6] <- spec[3,1:6] + get_spec(count1[idx3,3:4])
    spec[3,7] <- spec[3,7] + length(grep("indel",count1$Type[idx3]))
  }

  if (length(grep("P",tissues))>=2 && patient1 %in% FF) {
    idx <- grep("P",tissues)
    idx1 <- count1$Primary_average_CCF_adj >= 0.6
    idx2 <- count1$Primary_average_CCF_adj < 0.6 & rowSums(present[,idx]) > 1
    idx3 <- count1$Primary_average_CCF_adj < 0.6 & 
      rowSums(present[,idx]) == 1 &
	rowSums(present[,-idx]) == 0
    
    spec[7,1:6] <- spec[7,1:6] + get_spec(count1[idx1,3:4])
    spec[7,7] <- spec[7,7] + length(grep("indel",count1$Type[idx1]))
    spec[8,1:6] <- spec[8,1:6] + get_spec(count1[idx2,3:4])
    spec[8,7] <- spec[8,7] + length(grep("indel",count1$Type[idx2]))
    spec[9,1:6] <- spec[9,1:6] + get_spec(count1[idx3,3:4])
    spec[9,7] <- spec[9,7] + length(grep("indel",count1$Type[idx3]))
  }

  if (length(grep("BM",tissues))>=2 && patient1 %in% FFPE) {
    idx <- grep("BM",tissues)
    idx1 <- count1$BrainMet_average_CCF_adj >= 0.6
    idx2 <- count1$BrainMet_average_CCF_adj < 0.6 & rowSums(present[,idx]) > 1
    idx3 <- count1$BrainMet_average_CCF_adj < 0.6 & 
      rowSums(present[,idx]) == 1 &
	rowSums(present[,-idx]) == 0

    spec[4,1:6] <- spec[4,1:6] + get_spec(count1[idx1,3:4])
    spec[4,7] <- spec[4,7] + length(grep("indel",count1$Type[idx1]))
    spec[5,1:6] <- spec[5,1:6] + get_spec(count1[idx2,3:4])
    spec[5,7] <- spec[5,7] + length(grep("indel",count1$Type[idx2]))
    spec[6,1:6] <- spec[6,1:6] + get_spec(count1[idx3,3:4])
    spec[6,7] <- spec[6,7] + length(grep("indel",count1$Type[idx3]))
  }

  if (length(grep("LI",tissues))>=2 && patient1 %in% FF) {
    idx <- grep("LI",tissues)
    idx1 <- count1$Liver_average_CCF_adj >= 0.6
    idx2 <- count1$Liver_average_CCF_adj < 0.6 & rowSums(present[,idx]) > 1
    idx3 <- count1$Liver_average_CCF_adj < 0.6 & 
      rowSums(present[,idx]) == 1 &
	rowSums(present[,-idx]) == 0

    spec[10,1:6] <- spec[10,1:6] + get_spec(count1[idx1,3:4])
    spec[10,7] <- spec[10,7] + length(grep("indel",count1$Type[idx1]))
    spec[11,1:6] <- spec[11,1:6] + get_spec(count1[idx2,3:4])
    spec[11,7] <- spec[11,7] + length(grep("indel",count1$Type[idx2]))
    spec[12,1:6] <- spec[12,1:6] + get_spec(count1[idx3,3:4])
    spec[12,7] <- spec[12,7] + length(grep("indel",count1$Type[idx3]))
  }
}


for (i in 1:length(patients)) {
  patient1 <- patients[i]
  print(patient1)
  fn1 <- paste(dir1,patient1,"_MuTectSNV_Indel_Coding_CHAT_Filtered.txt",sep="")
  count1 <- read.delim(fn1,as.is=TRUE)
  count1 <- count1[count1$CtoT_filter == "Pass",]
  sns <- sub("_mutect","",grep("mutect",names(count1),value=TRUE))
  if (!all(sns %in% samples)) print(sns[! sns %in% samples])
  sns2 <- samplelabels[match(sns,samples)]
  tissues <- sapply(strsplit(sns2,"_"),"[",2)
  
  mutect <- count1[,grep("mutect",names(count1))]
  freqs <- count1[,grep("freq",names(count1))[-1]][,1:length(sns)]
  cover <- count1[,grep("cover",names(count1))[-1]]  
  reads <- round(cover*freqs)
  present <- (freqs >= 0.05 & reads >= 3) | freqs >= 0.1

  if (length(grep("P",tissues))>=2 && patient1 %in% FFPE) {
    idx <- grep("P",tissues)
    idx1 <- count1$Primary_average_CCF_adj >= 0.6
    idx2 <- count1$Primary_average_CCF_adj < 0.6 & rowSums(present[,idx]) > 1
    idx3 <- count1$Primary_average_CCF_adj < 0.6 & 
      rowSums(present[,idx]) == 1 &
	rowSums(present[,-idx]) == 0
    
    spec_f[1,1:6] <- spec_f[1,1:6] + get_spec(count1[idx1,3:4])
    spec_f[1,7] <- spec_f[1,7] + length(grep("indel",count1$Type[idx1]))
    spec_f[2,1:6] <- spec_f[2,1:6] + get_spec(count1[idx2,3:4])
    spec_f[2,7] <- spec_f[2,7] + length(grep("indel",count1$Type[idx2]))
    spec_f[3,1:6] <- spec_f[3,1:6] + get_spec(count1[idx3,3:4])
    spec_f[3,7] <- spec_f[3,7] + length(grep("indel",count1$Type[idx3]))
  }

  if (length(grep("P",tissues))>=2 && patient1 %in% FF) {
    idx <- grep("P",tissues)
    idx1 <- count1$Primary_average_CCF_adj >= 0.6
    idx2 <- count1$Primary_average_CCF_adj < 0.6 & rowSums(present[,idx]) > 1
    idx3 <- count1$Primary_average_CCF_adj < 0.6 & 
      rowSums(present[,idx]) == 1 &
	rowSums(present[,-idx]) == 0
    
    spec_f[7,1:6] <- spec_f[7,1:6] + get_spec(count1[idx1,3:4])
    spec_f[7,7] <- spec_f[7,7] + length(grep("indel",count1$Type[idx1]))
    spec_f[8,1:6] <- spec_f[8,1:6] + get_spec(count1[idx2,3:4])
    spec_f[8,7] <- spec_f[8,7] + length(grep("indel",count1$Type[idx2]))
    spec_f[9,1:6] <- spec_f[9,1:6] + get_spec(count1[idx3,3:4])
    spec_f[9,7] <- spec_f[9,7] + length(grep("indel",count1$Type[idx3]))
  }

  if (length(grep("BM",tissues))>=2 && patient1 %in% FFPE) {
    idx <- grep("BM",tissues)
    idx1 <- count1$BrainMet_average_CCF_adj >= 0.6
    idx2 <- count1$BrainMet_average_CCF_adj < 0.6 & rowSums(present[,idx]) > 1
    idx3 <- count1$BrainMet_average_CCF_adj < 0.6 & 
      rowSums(present[,idx]) == 1 &
	rowSums(present[,-idx]) == 0

    spec_f[4,1:6] <- spec_f[4,1:6] + get_spec(count1[idx1,3:4])
    spec_f[4,7] <- spec_f[4,7] + length(grep("indel",count1$Type[idx1]))
    spec_f[5,1:6] <- spec_f[5,1:6] + get_spec(count1[idx2,3:4])
    spec_f[5,7] <- spec_f[5,7] + length(grep("indel",count1$Type[idx2]))
    spec_f[6,1:6] <- spec_f[6,1:6] + get_spec(count1[idx3,3:4])
    spec_f[6,7] <- spec_f[6,7] + length(grep("indel",count1$Type[idx3]))
  }

  if (length(grep("LI",tissues))>=2 && patient1 %in% FF) {
    idx <- grep("LI",tissues)
    idx1 <- count1$Liver_average_CCF_adj >= 0.6
    idx2 <- count1$Liver_average_CCF_adj < 0.6 & rowSums(present[,idx]) > 1
    idx3 <- count1$Liver_average_CCF_adj < 0.6 & 
      rowSums(present[,idx]) == 1 &
	rowSums(present[,-idx]) == 0

    spec_f[10,1:6] <- spec_f[10,1:6] + get_spec(count1[idx1,3:4])
    spec_f[10,7] <- spec_f[10,7] + length(grep("indel",count1$Type[idx1]))
    spec_f[11,1:6] <- spec_f[11,1:6] + get_spec(count1[idx2,3:4])
    spec_f[11,7] <- spec_f[11,7] + length(grep("indel",count1$Type[idx2]))
    spec_f[12,1:6] <- spec_f[12,1:6] + get_spec(count1[idx3,3:4])
    spec_f[12,7] <- spec_f[12,7] + length(grep("indel",count1$Type[idx3]))
  }
}

pdf("mCRC_All_SNV_Indel_Group_Spectrum.pdf",width=6,height=6)
par(mar=c(8,4,3,1))

spec2 <- spec/rowSums(spec)
idx <- c("C|G->T|A","C|G->A|T","C|G->G|C","A|T->G|C","A|T->C|G","A|T->T|A","Indel")
spec2 <- spec2[,idx]
spec3 <- spec2[,-7] / rowSums(spec2[,-7])

barplot(t(spec2),width=0.8,space=0.25,col=2:8,xlim=c(0,nrow(spec2)+4),
	legend.text=TRUE,las=2,cex.names=1,
	ylab="Proportion",main="Coding SNVs and indels",
	args.legend=list("toprigh",cex=0.8))
barplot(t(spec3),width=0.8,space=0.25,col=2:7,xlim=c(0,nrow(spec2)+4),
	legend.text=TRUE,las=2,cex.names=1,
	ylab="Proportion",main="Coding SNVs",
	args.legend=list("toprigh",cex=0.8))

dev.off()

pdf("mCRC_All_Filtered_SNV_Indel_Group_Spectrum.pdf",width=6,height=6)
par(mar=c(8,4,3,1))

spec2 <- spec_f/rowSums(spec_f)
idx <- c("C|G->T|A","C|G->A|T","C|G->G|C","A|T->G|C","A|T->C|G","A|T->T|A","Indel")
spec2 <- spec2[,idx]
spec3 <- spec2[,-7] / rowSums(spec2[,-7])

barplot(t(spec2),width=0.8,space=0.25,col=2:8,xlim=c(0,nrow(spec2)+4),
	legend.text=TRUE,las=2,cex.names=1,
	ylab="Proportion",main="Coding SNVs and indels",
	args.legend=list("toprigh",cex=0.8))
barplot(t(spec3),width=0.8,space=0.25,col=2:7,xlim=c(0,nrow(spec2)+4),
	legend.text=TRUE,las=2,cex.names=1,
	ylab="Proportion",main="Coding SNVs",
	args.legend=list("toprigh",cex=0.8))

dev.off()
