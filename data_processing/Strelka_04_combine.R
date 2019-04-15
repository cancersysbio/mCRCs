#!/srv/gsfs0/software/R/3.2.2/bin/Rscript

##Combine snv and indel calls

inputpar <- commandArgs(TRUE)
if (length(inputpar) < 1 ) stop("Wrong number of input arguments: <patient_name>")

patient <- inputpar[1]
print(patient)

snv_file <- paste0(patient,"_MuTectSNV_Coding.txt")
indel_file <- paste0(patient,"_Strelka_Indel_Coding.txt")

snv <- read.delim(snv_file,as.is=TRUE)
indel <- read.delim(indel_file,as.is=TRUE)
  
##only keep start positions for indels and change names in indel to match names in snv
indel$End <- NULL
names(indel)[2] <- "Pos"
names(indel) <- sub("strelka","mutect",names(indel))

snv2 <- rbind(snv,indel)
snv2 <- snv2[order(match(snv2$Chr,c(1:22,"X","Y")),snv2$Pos),]

write.table(snv2,paste0(patient,"_MuTectSNV_Indel_Coding.txt"),
	    quote=FALSE,row.names=FALSE,sep="\t")
