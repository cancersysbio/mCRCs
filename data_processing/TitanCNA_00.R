#!/srv/gsfs0/software/R/R-3.1.1/bin/Rscript

##Use multiple tumor samples to filter out unreliable germline HET SNPs

library(TitanCNA)

inputpar <- commandArgs(TRUE)
if (length(inputpar) < 2) stop("Wrong number of input arguments: <normal.bam> <tumor.bam.list>")

normalbam <- inputpar[1]
tumorbams <- scan(inputpar[2],"c")

sn1 <- strsplit(basename(normalbam),".",fixed=TRUE)[[1]][1]

print(sn1)
print(bamfile)

normalvcf <- paste(sn1,".titancna.vcf",sep="")

fn <- paste(sn1,".alleleCounts.txt",sep="")
tmp <- extractAlleleReadCounts(bamFile=bamfile, 
			       bamIndex=paste(bamfile,".bai",sep=""),
			       positions=normalvcf, outputFilename = fn, 
			       pileupParam = PileupParam(max_depth=2000))
normal <- read.delim(fn,as.is=TRUE)

for (bam1 in tumorbams) {
  print(bam1)

  sn1 <- strsplit(basename(bam1),".",fixed=TRUE)[[1]][1]
  fn <- paste(sn1,".alleleCounts.txt",sep="")
  tmp <- extractAlleleReadCounts(bamFile=bam1, 
				 bamIndex=paste(bam1,".bai",sep=""),
				 positions=normalvcf, outputFilename = fn, 
				 pileupParam = PileupParam(max_depth=2000))
}

tumor <- NULL  
for (bam1 in tumorbams) {
  sn1 <- strsplit(basename(bam1),".",fixed=TRUE)[[1]][1]
  fn <- paste(sn1,".alleleCounts.txt",sep="")
  if (is.null(tumor)) {
    tumor <- read.delim(fn,as.is=TRUE)
    if (! all(tumor$chr == normal$chr & tumor$position == normal$position &
		tumor$ref == normal$ref & tumor$Nref == normal$Nref))
      stop("Not all equal")      
    }
  else {
    tumor2 <- read.delim(fn,as.is=TRUE)
    if (! all(tumor$chr == tumor2$chr & tumor$position == tumor2$position &
		tumor$ref == tumor2$ref & tumor$Nref == tumor2$Nref))
      stop("Not all equal")
    tumor$refCount <- tumor$refCount + tumor2$refCount
    tumor$NrefCount <- tumor$NrefCount + tumor2$NrefCount
  }
}

idx <- which((tumor$refCount/(tumor$refCount + tumor$NrefCount) >= 0.1 & 
	  tumor$NrefCount/(tumor$refCount + tumor$NrefCount) >= 0.1) | 
	    (tumor$refCount + tumor$NrefCount > 0 &
	       (normal$refCount + normal$NrefCount >= 10 & 
		  normal$refCount/(normal$refCount + normal$NrefCount) >= 0.3 &
		    normal$NrefCount/(normal$refCount + normal$NrefCount) >= 0.3)))

for (bam1 in tumorbams) {
  sn1 <- strsplit(basename(bam1),".",fixed=TRUE)[[1]][1]
  fn <- paste(sn1,".alleleCounts.txt",sep="")
  tumor <- read.delim(fn,as.is=TRUE)
  tumor <- tumor[idx,]
  write.table(tumor,file=fn,sep="\t",quote=FALSE,row.names=FALSE)
}
