#!/srv/gsfs0/software/R/R-3.1.1/bin/Rscript

##Allow two normal bam files, separated by "|", the first from matched normal sample used for calling germline heterozygous SNVs, the second from possibly pooled normal samples used for calculating depth ratios at those SNV loci

##Only use chromosome 1:22 since only het snp loci in normal are used 
##max copy number = 8
##no of clones: 1-3
##use TitanCNA 1.5.7
##using input ploidy as initial value

inputpar <- commandArgs(TRUE)
if (length(inputpar) < 3 ) stop("Wrong number of input arguments: <ploidy> <normal.bam> <tumor.bam> <target.bed>")

input_ploidy <- as.numeric(inputpar[1])
normalbam <- inputpar[2]
tumorbam <- inputpar[3]
targetbed <- NULL
if (length(inputpar) == 4) {
  if (!file.exists(inputpar[4]))
    stop("File ",inputpar[4]," doesn't exist")  
  targetbed <- inputpar[4]
}

normalbam1 <- strsplit(normalbam,"|",fixed=TRUE)[[1]][1]
normalbam2 <- strsplit(normalbam,"|",fixed=TRUE)[[1]][2]
if (is.na(normalbam2))
  normalbam2 <- normalbam1

if (!file.exists(tumorbam))
    stop("File ",tumorbam," doesn't exist")

library(TitanCNA)

chrs <- 1:22

targetseg <- NULL
if (!is.null(targetbed)) {
  targetseg <- read.table(targetbed,,as.is=TRUE)[,1:3]
}

sn <- strsplit(basename(tumorbam),".",fixed=TRUE)[[1]][1]
sn0 <- strsplit(basename(normalbam1),".",fixed=TRUE)[[1]][1]
normalvcf <- paste(sn0,".titancna.vcf",sep="")

fn <- paste(sn,".alleleCounts.txt",sep="")
if ( ! file.exists(fn)) {
  tmp <- extractAlleleReadCounts(bamFile=tumorbam, 
				 bamIndex=paste(tumorbam,".bai",sep=""),
				 positions=normalvcf, outputFilename = fn, 
				 pileupParam = PileupParam(max_depth=1000,min_nucleotide_depth=1))
}
snpData <- loadAlleleCounts(fn)

tumWig <- paste(sn,".readcount.wig",sep="")
sn02 <- strsplit(basename(normalbam2),".",fixed=TRUE)[[1]][1]
normWig <- paste(sn02,".readcount.wig",sep="")
gcWig <- "ref_genomes/hg19.gc.wig"
mapWig <- "ref_genomes/hg19.map.wig"
cnData <- correctReadDepth(tumWig,normWig,gcWig,mapWig,targetedSequence=targetseg)
cnData <- cnData[!is.na(cnData$logR),]

logR <- getPositionOverlap(snpData$chr, snpData$posn, cnData)
snpData$logR <- log(2^logR)  #transform the log ratio to natural logs
snpData <- filterData(snpData, chrs = chrs, minDepth = 10, maxDepth = 1000,
		      positionList = NULL)

titancnaresults <- vector("list",9)
n0 <- rep(c(0.05,0.3,0.5),3)
##n0 <- rep(c(0.3,0.5,0.7),3)
nc <- rep(1:3,each=3)

for (j in 1:length(titancnaresults)) { 
  numClusters <- nc[j]
  params <- loadDefaultParameters(copyNumber = 8, numberClonalClusters = numClusters, data=snpData)
  if (!is.null(targetseg)) {
    K <- length(params$genotypeParams$alphaKHyper)
    params$genotypeParams$alphaKHyper <- rep(1000, K)
  }
  params$normalParams$n_0 <- n0[j]
  params$ploidyParams$phi_0 <- input_ploidy
  
  ##allele frequencies for balanced states (allow higher AFs?)
  params$genotypeParams$rt[c(1,4,9,16,25)] <- min(0.6, params$genotypeParams$rt[4])
  
  convergeParams <- runEMclonalCN(snpData, gParams = params$genotypeParams,
				  nParams = params$normalParams,
				  pParams = params$ploidyParams,
				  sParams = params$cellPrevParams,
				  maxiter = 20, maxiterUpdate = 1500,
				  useOutlierState = FALSE, txnExpLen = 1e15,
				  txnZstrength = 5e5,
				  normalEstimateMethod = "map",
				  estimateS = TRUE, estimatePloidy = TRUE)
  optimalPath <- viterbiClonalCN(snpData, convergeParams)
  if (length(unique(optimalPath)) == 1) next
  results <- outputTitanResults(snpData, convergeParams, optimalPath,
				filename = NULL, posteriorProbs = FALSE,
				subcloneProfiles = TRUE)
  ploidy <- tail(convergeParams$phi, 1)
  
  titancnaresults[[j]] <- list(S_DbwIndex=computeSDbwIndex(results)$S_DbwIndex,
			       results=results,convergeParams=convergeParams)
}
save(titancnaresults,file=paste(sn,".TitanCNA.RData",sep=""))

pdf(paste(sn,".TitanCNA.pdf",sep=""),height=5,width=10)

for (j in 1:length(titancnaresults)) {
  if (is.null(titancnaresults[[j]])) next
  
  SD <- round(titancnaresults[[j]]$S_DbwIndex,3)
  convergeParams <- titancnaresults[[j]]$convergeParams
  results <- titancnaresults[[j]]$results
  nclones <- nrow(convergeParams$s)
  ploidy <- round(tail(convergeParams$phi, 1),2)
  meandepth <- round(mean(as.numeric(results$Depth)),2)
  npoints <- nrow(results)
  n <- round(tail(convergeParams$n,1),2)
  s <- round(convergeParams$s[1,ncol(convergeParams$s)],2)
  ploidy2 <- ploidy * (1 - n) + 2 * n
  
  layout(matrix(c(1,2,3,3),nrow=2))
  par(pty="m")
  par(mar=c(4,4,2,1))
  plotCNlogRByChr(results, chr = NULL, ploidy = ploidy2, 
		  ylim = c(-2, 2), cex = 0.25, main = sn,
		  xlab=paste(" SD=",SD," n=",n," s=",s," nc=",nclones," pl=",ploidy," np=",npoints," md=",meandepth,sep=""))

  par(mar=c(4,4,2,1))
  plotAllelicRatio(results, chr = NULL, ylim = c(0, 1), cex = 0.25,
		   xlab = "", main = "")
  
  allstate <- paste(results$Chr,results$TITANstate,results$ClonalCluster)
  changepoints <- c(1,which(allstate[-1] != allstate[-length(allstate)])+1)
  segments <- results[changepoints,c(1,2,2,6,7,8,9,10,11,12)]
  names(segments)[2:3] <- c("Position1","Position2")
  segments$Position2 <- results$Position[c(changepoints[-1]-1,length(allstate))]
  segments[[2]] <- as.numeric(segments[[2]])
  segments[[3]] <- as.numeric(segments[[3]])
  segments[[4]] <- as.numeric(segments[[4]])
  segments[[5]] <- as.numeric(segments[[5]])
  segments[[6]] <- as.numeric(segments[[6]])
  segments[[7]] <- as.numeric(segments[[7]])
  segments[[10]] <- as.numeric(segments[[10]])
  segments$NumMarker <- diff(c(changepoints,length(allstate)+1))
  for (k in 1:nrow(segments)) {
    af <- as.numeric(results$AllelicRatio[changepoints[k]:(changepoints[k]+segments$NumMarker[k]-1)])
    segments$AllelicRatio[k] <- mean(pmax(af,1-af))
    lr <- as.numeric(results$LogRatio[changepoints[k]:(changepoints[k]+segments$NumMarker[k]-1)])
    segments$LogRatio[k] <- mean(lr)
  }
  
  ylim1 <- quantile(rep(segments$LogRatio,segments$NumMarker),c(0.0001,0.9999))
  par(pty="s")
  smoothScatter(rep(segments$AllelicRatio,segments$NumMarker),
		rep(segments$LogRatio,segments$NumMarker),
		xlim=c(0.5,1),ylim=ylim1,
		main = paste(sn,j), xlab="Allelic ratio",ylab="Log ratio")
  abline(h=log2(2/ploidy2),lty=2)
}
dev.off()
