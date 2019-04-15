##Use midpoints between markers as boundaries of segments
titancna2seg <- function(titanresult,titanparams,adjust_end=TRUE) {

  major_cn_code <- c(0,1,2,1,3,2,4,3,2,5,4,3,6,5,4,3,7,6,5,4,8,7,6,5,4)
  minor_cn_code <- c(0,0,0,1,0,1,0,1,2,0,1,2,0,1,2,3,0,1,2,3,0,1,2,3,4)

  chr_sizse <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)

  ploidy <- round(tail(titanparams$phi,1),3)
  n <- round(tail(titanparams$n,1),3)

  titanresult$Position <- as.integer(titanresult$Position)
  titanresult$LogRatio <- as.numeric(titanresult$LogRatio)
  titanresult$AllelicRatio <- as.numeric(titanresult$AllelicRatio)
  titanresult$CopyNumber <- as.numeric(titanresult$CopyNumber)
  titanresult$CellularPrevalence[titanresult$CellularPrevalence == "NA"] <- NA
  titanresult$CellularPrevalence <- as.numeric(titanresult$CellularPrevalence)
  titanresult$ClonalCluster[is.na(titanresult$ClonalCluster)] <- 0
  cp2 <- c(which(titanresult$TITANstate[-1] != titanresult$TITANstate[-nrow(titanresult)] | 
		   titanresult$Chr[-1] != titanresult$Chr[-nrow(titanresult)] | 
		     titanresult$ClonalCluster[-1] !=  titanresult$ClonalCluster[-nrow(titanresult)]),
	   nrow(titanresult))
  cp1 <- c(1,cp2[-length(cp2)]+1)

  cnv <- data.frame(chrom=titanresult$Chr[cp1],
                    loc.start=titanresult$Position[cp1],
                    loc.end=titanresult$Position[cp2],
		    num.mark=cp2-cp1+1,
                    seg.mean=titanresult$LogRatio[cp1],
		    copynumber=titanresult$CopyNumber[cp1],
		    minor_cn=minor_cn_code[as.integer(titanresult$TITANstate[cp1])+1],
		    major_cn=major_cn_code[as.integer(titanresult$TITANstate[cp1])+1],
		    allelicratio=titanresult$AllelicRatio[cp1],
		    LOHcall=titanresult$TITANcall[cp1],
		    cellularprevalence=titanresult$CellularPrevalence[cp1],
		    ploidy=ploidy,
		    normalproportion=n)
  for (j in 1:length(cp1)) {
    cnv$seg.mean[j] <- mean(titanresult$LogRatio[cp1[j]:cp2[j]])
    cnv$allelicratio[j] <- mean(0.5+abs(0.5-titanresult$AllelicRatio[cp1[j]:cp2[j]]))
    if (j < length(cp1)) {
      if (titanresult$Chr[cp2[j]] == titanresult$Chr[cp1[j+1]]) {
	cnv$loc.end[j] <- round((cnv$loc.end[j]+cnv$loc.start[j+1])/2)
      }
    }
    if (j > 1) {
      if (titanresult$Chr[cp1[j]] == titanresult$Chr[cp2[j-1]]) {
	cnv$loc.start[j] <- cnv$loc.end[j-1]+1
      }
    }
  }
  ##Change end points of segments at chromosome ends
  if (adjust_end) 
    for (i in unique(cnv$chrom)) {
      cnv$loc.start[head(which(cnv$chrom == i),1)] <- 1
      cnv$loc.end[tail(which(cnv$chrom == i),1)] <- 
	chr_sizes[as.integer(as.character(i))]
    }
  cnv$logcopynumberratio <- log2(((cnv$copynumber - 2) * cnv$cellularprevalence + 2)/2)
  cnv$logcopynumberratio[is.na(cnv$logcopynumberratio)] <- 0
  ##cnv$seg.mean.2 <- cnv$seg.mean
  for (i in c(5,7,12)) {
    cnv[[i]] <- round(cnv[[i]],3)
  }
  return(cnv)
}

