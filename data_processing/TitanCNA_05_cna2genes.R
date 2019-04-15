##obtain per genes values from cna segmentation results

library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

cna2genes <- function(cna,genelist) {

  ##cna: dataframe with chr,start,end as the first 3 columns followed by other values
  ##genelist: vector of gene names, "" for all genes

  cna_gr <- GRanges(seqnames=paste("chr",cna[[1]],sep=""),
                    ranges=IRanges(start=cna[[2]],end=cna[[3]]),
                    strand=rep("*",nrow(cna)),
		    cna[,-(1:3)])
  
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  allgenes <- genes(txdb)
  
  if (!identical(genelist,"")) {
    geneids <- select(org.Hs.eg.db,keys=genelist,keytype="SYMBOL",columns="ENTREZID")$ENTREZID
    geneidx <- match(geneids, names(allgenes))
    genelist <- genelist[!is.na(geneidx)]
    geneidx <- geneidx[!is.na(geneidx)]
    genes_gr <- allgenes[geneidx]
  }
  else {
    allgenes <- sort(allgenes)
    geneids <- names(allgenes)
    genelist <- select(org.Hs.eg.db,keys=geneids,columns="SYMBOL",keytype="ENTREZID")$SYMBOL
    genes_gr <- allgenes[!is.na(genelist)]
    genelist <- genelist[!is.na(genelist)]
  }
  strand(genes_gr) <- "*"
  idx <- which(seqnames(genes_gr) %in% paste("chr",c(1:22,"X","Y"),sep=""))
  genes_gr <- genes_gr[idx]
  genelist <- genelist[idx]

  cnahits <- findOverlaps(genes_gr,cna_gr)
  geneidx <- queryHits(cnahits)
  genecounthits <- countQueryHits(cnahits)
  cnaidx <- subjectHits(cnahits)

  genes1hit <- which(genecounthits==1)
  genes1hit_geneidx <- geneidx[geneidx %in% genes1hit]
  genes1hit_cnaidx <- cnaidx[geneidx %in% genes1hit]
  genes_cna <- cbind(gene=genelist[genes1hit_geneidx],cna[genes1hit_cnaidx,])
  genes1hit_gr <- pintersect(genes_gr[genes1hit_geneidx],cna_gr[genes1hit_cnaidx])
  genes_cna[,3] <- start(genes1hit_gr)
  genes_cna[,4] <- end(genes1hit_gr)

  genes2hit <- which(genecounthits>1)
  if (length(genes2hit) > 0) {
    for (i in genes2hit) {
      geneidx1 <- which(geneidx == i)
      for (j in geneidx1) {
	overlap_gr <- intersect(genes_gr[geneidx[j]],cna_gr[cnaidx[j]])
	overlap_cna <- cbind(gene=genelist[geneidx[j]],cna[cnaidx[j],])
	overlap_cna[,3] <- start(overlap_gr)
	overlap_cna[,4] <- end(overlap_gr)
	genes_cna <- rbind(genes_cna,overlap_cna,make.row.names=FALSE)
      }
    }
  }

  genes0hit <- which(genecounthits==0)
  if (length(genes0hit) > 0) {
    genes_cna0 <- cbind(gene=genelist[genes0hit],cna[1,])
    genes_cna0[,2] <- as.character(seqnames(genes_gr)[genes0hit])
    genes_cna0[,3] <- start(genes_gr)[genes0hit]
    genes_cna0[,4] <- end(genes_gr)[genes0hit]
    if (ncol(genes_cna0) >= 5)
      for (i in 5:ncol(genes_cna0)) {
	genes_cna0[,i] <- NA
      }
    genes_cna <- rbind(genes_cna,genes_cna0,make.row.names=FALSE)
  }

  genes_cna <- genes_cna[ order(match(genes_cna[,1],genelist)),]
  genes_cna[,2] <- sub("chr","",genes_cna[,2])
  rownames(genes_cna) <- NULL

  return(genes_cna)
}
