##modified function to use results from TitanCNA
##Weighted average for SNVs not passing probalility cutoff value

getSampleCCF2 <- function (vcf, dd, thr.coverage = 1, upper.cov = 1.0) 
{
    ##Probability cutoff
    prob_cutoff <- 0.95

    s <- vcf$Tumor_cover
    f <- vcf$Tumor_freq
    b <- round(s*f)
    a <- s - b

    s0 <- vcf$Normal_cover
    f0 <- vcf$Normal_freq
    b0 <- round(s0*f0)
    a0 <- s0 - b0

    thr.upper <- quantile(s, upper.cov)
    vv.cov <- which(s >= thr.coverage) ##( & s0 >= thr.coverage)

    sam.dd <- dd
    MutInfo <- c()
    n <- nrow(sam.dd)
    vv.int <- vv.cov

    cc <- seq(0.01, 1, by = 0.01)
    computeSD <- function(N, S, f, cc) {
        M1list <- c()
        M2list <- c()
        MLElist <- c()
        for (ii in 1:length(N)) {
            PF <- sum(dbinom(S[ii], N[ii], f), na.rm = TRUE)
            M1 <- sum(dbinom(S[ii], N[ii], f) * cc, na.rm = TRUE)/PF
            M2 <- sum(dbinom(S[ii], N[ii], f) * cc^2, na.rm = TRUE)/PF
            M1list <- c(M1list, M1)
            M2list <- c(M2list, M2)
            MLElist <- c(MLElist, cc[which.max(dbinom(S[ii], 
                N[ii], f))])
        }
        return(list(M1 = MLElist, SD = sqrt(M2list - M1list^2)))
    }
    tmp.vcf <- vcf[vv.int, ]
    for (k in 1:n) {
        vv <- which(tmp.vcf[, 1] == sam.dd[k, 1] & tmp.vcf[, 
            2] >= sam.dd[k, 2] & tmp.vcf[, 2] <= sam.dd[k, 3])
        if (length(vv) == 0) 
            next
        nb <- sam.dd[k, 9]
        nt <- sam.dd[k, 10]
        tmp.b <- b[vv.int[vv]]
        tmp.s <- s[vv.int[vv]]
        tmp.a <- tmp.s - tmp.b
        tmp.b0 <- b0[vv.int[vv]]
        tmp.s0 <- s0[vv.int[vv]]
        nv <- length(vv)
        if (is.na(sam.dd[k, 8])) {
            INFO <- data.frame(CCF=rep(NA,nv), 
			       CCF_std=rep(NA,nv),
			       sAGP=rep(NA, nv), 
			       n_minor=rep(NA, nv), 
			       n_total=rep(NA, nv),
			       Lineage=rep(NA, nv), 
			       stringsAsFactors=FALSE)
            tmp.dd <- tmp.vcf[vv, ]
            tmp.dd <- cbind(tmp.dd,INFO)
            MutInfo <- rbind(MutInfo, tmp.dd)
            next
        }
        p <- sam.dd[k, 8]
        if (p > 1) 
            p <- 1
        if (nb == 1 & nt == 2) {
            ff = cc/2
            Ms = computeSD(tmp.s, tmp.b, ff, cc)
            CCF <- Ms$M1
            SDs <- Ms$SD
            IsEarly <- rep("C", nv)
	}
        else if (nb == 0 & nt == 0) {
	    nc <- nt * p + 2 * (1 - p)
	    ff = cc[cc<=(1-p)]/nc
	    Ms = computeSD(tmp.s, tmp.b, ff, cc[cc<=(1-p)])
            CCF <- Ms$M1
            SDs <- Ms$SD
            IsEarly <- rep("C", nv)
        }	
        else if (nt == 1) {
            nc <- nt * p + 2 * (1 - p)
            ff <- cc/nc
            Ms <- computeSD(tmp.s, tmp.b, ff, cc)
            CCF <- Ms$M1
            SDs <- Ms$SD
            fh.ea <- (p * (nt - nb) + 1 - p)/nc
            fl.ea <- (p * (nt - nb))/nc
            fh.t <- p/nc
            fh.e <- (1 - p)/nc
            pEarly.a <- pbeta(rep(fh.ea, nv), tmp.b + 1, tmp.a + 
                1) - pbeta(rep(fl.ea, nv), tmp.b + 1, tmp.a + 
                1)
            pLate <- pbeta(rep(fh.t, nv), tmp.b + 1, tmp.a + 
                1)
            pEuploid <- pbeta(rep(fh.e, nv), tmp.b + 1, tmp.a + 
                1)
            Ptot <- pEarly.a + pLate + pEuploid
            cp.A <- pEarly.a/Ptot
            cp.CD <- 1 - cp.A
            cp.C <- pLate/Ptot
            cp.D <- pEuploid/Ptot
            cp.AC <- 1 - cp.D
            cp.AD <- 1 - cp.C
            vv.A <- which(cp.A >= prob_cutoff)
            vv.CD <- which(cp.CD >= prob_cutoff & cp.C < prob_cutoff & cp.D < 
                prob_cutoff)
            vv.C <- which(cp.C >= prob_cutoff)
            vv.D <- which(cp.D >= prob_cutoff)
            vv.AC <- which(cp.AC >= prob_cutoff & cp.A < prob_cutoff & cp.C < 
                prob_cutoff)
            vv.AD <- which(cp.AD >= prob_cutoff & cp.A < prob_cutoff & cp.D < 
                prob_cutoff)
            IsEarly <- rep("A1/B/C", nv)
            IsEarly[vv.A] <- "A1"
            IsEarly[vv.C] <- "B"
            IsEarly[vv.D] <- "C"
            IsEarly[vv.AC] <- "A1/B"
            IsEarly[vv.AD] <- "A1/C"
            IsEarly[vv.CD] <- "B/C"
        }
        else if (nb == 0 | nt == 2 * nb) {
            nc <- nt * p + 2 * (1 - p)
            fh.ea <- (p * (nt - nb) + 1 - p)/nc
            fl.ea <- (p * (nt - nb))/nc
            fh.t <- p/nc
            fh.e <- (1 - p)/nc
            pEarly.a <- pbeta(rep(fh.ea, nv), tmp.b + 1, tmp.a + 
                1) - pbeta(rep(fl.ea, nv), tmp.b + 1, tmp.a + 
                1)
            pLate <- pbeta(rep(fh.t, nv), tmp.b + 1, tmp.a + 
                1)
            pEuploid <- pbeta(rep(fh.e, nv), tmp.b + 1, tmp.a + 
                1)
            Ptot <- pEarly.a + pLate + pEuploid
            cpEarly.a <- pEarly.a/Ptot
            cpLate.eup <- 1 - cpEarly.a
            cpLate <- pLate/Ptot
            cpEup <- pEuploid/Ptot
            vv.early.a <- which(cpEarly.a >= prob_cutoff)
            vv.late.eup <- which(cpLate.eup >= prob_cutoff)
            vv.late <- which(cpLate >= prob_cutoff)
            vv.Eup <- which(cpEup >= prob_cutoff)
            z <- tmp.b/tmp.s
            CCF <- rep(NA, nv)
            SDs <- rep(NA, nv)
            ff.early = (cc[cc>=p] - p + (nt - nb) * p)/nc
            Ms.early = computeSD(tmp.s[vv.early.a], tmp.b[vv.early.a], 
                ff.early, cc[cc>=p])
            CCF[vv.early.a] <- Ms.early$M1
            SDs[vv.early.a] <- Ms.early$SD
            ff.late = cc/nc
            Ms.late = computeSD(tmp.s[c(vv.late, vv.late.eup)], 
                tmp.b[c(vv.late, vv.late.eup)], ff.late, cc)
            CCF[c(vv.late, vv.late.eup)] <- Ms.late$M1
            SDs[c(vv.late, vv.late.eup)] <- Ms.late$SD

	    ##Weighted average
	    idx <- which(is.na(CCF))
	    Ms.early <- computeSD(tmp.s[idx], tmp.b[idx], ff.early, cc[cc>=p])
	    Ms.late <- computeSD(tmp.s[idx], tmp.b[idx], ff.late, cc)
	    CCF[idx] <- Ms.early$M1 * cpEarly.a[idx] + 
	      Ms.late$M1 * cpLate.eup[idx]

            CCF <- round(CCF, 3)
            IsEarly <- rep("A1/B/C", nv)
            IsEarly[vv.early.a] <- "A1"
            IsEarly[vv.late.eup] <- "B/C"
            IsEarly[vv.late] <- "B"
            IsEarly[vv.Eup] <- "C"
        }
        else if (nb >= 1 & nt > 2) {
            nc <- nt * p + 2 * (1 - p)
            fh.ea <- (p * (nt - nb) + 1 - p)/nc
            fl.ea <- (p * (nt - nb))/nc
            fh.eb <- (nb * p + 1 - p)/nc
            fl.eb <- nb * p/nc
            fh.t <- p/nc
            fh.e <- (1 - p)/nc
            pEarly.a <- pbeta(rep(fh.ea, nv), tmp.b + 1, tmp.a + 
                1) - pbeta(rep(fl.ea, nv), tmp.b + 1, tmp.a + 
                1)
            pEarly.b <- pbeta(rep(fh.eb, nv), tmp.b + 1, tmp.a + 
                1) - pbeta(rep(fl.eb, nv), tmp.b + 1, tmp.a + 
                1)
            pLate <- pbeta(rep(fh.t, nv), tmp.b + 1, tmp.a + 
                1)
            pEuploid <- pbeta(rep(fh.e, nv), tmp.b + 1, tmp.a + 
                1)
            Ptot <- pEarly.a + pEarly.b + pLate + pEuploid
            cp.A <- pEarly.a/Ptot
            cp.B <- pEarly.b/Ptot
            cp.C <- pLate/Ptot
            cp.D <- pEuploid/Ptot
            cp.AB <- 1 - cp.C - cp.D
            cp.AC <- 1 - cp.B - cp.D
            cp.AD <- 1 - cp.B - cp.C
            cp.BC <- 1 - cp.A - cp.D
            cp.BD <- 1 - cp.A - cp.C
            cp.CD <- 1 - cp.A - cp.B
            cp.ABC <- 1 - cp.D
            cp.ABD <- 1 - cp.C
            cp.ACD <- 1 - cp.B
            cp.BCD <- 1 - cp.A
            vv.A <- which(cp.A >= prob_cutoff)
            vv.B <- which(cp.B >= prob_cutoff)
            vv.C <- which(cp.C >= prob_cutoff)
            vv.D <- which(cp.D >= prob_cutoff)
            vv.CD <- which(cp.CD >= prob_cutoff & cp.C < prob_cutoff & cp.D < 
                prob_cutoff)
            vv.AB <- which(cp.AB >= prob_cutoff & cp.A < prob_cutoff & cp.B < 
                prob_cutoff)
            vv.AC <- which(cp.AC >= prob_cutoff & cp.A < prob_cutoff & cp.C < 
                prob_cutoff)
            vv.AD <- which(cp.AD >= prob_cutoff & cp.A < prob_cutoff & cp.D < 
                prob_cutoff)
            vv.BC <- which(cp.BC >= prob_cutoff & cp.B < prob_cutoff & cp.C < 
                prob_cutoff)
            vv.BD <- which(cp.BD >= prob_cutoff & cp.B < prob_cutoff & cp.D < 
                prob_cutoff)
            vv.BCD <- which(cp.BCD >= prob_cutoff & cp.BC < prob_cutoff & cp.BD < 
                prob_cutoff & cp.CD < prob_cutoff & cp.B < prob_cutoff & cp.C < prob_cutoff & 
                cp.D < prob_cutoff)
            vv.ABC <- which(cp.ABC >= prob_cutoff & cp.BC < prob_cutoff & cp.AB < 
                prob_cutoff & cp.AC < prob_cutoff & cp.B < prob_cutoff & cp.C < prob_cutoff & 
                cp.A < prob_cutoff)
            vv.ABD <- which(cp.ABD >= prob_cutoff & cp.AB < prob_cutoff & cp.AD < 
                prob_cutoff & cp.BD < prob_cutoff & cp.B < prob_cutoff & cp.D < prob_cutoff & 
                cp.A < prob_cutoff)
            vv.ACD <- which(cp.ACD >= prob_cutoff & cp.AC < prob_cutoff & cp.AD < 
                prob_cutoff & cp.CD < prob_cutoff & cp.A < prob_cutoff & cp.D < prob_cutoff & 
                cp.C < prob_cutoff)
            z <- tmp.b/tmp.s
            CCF <- rep(NA, nv)
            SDs <- rep(NA, nv)
            ff.A <- (cc[cc>=p] - p + (nt - nb) * p)/nc
            Ms.A <- computeSD(tmp.s[vv.A], tmp.b[vv.A], ff.A, cc[cc>=p])
            CCF[vv.A] <- Ms.A$M1
            SDs[vv.A] <- Ms.A$SD
            ff.B <- (cc[cc>=p] - p + nb * p)/nc
            Ms.B <- computeSD(tmp.s[vv.B], tmp.b[vv.B], ff.B, cc[cc>=p])
            CCF[vv.B] <- Ms.B$M1
            SDs[vv.B] <- Ms.B$SD
            ff.C <- cc/nc
            Ms.C <- computeSD(tmp.s[c(vv.C, vv.D, vv.CD)], tmp.b[c(vv.C, 
                vv.D, vv.CD)], ff.C, cc)
            CCF[c(vv.C, vv.D, vv.CD)] <- Ms.C$M1
            SDs[c(vv.C, vv.D, vv.CD)] <- Ms.C$SD
            if (nb == 1) {
                ff.BCD <- cc/nc
                Ms.BCD <- computeSD(tmp.s[c(vv.BCD, vv.BC, vv.BD)], 
                  tmp.b[c(vv.BCD, vv.BC, vv.BD)], ff.BCD, cc)
                CCF[c(vv.BCD, vv.BC, vv.BD)] <- Ms.BCD$M1
                SDs[c(vv.BCD, vv.BC, vv.BD)] <- Ms.BCD$SD
            }
 	    
	    ##Weighted average
	    idx <- which(is.na(CCF))
            Ms.A <- computeSD(tmp.s[idx], tmp.b[idx], ff.A, cc[cc>=p])
	    Ms.B <- computeSD(tmp.s[idx], tmp.b[idx], ff.B, cc[cc>=p])
	    Ms.C <- computeSD(tmp.s[idx], tmp.b[idx], ff.C, cc)
	    CCF[idx] <- Ms.A$M1 * cp.A[idx] +  Ms.B$M1 * cp.B[idx] +
	      Ms.C$M1 * cp.CD[idx]

	    IsEarly <- rep("A1/A2/B/C", nv)
            IsEarly[vv.A] <- "A1"
            IsEarly[vv.BCD] <- "A2/B/C"
            IsEarly[vv.CD] <- "B/C"
            IsEarly[vv.B] <- "A2"
            IsEarly[vv.C] <- "B"
            IsEarly[vv.D] <- "C"
            IsEarly[vv.AB] <- "A1/A2"
            IsEarly[vv.AC] <- "A1/B"
            IsEarly[vv.AD] <- "A1/C"
            IsEarly[vv.BC] <- "A2/B"
            IsEarly[vv.BD] <- "A2/C"
            IsEarly[vv.ACD] <- "A1/B/C"
            IsEarly[vv.ABD] <- "A1/A2/C"
            IsEarly[vv.ABC] <- "A1/A2/B"
        }
        p <- round(p, 3)
        CCF <- round(CCF, 3)
        SDs <- round(SDs, 3)         
	INFO <- data.frame(CCF=CCF, 
			   CCF_std=SDs,
			   sAGP=rep(p,nv), 
			   n_minor=rep(nb, nv), 
			   n_total=rep(nt, nv),
			   Lineage=IsEarly, 
			   stringsAsFactors=FALSE)
        tmp.dd <- tmp.vcf[vv, ]
        tmp.dd <- cbind(tmp.dd,INFO)
        MutInfo <- rbind(MutInfo, tmp.dd)
    }
    if (is.null(MutInfo)) 
        return()
    return(MutInfo)
}
