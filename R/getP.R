gchisq_p <- function(s, lambda, tol=0){
    lambda <- lambda[is.finite(lambda)]
    if(length(lambda)==0 || !is.finite(s) || s < 0) return(NA_real_)
    mx <- max(lambda)
    lambda <- lambda[lambda > mx * tol]
    if(length(lambda)==0) return(NA_real_)
    o <- CompQuadForm::davies(q=s, lambda=lambda, acc=1e-16, lim=50000)
    p <- o$Qq
    if(!is.finite(p) || o$ifault != 0){
        p <- CompQuadForm::imhof(q=s, lambda=lambda, epsabs=1e-16, epsrel=1e-16, limit=50000)$Qq
    }
    pmin(pmax(p, .Machine$double.xmin), 1)
}

#' Dynamic association testing for genotype dosage matrix
#'
#' `getP()` performs dynamic association testing for each variant in a genotype
#' dosage matrix. For each variant, genotype dosages are centred by allele
#' frequency and tested against the dynamic GP basis using a mixture chi-square
#' distribution.
#'
#' This function requires a `DynamicPGS` object after [gpreg2()] or
#' [prep_assoc()]. Optionally, posterior effect estimates and inverse covariance
#' matrices can be stored for downstream dynamic PGS calculation.
#'
#' @param adata A `DynamicPGS` object prepared by [gpreg2()] or [prep_assoc()].
#' @param Gall Numeric genotype dosage matrix with variants in rows and
#'   individuals in columns. Row names should be variant IDs. Columns should
#'   correspond to the unique individuals in `adata`.
#' @param delta2g Numeric prior variance parameter for variant effects.
#' @param Beta Logical. If `TRUE`, stores posterior effect estimates for each
#'   variant in `adata$Beta`.
#' @param Sinv Logical. If `TRUE`, stores inverse covariance matrices for each
#'   variant in `adata$Sinv`. This is required for standard errors in
#'   [getDynamicPGS()].
#' @param ncore Integer. Number of CPU cores used for parallel computation.
#'
#' @return The input `DynamicPGS` object updated with `pval`. If requested, also
#'   contains `Beta` and `Sinv`.
#'
#' @examples
#' \dontrun{
#' adata <- getP(adata, Gall, Beta = TRUE, Sinv = TRUE, ncore = 8)
#' head(adata$pval)
#' }
#'
#' @export
getP <- function(adata, Gall, delta2g=0.01, Beta=F, Sinv=F, ncore=max(1, parallel::detectCores()-1)){
    
    if(is.null(adata$PhiKdty)){
        stop("Step 2 optimization has not been performed yet...aborted")
    }
    
    maxni = adata$max_family_size
    N = adata$N
    M = adata$M
    P = adata$P
    Q = adata$Q
    
    y=adata$y
    
    G0=adata$G0
    tX=adata$tX
    D=adata$D
    KdtG0base=adata$KdtG0base
    mapid=adata$mapid
    mapid2=adata$mapid2
    CinvKdtG0base=adata$CinvKdtG0base
    CinvB=adata$CinvB
    DinvBtCinv=adata$DinvBtCinv
    PhiXty=adata$PhiXty
    PhiKdty=adata$PhiKdty
    sigma2=adata$sigma2
    
    
    pval = do.call("rbind", parallel::mclapply(1:nrow(Gall), function(l){
        print(l)
        gl = Gall[l,] # rbinom(Nd,2,0.3)
        af = mean(gl)/2
        gl = (gl-2*af) #/sqrt(af*2*(1-af))
        
        G = G0*gl[mapid]
        
        XtG = t(tX)%*%G
        DinvXtG = Solve(D, XtG)
        KdtG     =     KdtG0base[,1:(M+1)]*(c(0,gl)[mapid2[,1]+1])
        CinvKdtG = CinvKdtG0base[,1:(M+1)]*(c(0,gl)[mapid2[,1]+1])
        if(maxni>1){for(j in 2:maxni){
            KdtG =         KdtG +     KdtG0base[,1:(M+1)+(M+1)*(j-1)]*(c(0,gl)[mapid2[,j]+1]) # Nd(M+1) x (M+1)
            CinvKdtG = CinvKdtG + CinvKdtG0base[,1:(M+1)+(M+1)*(j-1)]*(c(0,gl)[mapid2[,j]+1]) # Nd(M+1) x (M+1)
        }}
        DinvBtCinvKdtG = DinvBtCinv%*%KdtG # (Q+M) x (M+1)
        PhiZtG = rbind(DinvXtG - DinvBtCinvKdtG, -t(DinvBtCinv)%*%XtG + CinvKdtG + CinvB%*%DinvBtCinvKdtG)
        
        ytVinvG = t(y)%*%G - t(c(PhiXty,PhiKdty))%*%rbind(XtG,KdtG)
        GVinvG  = t(G)%*%G - t(rbind(XtG,KdtG))%*%PhiZtG
        GVinvG  = (GVinvG + t(GVinvG)) / 2
        s = sum(ytVinvG^2)/sigma2
        lambda = eigen(GVinvG, symmetric=TRUE, only.values=TRUE)$values
        p1 = try(gchisq_p(s/N, lambda/N))
        if(is.character(p1)){
            pval = 1
        }
        
        if(Beta|Sinv){
            Dg = diag(M+1)/delta2g + GVinvG
            bg = Solve(Dg, t(ytVinvG))
            if(Sinv){
                return(c(p1, bg, c(Dg)))
            }else{
                return(c(p1, bg))
            }
        }else{
            return(p1)
        }
    }, mc.cores=ncore))
    
    adata$pval = pval[,1]
    names(adata$pval) = rownames(Gall)
    if(Beta|Sinv){
        beta = pval[,2:(M+2),drop=F]
        rownames(beta) = rownames(Gall)
        adata$Beta = beta
        if(Sinv){
            sinv = pval[,(M+3):(M+2+(M+1)^2),drop=F]
            rownames(sinv) = rownames(Gall)
            adata$Sinv = sinv/sigma2
        }
    }
    adata
}
