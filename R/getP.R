gchisq_p <- function(s, lambda, tol=0, switch_p=1e-12){
    lambda <- lambda[is.finite(lambda) & lambda > 0]
    if(length(lambda)==0 || !is.finite(s) || s < 0) return(NA_real_)
    mx <- max(lambda); lambda <- lambda[lambda > mx*tol]
    if(length(lambda)==0) return(NA_real_)
    spa <- function(q, l){
        mu <- sum(l)
        if(q <= mu) return(NA_real_)
        K <- function(t) -0.5*sum(log1p(-2*t*l))
        K1 <- function(t) sum(l/(1-2*t*l))
        K2 <- function(t) 2*sum(l^2/(1-2*t*l)^2)
        tmax <- (1/(2*max(l)))*(1-1e-12)
        tt <- tryCatch(uniroot(function(t) K1(t)-q, c(0,tmax), tol=.Machine$double.eps^0.5)$root, error=function(e) NA_real_)
        if(!is.finite(tt)) return(NA_real_)
        r <- sqrt(2*(tt*q-K(tt))); v <- tt*sqrt(K2(tt))
        if(!is.finite(r) || !is.finite(v) || r <= 0 || v <= 0) return(NA_real_)
        a <- pnorm(-r, log.p=TRUE); b <- dnorm(r, log=TRUE)+log(abs(1/v-1/r)); sg <- sign(1/v-1/r)
        logp <- if(sg > 0){m <- max(a,b); m+log(exp(a-m)+exp(b-m))} else {m <- max(a,b); z <- exp(a-m)-exp(b-m); if(z <= 0) return(exp(a)); m+log(z)}
        exp(logp)
    }
    o <- CompQuadForm::davies(q=s, lambda=lambda, acc=1e-12, lim=50000)
    p <- o$Qq
    if(!is.finite(p) || p < 0 || p > 1 || o$ifault != 0){
        p <- tryCatch(CompQuadForm::imhof(q=s, lambda=lambda, epsabs=1e-12, epsrel=1e-12, limit=50000)$Qq, error=function(e) NA_real_)
    }
    if(!is.finite(p) || p < switch_p){
        p2 <- spa(s, lambda)
        if(is.finite(p2)) p <- p2
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
