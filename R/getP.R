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
                return(c(p1, af, bg, c(Dg)))
            }else{
                return(c(p1, af, bg))
            }
        }else{
            return(c(p1, af))
        }
    }, mc.cores=ncore))
    
    adata$pval = pval[,1]
    adata$allele_frequency = pval[,2]
    names(adata$pval) = rownames(Gall)
    names(adata$allele_frequency) = rownames(Gall)
    if(Beta|Sinv){
        beta = pval[,3:(M+3),drop=F]
        rownames(beta) = rownames(Gall)
        adata$Beta = beta
        if(Sinv){
            sinv = pval[,(M+4):(M+3+(M+1)^2),drop=F]
            rownames(sinv) = rownames(Gall)
            adata$Sinv = sinv/sigma2
        }
    }
    adata
}






#' Plot time-varying SNP effect sizes
#'
#' `plotEffectSize()` visualises the estimated dynamic effect size of a variant
#' across the support of the time or age variable. The variant-specific effect is
#' reconstructed from the Gaussian process basis and plotted with an approximate
#' 95% confidence interval.
#'
#' @param adata A `DynamicPGS` object after association mapping. The object must
#'   contain `Beta`, `Sinv`, `ta`, `rho`, and related Gaussian process
#'   components.
#' @param variant Integer or character. Variant to plot. If an integer is supplied,
#'   it is used as the row index of `adata$Beta`. If a character string is supplied,
#'   it is matched against `rownames(adata$Beta)`.
#' @param xstar Numeric vector. Points at which the effect size is evaluated. If
#'   `NULL`, a sequence over `adata$support_x` is used.
#' @param col Colour used for the effect-size curve and confidence band.
#' @param add Logical. If `TRUE`, add the curve to an existing plot.
#' @param ci Logical. If `TRUE`, draw an approximate 95% confidence interval.
#' @param lwd Line width for the effect-size curve.
#' @param main Main title of the plot. If `NULL`, the variant name or index is used.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#'
#' @return Invisibly returns a `data.frame` containing the evaluation points,
#'   estimated effect sizes, and confidence interval limits.
#'
#' @details
#' The dynamic effect size is computed as
#' \deqn{G_0(x) \hat{\beta},}
#' where `G0` is the Gaussian process basis evaluated at `xstar`, and
#' `\hat{\beta}` is the variant-specific coefficient vector stored in
#' `adata$Beta`. The uncertainty is approximated using the inverse precision
#' matrix stored in `adata$Sinv`.
#'
#' @examples
#' \dontrun{
#' plotEffectSize(adata, variant=1)
#' plotEffectSize(adata, variant="rs12345", col=2)
#' plotEffectSize(adata, variant=1, xstar=seq(0, 60, by=1))
#' }
#'
#' @export
plotEffectSize <- function(adata, variant, xstar=NULL, col=1, add=FALSE, ci=TRUE, lwd=2, main=NULL, xlab="x", ylab="Effect size", genoEffect=F, ...){
    
    col.geno = c("#19324B", "#32C8E9", "#F51A57")
    if(genoEffect){ci=F}
    
    if(is.null(xstar)) xstar <- Seq(adata$support_x)
    ta <- adata$ta
    rho <- adata$rho*adata$r_rho
    M <- adata$M
    Sinv <- adata$Sinv
    B <- adata$Beta
    
    if(is.null(B)){stop("No effect size estimate. Use getP(...,Beta=T,Sinv=T).")}

    if(is.character(variant)) variant <- match(variant, rownames(B))
    if(length(variant) != 1 || is.na(variant)) stop("variant must be a valid row index or row name of adata$Beta.")
    if(is.null(main)) main <- if(!is.null(rownames(B))) rownames(B)[variant] else paste0("variant ", variant)

    popavg = 0
    if(genoEffect){
        beta0 = c(tail(adata$PhiXty,M), adata$PhiXty[1])
        Knm <- getK(xstar, ta, adata$rho)
        Kmm <- getK(ta, ta, adata$rho)
        R <- chol(Kmm)
        tKnm <- t(forwardsolve(t(R), t(Knm)))
        G00 <- cbind(tKnm, 1)
        popavg = c(G00%*%beta0)
    }
    
    Knm <- getK(xstar, ta, rho)
    Kmm <- getK(ta, ta, rho)
    R <- chol(Kmm)
    tKnm <- t(forwardsolve(t(R), t(Knm)))
    G0 <- cbind(tKnm, 1)

    af = adata$allele_frequency[variant]
    b <- B[variant,]
    S <- matrix(Sinv[variant,], M+1, M+1)
    S <- (S + t(S)) / 2

    m <- as.numeric(G0 %*% b)
    V <- try(solve(S, t(G0)), silent=TRUE)
    if(inherits(V, "try-error")) V <- SolveJitter(S, t(G0))
    v <- colSums(V * t(G0))
    v <- pmax(v, 0)

    upper <- m + 1.96 * sqrt(v)
    lower <- m - 1.96 * sqrt(v)

    if(!add){
        ylim <- if(ci) range(c(lower, upper), na.rm=TRUE) else range(m+popavg, na.rm=TRUE)
        plot(xstar, m, type="n", ylim=ylim, main=main, xlab=xlab, ylab=ylab, axes=FALSE, ...)
        axis(2, las=2)
        axis(1)
        box()
        abline(h=0, lty=2)
    }

    if(genoEffect){
        lines(xstar, popavg + (0-af*2)*m, col=col.geno[1], lwd=lwd)
        lines(xstar, popavg + (1-af*2)*m, col=col.geno[2], lwd=lwd)
        lines(xstar, popavg + (2-af*2)*m, col=col.geno[3], lwd=lwd)
    }else{
        if(ci) polygon(c(xstar, rev(xstar)), c(upper, rev(lower)), col=Alpha(col), border=NA)
        lines(xstar, m, col=col, lwd=lwd)
    }
    invisible(data.frame(x=xstar, beta=m, lower=lower, upper=upper))
}
