#' Compute dynamic polygenic scores over a continuous index
#'
#' `getDynamicPGS()` computes dynamic polygenic scores at user-specified values
#' of the continuous index, such as age or time. It combines genotype dosages
#' with variant-level dynamic effect estimates obtained by [getP()].
#'
#' The genotype matrix should contain the variants used in `adata$Beta`. Missing
#' variants are ignored, and the function reports how many index variants were
#' found.
#'
#' @param adata A `DynamicPGS` object after [getP()] with `Beta`. If standard
#'   errors are required, `Sinv` should also be present.
#' @param Gall Numeric genotype dosage matrix with variants in rows and samples
#'   in columns.
#' @param xstar Numeric vector of target index values at which dynamic PGS should
#'   be evaluated.
#' @param af Optional named numeric vector of allele frequencies for variants in
#'   `Gall`. If `NULL`, allele frequencies are estimated as row means divided by
#'   two.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{xstar}{Target index values.}
#'   \item{avg}{Estimated average trajectory at `xstar`.}
#'   \item{E}{Dynamic PGS matrix evaluated at `xstar`. Rows correspond to
#'     `xstar`; columns correspond to samples.}
#'   \item{SE}{Approximate standard errors of dynamic PGS values.}
#'   \item{sigma2}{Residual variance estimate from the fitted model.}
#' }
#'
#' @examples
#' \dontrun{
#' adata <- getDynamicPGS(adata, Gall = Gall, xstar = 0:60)
#' # PGS curve with standard error for the first subject
#' plot(adata, 1, Prediction=F)
#' # prediction interval for the first subject
#' plot(adata, 1, Prediction=T)
#' }
#'
#' @export
getDynamicPGS = function(adata, Gall, xstar=NULL, af=adata$af){
    
    if(is.null(af)){ af = rowMeans(Gall,na.rm=T)/2; names(af)=rownames(Gall) }
    if(is.null(xstar)){ xstar = Seq(adata$support_x)}
    
    ta = adata$ta
    rho = adata$rho*adata$r_rho
    M = adata$M
    Sinv = adata$Sinv
    B = adata$Beta
    L0 = nrow(B)
    beta0 = c(tail(adata$PhiXty,M), adata$PhiXty[1])
    
    com = intersect(rownames(Gall),rownames(B))
    af = af[com]
    B = B[com,,drop=F]
    Gall = Gall[com,,drop=F]
    L = nrow(B)
    cat(paste(nrow(Gall)," variants out of ",L0," proxy variants were found...\n",sep=""))
    
    Knm = getK(xstar, ta, adata$rho)
    Kmm = getK(ta, ta, adata$rho)
    R = chol(Kmm)
    tKnm = t(forwardsolve(t(R),t(Knm)))
    G00 = cbind(tKnm, 1)
    
    Knm = getK(xstar, ta, rho)
    Kmm = getK(ta, ta, rho)
    R = chol(Kmm)
    tKnm = t(forwardsolve(t(R),t(Knm)))
    G0 = cbind(tKnm, 1)
    
    cat("Computing dynamic PGS at x*\n")
    E=V=0
    for(l in 1:L){
        gl <- Gall[l,]
        gl = gl - af[l]*2
        E=E+G0%*%B[l,]%*%t(gl)
        V=V+colSums(solve(matrix(Sinv[l,],M+1),t(G0))*t(G0))%*%t(gl^2)
    }
    adata$xstar=xstar
    adata$pop_avg_xstar=G00%*%beta0
    adata$PGS=E
    adata$PGS_SE=sqrt(V)
    
    adata
}



#' Plot DynamicPGS object
#'
#' @param x A DynamicPGS object.
#' @param i Integer. Index of the individual/sample to plot.
#' @param Prediction Logical. If TRUE, plot predicted phenotype including the population average.
#' @param col Line colour.
#' @param add Logical. If TRUE, add to the existing plot.
#' @param xlab,ylab Axis labels.
#' @param lwd Line width.
#' @param ci Logical. If TRUE, draw 95% confidence interval.
#' @param ... Additional arguments passed to plot().
#'
#' @export
#' @method plot DynamicPGS
plot.DynamicPGS <- function(x, i=NULL, Prediction=FALSE, col=1, add=FALSE, xlab="x", ylab=NULL, lwd=2, ci=TRUE, ...) {
    if(is.null(x$xstar)) stop("x$xstar is missing.")
    if(is.null(x$PGS)) stop("x$PGS is missing.")
    if(is.null(x$PGS_SE)) stop("x$PGS_SE is missing.")
    if(is.null(i)) i <- 1

    xstar <- x$xstar

    if(Prediction){
        if(is.null(x$pop_avg)) stop("x$pop_avg is missing.")
        if(is.null(x$sigma2)) stop("x$sigma2 is missing.")
        pred <- as.numeric(x$pop_avg + x$PGS[,i])
        se <- sqrt(x$sigma2 + x$PGS_SE[,i]^2)
        if(is.null(ylab)) ylab <- "Predicted phenotype"
    }else{
        pred <- as.numeric(x$PGS[,i])
        se <- as.numeric(x$PGS_SE[,i])
        if(is.null(ylab)) ylab <- "Dynamic PGS"
    }

    upper <- pred + 1.96 * se
    lower <- pred - 1.96 * se

    if(!add){
        ylim <- if(ci) range(c(lower, upper), na.rm=TRUE) else range(pred, na.rm=TRUE)
        plot(xstar, pred, type="n", xlab=xlab, ylab=ylab, ylim=ylim, ...)
    }
    if(ci){
        polygon(c(xstar, rev(xstar)), c(upper, rev(lower)), col=Alpha(col), border=NA)
    }
    lines(xstar, pred, col=col, lwd=lwd)
    invisible(x)
}


#' Get a public DynamicPGS object
#'
#' Remove individual-level data from a DynamicPGS object and retain only model-level
#' parameters required for public use.
#'
#' @param adata A DynamicPGS object.
#'
#' @return A DynamicPGS_public object containing only `ta`, `rho`, `M`, `Sinv`,
#'   `Beta`, `sigma2`, and `PhiXty`.
#' @export
getPublicData = function(adata){
    keep0 = c("ta","rho","M","Sinv","Beta","sigma2","PhiXty","proxy","support_x","r_rho","allele_frequency")
    keep = intersect(keep0, names(adata))
    out = adata[keep]
    missing = setdiff(keep0, names(adata))
    if(length(missing)>0) warning("Missing fields: ", paste(missing, collapse=", "))
    class(out) = c("PublicDynamicPGS", "DynamicPGS")
    out
}

print.PublicDynamicPGS = function(x, ...){
    cat("Public DynamicPGS object\n")
    if(!is.null(x$M)) cat("Number of inducing points: ", x$M, "\n", sep="")
    if(!is.null(x$rho)) cat("rho: ", x$rho, "\n", sep="")
    if(!is.null(x$ta)) cat("Support: ", min(x$support_x), " to ", max(x$support_x), "\n", sep="")
    if(!is.null(x$sigma2)) cat("sigma2: ", x$sigma2, "\n", sep="")
    if(!is.null(x$proxy)) cat("N lead vars: ", length(unique(x$proxy$lead)), "\n", sep="")
    if(!is.null(x$proxy)) cat("N proxy vars: ", length(unique(x$proxy$proxy)), "\n", sep="")
    invisible(x)
}
