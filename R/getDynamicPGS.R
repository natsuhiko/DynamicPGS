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
#' @param xstar Numeric vector of target index values at which dynamic PGS should
#'   be evaluated.
#' @param Gall Numeric genotype dosage matrix with variants in rows and samples
#'   in columns.
#' @param af Optional named numeric vector of allele frequencies for variants in
#'   `Gall`. If `NULL`, allele frequencies are estimated as row means divided by
#'   two.
#' @param adata A `DynamicPGS` object after [getP()] with `Beta`. If standard
#'   errors are required, `Sinv` should also be present.
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
#' dpgs <- getDynamicPGS(xstar = 0:60, Gall = Gall, adata = adata)
#' matplot(dpgs$xstar, dpgs$E, type = "l")
#' }
#'
#' @export
getDynamicPGS = function(xstar=0:54, Gall, af=NULL, adata){
    
    if(is.null(af)){ af = rowMeans(Gall,na.rm=T)/2 }
    
    ta = adata$ta
    rho = adata$rho
    M = adata$M
    Sinv = adata$Sinv
    B = adata$Beta
    L0 = nrow(B)
    sigma2 = adata$sigma2
    beta0 = c(tail(adata$PhiXty,M), adata$PhiXty[1])
    
    com = intersect(rownames(Gall),rownames(B))
    af = af[match(com,rownames(Gall))]
    B = B[com,]
    Gall = Gall[com,]
    L = nrow(B)
    cat(paste(nrow(Gall)," variants out of ",L0," index variants were found...\n",sep=""))
    
    Knm = getK(xstar, ta, rho)
    Kmm = getK(ta, ta, rho)
    R = chol(Kmm)
    tKnm = t(forwardsolve(t(R),t(Knm)))
    G0 = cbind(tKnm, 1)
    
    cat("Computing dynamic PRS at x*\n")
    E=V=0
    for(l in 1:L){
        gl <- Gall[l,]
        gl = gl - af[l]*2
        E=E+G0%*%B[l,]%*%t(gl)
        V=V+sigma2*colSums(solve(matrix(Sinv[l,],M+1),t(G0))*t(G0))%*%t(gl^2)
    }
    list(xstar=xstar, avg=G0%*%beta0, E=E, SE=sqrt(V), sigma2=sigma2)
}
