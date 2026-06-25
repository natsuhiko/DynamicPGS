#' Prepare intermediate matrices for dynamic association mapping
#'
#' `prep_assoc()` computes and stores the linear algebra quantities required for
#' variant-level dynamic association testing. It is usually called automatically
#' by [gpreg2()], but can also be called explicitly when model parameters are
#' modified or recalibrated.
#'
#' The function constructs GP basis matrices for genotype effects, handles
#' family-block Cholesky factors, and stores quantities used repeatedly by
#' [getP()] to avoid recomputing expensive matrix operations for each variant.
#'
#' @param adata A `DynamicPGS` object after [gpreg1()] and [gpreg2()], containing
#'   `rho` and `delta2d`.
#' @param r_rho Numeric multiplier applied to the GP length-scale used for
#'   genotype effects.
#' @param r_delta2d Numeric multiplier applied to the individual-level variance
#'   components used in the genotype-effect model.
#' @param ncore Integer. Number of CPU cores used for parallel computation.
#'
#' @return The input `DynamicPGS` object updated with association-mapping
#'   intermediate quantities, including `G0`, `tX`, `D`, `KdtG0base`,
#'   `CinvKdtG0base`, `CinvB`, `DinvBtCinv`, `PhiXty`, `PhiKdty`, `sigma2`,
#'   `sigma0`, and individual mapping indices.
#'
#' @examples
#' \dontrun{
#' adata <- prep_assoc(adata, r_rho = 1, r_delta2d = 1, ncore = 4)
#' }
#'
#' @export
prep_assoc = function(adata, r_rho=1, r_delta2d=1, ncore=max(1, parallel::detectCores()-1)){
    
    if(is.null(adata$rho)){
        stop("Step 1 optimization has not been performed yet...aborted")
    }
    if(is.null(adata$delta2d)){
        stop("Step 2 optimization has not been performed yet...aborted")
    }
    
    y = adata$y
    x = adata$x
    iid = adata$iid
    uiid = unique(iid)
    fid = adata$fid
    X = adata$X
    delta2 = adata$delta2
    delta2d= adata$delta2d
    rho = adata$rho
    ta = adata$ta
    nh = adata$nh
    Lmat = adata$Lmat[match(iid, unique(iid)),]
    if(sum(Lmat$IID == iid) != length(iid)){
        return("Lmat is inappropriate!")
    }else{
        Lmat = as.matrix(Lmat[,-c(1:2),drop=F])
    }

    N = length(y)
    M = length(ta)
    Nd = length(table(iid))
    Nf = length(table(fid))

    if(Verbose){
        cat("Nd: "); cat(Nd); cat("\n")
        cat("Nf: "); cat(Nf); cat("\n")
        cat("N: "); cat(N); cat("\n")
        cat("M: "); cat(M); cat("\n")
    }
    P = length(nh)
    Q = ncol(X)
    if(Verbose){ cat("Q: "); cat(Q); cat("\n") }

    Knm = getK(x, ta, rho)
    Kmm = getK(ta, ta, rho)
    R = chol(Kmm)
    tKnm = t(forwardsolve(t(R), t(Knm)))
    tX = cbind(X, tKnm)
    XtX = crossprod(tX)
    Xty = crossprod(tX, y)

    mapid = match(iid, uiid)
    maxni = adata$max_family_size

    #
    # GP for SNP genotypes
    #
    Knm2 = getK(x, ta, rho*r_rho)
    Kmm2 = getK(ta, ta, rho*r_rho)
    R2 = chol(Kmm2)
    tKnm2 = t(forwardsolve(t(R2), t(Knm2)))
    G0 = cbind(tKnm2, 1)

    one_id_fun_map = function(ii){
        uii = unique(iid[ii])
        tx = tX[ii,,drop=FALSE]
        n1 = length(uii)
        y1 = y[ii]
        g0 = G0[ii,,drop=FALSE]
        tk = C1prior = NULL

        for(j in 1:n1){
            tk = cbind(tk, cbind(tKnm[ii,,drop=FALSE], 1) * Lmat[ii,j])
            C1prior = dbind(C1prior, diag(rep(1.0 / r_delta2d / delta2d, c(M,1))))
        }

        K = (M + 1) * n1
        C1 = C1prior + crossprod(tk)
        ldc1 = getLogDet(C1)
        B1 = crossprod(tk, tx)
        tky = crossprod(tk, y1)

        tkg = NULL
        for(j in uii){
            rows_j = iid[ii] == j
            tkg = cbind(tkg, crossprod(tk[rows_j,,drop=FALSE], g0[rows_j,,drop=FALSE]))
        }

        if(maxni > n1){
            tkg = cbind(tkg, matrix(0, K, (maxni-n1)*(M+1)))
        }
        
        res = Solve(C1, cbind(B1, tky, tkg))
        
        # calculating standard errors (Cinv + CinvBDinvBtCinv) for each indiv
        C1inv = Solve(C1,diag(nrow(C1)))
        C1invi = NULL
        for(j in 1:n1){
            C1invi = rbind(C1invi, C1inv[((j-1)*(M+1)+1):(j*(M+1)),((j-1)*(M+1)+1):(j*(M+1))])
        }

        cbind(
            B1,
            tky,
            tkg,
            res,
            matrix(1,nrow(B1),1) %*% matrix(c(match(uii,uiid), rep(0,maxni-n1)), nrow=1),
            C1invi
        )
    }
    ufid = unique(fid)
    fid_fac = factor(fid, levels=ufid)
    idx_list = split(seq_len(N), fid_fac)
    CinvB = do.call("rbind", parallel::mclapply(idx_list, one_id_fun_map, mc.cores=ncore))
    A = XtX + diag(1/delta2[rep(1:P,nh)])
    B = CinvB[,1:(Q+M)]
    Kdty = CinvB[,Q+M+1]
    KdtG0base = CinvB[,(Q+M+1+1):(Q+M+1+(M+1)*maxni)]

    CinvKdty = CinvB[,(Q+M+1+(M+1)*maxni+Q+M+1)]
    CinvKdtG0base = CinvB[,(Q+M+1+(M+1)*maxni+(Q+M)+1+1):(Q+M+1+(M+1)*maxni+(Q+M)+1+(M+1)*maxni)]
    mapid2 = CinvB[,(Q+M+1+(M+1)*maxni+(Q+M)+1+(M+1)*maxni+1):(Q+M+1+(M+1)*maxni+(Q+M)+1+(M+1)*maxni+maxni),drop=FALSE]
    Cinvi = CinvB[,(Q+M+1+(M+1)*maxni+(Q+M)+1+(M+1)*maxni+maxni+1):(Q+M+1+(M+1)*maxni+(Q+M)+1+(M+1)*maxni+maxni+(M+1))]
    CinvB = CinvB[,(Q+M+1+(M+1)*maxni+1):(Q+M+1+(M+1)*maxni+(Q+M))]

    D = A - t(B)%*%CinvB
    DinvBtCinv = Solve(D,t(CinvB))
    DinvBtCinvKdty = DinvBtCinv%*%Kdty
    
    Phii = do.call("rbind", parallel::mclapply(1:Nd, function(ii){
        c(Cinvi[((ii-1)*(M+1)+1):(ii*(M+1)),]+CinvB[((ii-1)*(M+1)+1):(ii*(M+1)),]%*%DinvBtCinv[,((ii-1)*(M+1)+1):(ii*(M+1))])
    }, mc.cores=ncore)) # Nd x (M+1)*(M+1) # don't forget to multiply s2

    PhiXty = c(Solve(D, Xty) - DinvBtCinvKdty)
    PhiKdty = c(- t(DinvBtCinv)%*%Xty + CinvKdty + CinvB%*%DinvBtCinvKdty)

    sigma2 = (sum(y^2)-sum(c(Xty,Kdty)*c(PhiXty,PhiKdty)))/N
    sigma0 = (sum(y^2)-sum(c(Xty)*c(Solve(A, Xty))))/N

    adata$delta2 = delta2
    adata$delta2d = delta2d

    adata$G0 = G0
    adata$tX = tX
    adata$D = D
    adata$KdtG0base = KdtG0base
    adata$mapid = mapid
    adata$mapid2 = mapid2
    adata$CinvKdtG0base = CinvKdtG0base
    adata$CinvB = CinvB
    adata$DinvBtCinv = DinvBtCinv
    adata$PhiXty = PhiXty
    adata$PhiKdty = PhiKdty
    adata$sigma2 = sigma2
    adata$sigma0 = sigma0
    adata$Phii = Phii

    adata
}

