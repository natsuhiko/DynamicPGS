
gpreg2 = function(adata, delta2d0=c(0.94,1.7), Verbose=F, Plot=F, ncore=max(1, parallel::detectCores()-1), MAXITR=100){
    library(parallel)
    
    if(is.null(adata$rho)){
        stop("Step 1 optimization has not been performed yet...aborted")
    }

    y = adata$y
    x = adata$x
    iid = adata$iid
    uiid = unique(iid)
    fid = adata$fid
    X = adata$X
    delta2 = adata$delta2
    rho = adata$rho
    ta = adata$ta
    nh = adata$nh
    Lmat = adata$Lmat[match(iid, unique(iid)),]
    if(sum(Lmat$IID == iid) != length(iid)){
        return("Lmat is inappropriate!")
    }else{
        Lmat = as.matrix(Lmat[,-c(1:2)])
    }

    N = length(y)
    M = length(ta)
    Nd = length(table(iid))
    Nf = length(table(fid))

    cat("Nd: "); cat(Nd); cat("\n")
    cat("Nf: "); cat(Nf); cat("\n")
    cat("N: "); cat(N); cat("\n")
    cat("M: "); cat(M); cat("\n")

    delta2d = delta2d0
    P = length(nh)
    Q = ncol(X)
    cat("Q: "); cat(Q); cat("\n")

    Knm = getK(x, ta, rho)
    Kmm = getK(ta, ta, rho)
    R = chol(Kmm)
    tKnm = t(forwardsolve(t(R), t(Knm)))
    tX = cbind(X, tKnm)
    XtX = crossprod(tX)
    Xty = crossprod(tX, y)
    titsias = max(0, N - sum(tKnm^2))
    y2 = sum(y^2)

    bfgs_state = NULL
    lball = NULL
    idx_list = split(seq_len(N), fid)

    one_id_fun = function(ii, delta2din){
        tx = tX[ii,,drop=FALSE]
        n1 = length(unique(iid[ii]))
        tk = C1prior = NULL

        for(j in 1:n1){
            tk = cbind(tk, cbind(tKnm[ii,,drop=FALSE], 1) * Lmat[ii,j])
            C1prior = dbind(C1prior, diag(rep(1 / delta2din, c(M,1))))
        }

        K = (M + 1) * n1
        y1 = y[ii]
        C1 = C1prior + crossprod(tk)
        ldc1 = getLogDet(C1)
        B1 = crossprod(tk, tx)
        tky = crossprod(tk, y1)
        res = Solve(C1, cbind(B1, diag(K), tky))

        cbind(
            res[,1:(Q+M),drop=FALSE],
            diag(res[,(Q+M+1):(Q+M+K),drop=FALSE]),
            res[,Q+M+K+1],
            c(ldc1, rep(0, K-1)),
            B1,
            tky
        )
    }

    eval_state = function(delta2in, delta2din){
        if(any(!is.finite(delta2in)) || any(delta2in <= 0)) return(NULL)
        if(any(!is.finite(delta2din)) || any(delta2din <= 0)) return(NULL)

        A1 = XtX + diag(1 / delta2in[rep(1:P, nh)])
        A1 = (A1 + t(A1)) / 2

        Ctmp = try(
            do.call("rbind", parallel::mclapply(idx_list, one_id_fun, delta2din=delta2din, mc.cores=ncore)),
            silent=TRUE
        )
        if(inherits(Ctmp, "try-error")) return(NULL)

        diagCinv1 = Ctmp[,Q+M+1]
        CinvKdty1 = Ctmp[,Q+M+2]
        logDetC1 = Ctmp[,Q+M+3]
        B1 = Ctmp[,(Q+M+3)+(1:(Q+M)),drop=FALSE]
        Kdty1 = Ctmp[,(Q+M+3)+(Q+M)+1]
        CinvB1 = Ctmp[,1:(Q+M),drop=FALSE]

        out = try({
            D1 = A1 - crossprod(B1, CinvB1)
            D1 = (D1 + t(D1)) / 2
            DinvBtCinv1 = Solve(D1, t(CinvB1))
            DinvBtCinvKdty1 = DinvBtCinv1 %*% Kdty1
            PhiXty1 = c(Solve(D1, Xty) - DinvBtCinvKdty1)
            PhiKdty1 = c(-t(DinvBtCinv1) %*% Xty + CinvKdty1 + CinvB1 %*% DinvBtCinvKdty1)
            sigma21 = (y2 - sum(c(Xty, Kdty1) * c(PhiXty1, PhiKdty1))) / N
            if(!is.finite(sigma21) || sigma21 <= 0) stop("bad sigma2")

            logDetD1 = getLogDet(D1)
            if(!is.finite(logDetD1)) stop("bad logDetD")

            lb1 = -N * log(sigma21) / 2 -
                sum(nh * log(delta2in)) / 2 -
                sum(Nd * c(M,1) * log(delta2din)) / 2 -
                sum(logDetC1) / 2 -
                logDetD1 / 2 -
                (delta2in[P] + delta2din[1]) * titsias / 2

            list(
                lb=lb1,
                A=A1,
                diagCinv=diagCinv1,
                CinvKdty=CinvKdty1,
                logDetC=logDetC1,
                B=B1,
                Kdty=Kdty1,
                CinvB=CinvB1,
                D=D1,
                DinvBtCinv=DinvBtCinv1,
                PhiXty=PhiXty1,
                PhiKdty=PhiKdty1,
                sigma2=sigma21
            )
        }, silent=TRUE)

        if(inherits(out, "try-error")) return(NULL)
        if(!is.finite(out$lb)) return(NULL)
        out
    }

    st = eval_state(delta2, delta2d)
    if(is.null(st)) stop("Initial state failed.")

    A = st$A
    diagCinv = st$diagCinv
    CinvKdty = st$CinvKdty
    logDetC = st$logDetC
    B = st$B
    Kdty = st$Kdty
    CinvB = st$CinvB
    D = st$D
    DinvBtCinv = st$DinvBtCinv
    PhiXty = st$PhiXty
    PhiKdty = st$PhiKdty
    sigma2 = st$sigma2
    lball = st$lb
    fail_count = 0

    for(itr in 1:MAXITR){
        if(Plot){
            par(mfcol=c(1,2))
            plot(rev(lball[1:min(length(lball),100)]), xlab="last 100 iterations", ylab="lower bound")
            x0 = Seq(range(x, na.rm=T))
            Knm0 = getK(x0, ta, rho)
            tKnm0 = t(forwardsolve(t(R), t(Knm0)))
            y0 = tKnm0 %*% tail(PhiXty, M) + PhiXty[1]
            boxplot(y ~ x, at=x0)
            lines(x0, y0, col=2, lwd=3)
        }

        cat("iteration: "); cat(itr); cat(" | lower bound = "); cat(lball[1]); cat("\n")

        if(itr > 20 && length(lball) >= 10){
            if((lball[1] - lball[10]) / 10 < 1e-3){
                break
            }
        }

        Dinv = Solve(D)

        dKd = rep(rep(delta2d, c(M,1)), Nd)
        dinv2 = 1 / dKd
        diagKdtVinvKd = dinv2 * (dKd - (diagCinv + colSums(DinvBtCinv * t(CinvB)))) * dinv2
        KdVinvy2 = c(PhiKdty * dinv2)^2

        g2 = rowSums(matrix(KdVinvy2 / sigma2 / 2 - diagKdtVinvKd / 2, M+1))
        g2 = c(sum(g2[1:M]), g2[M+1])
        g2[1] = g2[1] - titsias / 2
        g2 = g2 * delta2d

        theta = log(delta2d)
        grad_max = g2

        tmp = BFGS(theta, grad_max, state=bfgs_state, init_scale=1e-3, max_step=2)
        ss = tmp$direction
        bfgs_state_candidate = tmp$state

        dirderiv = sum(grad_max * ss)
        if(!is.finite(dirderiv) || dirderiv <= 0){
            if(Verbose){
                cat("BFGS direction is not ascent; use gradient direction.\n")
                cat("directional derivative =", dirderiv, "\n")
            }
            ss = grad_max
            ss = ss / max(1, max(abs(ss))) * 0.05
            bfgs_state_candidate = NULL
        }

        ss[ss > 0.1] = 0.1
        ss[ss < -0.1] = -0.1

        accepted = FALSE

        for(sign_dir in c(1, -1)){
            for(r in 0:20){
                delta2d1 = exp(log(delta2d) + sign_dir * ss / 2^r)

                st1 = eval_state(delta2, delta2d1)
                if(is.null(st1)) next

                lb1 = st1$lb

                if(Verbose){
                    cat("sign =", sign_dir, "r =", r, "lb1 =", lb1, "diff =", lb1 - lball[1], "\n")
                }

                if(is.finite(lb1) && lb1 > lball[1] - 1e-8){
                    delta2d = delta2d1

                    A = st1$A
                    diagCinv = st1$diagCinv
                    CinvKdty = st1$CinvKdty
                    logDetC = st1$logDetC
                    B = st1$B
                    Kdty = st1$Kdty
                    CinvB = st1$CinvB
                    D = st1$D
                    DinvBtCinv = st1$DinvBtCinv
                    PhiXty = st1$PhiXty
                    PhiKdty = st1$PhiKdty
                    sigma2 = st1$sigma2

                    lball = c(lb1, lball)

                    if(sign_dir == 1){
                        bfgs_state = bfgs_state_candidate
                    }else{
                        bfgs_state = NULL
                    }

                    accepted = TRUE
                    fail_count = 0

                    if(Verbose) print(delta2d)
                    break
                }
            }
            if(accepted) break
        }

        if(!accepted){
            fail_count = fail_count + 1
            bfgs_state = NULL
            grad_norm = max(abs(grad_max))

            if(Verbose){
                cat("No acceptable step at iteration", itr, "\n")
                cat("max abs gradient =", grad_norm, "\n")
                cat("fail_count =", fail_count, "\n")
            }

            if(grad_norm < 1e-4){
                if(Verbose) cat("Treat as converged because gradient is small.\n")
                break
            }

            if(fail_count >= 3){
                if(Verbose) cat("Stopping because no acceptable step was found repeatedly.\n")
                break
            }

            next
        }
    }
    adata$delta2d = delta2d
    adata$lball_d = lball
    adata = prep_assoc(adata)
    adata
    
}
