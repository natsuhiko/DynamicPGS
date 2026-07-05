#' Step 1 optimisation of the population-level Gaussian-process regression model
#'
#' `GPReg1()` fits the first-stage Gaussian-process regression model for the
#' longitudinal phenotype. This step estimates the population-level smooth
#' trajectory, covariate variance components, the GP length-scale parameter
#' `rho`, and the residual variance.
#'
#' This function should be run after [getData()]. The fitted object is required
#' by [GPReg2()] and downstream association testing.
#'
#' @param adata A `DynamicPGS` object returned by [getData()].
#' @param rho_a Shape-like hyperparameter for the inverse-gamma prior on `rho`.
#'   Internally scaled by the number of individuals.
#' @param rho_b Scale-like hyperparameter for the inverse-gamma prior on `rho`.
#'   Internally scaled by the number of individuals.
#' @param Verbose Logical. If `TRUE`, prints optimisation details.
#' @param Plot Logical. If `TRUE`, plots the lower bound and fitted trajectory
#'   during optimisation.
#' @param MAXITR Integer. Maximum number of optimisation iterations.
#'
#' @return The input `DynamicPGS` object updated with first-stage estimates,
#'   including `rho`, `delta2`, and `lball`.
#'
#' @examples
#' \dontrun{
#' adata <- getData("phenotype.tsv")
#' adata <- GPReg1(adata, Verbose = TRUE)
#' }
#'
#' @export
GPReg1 = function(adata, rho_a=0.01, rho_b=5, Verbose=F, Plot=F, MAXITR=100){

    rho_min=1e-3
    rho_max=3000
    delta_min=1e-8
    delta_max=1000
    jitter=1e-8
    tp = 1 # Titsias correction

    y = adata$y
    x = adata$x
    X = adata$X
    delta2 = adata$delta2
    rho = adata$rho
    ta = adata$ta
    nh = adata$nh

    N = length(y)
    M = length(ta)
    Nd = adata$Nd
    
    rho_a = rho_a*Nd
    rho_b = rho_b*Nd
    
    delta2 = c(delta2, 10)
    nh = c(nh, Context=M)
    P = length(nh)

    if(Verbose) print(nh)

    rho = rho_b / (rho_a + 1)

    cntxt = rep(0:1, c(ncol(X), M)) == 1
    y2 = sum(y^2)

    bfgs_state = NULL
    lball = NULL

    calc_state = function(delta2, rho) {
        if(any(!is.finite(delta2)) || any(delta2 <= 0)) return(NULL)
        if(!is.finite(rho) || rho <= 0) return(NULL)

        delta2 = pmin(pmax(delta2, delta_min), delta_max)
        rho = min(max(rho, rho_min), rho_max)

        Knm = getK(x, ta, rho)
        Kmm = getK(ta, ta, rho)
        Kmm = (Kmm + t(Kmm)) / 2
        Kmm = Kmm # + diag(jitter, nrow(Kmm))

        R = try(chol(Kmm), silent=TRUE)
        if(inherits(R, "try-error")) return(NULL)

        tKnm = t(forwardsolve(t(R), t(Knm)))
        tX = cbind(X, tKnm)
        XtX = crossprod(tX)
        Xty = crossprod(tX, y)
        titsias = tp * (N - sum(tKnm^2))

        A = XtX + diag(1 / delta2[rep(1:P, nh)])
        A = (A + t(A)) / 2

        PhiXty = try(SolveJitter(A, Xty), silent=TRUE)
        if(inherits(PhiXty, "try-error")) return(NULL)

        Phi = try(SolveJitter(A), silent=TRUE)
        if(inherits(Phi, "try-error")) return(NULL)

        sigma2 = (y2 - sum(Xty * PhiXty)) / N
        if(!is.finite(sigma2) || sigma2 <= 0) return(NULL)

        logdetA = getLogDet(A)
        if(!is.finite(logdetA)) return(NULL)

        log_prior_rho = -(rho_a + 1) * log(rho) - rho_b / rho

        lb = - N * log(sigma2) / 2 - sum(nh * log(delta2)) / 2 - logdetA / 2 - delta2[P] * titsias / 2 + log_prior_rho

        if(!is.finite(lb)) return(NULL)

        list(lb=lb,delta2=delta2,rho=rho,Knm=Knm,Kmm=Kmm,R=R,tKnm=tKnm,tX=tX,XtX=XtX,Xty=Xty,A=A,PhiXty=PhiXty,Phi=Phi,sigma2=sigma2,titsias=titsias)
    }

    # initial state
    st = calc_state(delta2, rho)
    if(is.null(st)) stop("Initial state failed.")

    delta2 = st$delta2
    rho = st$rho
    Knm = st$Knm
    Kmm = st$Kmm
    R = st$R
    tKnm = st$tKnm
    tX = st$tX
    XtX = st$XtX
    Xty = st$Xty
    A = st$A
    PhiXty = st$PhiXty
    Phi = st$Phi
    sigma2 = st$sigma2
    titsias = st$titsias
    lball = st$lb
    fail_count = 0

    for(itr in 1:MAXITR){

        if(Plot){
            par(mfcol=c(1,2))
            plot(rev(lball[1:min(length(lball), 100)]),xlab="last 100 iterations",ylab="lower bound")
            x0 = Seq(range(x, na.rm=TRUE))
            Knm0 = getK(x0, ta, rho)
            tKnm0 = t(forwardsolve(t(R), t(Knm0)))
            y0 = tKnm0 %*% tail(PhiXty, M) + PhiXty[1]
            boxplot(y ~ x, at=x0)
            lines(x0, y0, col=2, lwd=3)
        }
        if(Verbose){cat("iteration:", itr, "| lower bound =", lball[1], "| rho =", rho, "\n")}

        # convergence check
        if(itr > 10){ if((lball[1] - lball[10]) / 10 < 1e-5){ break } }

        # gradient for delta2
        dinv = 1 / delta2[rep(1:P, nh)]
        diagXtVinvX = dinv * (delta2[rep(1:P, nh)] - diag(Phi)) * dinv
        XtVinvy2 = c(PhiXty * dinv)^2

        g = unlist(lapply(split(XtVinvy2 / sigma2 / 2 - diagXtVinvX / 2, rep(1:P, nh)),sum))

        g[P] = g[P] - titsias / 2
        g = g * delta2

        # gradient for rho
        pKmm = getKprime(ta, ta, rho)
        pKnm = getKprime(x, ta, rho)

        tXrho = cbind(X, Knm)

        Phiinvrho =
            crossprod(tXrho) +
            dbind(
                diag(1 / delta2[rep(seq(P-1), nh[-P])]),
                Kmm / delta2[P]
            )

        Phiinvrho = (Phiinvrho + t(Phiinvrho)) / 2

        Phirho = SolveJitter(Phiinvrho)
        PhiXtrho = SolveJitter(Phiinvrho, t(tXrho))
        PhiXtyrho = c(SolveJitter(Phiinvrho, crossprod(tXrho, y)))

        Kmm_inv = SolveJitter(Kmm)
        KmminvKmn = Kmm_inv %*% t(Knm)

        grho =
            - sum(t(PhiXtrho)[, cntxt] * pKnm) +
            sum((delta2[P] * Kmm_inv - Phirho[cntxt, cntxt]) *
                    pKmm / delta2[P]) / 2

        grho =
            grho +
            sum(
                PhiXtyrho[cntxt] *
                    (
                        c(crossprod(y, pKnm)) -
                        c((t(PhiXtyrho) %*% t(tXrho)) %*% pKnm)
                    )
            ) / sigma2 -
            sum(t(pKmm / delta2[P] * PhiXtyrho[cntxt]) *
                    PhiXtyrho[cntxt]) / sigma2 / 2

        grho =
            grho +
            tp * delta2[P] * sum(KmminvKmn * t(pKnm)) -
            tp * delta2[P] *
                sum((KmminvKmn %*% t(KmminvKmn)) * pKmm) / 2

        # Convert d/d rho to d/d log(rho),
        # then add IG prior gradient:
        # d log p(rho) / d log rho = -(rho_a+1) + rho_b/rho
        grho = grho * rho - (rho_a + 1) + rho_b / rho

        theta = log(c(delta2, rho))
        grad_max = c(g, grho)

        tmp = BFGS(
            theta,
            grad_max,
            state=bfgs_state,
            init_scale=1e-3,
            max_step=2
        )

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

        # conservative step clipping on log scale
        ss[ss > 0.1] = 0.1
        ss[ss < -0.1] = -0.1

        accepted = FALSE

        for(sign_dir in c(1, -1)){
            for(r in 0:20){
                theta1 = theta + sign_dir * ss / 2^r

                theta1[1:P] =
                    pmin(pmax(theta1[1:P], log(delta_min)), log(delta_max))

                theta1[P+1] =
                    pmin(pmax(theta1[P+1], log(rho_min)), log(rho_max))

                delta21 = exp(theta1[1:P])
                rho1 = exp(theta1[P+1])

                st1 = calc_state(delta21, rho1)

                if(is.null(st1)){
                    if(Verbose) cat("sign =", sign_dir, "r =", r, "candidate failed\n")
                    next
                }

                lb1 = st1$lb

                if(Verbose){
                    cat("sign =", sign_dir,
                        "r =", r,
                        "| lb1 =", lb1,
                        "| diff =", lb1 - lball[1],
                        "| rho1 =", rho1, "\n")
                }

                if(is.finite(lb1) && lb1 > lball[1] - 1e-8){
                    delta2 = st1$delta2
                    rho = st1$rho
                    Knm = st1$Knm
                    Kmm = st1$Kmm
                    R = st1$R
                    tKnm = st1$tKnm
                    tX = st1$tX
                    XtX = st1$XtX
                    Xty = st1$Xty
                    A = st1$A
                    PhiXty = st1$PhiXty
                    Phi = st1$Phi
                    sigma2 = st1$sigma2
                    titsias = st1$titsias

                    lball = c(lb1, lball)

                    if(sign_dir == 1){
                        bfgs_state = bfgs_state_candidate
                    }else{
                        bfgs_state = NULL
                    }

                    accepted = TRUE
                    fail_count = 0

                    if(Verbose) print(c(delta2, rho))
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
                cat("No acceptable step found at iteration", itr, "\n")
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

    adata$P = P
    adata$nh = nh
    adata$rho = rho
    adata$A = A
    adata$sigma2 = sigma2
    adata$delta2 = delta2
    adata$lball = lball
    adata$PhiXty = PhiXty

    return(adata)
}
