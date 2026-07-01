getZd <- function(did){
    did <- factor(did, levels=unique(did))
    Matrix::sparseMatrix(
        i = seq_along(did),
        j = as.integer(did),
        x = 1,
        dims = c(length(did), nlevels(did))
    )
}

getKnmd <- function(Knm, did){
    Knm <- cbind(1, Knm)
    did <- factor(did, levels=unique(did))
    n <- nrow(Knm)
    M <- ncol(Knm)
    Nd <- nlevels(did)
    Matrix::sparseMatrix(
        i = rep(seq_len(n), M),
        j = as.integer(did) + rep(seq(0, by=Nd, length.out=M), each=n),
        x = as.double(Knm),
        dims = c(n, Nd*M)
    )
}



getK = function(x, ta, rho){
    X = outer(x, ta, "-")
    exp(-X^2/rho)
}

getKprime = function(x, ta, rho){
    X = outer(x, ta, "-")
    exp(-X^2/rho)*(X^2)/rho/rho
}

Seq=function(x){seq(x[1],x[2])}

BFGS = function(theta, grad_max, state = NULL,
                     init_scale = 1e-3,
                     eps = 1e-10,
                     max_step = 2) {
    # theta: log(c(delta2, rho))
    # grad_max: d lower_bound / d theta
    # internally use grad_min = -grad_max

    grad_min = -grad_max
    p = length(theta)

    if (is.null(state)) {
        H = diag(init_scale, p)
        direction = -H %*% grad_min  # = H %*% grad_max

        state = list(
            theta_prev = theta,
            grad_prev = grad_min,
            H = H
        )

        direction = as.numeric(direction)
        direction = pmax(pmin(direction, max_step), -max_step)

        return(list(direction = direction, state = state))
    }

    s = theta - state$theta_prev
    y = grad_min - state$grad_prev

    H = state$H
    sy = sum(s * y)

    if (is.finite(sy) && sy > eps) {
        rho_bfgs = 1 / sy
        I = diag(p)

        V = I - rho_bfgs * tcrossprod(s, y)
        H = V %*% H %*% t(V) + rho_bfgs * tcrossprod(s, s)

        H = (H + t(H)) / 2
    } else {
        # curvature condition failed
        # reset inverse Hessian
        H = diag(init_scale, p)
    }

    direction = -H %*% grad_min  # minimization descent = lower-bound ascent
    direction = as.numeric(direction)

    # avoid huge log-scale jumps
    direction = pmax(pmin(direction, max_step), -max_step)

    state = list(
        theta_prev = theta,
        grad_prev = grad_min,
        H = H
    )

    list(direction = direction, state = state)
}

LBFGS = function(X, G, eps = 1e-12) {
    M = ncol(X)
    if (M < 2) { return(-G[,1]) }
    S = X[, 1:(M-1), drop=FALSE] - X[, 2:M, drop=FALSE]
    Y = G[, 1:(M-1), drop=FALSE] - G[, 2:M, drop=FALSE]
    q = G[,1]
    a = rep(0, M-1)
    flag = rep(TRUE, M-1)
    sy = rep(NA_real_, M-1)
    for (i in 1:(M-1)) {
        sy[i] = sum(S[,i] * Y[,i])
        if (is.na(sy[i]) || sy[i] <= eps) {
            flag[i] = FALSE
            next
        }
        a[i] = sum(S[,i] * q) / sy[i]
        q = q - a[i] * Y[,i]
    }
    effi = rev(seq_len(M-1)[flag])
    if (length(effi) == 0) { return(-G[,1]) }
    i0 = min(effi)
    yy = sum(Y[,i0]^2)
    gammak = sy[i0] / yy

    if (is.na(gammak) || gammak <= 0 || yy <= eps) { gammak = 1 }
    z = gammak * q
    for (i in effi) {
        bi = sum(Y[,i] * z) / sy[i]
        z = z + S[,i] * (a[i] - bi)
    }
    return(-z)
}


Solve = function(x,y=diag(nrow(x))){
    r  = chol((x+t(x))/2)
    ra = forwardsolve(t(r),y)
    backsolve(r, ra)
}

dbind=function(...){
    dims=NULL

    lis = list(...)
    lis = lis[!unlist(lapply(lis,is.null))]
    
    
    N = length(lis)
    for(i in 1:length(lis)){
        x = dim(lis[[i]])
        if(is.null(x)){
            dims = rbind(dims, c(length(lis[[i]]),1))
        }else{
            dims = rbind(dims, x)
        }
    }
    A = array(0, apply(dims,2,sum))
    cdims = apply(rbind(c(0,0),dims),2,cumsum)
    for(i in 1:N){
        if(length(lis[[i]])>0){ A[(cdims[i,1]+1):cdims[i+1,1], (cdims[i,2]+1):cdims[i+1,2]]= lis[[i]] }
    }
    A
}

getLogDet=function(x){
    x=(x+t(x))/2;
    eval = try(eigen(x,T)[[1]]);
    if(is.character(eval)){
        return(NA)
    }
    sum(log(eval))
}

SolveJitter = function(A, B = diag(nrow(A)), jitter = 1e-8, max_try = 8) {
    A = (A + t(A)) / 2
    scaleA = mean(diag(A))
    if(!is.finite(scaleA) || scaleA <= 0) scaleA = 1

    for(k in 0:max_try) {
        jj = jitter * 10^k * scaleA
        R = try(chol(A + diag(jj, nrow(A))), silent = TRUE)
        if(!inherits(R, "try-error")) {
            return(backsolve(R, forwardsolve(t(R), B)))
        }
    }

    stop("SolveJitter failed")
}





Alpha=function(col, alpha=0.1*255){r=col2rgb(col); rgb(r[1,]/255,r[2,]/255,r[3,]/255,alpha/255)}


Seq <- function(x) seq(x[1], x[2], by=1)
