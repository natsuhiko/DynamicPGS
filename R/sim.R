getMockData = function(adata, Nd=1000){
    
    k = simulate_king(Nd)
    uiid = k$iid
    x = simulate_x_from_xlist(split(adata$x, adata$iid), Nd, iid=uiid)
    N = nrow(x)
    
    pc = data.frame(pc1=rnorm(Nd),pc2=rnorm(Nd))[match(x$IID,uiid),]
    sex = c("F","M")[(runif(Nd)>0.5)+1][match(x$IID,uiid)]
    
    M = adata$M
    ta = adata$ta
    rho = adata$rho
    beta0 = c(tail(adata$PhiXty,M), adata$PhiXty[1])
    Knm = getK(x$x, ta, rho)
    Kmm = getK(ta, ta, rho)
    R = chol(Kmm)
    tKnm = cbind(t(forwardsolve(t(R),t(Knm))),1)
    
    x=cbind(x,y=rep(NA,N))
    getData(x, Cov=cbind(pc,sex), king=k$king, inducing_points=ta)
}

simulate_king <- function(N=1000, max_3rd_degree_size=3, seed=1){
    set.seed(seed)
    iid <- sprintf("ID%06d", seq_len(N))
    fid <- sprintf("F%06d", seq_len(N))

    kin_range <- list(
        PO  = c(0.1768, 0.3536),   # parent-offspring / duplicate-ish upper avoided
        FS  = c(0.1768, 0.3536),   # full siblings
        S2  = c(0.0884, 0.1768),   # second-degree
        S3  = c(0.0442, 0.0884)    # third-degree
    )

    res <- list()
    used <- rep(FALSE, N)
    k <- 1

    while(any(!used)){
        avail <- which(!used)
        m <- sample(seq_len(max_3rd_degree_size), 1)
        m <- min(m, length(avail))
        members <- sample(avail, m)
        used[members] <- TRUE

        if(m >= 2){
            pairs <- t(combn(members, 2))
            keep <- runif(nrow(pairs)) < runif(1, 0.4, 0.9)
            pairs <- pairs[keep,,drop=FALSE]

            if(nrow(pairs) > 0){
                deg <- sample(c("PO","FS","S2","S3"), nrow(pairs), replace=TRUE, prob=c(0.15,0.20,0.30,0.35))
                kin <- mapply(function(d) runif(1, kin_range[[d]][1], kin_range[[d]][2]), deg)
                ibs0 <- ifelse(deg=="PO", 0,
                         ifelse(deg=="FS", runif(length(deg),0.001,0.01),
                         ifelse(deg=="S2", runif(length(deg),0.005,0.03),
                                        runif(length(deg),0.01,0.05))))
                hethet <- runif(length(deg), 0.08, 0.18)
                ibs0_seg <- runif(length(deg), 0, 0.02)

                res[[k]] <- data.frame(
                    FID1=fid[pairs[,1]], ID1=iid[pairs[,1]],
                    FID2=fid[pairs[,2]], ID2=iid[pairs[,2]],
                    N_SNP=sample(500000:700000, nrow(pairs), replace=TRUE),
                    HetHet=hethet,
                    IBS0=ibs0,
                    Kinship=kin,
                    Error=0,
                    stringsAsFactors=FALSE
                )
                k <- k+1
            }
        }
    }

    kin0 <- if(length(res)) do.call(rbind, res) else data.frame()
    kin0 <- kin0[kin0$Kinship >= 0.0442,,drop=FALSE]
    rownames(kin0) <- NULL
    list(king=kin0, iid=iid)
}

simulate_x_from_xlist <- function(xlist, N=1000, iid=NULL, seed=1, jitter=0, round_digits=1){
    set.seed(seed)
    xlist <- lapply(xlist, function(z) sort(round(as.numeric(z[is.finite(z)]), round_digits)))
    xlist <- xlist[lengths(xlist) > 0]
    if(length(xlist)==0) stop("xlist has no valid month-age vectors")
    if(is.null(iid)) iid <- sprintf("SIM%06d", seq_len(N))
    if(length(iid) != N) stop("length(iid) must be equal to N")
    out <- lapply(iid, function(id){
        x <- sample(xlist, 1)[[1]]
        if(jitter > 0) x <- pmax(0, x + rnorm(length(x), 0, jitter))
        data.frame(IID=id, x=x)
    })
    out <- do.call(rbind, out); 
    rownames(out) <- NULL;
    out
}
