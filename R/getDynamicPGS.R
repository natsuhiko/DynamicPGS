
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
