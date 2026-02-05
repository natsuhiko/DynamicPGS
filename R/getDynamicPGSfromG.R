# VID = Locus:CHR:POS:REF:ALT
getDynamicPGSfromG <-
function(xstar=0:54, G, VID, ta=0:11*5, rho=500){
    load("Data/beta0sigma2.Rbin")
    load("Data/BetaMatrix.Rbin")
    load("Data/StderrBetaMatrix.Rbin")
    load("Data/proxy.Rbin")

    source("R/getK.R")

    map = match(VID,paste(proxy[,1],proxy[,2],proxy[,3],proxy[,5],proxy[,6],sep=":"))
    proxy = proxy[map,]
    B = B[map,]
    Sinv = Sinv[map,]
    beta0 = beta0sigma2[[1]]
    sigma2 = beta0sigma2[[2]]
     
    L = length(unique(proxy[,1]))
    M = length(ta)
    Knm = getK(xstar, ta, rho)
    Kmm = getK(ta, ta, rho)
    R = chol(Kmm)
    tKnm = t(forwardsolve(t(R),t(Knm)))
    G0 = cbind(tKnm, 1) # N x M+1
    
    cat("Computing dynamic PRS at x*\n")
    E=V=0
    for(l in 1:L){
        gl <- G[l,]
        if(sum(is.na(gl))==0){
            gl = gl - mean(gl)
            E=E+G0%*%B[l,]%*%t(gl)
            V=V+sigma2*colSums(solve(matrix(Sinv[l,],M+1),t(G0))*t(G0))%*%t(gl^2)
        }
    }
    list(xstar=xstar, avg=G0%*%beta0, E=E, SE=sqrt(V), sigma2=sigma2)
}
