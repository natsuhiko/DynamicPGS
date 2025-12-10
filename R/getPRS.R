
getDynamicPRS = function(xstar=0:54, PGEN="/path/to/your/pgen/dir/", ta=0:11*5, rho=500){
    load("Data/beta0sigma2.Rbin")
    load("Data/BetaMatrix.Rbin")
    load("Data/StderrBetaMatrix.Rbin")
    load("Data/proxy.Rbin")

    source("R/getK.R")

    library(pgenlibr)
    
    L=length(unique(proxy[,1]))
    vexist=rep(F,nrow(proxy))
    vid=rep(NA,nrow(proxy))
    for(CHR in c(1:22,"X")){
	cat(paste("Searching proxy variants on chromosome",CHR,"...\n"))
	if(file.exists(paste(PGEN,"/chr",CHR,".pvar",sep=""))){
            pvar = read.table(paste(PGEN,"/chr",CHR,".pvar",sep=""),header=F,stringsAsFactors=F,comment="#")
            names(pvar) = c("CHROM", "POS", "ID", "REF", "ALT", "INFO")
    	    vexist = vexist | (paste(proxy[,2],proxy[,3],proxy[,5],proxy[,6])%in%paste(pvar[,1],pvar[,2],pvar[,4],pvar[,5]))
    	    tmp = pvar$ID[match(paste(proxy[,2],proxy[,3],proxy[,5],proxy[,6]),paste(pvar[,1],pvar[,2],pvar[,4],pvar[,5]))]
	    vid[!is.na(tmp)] = tmp[!is.na(tmp)]
	}
    }
    map = match(unique(proxy[vexist,1]),proxy[vexist,1])
    proxy = proxy[vexist,][map,]
    B = B[vexist,][map,]
    Sinv = Sinv[vexist,][map,]
    vid = vid[vexist][map]
    beta0 = beta0sigma2[[1]]
    sigma2 = beta0sigma2[[2]]
    cat(paste(nrow(proxy),"proxy variants out of 322 index variants were found...\n")) 
     
    L = length(unique(proxy[,1]))
    M = length(ta)
    Knm = getK(xstar, ta, rho)
    Kmm = getK(ta, ta, rho)
    R = chol(Kmm)
    tKnm = t(forwardsolve(t(R),t(Knm)))
    G0 = cbind(tKnm, 1) # N x M+1

    cat("Computing dynamic PRS at x*\n")
    psam <- read.table(paste(PGEN,"/chr22.psam",sep=""),header=F)
    E=V=0
    for(l in 1:L){
        CHR=gsub("chr","",proxy[l,2])
        pvar <- pgenlibr::NewPvar(paste(PGEN,"/chr",CHR,".pvar",sep=""))
        pgen <- pgenlibr::NewPgen(paste(PGEN,"/chr",CHR,".pgen",sep=""),pvar=pvar)
        gl <- pgenlibr::ReadList(pgen, GetVariantsById(pvar, vid[l]))
        gl = gl - mean(gl)
        E=E+G0%*%B[l,]%*%t(gl)
        V=V+sigma2*colSums(solve(matrix(Sinv[l,],M+1),t(G0))*t(G0))%*%t(gl^2)
    }
    list(xstar=xstar, avg=G0%*%beta0, E=E, V=V)
}
