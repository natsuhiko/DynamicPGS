plotROCandPRC <-
function(pval, flag, names=NULL, col=NULL, lwd=2, main="P-value classifier", add_auc=TRUE, legend_pos="bottomright", roc_xlim=c(0,1), roc_ylim=c(0,1)){
    if(is.vector(pval)) pval <- matrix(pval, ncol=1)
    if(is.data.frame(pval)) pval <- as.matrix(pval)
    if(!is.matrix(pval)) stop("pval must be a vector, matrix, or data.frame.")
    if(length(flag) != nrow(pval)) stop("length(flag) must match nrow(pval).")
    flag0 <- as.integer(flag); K <- ncol(pval)
    if(is.null(names)) names <- colnames(pval)
    if(is.null(names)) names <- paste0("model", seq_len(K))
    if(is.null(col)) col <- seq_len(K)
    col <- rep(col, length.out=K)
    
    get_curve <- function(p, flag){
        ok <- is.finite(p) & !is.na(flag) & p >= 0 & p <= 1 & flag %in% c(0,1)
        p <- p[ok]; flag <- flag[ok]
        P <- sum(flag == 1); N <- sum(flag == 0)
        if(length(p) < 2 || P == 0 || N == 0) return(NULL)
        score <- -log10(pmax(p, .Machine$double.xmin))
        d <- data.frame(score=score, pos=flag==1, neg=flag==0)
        g <- aggregate(cbind(pos,neg) ~ score, d, sum)
        g <- g[order(g$score, decreasing=TRUE),]
        TP <- cumsum(g$pos); FP <- cumsum(g$neg)
        tpr <- TP/P; fpr <- FP/N
        recall <- tpr; precision <- TP/pmax(TP+FP,1)
        roc <- data.frame(fpr=c(0,fpr,1), tpr=c(0,tpr,1))
        pr <- data.frame(recall=c(0,recall), precision=c(1,precision))
        auc_roc <- sum(diff(roc$fpr) * (head(roc$tpr,-1) + tail(roc$tpr,-1))/2)
        auc_pr <- sum(diff(pr$recall) * tail(pr$precision,-1))
        list(roc=roc, pr=pr, auc_roc=auc_roc, auc_pr=auc_pr, n=length(p), n_pos=P, n_neg=N)
    }
    
    curves <- lapply(seq_len(K), function(j) get_curve(pval[,j], flag0))
    names(curves) <- names
    valid <- !sapply(curves, is.null)
    if(!any(valid)) stop("No valid p-value column.")
    
    op <- par(mfrow=c(1,2),mgp=c(2,0.5,0),mar=c(4,4,1,1),cex=1,cex.lab=1,cex.axis=1); on.exit(par(op), add=TRUE)
    
    plot(0,0,type="n",xlim=roc_xlim,ylim=roc_ylim,xlab="False positive rate",ylab="True positive rate",main="",axes=F);axis(1);axis(2,las=2)
    abline(0,1,lty=2)
    for(j in which(valid)) lines(curves[[j]]$roc$fpr, curves[[j]]$roc$tpr, col=col[j], lwd=lwd)
    if(add_auc){
        leg <- sprintf("%s: AUC %.4f", names[valid], sapply(curves[valid], `[[`, "auc_roc"))
        legend(legend_pos, legend=leg, col=col[valid], lwd=lwd, bty="n", cex=0.8)
    } else legend(legend_pos, legend=names[valid], col=col[valid], lwd=lwd, bty="n", cex=0.8)
    
    prev <- mean(flag0 %in% 1, na.rm=TRUE)
    plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="Recall",ylab="Precision",main="",axes=F);axis(1);axis(2,las=2)
    abline(h=prev,lty=2)
    for(j in which(valid)) lines(curves[[j]]$pr$recall, curves[[j]]$pr$precision, col=col[j], lwd=lwd)
    if(add_auc){
        leg <- sprintf("%s: AUPRC %.4f", names[valid], sapply(curves[valid], `[[`, "auc_pr"))
        legend("topright", legend=leg, col=col[valid], lwd=lwd, bty="n", cex=0.8)
    } else legend("topright", legend=names[valid], col=col[valid], lwd=lwd, bty="n", cex=0.8)
    
    summary <- data.frame(name=names[valid], n=sapply(curves[valid], `[[`, "n"), n_pos=sapply(curves[valid], `[[`, "n_pos"), n_neg=sapply(curves[valid], `[[`, "n_neg"), auc_roc=sapply(curves[valid], `[[`, "auc_roc"), auc_pr=sapply(curves[valid], `[[`, "auc_pr"), row.names=NULL)
    invisible(list(summary=summary, curves=curves))
}
