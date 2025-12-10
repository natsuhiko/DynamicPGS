getLogDet <-
function(x){
    x=(x+t(x))/2;
    eval = try(eigen(x,T)[[1]]);
    if(is.character(eval)){
        return(NA)
    }
    sum(log(eval))
}
getK <-
function(x, ta, rho){
    X = outer(x, ta, "-")
    exp(-X^2/rho)
}
Solve <-
function(x,y=diag(nrow(x))){
    r  = chol((x+t(x))/2)
    ra = forwardsolve(t(r),y)
    backsolve(r, ra)
}
