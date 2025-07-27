
PEVMAX0 <-
  function(Train,Test, P, lambda=1e-5, C=NULL){
    PTrain<-P[rownames(P)%in%Train,]
    if (!is.null(C)){
    PEVmean<-max(diag(C%*%PTrain%*%solve(crossprod(PTrain)+lambda*diag(ncol(PTrain)),t(C%*%PTrain))))
    } else {PEVmean<-max(diag(PTrain%*%solve(crossprod(PTrain)+lambda*diag(ncol(PTrain)),t(PTrain))))
    }
    return(PEVmean)
  }