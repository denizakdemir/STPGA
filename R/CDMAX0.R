

CDMAX0 <-
  function(Train,Test, P, lambda=1e-5, C=NULL){
    PTrain<-P[rownames(P)%in%Train,]
    if (!is.null(C)){ PTrain<-C%*%PTrain}
    CDmean<-max(diag(PTrain%*%solve(crossprod(PTrain)+lambda*diag(ncol(PTrain)),t(PTrain)))/diag(tcrossprod(PTrain)))
    return(CDmean)
  }
