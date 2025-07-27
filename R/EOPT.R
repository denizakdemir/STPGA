


EOPT <-
  function(Train,Test=NULL, P, lambda=1e-5,C=NULL){
    PTrain<-P[rownames(P)%in%Train,]
    if (!is.null(C)){
    svdD<-svd(C%*%solve(crossprod(PTrain)+lambda*diag(ncol(PTrain)), t(C)), nu=0,nv=0)
    } else {
      svdD<-svd(crossprod(PTrain)+lambda*diag(ncol(PTrain)), nu=0,nv=0)
      
    }
    return(max(1/svdD$d))
  }