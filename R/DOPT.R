


DOPT <-
  function(Train,Test=NULL, P, lambda=1e-5, C=NULL){
    PTrain<-P[rownames(P)%in%Train,]
    if (!is.null(C)){
    D<-determinant(C%*%solve(crossprod(PTrain)+lambda*diag(ncol(PTrain)), t(C)), logarithm=TRUE)$modulus
        } else {
    D<--determinant(crossprod(PTrain)+lambda*diag(ncol(PTrain)), logarithm=TRUE)$modulus
      
    }
    return(D)
  }