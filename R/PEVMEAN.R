PEVMEAN <-
  function(Train,Test, P, lambda=1e-5, C=NULL){
    PTrain<-P[rownames(P)%in%Train,]
    PTest<-P[rownames(P)%in%Test,]
    if (length(Test)==1){PTest=matrix(PTest, nrow=1)}
    if (!is.null(C)){ PTest<-C%*%PTest}
    PEVmean<-mean(diag(PTest%*%solve(crossprod(PTrain)+lambda*diag(ncol(PTrain)),t(PTest))))
    return(PEVmean)
  }



