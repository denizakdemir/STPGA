

GAUSSMEANMM <-
  function(Train,Test, Kinv, K, lambda=1e-5,C=NULL, Vg=NULL, Ve=NULL){
    ntrain<-length(Train)
    M<-matrix(-1/ntrain,ntrain,ntrain)+diag(ntrain)
    namesinPforfactorlevels<-rownames(Kinv)
    factorTrain<-factor(as.character(Train), levels=namesinPforfactorlevels)
    factorTest<-factor(as.character(Test), levels=namesinPforfactorlevels)
    Ztrain<-as.matrix(model.matrix(~factorTrain-1))
    Ztest<-as.matrix(model.matrix(~factorTest-1))
    if (!is.null(C)){Ztest<-C%*%Ztest}
    Dmat<-Ztest%*%(K)%*%t(Ztest)-Ztest%*%(K)%*%t(Ztrain)%*%solve(Ztrain%*%K%*%t(Ztrain)+lambda*diag(nrow(Ztrain)))%*%Ztrain%*%(K)%*%t(Ztest)
    return(mean(diag(Dmat)))
  }

