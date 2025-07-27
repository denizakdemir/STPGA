

PEVMEANMM <-
  function(Train,Test, Kinv,K, lambda=1e-5, C=NULL, Vg=NULL, Ve=NULL){
    ntrain<-length(Train)
    M<-matrix(-1/ntrain,ntrain,ntrain)+diag(ntrain)
    X<-matrix(1, ncol=1, nrow=ntrain)
    namesinPforfactorlevels<-rownames(Kinv)
    factorTrain<-factor(as.character(Train), levels=namesinPforfactorlevels)
    factorTest<-factor(as.character(Test), levels=namesinPforfactorlevels)
    Ztrain<-as.matrix(model.matrix(~factorTrain-1))
    Ztest<-as.matrix(model.matrix(~factorTest-1))
   if ((is.null(Vg)||is.null(Ve))){
     if (!is.null(C)){Ztest<-C%*%Ztest}
     PEVMAT<-Ztest%*%solve(crossprod(Ztrain,M%*%Ztrain)+lambda*Kinv)%*%t(Ztest)
   } else {
     bigX<-kronecker(diag(ncol(Vg)), X)
     bigZtrain<-kronecker(diag(ncol(Vg)), Ztrain)
     bigZtest<-kronecker(diag(ncol(Vg)), Ztest)
     if (!is.null(C)){bigZtest<-C%*%bigZtest}
     Ginv<-kronecker(solve(Vg), Kinv)
     Rinv<-kronecker(solve(Ve), diag(ntrain))
     XtRinvX<-t(bigX)%*%Rinv%*%bigX
     ZtRinvX<-t(bigZtrain)%*%Rinv%*%bigX
     ZtRinvZ<-t(bigZtrain)%*%Rinv%*%bigZtrain
     PEVMAT<- (bigZtest%*%solve(ZtRinvZ+Ginv-ZtRinvX%*%solve(XtRinvX)%*%t(ZtRinvX))%*%t(bigZtest))
   }
    return(mean(diag(PEVMAT)))
  }

