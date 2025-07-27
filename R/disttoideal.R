disttoideal<-function(X){
  Xc<-apply(X,2,FUN=function(x){(x-min(x))/(max(x)-min(x))})
  return(as.matrix(dist(rbind(matrix(0,ncol=ncol(Xc), nrow=1),Xc)))[-1,1])
}