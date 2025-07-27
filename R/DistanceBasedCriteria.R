dist_to_test<-function(Train, Test, Dst, lambda, C){
  return(max(Dst[rownames(Dst)%in%Train,colnames(Dst)%in%Test]))
}

dist_to_test2<-function(Train, Test, Dst, lambda, C){
  return(mean(Dst[rownames(Dst)%in%Train,colnames(Dst)%in%Test]))
}

neg_dist_in_train<-function(Train, Test, Dst, lambda, C){
  Dt<-Dst[rownames(Dst)%in%Train,colnames(Dst)%in%Train]
  return(-min(Dt[lower.tri(Dt, diag=FALSE)]))
}

neg_dist_in_train2<-function(Train, Test, Dst, lambda, C){
  Dt<-Dst[rownames(Dst)%in%Train,colnames(Dst)%in%Train]
  return(-mean(Dt[lower.tri(Dt, diag=FALSE)]))
}