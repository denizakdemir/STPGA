makeonecross <-
  function(x1,x2,Candidates,mutprob=.5, mutintensity=2){
    n1<-length(unlist(x1))
    n2<-length(unlist(x2))
    n<-min(c(n1,n2))
    x1x2<-union(unlist(x1),unlist(x2))
    cross<-sample(x1x2,n, replace=F)
    #x1x2<-c(unlist(x1),unlist(x2))
    #cross<-c()
    #while(length(unique(cross))<n){cross<-c(cross,sample(x1x2,n, replace=T))}
    #cross<-unique(cross)[1:n]
    randnum<-runif(1)
    if (randnum<mutprob){
      setdiffres<-setdiff(Candidates,cross)
      ntoreplace<-min(c(rpois(1,min(mutintensity, n)),n,length(setdiffres)))
      cross[sample(1:n,ntoreplace)]<-sample(setdiffres,ntoreplace)
    }
    return(sort(cross))
  }
