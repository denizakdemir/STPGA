GenAlgForSubsetSelectionMO<-function(Pcs=NULL,Dist=NULL, Kernel=NULL, Candidates,Test,ntoselect,selectionstats,selectionstatstypes,plotdirections, npopGA=100,  mutprob=.8, mutintensity=1, nitGA=500,lambda=1e-6, plotiters=FALSE, mc.cores=1, InitPop=NULL, C=NULL, axes.labels=NULL){
    if(is.null(C)){
    C1<-diag(length(Test))
    C2<-diag(ncol(Pcs))
    }
  
  if (!is.null(Kernel)){Kinv<-chol2inv(chol(Kernel))}
  nstats<-length(selectionstats)
  GAobjfunc<-function(Train){
    outlist<-lapply(selectionstats,FUN = function(x) {
      if (x%in%c("DOPT", "AOPT", "EOPT")){C<-C2} else{C<-C1}
      if (selectionstatstypes[which(selectionstats%in%x)]%in%"Pcs"){
      
        return(do.call(x, list(Train, Test, Pcs, lambda, C)))
      } else if (selectionstatstypes[which(selectionstats%in%x)]%in%"Kernel"){
      
        return(do.call(x, list(Train, Test, Kinv,Kernel, lambda, C)))
      } else if (selectionstatstypes[which(selectionstats%in%x)]%in%"Dist"){
      
        return(do.call(x, list(Train, Test, Dist, lambda, C)))
      }
    })
    return(outlist)
  }
  
  
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
  
  GenerateCrossesfromElites<-
    function(Elites, Candidates, npopGA, mutprob,mc.cores=1, mutintensity=1, memoryfortabu=NULL){
      
      newcrosses<-mclapply(1:npopGA, FUN=function(x){
        x1<-Elites[[sample(1:length(Elites),1)]]
        x2<-Elites[[sample(1:length(Elites),1)]]
        out<-makeonecross(x1=x1,x2=x2,Candidates=Candidates,mutprob=mutprob, mutintensity=mutintensity)
        return(out)
      }, mc.cores=mc.cores,  mc.set.seed = T)
      return(newcrosses)
    }
  
  
  GAfunc<-function (ntoselect, npopGA,  mutprob, niterations) 
  {

    lengthfrontier<-0
    while(lengthfrontier<5){
    InitPop <- mclapply(1:npopGA, function(x) {
      return(sample(Candidates, ntoselect, replace=FALSE))
    }, mc.cores=mc.cores)
    
    InitPopFuncValues <- matrix(as.numeric(unlist(mclapply(InitPop,FUN = function(x) {
      GAobjfunc(x)}, mc.cores = mc.cores, mc.preschedule = T))), ncol=  nstats, byrow=TRUE)
    
    bestsolsfortraits<-apply(InitPopFuncValues, 2, which.min)
    
    frontier3<- which(!emoa::is_dominated(t(InitPopFuncValues[,1:nstats])))
    lengthfrontier<-length(frontier3)
    }
    xy.f <- matrix(InitPopFuncValues[frontier3, ],nrow=length(frontier3))
    colnames(xy.f)<-selectionstats
     if (plotiters){
       rbPal <- grDevices::colorRampPalette(c("light green", "yellow", "orange", "red")[4:1])
       xy.f2<-rbind(xy.f,InitPopFuncValues)
       PCH2<-c(rep(17,nrow(xy.f)),rep(1,nrow(InitPopFuncValues)))
       CEX2<-c(rep(1,nrow(xy.f)),rep(.3,nrow(InitPopFuncValues)))
       
       #This adds a column of color values
       # based on the y values
       
       colsval1<-disttoideal(xy.f)
       colsval2<-disttoideal(xy.f2)[-(1:nrow(xy.f))]
       colsval<-c(colsval1,(1+max(colsval1))*colsval2)
       PCH2[which.min(colsval)]<-13
       CEX2[which.min(colsval)]<-2
       Col1 <- rbPal(length(colsval1))[as.numeric(cut(colsval1,breaks = 20))]
       Col2 <- rep("#694002",nrow(xy.f2)-nrow(xy.f))
       Col<-c(Col1,Col2)
       Col<-scales::alpha(Col, 0.4)
       Col[which.min(colsval)]<-"#1e1e1d"
        if (ncol(xy.f)==2){
      plot(plotdirections[1]* xy.f2[,c(1)],plotdirections[2]* xy.f2[,c(2)], xlab=selectionstats[1],ylab=selectionstats[2], col=Col, pch=PCH2,cex=CEX2)
      } 
      if (ncol(xy.f2)==3){
        if (!is.null(axes.labels)){
          xlab=axes.labels[1]
          ylab=axes.labels[2]
          zlab=axes.labels[3]
          scatterplot3d::scatterplot3d(xy.f2*matrix(plotdirections, ncol=3,nrow=nrow(xy.f2), byrow=T),highlight.3d=FALSE, color=Col, pch=PCH2,cex.symbols=CEX2,grid=F, box=T, xlab=xlab,ylab=ylab,zlab=zlab)
          
        } else {
          scatterplot3d::scatterplot3d(xy.f2*matrix(plotdirections, ncol=3,nrow=nrow(xy.f2), byrow=T),highlight.3d=FALSE, color=Col, pch=PCH2,cex.symbols=CEX2,grid=F, box=T)
          
        }
        }
      if (ncol(xy.f2)>3){
        pairs(xy.f2*matrix(plotdirections, ncol=ncol(xy.f2),nrow=nrow(xy.f2), byrow=T), col=Col, pch=PCH2,cex=CEX2)
      }
      }
    
    
    #
    # Visualization.
    #
    # plot(xy.f[,c(1,2)], xlab="X", ylab="Y", pch=209, main="Quasiconvex Hull", col="red")
    # points(InitPopFuncValues[,c(1,2)], col="blue")
    
    ElitePop <- mclapply(frontier3, FUN = function(x){
      return(InitPop[[x]])
    }, mc.cores = mc.cores, mc.preschedule = T)
    ElitePopFuncValues <- InitPopFuncValues[frontier3,]
    
    for (iters in 1:niterations) {
      lengthfrontier<-0
      while(lengthfrontier<5){
      CurrentPop <- GenerateCrossesfromElites(Elites = ElitePop, 
                                              Candidates = Candidates, npopGA = npopGA, mutprob = mutprob)
      CurrentPop<-c(CurrentPop, ElitePop[1])
      CurrentPopFuncValues <- matrix(as.numeric(unlist(mclapply(CurrentPop,
                                                                FUN = function(x) {
                                                                  GAobjfunc(x)
                                                                }, mc.cores=mc.cores, mc.preschedule = T))), ncol=  nstats, byrow=TRUE)
      
      
      CurrentPop <- c(CurrentPop,ElitePop)
      notduplicated<-!duplicated(CurrentPop)
      CurrentPop<-CurrentPop[notduplicated]
      CurrentPopFuncValues<-rbind(CurrentPopFuncValues,ElitePopFuncValues)
      CurrentPopFuncValues<-CurrentPopFuncValues[notduplicated,]
      
      frontier3<- which(!emoa::is_dominated(t(CurrentPopFuncValues[,1:nstats])))
      lengthfrontier<-length(frontier3)
      }
      
      xy.f <- matrix(CurrentPopFuncValues[frontier3, ], nrow=length(frontier3))
      #
      # Visualization.
      #
      colnames(xy.f)<-selectionstats
      
      if (plotiters){
        rbPal <- grDevices::colorRampPalette(c("light green", "yellow", "orange", "red")[4:1])
        xy.f2<-rbind(xy.f,InitPopFuncValues)
        PCH2<-c(rep(17,nrow(xy.f)),rep(1,nrow(InitPopFuncValues)))
        CEX2<-c(rep(1,nrow(xy.f)),rep(.3,nrow(InitPopFuncValues)))
        
        #This adds a column of color values
        # based on the y values
        
        colsval1<-disttoideal(xy.f)
        colsval2<-disttoideal(xy.f2)[-(1:nrow(xy.f))]
        colsval<-c(colsval1,(1+max(colsval1))*colsval2)
        PCH2[which.min(colsval)]<-13
        CEX2[which.min(colsval)]<-2
        Col1 <- rbPal(length(colsval1))[as.numeric(cut(colsval1,breaks = 20))]
        Col2 <- rep("#694002",nrow(xy.f2)-nrow(xy.f))
        Col<-c(Col1,Col2)
        Col<-scales::alpha(Col, 0.4)
        Col[which.min(colsval)]<-"#1e1e1d"
        if (ncol(xy.f)==2){
          plot(xy.f2[,c(1)],xy.f2[,c(2)], xlab=selectionstats[1],ylab=selectionstats[2], col=Col, pch=PCH2,cex=CEX2)
        } 
        if (ncol(xy.f2)==3){
          if (!is.null(axes.labels)){
            xlab=axes.labels[1]
            ylab=axes.labels[2]
            zlab=axes.labels[3]
            scatterplot3d::scatterplot3d(xy.f2*matrix(plotdirections, ncol=3,nrow=nrow(xy.f2), byrow=T),highlight.3d=FALSE, color=Col, pch=PCH2,cex.symbols=CEX2,grid=F, box=T, xlab=xlab,ylab=ylab,zlab=zlab)
            
          } else {
            scatterplot3d::scatterplot3d(xy.f2*matrix(plotdirections, ncol=3,nrow=nrow(xy.f2), byrow=T),highlight.3d=FALSE, color=Col, pch=PCH2,cex.symbols=CEX2,grid=F, box=T)
            
          }
        }
        if (ncol(xy.f2)>3){
          pairs(xy.f2*matrix(plotdirections, ncol=ncol(xy.f2),nrow=nrow(xy.f2), byrow=T), col=Col, pch=PCH2,cex=CEX2)
        }
      }
      
      #points3d(CurrentPopFuncValues[,c(1,2,3)], col="blue")
      
      ElitePop <- mclapply(frontier3, FUN = function(x){
        return(CurrentPop[[x]])
      }, mc.cores = mc.cores, mc.preschedule = T)
      ElitePopFuncValues <- CurrentPopFuncValues[frontier3,]
    }
    
    return(list(ElitePop, ElitePopFuncValues))
  }
  
  outGA<-GAfunc(ntoselect=ntoselect, npopGA=npopGA, mutprob=mutprob, 
                niterations=nitGA) 
  return(outGA)
  
}
