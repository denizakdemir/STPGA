GenerateCrossesfromElites<-
function(Elites, Candidates, npop, mutprob,mc.cores=1, mutintensity=1, memoryfortabu=NULL){
  checkmemoryfortabu<-function(out){
    n=length(out)
    outTF<-F
    for (i in 1:length(memoryfortabu)){
    outTF=Reduce("||", lapply(memoryfortabu[[i]], function(x){return(n==length(union(x,out)))}))
    if (outTF){
     output=outTF
     break
     } else {output=F}
    }
    return(output)
  }
  
  checkmemoryfortabu2<-function(out){
    n=length(out)
    outTF<-F
   for (i in 1:length(memoryfortabu)){
    for (j in 1:length(memoryfortabu[[i]])){
      if (union(memoryfortabu[[i]][[j]],out)==n){outTF=TRUE}
      if (outTF){break}
    }
     if (outTF){break}
   }
     return(outTF)
  }
newcrosses<-mclapply(1:npop, FUN=function(x){
if (is.null(memoryfortabu)){
  x1<-Elites[[sample(1:length(Elites),1)]]
  x2<-Elites[[sample(1:length(Elites),1)]]
  out<-makeonecross(x1=x1,x2=x2,Candidates=Candidates,mutprob=mutprob, mutintensity=mutintensity)
  return(out)
} else {
  inmemorytabu=T
  x1<-Elites[[sample(1:length(Elites),1)]]
  mutintensitytabu=mutintensity
  mutprobtabu=mutprob
  while (inmemorytabu){
    x2<-Elites[[sample(1:length(Elites),1)]]
  out<-makeonecross(x1=x1,x2=x2,Candidates=Candidates,mutprob=mutprobtabu, mutintensity=mutintensitytabu)  
  mutintensitytabu=min(mutintensitytabu*1.1,floor(length(out)/5))
  mutprobtabu=min(mutprobtabu*1.1,1)
  if (!checkmemoryfortabu(out)){inmemorytabu=F}
 
  }
  return(out)
  }
}, mc.cores=mc.cores,  mc.set.seed = T)
	return(newcrosses)
}
