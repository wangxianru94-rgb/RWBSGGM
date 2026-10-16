standardize<- function(X){
  n<- nrow(X)
  d<- ncol(X)
  X.mean<- matrix(rep(apply(X,2,mean),n),nrow=n,byrow=T)
  X.sd<- apply(X,2,sd)
  X.stand<- (X-X.mean)/X.sd
  X.stand
  
}