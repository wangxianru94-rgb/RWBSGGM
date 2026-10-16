#==================================================================
block.precision.matrix=function(dd,alpha){
alpha=alpha
d=dd
k=floor(d/5)
NN=seq(1,d,5)
NN2=c(1:d)
if(alpha>0){ 
  location.change=(floor(d*(1-alpha)+1):d)
  NN2[location.change]=sample(location.change,length(location.change),replace=F)
}
A=matrix(0,nrow=d,ncol=d)
for(i in 1:k){
  for(j in NN[i]:(NN[i]+4)){
    for(l in j:(NN[i]+4)){
      A[NN2[j],NN2[l]]=0.8
    }
  }
}
A=A+t(A)
diag(A)=1
A
}

#======================================================
block.precision.matrix2=function(dd,alpha){
  alpha=alpha
  d=dd
  k=floor(d/5)
  NN=seq(1,d,5)
  NN2=c(1:d)
  if(alpha>0){ 
    location.change=(d-alpha):d
    NN2[location.change]=sample(location.change,length(location.change),replace=F)
  }
  A=matrix(0,nrow=d,ncol=d)
  for(i in 1:k){
    for(j in NN[i]:(NN[i]+4)){
      for(l in j:(NN[i]+4)){
        A[NN2[j],NN2[l]]=0.8
      }
    }
  }
  A=A+t(A)
  diag(A)=1
  A
}

#=======================================================
chain.graph2=function(p,alpha){
  # alpha denotes the difference between the two graphs
  TT=runif(p-1,0.5,1) 
  PI=c(1:p)
  location.change=(p-alpha):p
  PI[location.change]=sample(location.change,length(location.change),replace=F)
  SIGMA1=matrix(0,nrow=p,ncol=p)
  SIGMA2=matrix(0,nrow=p,ncol=p)
  
  for (i in 1:(p-1)){
    for(j in (i+1):p){
      a=max(i,j) 
      b=min(i,j) 
      SIGMA1[i,j]=exp(-abs(sum(TT[(b):(a-1)]))/5)
    }
  }
  
  for (i in 1:(p-1)){
    for(j in (i+1):p){
      a=max(PI[i],PI[j]) 
      b=min(PI[i],PI[j]) 
      SIGMA2[i,j]=exp(-abs(sum(TT[(b):(a-1)]))/5)
    }
  }
  
  
  SIGMA1=SIGMA1+t(SIGMA1)
  diag(SIGMA1)=rep(1,p)
  OMEGA1=solve(SIGMA1)
  OMEGA.VECTOR1=as.vector(OMEGA1)
  OMEGA.VECTOR1[abs(OMEGA.VECTOR1)<=0.001]=0
  PRECISION1=matrix(OMEGA.VECTOR1,nrow=p)
  PRECISION11=diag(diag(PRECISION1)^{-1/2})%*%PRECISION1%*%diag(diag(PRECISION1)^{-1/2})
  
  SIGMA2=SIGMA2+t(SIGMA2)
  diag(SIGMA2)=rep(1,p)
  OMEGA2=solve(SIGMA2)
  OMEGA.VECTOR2=as.vector(OMEGA2)
  OMEGA.VECTOR2[abs(OMEGA.VECTOR2)<=0.001]=0
  PRECISION2=matrix(OMEGA.VECTOR2,nrow=p)
  PRECISION22=diag(diag(PRECISION2)^{-1/2})%*%PRECISION2%*%diag(diag(PRECISION2)^{-1/2})
  
  res=list(graph1=PRECISION11,graph2=PRECISION22)
  res
} 


#==================================================================
chain.graph=function(p,alpha){
  # alpha denotes the difference between the two graphs
  TT=runif(p-1,0.5,1) 
  PI=c(1:p)
  location.change=(floor(p*(1-alpha)+1):p)
  PI[location.change]=sample(location.change,length(location.change),replace=F)
  SIGMA1=matrix(0,nrow=p,ncol=p)
  SIGMA2=matrix(0,nrow=p,ncol=p)
  
  for (i in 1:(p-1)){
    for(j in (i+1):p){
      a=max(i,j) 
      b=min(i,j) 
      SIGMA1[i,j]=exp(-abs(sum(TT[(b):(a-1)]))/5)
    }
  }
  
  for (i in 1:(p-1)){
    for(j in (i+1):p){
      a=max(PI[i],PI[j]) 
      b=min(PI[i],PI[j]) 
      SIGMA2[i,j]=exp(-abs(sum(TT[(b):(a-1)]))/5)
    }
  }
  
  
  SIGMA1=SIGMA1+t(SIGMA1)
  diag(SIGMA1)=rep(1,p)
  OMEGA1=solve(SIGMA1)
  OMEGA.VECTOR1=as.vector(OMEGA1)
  OMEGA.VECTOR1[abs(OMEGA.VECTOR1)<=0.001]=0
  PRECISION1=matrix(OMEGA.VECTOR1,nrow=p)
  PRECISION11=diag(diag(PRECISION1)^{-1/2})%*%PRECISION1%*%diag(diag(PRECISION1)^{-1/2})
  
  SIGMA2=SIGMA2+t(SIGMA2)
  diag(SIGMA2)=rep(1,p)
  OMEGA2=solve(SIGMA2)
  OMEGA.VECTOR2=as.vector(OMEGA2)
  OMEGA.VECTOR2[abs(OMEGA.VECTOR2)<=0.001]=0
  PRECISION2=matrix(OMEGA.VECTOR2,nrow=p)
  PRECISION22=diag(diag(PRECISION2)^{-1/2})%*%PRECISION2%*%diag(diag(PRECISION2)^{-1/2})
  
  res=list(graph1=PRECISION11,graph2=PRECISION22)
  res
} 

#==================================================================================

recovery=function(thr1,thr2,precision1,precision2,precision1.head,precision2.head){
  p=ncol(precision1)
  best.thr1=thr1
  adjacent1.true=matrix(as.numeric(abs(precision1)>0),nrow=p)       
  edge1.true=(sum(adjacent1.true)-p)/2                             
  temp.adjacent1=matrix(as.numeric(abs(precision1.head)>best.thr1),nrow=p)
  temp.adjacent1=temp.adjacent1+t(temp.adjacent1)
  temp.adjacent1=as.vector(temp.adjacent1)
  temp.adjacent1[which(temp.adjacent1==2)]=1
  adjacent1.head=matrix(temp.adjacent1,nrow=p)
  diag(adjacent1.head)=rep(1,p)
  edge1.head=(sum(adjacent1.head)-p)/2
  true1.estimate=(length(which( adjacent1.true==1& adjacent1.head==1 ))-p)/2
  precision1.score=true1.estimate/edge1.head
  recall1.score=true1.estimate/edge1.true
  F.score1=2*precision1.score*recall1.score/(precision1.score+recall1.score)
  
  #========recovery 2================================================
  best.thr2=thr2
  adjacent2.true=matrix(as.numeric(abs(precision2)>0),nrow=p)  
  edge2.true=(sum(adjacent2.true)-p)/2     
  
  temp.adjacent2=matrix(as.numeric(abs(precision2.head)>best.thr2),nrow=p)
  temp.adjacent2=temp.adjacent2+t(temp.adjacent2)
  temp.adjacent2=as.vector(temp.adjacent2)
  temp.adjacent2[which(temp.adjacent2==2)]=1
  adjacent2.head=matrix(temp.adjacent2,nrow=p)
  diag(adjacent2.head)=rep(1,p)
  edge2.head=(sum(adjacent2.head)-p)/2
  true2.estimate=(length(which( adjacent2.true==1& adjacent2.head==1 ))-p)/2
  precision2.score=true2.estimate/edge2.head
  recall2.score=true2.estimate/edge2.true
  F.score2=2*precision2.score*recall2.score/(precision2.score+recall2.score)
  res=c(precision1.score,recall1.score,F.score1,precision2.score,recall2.score,F.score2)
  res
}

#==============================================================================
cusum=function(x){
  n=length(x)
  cus=numeric()
  for( i in 1:(n-1)){
    cus[i]=sqrt(n)*sqrt(i/n)*sqrt(1-i/n)*  (mean(x[1:i])-mean(x[(i+1):n]))
  }
  cus
}





