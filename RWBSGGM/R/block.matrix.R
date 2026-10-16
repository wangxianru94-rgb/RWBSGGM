block.matrix=function(dd){
d=dd
k=floor(d/5)
NN=seq(1,d,5)
A=matrix(0,nrow=d,ncol=d)
for(i in 1:k){
  for(j in NN[i]:(NN[i]+4)){
    for(l in j:(NN[i]+4)){
      A[j,l]=0.6
    }
  }
}
A=A+t(A)
diag(A)=runif(d,1,2)
A
}
















