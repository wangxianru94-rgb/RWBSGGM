#' generate band type covarianc matrix
#'
#' @author LIUBIN
#' @param dd dimension
#' @return d*d matrix
#' @export

band.matrix=function(dd){
  d=dd
  A=matrix(0,nrow=d,ncol=d)
  for(i in 1:d){
    for(j in 1:i){
      A[i,j]=0.8^abs(i-j)
    }
  }
  A=A+t(A)-diag(rep(1,d))
  A
}  