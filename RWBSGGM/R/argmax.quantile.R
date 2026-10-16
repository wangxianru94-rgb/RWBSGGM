#library(pracma)  # 用于数值积分和其他数学运算

# 定义生成 [0, ∞) 上布朗运动的函数
generate_brownian_motion_positive <- function(t_max, n) {
  dt <- t_max / n  # 计算时间步长
  W <- c(0, cumsum(rnorm(n, mean = 0, sd = sqrt(dt))))  # 生成 [0, ∞) 上的布朗运动过程
  t_values <- seq(0, t_max, length.out = n + 1)  # 生成时间序列
  return(data.frame(t = t_values, W = W))
}



brownian_motion.full.generate=function(s1,s2,n){
  w1=generate_brownian_motion_positive (s2,n)
  w2=generate_brownian_motion_positive (-s1,n)
  w2$t=-1*rev(w2$t)
  w2$W=-1*rev(w2$W)
  W_full <- rbind(w2, w1[-1, ])
  return(W_full)
}


Z_process <- function(sigma1,sigma1_star,sigma2,sigma2_star,s_values=seq(-20,20,length.out =2*1000+1)) {
  # 计算常数
  
  
  # 生成布朗运动
  t_range <- range(s_values)
  W1 <- brownian_motion.full.generate(t_range[1], t_range[2], (length(s_values)-1)/2)
  Z <- numeric(length(s_values))  # 初始化Z(s)的值
  for (i in seq_along(s_values)) {  # 对每个s值进行迭代
    s <- s_values[i]
    if (s > 0) {
      Z[i] <-  2*sigma1* W1$W[i]-sigma1_star^2*abs(s)  # 计算s > 0时的Z(s)
    } else if (s < 0) {
      Z[i] <- 2*sigma2* W1$W[i]-sigma2_star^2*abs(s)  # 计算s < 0时的Z(s)
    } else {
      Z[i] <- 0  # s = 0时Z(s) = 0
    }
  }
  return(Z)
}

argmax.quantile=function(sigma1,sigma1_star,sigma2,sigma2_star,
                          s_values=seq(-20,20,length.out=2*1000+1),
                          nsim=2000){
  argmax.sim=numeric(nsim)
  for(i in seq_len(nsim)){
    Z=Z_process(sigma1,sigma1_star,sigma2,sigma2_star,s_values=s_values)
    argmax.sim[i]=s_values[which.max(Z)]
  }
  quantile.argmax=quantile(argmax.sim,c(0.05,0.95))
  return(quantile.argmax)
}









