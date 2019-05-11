install.packages("mirt", dependencies = TRUE)
library(mirt)

# the other packages
library(latex2exp) # for latex coding
library(tidyverse) # for data transformation and visualization

# Log likelihood function for theta graph----
# P(theta) in two-parameter logisticmodel
ptheta <- function(theta,a,b,c,D=1.702){
  c+(1-c)/(1+exp(-D*a*(theta-b)))
}

iif <- function(theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  I <- D^2*a^2*(1-p)*(p-c)^2/((1-c)^2*p)
}

# 対数尤度
LL <- function(u,theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  sum(u*log(p)+(1-u)*log(1-p),na.rm = T)
}

# 一階偏微分
fpd <- function(xi,theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  D*sum(a*(xi - p)*(p-c)/(p*(1-c)), na.rm = T)
}

# 二階偏微分
spd <- function(xi,theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  D^2*sum(a^2*(p-c)*(xi*c-p^2)*(1-p)/(p^2*(1-c)^2),na.rm = T)
}

# テスト情報量（尤度関数の二階偏微分の負の期待値）
pitheta <- function(theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  I <- D^2*sum(a^2*(1-p)*(p-c)^2/((1-c)^2*p), na.rm = T)
}

line_ic <- function(y,a,x,X){
  #傾きとx,yの値から切片を求めて直線の値yを再計算する関数
  b <- y-(a*x)
  a*X+b
}

line_ic2 <- function(y,a,x,X){
  #傾きとx,yの値から切片を求めて直線の値yを再計算する関数
  a <- -a
  b <- y-(a*x)
  a*X+b
} 

nr_fun <- function(t0,fp,fpp){
  # ニュートンラフソンの更新した値を計算する関数
  t0-fp/fpp
}

fs_fun <- function(t0,fp,I){
  # フィッシャースコアリングの更新した値を計算する関数
  t0+fp/I
}


# response pattern"
set.seed(123)
u <- sample(c(0,1),30, replace = T)
# item parameters
# 2PLM
a <- rlnorm(30,meanlog = -0.5,sdlog = 0.5)
b <- rnorm(30)
c <- rep(0,30)

dat_L <- data.frame(x=seq(-4,4,length.out = 101))
dat_L$LL <- apply(matrix(dat_L$x), 1, LL, u=u,a=a,b=b,c=c,D=1.702) 
dat_L$fpd <- apply(matrix(dat_L$x), 1, fpd, xi=u,a=a,b=b,c=c,D=1.702) 
dat_L$spd <- apply(matrix(dat_L$x), 1, spd, xi=u,a=a,b=b,c=c,D=1.702) 
dat_L$I <- -apply(matrix(dat_L$x), 1, pitheta, a=a,b=b,c=c,D=1.702) 

LL_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=LL))+
  labs(x=TeX("$\\theta$"), y=TeX("$lnL(\\theta)$"))

fpd_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=fpd))+
  labs(x=TeX("$\\theta$"), y=TeX("$lnL'(\\theta)$"))

spd_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=spd))+
  labs(x=TeX("$\\theta$"), y=TeX("$lnL''(\\theta)$"))

info_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=I))+
  labs(x=TeX("$\\theta$"), y=TeX("$I(\\theta)$"))

LL_g

fpd_g

spd_g

info_g


# まずは-2から
NR_g <- fpd_g+
  stat_function(fun=line_ic, args = list(y=fpd(u,-2,a,b,c,D=1.702), a=spd(u,-2,a,b,c,D=1.702), x=-2), aes(colour=1))
t2 <- nr_fun(-2,fpd(u,-2,a,b,c,D=1.702),spd(u,-2,a,b,c,D=1.702))

# 自動作図
NR_g <- fpd_g
t0 <- -2
t <- 0
conv <- T
while(conv){
  t <- t+1
  slope <- spd(u,t0,a,b,c,D=1.702)
  NR_g <- NR_g+
    stat_function(fun=line_ic, args = list(y=fpd(u,t0,a,b,c,D=1.702), a=slope, x=t0),
                  aes_q(colour=sprintf("iterations=%d, theta=%.3f",t,t0)))
  t1 <- nr_fun(t0,fpd(u,t0,a,b,c,D=1.702),spd(u,t0,a,b,c,D=1.702))
  if(abs(t1-t0)<0.001) conv <- FALSE
  t0 <- t1
}
NR_g <- NR_g+labs(x=TeX("$\\theta$"), y=TeX("$lnL'(\\theta)$"), colour=TeX("iterations"))

FS_g <- fpd_g
t0 <- -2
t <- 0
conv <- T
while(conv){
  t <- t+1
  I <- pitheta(t0,a,b,c,D=1.702)
  FS_g <- FS_g+
    stat_function(fun=line_ic2, args = list(y=fpd(u,t0,a,b,c,D=1.702), a=I, x=t0),
                  aes_q(colour=sprintf("iterations=%d, theta=%.3f",t,t0)))
  t1 <- fs_fun(t0,fpd(u,t0,a,b,c,D=1.702),pitheta(t0,a,b,c,D=1.702))
  if(abs(t1-t0)<0.001) conv <- FALSE
  t0 <- t1
}
FS_g <- FS_g+labs(x=TeX("$\\theta$"), y=TeX("$lnL'(\\theta)$"), colour=TeX("iterations"))

NR_g

FS_g


# 3PLM
set.seed(123)
a <- rlnorm(30,meanlog = -0.5,sdlog = 0.5)
b <- rnorm(30)
c <- rbeta(30,2,10)

dat_L <- data.frame(x=seq(-4,4,length.out = 101))
dat_L$LL <- apply(matrix(dat_L$x), 1, LL, u=u,a=a,b=b,c=c,D=1.702) 
dat_L$fpd <- apply(matrix(dat_L$x), 1, fpd, xi=u,a=a,b=b,c=c,D=1.702) 
dat_L$spd <- apply(matrix(dat_L$x), 1, spd, xi=u,a=a,b=b,c=c,D=1.702) 
dat_L$I <- -apply(matrix(dat_L$x), 1, pitheta, a=a,b=b,c=c,D=1.702) 

LL_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=LL))+
  labs(x=TeX("$\\theta$"), y=TeX("$lnL(\\theta)$"))

fpd_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=fpd))+
  labs(x=TeX("$\\theta$"), y=TeX("$lnL'(\\theta)$"))

spd_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=spd))+
  labs(x=TeX("$\\theta$"), y=TeX("$lnL''(\\theta)$"))

info_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=I))+
  labs(x=TeX("$\\theta$"), y=TeX("$I(\\theta)$"))

LL_g

fpd_g

spd_g

info_g

# 自動作図
NR_g <- fpd_g
t0 <- 0
t <- 0
conv <- T
while(conv){
  t <- t+1
  slope <- spd(u,t0,a,b,c,D=1.702)
  NR_g <- NR_g+
    stat_function(fun=line_ic, args = list(y=fpd(u,t0,a,b,c,D=1.702), a=slope, x=t0),
                  aes_q(colour=sprintf("iterations=%d, theta=%.3f",t,t0)))
  t1 <- nr_fun(t0,fpd(u,t0,a,b,c,D=1.702),spd(u,t0,a,b,c,D=1.702))
  if(abs(t1-t0)<0.001) conv <- FALSE
  t0 <- t1
}
NR_g <- NR_g+labs(x=TeX("$\\theta$"), y=TeX("$lnL'(\\theta)$"), colour=TeX("iterations"))

FS_g <- fpd_g
t0 <- 0
t <- 0
conv <- T
while(conv){
  t <- t+1
  I <- pitheta(t0,a,b,c,D=1.702)
  FS_g <- FS_g+
    stat_function(fun=line_ic2, args = list(y=fpd(u,t0,a,b,c,D=1.702), a=I, x=t0),
                  aes_q(colour=sprintf("iterations=%d, theta=%.3f",t,t0)))
  t1 <- fs_fun(t0,fpd(u,t0,a,b,c,D=1.702),pitheta(t0,a,b,c,D=1.702))
  if(abs(t1-t0)<0.001) conv <- FALSE
  t0 <- t1
}
FS_g <- FS_g+labs(x=TeX("$\\theta$"), y=TeX("$lnL'(\\theta)$"), colour=TeX("iterations"))

NR_g
FS_g


# MAP

# 尤度関数
LL_b <- function(u,theta,a,b,c,mu,sigma,D){
  p <- ptheta(theta,a,b,c,D)
  sum(log(p)*u+log(1-p)*(1-u)-0.5*((theta-mu)/sigma)^2,na.rm = T)
}

# 一階偏微分
fpdLPD <- function(xi,theta,a,b,c,mu,sigma,D){
  p <- ptheta(theta,a,b,c,D)
  D*sum(a*(xi - p)*(p-c)/(p*(1-c)), na.rm = T)-1/sigma^2*(theta-mu)
}

# 二階偏微分
spdLPD <- function(xi,theta,a,b,c,sigma,D){
  p <- ptheta(theta,a,b,c,D)
  D^2*sum(a^2*(p-c)*(xi*c-p^2)*(1-p)/(p^2*(1-c)^2),na.rm = T) - 1/sigma
}

set.seed(123)
a <- rlnorm(30,meanlog = -0.5, sdlog = 0.5)
b <- rnorm(30)
c <- rbeta(30,2,10)

dat_L <- data.frame(x=seq(-4,4,length.out = 101))
dat_L$LL <- apply(matrix(dat_L$x), 1, LL, u=u,a=a,b=b,c=c,D=1.702) 
dat_L$fpd <- apply(matrix(dat_L$x), 1, fpd, xi=u,a=a,b=b,c=c,D=1.702) 
dat_L$spd <- apply(matrix(dat_L$x), 1, spd, xi=u,a=a,b=b,c=c,D=1.702) 
dat_L$LL_B <- apply(matrix(dat_L$x), 1, LL_b, u=u,a=a,b=b,c=c,D=1.702,mu=0,sigma=1) 
dat_L$fpdLPD <- apply(matrix(dat_L$x), 1, fpdLPD, xi=u,a=a,b=b,c=c,D=1.702,mu=0,sigma=1) 
dat_L$spdLPD <- apply(matrix(dat_L$x), 1, spdLPD, xi=u,a=a,b=b,c=c,D=1.702,sigma=1) 

MAP_LL_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=LL,colour="MLE"))+
  geom_line(aes(y=LL_B, colour="MAP"))+
  labs(x=TeX("$\\theta$"),y="",colour="mathod")

MAP_fpd_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=fpd,colour="MLE"))+
  geom_line(aes(y=fpdLPD, colour="MAP"))+
  labs(x=TeX("$\\theta$"),y="",colour="mathod")

MAP_spd_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=spd,colour="MLE"))+
  geom_line(aes(y=spdLPD, colour="MAP"))+
  labs(x=TeX("$\\theta$"),y="",colour="mathod")

# 事前分布のパラメタを変化させた場合の尤度方程式の変化
dat_L <- data.frame(x=seq(-4,4,length.out = 101))
dat_L$fpd <- apply(matrix(dat_L$x), 1, fpd, xi=u,a=a,b=b,c=c,D=1.702) 
dat_L$fpdLPD1 <- apply(matrix(dat_L$x), 1, fpdLPD, xi=u,a=a,b=b,c=c,D=1.702,mu=0,sigma=1) 
dat_L$fpdLPD2 <- apply(matrix(dat_L$x), 1, fpdLPD, xi=u,a=a,b=b,c=c,D=1.702,mu=5,sigma=1) 
dat_L$fpdLPD3 <- apply(matrix(dat_L$x), 1, fpdLPD, xi=u,a=a,b=b,c=c,D=1.702,mu=0,sigma=0.8) 

MAP_fpd2_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=fpd,colour="MLE"))+
  geom_line(aes(y=fpdLPD1, colour="N(0,1)"))+
  geom_line(aes(y=fpdLPD2, colour="N(5,1)"))+
  geom_line(aes(y=fpdLPD3, colour="N(0,0.64)"))+
  labs(x=TeX("$\\theta$"),y=TeX("$ln\\L'(\\theta)+lnp'(\\theta)$"),colour="condition")

MAP_fpd2_g

# test data


# Log likelihood function for item parameters----

# IRT analysis via mirt package

# test data
