# mirt package----
install.packages("mirt", dependencies = TRUE)
library(mirt)

# tutorial site----
# https://github.com/philchalmers/mirt/wiki

# the other packages----
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


# Expected Log Complete Likelihood function for item parameters----

# likelihood
L <- function(u, theta, a, b){
  # c とDは固定
  p <- ptheta(theta = theta, a = a, b = b, c = 0, D = 1)
  prod(p^u*(1-p)^(1-u))
}
# expected complete log likelihood
ELL <- function(r, N, a, b, theta){
  p <- ptheta(theta = theta, a = a, b = b, c = 0, D = 1)
  sum(r*log(p)+(N-r)*log(1-p))
}

# distribution for marginal
node <- seq(-4, 4, length.out = 31)
weight <- dnorm(node)/sum(dnorm(node)) # normalization

# data
set.seed(0204)
at <- rlnorm(30, 0.2, 0.5) # non negetive real
bt <- rnorm(30, 0, 1.5) # real
tt <- rnorm(5000, 0, 1) # real 
dat1 <- simdata(a = at, d = -bt*at, Theta = as.matrix(tt), itemtype = "dich")

# Estep
Lim <- matrix(nrow = 5000, ncol = 31) # subjects * nodes
Gim <- matrix(nrow = 5000, ncol = 31) # subjects * nodes
for(i in 1:5000){
  for(m in 1:31) Lim[i,m] <- L(dat1[i,], node[m], at, bt)
}

for(i in 1:5000){
  Gim[i,] <- Lim[i,]*weight / sum(Lim[i,]*weight)
}

Nm <- colSums(Gim) # 各ノードごとの期待受検者度数
rjm <- t(dat1) %*% Gim # 各ノードごとの期待正答受検者度数

# 等高線

# item 1
Estep_Likelihood <- matrix(nrow = 31, ncol = 31)
ai <-  seq(0, 3, length.out = 31)
bi <- seq(-1, 4, length.out = 31)
for(i in 1:31){
  for(j in 1:31){
    Estep_Likelihood[i, j] <- ELL(rjm[1,], Nm, ai[i], bi[j], node)
  }
}

colnames(Estep_Likelihood) <- bi
Estep_Likelihood <- cbind(a = ai, Estep_Likelihood)
Estep_Likelihood %<>% 
  as_tibble() %>% 
  gather(key = "b", value = "ELL", -a) %>% 
  map_df(as.numeric)
# ggplot
Estep_Likelihood %>% ggplot(aes(x = b, y = a, z = ELL))+
  geom_contour(bins = 100)+
  labs(title = paste("a =", round(at[1], digits = 5), ", b =", round(bt[1], digits = 5)))

# item 2
Estep_Likelihood <- matrix(nrow = 31, ncol = 31)
ai <-  seq(0, 0.5, length.out = 31)
bi <- seq(-1, 3, length.out = 31)
for(i in 1:31){
  for(j in 1:31){
    Estep_Likelihood[i, j] <- ELL(rjm[2,], Nm, ai[i], bi[j], node)
  }
}

colnames(Estep_Likelihood) <- bi
Estep_Likelihood <- cbind(a = ai, Estep_Likelihood)
Estep_Likelihood %<>% 
  as_tibble() %>% 
  gather(key = "b", value = "ELL", -a) %>% 
  map_df(as.numeric)
# ggplot
Estep_Likelihood %>% ggplot(aes(x = b, y = a, z = ELL))+
  geom_contour(bins = 100)+
  labs(title = paste("a =", round(at[2], digits = 5), ", b =", round(bt[2], digits = 5)))

# example plot for multigroup distribution----
tibble(theta = c(-4:4)) %>% ggplot(aes(x = theta))+
  stat_function(fun = dnorm, args = list(mean = -1, sd = 1), colour = 2)+
  stat_function(fun = dnorm, args = list(mean = 1, sd = 1), colour = 3)+
  stat_function(fun = function(x, arg1, arg2) 0.5*dnorm(x, arg1[1], arg1[2])+0.5*dnorm(x, arg2[1], arg2[2]) , args = list(arg1 = c(-1,1), arg2 = c(1, 1)))

# Log likelihood function for item parameters----


# Single group IRT analysis via mirt package----

# simulation data
set.seed(0204)
at <- rlnorm(30, 0.2, 0.5) # non negetive real
bt <- rnorm(30, 0, 1.5) # real
tt <- rnorm(5000, 0, 1) # real 
dat1 <- simdata(a = at, d = -bt*at, Theta = as.matrix(tt), itemtype = "dich") # b = -d/a then d = -a*b

# model object
mod1 <- mirt.model('Factor1 = 1-30') # One factor

# estimation
fit1 <- mirt(data = dat1, model = mod1, itemtype = "2PL", SE = TRUE, SE.type = 'complete') # mirt(data = dat, model = 1, itemtype = "2PL", ) でもOK

# extract parameter and fit index
par1 <- coef(fit1, IRTpars = T, simplify = T, printSE = T)
ml1 <- fscores(fit1, method = "ML")

# trace plot
plot(fit1, type = "trace")
itemplot(fit1, 1)

# fit
# (total) test fit index
M2(fit1)

# item fit index
# infit statistics can be estimated only if model is Rasch.
itemfit(fit1, fit_stats = "S_X2")
itemfit(fit1, fit_stats = "G2")
# empirical item trace line plot
itemfit(fit1, empirical.plot = c(3))

# person fit index
personfit(fit1, method = "MAP") # this values follow standard normal distribution (empirically)

# expected total score (test characteristic curve)
plot(fit1)

# test information function
testinfo(fit1, sort(ml1)) %>% plot(x = sort(ml1), type = "l", xlim = c(-6,6))
plot(fit1, type = "info")
##
#  dentype = 'Davidian-4'
# plot(fit1, type = "Davidian") # only if estimated via Davidian curve IRT


# Multiple group IRT analysis via mirt packagte----
set.seed(19930828)
at1 <- rlnorm(45, 0.2, 0.5) # non negetive real
bt1 <- rnorm(45, 0, 1) # real
# group 1
tt1 <- rnorm(5000, 0, 1) # real 
tmp1 <- simdata(a = at1, d = -bt1*at1, Theta = as.matrix(tt1), itemtype = "dich") %>% as.data.frame()
tmp1[,31:45] <- NA
# at1 <- at1[1:30]
# bt1 <- bt1[1:30]

set.seed(19930829)
at2 <- c( at1[1:30], rlnorm(15, 0.2, 0.5) )# non negetive real
bt2 <- c( bt1[1:30], rnorm(15, 1, 1) ) # real
# group 2
tt2 <- rnorm(5000, 1, 1) # real 
tmp2 <- simdata(a = at2, d = -bt2*at2, Theta = as.matrix(tt2), itemtype = "dich")  %>% as.data.frame()
tmp2[,1:15] <- NA

# join
dat2 <- rbind(tmp1, tmp2)

# group index vector (must be character vector)
grp2 <- c(rep("G1",5000), rep("G2",5000))

# model
mod2 <- mirt.model('Factor1 = 1-45') # One factor

# constraint
const <- c("free_var", "free_means", colnames(dat2))

# estimate
fit2 <- multipleGroup(dat2, mod2, group = grp2, invariance = const, dentype = "empiricalhist",
                      itemtype = "2PL", empiricalhist = T, accelerate = "squarem")

# 
coef(fit2, simplify = T, IRTpars = T)
plot(fit2, type = 'empiricalhist', npts = 60)

# grp2 <- c(rep(1,5000), rep(2,5000))
# irtfun2::estip2(cbind(grp2, dat2), IDc = 0, Gc = 1, fc = 2, D = 1.0)
