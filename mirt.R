# mirt package
library(mirt)
library(irtfun2)

# mirt function run test
dat1 <- sim_data_1[,-1]
mod1 <- mirt.model('G = 1-30')
fit1 <- mirt(dat1, mod1, itemtype = "2PL")
coef(fit1, simplyfy = T)

dat2 <- dat_2[, -1]
mod2 <- mirt.model('g = 1-5')
fit2 <- mirt(dat2, mod2, itemtype = "graded")
fscores(fit2)
personfit(fit2)


# plausible values check
# compare irtfun2 and mirt
res3 <- estip2(sim_data_2, fc = 2, D = 1.0)
pv3 <- estheta(sim_data_2, param = res3$para, est = "PVs", D = 1.0, fc = 2, method = "NR", sampling_engine = "rejection_Cpp", gc = 0)
test <- estheta(sim_data_2, param = res3$para, est = "MAP", method = "NR", D = 1.0, fc = 2, gc = 0)
test2 <- estheta(sim_data_2, param = res3$para, est = "MAP", method = "Brent", D = 1.0, fc = 2, gc = 0)
test$`MAPmean&sd`
test2$`MAPmean&sd`

stat_pv(pv3$PVs_only)

dat4 <- sim_data_2[, -1]
mod1
fit4 <- mirt(dat4, mod1, itemtype = "2PL")
coef(fit4, IRTpars = T, simplify = T)
pv4 <- fscores(fit4, plausible.draws = 10)
lapply(pv4, mean)
lapply(pv4, sd)

# validity check
# for easy estimation
sim_data_2$ID <- formatC(c(1:3000), width = 5, flag = 0)
write.table(sim_data_2, file = "sim_dat_2.dat", row.names = F, col.names = F, sep = "")
pv5 <- read.csv(file = "C:/Users/bc0089985/Downloads/EasyEstimation/sim_dat_2ThetaPV.csv")
lapply(pv5[,-1], mean)
lapply(pv5[,-1], sd)


# Multi Group estimation in mirt function
# data
dat6 <- sim_dat_st[,c(-1, -2)]
mod6 <- mirt.model('F = 1-70')
grp6 <- as.vector(apply(matrix(sim_dat_st$grade, ncol = 1), 1, paste0, "D"))
# grp6 <- c(rep("G1", 500), rep("G2", 500), rep("G3", 500), rep("G4", 500), rep("G5", 500))
fit6 <- multipleGroup(dat6, mod6, grp6, invariance=c('free_var','free_means', colnames(dat6)), 
                      rotate = "none", itemtype = "2PL", quadpts = 31)
param <- coef(fit6, IRTpars = T)

estip2(sim_dat_st, Gc = 2, D = 1)
