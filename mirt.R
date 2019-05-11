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


# Multi Group estimation in mirt function
# data
dat6 <- sim_dat_st[,c(-1, -2)]
mod6 <- mirt.model('F = 1-70')
grp6 <- as.vector(apply(matrix(sim_dat_st$grade, ncol = 1), 1, paste0, "D"))
# grp6 <- c(rep("G1", 500), rep("G2", 500), rep("G3", 500), rep("G4", 500), rep("G5", 500))
fit6 <- multipleGroup(dat6, mod6, grp6, invariance=c('free_var','free_means', colnames(dat6)), 
                      rotate = "none", itemtype = "2PL", quadpts = 31)
param <- coef(fit6, IRTpars = T, simplify = T)

