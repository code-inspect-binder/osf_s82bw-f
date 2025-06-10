
##########################################
## SIMULATION:
##
## Test the estimation of continuous
## predictor variables in the latent-trait
## MPT approach
##########################################

# Install and load TreeBUGS
install.packages("TreeBUGS")
library(TreeBUGS)
# setwd("C:/mpt/")



##########################################
## Data-generating values

# participants:
N <- 50
# standardized regression coefficients
b <- c(d1 = .5, D1 = -.3)
# covariate scaling (should not matter):
covMean = 100
covSD = 15
# latent probit-means:
mu <-  c(a = .3, b = -.1, d1 = .6,  D1 = .3)
# latent SD:
sig <- c(a = .6, b = .5, d1 = 1, D1 = .2)
# number of 2HTSM items:
numItems <- c(E=16, U=16, N=32)



##########################################
#### data-generating function

replication <- function(N, b, mu, sig, numItems, covMean=0, covSD=1){

  library(MASS)
  theta <- mvrnorm(N, mu = mu, Sigma = diag(4)*sig^2)
  covData <- data.frame(cov = rnorm(N, covMean, covSD))

  # latent-probit regression
  theta[,"d1"] <- b["d1"]*(covData$cov-covMean)/covSD  + theta[,"d1"]
  theta[,"D1"] <- b["D1"]*(covData$cov-covMean)/covSD  + theta[,"D1"]

  # generate MPT data:
  gendat <- genMPT(theta = pnorm(theta),
                   restrictions = "model/restrictions.txt",
                   numItems = numItems,
                   eqnfile="model/2htsm.eqn")

  # fit model
  sim.pred <- traitMPT(eqnfile="model/2htsm.eqn",
                       data = gendat,
                       restrictions = "model/restrictions.txt",
                       covData = covData,
                       predStructure = list("D1 d1 ; cov"),
                       n.chain = 8, n.iter = 20000, n.adapt = 5000,
                       n.burnin = 5000, n.thin = 10,silent.jags=TRUE )
  ###### check convergence manually in several replications!
  # plot(sim.pred)
  # plot(sim.pred, "slope")
  # plot(sim.pred, "sigma")
  # plot(sim.pred, "rho")
  # summary(sim.pred)
  sim.pred$runjags <- NULL  # do not save MCMC samples
  sim.pred
}



##########################################
#### Run simulation

# number of replications:
M <- 500
print(t0 <- Sys.time())
sim <- list()
for(m in 1:M){
  try(sim[[m]] <- replication(N, b, mu, sig, numItems = numItems,
                              covMean = covMean, covSD = covSD))
  save(sim, file="sim_tmp.RData")
  cat("\n", m, "/", M, "\n")
}
t1 <- Sys.time()
print(t1-t0)


##########################################
#### extract relevant results

# select convergent replications
rep <- sapply(sim, function(xx){
  all(xx$mcmc.summ[,"Rhat"] < 1.05, na.rm = TRUE)
})
sel <- c("Mean","2.5%","97.5%")
res <- do.call("rbind", sapply(sim[rep], function(xx){
  x <- xx$summary$groupParameters
  list(x$slope[,sel]*sd(xx$mptInfo$covData$cov), # get standardized slope!
       x$mean[,sel],
       x$sigma[,sel])
}))
res <- data.frame(res, Parameter = rownames(res))

res$True <- c(rev(b), pnorm(mu), sig)
res$CI <- with(res, True >= X2.5. & True <= X97.5.)
res$overZero <- with(res, (0 <= X2.5. & 0 <= X97.5.) | 0 >= X2.5. & 0 >= X97.5.)
res$AbsDiff <- with(res, abs(Mean - True))
head(res)


##########################################
#### Summarize results

library(dplyr)
tab <- res %>%
  group_by(Parameter) %>%
  summarise(
    True = mean(True),
    Mean = mean(Mean),
    mean2.5 = mean(X2.5.),
    mean97.5 = mean(X97.5.),
    meanAbsBias=mean(AbsDiff),
    Percent.overlap = mean(CI)*100,
    Percent.zero = mean(overZero)*100)
tab[c(5:8,1:4,9:10),] %>% print.data.frame(digits = 2)


### save results:
# write.csv(tab[c(5:8,1:4,9:10),], row.names = F, quote=F, file="sim_res.csv")
# save(sim, res, file="sim_results_M=500.RData")
