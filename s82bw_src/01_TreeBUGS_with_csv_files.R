##################################
## Supplementary to:
## Heck, Arndold, & Arndold (under revision).
## TreeBUGS: An R Package for Hierarchical
## Multinomial-Processing-Tree Modeling
##
## This script shows how to use TreeBUGS with
## external files (EQN for defining MPT models,
## CSV for data)
##################################


##################################
##### Setup

# Change the working directory, e.g.:
# setwd("C:/mpt")

# Install and load TreeBUGS from CRAN:
install.packages("TreeBUGS")
# Install newest version from GitHub:
# install.packages("devtools")
# library(devtools)
# install_github("denis-arnold/TreeBUGS", build_vignettes = TRUE)

# load TreeBUGS
library(TreeBUGS)

# Access help files:
?TreeBUGS
?traitMPT
?betaMPT
vignette("TreeBUGS_1_intro")
vignette("TreeBUGS_2_extended")



###################################
##### Testing participant heterogeneity

# chi-square test:
testHetChi(freq = "data/data_retrieval.csv",
           tree = c("E","E","E", "U","U","U", "N","N","N") )
testHetChi(freq = "data/data_encoding.csv",
           tree = c("E","E","E", "U","U","U", "N","N","N") )

# permutation test:
test.het <- testHetPerm(data = "data/data_retrieval_long.csv",
                        rep = 10000,
                        tree = list(c("EE","EU","EN"),
                                    c("UE","UU","UN"),
                                    c("NE","NU","NN")))
test.het[2:3]

# assess heterogeneity graphically
plotFreq("data/data_retrieval.csv", eqn = "model/2htsm.eqn")
plotFreq("data/data_encoding.csv", eqn = "model/2htsm.eqn")

# one line per individual
plotFreq("data/data_retrieval.csv", eqn = "model/2htsm.eqn",
         freq = FALSE, boxplot = FALSE)



###################################
##### Model fitting

# Fitting a latent-trait MPT with defaults (not recommended):
m.retrieval.default  <- traitMPT(eqnfile="model/2htsm.eqn",
                                 data = "data/data_retrieval.csv",
                                 restrictions = "model/restrictions.txt")
summary(m.retrieval.default)

# Fitting a latent-trait MPT with adjusted arguments:
m.retrieval <- traitMPT(eqnfile="model/2htsm.eqn",
                        data = "data/data_retrieval.csv",
                        restrictions = "model/restrictions.txt",
                        covData = "data/age_retrieval.csv",
                        modelfilename = "results/2htsm_traitMPT.jags",
                        transformedParameters = list("deltaDd=D1-d1"),
                        parEstFile = "results/results_retrieval_traitMPT.txt",
                        n.chain = 4, n.iter = 50000, n.adapt = 10000,
                        n.burnin = 10000, n.thin = 10,
                        ppp = 5000, dic = TRUE)
summary(m.retrieval)

# Fitting a beta-MPT model

m.retrieval.beta <- betaMPT(eqnfile="model/2htsm.eqn",
                            data = "data/data_retrieval.csv",
                            restrictions = "model/restrictions.txt",
                            modelfilename = "results/2htsm_betaMPT.jags",
                            transformedParameters = list("deltaDd=D1-d1"),
                            parEstFile = "results/results_retrieval_betaMPT.txt",
                            n.chain = 4, n.iter = 50000, n.adapt = 10000,
                            n.burnin = 10000, n.thin = 10,
                            ppp = 5000, dic = TRUE
                            )
summary(m.retrieval.beta)



###################################
##### Monitoring convergence graphically

?plot.traitMPT
plot(m.retrieval, parameter = "mean")
plot(m.retrieval, parameter = "mean", type = "acf")
plot(m.retrieval, parameter = "mean", type = "gelman")

# other parameters:
plot(m.retrieval, parameter = "sigma")
plot(m.retrieval, parameter = "rho")

# continue MCMC sampling
m.retrieval2 <- extendMPT(m.retrieval.default,
                          n.iter = 5000, n.adapt = 3000)
summary(m.retrieval2)



###################################
##### Check and adjust prior distributions

##### beta-MPT
# (A) zeros-trick (Smith & Batchelder, 2010):
plotPrior(list(alpha="zero",beta="zero"))
# (B) TreeBUGS default
plotPrior(prior = list(alpha = "dgamma(1,.1)",
                       beta = "dgamma(1,.1)"))

###### traitMPT
# (A) TreeBUGS default
plotPrior(prior = list(mu = "dnorm(0,1)",
                       xi = "dunif(0,10)",
                       V=diag(2), df=3))
# (B) higher precision=4 for mean "mu"
plotPrior(prior = list(mu = c("dnorm(0,1)", "dnorm(0,4)"),
                       xi = "dunif(0,10)",
                       V=diag(2), df=3))

# Adjusting priors for group-level parameters
m.retrieval3  <- traitMPT(eqnfile="model/2htsm.eqn",
                          data = "data/data_retrieval.csv",
                          restrictions = "model/restrictions.txt",
                          mu = c( a="dnorm(0,4)", b="dnorm(0,4)",
                                  d1="dnorm(0,1)", D1="dnorm(0,1)"))
summary(m.retrieval3)



###################################
##### Testing model fit

# graphically:
plotFit(m.retrieval)
plotFit(m.retrieval, stat = "cov")

# posterior predictive p-values:
ppp.retrieval <- PPP(m.retrieval, M=1000)
ppp.retrieval

### posterior-predictive samples
# (A) for participants in data set
posteriorPredictive(m.retrieval, M=1, nCPU = 1)
# (B) for a novel participant (samples MPT parameters for person first)
posteriorPredictive(m.retrieval, M=5,
                    numItems=c(E=16,N=32,U=16), nCPU = 1)

posteriorPredictive(m.retrieval, M=5,
                    numItems=c(E=16,N=32,U=16), nCPU = 1)



###################################
##### Plotting and extracting parameters

# plot parameters:
plotParam(m.retrieval)

# plot estimated hierarchical distribution against plausible values
plotDistribution(m.retrieval, scale = "latent")

# compare prior to posterior density
plotPriorPost(m.retrieval)

# extract parameter estimates:
getParam(m.retrieval, parameter = "theta", stat = "mean")
getParam(m.retrieval, parameter = "mean", stat = "summary")
getParam(m.retrieval, parameter = "rho")



###################################
##### Within-subject tests

# Example: two-high threshold model (2HTM, included in TreeBUGS)
htm <- system.file("MPTmodels/2htm.eqn", package="TreeBUGS")

# check model equations and parameters
readEQN(htm, paramOrder=TRUE)

# create EQN file for within-subject manipulations
withinSubjectEQN(htm,
                 labels = c("high","low"),  # factor labels
                 constant=c("g"))           # constant parameters

# compute differences/ratios of MPT parameters within participants
transpar <- transformedParameters(fittedModel = m.retrieval,
                                  transformedParameters=list("deltaDd=D1-d1"),
                                  level = "individual")
summary(transpar[,1:5])$quantiles



###################################
##### Between-subject tests

# fit latent-trait MPT for second between-subject condition (encoding)
m.encoding <- traitMPT(eqnfile="model/2htsm.eqn",
                       data = "data/data_encoding.csv",
                       restrictions = "model/restrictions.txt",
                       modelfilename = "results/2htsm.jags",
                       covData = "data/age_encoding.csv",
                       transformedParameters = list("deltaDd=D1-d1"),
                       parEstFile = "results/results_encoding.txt",
                       n.chain = 4, n.iter = 50000, n.adapt = 10000,
                       n.burnin = 10000, n.thin = 10,
                       ppp = 5000, dic = TRUE)
summary(m.encoding)

# get (A) credibility intervals and (B) Bayesian p-values
# for the difference in parameters D and d:
betweenSubjectMPT(m.retrieval, m.encoding, par1 = "D1")
betweenSubjectMPT(m.retrieval, m.encoding, par1 = "d1")



###################################
##### Including continuous predictors

m.predictor <- traitMPT(eqnfile="model/2htsm.eqn",
                        data = "data/data_retrieval.csv",
                        restrictions = "model/restrictions.txt",
                        modelfilename = "results/2htsm_predictor.jags",
                        covData = "data/pc_retrieval.csv",
                        predStructure = list("b ; pc"),
                        parEstFile = "results/results_predictor.txt",
                        n.chain = 4, n.iter = 70000, n.adapt = 20000,
                        n.burnin = 15000, n.thin = 15)
summary(m.predictor)
# Note: unstandardized regression weights are reported!
#       (the latent-trait regression equation is:
#        theta_i = Phi( mu + slope*covariate_i + delta_i)
#     [MPT param.]    [mean]   [regression]    [person effect]



###################################
##### Including discrete predictors for between-subject manipulations

m.both.conditions <- traitMPT(eqnfile="model/2htsm.eqn",
                              data = "data/data_both.csv",
                              restrictions = "model/restrictions.txt",
                              modelfilename = "results/2htsm_between.jags",
                              covData = "data/group.csv",
                              predStructure = list("a b d1 D1 ; Group"),
                              predType = c("f"),  # fixed effects
                              parEstFile = "results/results_both.txt",
                              n.chain = 4, n.iter = 70000, n.adapt = 20000,
                              n.burnin = 15000, n.thin = 15)
summary(m.both.conditions)
round(getGroupMeans(m.both.conditions), 3)
