##################################
## Supplementary to:
## Heck, Arndold, & Arndold (under revision).
## TreeBUGS: An R Package for Hierarchical
## Multinomial-Processing-Tree Modeling
##
## This script shows how to use TreeBUGS within R
## without any external files.
##################################


##################################
##### Setup

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



##################################
##### Handle MPT models, constraints, and data in R

# load data by Arnold et al. (2013) from package
?arnold2013
data(arnold2013)
head(arnold2013)
# get frequencies:
arnold.retrieval <- subset(arnold2013, group == "retrieval",
                           select=5:13)
arnold.encoding <- subset(arnold2013, group == "encoding",
                          select=5:13)

# define MPT model in R (cf. "model/2htsm.eqn")
htsm.eqn <- "###  2HTSM ###
E     EE     D1*d1
E     EE     D1*(1-d1)*a
E     EU     D1*(1-d1)*(1-a)
E     EE     (1-D1)*b*g
E     EU     (1-D1)*b*(1-g)
E     EN     (1-D1)*(1-b)
U     UU     D2*d2
U     UE     D2*(1-d2)*a
U     UU     D2*(1-d2)*(1-a)
U     UE     (1-D2)*b*g
U     UU     (1-D2)*b*(1-g)
U     UN     (1-D2)*(1-b)
N     NN     D3
N     NE     (1-D3)*b*g
N     NU     (1-D3)*b*(1-g)
N     NN     (1-D3)*(1-b)
"

# MPT equality constraints in R (cf. "model/restrictions.txt")
htsm.restr <- list("D1=D2=D3", "d1=d2", "g=a")

# check EQN input
readEQN(file = htsm.eqn, restrictions = htsm.restr, paramOrder = TRUE)

# [conversion of MPTinR model files into EQN files for TreeBUGS:]
#    library(MPTinR)
#    make.eqn(MPTinR.filename, eqn.filename)
#    readEQN(eqn.filename)


###################################
##### Testing participant heterogeneity

# chi-square test:
testHetChi(freq = arnold.retrieval,
           tree = c("E","E","E", "U","U","U", "N","N","N") )
testHetChi(freq = arnold.encoding,
           tree = c("E","E","E", "U","U","U", "N","N","N") )

### permutation test: (data not included in TreeBUGS)
# arnold.r.long <- read.csv("data/data_retrieval_long.csv")
# test.het <- testHetPerm(data = arnold.r.long,
#                         rep = 10000,
#                         tree = list(c("EE","EU","EN"),
#                                     c("UE","UU","UN"),
#                                     c("NE","NU","NN")))
# test.het[2:3]

# assess heterogeneity graphically
plotFreq(arnold.retrieval, eqn = htsm.eqn)
plotFreq(arnold.encoding, eqn = htsm.eqn)

# one line per individual
plotFreq(arnold.retrieval, eqn = htsm.eqn,
         freq = FALSE, boxplot = FALSE)



###################################
##### Model fitting

# Fitting a latent-trait MPT with defaults (not recommended):
m.retrieval.default  <- traitMPT(eqnfile=htsm.eqn,
                                 data = arnold.retrieval,
                                 restrictions = htsm.restr)
summary(m.retrieval.default)

# Fitting a latent-trait MPT with adjusted arguments:
m.retrieval <- traitMPT(eqnfile=htsm.eqn,
                        data = arnold.retrieval,
                        restrictions = htsm.restr,
                        covData = subset(arnold2013,
                                         group=="retrieval", "age"),
                        # modelfilename = "results/2htsm_traitMPT.jags",
                        transformedParameters = list("deltaDd=D1-d1"),
                        # parEstFile="results/results_retrieval_traitMPT.txt",
                        n.chain = 4, n.iter = 50000, n.adapt = 10000,
                        n.burnin = 10000, n.thin = 10,
                        ppp = 5000, dic = TRUE)
summary(m.retrieval)

# Fitting a beta-MPT model

m.retrieval.beta <- betaMPT(eqnfile=htsm.eqn,
                            data = arnold.retrieval,
                            restrictions = htsm.restr,
                            transformedParameters = list("deltaDd=D1-d1"),
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
m.retrieval3  <- traitMPT(eqnfile=htsm.eqn,
                          data = arnold.retrieval,
                          restrictions = htsm.restr,
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
m.encoding <- traitMPT(eqnfile=htsm.eqn,
                       data = arnold.encoding,
                       restrictions = htsm.restr,
                       covData = subset(arnold2013,
                                        group=="encoding", "age"),
                       transformedParameters = list("deltaDd=D1-d1"),
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

m.predictor <- traitMPT(eqnfile=htsm.eqn,
                        data = arnold.retrieval,
                        restrictions = htsm.restr,
                        covData = subset(arnold2013,
                                         group=="retrieval", "pc"),
                        predStructure = list("b ; pc"),
                        n.chain = 4, n.iter = 70000, n.adapt = 20000,
                        n.burnin = 15000, n.thin = 15)
summary(m.predictor)
# Note: unstandardized regression weights are reported!
#       (the latent-trait regression equation is:
#        theta_i = Phi( mu + slope*covariate_i + delta_i)
#     [MPT param.]    [mean]   [regression]    [person effect]



###################################
##### Including discrete predictors for between-subject manipulations

m.both.conditions <- traitMPT(eqnfile=htsm.eqn,
                              data = arnold2013[,5:13], # only frequencies
                              restrictions = htsm.restr,
                              covData = subset(arnold2013, select="group"),
                              predStructure = list("a b d1 D1 ; group"),
                              predType = c("f"),  # fixed effects
                              n.chain = 4, n.iter = 70000, n.adapt = 20000,
                              n.burnin = 15000, n.thin = 15)
summary(m.both.conditions)
round(getGroupMeans(m.both.conditions), 3)
