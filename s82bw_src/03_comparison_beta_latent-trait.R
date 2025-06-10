
# install.packages("TreeBUGS")
# install.packages("dplyr")
library(TreeBUGS)
library(dplyr)


###################################
##### Fitting a latent-trait MPT model

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


###################################
##### Fitting a beta-MPT model

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

### comparison of parameter estimates
cors <- diag(cor(t(m.retrieval$summary$individParameters[,,"Mean"]),
                 t(m.retrieval.beta$summary$individParameters[,,"Mean"])))
absDiff <- rowMeans(abs(m.retrieval$summary$individParameters[,,"Mean"]-
                          m.retrieval.beta$summary$individParameters[,,"Mean"]))
mean_sd <- probitInverse(fittedModel = m.retrieval)
summ <- summarizeMCMC(mean_sd)
tab <- cbind(summ[seq(1,8,2),1:2],
             summ[seq(2,8,2),1:2],
             m.retrieval.beta$summary$groupParameters$mean[,1:2],
             m.retrieval.beta$summary$groupParameters$SD[,1:2],
             Cor=cors, absoluteDiff=absDiff)
colnames(tab)[c(1,3,5,7)] <- c("trait_mean","trait_sig","beta_mean","beta_SD")
round(tab,3)
write.csv(tab, file="comparison.csv", row.names = T, quote=F)

plotFit(m.retrieval.beta)
m.retrieval$summ$fitStatistics$overall[c(3,6)]
m.retrieval.beta$summ$fitStatistics$overall[c(3,6)]
m.retrieval$summary$dic
m.retrieval.beta$summary$dic
