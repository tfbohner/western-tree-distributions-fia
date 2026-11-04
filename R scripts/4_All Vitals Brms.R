## 4_All Vitals Brms
## brms versions of the occurrence and vital-rate models, used to cross-check the Stan code. Run after 4_All Vitals Stan Models.R loads the data.
## Author: Teresa Bohner
## Project: where western trees grow vs. where they thrive (see README)

## occurrence----
# data for brms
occurrence$swb <- scale(occurrence$ptcoeff)
prior <- c(prior_string("normal(0,10)", class = "b"))

pa.brms <- brm(pres~ -1 + intercept + swb + I(swb^2), data=occurrence, family=bernoulli(), 
               prior=prior, chains = 3, cores=3)
r2_all[i] <- bayes_R2(pa.brms)[1]
saveRDS(pa.brms, paste("Ch 1/saved models/pa models/pa.fit",i,"b.rds", sep=""))

## abundance----
# data for brms
abundance$swb <- scale(abundance$ptcoeff)
prior <- c(prior_string("normal(0,10)", class = "b"))

nb.brms <- brm(totabun~ -1 + intercept + swb + I(swb^2), data=abundance, family=negbinomial(),
               prior=prior, chains=3, cores=3)

r2_all[i] <- bayes_R2(nb.brms)[1]
saveRDS(nb.brms, paste("Ch 1/saved models/abun models/abun.fit",i,"b.rds", sep=""))


## growth ----
# data for brms
brms.dat <- dplyr::select(trees, c(gensp, baicm, plot.uniqueID)) %>%
  bind_cols(as.data.frame(cov.mat))
names(brms.dat)[4:8] <- c("dia","swb", "avgba", "adcount", "temp")
brms.dat$plot <- as.numeric(droplevels(as.factor(brms.dat$plot.uniqueID)))

grow.brms1 <- brm(baicm~ -1 + intercept + (1|plot) + dia + I(dia^2) +
                    swb + I(swb^2) +
                    temp + I(temp^2) +
                    adcount + avgba, 
                  family=gaussian(), 
                  data=brms.dat, chains=3, cores=3, prior = c(prior_(~normal(0,1), class=~sd),
                                                              prior_(~normal(0,10), class=~b),
                                                              prior_(~normal(0,10), class=~sigma)))

r2_all[i] <- bayes_R2(grow.brms1)[1]


saveRDS(grow.brms1, paste("Ch 1/saved models/grow models/grow.fit",i,"b.rds", sep=""))

## survival ----

brms.dat <- dplyr::select(surv, c(gensp, surv, plot.uniqueID)) %>%
  bind_cols(as.data.frame(cov.mat))
names(brms.dat)[4:8] <- c("dia","swb", "avgba", "adcount", "temp")
brms.dat$plot <- as.numeric(droplevels(as.factor(brms.dat$plot.uniqueID)))


fit.bern <- brm(bf(surv~ -1 + intercept + (1|plot) + dia + swb + avgba + adcount + temp + I(dia^2) + I(swb^2) + I(temp^2)), 
                data=brms.dat, family=bernoulli(), cores=4)

r2_all[i] <- bayes_R2(fit.bern)[1]


saveRDS(fit.bern, paste("Ch 1/saved models/surv models/surv.fit",i,"b.rds", sep=""))

## recruitment----
# data for brms
brms.dat <- dplyr::select(recruits, c(spid, N.sdlgs)) %>%
  cbind(cov.mat)
names(brms.dat)[3:6] <- c("swb", "avgba", "adcount", "temp")


r2_all <- rep(NA, length=max(brms.dat$spid))
for(i in 1:max(brms.dat$spid)){
  fit.zinb <- brm(bf(N.sdlgs ~ -1 + intercept + adcount + avgba + swb + temp + I(swb^2) + I(temp^2),zi~1),
                  data = filter(brms.dat, spid==i), family = zero_inflated_negbinomial())
  r2_all[i] <- bayes_R2(fit.zinb)[1]
}

r2_all <- cbind(data.frame("spid"=spptable$spid,
                           "gensp"=spptable$gensp), r2_all) %>%
  rename("r2_all"="x")
