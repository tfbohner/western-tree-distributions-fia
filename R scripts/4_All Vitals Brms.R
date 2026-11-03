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
# get_prior(totabun~ -1 + intercept + swb + I(swb^2), data=abundance, family=negbinomial())
prior <- c(prior_string("normal(0,10)", class = "b"))

# make_stancode(totabun~ -1 + intercept + swb + I(swb^2), data=abundance, family=negbinomial())
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

# get_prior(baicm~ -1 + intercept + (1|plot) + dia + I(dia^2) +
#             swb + I(swb^2) +
#             temp + I(temp^2) +
#             adcount + avgba, 
#           family=gaussian(), 
#           data=brms.dat)
# 
# make_stancode(baicm~ -1 + intercept + (1|plot) + dia + I(dia^2) +
#             swb + I(swb^2) +
#             temp + I(temp^2) +
#             adcount + avgba, 
#           family=gaussian(), 
#           data=brms.dat, prior = c(prior_(~normal(0,1), class=~sd),
#                                    prior_(~normal(0,10), class=~sigma)))
r2_all[i] <- bayes_R2(grow.brms1)[1]

# plot((fixef(grow.brms1)[1] + fixef(grow.brms1)[4]*swb_pred + fixef(grow.brms1)[5]*swb_pred^2),
#      (summary(grow.fit, pars="a_sp")$summary[1] +
#         summary(grow.fit, pars="b_swb")$summary[1]*swb_pred +
#         summary(grow.fit, pars="b2_swb")$summary[1]*swb_pred^2))
# abline(0,1, col="red")

saveRDS(grow.brms1, paste("Ch 1/saved models/grow models/grow.fit",i,"b.rds", sep=""))

## survival ----

brms.dat <- dplyr::select(surv, c(gensp, surv, plot.uniqueID)) %>%
  bind_cols(as.data.frame(cov.mat))
names(brms.dat)[4:8] <- c("dia","swb", "avgba", "adcount", "temp")
brms.dat$plot <- as.numeric(droplevels(as.factor(brms.dat$plot.uniqueID)))

# make_stancode(bf(surv~ -1+ (1|plot) + dia + swb + avgba + adcount + temp + I(dia^2) + I(swb^2) + I(temp^2)), 
#               data=brms.dat, family=bernoulli())
# get_prior(bf(surv~ -1+ (1|plot) + dia + swb + avgba + adcount + temp + I(dia^2) + I(swb^2) + I(temp^2)), 
#               data=brms.dat, family=bernoulli())

fit.bern <- brm(bf(surv~ -1 + intercept + (1|plot) + dia + swb + avgba + adcount + temp + I(dia^2) + I(swb^2) + I(temp^2)), 
                data=brms.dat, family=bernoulli(), cores=4)

r2_all[i] <- bayes_R2(fit.bern)[1]

# plot(inv.logit(fixef(fit.bern)[1] + fixef(fit.bern)[3]*swb_pred + fixef(fit.bern)[8]*swb_pred^2),
#      inv.logit(pr$`50%`[pr$spid==i]))

saveRDS(fit.bern, paste("Ch 1/saved models/surv models/surv.fit",i,"b.rds", sep=""))

## recruitment----
# data for brms
brms.dat <- dplyr::select(recruits, c(spid, N.sdlgs)) %>%
  cbind(cov.mat)
names(brms.dat)[3:6] <- c("swb", "avgba", "adcount", "temp")

# make_stancode(bf(N.sdlgs ~ -1 + intercept + adcount + avgba + swb + temp + I(swb^2) + I(temp^2),zi~1),
#               data = filter(brms.dat, spid==1), family = zero_inflated_negbinomial())
# get_prior(bf(N.sdlgs ~ -1 + intercept + adcount + avgba + swb + temp + I(swb^2) + I(temp^2),zi~1),
#           data = filter(brms.dat, spid==1), family = zero_inflated_negbinomial())

r2_all <- rep(NA, length=max(brms.dat$spid))
for(i in 1:max(brms.dat$spid)){
  fit.zinb <- brm(bf(N.sdlgs ~ -1 + intercept + adcount + avgba + swb + temp + I(swb^2) + I(temp^2),zi~1),
                  data = filter(brms.dat, spid==i), family = zero_inflated_negbinomial())
  r2_all[i] <- bayes_R2(fit.zinb)[1]
}

r2_all <- cbind(data.frame("spid"=spptable$spid,
                           "gensp"=spptable$gensp), r2_all) %>%
  rename("r2_all"="x")
