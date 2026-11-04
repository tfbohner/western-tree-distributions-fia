## 6_Pres abs with temp comparison
## Sensitivity check: occurrence models with soil water balance alone vs. soil water balance plus temperature.
## Author: Teresa Bohner
## Project: where western trees grow vs. where they thrive (see README)

library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(boot)
library(shinystan)
library(bayesplot)
library(KernSmooth)
library(outliers)
library(loo)
library(brms)
library(ggpmisc)

## Paths below follow the original project folder layout ("Ch 1/..."); point them at your data

#################### Occurence data manipulation----
rm(list = ls())
occurrence <- read.csv("Ch 1/data for stan/occurrence.new.csv")
spptable <- read.csv('Ch 1/data for stan/spptable.new.csv')
spptable <- spptable[,c(1:2)]
names(spptable)[2] <- "gensp"
spptable$gensp <- as.character(spptable$gensp)

length(unique(occurrence$plot.uniqueID))

tableinfo_occ <- dplyr::group_by(occurrence, gensp) %>%
  dplyr::summarize(nplots_p = sum(pres, na.rm=T))


N_spp <- length(unique(occurrence$spid))
maxmin <- dplyr::group_by(occurrence, spid)%>%
  dplyr::summarise(min_swb = min(ptcoeff, na.rm=T),
            max_swb = max(ptcoeff, na.rm=T))
theta_all <- matrix(NA, nrow = 200*N_spp, ncol=6)
theta_all <- data.frame(theta_all)

maxminopt_all_med <- data.frame(matrix(NA, nrow=N_spp, ncol=17))
r2_all <- rep(NA, length=N_spp)

occurrence1 <- occurrence
maxmin1 <- maxmin

## Occurrence loop----
for(i in 1:N_spp){
  occurrence <- occurrence1[occurrence1$spid==i,]
  maxmin <- maxmin1[maxmin1$spid==i,]
  
  swb.raw <- occurrence$ptcoeff
  swb.mean <- mean(swb.raw)
  swb.sd <- sd(swb.raw)
  swb <- as.vector(scale(swb.raw))
  
  temp.raw <- occurrence$temp
  temp.mean <- mean(temp.raw)
  temp.sd <- sd(temp.raw)
  temp <- as.vector(scale(temp.raw))
  
  min_swb <- (maxmin$min_swb-swb.mean)/swb.sd
  max_swb <- (maxmin$max_swb-swb.mean)/swb.sd
  
  N <- dim(occurrence)[1]
  N_pred <- 200
  N_plots <- length(unique(occurrence$plot.uniqueID))
  plot <- as.vector(as.integer(factor(occurrence$plot.uniqueID)))
  pres <- as.vector(occurrence$pres)
  swb <- swb
  temp <- temp
  min_swb <- min_swb
  max_swb <- max_swb
  swb_pred <- seq(min(swb),max(swb), length.out=N_pred)

  pa_data <- list(N=N, N_plots=N_plots, plot=plot, N_pred=N_pred,
                  pres=pres, swb=swb, swb_pred=swb_pred, temp=temp, 
                  min_swb=min_swb, max_swb=max_swb)

  pars=c("a","b_swb","b2_swb", "b_temp","b2_temp",
         "opt_swb","thetamin_swb","thetamax_swb","thetaopt_swb",
         "theta_pred") # "y_pred","log_lik"

  ## Occurence models----
  rstan_options(auto_write = TRUE) 
  pa.compiled <- stan_model(file="Stan/presence_absence/binomial_swb_temp_spp.stan")
  pa.fit<- sampling(
    object = pa.compiled
    , data=pa_data
    , chains=3
    , iter=2000
    , warmup=1000
    , thin=1
    , refresh=100
    , pars=pars
    , cores=3
    # , control= list(adapt_delta=0.9, max_treedepth=18)
  )

  
  ## Occurrence diagnostics----
  
  ## Occurrence output----
  theta <- as.data.frame(summary(pa.fit, pars="theta_pred")$summary[, c(1,4,6,8)])
  spid <- rep(i, N_pred)
  theta$spid <- spid
  theta <- merge(theta, spptable, by='spid')
  
  theta_all[c((i*200-200+1):(i*200)),] <- theta
  names(theta_all) <- names(theta)
  
  # output data for VR modeling
  # median and intervals
  maxminopt <- maxmin
  names(maxminopt) <- c('spid','min_swb_p','max_swb_p')
  opt_swb <- summary(pa.fit, pars="opt_swb")$summary[, 6]
  opt_swb_2.5 <- summary(pa.fit, pars="opt_swb")$summary[, 4]
  opt_swb_25 <- summary(pa.fit, pars="opt_swb")$summary[, 5]
  opt_swb_75 <- summary(pa.fit, pars="opt_swb")$summary[, 7]
  opt_swb_97.5 <- summary(pa.fit, pars="opt_swb")$summary[, 8]
  maxminopt$peak_swb_p <- swb.mean + opt_swb*swb.sd
  maxminopt$peak_swb_p_2.5 <- swb.mean + opt_swb_2.5*swb.sd
  maxminopt$peak_swb_p_25 <- swb.mean + opt_swb_25*swb.sd
  maxminopt$peak_swb_p_75 <- swb.mean + opt_swb_75*swb.sd
  maxminopt$peak_swb_p_97.5 <- swb.mean + opt_swb_97.5*swb.sd
  maxminopt$thetamin_swb_p <- inv.logit(summary(pa.fit, pars="thetamin_swb")$summary[, 6])
  maxminopt$thetamax_swb_p <- inv.logit(summary(pa.fit, pars="thetamax_swb")$summary[, 6])
  maxminopt$thetapeak_swb_p <- inv.logit(summary(pa.fit, pars="thetaopt_swb")$summary[, 6])
  maxminopt$peak_swb_p <- ifelse(maxminopt$min_swb_p<=maxminopt$peak_swb_p&maxminopt$peak_swb_p<=maxminopt$max_swb_p, maxminopt$peak_swb_p,NA)
  maxminopt$thetapeak_swb_p <- ifelse(!is.na(maxminopt$peak_swb_p), maxminopt$thetapeak_swb_p,NA)
  maxminopt$opt_swb_p <- ifelse(!is.na(maxminopt$thetapeak_swb_p),ifelse(maxminopt$thetapeak_swb_p>maxminopt$thetamin_swb_p & maxminopt$thetapeak_swb_p>maxminopt$thetamax_swb_p,
                                                                         maxminopt$peak_swb_p,
                                                                         ifelse(maxminopt$thetamin_swb_p>maxminopt$thetamax_swb_p, 
                                                                                maxminopt$min_swb_p, maxminopt$max_swb_p)),
                                ifelse(maxminopt$thetamin_swb_p>maxminopt$thetamax_swb_p, maxminopt$min_swb_p, maxminopt$max_swb_p))
  maxminopt$opt_swb_p_2.5 <- ifelse(!is.na(maxminopt$thetapeak_swb_p),ifelse(maxminopt$thetapeak_swb_p>maxminopt$thetamin_swb_p & maxminopt$thetapeak_swb_p>maxminopt$thetamax_swb_p,
                                                                            maxminopt$peak_swb_p_2.5,
                                                                            ifelse(maxminopt$thetamin_swb_p>maxminopt$thetamax_swb_p, 
                                                                                   maxminopt$min_swb_p, maxminopt$max_swb_p)),
                                   ifelse(maxminopt$thetamin_swb_p>maxminopt$thetamax_swb_p, maxminopt$min_swb_p, maxminopt$max_swb_p))
  maxminopt$opt_swb_p_25 <- ifelse(!is.na(maxminopt$thetapeak_swb_p),ifelse(maxminopt$thetapeak_swb_p>maxminopt$thetamin_swb_p & maxminopt$thetapeak_swb_p>maxminopt$thetamax_swb_p,
                                                                            maxminopt$peak_swb_p_25,
                                                                            ifelse(maxminopt$thetamin_swb_p>maxminopt$thetamax_swb_p, 
                                                                                   maxminopt$min_swb_p, maxminopt$max_swb_p)),
                                   ifelse(maxminopt$thetamin_swb_p>maxminopt$thetamax_swb_p, maxminopt$min_swb_p, maxminopt$max_swb_p))
  maxminopt$opt_swb_p_75 <- ifelse(!is.na(maxminopt$thetapeak_swb_p),ifelse(maxminopt$thetapeak_swb_p>maxminopt$thetamin_swb_p & maxminopt$thetapeak_swb_p>maxminopt$thetamax_swb_p,
                                                                            maxminopt$peak_swb_p_75,
                                                                            ifelse(maxminopt$thetamin_swb_p>maxminopt$thetamax_swb_p, 
                                                                                   maxminopt$min_swb_p, maxminopt$max_swb_p)),
                                   ifelse(maxminopt$thetamin_swb_p>maxminopt$thetamax_swb_p, maxminopt$min_swb_p, maxminopt$max_swb_p))
  maxminopt$opt_swb_p_97.5 <- ifelse(!is.na(maxminopt$thetapeak_swb_p),ifelse(maxminopt$thetapeak_swb_p>maxminopt$thetamin_swb_p & maxminopt$thetapeak_swb_p>maxminopt$thetamax_swb_p,
                                                                            maxminopt$peak_swb_p_97.5,
                                                                            ifelse(maxminopt$thetamin_swb_p>maxminopt$thetamax_swb_p, 
                                                                                   maxminopt$min_swb_p, maxminopt$max_swb_p)),
                                   ifelse(maxminopt$thetamin_swb_p>maxminopt$thetamax_swb_p, maxminopt$min_swb_p, maxminopt$max_swb_p))
  
  maxminopt$thetaopt_swb_p <- ifelse(!is.na(maxminopt$thetapeak_swb_p),ifelse(maxminopt$thetapeak_swb_p>maxminopt$thetamin_swb_p & maxminopt$thetapeak_swb_p>maxminopt$thetamax_swb_p,
                                                                              maxminopt$thetapeak_swb_p,
                                                                              ifelse(maxminopt$thetamin_swb_p>maxminopt$thetamax_swb_p, 
                                                                                     maxminopt$thetamin_swb_p, maxminopt$thetamax_swb_p)),
                                     ifelse(maxminopt$thetamin_swb_p>maxminopt$thetamax_swb_p, maxminopt$thetamin_swb_p, maxminopt$thetamax_swb_p))
  
  maxminopt_all_med[i,] <- maxminopt
  names(maxminopt_all_med) <- names(maxminopt)
  
  print(paste("spp[",i,"] complete"), sep="")
  
}

maxminopt_p <- maxminopt_all_med
theta <- theta_all
occurrence <- occurrence1

theta2 <- theta
maxminopt_p2 <- maxminopt_p

spptable2 <- dplyr::select(theta, c(spid, gensp)) %>% group_by(spid, gensp) %>% summarize()
maxminopt_p2 <- merge(maxminopt_p2, spptable2, by="spid")
maxminopt_p <- merge(maxminopt_p, spptable, by="spid")

## Occurrence figures ----
trans_gray <- adjustcolor(col = "gray",alpha.f = 0.4)


load("Ch 1/saved models/distribution modeling swb.RData")

pdf("Output pdfs/pres by species_partialvfull.pdf")
par(mfrow=c(2,2))
for(i in 1:N_spp){
  xdum <- seq(min(occurrence$ptcoeff[occurrence$spid==i],na.rm=T),max(occurrence$ptcoeff[occurrence$spid==i],na.rm=T), length.out=N_pred)
  xmean <- mean(occurrence$ptcoeff[occurrence$spid==i],na.rm=T)
  xsd <- sd(occurrence$ptcoeff[occurrence$spid==i],na.rm=T)
  
  ncuts<-20
  chopsize<-cut(occurrence$ptcoeff[occurrence$spid==i],ncuts)
  
  binmean.size<-as.numeric(sapply(split(occurrence$ptcoeff[occurrence$spid==i],chopsize),mean, na.rm=T))
  binmean.occur<-as.numeric(sapply(split(occurrence$pres[occurrence$spid==i],chopsize),mean, na.rm=T))
  plot(occurrence$ptcoeff[occurrence$spid==i], occurrence$pres[occurrence$spid==i], type="n", 
       main=spptable[i,2],
       xlab="swb gradient", 
       ylab="pr(occurrence)")
  points(binmean.size, binmean.occur)
  mtext("wet", side=1, line=2, adj=1)
  mtext("dry", side=1, line=2, adj=0)
  
  polygon(c(xdum,rev(xdum)),c(inv.logit(theta$`97.5%`[theta$spid==i]),
                              rev(inv.logit(theta$`2.5%`[theta$spid==i]))),col=trans_gray, border=NA)
  lines(xdum, inv.logit(theta$`50%`[theta$spid==i]), col=4)
  lines(xdum, inv.logit(theta2$`50%`[theta$spid==i]), col=2)

  abline(v=maxminopt_all_med$opt_swb_p[maxminopt_all_med$spid==i], col=4, lty=2)
  abline(v=maxminopt_p2$opt_swb_p[maxminopt_p2$spid==i], col=2, lty=2)
  
  legend("topright",
         legend=c("without temperature","with temperature"),
         lty=c(1,1), col=c("blue","red"), cex=0.85)
  

}

dev.off()

comparison <- merge(dplyr::select(maxminopt_p2, c(gensp, opt_swb_p)), dplyr::select(maxminopt_p, c(gensp, opt_swb_p)), by='gensp') %>%
  rename(opt_swb_temp=opt_swb_p.x, opt_swb=opt_swb_p.y)
pdf("Output pdfs/jan19 abundance figures/occurrence comparision.pdf", width = 5, height = 5)
ggplot(comparison, aes(opt_swb, opt_swb_temp)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method='lm', se=FALSE, color='black') +
  # geom_abline(intercept=0,slope=1, color="red") +
  ylab("SWB optimium with temperature") + xlab("SWB optimum without temperature") +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.05, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = .15, label.x.npc = 0.85, size = 4) +
  ggtitle("Occurrence models")

dev.off()

plot(maxminopt_p$opt_swb_p, maxminopt_p2$opt_swb_p,
     xlab="SWB optimum without temperature",
     ylab="SWB optimum including temperature")
abline(0,1)


wilcox.test(comparison$opt_swb_temp, comparison$opt_swb)

save.image('Ch 1/saved models/distribution modeling comparison.RData', compress=T)


## Abundance data manipulation----
## read in data
rm(list = ls())
abundance <- read.csv("Ch 1/data for stan/abundance.new.csv")
spptable <- read.csv('Ch 1/data for stan/spptable.new.csv')
spptable <- spptable[,c(1:2)]
names(spptable)[2] <- "gensp"
spptable$gensp <- as.character(spptable$gensp)
maxminopt <- read.csv('Ch 1/data for stan/distromaxminopt_swb_temp.csv')
maxminopt <- merge(maxminopt, spptable, by="gensp")
theta <- read.csv('Ch 1/data for stan/theta_swb_temp.csv')
theta$seq <- as.numeric(rownames(theta))
theta <- merge(theta, spptable, by="gensp")
theta <- theta[order(theta$seq),]
swb_axes_p <- read.csv('Ch 1/data for stan/swb_axes_p_temp.csv')
swb_axes_p <- merge(swb_axes_p, spptable, by="gensp")

abundance$avgba <- abundance$totba/abundance$totabun

maxmin <- dplyr::group_by(abundance, spid)%>%
  dplyr::summarise(min_swb = min(ptcoeff, na.rm=T),
                   max_swb = max(ptcoeff, na.rm=T))

N_spp <- length(unique(abundance$spid))

eta_all <- matrix(NA, nrow = 200*N_spp, ncol=6)
eta_all <- data.frame(eta_all)

med_SWBs_all <- data.frame(matrix(NA, nrow=N_spp, ncol=22))
r2_all <- rep(NA, length=N_spp)

abund1 <- abundance
maxmin1 <- maxmin
maxminopt1 <- maxminopt


spids <- dplyr::group_by(abund1, gensp) %>%
  dplyr::summarise(spid = unique(spid))
options(tibble.print_max = Inf)
spids

#### Abundance loop----
for(i in 34:N_spp){
  abundance <- abund1[abund1$spid==i,]
  maxmin <- maxmin1[maxmin1$spid==i,]
  maxminopt <- maxminopt1[maxminopt1$spid==i,]
  
  ## scale covariates 
  cov.mat.raw <- cbind(abundance$ptcoeff, abundance$temp, abundance$avgba)## data for Stan , trees$dist_nowCWD, trees$dist_histCWD
  cov.mat <- apply(cov.mat.raw, 2, scale)
  cov.means <- colMeans(cov.mat.raw, na.rm=T)
  cov.sds <- apply(cov.mat.raw, 2, sd, na.rm=T)
  
  min_swb <- (maxmin$min_swb-cov.means[1])/cov.sds[1]
  max_swb <- (maxmin$max_swb-cov.means[1])/cov.sds[1]
  
  ## put opt on same swb scale
  opt_swb_p <- (maxminopt$opt_swb_p-cov.means[1])/cov.sds[1]
  
  # data for stan
  N <- dim(abundance)[1]
  N_pred <- 200
  N_plts <- length(unique(trees$plot.uniqueID))
  plot <- as.numeric(droplevels(as.factor(trees$plot.uniqueID)))
  spp <- abundance$spid
  
  abun <- as.vector(abundance$totabun)
  
  swb <- cov.mat[,1]
  temp <- cov.mat[,2]
  swb_pred <- seq(min(cov.mat[,1],na.rm=T),max(cov.mat[,1],na.rm=T), length.out=N_pred)
  min_swb <- min_swb
  max_swb <- max_swb
  pres_swb_dat <- opt_swb_p
  avgba <- cov.mat[,3]
  
  abun_data <- list(N=N, N_pred=N_pred,
                    abun=abun, swb=swb, swb_pred=swb_pred, 
                    temp=temp, #avgba=avgba,
                    pres_swb_dat=pres_swb_dat, min_swb=min_swb, max_swb=max_swb)
  
  pars=c("a","b_swb","b2_swb","b_temp","b2_temp","phi"
         ,"wet_swb","dry_swb","pres_swb","opt_swb", "par_swb"
         ,"n_hat_wet","n_hat_dry","n_hat_opt", "n_hat_pres"
         ,"eta_pred") # ,'y_rep' ,"b_ba",
  
  ## Abundance models---- 
  rstan_options(auto_write = TRUE) 
  nb.compiled <- stan_model(file="Stan/abundance/abundance_new_spp_nb_temp.stan")
  
  nb.fit<- sampling(
    object = nb.compiled
    , data=abun_data
    , chains=3
    , iter=1000
    , warmup=500
    , thin=1
    , refresh=100
    , pars=pars
    , cores=3
    , control= list(adapt_delta=0.9)
  )
  
  # data for brms
  # # get_prior(totabun~ -1 + intercept + swb + I(swb^2), data=abundance, family=negbinomial())
  # 
  # # make_stancode(totabun~ -1 + intercept + swb + I(swb^2), data=abundance, family=negbinomial())
  # nb.brms <- brm(totabun~ -1 + intercept + swb + I(swb^2), data=abundance, family=negbinomial(),
  #                prior=prior, chains=3, cores=3)
  # 
  
  # ## Abundance diagnostics----
  
  
  ## Abundance output----
  # predicted curve
  eta <- as.data.frame(summary(nb.fit, pars="eta_pred")$summary[, c(1,4,6,8)])
  spid <- rep(i, N_pred)
  eta$spid <- spid
  eta <- merge(eta, spptable, by='spid')
  
  if(i==1) eta_all <- eta
  if(i>1) eta_all <- rbind(eta_all, eta)
  
  ## betas
  betas <- as.data.frame(summary(nb.fit, pars=c("b_swb","b2_swb"),probs = c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))$summary[, c(1:10, 12)])
  betas$param <- rownames(betas)
  
  betas$spid <- i
  betas <- merge(betas, spptable, by='spid')
  
  if(i==1) betas_all <- betas
  if(i>1) betas_all <- rbind(betas_all, betas)
  
  ## SWB and subsequent g-hat values
  abun_swb <- as.data.frame(summary(nb.fit, pars="opt_swb")$summary[, c(1,4,5,6,7,8)])
  abun_swb_unscaled <- cov.means[1] + abun_swb*cov.sds[1]
  par_swb <- as.data.frame(summary(nb.fit, pars="par_swb")$summary[,1])
  par_swb_unscaled <- cov.means[1] + par_swb*cov.sds[1]
  pres_swb <- as.data.frame(summary(nb.fit, pars="pres_swb")$summary[,1])
  pres_swb_unscaled <- cov.means[1] + pres_swb*cov.sds[1]
  wet_swb <- as.data.frame(summary(nb.fit, pars="wet_swb")$summary[,1])
  wet_swb_unscaled <- cov.means[1] + wet_swb*cov.sds[1]
  dry_swb <- as.data.frame(summary(nb.fit, pars="dry_swb")$summary[,1])
  dry_swb_unscaled <- cov.means[1] + dry_swb*cov.sds[1]
  
  n_hat_abun <- as.data.frame(summary(nb.fit, pars="n_hat_opt")$summary[, c(1,4,6,8)])
  n_hat_pres <- as.data.frame(summary(nb.fit, pars="n_hat_pres")$summary[, c(1,4,6,8)])
  n_hat_wet <- as.data.frame(summary(nb.fit, pars="n_hat_wet")$summary[, c(1,4,6,8)])
  n_hat_dry <- as.data.frame(summary(nb.fit, pars="n_hat_dry")$summary[, c(1,4,6,8)])
  
  # Medians
  SWBs <- as.data.frame(cbind(as.list(abun_swb_unscaled[2:6,1]), pres_swb_unscaled, wet_swb_unscaled, dry_swb_unscaled, 
                              n_hat_abun[3,1], n_hat_pres[3,1], n_hat_wet[3,1], n_hat_dry[3,1]))
  names(SWBs) <- c('abun_swb_2.5','abun_swb_25','abun_swb','abun_swb_75','abun_swb_97.5','pres_swb','wet_swb','dry_swb','n_hat_abun', 'n_hat_pres','n_hat_wet','n_hat_dry')
  SWBs$spid <- i
  
  # Choose greater: parabola or endpoint optimum
  SWBs$abun_swb <- ifelse(SWBs$wet_swb>=SWBs$abun_swb&SWBs$abun_swb>=SWBs$dry_swb, SWBs$abun_swb,NA)
  SWBs$n_hat_abun <- ifelse(!is.na(SWBs$abun_swb), SWBs$n_hat_abun,NA)
  SWBs$opt_swb <- ifelse(!is.na(SWBs$n_hat_abun),ifelse(SWBs$n_hat_abun>SWBs$n_hat_wet & SWBs$n_hat_abun>SWBs$n_hat_dry
                                                        , SWBs$abun_swb, 
                                                        ifelse(SWBs$n_hat_dry>SWBs$n_hat_wet 
                                                               , SWBs$dry_swb, SWBs$wet_swb)), 
                         ifelse(is.na(SWBs$n_hat_abun),ifelse(SWBs$n_hat_dry>SWBs$n_hat_wet 
                                                              , SWBs$dryswb, SWBs$wetswb)))
  
  SWBs$n_hat_opt <- ifelse(!is.na(SWBs$n_hat_abun),ifelse(SWBs$n_hat_abun>SWBs$n_hat_wet & SWBs$n_hat_abun>SWBs$n_hat_dry
                                                          , SWBs$n_hat_abun, 
                                                          ifelse(SWBs$n_hat_dry>SWBs$n_hat_wet 
                                                                 , SWBs$n_hat_dry, SWBs$n_hat_wet)), 
                           ifelse(is.na(SWBs$n_hat_abun),ifelse(SWBs$n_hat_dry>SWBs$n_hat_wet 
                                                                , SWBs$n_hat_dry, SWBs$n_hat_wet)))
  
  ## Choose minimum for sensitivity
  SWBs$min_swb <- ifelse(!is.na(SWBs$n_hat_abun),ifelse(SWBs$n_hat_abun<SWBs$n_hat_wet & SWBs$n_hat_abun<SWBs$n_hat_dry
                                                        , SWBs$abun_swb, 
                                                        ifelse(SWBs$n_hat_dry<SWBs$n_hat_wet 
                                                               , SWBs$dry_swb, SWBs$wet_swb)), 
                         ifelse(is.na(SWBs$n_hat_abun),ifelse(SWBs$n_hat_dry<SWBs$n_hat_wet 
                                                              , SWBs$dryswb, SWBs$wetswb)))
  
  SWBs$n_hat_min <- ifelse(!is.na(SWBs$n_hat_abun),ifelse(SWBs$n_hat_abun<SWBs$n_hat_wet & SWBs$n_hat_abun<SWBs$n_hat_dry
                                                          , SWBs$n_hat_abun, 
                                                          ifelse(SWBs$n_hat_dry<SWBs$n_hat_wet 
                                                                 , SWBs$n_hat_dry, SWBs$n_hat_wet)), 
                           ifelse(is.na(SWBs$n_hat_abun),ifelse(SWBs$n_hat_dry<SWBs$n_hat_wet 
                                                                , SWBs$n_hat_dry, SWBs$n_hat_wet)))
  # Calculate mismatch and cost
  SWBs$mismatch_p <- SWBs$opt_swb - SWBs$pres_swb
  
  SWBs$cost_p <- SWBs$n_hat_opt - SWBs$n_hat_pres
  SWBs$cost_p <- ifelse(SWBs$cost_p<0, 0, SWBs$cost_p)

  
  SWBs$range <- SWBs$n_hat_opt - SWBs$n_hat_min
  SWBs$sens <- SWBs$range/(SWBs$wet_swb - SWBs$dry_swb)
  
  SWBs <- merge(SWBs, spptable, by='spid')
  
  med_SWBs_all[i,] <- SWBs
  names(med_SWBs_all) <- names(SWBs)
  
  print(paste("spp[",i,"] complete"), sep="")
  
}


## Abundance figures----
abunswbs_2 <- med_SWBs_all
eta2 <- eta_all

load("Ch 1/saved models/abundance loop_temp.RData")

trans_gray <- adjustcolor(col = "gray",alpha.f = 0.4)


pdf("Output pdfs/abun by species_fullvpartial.pdf")
par(mfrow=c(2,2))
for(i in 1:N_spp){
  xdum <- seq(min(abund1$ptcoeff[abund1$spid==i],na.rm=T),max(abund1$ptcoeff[abund1$spid==i],na.rm=T), length.out=N_pred)
  xmean <- mean(abund1$ptcoeff[abund1$spid==i],na.rm=T)
  xsd <- sd(abund1$ptcoeff[abund1$spid==i],na.rm=T)
  ncuts<-20
  chopsize<-cut(abund1$ptcoeff[abund1$spid==i],ncuts)
  
  ###start here
  
  binmean.size<-as.numeric(sapply(split(abund1$ptcoeff[abund1$spid==i],chopsize),mean, na.rm=T))
  binmean.abun<-as.numeric(sapply(split(abund1$totabun[abund1$spid==i],chopsize),mean, na.rm=T))
  plot(abund1$ptcoeff[abund1$spid==i], abund1$totabun[abund1$spid==i], type="n", ylim=c(min(binmean.abun, na.rm=T)-5, max(binmean.abun, na.rm=T)+10),
       main=spptable[i,2],
       xlab="swb gradient", 
       ylab="Abundance")
  points(binmean.size, binmean.abun)
  mtext("wet", side=1, line=2, adj=1)
  mtext("dry", side=1, line=2, adj=0)
  
  polygon(c(xdum,rev(xdum)),c(eta_all$`97.5%`[eta_all$spid==i],rev(eta_all$`2.5%`[eta_all$spid==i])),col=trans_gray, border=NA)
  lines(xdum, eta_all$`50%`[eta_all$spid==i], col=4)
  lines(xdum, eta2$`50%`[eta2$spid==i], col=2)
  abline(v=abunswbs$opt_swb[abunswbs$spid==i], col=4, lty=2)
  abline(v=abunswbs_2$opt_swb[abunswbs_2$spid==i], col=2, lty=2)
  
  legend("topright",
         legend=c("without temperature","with temperature"),
         lty=c(1,1), col=c("blue","red"), cex=0.85)
  
}
dev.off()

comparison <- merge(dplyr::select(abunswbs_2, c(gensp, opt_swb)), dplyr::select(abunswbs, c(gensp, opt_swb)), by='gensp') %>%
  rename(opt_swb_temp=opt_swb.x, opt_swb=opt_swb.y)
pdf("Output pdfs/jan19 abundance figures/abundance comparision.pdf", width = 5, height = 5)
ggplot(comparison, aes(opt_swb, opt_swb_temp)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method='lm', se=FALSE, color='black') +
  # geom_abline(intercept=0,slope=1, color="red") +
  ylab("SWB optimium with temperature") + xlab("SWB optimum without temperature") +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.05, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = .15, label.x.npc = 0.85, size = 4) +
  ggtitle("Density models")
dev.off()

wilcox.test(comparison$opt_swb_temp, comparison$opt_swb)

save.image('Ch 1/saved models/abundance modeling comparison.RData', compress=T)

### Growth models----
rm(list = ls())
load("Ch 1/saved models/growth loop_temp.RData")
## for trying out different distributions
i=1

trees <- trees1[trees1$spid==i,]

## scale covariates 
cov.mat.raw <- cbind(trees$dia, trees$ptcoeff, trees$avgba, trees$adcount, trees$temp)## data for Stan , trees$dist_nowCWD, trees$dist_histCWD
cov.mat <- apply(cov.mat.raw, 2, scale)
cov.means <- colMeans(cov.mat.raw, na.rm=T)
cov.sds <- apply(cov.mat.raw, 2, sd, na.rm=T)

# data for brms
brms.dat <- dplyr::select(trees, c(gensp, baicm, plot.uniqueID)) %>%
  bind_cols(as.data.frame(cov.mat.raw)) 
#%>%
names(brms.dat)[4:8] <- c("dia","swb", "avgba", "adcount", "temp")
brms.dat$plot <- as.numeric(droplevels(as.factor(brms.dat$plot.uniqueID)))

constant <- abs(min(brms.dat$baicm )) + 0.01
## add constant to make bai all positive (for positive distributions)
brms.dat$posbai <- brms.dat$baicm +constant
brms.dat$logbai <- log(brms.dat$posbai)

boxcox(brms.dat$posbai, optimize = T)

powerTransform(brms.dat$baicm, family='bcnPower')
brms.dat$bcpos <- bcPower(brms.dat$posbai, lambda=.27)


mod <- lm(posbai + 1 ~ dia + swb, data=brms.dat)
mod.aux <- update(mod, . ~ . + boxCoxVariable(posbai + 1))
summary(mod.aux)


## box cox transformations with negatives
brms.dat$bcbai <- bcnPower(brms.dat$baicm, lambda=0.21, gamma=0.1) ## don't exactly know how to chooose lambda and gamma yet
brms.dat$yjbai <- yjPower(brms.dat$baicm, 0)
hist(brms.dat$bcbai)
hist(brms.dat$yjbai) ## this method doesn't seem to help

## Normal error distributions----
## all data, no transormations
grow.brms.norm <- brm(bf(baicm~ -1 + intercept + (1|plot) + dia + I(dia^2) +
                    swb + I(swb^2) +
                    temp + I(temp^2) +
                    adcount + avgba), 
                  family=gaussian()
                 , data=brms.dat
                 , chains=3, cores=3
                 )

## all data, boxcox transformation
grow.brms.norm.bc <- brm(bf(bcbai~ -1 + intercept + (1|plot) + dia + I(dia^2) +
                           swb + I(swb^2) +
                           temp + I(temp^2) +
                           adcount + avgba), 
                      family=gaussian()
                      , data=brms.dat
                      , chains=3, cores=3
)


grow.brms.norm.bcp <- brm(bf(bcpos~ -1 + intercept + (1|plot) + dia + I(dia^2) +
                              swb + I(swb^2) +
                              temp + I(temp^2) +
                              adcount + avgba), 
                         family=gaussian()
                         , data=brms.dat
                         , chains=3, cores=3
)

## all data, no transformation, sigma varies with diameter
  ## this one actually looks pretty good now that I've changed over to standardized residuals (recommended)
stanvars <- stanvar(scode = "  
                    real par_swb;
    par_swb = -b[1]/(2*b[2]);", 
                    block = "genquant")
grow.brms.normv <- brm(bf(baicm~ -1 + intercept + (1|plot) + swb + I(swb^2) +
                            dia + I(dia^2) +
                           temp + I(temp^2) +
                           adcount + avgba, sigma~dia), 
                      family=gaussian()
                      , data=brms.dat
                      , chains=3, cores=3
)

grow.1 <- brm(bf(baicm~ (1|plot) + swb + I(swb^2) +
                   dia + I(dia^2) +
                   temp + I(temp^2) +
                   adcount + avgba, sigma~dia), 
                       family=gaussian()
                       , data=brms.dat
                       , chains=3, cores=3
)


list_of_data <- marginal_effects(grow.1, resolution=1000)
par <- list_of_data$swb[which.max(list_of_data$swb$estimate__),1:7]
par2 <- par
par2$swb <- maxminopt$opt_swb_p
pars <- rbind(par, par2)
fitted(grow.1, newdata=pars)
draws.opt <- fitted(grow.1, newdata=pars, summary=F)
mean(draws.opt[,1]-draws.opt[,2])
pred <- predict(grow.1, newdata = list_of_data$swb, summary = T)
fit <- fitted(grow.1, newdata = list_of_data$swb, summary = T)

## all data, log + constant transformation
grow.brms.norml <- brm(bf(logbai~ -1 + intercept + (1|plot) + dia + I(dia^2) +
                            swb + I(swb^2) +
                            temp + I(temp^2) +
                            adcount + avgba), 
                       family=gaussian()
                       , data=brms.dat
                       , chains=3, cores=3
)

## Exponential modified gaussian error distribution----
## all data, no transformation
grow.brms.exg <- brm(bf(baicm~ -1 + intercept + (1|plot) + dia + I(dia^2) +
                           swb + I(swb^2) +
                           temp + I(temp^2) +
                           adcount + avgba), 
                      family=exgaussian()
                      , data=brms.dat
                      , chains=3, cores=3
)

## Skewed normal error distribution----
grow.brms.sn <- brm(bf(baicm~ -1 + intercept + (1|plot) + dia + I(dia^2) +
                          swb + I(swb^2) +
                          temp + I(temp^2) +
                          adcount + avgba), 
                     family=skew_normal()
                     , data=brms.dat
                     , chains=3, cores=3
)
grow.brms.sn.bc <- brm(bf(bcbai~ -1 + intercept + (1|plot) + dia + I(dia^2) +
                         swb + I(swb^2) +
                         temp + I(temp^2) +
                         adcount + avgba), 
                    family=skew_normal()
                    , data=brms.dat
                    , chains=3, cores=3
)
# with increasing variance
grow.brms.snv <- brm(bf(baicm~ -1 + intercept + (1|plot) + dia + I(dia^2) +
                         swb + I(swb^2) +
                         temp + I(temp^2) +
                         adcount + avgba, sigma~dia), 
                    family=skew_normal()
                    , data=brms.dat
                    , chains=3, cores=3
)


## Lognormal error distribution----
## all data + constant
grow.brms.lnorm <- brm(bf(posbai~ -1 + intercept + (1|plot) + dia + I(dia^2) +
                           swb + I(swb^2) +
                           temp + I(temp^2) +
                           adcount + avgba), 
                      family=lognormal()
                      , data=brms.dat
                      , chains=3, cores=3
)

## Mixture model----
## can we model this as a mixture of two gaussian processes? Not sure how to interpret this anyhow...and seems complicated
# prior <- c(
#   prior(normal(0, 0.5), sigma1),
# )
# 
# grow.brms.mix <- brm(bf(baicm~ (1|plot) + dia + I(dia^2) +
#                                 swb + I(swb^2) +
#                                 temp + I(temp^2) +
#                                 adcount + avgba), prior=prior,
#                           family=mix, data=brms.dat, chains=3, cores=3)

## Also tried: square root transformations, student t, hurdle models, and gams none seemed promising
l1 <- loo(grow.brms.norm)
l2 <- loo(grow.brms.norm.bc)
l3 <- loo(grow.brms.normv)
l4 <- loo(grow.brms.sn)
l5 <- loo(grow.brms.snv)
l6 <- loo(grow.brms.sn.bc)

compare(l2, l5, l6) 
loo_model_weights(list(l6, l7), method = "pseudobma")

## Change object for plotting----
grow.brms <- grow.brms.snv # *
grow.brms <- grow.brms.normv # *pref
grow.brms <- grow.brms.norm.bc
grow.brms <- grow.brms.norm   # Jeff: interesting that ordinary and pearson residuals are close to identical in normal model while show bid difference in normv


## posterior predictive distributions
pp_check(grow.brms) #+ xlim(-20, 150)
pp_check(grow.brms, type="scatter_avg")

## marginal effects
marginal_effects(grow.brms)

## residuals vs fitted value
fit <- fitted(grow.brms, summary=TRUE)

res <- residuals(grow.brms, summary = TRUE, method="fitted",type = "pearson")
plot(fit[,1], res[,1], xlab="fitted value", ylab="pearson residual", main=spptable$gensp[spptable$spid==i]) 
abline(0,0, col="red")

res.ord <- residuals(grow.brms, summary = TRUE, method="fitted",type = "ordinary")
plot(fit[,1], res.ord[,1], xlab="fitted value", ylab="ordinary residual", main=spptable$gensp[spptable$spid==i]) 
abline(0,0, col="red")

## distribution of residuals
par(mfrow=c(1,3))
hist(res[,1],main='pearson', 40)
hist(res.ord[,1],main='ordinary', 70)
plot(res[,1],res.ord[,1],xlab='pearson',ylab='ordinary')

## residuals vs raw data
plot(brms.dat$baicm, res[,1])
abline(0,0, col='red')

## residuals vs predictors
plot(brms.dat$dia, res[,1], xlab="diameter", ylab="residual", main=spptable$gensp[spptable$spid==i])
abline(0,0, col="red")

plot(brms.dat$swb, res[,1], xlab="swb", ylab="residual", main=spptable$gensp[spptable$spid==i])
abline(0,0, col="red")

plot(brms.dat$temp, res[,1], xlab="temperature", ylab="residual", main=spptable$gensp[spptable$spid==i])
abline(0,0, col="red")

plot(brms.dat$adcount, res[,1], xlab="stem density", ylab="residual", main=spptable$gensp[spptable$spid==i])
abline(0,0, col="red")

plot(brms.dat$avgba, res[,1], xlab="average stem size", ylab="residual", main=spptable$gensp[spptable$spid==i])
abline(0,0, col="red")


#### All species pdfs----
## normal model----
pdf("Output pdfs/growth_residuals.pdf")
for(i in 1:N_spp) {
  trees <- trees1[trees1$spid==i,]
  
  ## scale covariates 
  cov.mat.raw <- cbind(trees$dia, trees$ptcoeff, trees$avgba, trees$adcount, trees$temp)## data for Stan , trees$dist_nowCWD, trees$dist_histCWD
  cov.mat <- apply(cov.mat.raw, 2, scale)
  cov.means <- colMeans(cov.mat.raw, na.rm=T)
  cov.sds <- apply(cov.mat.raw, 2, sd, na.rm=T)
  
  # data for brms
  brms.dat <- dplyr::select(trees, c(gensp, baicm, plot.uniqueID)) %>%
    bind_cols(as.data.frame(cov.mat))
  names(brms.dat)[4:8] <- c("dia","swb", "avgba", "adcount", "temp")
  brms.dat$plot <- as.numeric(droplevels(as.factor(brms.dat$plot.uniqueID)))
  
  grow.brms1 <- readRDS(paste("Ch 1/saved models/grow models/grow.fit",i,"b.rds", sep=""))
  
  res <- residuals(grow.brms1, summary = TRUE)
  fit <- fitted(grow.brms1, summary=TRUE)
  plot(fit[,1], res[,1], xlab="fitted value", ylab="residual", main=spptable$gensp[spptable$spid==i]) 
  abline(0,0, col="red")
  
  pp_check(grow.brms1)

}

dev.off()

## sigma~dia model----
pdf("Output pdfs/sigmadia_growth_residuals.pdf")
for(i in 1:N_spp) {
  trees <- trees1[trees1$spid==i,]
  
  ## scale covariates 
  cov.mat.raw <- cbind(trees$dia, trees$ptcoeff, trees$avgba, trees$adcount, trees$temp)## data for Stan , trees$dist_nowCWD, trees$dist_histCWD
  cov.mat <- apply(cov.mat.raw, 2, scale)
  cov.means <- colMeans(cov.mat.raw, na.rm=T)
  cov.sds <- apply(cov.mat.raw, 2, sd, na.rm=T)
  
  # data for brms
  brms.dat <- dplyr::select(trees, c(gensp, baicm, plot.uniqueID)) %>%
    bind_cols(as.data.frame(cov.mat))
  names(brms.dat)[4:8] <- c("dia","swb", "avgba", "adcount", "temp")
  brms.dat$plot <- as.numeric(droplevels(as.factor(brms.dat$plot.uniqueID)))
  
  grow.brms.normv <- brm(bf(baicm~ -1 + intercept + (1|plot) + dia + I(dia^2) +
                              swb + I(swb^2) +
                              temp + I(temp^2) +
                              adcount + avgba, sigma~dia), 
                         family=gaussian()
                         , data=brms.dat
                         , chains=3, cores=3
  )
  
  res <- residuals(grow.brms.normv, summary = TRUE)
  fit <- fitted(grow.brms.normv, summary=TRUE)
  plot(fit[,1], res[,1], xlab="fitted value", ylab="residual", main=spptable$gensp[spptable$spid==i]) 
  abline(0,0, col="red")
  
  pp_check(grow.brms.normv)
  
}

dev.off()

## log transform----
pdf("Output pdfs/growth_residuals_logmod.pdf")
for(i in 1:N_spp) {
  trees <- trees1[trees1$spid==i,]
  
  ## scale covariates 
  cov.mat.raw <- cbind(trees$dia, trees$ptcoeff, trees$avgba, trees$adcount, trees$temp)## data for Stan , trees$dist_nowCWD, trees$dist_histCWD
  cov.mat <- apply(cov.mat.raw, 2, scale)
  cov.means <- colMeans(cov.mat.raw, na.rm=T)
  cov.sds <- apply(cov.mat.raw, 2, sd, na.rm=T)
  
  # data for brms
  brms.dat <- dplyr::select(trees, c(gensp, baicm, plot.uniqueID)) %>%
    bind_cols(as.data.frame(cov.mat))
  names(brms.dat)[4:8] <- c("dia","swb", "avgba", "adcount", "temp")
  brms.dat$plot <- as.numeric(droplevels(as.factor(brms.dat$plot.uniqueID)))
  
  constant <- abs(min(brms.dat$baicm)) +1
  brms.dat$logbai <- log(brms.dat$baicm + constant)
  
  grow.brms2 <- brm(logbai~ -1 + intercept + (1|plot) + dia + I(dia^2) +
                      swb + I(swb^2) +
                      temp + I(temp^2) +
                      adcount + avgba, 
                    family=gaussian(), 
                    data=brms.dat, chains=3, cores=3, prior = c(prior_(~normal(0,1), class=~sd),
                                                                prior_(~normal(0,10), class=~b),
                                                                prior_(~normal(0,10), class=~sigma)))
  
  saveRDS(grow.brms2, paste("Ch 1/saved models/grow models/grow.fit.log",i,"b.rds", sep=""))
  
  res <- residuals(grow.brms2, summary = TRUE)
  fit <- fitted(grow.brms2, summary=TRUE)
  plot(fit[,1], res[,1], xlab="fitted value", ylab="residual", main=spptable$gensp[spptable$spid==i]) 
  abline(0,0, col="red")
  
}

dev.off()

## subset and log transform----
pdf("Output pdfs/growth_residuals_subset_logmod.pdf")
for(i in 1:N_spp) {
  trees <- trees1[trees1$spid==i,]
  
  ## scale covariates 
  cov.mat.raw <- cbind(trees$dia, trees$ptcoeff, trees$avgba, trees$adcount, trees$temp)## data for Stan , trees$dist_nowCWD, trees$dist_histCWD
  cov.mat <- apply(cov.mat.raw, 2, scale)
  cov.means <- colMeans(cov.mat.raw, na.rm=T)
  cov.sds <- apply(cov.mat.raw, 2, sd, na.rm=T)
  
  # data for brms
  brms.dat <- dplyr::select(trees, c(gensp, baicm, plot.uniqueID)) %>%
    bind_cols(as.data.frame(cov.mat)) %>%
    subset(baicm>0)
  names(brms.dat)[4:8] <- c("dia","swb", "avgba", "adcount", "temp")
  
  brms.dat$plot <- as.numeric(droplevels(as.factor(brms.dat$plot.uniqueID)))
  
  brms.dat$logbai <- log(brms.dat$baicm)
  
  grow.brms3 <- brm(bf(logbai~ -1 + intercept + (1|plot) + dia + I(dia^2) +
                         swb + I(swb^2) +
                         temp + I(temp^2) +
                         adcount + avgba), 
                       family=gaussian(), 
                       data=brms.dat, chains=3, cores=3, 
                    prior = c(prior_(~normal(0,1), class=~sd),
                              prior_(~normal(0,10), class=~b)
                              , prior_(~normal(0,10), class=~sigma))
                    )
  
  saveRDS(grow.brms3, paste("Ch 1/saved models/grow models/grow.fit.subset.log",i,"b.rds", sep=""))
  
  res <- residuals(grow.brms3, summary = TRUE)
  fit <- fitted(grow.brms3, summary=TRUE)
  plot(fit[,1], res[,1], xlab="fitted value", ylab="residual", main=spptable$gensp[spptable$spid==i]) 
  abline(0,0, col="red")

  
}

dev.off()

## simulation to check residuals----
x1 <- rnorm(1000,11,.5)
x2<-rnorm(1000,5,95)

b1 <- 0.1
b2 <- 0.9
sigma <- 0.4

eps <- rnorm(x1,1000,sigma)
y <- b1*x1  + b2*x2  +eps
Y<-as.matrix(y)
X<-cbind(as.matrix(x1),as.matrix(x2))

dat <- as.data.frame(cbind(Y, X))
names(dat) <- c("y", "x1", "x2")

mod <- brm(y~x1 +x2, data=dat, family=gaussian(), chains = 3, cores=3)

pp_check(mod, type="scatter_avg", nsamples=100)
res <- residuals(mod, summary = TRUE)
fit <- predict(fitted, summary=TRUE)
plot(fit[,1], res[,1]) 
abline(0,0, col="red")


