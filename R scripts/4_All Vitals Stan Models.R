## Author: Teresa Bohner
## Purpose: Chapter 1 analyses: Occurrence, abundance (density), and vital rate modeling
## Date Modified: 5-8-2019
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

data.path <- "Processed Data/Analysis Ready Data/"
mod.path <- "Stan Models/"
save.path <- "Stan Models/Saved Models/"

#################### Occurence data manipulation----
rm(list = ls())
occurrence <- read.csv(paste0(data.path, "occurrence.new.csv"))
spptable <- read.csv(paste0(data.path, 'spptable.new.csv'))
spptable <- spptable[,c(1:2)]
names(spptable)[2] <- "gensp"
spptable$gensp <- as.character(spptable$gensp)
# abundance <- rename(abundance, "cwd"="ptcoeff")
# occurrence <- rename(occurrence, "cwd"="ptcoeff")

length(unique(occurrence$plot.uniqueID))

tableinfo_occ <- dplyr::group_by(occurrence, gensp) %>%
  dplyr::summarize(nplots_p = sum(pres, na.rm=T))

write.csv(tableinfo_occ, paste0(data.path, 'tableinfo.csv'), row.names = F)

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
  
  min_swb <- (maxmin$min_swb-swb.mean)/swb.sd
  max_swb <- (maxmin$max_swb-swb.mean)/swb.sd
  
  N <- dim(occurrence)[1]
  N_pred <- 200
  N_plots <- length(unique(occurrence$plot.uniqueID))
  plot <- as.vector(as.integer(factor(occurrence$plot.uniqueID)))
  pres <- as.vector(occurrence$pres)
  swb <- swb
  min_swb <- min_swb
  max_swb <- max_swb
  swb_pred <- seq(min(swb),max(swb), length.out=N_pred)

  pa_data <- list(N=N, N_plots=N_plots, plot=plot, N_pred=N_pred,
                  pres=pres, swb=swb, swb_pred=swb_pred,
                  min_swb=min_swb, max_swb=max_swb)

  pars=c("a","b_swb","b2_swb",
         "opt_swb","thetamin_swb","thetamax_swb","thetaopt_swb",
         "theta_pred") # "y_pred","log_lik"

  ## Occurence models----
  rstan_options(auto_write = TRUE) 
  pa.compiled <- stan_model(file=paste0(mod.path, "binomial_swb_spp.stan"))
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
  # # extract pareto-k
  # log_lik1 <- extract_log_lik(pa.fit, merge_chains = FALSE)
  # rel_n_eff <- relative_eff(exp(log_lik1))
  # loo1 <- loo(log_lik1, r_eff = rel_n_eff, cores = 2, save_psis = TRUE)
  # # plot(loo1, label_points = TRUE)
  # pareto_k_all$N_0.7[i] <- length(which(loo1$diagnostics$pareto_k>0.7))
  # pareto_k_all$N_1[i] <- length(which(loo1$diagnostics$pareto_k>1))
  # pareto_k_all$N_tot[i] <- length(loo1$diagnostics$pareto_k)
  
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
  saveRDS(pa.fit, paste(save.path, "pa models/pa.fit",i,".rds", sep=""))

  
}

maxminopt_p <- maxminopt_all_med
theta <- theta_all
occurrence <- occurrence1

## Occurrence figures ----
trans_gray <- adjustcolor(col = "gray",alpha.f = 0.4)

pdf("Output pdfs/pres by species.pdf")
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
  # points(swb[spp==i], jitter(pres[spp==i], 0.1))
  mtext("wet", side=1, line=2, adj=1)
  mtext("dry", side=1, line=2, adj=0)
  
  polygon(c(xdum,rev(xdum)),c(inv.logit(theta$`97.5%`[theta$spid==i]),
                              rev(inv.logit(theta$`2.5%`[theta$spid==i]))),col=trans_gray, border=NA)
  lines(xdum, inv.logit(theta$`50%`[theta$spid==i]), col=4)
  # lines(xdum, inv.logit(fixef(pa.brms)[1] + fixef(pa.brms)[2]*swb_pred + fixef(pa.brms)[3]*swb_pred^2))
  # 
  # plot(inv.logit(fixef(pa.brms)[1] + fixef(pa.brms)[2]*swb_pred + fixef(pa.brms)[3]*swb_pred^2),
  #      inv.logit(theta$`50%`[theta$spid==i]))

  abline(v=(maxminopt_all_mean$opt_swb_p[maxminopt_all_mean$spid==i]-xmean)/xsd, col=2, lty=2)
  abline(v=maxminopt_all_med$opt_swb_p[maxminopt_all_med$spid==i], col=4, lty=2)
  abline(v=maxminopt_all_med$opt_swb_p_25[maxminopt_all_med$spid==i], col="gray", lty=2)
  abline(v=maxminopt_all_med$opt_swb_p_75[maxminopt_all_med$spid==i], col="gray", lty=2)
  

}

dev.off()

## Occurrence data export----
# merge abundance and occurrence results
maxminopt <- maxminopt_p
maxminopt <- merge(maxminopt, spptable, by='spid')
maxminopt <- bind_cols(maxminopt, as.data.frame(r2_all))

swb_pred_p <- matrix(data=NA, ncol = N_spp, nrow=200)
for(i in 1:N_spp){
  swb_pred_p[1:200,i] <- seq(min(occurrence$ptcoeff[occurrence$spid==i],na.rm=T),max(occurrence$ptcoeff[occurrence$spid==i],na.rm=T), length.out=N_pred)
}
swb_pred_p <- as.data.frame(swb_pred_p)
names(swb_pred_p) <- spptable$gensp
swb_axes_p <- melt(swb_pred_p, variable.name="gensp", value.name="swb_pred")


# save csv's for vital rate modeling
write.csv(dplyr::select(maxminopt, -spid), paste0(data.path, 'distromaxminopt_swb.csv'), row.names = F)
write.csv(dplyr::select(theta, -spid), paste0(data.path, 'theta_swb.csv'), row.names = F)
write.csv(swb_axes_p, paste0(data.path, 'swb_axes_p.csv'), row.names = F)

save.image(paste0(save.path, 'distribution modeling swb.RData'), compress=T)


## Abundance data manipulation----
## read in data
rm(list = ls())
abundance <- read.csv(paste0(data.path, "abundance.new.csv"))
spptable <- read.csv(paste0(data.path, 'spptable.new.csv'))
spptable <- spptable[,c(1:2)]
names(spptable)[2] <- "gensp"
spptable$gensp <- as.character(spptable$gensp)
maxminopt <- read.csv(paste0(data.path, 'distromaxminopt_swb.csv'))
maxminopt <- merge(maxminopt, spptable, by="gensp")
theta <- read.csv(paste0(data.path, 'theta_swb.csv'))
theta$seq <- as.numeric(rownames(theta))
theta <- merge(theta, spptable, by="gensp")
theta <- theta[order(theta$seq),]
swb_axes_p <- read.csv(paste0(data.path, 'swb_axes_p.csv'))
swb_axes_p <- merge(swb_axes_p, spptable, by="gensp")

abundance$avgba <- abundance$totba/abundance$totabun

maxmin <- dplyr::group_by(abundance, spid)%>%
  dplyr::summarise(min_swb = min(ptcoeff, na.rm=T),
                   max_swb = max(ptcoeff, na.rm=T))

N_spp <- length(unique(abundance$spid))

eta_all <- matrix(NA, nrow = 200*N_spp, ncol=6)
eta_all <- data.frame(eta_all)

med_SWBs_all <- data.frame(matrix(NA, nrow=N_spp, ncol=30))
r2_all <- rep(NA, length=N_spp)
# pareto_k_all <- data.frame(matrix(NA, nrow=N_spp, ncol=3))
# pareto_k_all <- cbind(spptable, pareto_k_all)
# names(pareto_k_all)[3:5] <- c("N_0.7", "N_1", "N_tot")

abund1 <- abundance
maxmin1 <- maxmin
maxminopt1 <- maxminopt


spids <- dplyr::group_by(abund1, gensp) %>%
  dplyr::summarise(spid = unique(spid))
options(tibble.print_max = Inf)
spids

#### Abundance loop----
for(i in 1:N_spp){
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
                    # temp=temp, avgba=avgba,
                    pres_swb_dat=pres_swb_dat, min_swb=min_swb, max_swb=max_swb)
  
  pars=c("a","b_swb","b2_swb","phi"
         ,"wet_swb","dry_swb","pres_swb","opt_swb", "par_swb"
         ,"n_hat_wet","n_hat_dry","n_hat_opt", "n_hat_pres"
         ,"eta_pred") # ,'y_rep' "b_temp","b2_temp","b_ba",
  
  ## Abundance models---- 
  rstan_options(auto_write = TRUE) 
  nb.compiled <- stan_model(file=paste0(mod.path, "abundance_new_spp_nb.stan"))
  
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
  
  saveRDS(nb.fit, paste(save.path, "abun models/abun.fit",i,".rds", sep=""))

  # ## Abundance diagnostics----
  # posterior <- as.array(nb.fit)
  # color_scheme_set("viridis")
  # pdf(paste("Output pdfs/abun regression coeff/abun regression coeff", i, ".pdf", sep=""))
  # # print(plot(loo1, label_points = TRUE))
  # print(mcmc_intervals(posterior, pars=c("a","b_swb","b_ba","b_temp","b2_swb","b2_temp", "phi")) +
  #         ggplot2::labs(title = spptable$gensp[i]))
  # print(mcmc_combo(posterior, pars=c("a","b_swb","b_ba","b_temp")))
  # print(mcmc_acf_bar(posterior, pars=c("a","b_swb","b_ba","b_temp")))
  # print(mcmc_combo(posterior, pars=c("b2_swb","b2_temp","phi","lp__")))
  # print(mcmc_acf_bar(posterior, pars=c("b2_swb","b2_temp","phi","lp__")))
  # print(mcmc_rhat_hist(rhat(nb.fit)[1:7]))
  # print(mcmc_neff_hist(neff_ratio(nb.fit)[1:7]))
  # dev.off()
  
  posterior <- as.array(nb.fit)
  color_scheme_set("viridis")
  pdf(paste("Output pdfs/abun regression coeff/abun regression coeff", i, ".pdf", sep=""))
  # print(plot(loo1, label_points = TRUE))
  print(mcmc_intervals(posterior, pars=c("a","b_swb","b2_swb","phi")) +
          ggplot2::labs(title = spptable$gensp[i]))
  print(mcmc_combo(posterior, pars=c("a","b_swb","b2_swb","phi")))
  print(mcmc_acf_bar(posterior, pars=c("a","b_swb","b2_swb","phi")))
  print(mcmc_rhat_hist(rhat(nb.fit)[1:4]))
  # print(mcmc_neff_hist(neff_ratio(nb.fit)[1:4]))
  dev.off()
  
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
abunswbs <- med_SWBs_all
trans_gray <- adjustcolor(col = "gray",alpha.f = 0.4)

pdf("Output pdfs/abun by species_simple.pdf")
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
  # points(swb[spp==i], jitter(pres[spp==i], 0.1))
  mtext("wet", side=1, line=2, adj=1)
  mtext("dry", side=1, line=2, adj=0)
  
  polygon(c(xdum,rev(xdum)),c(eta_all$`97.5%`[eta_all$spid==i],rev(eta_all$`2.5%`[eta_all$spid==i])),col=trans_gray, border=NA)
  lines(xdum, eta_all$`50%`[eta_all$spid==i], col=4)
  # abline(v=med_SWBs_all$wet_swb[med_SWBs_all$spid==i], lwd=2)
  # abline(v=med_SWBs_all$dry_swb[med_SWBs_all$spid==i], lwd=2)
  abline(v=med_SWBs_all$opt_swb[med_SWBs_all$spid==i], col=4, lty=2)
  
}
dev.off()

pdf("Output pdfs/jan19 abundance figures/abundance and occurrence median_simple.pdf")
for(i in 1:N_spp){
  xdum <- seq(min(abund1$ptcoeff[abund1$spid==i],na.rm=T),max(abund1$ptcoeff[abund1$spid==i],na.rm=T), length.out=N_pred)
  
  ydum <- seq(0,(max(eta_all$`50%`[eta_all$spid==i])+10), length.out=length(xdum))
  
  x2dum <- sort(swb_axes_p$swb_pred[swb_axes_p$spid==i])
  y2dum <- seq(0,1, length.out=length(xdum))
  lowlim <- ifelse(med_SWBs_all$pres_swb[med_SWBs_all$spid==i]<med_SWBs_all$dry_swb[med_SWBs_all$spid==i], 
                   med_SWBs_all$pres_swb[med_SWBs_all$spid==i],
                   med_SWBs_all$dry_swb[med_SWBs_all$spid==i])
  hilim <- ifelse(med_SWBs_all$pres_swb[med_SWBs_all$spid==i]>med_SWBs_all$wet_swb[med_SWBs_all$spid==i], 
                  med_SWBs_all$pres_swb[med_SWBs_all$spid==i],
                  med_SWBs_all$wet_swb[med_SWBs_all$spid==i])
  
  par(mar = c(5,5,2,5))
  plot(x2dum, y2dum,type="n", axes=F, xlab=NA, ylab=NA
       , xlim=c(lowlim, hilim)) #
  axis(side = 4)
  mtext(side = 4, line = 3, 'Occurrence')
  # occurrence curve
  polygon(c(x2dum,rev(x2dum)),c(inv.logit(theta$X2.5.[theta$spid==i]),rev(inv.logit(theta$X97.5.[theta$spid==i]))),col=trans_gray, border=NA)
  lines(x2dum, inv.logit(theta$X50.[theta$spid==i]), col=2)
  
  par(new=T)  
  plot(xdum, ydum, type="n", main=spptable[i,2]
       , xlim=c(lowlim, hilim)
       , xlab="swb gradient", ylab="Abundance") #
  # growth curve
  polygon(c(xdum,rev(xdum)),c(eta_all$`97.5%`[eta_all$spid==i],rev(eta_all$`2.5%`[eta_all$spid==i])),col=trans_gray, border=NA)
  lines(xdum, eta_all$`50%`[eta_all$spid==i], col=1)
  
  abline(v=med_SWBs_all$opt_swb[med_SWBs_all$spid==i], col=1, lty=0)
  segments(x0=med_SWBs_all$opt_swb[med_SWBs_all$spid==i], y0=-10,
           x1=med_SWBs_all$opt_swb[med_SWBs_all$spid==i],y1=med_SWBs_all$n_hat_opt[med_SWBs_all$spid==i], col=1, lty=2)
  abline(h=med_SWBs_all$g_hat_opt[med_SWBs_all$spid==i], col=1, lty=0)
  segments(x0=0, y0=med_SWBs_all$n_hat_opt[med_SWBs_all$spid==i],
           x1=med_SWBs_all$opt_swb[med_SWBs_all$spid==i],y1=med_SWBs_all$n_hat_opt[med_SWBs_all$spid==i], col=1, lty=2)
  abline(h=med_SWBs_all$n_hat_pres[med_SWBs_all$spid==i], col=2, lty=0)
  segments(x0=0, y0=med_SWBs_all$n_hat_pres[med_SWBs_all$spid==i],
           x1=med_SWBs_all$pres_swb[med_SWBs_all$spid==i],y1=med_SWBs_all$n_hat_pres[med_SWBs_all$spid==i], col=2, lty=3)
  
  abline(v=med_SWBs_all$pres_swb[med_SWBs_all$spid==i], col=2, lty=0)
  segments(x0=med_SWBs_all$pres_swb[med_SWBs_all$spid==i], y0=-10,
           x1=med_SWBs_all$pres_swb[med_SWBs_all$spid==i],y1=med_SWBs_all$n_hat_pres[med_SWBs_all$spid==i], col=2, lty=3)
  
  legend("topright",
         legend=c("abundance","occurrence"),
         lty=c(1,1), col=c("black","red"), cex=0.85)
  
}
dev.off()

## Abundance data export----
# pareto_k_all$k_0.7_percent <- pareto_k_all$N_0.7/pareto_k_all$N_tot*100
# pareto_k_all <- pareto_k_all[order(-pareto_k_all$k_0.7_percent),]
# write.csv(pareto_k_all, 'grow pareto k.csv', row.names=F)
med_SWBs_all <- bind_cols(med_SWBs_all, as.data.frame(r2_all))
write.csv(betas_all, paste0(data.path, "abun_betas_swb_temp.csv"), row.names=F)
write.csv(med_SWBs_all, paste0(data.path, "abunswbs_temp.csv"), row.names = F)
save.image(paste0(save.path, "abundance loop_temp.RData"))



## Growth data manipulation----
# *Clear workspace to start 
rm(list = ls())
# read in data
spptable <- na.omit(read.csv(paste0(data.path, 'spptable.new.csv'))[,c(1,4)])
names(spptable)[2] <- "gensp"
spptable$gensp <- as.character(spptable$gensp)
maxminopt <- read.csv(paste0(data.path, 'distromaxminopt_swb.csv'))
maxminopt <- merge(maxminopt, spptable, by="gensp")
theta <- read.csv(paste0(data.path, 'theta_swb.csv'))
theta$seq <- as.numeric(rownames(theta))
theta <- merge(theta, spptable, by="gensp")
theta <- theta[order(theta$seq),]
swb_axes_p <- read.csv(paste0(data.path, 'swb_axes_p.csv'))
swb_axes_p <- merge(swb_axes_p, spptable, by="gensp")
trees <- read.csv(paste0(data.path, 'trees.new.csv'))

# trees <- rename(trees, "cwd"="ptcoeff")
tableinfo <- read.csv(paste0(data.path, 'tableinfo.csv'))

# trees <-trees[trees$bai> -0.03&trees$bai< 0.03,]

hist(trees$bai)
trees$baicm <- trees$bai*10000
trees$baratio <- trees$ba/trees$avgba

trees2 <- group_by(trees, gensp) %>%
  mutate(mean = mean(baicm, na.rm=T),
            sd = sd(baicm, na.rm=T),
            upper = sd(baicm, na.rm=T)*3 + mean(baicm, na.rm=T),
            lower = mean(baicm, na.rm=T) - sd(baicm, na.rm=T)*3,
            chiscore = scores(baicm, type="chisq", prob=0.999),
            chiscore2 = scores(ptcoeff, type="chisq", prob=0.999))

trees2$outlier <- ifelse(trees2$baicm<trees2$upper&trees2$baicm>trees2$lower, "no", "yes")

length(which(trees2$chiscore=="TRUE"))#/nrow(trees2)
length(which(trees2$chiscore2=="TRUE"))#/nrow(trees2)
length(which(trees2$chiscore=="TRUE"&trees2$chiscore2=="TRUE"))
length(which(trees2$outlier=="yes"))#/nrow(trees2)

trees <- subset(trees2, chiscore!="TRUE")

tableinfo_grow <- dplyr::group_by(trees, gensp) %>%
  dplyr::summarize(ntrees_g = length(unique(tree.uniqueID)),
                   nplots_g = length(unique(plot.uniqueID)))
tableinfo <- merge(tableinfo, tableinfo_grow, by="gensp", all=T)
write.csv(tableinfo, paste0(data.path, 'tableinfo.csv'), row.names = F)

maxgrow <- group_by(trees, gensp) %>%
  summarize(max = max(baicm, na.rm=T),
            mean = mean(baicm, na.rm=T))
write.table(maxgrow, paste0(data.path, "max_grow_obs.csv"), row.names = F)

maxmin <- dplyr::group_by(trees, spid)%>%
  dplyr::summarise(min_swb = min(ptcoeff, na.rm=T),
            max_swb = max(ptcoeff, na.rm=T))
N_spp <- length(unique(trees$spid))

mu_all <- matrix(NA, nrow = 200*N_spp, ncol=6)
mu_all <- data.frame(mu_all)

med_SWBs_all <- data.frame(matrix(NA, nrow=N_spp, ncol=26))
r2_all <- rep(NA, length=N_spp)

# pareto_k_all <- data.frame(matrix(NA, nrow=N_spp, ncol=3))
# pareto_k_all <- cbind(spptable, pareto_k_all)
# names(pareto_k_all)[3:5] <- c("N_0.7", "N_1", "N_tot")

trees1 <- trees
maxmin1 <- maxmin
maxminopt1 <- maxminopt

spids <- dplyr::group_by(trees1, gensp) %>%
  dplyr::summarise(spid = unique(spid))
options(tibble.print_max = Inf)
spids

pdf("Output pdfs/grow_hist.pdf")
for(i in 1:N_spp){
  trees <- trees1[trees1$spid==i,]
  par(mfrow=c(1,2))
  hist(trees$bai, main = spids$gensp[i])
  hist(log(trees$bai), main = paste("log(", spids$gensp[i], ")", collapse=""))
}
dev.off()
## Growth loop----
for(i in 3:N_spp){
  trees <- trees1[trees1$spid==i,]
  maxmin <- maxmin1[maxmin1$spid==i,]
  maxminopt <- maxminopt1[maxminopt1$spid==i,]
  
  ## scale covariates 
  cov.mat.raw <- cbind(trees$dia, trees$ptcoeff, trees$avgba, trees$adcount, trees$temp)## data for Stan , trees$dist_nowCWD, trees$dist_histCWD
  cov.mat <- apply(cov.mat.raw, 2, scale)
  cov.means <- colMeans(cov.mat.raw, na.rm=T)
  cov.sds <- apply(cov.mat.raw, 2, sd, na.rm=T)
  
  
  min_swb <- (maxmin$min_swb-cov.means[2])/cov.sds[2]
  max_swb <- (maxmin$max_swb-cov.means[2])/cov.sds[2]
  
  ## put opt on same swb scale
  opt_swb_p <- (maxminopt$opt_swb_p-cov.means[2])/cov.sds[2]
  
  # data for stan
  N <- dim(trees)[1]
  N_pred <- 200
  spp <- trees$spid
  y <- as.vector(trees$baicm)
  dia <- cov.mat[,1]
  swb <- cov.mat[,2]
  temp <- cov.mat[,5]
  swb_pred <- seq(min(cov.mat[,2],na.rm=T),max(cov.mat[,2],na.rm=T), length.out=N_pred)
  min_swb <- min_swb
  max_swb <- max_swb
  pres_swb_dat <- opt_swb_p
  avgba <- cov.mat[,3]
  adcount <- cov.mat[,4]
  N_plts <- length(unique(trees$plot.uniqueID))
  plot <- as.numeric(droplevels(as.factor(trees$plot.uniqueID)))
  
  grow_data <- list(N=N, min_swb=min_swb, max_swb=max_swb, pres_swb_dat=pres_swb_dat,
                    y=y, dia=dia, swb=swb, temp, swb_pred=swb_pred, 
                    avgba=avgba, adcount=adcount, N_plts=N_plts, plot=plot) 
  
  pars=c("a_sp","b_dia","b_swb","b_ba","b_adcount","b_temp","b2_dia","b2_swb","b2_temp"
         ,"sigma_dia","sigma_c","wet_swb","dry_swb","pres_swb","opt_swb", "par_swb"
         ,"g_hat_wet","g_hat_dry","g_hat_opt", "g_hat_pres"
         , "mu_pred") #"log_lik"
  
  ## Growth models----
  rstan_options(auto_write = TRUE) 
  grow.compiled <- stan_model(file=paste0(mod.path, "grow_ind_swb_temp_spp_var.stan"))
  
  grow.fit<- sampling(
    object = grow.compiled
    , data=grow_data
    , chains=3
    , iter=1000
    , warmup=500
    , thin=1
    , refresh=100
    , pars=pars
    , cores=3
    , save_warmup=FALSE
    , control= list(max_treedepth=11, adapt_delta=0.99)
  )
  
  
  saveRDS(grow.fit, paste(save.path, "grow models/grow.fit",i,".rds", sep=""))
  
  ## Growth diagnostics----
  # # extract pareto-k
  # log_lik1 <- extract_log_lik(grow.fit, merge_chains = FALSE)
  # rel_n_eff <- relative_eff(exp(log_lik1))
  # loo1 <- loo(log_lik1, r_eff = rel_n_eff, cores = 2, save_psis = TRUE)
  # # plot(loo1, label_points = TRUE)
  # pareto_k_all$N_0.7[i] <- length(which(loo1$diagnostics$pareto_k>0.7))
  # pareto_k_all$N_1[i] <- length(which(loo1$diagnostics$pareto_k>1))
  # pareto_k_all$N_tot[i] <- length(loo1$diagnostics$pareto_k)

  posterior <- as.array(grow.fit)
  color_scheme_set("viridis")
  pdf(paste("Output pdfs/growth regression coeff/grow regression coeff", i, ".pdf", sep=""))
  # print(plot(loo1, label_points = TRUE))
  print(mcmc_intervals(posterior, pars=c("a_sp","b_dia","b_swb","b_ba","b_adcount","b_temp", "b2_dia","b2_swb","b2_temp","sigma_dia","sigma_c")) +
          ggplot2::labs(title = spptable$gensp[i]))
  print(mcmc_combo(posterior, pars=c("a_sp","b_dia","b_swb","b_ba","b_adcount","b_temp")))
  print(mcmc_acf_bar(posterior, pars=c("a_sp","b_dia","b_swb","b_ba","b_adcount","b_temp")))
  print(mcmc_combo(posterior, pars=c("b2_dia","b2_swb","b2_temp","sigma_dia","sigma_c","lp__")))
  print(mcmc_acf_bar(posterior, pars=c("b2_dia","b2_swb","b2_temp","sigma_dia","sigma_c","lp__")))
  print(mcmc_rhat_hist(rhat(grow.fit)[2:12]))
  # print(mcmc_neff_hist(neff_ratio(grow.fit)[2:12]))
  dev.off()
  
  ## Growth output----
  # predicted curve
  mu <- as.data.frame(summary(grow.fit, pars="mu_pred")$summary[, c(1,4,6,8)])
  spid <- rep(i, N_pred)
  mu$spid <- spid
  mu <- merge(mu, spptable, by='spid')
  
  if(i==1) mu_all <- mu
  if(i>1) mu_all <- rbind(mu_all, mu)
  # mu_all[c((i*200-200+1):(i*200)),] <- mu
  
  ## betas
  betas <- as.data.frame(summary(grow.fit, pars=c("b_dia","b_swb","b_ba","b_adcount","b_temp","b2_dia","b2_swb","b2_temp"),probs = c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))$summary[, c(1:10, 12)])
  betas$param <- rownames(betas)
  
  betas$spid <- i
  betas <- merge(betas, spptable, by='spid')
  
  if(i==1) betas_all <- betas
  if(i>1) betas_all <- rbind(betas_all, betas)
  betas_all[c((i*8-8+1):(i*8)),] <- betas
  
  ## SWB and subsequent g-hat values
  grow_swb <- as.data.frame(summary(grow.fit, pars="opt_swb")$summary[, c(1,4,5,6,7,8)])
  grow_swb_unscaled <- cov.means[2] + grow_swb*cov.sds[2]
  par_swb <- as.data.frame(summary(grow.fit, pars="par_swb")$summary[,1])
  par_swb_unscaled <- cov.means[2] + par_swb*cov.sds[2]
  pres_swb <- as.data.frame(summary(grow.fit, pars="pres_swb")$summary[,1])
  pres_swb_unscaled <- cov.means[2] + pres_swb*cov.sds[2]
  wet_swb <- as.data.frame(summary(grow.fit, pars="wet_swb")$summary[,1])
  wet_swb_unscaled <- cov.means[2] + wet_swb*cov.sds[2]
  dry_swb <- as.data.frame(summary(grow.fit, pars="dry_swb")$summary[,1])
  dry_swb_unscaled <- cov.means[2] + dry_swb*cov.sds[2]
  
  g_hat_grow <- as.data.frame(summary(grow.fit, pars="g_hat_opt")$summary[, c(1,4,6,8)])
  g_hat_pres <- as.data.frame(summary(grow.fit, pars="g_hat_pres")$summary[, c(1,4,6,8)])
  g_hat_wet <- as.data.frame(summary(grow.fit, pars="g_hat_wet")$summary[, c(1,4,6,8)])
  g_hat_dry <- as.data.frame(summary(grow.fit, pars="g_hat_dry")$summary[, c(1,4,6,8)])
  
  # Medians
  SWBs <- as.data.frame(cbind(as.list(grow_swb_unscaled[2:6,1]), pres_swb_unscaled, wet_swb_unscaled, dry_swb_unscaled, 
                              g_hat_grow[3,1], g_hat_pres[3,1], g_hat_wet[3,1], g_hat_dry[3,1]))
  names(SWBs) <- c('grow_swb_2.5','grow_swb_25','grow_swb','grow_swb_75','grow_swb_97.5','pres_swb','wet_swb','dry_swb','g_hat_grow', 'g_hat_pres','g_hat_wet','g_hat_dry')
  SWBs$spid <- i
  
  # Choose greater: parabola or endpoint optimum
  SWBs$grow_swb <- ifelse(SWBs$wet_swb>=SWBs$grow_swb&SWBs$grow_swb>=SWBs$dry_swb, SWBs$grow_swb,NA)
  SWBs$g_hat_grow <- ifelse(!is.na(SWBs$grow_swb), SWBs$g_hat_grow,NA)
  SWBs$opt_swb <- ifelse(!is.na(SWBs$g_hat_grow),ifelse(SWBs$g_hat_grow>SWBs$g_hat_wet & SWBs$g_hat_grow>SWBs$g_hat_dry
                                                        , SWBs$grow_swb, 
                                                        ifelse(SWBs$g_hat_dry>SWBs$g_hat_wet 
                                                               , SWBs$dry_swb, SWBs$wet_swb)), 
                         ifelse(is.na(SWBs$g_hat_grow),ifelse(SWBs$g_hat_dry>SWBs$g_hat_wet 
                                                              , SWBs$dryswb, SWBs$wetswb)))
  SWBs$opt_swb_2.5 <- ifelse(!is.na(SWBs$g_hat_grow),ifelse(SWBs$g_hat_grow>SWBs$g_hat_wet & SWBs$g_hat_grow>SWBs$g_hat_dry
                                                            , SWBs$grow_swb_2.5, 
                                                            ifelse(SWBs$g_hat_dry>SWBs$g_hat_wet 
                                                                   , SWBs$dry_swb, SWBs$wet_swb)), 
                             ifelse(is.na(SWBs$g_hat_grow),ifelse(SWBs$g_hat_dry>SWBs$g_hat_wet 
                                                                  , SWBs$dryswb, SWBs$wetswb)))
  SWBs$opt_swb_25 <- ifelse(!is.na(SWBs$g_hat_grow),ifelse(SWBs$g_hat_grow>SWBs$g_hat_wet & SWBs$g_hat_grow>SWBs$g_hat_dry
                                                           , SWBs$grow_swb_25, 
                                                           ifelse(SWBs$g_hat_dry>SWBs$g_hat_wet 
                                                                  , SWBs$dry_swb, SWBs$wet_swb)), 
                            ifelse(is.na(SWBs$g_hat_grow),ifelse(SWBs$g_hat_dry>SWBs$g_hat_wet 
                                                                 , SWBs$dryswb, SWBs$wetswb)))
  SWBs$opt_swb_75 <- ifelse(!is.na(SWBs$g_hat_grow),ifelse(SWBs$g_hat_grow>SWBs$g_hat_wet & SWBs$g_hat_grow>SWBs$g_hat_dry
                                                           , SWBs$grow_swb_75, 
                                                           ifelse(SWBs$g_hat_dry>SWBs$g_hat_wet 
                                                                  , SWBs$dry_swb, SWBs$wet_swb)), 
                            ifelse(is.na(SWBs$g_hat_grow),ifelse(SWBs$g_hat_dry>SWBs$g_hat_wet 
                                                                 , SWBs$dryswb, SWBs$wetswb)))
  SWBs$opt_swb_97.5 <- ifelse(!is.na(SWBs$g_hat_grow),ifelse(SWBs$g_hat_grow>SWBs$g_hat_wet & SWBs$g_hat_grow>SWBs$g_hat_dry
                                                             , SWBs$grow_swb_97.5, 
                                                             ifelse(SWBs$g_hat_dry>SWBs$g_hat_wet 
                                                                    , SWBs$dry_swb, SWBs$wet_swb)), 
                              ifelse(is.na(SWBs$g_hat_grow),ifelse(SWBs$g_hat_dry>SWBs$g_hat_wet 
                                                                   , SWBs$dryswb, SWBs$wetswb)))
  
  SWBs$g_hat_opt <- ifelse(!is.na(SWBs$g_hat_grow),ifelse(SWBs$g_hat_grow>SWBs$g_hat_wet & SWBs$g_hat_grow>SWBs$g_hat_dry
                                                          , SWBs$g_hat_grow, 
                                                          ifelse(SWBs$g_hat_dry>SWBs$g_hat_wet 
                                                                 , SWBs$g_hat_dry, SWBs$g_hat_wet)), 
                           ifelse(is.na(SWBs$g_hat_grow),ifelse(SWBs$g_hat_dry>SWBs$g_hat_wet 
                                                                , SWBs$g_hat_dry, SWBs$g_hat_wet)))
  
  ## Choose minimum for sensitivity
  SWBs$min_swb <- ifelse(!is.na(SWBs$g_hat_grow),ifelse(SWBs$g_hat_grow<SWBs$g_hat_wet & SWBs$g_hat_grow<SWBs$g_hat_dry
                                                        , SWBs$grow_swb, 
                                                        ifelse(SWBs$g_hat_dry<SWBs$g_hat_wet 
                                                               , SWBs$dry_swb, SWBs$wet_swb)), 
                         ifelse(is.na(SWBs$g_hat_grow),ifelse(SWBs$g_hat_dry<SWBs$g_hat_wet 
                                                              , SWBs$dryswb, SWBs$wetswb)))
  
  SWBs$g_hat_min <- ifelse(!is.na(SWBs$g_hat_grow),ifelse(SWBs$g_hat_grow<SWBs$g_hat_wet & SWBs$g_hat_grow<SWBs$g_hat_dry
                                                          , SWBs$g_hat_grow, 
                                                          ifelse(SWBs$g_hat_dry<SWBs$g_hat_wet 
                                                                 , SWBs$g_hat_dry, SWBs$g_hat_wet)), 
                           ifelse(is.na(SWBs$g_hat_grow),ifelse(SWBs$g_hat_dry<SWBs$g_hat_wet 
                                                                , SWBs$g_hat_dry, SWBs$g_hat_wet)))
  # Calculate mismatch and cost
  SWBs$mismatch_p <- SWBs$opt_swb - SWBs$pres_swb
  
  SWBs$cost_p <- SWBs$g_hat_opt - SWBs$g_hat_pres
  SWBs$cost_p <- ifelse(SWBs$cost_p<0, 0, SWBs$cost_p)
  
  SWBs$range <- SWBs$g_hat_opt - SWBs$g_hat_min
  SWBs$sens <- SWBs$range/(SWBs$wet_swb - SWBs$dry_swb)
  
  SWBs <- merge(SWBs, spptable, by='spid')
  
  med_SWBs_all[i,] <- SWBs
  names(med_SWBs_all) <- names(SWBs)
  
  print(paste("spp[",i,"] complete"), sep="")
}

growswbs <- med_SWBs_all

## Growth figures----
trans_gray <- adjustcolor(col = "gray",alpha.f = 0.4)

pdf("Output pdfs/grow by species.pdf")
for(i in 1:N_spp){
  xdum <- seq(min(trees1$ptcoeff[trees1$spid==i],na.rm=T),max(trees1$ptcoeff[trees1$spid==i],na.rm=T), length.out=N_pred)
  xmean <- mean(trees1$ptcoeff[trees1$spid==i],na.rm=T)
  xsd <- sd(trees1$ptcoeff[trees1$spid==i],na.rm=T)
  ncuts<-20
  chopsize<-cut(trees1$ptcoeff[trees1$spid==i],ncuts)
  
  ###start here
  
  binmean.size<-as.numeric(sapply(split(trees1$ptcoeff[trees1$spid==i],chopsize),mean, na.rm=T))
  binmean.grow<-as.numeric(sapply(split(trees1$baicm[trees1$spid==i],chopsize),mean, na.rm=T))
  plot(trees1$ptcoeff[trees1$spid==i], trees1$baicm[trees1$spid==i], type="n", ylim=c(min(binmean.grow, na.rm=T)-5, max(binmean.grow, na.rm=T)+10),
       main=spptable[i,2],
       xlab="swb gradient", 
       ylab="BAI (cm^2)")
  points(binmean.size, binmean.grow)
  # points(swb[spp==i], jitter(pres[spp==i], 0.1))
  mtext("wet", side=1, line=2, adj=1)
  mtext("dry", side=1, line=2, adj=0)
  
  polygon(c(xdum,rev(xdum)),c(mu_all$`97.5%`[mu_all$spid==i],rev(mu_all$`2.5%`[mu_all$spid==i])),col=trans_gray, border=NA)
  lines(xdum, mu_all$`50%`[mu_all$spid==i], col=4)
  # abline(v=med_SWBs_all$wet_swb[med_SWBs_all$spid==i], lwd=2)
  # abline(v=med_SWBs_all$dry_swb[med_SWBs_all$spid==i], lwd=2)
  abline(v=med_SWBs_all$opt_swb[med_SWBs_all$spid==i], col=4, lty=2)
  
}
dev.off()

trans_gray <- adjustcolor(col = "gray",alpha.f = 0.4)


pdf("Output pdfs/growth and occurrence median.pdf")
for(i in 1:N_spp){
  xdum <- seq(min(trees1$ptcoeff[trees1$spid==i],na.rm=T),max(trees1$ptcoeff[trees1$spid==i],na.rm=T), length.out=N_pred)
  
  ydum <- seq((min(mu_all$`50%`[mu_all$spid==i])-10),(max(mu_all$`50%`[mu_all$spid==i])+10), length.out=length(xdum))
  
  x2dum <- sort(swb_axes_p$swb_pred[swb_axes_p$spid==i])
  # y2dum <- seq(0,(max(eta$mean[eta$spid==i])+10), length.out=length(xdum))
  y2dum <- seq(0,1, length.out=length(xdum))
  lowlim <- ifelse(med_SWBs_all$pres_swb[med_SWBs_all$spid==i]<med_SWBs_all$dry_swb[med_SWBs_all$spid==i], 
                   med_SWBs_all$pres_swb[med_SWBs_all$spid==i],
                   med_SWBs_all$dry_swb[med_SWBs_all$spid==i])
  hilim <- ifelse(med_SWBs_all$pres_swb[med_SWBs_all$spid==i]>med_SWBs_all$wet_swb[med_SWBs_all$spid==i], 
                  med_SWBs_all$pres_swb[med_SWBs_all$spid==i],
                  med_SWBs_all$wet_swb[med_SWBs_all$spid==i])
  
  par(mar = c(5,5,2,5))
  plot(x2dum, y2dum,type="n", axes=F, xlab=NA, ylab=NA
       , xlim=c(lowlim, hilim)) #
  axis(side = 4)
  mtext(side = 4, line = 3, 'Occurrence')
  # occurrence curve
  polygon(c(x2dum,rev(x2dum)),c(inv.logit(theta$X2.5.[theta$spid==i]),rev(inv.logit(theta$X97.5.[theta$spid==i]))),col=trans_gray, border=NA)
  lines(x2dum, inv.logit(theta$X50.[theta$spid==i]), col=2)
  
  par(new=T)  
  plot(xdum, ydum, type="n", main=spptable[i,2]
       , xlim=c(lowlim, hilim)
       , xlab="swb gradient", ylab=expression("BAI" ~(cm^{2}))) #
  # growth curve
  polygon(c(xdum,rev(xdum)),c(mu_all$`97.5%`[mu_all$spid==i],rev(mu_all$`2.5%`[mu_all$spid==i])),col=trans_gray, border=NA)
  lines(xdum, mu_all$`50%`[mu_all$spid==i], col=1)
  
  abline(v=med_SWBs_all$opt_swb[med_SWBs_all$spid==i], col=1, lty=0)
  segments(x0=med_SWBs_all$opt_swb[med_SWBs_all$spid==i], y0=-21,x1=med_SWBs_all$opt_swb[med_SWBs_all$spid==i],y1=med_SWBs_all$g_hat_opt[med_SWBs_all$spid==i], col=1, lty=2)
  abline(h=med_SWBs_all$g_hat_opt[med_SWBs_all$spid==i], col=1, lty=0)
  segments(x0=-50, y0=med_SWBs_all$g_hat_opt[med_SWBs_all$spid==i],x1=med_SWBs_all$opt_swb[med_SWBs_all$spid==i],y1=med_SWBs_all$g_hat_opt[med_SWBs_all$spid==i], col=1, lty=2)
  abline(h=med_SWBs_all$g_hat_pres[med_SWBs_all$spid==i], col=2, lty=0)
  segments(x0=-50, y0=med_SWBs_all$g_hat_pres[med_SWBs_all$spid==i],x1=med_SWBs_all$pres_swb[med_SWBs_all$spid==i],y1=med_SWBs_all$g_hat_pres[med_SWBs_all$spid==i], col=2, lty=3)
  
  abline(v=med_SWBs_all$pres_swb[med_SWBs_all$spid==i], col=2, lty=0)
  segments(x0=med_SWBs_all$pres_swb[med_SWBs_all$spid==i], y0=-21,x1=med_SWBs_all$pres_swb[med_SWBs_all$spid==i],y1=med_SWBs_all$g_hat_pres[med_SWBs_all$spid==i], col=2, lty=3)
  
  legend("topright",
         legend=c("growth","occurrence"),
         lty=c(1,1), col=c("black","red"), cex=0.85)
  
}
dev.off()

## Growth data export----
# pareto_k_all$k_0.7_percent <- pareto_k_all$N_0.7/pareto_k_all$N_tot*100
# pareto_k_all <- pareto_k_all[order(-pareto_k_all$k_0.7_percent),]
# write.csv(pareto_k_all, 'grow pareto k.csv', row.names=F)
write.csv(betas_all, paste0(data.path, "grow_betas_swb_temp.csv"), row.names=F)

med_SWBs_all <- bind_cols(med_SWBs_all, as.data.frame(r2_all))
write.csv(med_SWBs_all, paste0(data.path, "growswbs_temp.csv"), row.names = F)

## all models are saved as RDS separately
remove(grow.fit)
save.image(paste0(save.path, "growth loop_temp.RData"))


## Survival data manipulation----
# *Clear workspace to start 
rm(list = ls())
# read in data
spptable <- na.omit(read.csv(paste0(data.path, 'spptable.new.csv'))[,c(1,5)])
names(spptable)[2] <- "gensp"
# spptable <- subset(spptable, gensp=="Pinus monophylla")
# spptable$spid <- 1
spptable$gensp <- as.character(spptable$gensp)
#spptable <- spptable[1:5,]
maxminopt <- read.csv(paste0(data.path, 'distromaxminopt_swb.csv'))
maxminopt <- merge(maxminopt, spptable, by="gensp")
theta <- read.csv(paste0(data.path, 'theta_swb.csv'))
theta$seq <- as.numeric(rownames(theta))
theta <- merge(theta, spptable, by="gensp")
theta <- theta[order(theta$seq),]
swb_axes_p <- read.csv(paste0(data.path, 'swb_axes_p.csv'))
swb_axes_p <- merge(swb_axes_p, spptable, by="gensp")
growswbs <- read.csv(paste0(data.path, 'growswbs_temp.csv'))
surv <- read.csv(paste0(data.path, 'surv.new.csv'))
tableinfo <- read.csv(paste0(data.path, 'tableinfo.csv'))

surv <- droplevels(subset(surv, gensp%in%spptable$gensp))

## If using a different subset of species, need to reset the IDs 

surv$spid <- as.numeric(as.factor(surv$gensp))

tableinfo_surv <- dplyr::group_by(surv, gensp) %>%
  dplyr::summarize(ntrees_s = length(unique(tree.uniqueID)),
            nplots_s = length(unique(plot.uniqueID)))
tableinfo <- merge(tableinfo, tableinfo_surv, by="gensp", all=T)
write.csv(tableinfo, paste0(data.path, 'tableinfo.csv', row.names = F))

surv2 <- group_by(surv, gensp) %>%
  mutate(mean = mean(ptcoeff, na.rm=T),
         sd = sd(ptcoeff, na.rm=T),
         upper = sd(ptcoeff, na.rm=T)*3 + mean(ptcoeff, na.rm=T),
         lower = mean(ptcoeff, na.rm=T) - sd(ptcoeff, na.rm=T)*3,
         chiscore = scores(ptcoeff, type="chisq", prob=0.999))

surv2$outlier <- ifelse(surv2$ptcoeff<surv2$upper&surv2$ptcoeff>surv2$lower, "no", "yes")

length(which(surv2$chiscore=="TRUE"))#/nrow(surv2)
length(which(surv2$chiscore=="TRUE"))/nrow(surv2)
length(which(surv2$outlier=="yes"))#/nrow(surv2)

surv3 <- subset(surv2, outlier=="yes")
# surv <- surv2

maxmin <- dplyr::group_by(surv, spid)%>%
  dplyr::summarise(min_swb = min(ptcoeff, na.rm=T),
                   max_swb = max(ptcoeff, na.rm=T))

N_spp <- length(unique(surv$spid))

pr_all <- matrix(NA, nrow = 200*N_spp, ncol=6)
pr_all <- data.frame(pr_all)

med_SWBs_all <- data.frame(matrix(NA, nrow=N_spp, ncol=26))

r2_all <- rep(NA, length=N_spp)

# pareto_k_all <- data.frame(matrix(NA, nrow=N_spp, ncol=3))
# pareto_k_all <- cbind(spptable, pareto_k_all)
# names(pareto_k_all)[3:5] <- c("N_0.7", "N_1", "N_tot")

surv1 <- surv
maxmin1 <- maxmin
maxminopt1 <- maxminopt

## Survival loop----
for(i in 1:N_spp){
  surv <- surv1[surv1$spid==i,]
  maxmin <- maxmin1[maxmin1$spid==i,]
  maxminopt <- maxminopt1[maxminopt1$spid==i,]

  ## scale covariates
  cov.mat.raw <- cbind(surv$dia, surv$ptcoeff, surv$avgba, surv$adcount, surv$temp)## data for Stan , surv$dist_nowCWD, surv$dist_histCWD
  cov.mat <- apply(cov.mat.raw, 2, scale)
  cov.means <- colMeans(cov.mat.raw, na.rm=T)
  cov.sds <- apply(cov.mat.raw, 2, sd, na.rm=T)

  min_swb <- (maxmin$min_swb-cov.means[2])/cov.sds[2]
  max_swb <- (maxmin$max_swb-cov.means[2])/cov.sds[2]

  ## put opt on same SWB scale
  opt_swb_p <- (maxminopt$opt_swb_p-cov.means[2])/cov.sds[2]

  # data for stan
  N <- dim(surv)[1]
  N_pred <- 200
  survival <- as.vector(surv$surv)
  dia <- cov.mat[,1]
  swb <- cov.mat[,2]
  swb_pred <- seq(min(cov.mat[,2],na.rm=T),max(cov.mat[,2],na.rm=T), length.out=N_pred)
  min_swb <- min_swb
  max_swb <- max_swb
  pres_swb_dat <- opt_swb_p
  avgba <- cov.mat[,3]
  adcount <- cov.mat[,4]
  temp <- cov.mat[,5]

  plot <- as.numeric(droplevels(as.factor(surv$plot.uniqueID)))
  N_plts <- length(unique(plot))


  surv_data <- list(N=N, N_plts=N_plts, plot=plot, min_swb=min_swb, max_swb=max_swb, pres_swb_dat=pres_swb_dat,
                    surv=survival, dia=dia, swb=swb, swb_pred=swb_pred, temp=temp,
                    avgba=avgba, adcount=adcount)

  ## Survival models----
  rstan_options(auto_write = TRUE)
  surv.compiled <- stan_model(file=paste0(mod.path, "surv_swb_temp_spp.stan"))
  
  pars=c("a_sp","b_dia","b_swb","b_ba","b_adcount","b2_swb", "b2_dia", "b_temp", "b2_temp"
         ,"sigma_c","wet_swb","dry_swb","pres_swb","par_swb", "opt_swb"
         ,"s_hat_wet","s_hat_dry","s_hat_opt", "s_hat_pres"
         , "pr_pred")#"pr_hat","log_lik"
  
  surv.fit<- sampling(
    object = surv.compiled
    , data=surv_data
    , chains=3
    , iter=2000
    , warmup=1000
    , thin=1
    , refresh=100
    , pars=pars
    , cores=3
    , save_warmup = FALSE
    , control= list(adapt_delta=0.95) #max_treedepth=12,
  )

  saveRDS(surv.fit, paste(save.path, "surv models/surv.fit",i,".rds", sep=""))
  
  ## Survival diagnostics----
  # # extract pareto-k
  # log_lik1 <- extract_log_lik(surv.fit, merge_chains = FALSE)
  # rel_n_eff <- relative_eff(exp(log_lik1))
  # loo1 <- loo(log_lik1, r_eff = rel_n_eff, cores = 2, save_psis = TRUE)
  # # plot(loo1, label_points = TRUE)
  # pareto_k_all$N_0.7[i] <- length(which(loo1$diagnostics$pareto_k>0.7))
  # pareto_k_all$N_1[i] <- length(which(loo1$diagnostics$pareto_k>1))
  # pareto_k_all$N_tot[i] <- length(loo1$diagnostics$pareto_k)
  posterior <- as.array(surv.fit)
  color_scheme_set("viridis")
  pdf(paste("Output pdfs/surv regression coeff/surv regression coeff", i, ".pdf", sep=""))
  # print(plot(loo1, label_points = TRUE))
  print(mcmc_intervals(posterior, pars=c("a_sp","b_dia","b_swb","b_ba","b_adcount","b_temp", "b2_dia","b2_swb","b2_temp","sigma_c")) +
          ggplot2::labs(title = spptable$gensp[i]))
  print(mcmc_combo(posterior, pars=c("a_sp","b_dia","b_swb","b_ba","b_adcount","b_temp")))
  print(mcmc_acf_bar(posterior, pars=c("a_sp","b_dia","b_swb","b_ba","b_adcount","b_temp")))
  print(mcmc_combo(posterior, pars=c("b2_dia","b2_swb","b2_temp","sigma_c","lp__")))
  print(mcmc_acf_bar(posterior, pars=c("b2_dia","b2_swb","b2_temp","sigma_c","lp__")))
  print(mcmc_rhat_hist(rhat(surv.fit)[1:10]))
  # print(mcmc_neff_hist(neff_ratio(surv.fit)[1:10]))
  dev.off()

  ## Survival output----
  # predicted curve
  pr <- as.data.frame(summary(surv.fit, pars="pr_pred")$summary[, c(1,4,6,8)])
  spid <- rep(i, N_pred)
  pr$spid <- spid
  pr <- merge(pr, spptable, by='spid')

  if(i==1) pr_all <- pr
  if(i>1) pr_all <- rbind(pr_all, pr)

  ## betas
  betas <- as.data.frame(summary(surv.fit, pars=c("b_dia","b_swb","b_ba","b_adcount","b_temp","b2_dia","b2_swb","b2_temp"),probs = c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))$summary[, c(1:10, 12)])
  betas$param <- rownames(betas)

  betas$spid <- i
  betas <- merge(betas, spptable, by='spid')

  if(i==1) betas_all <- betas
  if(i>1) betas_all <- rbind(betas_all, betas)

  ## SWB and subsequent s-hat values
  surv_swb <- as.data.frame(summary(surv.fit, pars="opt_swb")$summary[, c(1,4,5,6,7,8)])
  surv_swb_unscaled <- cov.means[2] + surv_swb*cov.sds[2]
  pres_swb <- summary(surv.fit, pars="pres_swb")$summary[,1]
  pres_swb_unscaled <- cov.means[2] + pres_swb*cov.sds[2]
  wet_swb <- summary(surv.fit, pars="wet_swb")$summary[,1]
  wet_swb_unscaled <- cov.means[2] + wet_swb*cov.sds[2]
  dry_swb <- summary(surv.fit, pars="dry_swb")$summary[,1]
  dry_swb_unscaled <- cov.means[2] + dry_swb*cov.sds[2]

  s_hat_surv <- inv.logit(summary(surv.fit, pars="s_hat_opt")$summary[, c(1,4,6,8)])
  s_hat_pres <- inv.logit(summary(surv.fit, pars="s_hat_pres")$summary[, c(1,4,6,8)])
  s_hat_wet <- inv.logit(summary(surv.fit, pars="s_hat_wet")$summary[, c(1,4,6,8)])
  s_hat_dry <- inv.logit(summary(surv.fit, pars="s_hat_dry")$summary[, c(1,4,6,8)])

  # Medians
  SWBs <- as.data.frame(t(c(surv_swb_unscaled[2:6,1], pres_swb_unscaled, wet_swb_unscaled, dry_swb_unscaled,
                              s_hat_surv[3], s_hat_pres[3], s_hat_wet[3], s_hat_dry[3])))
  names(SWBs) <- c('surv_swb_2.5','surv_swb_25','surv_swb','surv_swb_75','surv_swb_97.5','pres_swb','wet_swb','dry_swb','s_hat_surv','s_hat_pres','s_hat_wet','s_hat_dry')
  SWBs$spid <- i

  # Choose greater: parabola or endpoint optimum
  SWBs$surv_swb <- ifelse(SWBs$wet_swb>=SWBs$surv_swb&SWBs$surv_swb>=SWBs$dry_swb, SWBs$surv_swb,NA)
  SWBs$s_hat_surv <- ifelse(!is.na(SWBs$surv_swb), SWBs$s_hat_surv,NA)
  SWBs$opt_swb <- ifelse(!is.na(SWBs$s_hat_surv),ifelse(SWBs$s_hat_surv>SWBs$s_hat_wet & SWBs$s_hat_surv>SWBs$s_hat_dry
                                                        , SWBs$surv_swb,
                                                        ifelse(SWBs$s_hat_dry>SWBs$s_hat_wet
                                                               , SWBs$dry_swb, SWBs$wet_swb)),
                         ifelse(is.na(SWBs$s_hat_surv),ifelse(SWBs$s_hat_dry>SWBs$s_hat_wet
                                                              , SWBs$dry_swb, SWBs$wet_swb)))
  SWBs$opt_swb_2.5 <- ifelse(!is.na(SWBs$s_hat_surv),ifelse(SWBs$s_hat_surv>SWBs$s_hat_wet & SWBs$s_hat_surv>SWBs$s_hat_dry
                                                            , SWBs$surv_swb_2.5,
                                                            ifelse(SWBs$s_hat_dry>SWBs$s_hat_wet
                                                                   , SWBs$dry_swb, SWBs$wet_swb)),
                             ifelse(is.na(SWBs$s_hat_surv),ifelse(SWBs$s_hat_dry>SWBs$s_hat_wet
                                                                  , SWBs$dry_swb, SWBs$wet_swb)))
  SWBs$opt_swb_25 <- ifelse(!is.na(SWBs$s_hat_surv),ifelse(SWBs$s_hat_surv>SWBs$s_hat_wet & SWBs$s_hat_surv>SWBs$s_hat_dry
                                                           , SWBs$surv_swb_25,
                                                           ifelse(SWBs$s_hat_dry>SWBs$s_hat_wet
                                                                  , SWBs$dry_swb, SWBs$wet_swb)),
                            ifelse(is.na(SWBs$s_hat_surv),ifelse(SWBs$s_hat_dry>SWBs$s_hat_wet
                                                                 , SWBs$dry_swb, SWBs$wet_swb)))
  SWBs$opt_swb_75 <- ifelse(!is.na(SWBs$s_hat_surv),ifelse(SWBs$s_hat_surv>SWBs$s_hat_wet & SWBs$s_hat_surv>SWBs$s_hat_dry
                                                           , SWBs$surv_swb_75,
                                                           ifelse(SWBs$s_hat_dry>SWBs$s_hat_wet
                                                                  , SWBs$dry_swb, SWBs$wet_swb)),
                            ifelse(is.na(SWBs$s_hat_surv),ifelse(SWBs$s_hat_dry>SWBs$s_hat_wet
                                                                 , SWBs$dry_swb, SWBs$wet_swb)))

  SWBs$opt_swb_97.5 <- ifelse(!is.na(SWBs$s_hat_surv),ifelse(SWBs$s_hat_surv>SWBs$s_hat_wet & SWBs$s_hat_surv>SWBs$s_hat_dry
                                                             , SWBs$surv_swb_97.5,
                                                             ifelse(SWBs$s_hat_dry>SWBs$s_hat_wet
                                                                    , SWBs$dry_swb, SWBs$wet_swb)),
                              ifelse(is.na(SWBs$s_hat_surv),ifelse(SWBs$s_hat_dry>SWBs$s_hat_wet
                                                                   , SWBs$dry_swb, SWBs$wet_swb)))

  SWBs$s_hat_opt <- ifelse(!is.na(SWBs$s_hat_surv),ifelse(SWBs$s_hat_surv>SWBs$s_hat_wet & SWBs$s_hat_surv>SWBs$s_hat_dry
                                                          , SWBs$s_hat_surv,
                                                          ifelse(SWBs$s_hat_dry>SWBs$s_hat_wet
                                                                 , SWBs$s_hat_dry, SWBs$s_hat_wet)),
                           ifelse(is.na(SWBs$s_hat_surv),ifelse(SWBs$s_hat_dry>SWBs$s_hat_wet
                                                                , SWBs$s_hat_dry, SWBs$s_hat_wet)))
  # Choose minimum
  SWBs$min_swb <- ifelse(!is.na(SWBs$s_hat_surv),ifelse(SWBs$s_hat_surv<SWBs$s_hat_wet & SWBs$s_hat_surv<SWBs$s_hat_dry
                                                        , SWBs$surv_swb,
                                                        ifelse(SWBs$s_hat_dry<SWBs$s_hat_wet
                                                               , SWBs$dry_swb, SWBs$wet_swb)),
                         ifelse(is.na(SWBs$s_hat_surv),ifelse(SWBs$s_hat_dry<SWBs$s_hat_wet
                                                              , SWBs$dry_swb, SWBs$wet_swb)))

  SWBs$s_hat_min <- ifelse(!is.na(SWBs$s_hat_surv),ifelse(SWBs$s_hat_surv<SWBs$s_hat_wet & SWBs$s_hat_surv<SWBs$s_hat_dry
                                                          , SWBs$s_hat_surv,
                                                          ifelse(SWBs$s_hat_dry<SWBs$s_hat_wet
                                                                 , SWBs$s_hat_dry, SWBs$s_hat_wet)),
                           ifelse(is.na(SWBs$s_hat_surv),ifelse(SWBs$s_hat_dry<SWBs$s_hat_wet
                                                                , SWBs$s_hat_dry, SWBs$s_hat_wet)))
  # Calculate mismatch and cost
  SWBs$mismatch_p <- SWBs$opt_swb - SWBs$pres_swb
  SWBs$cost_p <- SWBs$s_hat_opt - SWBs$s_hat_pres
  SWBs$cost_p <- ifelse(SWBs$cost_p<0,0,SWBs$cost_p)

  SWBs$range <- SWBs$s_hat_opt - SWBs$s_hat_min
  SWBs$sens <- SWBs$range/(SWBs$wet_swb - SWBs$dry_swb)

  med_SWBs <- merge(SWBs, spptable, by='spid')
  med_SWBs_all[i,] <- med_SWBs
  names(med_SWBs_all) <- names(med_SWBs)

  print(paste("spp[",i,"] complete"), sep="")
  }

survswbs <- med_SWBs_all


## Survival figures----
trans_gray <- adjustcolor(col = "gray",alpha.f = 0.4)

pdf("Output pdfs/surv by species.pdf")
for(i in 1:N_spp){
  xdum <- seq(min(surv1$ptcoeff[surv1$spid==i],na.rm=T),max(surv1$ptcoeff[surv1$spid==i],na.rm=T), length.out=N_pred)
  xmean <- mean(surv1$ptcoeff[surv1$spid==i],na.rm=T)
  xsd <- sd(surv1$ptcoeff[surv1$spid==i],na.rm=T)
  ncuts<-20
  chopsize<-cut(surv1$ptcoeff[surv1$spid==i],ncuts)
  
  ###start here
  
  binmean.size<-as.numeric(sapply(split(surv1$ptcoeff[surv1$spid==i],chopsize),mean, na.rm=T))
  binmean.surv<-as.numeric(sapply(split(surv1$surv[surv1$spid==i],chopsize),mean, na.rm=T))
  plot(surv1$ptcoeff[surv1$spid==i], surv1$surv[surv1$spid==i], type="n", ylim=c(0, 1),
       main=spptable[i,2],
       xlab="swb gradient", 
       ylab="survival")
  points(binmean.size, binmean.surv)
  # points(swb[spp==i], jitter(pres[spp==i], 0.1))
  mtext("wet", side=1, line=2, adj=1)
  mtext("dry", side=1, line=2, adj=0)
  
  polygon(c(xdum,rev(xdum)),c(inv.logit(pr_all$`97.5%`[pr_all$spid==i]),rev(inv.logit(pr_all$`2.5%`[pr_all$spid==i]))),col=trans_gray, border=NA)
  lines(xdum, inv.logit(pr_all$`50%`[pr_all$spid==i]), col=4)
  # abline(v=med_SWBs_all$wet_swb[med_SWBs_all$spid==i], lwd=2)
  # abline(v=med_SWBs_all$dry_swb[med_SWBs_all$spid==i], lwd=2)
  abline(v=med_SWBs_all$opt_swb[med_SWBs_all$spid==i], col=4, lty=2)
  # lines(xdum, inv.logit(fixef(fit.bern)[1] + fixef(fit.bern)[3]*swb_pred + fixef(fit.bern)[8]*swb_pred^2))

}
dev.off()

pdf("Output pdfs/surv and occurrence.pdf")
for(i in 1:N_spp){
  xdum <- seq(min(surv1$ptcoeff[surv1$spid==i],na.rm=T),max(surv1$ptcoeff[surv1$spid==i],na.rm=T), length.out=N_pred)
  x2dum <- sort(swb_axes_p$swb_pred[swb_axes_p$spid==i])
  y2dum <- rep(seq(0,max(1), length.out=length(xdum)),N_spp)
  y2dum <- matrix(y2dum, length(xdum), N_spp)
  
  par(mar = c(5,5,2,5))
  plot(xdum, y2dum[,i], type="n", axes=F, xlab=NA, ylab=NA, 
       xlim=c(survswbs$dry_swb[survswbs$spid==i],survswbs$wet_swb[survswbs$spid==i]))
  axis(side = 4)
  mtext(side = 4, line = 3, 'Occurrence')
  # occurrence curve
  polygon(c(x2dum,rev(x2dum)),c(inv.logit(theta$X2.5.[theta$spid==i]),rev(inv.logit(theta$X97.5.[theta$spid==i]))),col=trans_gray, border=NA)
  lines(x2dum, inv.logit(theta$X50.[theta$spid==i]), col=2)
  
  par(new=T)  
  plot(surv1$ptcoeff[surv1$spid==i], surv1$surv[surv1$spid==i], type="n", ylim=c(0,1), main=spptable[i,2],xlab="swb gradient", ylab="p(Survival)",
       xlim=c(survswbs$dry_swb[survswbs$spid==i],survswbs$wet_swb[survswbs$spid==i]))
  # survival curve
  polygon(c(xdum,rev(xdum)),c(inv.logit(pr_all$`97.5%`[pr_all$spid==i]),rev(inv.logit(pr_all$`2.5%`[pr_all$spid==i]))),col=trans_gray, border=NA)
  lines(xdum, inv.logit(pr_all$`50%`[pr_all$spid==i]), col=1)
  
  abline(v=survswbs$opt_swb[survswbs$spid==i], col=1, lty=0)
  segments(x0=survswbs$opt_swb[survswbs$spid==i], y0=0,x1=survswbs$opt_swb[survswbs$spid==i],y1=(survswbs$s_hat_opt[survswbs$spid==i]), col=1, lty=2)
  abline(h=(survswbs$s_hat_opt[survswbs$spid==i]), col=1, lty=0)
  segments(x0=-50, y0=(survswbs$s_hat_opt[survswbs$spid==i]),x1=survswbs$opt_swb[survswbs$spid==i],y1=(survswbs$s_hat_opt[survswbs$spid==i]), col=1, lty=2)
  abline(h=(survswbs$s_hat_pres[survswbs$spid==i]), col=1, lty=0)
  segments(x0=-50, y0=(survswbs$s_hat_pres[survswbs$spid==i]),x1=survswbs$pres_swb[survswbs$spid==i],y1=(survswbs$s_hat_pres[survswbs$spid==i]), col=2, lty=2)
  
  abline(v=survswbs$pres_swb[survswbs$spid==i], col=2, lty=0)
  segments(x0=survswbs$pres_swb[survswbs$spid==i], y0=-11,x1=survswbs$pres_swb[survswbs$spid==i],y1=(survswbs$s_hat_pres[survswbs$spid==i]), col=2, lty=2)
  
  legend("topright",
         legend=c("survival","occurrence"),
         lty=c(1,1), col=c("black","red"), cex=0.85)
}
dev.off()

## Survival data export----
# pareto_k_all$k_0.7_percent <- pareto_k_all$N_0.7/pareto_k_all$N_tot*100
# pareto_k_all <- pareto_k_all[order(-pareto_k_all$k_0.7_percent),]
# write.csv(pareto_k_all, 'surv pareto k.csv', row.names=F)
med_SWBs_all <- bind_cols(med_SWBs_all, as.data.frame(r2_all))
write.csv(med_SWBs_all, paste0(data.path, 'survswbs_temp.csv'), row.names = F)
write.csv(betas_all, paste0(data.path, 'surv_betas_swb_temp.csv'), row.names=F)

## all models are saved as RDS separately
remove(surv.fit)
save.image(paste0(save.path, 'survival model_temp.RData'), compress=T)



## Recruitment data manipulation----
# *Clear workspace to start 
rm(list = ls())
# read in data
spptable <- na.omit(read.csv(paste0(data.path, 'spptable.new.csv'))[,c(1,6)])
names(spptable)[2] <- "gensp"
## subset for debugging
# spptable <- subset(spptable, spid%in%c(1:2))

maxminopt <- read.csv(paste0(data.path, 'distromaxminopt_swb.csv'))
maxminopt <- merge(maxminopt, spptable, by="gensp")
theta <- read.csv(paste0(data.path, 'theta_swb.csv'))
theta$seq <- as.numeric(rownames(theta))
theta <- merge(theta, spptable, by="gensp")
theta <- theta[order(theta$seq),]
swb_axes_p <- read.csv(paste0(data.path, 'swb_axes_p.csv'))
swb_axes_p <- merge(swb_axes_p, spptable, by="gensp")
growswbs <- read.csv(paste0(data.path, 'growswbs_temp.csv'))
survswbs <- read.csv(paste0(data.path, 'survswbs_temp.csv'))
recruits <- dplyr::select(read.csv(paste0(data.path, 'recruits.new.csv')), -spid)
# recruits <- rename(recruits, "cwd"="ptcoeff")
tableinfo <- read.csv(paste0(data.path, 'tableinfo.csv'))

recruits <- merge(recruits, spptable)

##is the problem repeat measurements?
recruits$spp_plot_code <- as.factor(paste(recruits$spid, recruits$plot.uniqueID, sep='.'))
recruits$spp_plot <- as.numeric(recruits$spp_plot_code)
recruits <- recruits[order(recruits$plotyear),]
recruits <- recruits[!duplicated(recruits$spp_plot),]

maxrec <- group_by(recruits, gensp) %>%
  summarize(max = max(N.sdlgs, na.rm=T),
            mean = mean(N.sdlgs[which(N.sdlgs>0)], na.rm=T))
write.table(maxrec, paste0(data.path, "max_rec_obs.csv"), row.names = F)

rec2 <- group_by(recruits, gensp) %>%
  mutate(mean = mean(ptcoeff, na.rm=T),
         sd = sd(ptcoeff, na.rm=T),
         upper = sd(ptcoeff, na.rm=T)*3 + mean(ptcoeff, na.rm=T),
         lower = mean(ptcoeff, na.rm=T) - sd(ptcoeff, na.rm=T)*3,
         chiscore = scores(ptcoeff, type="chisq", prob=0.999))

rec2$outlier <- ifelse(rec2$ptcoeff<rec2$upper&rec2$ptcoeff>rec2$lower, "no", "yes")

length(which(rec2$chiscore=="TRUE"))#/nrow(rec2)
length(which(rec2$chiscore=="TRUE"))/nrow(rec2)
length(which(rec2$outlier=="yes"))#/nrow(rec2)

# recruits <- subset(rec2, chiscore!="TRUE")

tableinfo_rec <- dplyr::group_by(recruits, gensp) %>%
  dplyr::summarize(nplots_r = length(unique(plot.uniqueID)),
            npres_r = sum(pres, na.rm=T))
tableinfo <- merge(tableinfo, tableinfo_rec, by="gensp", all=T)
write.csv(tableinfo, paste0(data.path, 'tableinfo.csv'), row.names = F)

## scale covariates 
cov.mat.raw <- cbind(recruits$ptcoeff, recruits$avgba, recruits$adcount, recruits$temp)## data for Stan , rec$dist_nowCWD, rec$dist_histCWD
cov.mat <- apply(cov.mat.raw, 2, scale)
cov.means <- colMeans(cov.mat.raw, na.rm=T)
cov.sds <- apply(cov.mat.raw, 2, sd, na.rm=T)

maxmin <- dplyr::group_by(recruits, spid)%>%
  dplyr::summarise(min_swb = min(ptcoeff, na.rm=T),
                   max_swb = max(ptcoeff, na.rm=T))
min_swb <- (maxmin$min_swb-cov.means[1])/cov.sds[1]
max_swb <- (maxmin$max_swb-cov.means[1])/cov.sds[1]

## put opt on same swb scale
opt_swb_p <- (maxminopt$opt_swb_p-cov.means[1])/cov.sds[1]

# data for stan
N <- dim(recruits)[1]
N_pred <- 200
N_spp <- length(unique(recruits$spid))
spp <- recruits$spid
sdlgs <- as.vector(recruits$N.sdlgs)
swb <- cov.mat[,1]
swb_pred <- seq(min(cov.mat[,1],na.rm=T),max(cov.mat[,1],na.rm=T), length.out=N_pred)
min_swb <- min_swb
max_swb <- max_swb
pres_swb_dat <- opt_swb_p
avgba <- cov.mat[,2]
adcount <- cov.mat[,3]
temp <- cov.mat[,4]

swb_pred_r <- cov.means[1] + swb_pred*cov.sds[1]

recruit_data <- list(
  N=N,  N_spp=N_spp
  ,swb=swb, adcount= adcount, avgba =  avgba, temp=temp
  ,sdlgs=sdlgs
  ,spp=spp
  ,swb_pred=swb_pred
  ,min_swb=min_swb, max_swb=max_swb, pres_swb_dat=pres_swb_dat
  ,N_pred=N_pred)

pars=c("b_swb","b_ba","b_adcount","b2_swb", "b_temp", "b2_temp"  
       ,"phi"
       ,"a_pres", "a_abund"
       ,"wet_swb","dry_swb","pres_swb","par_swb", "opt_swb"
       ,"r_hat_wet","r_hat_dry", "r_hat_pres", "r_hat_opt"
       , "mu_pred"
) #, "sdlgs_rep", "log_lik"

## Recruitment models----
recruit.compiled <- stan_model(file=paste0(mod.path, "recruit_plot_swb_temp.stan"))

recruit.fit<- sampling(
  object = recruit.compiled
  , data=recruit_data
  , chains=3
  , iter=1000
  , warmup=500
  , thin=1
  , refresh=100
  , pars=pars
  , cores=3
  , control= list(adapt_delta=0.9) #, max_treedepth=18
)


saveRDS(recruit.fit, paste0(save.path, "recruit_mod_stan.rds"))

## Recruitment diagnostics----
# log_lik1 <- extract_log_lik(recruit.fit, merge_chains = FALSE)
# rel_n_eff <- relative_eff(exp(log_lik1))
# loo1 <- loo(log_lik1, r_eff = rel_n_eff, cores = 2, save_psis = TRUE)
# # plot(loo1, label_points = TRUE)
# 
# pareto_k_all <- dplyr::select(recruits, c(spid, gensp, N.sdlgs, ptcoeff))
# 
# pareto_k_all$pareto_k <- loo1$diagnostics$pareto_k
# 
# pareto_k_spp <- group_by(pareto_k_all, spid) %>%
#   summarize(gensp = first(gensp),
#             N_0.7 = length(which(pareto_k>0.7)),
#             N_1 = length(which(pareto_k>1)),
#             N_tot = length(pareto_k),
#             k_0.7_percent = N_0.7/N_tot*100)

inits <- get_inits(recruit.fit)
inits_chain1 <- inits[[1]]
print(inits_chain1)

posterior <- as.array(recruit.fit)
print(head(posterior))

color_scheme_set("viridis")

pdf("Output pdfs/recruitparams.pdf")
for(i in 1:N_spp){
  pars1 <- c(paste("a_pres[", i,"]", sep=""),
             paste("a_abund[", i,"]", sep=""),
             paste("phi[", i,"]", sep=""))
             
  pars2 <- c(paste("b_swb[", i,"]", sep=""),
             paste("b_temp[", i,"]", sep=""),
             paste("b_ba[", i,"]", sep=""),
             paste("b_adcount[", i,"]", sep=""),
             paste("b2_temp[", i,"]", sep=""),
             paste("b2_swb[", i,"]", sep=""))
  parsall <- c(pars1, pars2)
  print(mcmc_intervals(posterior, pars=parsall))
  print(mcmc_combo(posterior, pars=pars1))
  print(mcmc_acf_bar(posterior, pars=pars1))
  print(mcmc_combo(posterior, pars=c(pars2, "lp__")))
  print(mcmc_acf_bar(posterior, pars=c(pars2, "lp__")))
  # print(mcmc_rhat_hist(rhat(grow.fit)[1:9]))
  # print(mcmc_neff_hist(neff_ratio(grow.fit)[1:9]))
}
dev.off()

## Recruitment output----
betas <- as.data.frame(summary(recruit.fit, pars=c("b_swb","b_temp","b_ba","b_adcount","b2_swb","b2_temp"),probs = c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))$summary[, c(1:10, 12)])
betas$param <- rownames(betas)
betas$spid <- rep(seq(1,N_spp, by=1),6)
betas <- merge(betas, spptable, by='spid')
betas$param <- sapply(strsplit(betas$param, '\\['), "[", 1 )

mu_rec <- as.data.frame(summary(recruit.fit, pars="mu_pred")$summary[, c(1,4,6,8)])
spid <- sort(rep(seq(1,N_spp,by=1),200))
mu_rec$spid <- spid
mu_rec <- merge(mu_rec, spptable, by="spid")

## SWB and subsequent g-hat values
rec_swb <- as.data.frame(summary(recruit.fit, pars="opt_swb")$summary[, c(1,4:8)])
rec_swb_unscaled <- cov.means[1] + rec_swb*cov.sds[1]
pres_swb <- summary(recruit.fit, pars="pres_swb")$summary[,1]
pres_swb_unscaled <- cov.means[1] + pres_swb*cov.sds[1]
wet_swb <- summary(recruit.fit, pars="wet_swb")$summary[,1]
wet_swb_unscaled <- cov.means[1] + wet_swb*cov.sds[1]
dry_swb <- summary(recruit.fit, pars="dry_swb")$summary[,1]
dry_swb_unscaled <- cov.means[1] + dry_swb*cov.sds[1]

r_hat_rec <- (summary(recruit.fit, pars="r_hat_opt")$summary[, c(1,4,6,8)])
r_hat_pres <- (summary(recruit.fit, pars="r_hat_pres")$summary[, c(1,4,6,8)])
r_hat_wet <- (summary(recruit.fit, pars="r_hat_wet")$summary[, c(1,4,6,8)])
r_hat_dry <- (summary(recruit.fit, pars="r_hat_dry")$summary[, c(1,4,6,8)])

# Medians
SWBs <- as.data.frame(cbind(rec_swb_unscaled[,2:6], pres_swb_unscaled, wet_swb_unscaled, dry_swb_unscaled, 
                            r_hat_rec[,3], r_hat_pres[,3], r_hat_wet[,3], r_hat_dry[,3]))
names(SWBs) <- c('rec_swb_2.5','rec_swb_25','rec_swb','rec_swb_75','rec_swb_97.5','pres_swb','wet_swb','dry_swb','r_hat_rec','r_hat_pres','r_hat_wet','r_hat_dry')
SWBs$spid <- seq(1,N_spp,by=1)

# Choose greater: parabola or endpoint optimum
SWBs$rec_swb <- ifelse(SWBs$wet_swb>=SWBs$rec_swb&SWBs$rec_swb>=SWBs$dry_swb, SWBs$rec_swb,NA)
SWBs$r_hat_rec <- ifelse(SWBs$wet_swb>=SWBs$rec_swb&SWBs$rec_swb>=SWBs$dry_swb, SWBs$r_hat_rec,NA)
SWBs$opt_swb <- ifelse(!is.na(SWBs$r_hat_rec),ifelse(SWBs$r_hat_rec>SWBs$r_hat_wet & SWBs$r_hat_rec>SWBs$r_hat_dry
                                                     , SWBs$rec_swb, 
                                                     ifelse(SWBs$r_hat_dry>SWBs$r_hat_wet 
                                                            , SWBs$dry_swb, SWBs$wet_swb)), 
                       ifelse(is.na(SWBs$r_hat_rec),ifelse(SWBs$r_hat_dry>SWBs$r_hat_wet 
                                                           , SWBs$dry_swb, SWBs$wet_swb)))
SWBs$opt_swb_2.5 <- ifelse(!is.na(SWBs$r_hat_rec),ifelse(SWBs$r_hat_rec>SWBs$r_hat_wet & SWBs$r_hat_rec>SWBs$r_hat_dry
                                                         , SWBs$rec_swb_2.5, 
                                                         ifelse(SWBs$r_hat_dry>SWBs$r_hat_wet 
                                                                , SWBs$dry_swb, SWBs$wet_swb)), 
                           ifelse(is.na(SWBs$r_hat_rec),ifelse(SWBs$r_hat_dry>SWBs$r_hat_wet 
                                                               , SWBs$dry_swb, SWBs$wet_swb)))

SWBs$opt_swb_25 <- ifelse(!is.na(SWBs$r_hat_rec),ifelse(SWBs$r_hat_rec>SWBs$r_hat_wet & SWBs$r_hat_rec>SWBs$r_hat_dry
                                                        , SWBs$rec_swb_25, 
                                                        ifelse(SWBs$r_hat_dry>SWBs$r_hat_wet 
                                                               , SWBs$dry_swb, SWBs$wet_swb)), 
                          ifelse(is.na(SWBs$r_hat_rec),ifelse(SWBs$r_hat_dry>SWBs$r_hat_wet 
                                                              , SWBs$dry_swb, SWBs$wet_swb)))

SWBs$opt_swb_75 <- ifelse(!is.na(SWBs$r_hat_rec),ifelse(SWBs$r_hat_rec>SWBs$r_hat_wet & SWBs$r_hat_rec>SWBs$r_hat_dry
                                                        , SWBs$rec_swb_75, 
                                                        ifelse(SWBs$r_hat_dry>SWBs$r_hat_wet 
                                                               , SWBs$dry_swb, SWBs$wet_swb)), 
                          ifelse(is.na(SWBs$r_hat_rec),ifelse(SWBs$r_hat_dry>SWBs$r_hat_wet 
                                                              , SWBs$dry_swb, SWBs$wet_swb)))

SWBs$opt_swb_97.5 <- ifelse(!is.na(SWBs$r_hat_rec),ifelse(SWBs$r_hat_rec>SWBs$r_hat_wet & SWBs$r_hat_rec>SWBs$r_hat_dry
                                                          , SWBs$rec_swb_97.5, 
                                                          ifelse(SWBs$r_hat_dry>SWBs$r_hat_wet 
                                                                 , SWBs$dry_swb, SWBs$wet_swb)), 
                            ifelse(is.na(SWBs$r_hat_rec),ifelse(SWBs$r_hat_dry>SWBs$r_hat_wet 
                                                                , SWBs$dry_swb, SWBs$wet_swb)))


SWBs$r_hat_opt <- ifelse(!is.na(SWBs$r_hat_rec),ifelse(SWBs$r_hat_rec>SWBs$r_hat_wet & SWBs$r_hat_rec>SWBs$r_hat_dry
                                                       , SWBs$r_hat_rec, 
                                                       ifelse(SWBs$r_hat_dry>SWBs$r_hat_wet 
                                                              , SWBs$r_hat_dry, SWBs$r_hat_wet)), 
                         ifelse(is.na(SWBs$r_hat_rec),ifelse(SWBs$r_hat_dry>SWBs$r_hat_wet 
                                                             , SWBs$r_hat_dry, SWBs$r_hat_wet)))

# Choose minimum
SWBs$min_swb <- ifelse(!is.na(SWBs$r_hat_rec),ifelse(SWBs$r_hat_rec<SWBs$r_hat_wet & SWBs$r_hat_rec<SWBs$r_hat_dry
                                                     , SWBs$rec_swb, 
                                                     ifelse(SWBs$r_hat_dry<SWBs$r_hat_wet 
                                                            , SWBs$dry_swb, SWBs$wet_swb)), 
                       ifelse(is.na(SWBs$r_hat_rec),ifelse(SWBs$r_hat_dry<SWBs$r_hat_wet 
                                                           , SWBs$dry_swb, SWBs$wet_swb)))

SWBs$r_hat_min <- ifelse(!is.na(SWBs$r_hat_rec),ifelse(SWBs$r_hat_rec<SWBs$r_hat_wet & SWBs$r_hat_rec<SWBs$r_hat_dry
                                                       , SWBs$r_hat_rec, 
                                                       ifelse(SWBs$r_hat_dry<SWBs$r_hat_wet 
                                                              , SWBs$r_hat_dry, SWBs$r_hat_wet)), 
                         ifelse(is.na(SWBs$r_hat_rec),ifelse(SWBs$r_hat_dry<SWBs$r_hat_wet 
                                                             , SWBs$r_hat_dry, SWBs$r_hat_wet)))

SWBs <-  merge(SWBs, spptable, by='spid')
SWBs <- merge(SWBs, maxrec[1:2], by="gensp")
SWBs$r_hat_opt <- ifelse(SWBs$r_hat_opt<SWBs$max, SWBs$r_hat_opt, SWBs$max)

SWBs$mismatch_p <- SWBs$opt_swb - SWBs$pres_swb
SWBs$cost_p <- SWBs$r_hat_opt - SWBs$r_hat_pres
SWBs$cost_p <- ifelse(SWBs$cost_p<0,0,SWBs$cost_p)

SWBs$range <- SWBs$r_hat_opt - SWBs$r_hat_min
SWBs$sens <- SWBs$range/(SWBs$wet_swb-SWBs$dry_swb)

med_SWBs <- SWBs

recswbs <- med_SWBs

## Recruitment figures----
swb_raw <- cov.means[1] + swb*cov.sds[1]
xdum <- swb_pred_r
ydum <- rep(seq(0, 50, length.out=length(xdum)),N_spp)
ydum <- matrix(ydum, length(xdum), N_spp)
trans_gray <- adjustcolor(col = "gray",alpha.f = 0.4)

pdf(file="Output pdfs/rec by species.pdf")
for(i in(1:N_spp)){
  
  ncuts<-10
  chopsize<-cut(swb_raw[spp==i],ncuts)
  binmean.size<-as.numeric(sapply(split(swb_raw[spp==i],chopsize),mean, na.rm=T))
  binmean.sdlgs<-as.numeric(sapply(split(sdlgs[spp==i],chopsize),mean, na.rm=T))
  
  yhi <- max(binmean.sdlgs, na.rm=T) +5
  plot(swb_raw[spp==i], sdlgs[spp==i], type="n", 
       main=spptable[i,2],
       xlab="swb gradient", 
       ylab="seedlings",
       ylim=c(0, yhi))
  points(binmean.size, binmean.sdlgs)
  mtext("wet", side=1, line=2, adj=1)
  mtext("dry", side=1, line=2, adj=0)

  polygon(c(xdum,rev(xdum)),c(mu_rec$`97.5%`[mu_rec$spid==i],rev(mu_rec$`2.5%`[mu_rec$spid==i])),col=trans_gray, border=NA)
  lines(xdum, mu_rec$`50%`[mu_rec$spid==i], col=4)
  abline(v=med_SWBs$wet_swb[med_SWBs$spid==i], lwd=2)
  abline(v=med_SWBs$dry_swb[med_SWBs$spid==i], lwd=2)
  abline(v=med_SWBs$opt_swb[med_SWBs$spid==i], col=2, lty=2)
  abline(v=med_SWBs$opt_swb[med_SWBs$spid==i], col=4, lty=2)
  
  legend("topright",
         legend=c("mean","limits"),
         lty=c(2,1), lwd=c(1,2), col=c("red","black"), cex=0.85)
  
}
dev.off()

y2dum <- rep(seq(0,max(1), length.out=length(xdum)),N_spp)
y2dum <- matrix(y2dum, length(xdum), N_spp)

pdf("Output pdfs/rec and occurrence.pdf")
for(i in 1:N_spp){
  x2dum <- sort(swb_axes_p$swb_pred[swb_axes_p$spid==i])
  ydum <- seq(0,(max(recswbs$r_hat_opt[recswbs$spid==i])+10), length.out=length(xdum))
  par(mar = c(5,5,2,5))
  plot(xdum, y2dum[,i], type="n", axes=F, xlab=NA, ylab=NA,
       xlim=c(recswbs$dry_swb[recswbs$spid==i],recswbs$wet_swb[recswbs$spid==i]))
  axis(side = 4)
  mtext(side = 4, line = 3, 'Occurrence')
  # occurrence curve
  polygon(c(x2dum,rev(x2dum)),c(inv.logit(theta$X2.5.[theta$spid==i]),rev(inv.logit(theta$X97.5.[theta$spid==i]))),col=trans_gray, border=NA)
  lines(x2dum, inv.logit(theta$X50.[theta$spid==i]), col=2)

  par(new=T)
  plot(xdum, ydum, type="n", main=spptable[i,2],xlab="swb gradient", ylab="Seedlings",
       xlim=c(recswbs$dry_swb[recswbs$spid==i],recswbs$wet_swb[recswbs$spid==i]))
  # recruitment curve
  polygon(c(xdum,rev(xdum)),c(mu_rec$`97.5%`[mu_rec$spid==i],rev(mu_rec$`2.5%`[mu_rec$spid==i])),col=trans_gray, border=NA)
  lines(xdum, mu_rec$`50%`[mu_rec$spid==i], col=1)

  abline(v=recswbs$opt_swb[recswbs$spid==i], col=1, lty=0)
  segments(x0=recswbs$opt_swb[recswbs$spid==i], y0=-11,x1=recswbs$opt_swb[recswbs$spid==i],y1=recswbs$r_hat_opt[recswbs$spid==i], col=1, lty=2)
  abline(h=recswbs$r_hat_opt[recswbs$spid==i], col=1, lty=0)
  segments(x0=-50, y0=recswbs$r_hat_opt[recswbs$spid==i],x1=recswbs$opt_swb[recswbs$spid==i],y1=recswbs$r_hat_opt[recswbs$spid==i], col=1, lty=2)
  abline(h=recswbs$r_hat_pres[recswbs$spid==i], col=1, lty=0)
  segments(x0=-50, y0=recswbs$r_hat_pres[recswbs$spid==i],x1=recswbs$pres_swb[recswbs$spid==i],y1=recswbs$r_hat_pres[recswbs$spid==i], col=2, lty=2)

  abline(v=recswbs$pres_swb[recswbs$spid==i], col=2, lty=0)
  segments(x0=recswbs$pres_swb[recswbs$spid==i], y0=-11,x1=recswbs$pres_swb[recswbs$spid==i],y1=recswbs$r_hat_pres[recswbs$spid==i], col=2, lty=2)

  legend("topright",
         legend=c("recruitment","occurrence"),
         lty=c(1,1), col=c("black","red"), cex=0.85)
}
dev.off()

## Recruitment data export----
recswbs <- merge(recswbs, r2_all, by=c("spid", 'gensp'))
write.csv(recswbs, paste0(data.path, 'recswbs_temp.csv'), row.names = F)
write.csv(r2_all, paste0(data.path, "rec_r2.csv"), row.names=F)
write.csv(betas, paste0(data.path, 'rec_betas_swb_temp.csv'), row.names = F)

save.image(paste0(save.path, 'repro model_temp.RData'), compress=T)

## allvital stuff----
recswbs <- read.csv(paste0(data.path, 'recswbs_temp.csv'))
growswbs <- read.csv(paste0(data.path, 'growswbs_temp.csv'))
survswbs <- read.csv(paste0(data.path, 'survswbs_temp.csv'))
spptable <- na.omit(read.csv(paste0(data.path, 'spptable.new.csv'))[,c(1,2)])
names(spptable)[2] <- "gensp"
maxminopt <- read.csv(paste0(data.path, 'distromaxminopt_swb.csv')) %>%
  rename(dry_swb_p = min_swb_p, wet_swb_p = max_swb_p)

surv_all <- merge(survswbs, maxminopt, by="gensp") %>%
rename(opt_swb_s = opt_swb, opt_swb_s_2.5 = opt_swb_2.5, opt_swb_s_25 = opt_swb_25, 
       opt_swb_s_75 = opt_swb_75, opt_swb_s_97.5 = opt_swb_97.5,
       dry_swb_s = dry_swb, wet_swb_s = wet_swb)

rec_all <- merge(recswbs, maxminopt, by="gensp") %>%
  rename(opt_swb_r = opt_swb, opt_swb_r_2.5 = opt_swb_2.5, opt_swb_r_25 = opt_swb_25, 
         opt_swb_r_75 = opt_swb_75, opt_swb_r_97.5 = opt_swb_97.5,
         dry_swb_r = dry_swb, wet_swb_r = wet_swb)

grow_all <- merge(growswbs, maxminopt, by="gensp") %>%
  rename(opt_swb_g = opt_swb, opt_swb_g_2.5 = opt_swb_2.5, opt_swb_g_25 = opt_swb_25, 
         opt_swb_g_75 = opt_swb_75, opt_swb_g_97.5 = opt_swb_97.5,
         dry_swb_g = dry_swb, wet_swb_g = wet_swb)

vr_all <- dplyr::select(surv_all, gensp, spid, wet_swb_s, dry_swb_s, opt_swb_s, opt_swb_s_2.5, 
                        opt_swb_s_25, opt_swb_s_75, opt_swb_s_97.5,
                        dry_swb_p, wet_swb_p, opt_swb_p, opt_swb_p_2.5, 
                        opt_swb_p_25, opt_swb_p_75, opt_swb_p_97.5) %>%
  merge(dplyr::select(rec_all, gensp, wet_swb_r, dry_swb_r, opt_swb_r, opt_swb_r_2.5, 
                      opt_swb_r_25, opt_swb_r_75, opt_swb_r_97.5), by="gensp", all.x=T) %>%
  merge(dplyr::select(grow_all, gensp, wet_swb_g, dry_swb_g, opt_swb_g, opt_swb_g_2.5, 
                      opt_swb_g_25, opt_swb_g_75, opt_swb_g_97.5), by="gensp", all.x=T)


vr_all2 <- melt(vr_all, id.vars=c("gensp", "spid"), measure.vars=c("opt_swb_p","opt_swb_s","opt_swb_g","opt_swb_r"))%>%
  rename(vital_rate=variable, median=value)

vr_all2$interval_2.5 <- melt(vr_all, id.vars="gensp", measure.vars=c("opt_swb_p_2.5","opt_swb_s_2.5","opt_swb_g_2.5","opt_swb_r_2.5"))$value
vr_all2$interval_25 <- melt(vr_all, id.vars="gensp", measure.vars=c("opt_swb_p_25","opt_swb_s_25","opt_swb_g_25","opt_swb_r_25"))$value
vr_all2$interval_75 <- melt(vr_all, id.vars="gensp", measure.vars=c("opt_swb_p_75","opt_swb_s_75","opt_swb_g_75","opt_swb_r_75"))$value
vr_all2$interval_97.5 <- melt(vr_all, id.vars="gensp", measure.vars=c("opt_swb_p_97.5","opt_swb_s_97.5","opt_swb_g_97.5","opt_swb_r_97.5"))$value
vr_all2$wet_swb <- melt(vr_all, id.vars="gensp", measure.vars=c("wet_swb_p","wet_swb_s","wet_swb_g","wet_swb_r"))$value
vr_all2$dry_swb <- melt(vr_all, id.vars="gensp", measure.vars=c("dry_swb_p","dry_swb_s","dry_swb_g","dry_swb_r"))$value


pdf("Output pdfs/all vitals SWB with 95 CI.pdf")
for(i in 1:10){
  vr_all_sub<- subset(vr_all2, spid==i*6|spid==i*6-1|spid==i*6-2|spid==i*6-3|spid==i*6-4|spid==i*6-5)
  
  print(ggplot(vr_all_sub, aes(vital_rate, median, fill=vital_rate)) +
          geom_segment(aes(x=vital_rate,y=interval_2.5, xend=vital_rate, yend=interval_97.5), colour="gray") +
          geom_segment(aes(x=vital_rate,y=interval_25, xend=vital_rate, yend=interval_75, colour=gensp)) +
          geom_point(aes(color=gensp)) + 
          coord_flip() +
          facet_wrap(~gensp, ncol = 1)+
          scale_x_discrete(labels=c("opt_swb_g" = "growth", "opt_swb_r" = "recruitment",
                                    "opt_swb_s" = "survival", "opt_swb_p" = "occurrence")) +
          ylab(label="median optimum SWB value") +
          xlab(label="")+
          theme(panel.grid.major.y = element_line(linetype="dotted", color="gray", size = 0.3)) +
          theme(legend.position = "none")
  )
}
dev.off()

vr_pres <- subset(vr_all2, vital_rate=="opt_swb_p")%>%
  dplyr::select(gensp, median)%>%
  rename("pres"="median")


vr_all3 <- merge(vr_all2, vr_pres, by="gensp")
vr_all3$sig <- ifelse(vr_all3$interval_2.5>vr_all3$pres&
                        vr_all3$interval_97.5>vr_all3$pres, 1,
                      ifelse(vr_all3$interval_2.5<vr_all3$pres&
                               vr_all3$interval_97.5<vr_all3$pres, 1,0))
vr_all3$sig[which(vr_all3$vital_rate=="opt_swb_p")] <- NA

pdf("Output pdfs/all vitals sig SWB with 95 CI.pdf")
for(i in 1:10){
  vr_all_sub<- subset(vr_all3, spid==i*6|spid==i*6-1|spid==i*6-2|spid==i*6-3|spid==i*6-4|spid==i*6-5)
  
  print(ggplot(vr_all_sub, aes(vital_rate, median, fill=vital_rate)) +
          geom_segment(aes(x=vital_rate,y=interval_2.5, xend=vital_rate, yend=interval_97.5, alpha=sig), colour="gray") +
          geom_segment(aes(x=vital_rate,y=interval_25, xend=vital_rate, yend=interval_75, colour=gensp, alpha=sig)) +
          geom_point(aes(color=gensp, alpha=sig)) + 
          scale_alpha_continuous(range=c(0.3,1))+
          coord_flip() +
          facet_wrap(~gensp, ncol = 1)+
          scale_x_discrete(labels=c("opt_swb_g" = "growth", "opt_swb_r" = "recruitment",
                                    "opt_swb_s" = "survival", "opt_swb_p" = "occurrence")) +
          ylab(label="median optimum SWB value") +
          geom_hline(aes(yintercept=pres), linetype="dotted") +
          xlab(label="")+
          theme(panel.grid.major.y = element_line(linetype="dotted", color="gray", size = 0.3)) +
          theme(legend.position = "none")
  )
}
dev.off()

vr_all4 <- dplyr::select(vr_all, gensp, opt_swb_p, opt_swb_s, opt_swb_g, opt_swb_r
                         , wet_swb_p, wet_swb_s, wet_swb_g, wet_swb_r
                         , dry_swb_p, dry_swb_s, dry_swb_g, dry_swb_r)%>%
  gather(var, val, -gensp) %>%
  separate(var, into=c("A","B","vr"), sep = '_') %>%
  unite("var", c("A","B")) %>%
  spread(var, val) %>%
  drop_na()
  

vr_all4$maximum <- ifelse(vr_all4$opt_swb==vr_all4$wet_swb, "wet",
                  ifelse(vr_all4$opt_swb==vr_all4$dry_swb, "dry", "intermediate"))

pdf("Output pdfs/all vr shape summary.pdf")
ggplot(vr_all4, aes(x=vr)) +
  geom_bar(aes(fill=maximum), stat="count")
dev.off()

