## 5_Trait phylo brms
## Phylogenetic regressions (brms) of mismatch and cost on traits, with LOO model selection, and the trait figure.
## Author: Teresa Bohner
## Project: where western trees grow vs. where they thrive (see README)

## packages----
library(tidyverse)
library(taxize)
library(brranching)
library(MCMCglmm)
library(brms)
library(bayesplot)
library(gridExtra)
library(ggtree)
library(tidybayes)
library(ggplot2)
library(ggstance)
library(cowplot)


## NCBI API key for taxize/brranching: set ENTREZ_KEY in ~/.Renviron (never in code)
getkey(service="entrez")

## 1: Read in data files ----
data.path <- "Processed Data/Analysis Ready Data/"
table.path <- "Output tables/"
fig.path <- "Output pdfs/"

all.data0 <- read.csv(paste0(data.path, "all_mismatch_and_traits.csv"))
all.data <- dplyr::select(all.data0, gensp,  vr,mismatch_p, abmis_p, cost_p, wood_density, sla_value, value_RootDepth, logSeedMass, value_Lifespan, MaxHeightM, med_plot_swb)

summary(all.data)

spplist <- as.character(unique(all.data$gensp))

spplist_phylocom = brranching::phylomatic_names(spplist, 
                                                format = "isubmit", 
                                                db = "ncbi")
species <- as.data.frame(cbind(spplist, spplist_phylocom))
names(species) <- c("gensp", "phylo")
species$phylo_gensp <- str_split_fixed(species$phylo, pattern="/", 3)[,3]

all.data <- merge(all.data, species, by="gensp")

tr2 <- brranching::phylomatic(spplist_phylocom, storedtree = "zanne2014", outformat="newick")

pdf(paste0(fig.path, "phlogenetic tree.pdf"), width = 6, height=14)
plot(tr2)
dev.off()
#NOTE: 2 taxa not matched: NA/pinaceae/abies/abies_shastensis/pinaceae/abies/abies_shastensis, 
#NA/pinaceae/larix/larix_lyallii/pinaceae/larix/larix_lyallii, ]

all.data <- filter(all.data, !gensp%in%c("Abies shastensis", "Larix lyallii"))

tr3 <- tr2
tr3$node.label <- "NA"
tr3$tip.label <- tolower(tr3$tip.label)

tr3 <- phytools::force.ultrametric(tr3)

inv.phylo <- MCMCglmm::inverseA(tr3, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

## Abundance model ----
abun.data <- filter(all.data, vr=="0abun")
head(abun.data)

for(i in 1:nrow(abun.data)){
  abun.data$tot_traits[i] <- length(which(!is.na(abun.data[i, 6:11])))
}

abun.data$scaleSWB <- scale(abun.data$med_plot_swb)
mod.a1 <- brm(mismatch_p~ scaleSWB, data=abun.data, control=list(adapt_delta=0.99))
mod.a2 <- brm(abmis_p~ scaleSWB, data=abun.data, family=lognormal(), control=list(adapt_delta=0.99))

X0 <- dplyr::select(abun.data, c(wood_density, sla_value, value_RootDepth, logSeedMass, value_Lifespan, MaxHeightM))
X <- scale(X0)

abun.data1 <- as.data.frame(cbind(dplyr::select(abun.data, c(phylo_gensp, mismatch_p)), X)) %>%
  rename(response=mismatch_p)
abun.data2 <- as.data.frame(cbind(dplyr::select(abun.data, c(phylo_gensp, abmis_p)), X)) %>%
  rename(response=abmis_p)

bform <- bf(response ~ (1|phylo_gensp) + mi(wood_density)  +  mi(sla_value) + mi(value_RootDepth) + mi(logSeedMass) + mi(value_Lifespan) + mi(MaxHeightM)) + 
  bf(wood_density | mi() ~ 1) + set_rescor(FALSE) +
  bf(sla_value | mi() ~ 1) + set_rescor(FALSE) +
  bf(value_RootDepth | mi() ~ 1) + set_rescor(FALSE) +
  bf(logSeedMass | mi() ~ 1) + set_rescor(FALSE) +
  bf(value_Lifespan | mi() ~ 1) + set_rescor(FALSE) +
  bf(MaxHeightM | mi() ~ 1) + set_rescor(FALSE)

fit_abun <- brm(bform, data = abun.data1, cov_ranef = list(phylo_gensp=A), control = list(adapt_delta=0.99, max_treedepth = 12))
fit_abun1 <- brm(bform, data = abun.data1, control = list(adapt_delta=0.99, max_treedepth = 12))

bform2 <- bf(response ~ (1|phylo_gensp) + mi(wood_density)  +  mi(sla_value) + mi(value_RootDepth) + 
               mi(logSeedMass) + mi(value_Lifespan) + mi(MaxHeightM), family=lognormal()) + 
  bf(wood_density | mi() ~ 1) + set_rescor(FALSE) +
  bf(sla_value | mi() ~ 1) + set_rescor(FALSE) + 
  bf(value_RootDepth | mi() ~ 1) + set_rescor(FALSE) + 
  bf(logSeedMass | mi() ~ 1) + set_rescor(FALSE) + 
  bf(value_Lifespan | mi() ~ 1) + set_rescor(FALSE) +
  bf(MaxHeightM | mi() ~ 1) + set_rescor(FALSE) 

fit_abun2 <- brm(bform2, data = abun.data2, cov_ranef = list(phylo_gensp=A), control = list(adapt_delta=0.99, max_treedepth = 12))
fit_abun3 <- brm(bform2, data = abun.data2, control = list(adapt_delta=0.99, max_treedepth = 14))

hyp <- "sd_phylo_gensp__response_Intercept^2 / (sd_phylo_gensp__response_Intercept^2 + sigma_response^2) = 0"
(hyp <- hypothesis(fit_abun, hyp, class = NULL))
plot(hyp)

hyp <- "sd_phylo_gensp__response_Intercept^2 / (sd_phylo_gensp__response_Intercept^2 + sigma_response^2) = 0"
(hyp <- hypothesis(fit_abun2, hyp, class = NULL))
plot(hyp)

lambdas <- data.frame("vr"=as.character("0abun"),
                      "lambda"=hyp$hypothesis$Estimate)
lambdas$vr <- as.character(lambdas$vr)
  
loo.0 <- loo(fit_abun, resp = "response")
loo.1 <- loo(fit_abun1, resp = "response")
compare(loo.0, loo.1) # negative favors first, postive favors second
loo_model_weights(list(loo.0, loo.1), method = "pseudobma")

loo.2 <- loo(fit_abun2, resp = "response")
loo.3 <- loo(fit_abun3, resp = "response")
compare(loo.2, loo.3)
loo_model_weights(list(loo.2, loo.3), method = "pseudobma")


## Growth model ----
growth.data <- filter(all.data, vr=="2growth")
head(growth.data)
growth.data$scaleSWB <- scale(growth.data$med_plot_swb)
mod.g1 <- brm(mismatch_p~ scaleSWB, data=growth.data, control=list(adapt_delta=0.99))
mod.g2 <- brm(abmis_p~ scaleSWB, data=growth.data, family=lognormal(), control=list(adapt_delta=0.99))


X0 <- dplyr::select(growth.data, c(wood_density, sla_value, value_RootDepth, logSeedMass, value_Lifespan, MaxHeightM))
X <- scale(X0)

growth.data1 <- as.data.frame(cbind(dplyr::select(growth.data, c(phylo_gensp, mismatch_p)), X)) %>%
  rename(response=mismatch_p)
growth.data2 <- as.data.frame(cbind(dplyr::select(growth.data, c(phylo_gensp, abmis_p)), X)) %>%
  rename(response=abmis_p)
  
fit_grow <- brm(bform, data = growth.data1, cov_ranef = list(phylo_gensp=A), control = list(adapt_delta=0.99, max_treedepth = 12))
fit_grow1 <- brm(bform, data = growth.data1, control = list(adapt_delta=0.99, max_treedepth = 12))

fit_grow2 <- brm(bform2, data = growth.data2, cov_ranef = list(phylo_gensp=A), control = list(adapt_delta=0.99, max_treedepth = 12))
fit_grow3 <- brm(bform2, data = growth.data2, control = list(adapt_delta=0.99, max_treedepth = 12))

hyp <- "sd_phylo_gensp__response_Intercept^2 / (sd_phylo_gensp__response_Intercept^2 + sigma_response^2) = 0"
(hyp <- hypothesis(fit_grow, hyp, class = NULL))
plot(hyp)

hyp <- "sd_phylo_gensp__response_Intercept^2 / (sd_phylo_gensp__response_Intercept^2 + sigma_response^2) = 0"
(hyp <- hypothesis(fit_grow2, hyp, class = NULL))
plot(hyp)

lambdas <- rbind(lambdas, c("2growth", hyp$hypothesis$Estimate))


loo.0 <- loo(fit_grow, resp = "response")
loo.1 <- loo(fit_grow1, resp = "response")
compare(loo.0, loo.1) 
loo_model_weights(list(loo.0, loo.1), method = "pseudobma")

loo.2 <- loo(fit_grow2, resp = "response")
loo.3 <- loo(fit_grow3, resp = "response")
compare(loo.2, loo.3) 
loo_model_weights(list(loo.2, loo.3), method = "pseudobma")


## Survival model ----
surv.data <- filter(all.data, vr=="1surv")
head(surv.data)

surv.data$scaleSWB <- scale(surv.data$med_plot_swb)
mod.s1 <- brm(mismatch_p~ scaleSWB, data=surv.data, control=list(adapt_delta=0.99))
mod.s2 <- brm(abmis_p~ scaleSWB, data=surv.data, family=lognormal(), control=list(adapt_delta=0.99))

X0 <- dplyr::select(surv.data, c(wood_density, sla_value, value_RootDepth, logSeedMass, value_Lifespan, MaxHeightM))
X <- scale(X0)

surv.data1 <- as.data.frame(cbind(dplyr::select(surv.data, c(phylo_gensp, mismatch_p)), X)) %>%
  rename(response=mismatch_p)

surv.data2 <- as.data.frame(cbind(dplyr::select(surv.data, c(phylo_gensp, abmis_p)), X)) %>%
  rename(response=abmis_p)

fit_surv <- brm(bform, data = surv.data1, cov_ranef = list(phylo_gensp=A), control = list(adapt_delta=0.99, max_treedepth = 12))
fit_surv1 <- brm(bform, data = surv.data1, control = list(adapt_delta=0.99, max_treedepth = 12))

fit_surv2 <- brm(bform2, data = surv.data2, cov_ranef = list(phylo_gensp=A), control = list(adapt_delta=0.99, max_treedepth = 12))
fit_surv3 <- brm(bform2, data = surv.data2, control = list(adapt_delta=0.99, max_treedepth = 12))

hyp <- "sd_phylo_gensp__response_Intercept^2 / (sd_phylo_gensp__response_Intercept^2 + sigma_response^2) = 0"
(hyp <- hypothesis(fit_surv, hyp, class = NULL))
plot(hyp)

hyp <- "sd_phylo_gensp__response_Intercept^2 / (sd_phylo_gensp__response_Intercept^2 + sigma_response^2) = 0"
(hyp <- hypothesis(fit_surv2, hyp, class = NULL))
plot(hyp)

lambdas <- rbind(lambdas, c("1surv", hyp$hypothesis$Estimate))

loo.0 <- loo(fit_surv, resp = "response")
loo.1 <- loo(fit_surv1, resp = "response")
compare(loo.0, loo.1) 
loo_model_weights(list(loo.0, loo.1), method = "pseudobma")

loo.2 <- loo(fit_surv2, resp = "response")
loo.3 <- loo(fit_surv3, resp = "response")
compare(loo.2, loo.3) 
loo_model_weights(list(loo.2, loo.3), method = "pseudobma")


## Recruitment model ----
rec.data <- filter(all.data, vr=="3rec")
head(rec.data)

rec.data$scaleSWB <- scale(rec.data$med_plot_swb)
mod.r1 <- brm(mismatch_p~ scaleSWB, data=rec.data, control=list(adapt_delta=0.99))
mod.r2 <- brm(abmis_p~ scaleSWB, data=rec.data, family=lognormal(), control=list(adapt_delta=0.99))

X0 <- dplyr::select(rec.data, c(wood_density, sla_value, value_RootDepth, logSeedMass, value_Lifespan, MaxHeightM))
X <- scale(X0)

rec.data1 <- as.data.frame(cbind(dplyr::select(rec.data, c(phylo_gensp, mismatch_p)), X)) %>%
  rename(response=mismatch_p)
rec.data2 <- as.data.frame(cbind(dplyr::select(rec.data, c(phylo_gensp, abmis_p)), X)) %>%
  rename(response=abmis_p)

fit_rec <- brm(bform, data = rec.data1, cov_ranef = list(phylo_gensp=A), control = list(adapt_delta=0.99, max_treedepth = 12))
fit_rec1 <- brm(bform, data = rec.data1, control = list(adapt_delta=0.99, max_treedepth = 12))

fit_rec2 <- brm(bform2, data = rec.data2, cov_ranef = list(phylo_gensp=A), control = list(adapt_delta=0.99, max_treedepth = 12))
fit_rec3 <- brm(bform2, data = rec.data2, control = list(adapt_delta=0.99, max_treedepth = 12))

hyp <- "sd_phylo_gensp__response_Intercept^2 / (sd_phylo_gensp__response_Intercept^2 + sigma_response^2) = 0"
(hyp <- hypothesis(fit_rec, hyp, class = NULL))
plot(hyp)

hyp <- "sd_phylo_gensp__response_Intercept^2 / (sd_phylo_gensp__response_Intercept^2 + sigma_response^2) = 0"
(hyp <- hypothesis(fit_rec2, hyp, class = NULL))
plot(hyp)

lambdas <- rbind(lambdas, c("3rec", hyp$hypothesis$Estimate))

loo.0 <- loo(fit_rec, resp = "response")
loo.1 <- loo(fit_rec1, resp = "response")
compare(loo.0, loo.1) 
loo_model_weights(list(loo.0, loo.1), method = "pseudobma")

loo.2 <- loo(fit_rec2, resp = "response")
loo.3 <- loo(fit_rec3, resp = "response")
compare(loo.2, loo.3) 
loo_model_weights(list(loo.2, loo.3), method = "pseudobma")

write.csv(lambdas, paste0(table.path, "phylo_lambda.csv"), row.names=F)

## Interval plots ----


abun1 <- fit_abun%>%
  gather_draws(c(bsp_response_miwood_density, bsp_response_misla_value, bsp_response_mivalue_RootDepth, bsp_response_milogSeedMass, bsp_response_mivalue_Lifespan, bsp_response_miMaxHeightM)) %>%
  median_qi(.width = c(.95, .9, .85)) %>%
  mutate(vr="Density")

grow1 <- fit_grow%>%
  gather_draws(c(bsp_response_miwood_density, bsp_response_misla_value, bsp_response_mivalue_RootDepth, bsp_response_milogSeedMass, bsp_response_mivalue_Lifespan, bsp_response_miMaxHeightM)) %>%
  median_qi(.width = c(.95, .9, .85)) %>%
  mutate(vr="Growth")

surv1 <- fit_surv%>%
  gather_draws(c(bsp_response_miwood_density, bsp_response_misla_value, bsp_response_mivalue_RootDepth, bsp_response_milogSeedMass, bsp_response_mivalue_Lifespan, bsp_response_miMaxHeightM)) %>%
  median_qi(.width = c(.95, .9, .85))%>%
  mutate(vr="Survival")

rec1 <- fit_rec%>%
  gather_draws(c(bsp_response_miwood_density, bsp_response_misla_value, bsp_response_mivalue_RootDepth, bsp_response_milogSeedMass, bsp_response_mivalue_Lifespan, bsp_response_miMaxHeightM)) %>%
  median_qi(.width = c(.95, .9, .85))%>%
  mutate(vr="Recruitment") 

all_vr_raw <- rbind(abun1, grow1, surv1, rec1) %>%
  mutate(sig95=ifelse(.width==0.95&.lower>0&.upper>0, "*", ifelse(.width==0.95&.lower<0&.upper<0,"*", "")),
         sig90=ifelse(.width==0.90&.lower>0&.upper>0, "*", ifelse(.width==0.90&.lower<0&.upper<0,"*", "")),
         sig85=ifelse(.width==0.85&.lower>0&.upper>0, "*", ifelse(.width==0.85&.lower<0&.upper<0,"*", "")))

plot_raw <- ggplot(all_vr_raw, aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointintervalh(aes(color=vr), size_range = c(0.5, 5), position = position_dodgev(-0.5)) +
  scale_color_manual(values = c("Density"="purple","Growth"="dark green", "Survival" = "blue", "Recruitment" = "red")) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels=c("log (Seed Mass)", "Maximum Height","SLA", "Lifespan", "Root Depth (min)", "Wood Density")) +
  ylab("") + xlab("Effect on Mismatch") +
  geom_text(aes(x=-26, label = sig95, color = vr, size=1.75), position = position_dodgev(-0.5)) +
  geom_text(aes(x=-25.12, label = sig90, color = vr, size=1.75), position = position_dodgev(-0.5)) +
  geom_text(aes(x=-24.24, label = sig85, color = vr, size=1.75), position = position_dodgev(-0.5)) +
  theme(legend.position = "none") +
  xlim(c(-26, 26)) +
  ggtitle("a)") +
  theme(axis.text=element_text(size=13, color="black"), axis.title=element_text(size=15), plot.title = element_text(size=20, face="bold"))

abun2 <- fit_abun2%>%
  gather_draws(c(bsp_response_miwood_density, bsp_response_misla_value, bsp_response_mivalue_RootDepth, bsp_response_milogSeedMass, bsp_response_mivalue_Lifespan, bsp_response_miMaxHeightM)) %>%
  median_qi(.width = c(.95, .9, .85)) %>%
  mutate(vr="Density") 

grow2 <- fit_grow2%>%
  gather_draws(c(bsp_response_miwood_density, bsp_response_misla_value, bsp_response_mivalue_RootDepth, bsp_response_milogSeedMass, bsp_response_mivalue_Lifespan, bsp_response_miMaxHeightM)) %>%
  median_qi(.width = c(.95, .9, .85)) %>%
  mutate(vr="Growth")

surv2 <- fit_surv2%>%
  gather_draws(c(bsp_response_miwood_density, bsp_response_misla_value, bsp_response_mivalue_RootDepth, bsp_response_milogSeedMass, bsp_response_mivalue_Lifespan, bsp_response_miMaxHeightM)) %>%
  median_qi(.width = c(.95, .9, .85))%>%
  mutate(vr="Survival")

rec2 <- fit_rec2%>%
  gather_draws(c(bsp_response_miwood_density, bsp_response_misla_value, bsp_response_mivalue_RootDepth, bsp_response_milogSeedMass, bsp_response_mivalue_Lifespan, bsp_response_miMaxHeightM)) %>%
  median_qi(.width = c(.95, .9, .85))%>%
  mutate(vr="Recruitment") 

all_vr_abs <- rbind(abun2, grow2, surv2, rec2) %>%
  mutate(sig95=ifelse(.width==0.95&.lower>0&.upper>0, "*", ifelse(.width==0.95&.lower<0&.upper<0,"*", "")),
         sig90=ifelse(.width==0.90&.lower>0&.upper>0, "*", ifelse(.width==0.90&.lower<0&.upper<0,"*", "")),
         sig85=ifelse(.width==0.85&.lower>0&.upper>0, "*", ifelse(.width==0.85&.lower<0&.upper<0,"*", "")))

plot_abs <- ggplot(all_vr_abs, aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointintervalh(aes(color=vr), size_range = c(0.5, 5), position = position_dodgev(-0.5)) +
  scale_color_manual(values = c("Density"="purple","Growth"="dark green", "Survival" = "blue", "Recruitment" = "red")) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels=c("log (Seed Mass)", "Maximum Height","SLA", "Lifespan", "Root Depth (min)", "Wood Density")) +
  ylab("") + xlab("Effect on Mismatch Magnitude") +
  geom_text(aes(x=-1.5, label = sig95, color = vr, size=1.75), position = position_dodgev(-0.5)) +
  geom_text(aes(x=-1.45, label = sig90, color = vr, size=1.75), position = position_dodgev(-0.5)) +
  geom_text(aes(x=-1.4, label = sig85, color = vr, size=1.75), position = position_dodgev(-0.5)) +
  xlim(c(-1.5,1.5)) +
  theme(legend.position = "none") +
  ggtitle("b)") +
  theme(axis.text=element_text(size=13, color="black"), axis.title=element_text(size=15), plot.title = element_text(size=20, face="bold"))

forlegend <- ggplot(all_vr_abs, aes (x = .value, fill=vr)) +
  geom_bar(aes(fill=vr)) +
  scale_fill_manual(name= "",
                    values = c("Density"="purple","Growth"="dark green", "Survival" = "blue", "Recruitment" = "red"))

legend <- get_legend(forlegend)
grid.arrange(legend)

pdf(paste0(fig.path, "fig4_intervals_brms_no_phylo.pdf"), width = 13, height = 5)
ggdraw(plot_grid(plot_grid(plot_raw, plot_abs, nrow=1, align='h'),
                 plot_grid(legend, nrow=1),
                 rel_widths=c(1, 0.13)))
dev.off()

abun3 <- mod.a1%>%
  gather_draws(b_scaleSWB) %>%
  median_qi(.width = c(.95, .9, .85)) %>%
  mutate(vr="Density")
surv3 <- mod.s1%>%
  gather_draws(b_scaleSWB) %>%
  median_qi(.width = c(.95, .9, .85)) %>%
  mutate(vr="Survival")
grow3 <- mod.g1%>%
  gather_draws(b_scaleSWB) %>%
  median_qi(.width = c(.95, .9, .85)) %>%
  mutate(vr="Growth")
rec3 <- mod.r1%>%
  gather_draws(b_scaleSWB) %>%
  median_qi(.width = c(.95, .9, .85)) %>%
  mutate(vr="Recruitment")
all_vr3 <- rbind(abun3, grow3, surv3, rec3) %>%
  mutate(sig95=ifelse(.width==0.95&.lower>0&.upper>0, "*", ifelse(.width==0.95&.lower<0&.upper<0,"*", "")),
         sig90=ifelse(.width==0.90&.lower>0&.upper>0, "*", ifelse(.width==0.90&.lower<0&.upper<0,"*", "")),
         sig85=ifelse(.width==0.85&.lower>0&.upper>0, "*", ifelse(.width==0.85&.lower<0&.upper<0,"*", "")))


abun4 <- mod.a2%>%
  gather_draws(b_scaleSWB) %>%
  median_qi(.width = c(.95, .9, .85)) %>%
  mutate(vr="Density")
surv4 <- mod.s2%>%
  gather_draws(b_scaleSWB) %>%
  median_qi(.width = c(.95, .9, .85)) %>%
  mutate(vr="Survival")
grow4 <- mod.g2%>%
  gather_draws(b_scaleSWB) %>%
  median_qi(.width = c(.95, .9, .85)) %>%
  mutate(vr="Growth")
rec4 <- mod.r2%>%
  gather_draws(b_scaleSWB) %>%
  median_qi(.width = c(.95, .9, .85)) %>%
  mutate(vr="Recruitment")
all_vr4 <- rbind(abun4, grow4, surv4, rec4) %>%
  mutate(sig95=ifelse(.width==0.95&.lower>0&.upper>0, "*", ifelse(.width==0.95&.lower<0&.upper<0,"*", "")),
         sig90=ifelse(.width==0.90&.lower>0&.upper>0, "*", ifelse(.width==0.90&.lower<0&.upper<0,"*", "")),
         sig85=ifelse(.width==0.85&.lower>0&.upper>0, "*", ifelse(.width==0.85&.lower<0&.upper<0,"*", "")))

all_vr_both <- bind_rows(all_vr3, all_vr_raw)
all_vr_both2 <- bind_rows(all_vr4, all_vr_abs)

plot_swb_raw <- ggplot(all_vr_both, aes(y = reorder(.variable, desc(.variable)), x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointintervalh(aes(color=vr), size_range = c(0.5, 5), position = position_dodgev(-0.5)) +
  scale_color_manual(values = c("Density"="purple","Growth"="dark green", "Survival" = "blue", "Recruitment" = "red")) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels=c("Wood Density", "Root Depth (min)", "Lifespan", "SLA", "Maximum Height", "log (Seed Mass)", "Median plot SWB")) +
  ylab("") + xlab("Effect on Mismatch") +
  geom_text(aes(x=-26, label = sig95, color = vr, size=1.75), position = position_dodgev(-0.5)) +
  geom_text(aes(x=-25.12, label = sig90, color = vr, size=1.75), position = position_dodgev(-0.5)) +
  geom_text(aes(x=-24.24, label = sig85, color = vr, size=1.75), position = position_dodgev(-0.5)) +
  theme(legend.position = "none") +
  xlim(c(-26, 26)) +
  # ggtitle("a)") +
  geom_hline(yintercept = 6.5) +
  theme(axis.text=element_text(size=13, color="black"), axis.title=element_text(size=15), plot.title = element_text(size=20, face="bold"))


plot_swb_abs <- ggplot(all_vr_both2, aes(y = reorder(.variable, desc(.variable)), x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointintervalh(aes(color=vr), size_range = c(0.5, 5), position = position_dodgev(-0.5)) +
  scale_color_manual(values = c("Density"="purple","Growth"="dark green", "Survival" = "blue", "Recruitment" = "red")) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels=c("Wood Density", "Root Depth (min)", "Lifespan", "SLA", "Maximum Height", "log (Seed Mass)", "Median plot SWB")) +
  ylab("") + xlab("Effect on Mismatch Magnitude") +
  geom_text(aes(x=-1.5, label = sig95, color = vr, size=1.75), position = position_dodgev(-0.5)) +
  geom_text(aes(x=-1.45, label = sig90, color = vr, size=1.75), position = position_dodgev(-0.5)) +
  geom_text(aes(x=-1.4, label = sig85, color = vr, size=1.75), position = position_dodgev(-0.5)) +
  xlim(c(-1.5,1.5)) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 6.5) +
  # ggtitle("a)") +
  theme(axis.text=element_text(size=13, color="black"), axis.title=element_text(size=15), plot.title = element_text(size=20, face="bold"))


pdf(paste0(fig.path, "fig4_intervals_brms_raw.pdf"), width = 8, height = 5)
ggdraw(plot_grid(plot_grid(plot_swb_raw, nrow=1, align='h'),
                 plot_grid(legend, nrow=1),
                 rel_widths=c(1, 0.3)))
dev.off()

pdf(paste0(fig.path, "fig4_intervals_brms_abs.pdf"), width = 8, height = 5)
ggdraw(plot_grid(plot_grid(plot_swb_abs, nrow=1, align='h'),
                 plot_grid(legend, nrow=1),
                 rel_widths=c(1, 0.3)))
dev.off()

