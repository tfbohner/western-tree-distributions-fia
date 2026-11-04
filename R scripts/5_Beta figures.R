## 5_Beta figures
## Figures of regression coefficients by species, traits and family.
## Author: Teresa Bohner
## Project: where western trees grow vs. where they thrive (see README)

library(bbmle)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(ggalt)
library(stringr)
library(reshape2)
library(dplyr)
library(doBy)
library(data.table)
library(lme4)
library(MASS)
library(ape)
library(knitr)
library(stargazer)
## not sure if we need all of these

source("R scripts/5_Beta summaries.R")

## ggplot labelling miscellainy
pred_names <- list(
  b1 = "Stem Count", 
  b2 = "Mean (BA)", 
  b3 = "SWB", 
  b4 = expression("SWB"^2), 
  b5 = "Temperature",
  b6 = expression("Temperature"^2),
  b7 = "Diameter", 
  b8 = expression("Diameter"^2)
)
pred_labeller <- function(variable,value){
  return(pred_names[value])
  
}

dens_names <- list(
  b1 = "SWB", 
  b2 = expression("SWB"^2)
)
dens_labeller <- function(variable,value){
  return(dens_names[value])
  
}


## Visualize regression coefficients----
pdf("Output pdfs/param grids.pdf", width = 8.5, height = 11)
ggplot(sbetas, aes(betas, reorder(gensp, -spid))) +
  geom_tile(aes(alpha=signif),fill="medium blue", color="black") + 
  guides(alpha=guide_legend(title="Proportion of posterior \nnot overlapping zero")) +
  coord_fixed(ratio=0.5)+
  theme_classic() +
  xlab(label="Survival") + ylab(label="") +
  scale_x_discrete(position = "top", labels=paramnames$predictor) +
  theme(axis.text.x = element_text(angle = 60, hjust = 0)) +
  theme(axis.text=element_text(size=11), axis.title.x = element_text(size=12, face="bold"))

ggplot(sbetas, aes(betas, reorder(gensp, -spid, fill=sign, alpha=signif))) +
  geom_tile(aesfill=interaction(as.factor(sign), as.factor(signif)), color="black") + 
  scale_fill_manual(
    values=c(alpha("red3", .1), alpha("red3",.7), alpha("red3",1),
             alpha("medium blue",.1), alpha("medium blue",.7), alpha("medium blue",1))) +
  # labels=c("Negative, NS", "Negative >90%", "Negative >95%",
  #          "Positive, NS", "Positive >90%", "Positive >95%"),
  # guides(alpha=guide_legend(title="Proportion of posterior \nnot overlapping zero")) +
  # guides(fill=guide_legend(title="Sign of median \nparameter estimate")) +
  coord_fixed(ratio=0.5)+
  theme_classic() +
  xlab(label="Survival") + ylab(label="") +
  scale_x_discrete(position = "top", labels=paramnames$predictor) +
  theme(axis.text.x = element_text(angle = 60, hjust = 0)) +
  theme(axis.text=element_text(size=11), axis.title.x = element_text(size=12, face="bold"))

ggplot(gbetas, aes(betas, reorder(gensp, -spid))) +
  geom_tile(aes(alpha=signif), fill="dark green", color="black") + 
  guides(alpha=guide_legend(title="Proportion of posterior \nnot overlapping zero")) +
  coord_fixed(ratio=0.5)+
  theme_classic() +
  xlab(label="Growth") + ylab(label="") +
  scale_x_discrete(position = "top", labels=paramnames$predictor) +
  theme(axis.text.x = element_text(angle = 60, hjust = 0)) +
  theme(axis.text=element_text(size=11), axis.title.x = element_text(size=12, face="bold"))
ggplot(gbetas, aes(betas, reorder(gensp, -spid))) +
  geom_tile(aes(alpha=signif, fill=sign), color="black") + 
  scale_fill_manual(values=c("red3", "medium blue")) +
  guides(alpha=guide_legend(title="Proportion of posterior \nnot overlapping zero")) +
  guides(fill=guide_legend(title="Sign of median \nparameter estimate")) +
  coord_fixed(ratio=0.5)+
  theme_classic() +
  xlab(label="Growth") + ylab(label="") +
  scale_x_discrete(position = "top", labels=paramnames$predictor) +
  theme(axis.text.x = element_text(angle = 60, hjust = 0)) +
  theme(axis.text=element_text(size=11), axis.title.x = element_text(size=12, face="bold"))

ggplot(rbetas, aes(betas, reorder(gensp, -spid))) +
  geom_tile(aes(alpha=signif), fill="red3", color="black") + 
  guides(alpha=guide_legend(title="Proportion of posterior \nnot overlapping zero")) +
  coord_fixed(ratio=0.5)+
  theme_classic() +
  xlab(label="Recruitment") + ylab(label="") +
  scale_x_discrete(position = "top", labels=paramnames2$predictor) +
  theme(axis.text.x = element_text(angle = 60, hjust = 0)) +
  theme(axis.text=element_text(size=11), axis.title.x = element_text(size=12, face="bold"))
ggplot(rbetas, aes(betas, reorder(gensp, -spid))) +
  geom_tile(aes(alpha=signif, fill=sign), color="black") + 
  scale_fill_manual(values=c("red3", "medium blue")) +
  guides(alpha=guide_legend(title="Proportion of posterior \nnot overlapping zero")) +
  guides(fill=guide_legend(title="Sign of median \nparameter estimate")) +
  coord_fixed(ratio=0.5)+
  theme_classic() +
  xlab(label="Recruitment") + ylab(label="") +
  scale_x_discrete(position = "top", labels=paramnames$predictor) +
  theme(axis.text.x = element_text(angle = 60, hjust = 0)) +
  theme(axis.text=element_text(size=11), axis.title.x = element_text(size=12, face="bold"))

dev.off()

pdf("Output pdfs/param histograms.pdf")
ggplot(allbetas, aes(X50., fill=param, alpha=signif)) +
  geom_histogram(color="black") +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  facet_grid(param~vr, scales="free") +
  theme_classic() +
  xlab("Median parameter estimate") + ylab("Number of species") +
  theme(legend.position = "none")

ggplot(allbetas, aes(X50., fill=param, alpha=signif)) +
  geom_histogram(color="black") +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  facet_wrap(param~vr, ncol=3, scales="free") +
  theme_classic() +
  xlab("Median parameter estimate") + ylab("Number of species") +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(legend.position = "none")

ggplot(allbetas, aes(X50., stat(count/3), fill=param, alpha=signif)) +
  geom_density(size=0.5, position="stack") +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  facet_grid(param~vr, scales="free") +
  # facet_wrap(~param+vr, ncol = 3, scales="free") +
  theme_classic() +
  xlab("Median parameter estimate") + ylab("Number of species") +
  theme(legend.position = "none") 

ggplot(allbetas, aes(X50., stat(count/3), fill=param, alpha=signif)) +
  geom_density(size=0.5, position="stack") +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  # facet_grid(param~vr, scales="free") +
  facet_wrap(~param+vr, ncol = 3, scales="free") +
  theme_classic() +
  xlab("Median parameter estimate") + ylab("Number of species") +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(legend.position = "none") 
dev.off()

pdf("Output pdfs/param intervals density.pdf", width =5, height=8.5)
ggplot(nbetas, aes(x=reorder(G.sp, -spid), y=X50.)) +
  geom_pointrange(aes(ymin=X5., ymax=X95., color=sign, alpha=signif90), shape=16) +
  scale_color_manual(values = c("red3", "medium blue")) +
  scale_alpha(range=c(0.2,1)) +
  xlab("") + ylab("Effect on Density") + 
  geom_hline(aes(yintercept=0), linetype="dashed") +
  coord_flip() +
  facet_wrap(~betas, ncol = 8, scales = "free_x", labeller = dens_labeller) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(size=8), axis.text.y = element_text(size=9, face="italic"))+
  theme(legend.position="none") 
dev.off()

pdf("Output pdfs/param intervals growth.pdf", width =11, height=8.5)
ggplot(gbetas, aes(x=reorder(G.sp, -spid), y=X50.)) +
  geom_pointrange(aes(ymin=X5., ymax=X95., color=sign, alpha=signif90), shape=16) +
  scale_color_manual(values = c("red3", "medium blue")) +
  scale_alpha(range=c(0.2,1)) +
  xlab("") + ylab("Effect on Growth") + 
  geom_hline(aes(yintercept=0), linetype="dashed") +
  coord_flip() +
  facet_wrap(~betas, ncol = 8, scales = "free_x", labeller = pred_labeller) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(size=8), axis.text.y = element_text(size=9, face="italic"))+
  theme(legend.position="none") 
dev.off()

pdf("Output pdfs/param intervals surv.pdf", width =11, height=8.5)
ggplot(sbetas, aes(x=reorder(G.sp, -spid), y=X50.)) +
  geom_pointrange(aes(ymin=X5., ymax=X95., color=sign, alpha=signif90), shape=16) +
  scale_color_manual(values = c("red3", "medium blue")) +
  scale_alpha(range=c(0.2,1)) +
  xlab("") + ylab("Effect on Survival") + 
  geom_hline(aes(yintercept=0), linetype="dashed") +
  coord_flip() +
  facet_wrap(~betas, ncol = 8, scales = "free_x", labeller = pred_labeller) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(size=8), axis.text.y = element_text(size=10, face="italic"))+
  theme(legend.position="none") 
dev.off()

pdf("Output pdfs/param intervals rec.pdf", width =11, height=8.5)
ggplot(rbetas, aes(x=reorder(G.sp, -spid), y=X50.)) +
  geom_pointrange(aes(ymin=X5., ymax=X95., color=sign, alpha=signif90), shape=16) +
  scale_color_manual(values = c("red3", "medium blue")) +
  scale_alpha(range=c(0.2,1)) +
  xlab("") + ylab("Effect on Recruitment") + 
  geom_hline(aes(yintercept=0), linetype="dashed") +
  coord_flip() +
  facet_wrap(~betas, ncol = 8, scales = "free_x", labeller = pred_labeller) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(axis.text.x = element_text(size=8), axis.text.y = element_text(size=10, face="italic"))+
  theme(legend.position="none") 
dev.off()

## regression coefficients with traits----
allbetas <- bind_rows(sbetas2, gbetas2, rbetas2)

## merge trait info
species_traits <- read.csv("Ch 1/data for stan/all_traits.csv")
allbetas <- merge(allbetas, species_traits, all.x=T)
allbetas$vr <- as.factor(allbetas$vr)

## betas by family----
pdf("Output pdfs/regression coeffs by family.pdf")
ggplot(allbetas, aes(x=family, y=X50._b_swb, fill=vr)) +
  geom_boxplot(position="dodge") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic() +
  ylab("Median parameter estimate") + xlab("Family") +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
  ggtitle("1st order SWB regression coefficient")

ggplot(allbetas, aes(x=family, y=X50._b2_swb, fill=vr)) +
  geom_boxplot(position="dodge") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic() +
  ylab("Median parameter estimate") + xlab("Family") +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
  ggtitle("2nd order SWB regression coefficient")

ggplot(allbetas, aes(x=family, y=X50._b_temp, fill=vr)) +
  geom_boxplot(position="dodge") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic() +
  ylab("Median parameter estimate") + xlab("Family") +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
  ggtitle("1st order temp regression coefficient")

ggplot(allbetas, aes(x=family, y=X50._b2_temp, fill=vr)) +
  geom_boxplot(position="dodge") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic() +
  ylab("Median parameter estimate") + xlab("Family") +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
  ggtitle("2nd order temp regression coefficient")

ggplot(allbetas, aes(x=family, y=X50._b_dia, fill=vr)) +
  geom_boxplot(position="dodge") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic() +
  ylab("Median parameter estimate") + xlab("Family") +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
  ggtitle("1st order diameter regression coefficient")

ggplot(allbetas, aes(x=family, y=X50._b2_dia, fill=vr)) +
  geom_boxplot(position="dodge") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic() +
  ylab("Median parameter estimate") + xlab("Family") +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
  ggtitle("2nd order diameter regression coefficient")

ggplot(allbetas, aes(x=family, y=X50._b_ba, fill=vr)) +
  geom_boxplot(position="dodge") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic() +
  ylab("Median parameter estimate") + xlab("Family") +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
  ggtitle("1st order mean(BA) regression coefficient")

ggplot(allbetas, aes(x=family, y=X50._b_adcount, fill=vr)) +
  geom_boxplot(position="dodge") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic() +
  ylab("Median parameter estimate") + xlab("Family") +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
  ggtitle("1st order adcount regression coefficient")

dev.off()

## betas by traits----
pdf("Output pdfs/b_SWB by traits.pdf")
ggplot(allbetas, aes(x=med_plot_swb, y=X50._b_swb)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE) +
  facet_wrap(~vr)

ggplot(allbetas, aes(x=med_plot_temp, y=X50._b_swb)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE) +
  facet_wrap(~vr)

ggplot(allbetas, aes(x=med_range_lat, y=X50._b_swb)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE) +
  facet_wrap(~vr)

ggplot(allbetas, aes(x=med_plot_elev, y=X50._b_swb)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE) +
  facet_wrap(~vr)

ggplot(allbetas, aes(x=wood_density, y=X50._b_swb)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE) +
  facet_wrap(~vr)

ggplot(allbetas, aes(x=Drought.tolerance, y=X50._b_swb)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE) +
  facet_wrap(~vr)

ggplot(allbetas, aes(x=Shade.tolerance, y=X50._b_swb)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE) +
  facet_wrap(~vr)

ggplot(allbetas, aes(x=Waterlogging.tolerance, y=X50._b_swb)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE) +
  facet_wrap(~vr)

dev.off()

