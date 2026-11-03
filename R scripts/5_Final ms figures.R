### Description: Final figures for FIA ms
### Author: Teresa Bohner
### Date modified: 2 Oct 2018

## packages----
library(tidyverse)
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggpmisc)
library(cowplot)

## data----
data.path <- "Processed Data/Analysis Ready Data/"

SWBs <- read.csv(paste0(data.path, "Call_mismatch_and_traits.csv"))


## Range overlap plot----
ggplot(SWBs, aes(overlap, abmis_p, color=vr)) +
  geom_point() +
  geom_smooth(method='lm', se = FALSE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text_npc', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, size = 4)

ggplot(SWBs, aes(overlap, mismatch_p, color=vr)) +
  geom_point() +
  geom_smooth(method='lm', se = FALSE) +
  geom_hline(yintercept=0, linetype="dashed") +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text_npc', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, size = 4)

ggplot(SWBs, aes(overlap_bien, abmis_p, color=vr)) +
  geom_point() +
  geom_smooth(method='lm', se = FALSE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text_npc', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, size = 4)

ggplot(SWBs, aes(overlap_bien, mismatch_p, color=vr)) +
  geom_point() +
  geom_smooth(method='lm', se = FALSE) +
  geom_hline(yintercept=0, linetype="dashed") +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text_npc', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, size = 4) +
  theme_classic()

p1 <- ggplot(SWBs, aes(overlap_max, abmis_p, color=vr)) +
  geom_point() +
  geom_smooth(method='lm', se = FALSE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text_npc', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y= "top", size = 4) +
  theme_classic() +
  ylab("Mismatch magnitude") + xlab("% range overlap, larger value") +
  # scale_color_discrete(name="",
  #                      breaks=c("0abun", "1surv", "2growth", "3rec"),
  #                      labels=c("Density", "Survival", "Growth", "Recruitment"))
  theme(legend.position = "none")
p1

p2 <- ggplot(SWBs, aes(overlap_max, mismatch_p, color=vr), labeller=vr_labeller) +
  geom_point() +
  geom_smooth(method='lm', se = FALSE,) +
  geom_hline(yintercept=0, linetype="dashed") +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text_npc', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, size = 4) +
  theme_classic() +
  ylab("Mismatch") + xlab("% range overlap, larger value") +
  # scale_color_discrete(name="",
  #                     breaks=c("0abun", "1surv", "2growth", "3rec"),
  #                     labels=c("Density", "Survival", "Growth", "Recruitment"))
  theme(legend.position = "none")
p2

forlegend <- ggplot(SWBs, aes(overlap_max, mismatch_p, color=vr)) +
  geom_line(aes(color=vr), size=1) +
  theme_classic() +
  scale_color_discrete(name="",
                       breaks=c("0abun", "1surv", "2growth", "3rec"),
                       labels=c("Density", "Survival", "Growth", "Recruitment"))
legend <- get_legend(forlegend)
grid.arrange(legend)

ggdraw(plot_grid(plot_grid(p1, p2, nrow=1, align='h'),
                 plot_grid(legend, nrow=1),
                 rel_widths=c(1, 0.13)))

pdf("Output pdfs/jan19 abundance figures/mismatch by range.pdf", width = 9, height = 4)
ggdraw(plot_grid(plot_grid(p1, p2, nrow=1, align='h'),
                 plot_grid(legend, nrow=1),
                 rel_widths=c(1, 0.15)))
dev.off()


Ntab <- read.csv(paste0(data.path,"tableinfo.csv"))

Ntablong <- Ntab %>%
  dplyr::select(c("gensp","nplots_p", "ntrees_g", "ntrees_s", "npres_r")) %>%
  rename("0abun"="nplots_p", "2growth"="ntrees_g", "1surv"="ntrees_s", "3rec"="npres_r") %>%
  gather("0abun", "2growth", "1surv", "3rec", key="vr", value="N")

SWBsplus <- merge(SWBs[,1:82], Ntablong, by=c("gensp", "vr"), all=T)
Nmismatchtest <- SWBsplus %>%
  data.table()%>%
  data.table::dcast(gensp ~vr,
                    value.var=c("mismatch_p",
                                "abmis_p",
                                "N", "r2_all"))
##N plots----
p1 <- ggplot(Nmismatchtest, aes(N_0abun, abmis_p_0abun)) +
  geom_point() +
  geom_smooth(method="lm") + ylab("Density mismatch magnitude") + xlab("N records for density model") +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.85, size = 4) +
  theme_classic() 
p1

p2 <- ggplot(Nmismatchtest, aes(N_1surv, abmis_p_1surv)) +
  geom_point() +
  geom_smooth(method="lm") + ylab("Survival mismatch magnitude") + xlab("N records for survival model") +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.85, size = 4) +
  theme_classic()
p2

p3 <-ggplot(Nmismatchtest, aes(N_2growth, abmis_p_2growth)) +
  geom_point() +
  geom_smooth(method="lm") + ylab("Growth mismatch magnitude") + xlab("N records for growth model") +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.85, size = 4) +
  theme_classic()
p3

p4 <- ggplot(Nmismatchtest, aes(N_3rec, abmis_p_3rec)) +
  geom_point() +
  geom_smooth(method="lm") + ylab("Recruitment mismatch magnitude") + xlab("N records for recruitment model") +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.85, size = 4) +
  theme_classic()
p4


p5 <- ggplot(Nmismatchtest, aes(N_0abun, mismatch_p_0abun)) +
  geom_point() +
  geom_smooth(method="lm") + ylab("Density mismatch") + xlab("N records for density model") +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.85, size = 4) +
  theme_classic() +
  geom_hline(yintercept=0, linetype="dashed")
p5

p6 <- ggplot(Nmismatchtest, aes(N_1surv, mismatch_p_1surv)) +
  geom_point() +
  geom_smooth(method="lm") + ylab("Survival mismatch") + xlab("N records for survival model") +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.85, size = 4) +
  theme_classic() +
  geom_hline(yintercept=0, linetype="dashed")
p6

p7 <-ggplot(Nmismatchtest, aes(N_2growth, mismatch_p_2growth)) +
  geom_point() +
  geom_smooth(method="lm") + ylab("Growth mismatch") + xlab("N records for growth model") +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.85, size = 4) +
  theme_classic() +
  geom_hline(yintercept=0, linetype="dashed")
p7

p8 <- ggplot(Nmismatchtest, aes(N_3rec, mismatch_p_3rec)) +
  geom_point() +
  geom_smooth(method="lm") + ylab("Recruitment mismatch") + xlab("N records for recruitment model") +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.85, size = 4) +
  theme_classic() +
  geom_hline(yintercept=0, linetype="dashed")
p8


pdf("Output pdfs/jan19 abundance figures/mismatch by N.pdf", width = 8, height = 8)
grid.arrange(p1, p2, p3, p4, nrow=2)
grid.arrange(p5, p6, p7, p8, nrow=2)
dev.off()

##R2 plots----
p1 <- ggplot(Nmismatchtest, aes(r2_all_0abun, abmis_p_0abun)) +
  geom_point() +
  geom_smooth(method="lm") + ylab("Density mismatch magnitude") + xlab(expression(paste("Density model fit, ", italic("R" ^"2"), sep=""))) +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.85, size = 4) +
  theme_classic() 
p1

p2 <- ggplot(Nmismatchtest, aes(r2_all_1surv, abmis_p_1surv)) +
  geom_point() +
  geom_smooth(method="lm") + ylab("Survival mismatch magnitude") + xlab(expression(paste("Survival model fit, ", italic("R" ^"2"), sep=""))) +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.85, size = 4) +
  theme_classic()
p2

p3 <-ggplot(Nmismatchtest, aes(r2_all_2growth, abmis_p_2growth)) +
  geom_point() +
  geom_smooth(method="lm") + ylab("Growth mismatch magnitude") + xlab(expression(paste("Growth model fit, ", italic("R" ^"2"), sep=""))) +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.85, size = 4) +
  theme_classic()
p3

p4 <- ggplot(Nmismatchtest, aes(r2_all_3rec, abmis_p_3rec)) +
  geom_point() +
  geom_smooth(method="lm") + ylab("Recruitment mismatch magnitude") + xlab(expression(paste("Recruitment model fit, ", italic("R" ^"2"), sep=""))) +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.85, size = 4) +
  theme_classic()
p4


p5 <- ggplot(Nmismatchtest, aes(r2_all_0abun, mismatch_p_0abun)) +
  geom_point() +
  geom_smooth(method="lm") + ylab("Density mismatch") + xlab(expression(paste("Density model fit, ", italic("R" ^"2"), sep=""))) +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.85, size = 4) +
  theme_classic() +
  geom_hline(yintercept=0, linetype="dashed")
p5

p6 <- ggplot(Nmismatchtest, aes(r2_all_1surv, mismatch_p_1surv)) +
  geom_point() +
  geom_smooth(method="lm") + ylab("Survival mismatch") + xlab(expression(paste("Survival model fit, ", italic("R" ^"2"), sep=""))) +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.85, size = 4) +
  theme_classic() +
  geom_hline(yintercept=0, linetype="dashed")
p6

p7 <-ggplot(Nmismatchtest, aes(r2_all_2growth, mismatch_p_2growth)) +
  geom_point() +
  geom_smooth(method="lm") + ylab("Growth mismatch") + xlab(expression(paste("Growth model fit, ", italic("R" ^"2"), sep=""))) +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.85, size = 4) +
  theme_classic() +
  geom_hline(yintercept=0, linetype="dashed")
p7

p8 <- ggplot(Nmismatchtest, aes(r2_all_3rec, mismatch_p_3rec)) +
  geom_point() +
  geom_smooth(method="lm") + ylab("Recruitment mismatch") + xlab(expression(paste("Recruitment model fit, ", italic("R" ^"2"), sep=""))) +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.85, size = 4) +
  theme_classic() +
  geom_hline(yintercept=0, linetype="dashed")
p8


pdf("Output pdfs/jan19 abundance figures/mismatch by R2.pdf", width = 8, height = 8)
grid.arrange(p1, p2, p3, p4, nrow=2)
grid.arrange(p5, p6, p7, p8, nrow=2)
dev.off()


## Fig 1 is conceptual. 

# wide format for bivariate plots
SWBwide <- SWBs %>%
  data.table()%>%
  data.table::dcast(gensp + family + gymno2~vr,
                    value.var=c("mismatch_p",
                                "abmis_p",
                                "cost_p", 
                                "scaled_cost_p"))
hist(SWBwide$abmis_p_0abun)
hist(SWBwide$abmis_p_1surv)
hist(SWBwide$abmis_p_2growth)
hist(SWBwide$abmis_p_3rec)

mean(SWBwide$abmis_p_0abun, na.rm=T)
mean(SWBwide$abmis_p_1surv, na.rm=T)
mean(SWBwide$abmis_p_2growth, na.rm=T)
mean(SWBwide$abmis_p_3rec, na.rm=T)

wilcox.test(x = SWBwide$abmis_p_0abun)
wilcox.test(x = SWBwide$abmis_p_1surv)
wilcox.test(x = SWBwide$abmis_p_2growth)
wilcox.test(x = SWBwide$abmis_p_3rec)

wilcox.test(x = SWBwide$mismatch_p_0abun)
wilcox.test(x = SWBwide$mismatch_p_1surv)
wilcox.test(x = SWBwide$mismatch_p_2growth)
wilcox.test(x = SWBwide$mismatch_p_3rec)
i=1

SWBs$lineint <- ifelse(SWBs$vr=="0abun", SWBs$mismatch_p, NA)
pdf("Output pdfs/jan19 abundance figures/within species comparison_simple.pdf")
for(i in 1:60) {
  print(
    ggplot(dplyr::filter(SWBs, spid==i), aes(x=mismatch_p, y=scaled_cost_p)) +
      geom_vline(aes(xintercept = lineint, color=vr), linetype="dashed") +
      geom_linerange(aes(ymin=0, ymax=scaled_cost_p, color=vr), size=4, alpha=0.75) +
      theme_classic() +
      xlab("Mismatch") + ylab("Scaled cost") +
      ggtitle(SWBs$gensp[SWBs$spid==i])
  )
}
dev.off()


  ## Fig S6----
p0 <-ggplot(SWBwide, aes(mismatch_p_0abun)) +
  geom_histogram(aes(y = ..count..),
                 binwidth=6,
                 colour="black", fill="dark gray") +
  geom_density(aes(y = ..density..*(60*6)), alpha=.4, fill="gray") +
  scale_x_continuous(name ="Density Mismatch", 
                     breaks = seq(-60,60, by=20), 
                     limits=c(-max(SWBwide$mismatch_p_0abun, na.rm=T), max(SWBwide$mismatch_p_0abun, na.rm=T))) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  ylab("Number of species") +
  ggtitle("a)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p0
p1 <-ggplot(SWBwide, aes(mismatch_p_1surv)) +
  geom_histogram(aes(y = ..count..),
                 binwidth=6,
                 colour="black", fill="dark gray") +
  geom_density(aes(y = ..density..*(60*6)), alpha=.4, fill="gray") +
  scale_x_continuous(name ="Survival Mismatch", 
                     breaks = seq(-60,60, by=20), 
                     limits=c(-max(SWBwide$abmis_p_1surv, na.rm=T), max(SWBwide$abmis_p_1surv, na.rm=T))) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  ggtitle("b)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))+
  theme(axis.title.y=element_blank())
p1
p2 <-ggplot(SWBwide, aes(mismatch_p_2growth)) +
  geom_histogram(aes(y=..count..),      # Histogram with density instead of count on y-axis
                 binwidth=6,
                 colour="black", fill="dark gray") +
  geom_density(aes(y=..density..*(50*6)), alpha=.4, fill="gray") +
  scale_x_continuous(name ="Growth Mismatch", 
                     breaks = seq(-60,60, by=20), 
                     limits=c(-max(SWBwide$abmis_p_1surv, na.rm=T), max(SWBwide$abmis_p_1surv, na.rm=T))) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  ggtitle("c)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold")) +
  theme(axis.title.y=element_blank())
p2
p3 <-ggplot(SWBwide, aes(mismatch_p_3rec)) +
  geom_histogram(aes(y=..count..),      # Histogram with density instead of count on y-axis
                 binwidth=6,
                 colour="black", fill="dark gray") +
  geom_density(aes(y=..density..*(55*6)), alpha=.4, fill="gray") +
  scale_x_continuous(name ="Recruitment Mismatch", 
                     breaks = seq(-60,60, by=20), 
                     limits=c(-max(SWBwide$abmis_p_1surv, na.rm=T), max(SWBwide$abmis_p_1surv, na.rm=T))) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  ggtitle("d)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold")) +
  theme(axis.title.y=element_blank())
p3

pdf("Output pdfs/jan19 abundance figures/Fig S6_abundance_simple.pdf", width = 12, height = 4)
grid.arrange(p0, p1, p2, p3, nrow=1)
dev.off()


## Fig 2----
## dens ~ Vr raw mismatch scatterplots

corr1 <- cor.test(SWBwide$mismatch_p_1surv, SWBwide$mismatch_p_0abun, method = 'spearman')
corr2 <- cor.test(SWBwide$mismatch_p_2growth, SWBwide$mismatch_p_0abun, method = 'spearman')
corr3 <- cor.test(SWBwide$mismatch_p_3rec, SWBwide$mismatch_p_0abun, method = 'spearman')

lm1 <- lm(SWBwide$mismatch_p_0abun ~ SWBwide$mismatch_p_1surv)
lm2 <- lm(SWBwide$mismatch_p_0abun ~ SWBwide$mismatch_p_2growth)
lm3 <- lm(SWBwide$mismatch_p_0abun ~ SWBwide$mismatch_p_3rec)

p1 <-ggplot(SWBwide, aes(mismatch_p_1surv, mismatch_p_0abun)) +
  geom_point(size=2.5, alpha=0.5) +
  # geom_smooth(method="lm", se=FALSE) +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, size = 4) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  scale_x_continuous(name ="Survival mismatch", 
                     breaks = seq(-60,60, by=20), 
                     limits=c(-max(SWBwide$abmis_p_1surv, na.rm=T), max(SWBwide$abmis_p_1surv, na.rm=T))) +
  scale_y_continuous(name ="Density mismatch", 
                     breaks = seq(-80,80, by=20), 
                     limits=c(-max(SWBwide$abmis_p_0abun, na.rm=T), max(SWBwide$abmis_p_0abun, na.rm=T))) +
  ggtitle("a)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p1
if(summary(lm1)$coefficients[2,4]<=0.05){
  p1 <- p1 + geom_abline(intercept = coef(lm1)[1], slope = coef(lm1)[2], size=1)
}   else{
  p1 <- p1 + geom_abline(intercept = coef(lm1)[1], slope = coef(lm1)[2], size=1, linetype="dashed", color="gray")
}
p1

p2 <- ggplot(SWBwide, aes(mismatch_p_2growth, mismatch_p_0abun)) +
  geom_point(size=2.5, alpha = 0.5) +
  # geom_smooth(method="lm", se=FALSE) +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, size = 4) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  scale_x_continuous(name ="Growth mismatch", 
                     breaks = seq(-40,40, by=20), 
                     limits=c(-max(SWBwide$abmis_p_2growth, na.rm=T), max(SWBwide$abmis_p_2growth, na.rm=T))) +
  scale_y_continuous(name ="Density mismatch", 
                     breaks = seq(-80,80, by=20), 
                     limits=c(-max(SWBwide$abmis_p_0abun, na.rm=T), max(SWBwide$abmis_p_0abun, na.rm=T))) +
  ggtitle("b)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold")) 

if(summary(lm2)$coefficients[2,4]<=0.05){
  p2 <- p2 + geom_abline(intercept = coef(lm2)[1], slope = coef(lm2)[2], size=1)
}   else{
  p2 <- p2 + geom_abline(intercept = coef(lm2)[1], slope = coef(lm2)[2], size=1, linetype="dashed", color="gray")
}
p2

p3 <- ggplot(SWBwide, aes(mismatch_p_3rec, mismatch_p_0abun)) +
  geom_point(size=2.5, alpha = 0.5) +
  # geom_jitter(width = 0.5, height = 0.5, size=2.5, alpha=0.5) +
  # geom_smooth(method="lm", se=FALSE) +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, size = 4) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  scale_x_continuous(name ="Recruitment mismatch", 
                     breaks = seq(-40,40, by=20), 
                     limits=c(-max(SWBwide$abmis_p_3rec, na.rm=T), max(SWBwide$abmis_p_3rec, na.rm=T))) +
  scale_y_continuous(name ="Density mismatch", 
                     breaks = seq(-80,80, by=20), 
                     limits=c(-max(SWBwide$abmis_p_0abun, na.rm=T), max(SWBwide$abmis_p_0abun, na.rm=T))) +
  ggtitle("c)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))

if(summary(lm3)$coefficients[2,4]<=0.05){
  p3 <- p3 + geom_abline(intercept = coef(lm3)[1], slope = coef(lm3)[2], size=1)
}   else{
  p3 <- p3 + geom_abline(intercept = coef(lm3)[1], slope = coef(lm3)[2], size=1, linetype="dashed", color="gray")
}

p3

pdf("Output pdfs/jan19 abundance figures/Fig 2_abundance_simple.pdf", width = 9, height = 3)
grid.arrange(p1, p2, p3, nrow=1)
dev.off()

## Fig S7----
## Vr ~ Vr raw mismatch scatterplots

corr1 <- cor.test(SWBwide$mismatch_p_1surv, SWBwide$mismatch_p_2growth, method = 'spearman')
corr2 <- cor.test(SWBwide$mismatch_p_2growth, SWBwide$mismatch_p_3rec, method = 'spearman')
corr3 <- cor.test(SWBwide$mismatch_p_3rec, SWBwide$mismatch_p_1surv, method = 'spearman')

lm1 <- lm(SWBwide$mismatch_p_2growth ~ SWBwide$mismatch_p_1surv)
lm2 <- lm(SWBwide$mismatch_p_3rec ~ SWBwide$mismatch_p_2growth)
lm3 <- lm(SWBwide$mismatch_p_1surv ~ SWBwide$mismatch_p_3rec)

p4 <-ggplot(SWBwide, aes(mismatch_p_1surv, mismatch_p_2growth)) +
  geom_point(size=2.5, alpha=0.5) +
  # geom_smooth(method="lm", se=FALSE) +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, size = 4) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  scale_x_continuous(name ="Survival mismatch", 
                     breaks = seq(-60,60, by=20), 
                     limits=c(-max(SWBwide$abmis_p_1surv, na.rm=T), max(SWBwide$abmis_p_1surv, na.rm=T))) +
  scale_y_continuous(name ="Growth mismatch", 
                     breaks = seq(-40,40, by=20), 
                     limits=c(-max(SWBwide$abmis_p_2growth, na.rm=T), max(SWBwide$abmis_p_2growth, na.rm=T))) +
  ggtitle("d)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))

if(summary(lm1)$coefficients[2,4]<=0.05){
  p4 <- p4 + geom_abline(intercept = coef(lm1)[1], slope = coef(lm1)[2], size=1)
}   else{
  p4 <- p4 + geom_abline(intercept = coef(lm1)[1], slope = coef(lm1)[2], size=1, linetype="dashed", color="gray")
}
p4

p5 <- ggplot(SWBwide, aes(mismatch_p_2growth, mismatch_p_3rec)) +
  geom_point(size=2.5, alpha = 0.5) +
  # geom_smooth(method="lm", se=FALSE) +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, size = 4) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  scale_x_continuous(name ="Growth mismatch", 
                     breaks = seq(-40,40, by=20), 
                     limits=c(-max(SWBwide$abmis_p_2growth, na.rm=T), max(SWBwide$abmis_p_2growth, na.rm=T))) +
  scale_y_continuous(name ="Recruitment mismatch", 
                     breaks = seq(-40,40, by=20), 
                     limits=c(-max(SWBwide$abmis_p_3rec, na.rm=T), max(SWBwide$abmis_p_3rec, na.rm=T))) +
  ggtitle("e)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))

if(summary(lm2)$coefficients[2,4]<=0.05){
  p5 <- p5 + geom_abline(intercept = coef(lm2)[1], slope = coef(lm2)[2], size=1)
}   else{
  p5 <- p5 + geom_abline(intercept = coef(lm2)[1], slope = coef(lm2)[2], size=1, linetype="dashed", color="gray")
}
p5

p6 <- ggplot(SWBwide, aes(mismatch_p_3rec, mismatch_p_1surv)) +
  geom_point(size=2.5, alpha = 0.5) +
  # geom_jitter(width = 0.5, height = 0.5, size=2.5, alpha=0.5) +
  # geom_smooth(method="lm", se=FALSE) +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, size = 4) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  scale_x_continuous(name ="Recruitment mismatch", 
                     breaks = seq(-40,40, by=20), 
                     limits=c(-max(SWBwide$abmis_p_3rec, na.rm=T), max(SWBwide$abmis_p_3rec, na.rm=T))) +
  scale_y_continuous(name ="Survival mismatch", 
                     breaks = seq(-60,60, by=20), 
                     limits=c(-max(SWBwide$abmis_p_1surv, na.rm=T), max(SWBwide$abmis_p_1surv, na.rm=T))) +
  ggtitle("f)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))

if(summary(lm3)$coefficients[2,4]<=0.05){
  p6 <- p6 + geom_abline(intercept = coef(lm3)[1], slope = coef(lm3)[2], size=1)
}   else{
  p6 <- p6 + geom_abline(intercept = coef(lm3)[1], slope = coef(lm3)[2], size=1, linetype="dashed", color="gray")
}
p6

pdf("Output pdfs/jan19 abundance figures/all performance corr.pdf", width = 9, height = 6)
grid.arrange(p1, p2, p3, p4, p5, p6, nrow=2)
dev.off()

## Fig 3----
## Mismatch ~ Cost raw value scatterplots

lm4 <- lm(cost_p_1surv ~ mismatch_p_1surv, data = SWBwide[which(SWBwide$mismatch_p_1surv<0),])
lm5 <- lm(cost_p_1surv ~ mismatch_p_1surv, data = SWBwide[which(SWBwide$mismatch_p_1surv>0),])
lm6 <- lm(cost_p_2growth ~ mismatch_p_2growth, data = SWBwide[which(SWBwide$mismatch_p_2growth<0),])
lm7 <- lm(cost_p_2growth ~ mismatch_p_2growth, data = SWBwide[which(SWBwide$mismatch_p_2growth>0),])
lm8 <- lm(cost_p_3rec ~ mismatch_p_3rec, data = SWBwide[which(SWBwide$mismatch_p_3rec<0),])
lm9 <- lm(cost_p_3rec ~ mismatch_p_3rec, data = SWBwide[which(SWBwide$mismatch_p_3rec>0),])
lm10 <- lm(cost_p_0abun ~ mismatch_p_0abun, data = SWBwide[which(SWBwide$mismatch_p_0abun<0),])
lm11 <- lm(cost_p_0abun ~ mismatch_p_0abun, data = SWBwide[which(SWBwide$mismatch_p_0abun>0),])

corr4 <- cor.test(SWBwide$mismatch_p_1surv[which(SWBwide$mismatch_p_1surv<0)], SWBwide$cost_p_1surv[which(SWBwide$mismatch_p_1surv<0)], method="spearman")
corr5 <- cor.test(SWBwide$mismatch_p_1surv[which(SWBwide$mismatch_p_1surv>0)], SWBwide$cost_p_1surv[which(SWBwide$mismatch_p_1surv>0)], method="spearman")
corr6 <- cor.test(SWBwide$mismatch_p_2growth[which(SWBwide$mismatch_p_2growth<0)], SWBwide$cost_p_2growth[which(SWBwide$mismatch_p_2growth<0)], method="spearman")
corr7 <- cor.test(SWBwide$mismatch_p_2growth[which(SWBwide$mismatch_p_2growth>0)], SWBwide$cost_p_2growth[which(SWBwide$mismatch_p_2growth>0)], method="spearman")
corr8 <- cor.test(SWBwide$mismatch_p_3rec[which(SWBwide$mismatch_p_3rec<0)], SWBwide$cost_p_3rec[which(SWBwide$mismatch_p_3rec<0)], method="spearman")
corr9 <- cor.test(SWBwide$mismatch_p_3rec[which(SWBwide$mismatch_p_3rec>0)], SWBwide$cost_p_3rec[which(SWBwide$mismatch_p_3rec>0)], method="spearman")
corr10 <- cor.test(SWBwide$mismatch_p_0abun[which(SWBwide$mismatch_p_0abun<0)], SWBwide$cost_p_0abun[which(SWBwide$mismatch_p_0abun<0)], method="spearman")
corr11 <- cor.test(SWBwide$mismatch_p_0abun[which(SWBwide$mismatch_p_0abun>0)], SWBwide$cost_p_0abun[which(SWBwide$mismatch_p_0abun>0)], method="spearman")

p0 <- ggplot(SWBwide, aes(mismatch_p_0abun, cost_p_0abun)) +
  geom_point(size=2.5, alpha=0.5) +
  # geom_smooth(data = SWBwide[which(SWBwide$mismatch_p_1surv<0),], method = "lm", se=F) +
  # geom_smooth(data = SWBwide[which(SWBwide$mismatch_p_1surv>0),], method = "lm", se=F) +
  stat_poly_eq(data = SWBwide[which(SWBwide$mismatch_p_0abun<0),],
               aes(label = ..rr.label..), label.y.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(data = SWBwide[which(SWBwide$mismatch_p_0abun<0),],
                  method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, size = 4) +
  stat_poly_eq(data = SWBwide[which(SWBwide$mismatch_p_0abun>0),],
               aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.9, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(data = SWBwide[which(SWBwide$mismatch_p_0abun>0),],
                  method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.9, size = 4) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  scale_x_continuous(name ="Mismatch", 
                     breaks = seq(-80,80, by=20), 
                     limits=c(-max(SWBwide$abmis_p_0abun, na.rm=T), max(SWBwide$abmis_p_0abun, na.rm=T))
  ) +
  scale_y_continuous(name ="Density cost, trees per plot", 
                     breaks = seq(0, 60, by=10), 
                     limits=c(-2, max(SWBwide$cost_p_0abun, na.rm=T))) +
  ggtitle("a)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p0
if(summary(lm10)$coefficients[2,4]<=0.05){
  p0 <- p0 + geom_segment(aes(x=0, xend=(-max(SWBwide$abmis_p_0abun, na.rm=T)), 
                   y=coef(lm10)[1], yend=(coef(lm10)[1]+coef(lm10)[2]*(-max(SWBwide$abmis_p_0abun, na.rm=T)))), 
               size=1)

}   else{
  p0 <- p0 + geom_segment(aes(x=0, xend=(-max(SWBwide$abmis_p_0abun, na.rm=T)), 
                              y=coef(lm10)[1], yend=coef(lm10)[1]+coef(lm10)[2]*(-max(SWBwide$abmis_p_0abun, na.rm=T))), 
                          size=1, linetype ="dashed", color="gray") 
}
p0
if(summary(lm11)$coefficients[2,4]<=0.05){
  p0 <- p0 + geom_segment(aes(x=0, xend=max(SWBwide$abmis_p_0abun, na.rm=T), 
                                y=coef(lm11)[1], yend=coef(lm11)[1]+coef(lm11)[2]*max(SWBwide$abmis_p_0abun, na.rm=T)), 
                            size=1) 
}   else{
  p0 <- p0 + geom_segment(aes(x=0, xend=max(SWBwide$abmis_p_0abun, na.rm=T), 
                                y=coef(lm11)[1], yend=coef(lm11)[1]+coef(lm11)[2]*max(SWBwide$abmis_p_0abun, na.rm=T)), 
                            size=1, linetype ="dashed", color="gray")
}
p0

p1 <- ggplot(SWBwide, aes(mismatch_p_1surv, cost_p_1surv)) +
  geom_point(size=2.5, alpha=0.5) +
  # geom_smooth(data = SWBwide[which(SWBwide$mismatch_p_1surv<0),], method = "lm", se=F) +
  # geom_smooth(data = SWBwide[which(SWBwide$mismatch_p_1surv>0),], method = "lm", se=F) +
  stat_poly_eq(data = SWBwide[which(SWBwide$mismatch_p_1surv<0),],
               aes(label = ..rr.label..), label.y.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(data = SWBwide[which(SWBwide$mismatch_p_1surv<0),],
                  method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, size = 4) +
  stat_poly_eq(data = SWBwide[which(SWBwide$mismatch_p_1surv>0),],
               aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.9, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(data = SWBwide[which(SWBwide$mismatch_p_1surv>0),],
                  method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.9, size = 4) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  scale_x_continuous(name ="Mismatch", 
                     breaks = seq(-60,60, by=20), 
                     limits=c(-max(SWBwide$abmis_p_1surv, na.rm=T), max(SWBwide$abmis_p_1surv, na.rm=T))
                     ) + #paste("Survival cost (", expression(italic(p)), "(Survival))", sep = ""),
  scale_y_continuous(name = expression(paste("Survival cost, ",italic("p"), "(Survival)", sep="")),
                     breaks = seq(0, 0.4, by=0.1), 
                     limits=c(-0.025, max(SWBwide$cost_p_1surv, na.rm=T))) +
  ggtitle("b)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p1
if(summary(lm4)$coefficients[2,4]<=0.05){
  p1 <- p1 + geom_segment(aes(x=0, xend=(-max(SWBwide$abmis_p_1surv, na.rm=T)), 
                   y=coef(lm4)[1], yend=coef(lm4)[1]+coef(lm4)[2]*(-max(SWBwide$abmis_p_1surv, na.rm=T))), 
               size=1)
}   else{
  p1 <- p1 + geom_segment(aes(x=0, xend=(-max(SWBwide$abmis_p_1surv, na.rm=T)), 
                              y=coef(lm4)[1], yend=coef(lm4)[1]+coef(lm4)[2]*(-max(SWBwide$abmis_p_1surv, na.rm=T))), 
                          size=1, linetype ="dashed", color="gray")

}
p1

if(summary(lm5)$coefficients[2,4]<=0.05){
  p1 <- p1 + geom_segment(aes(x=0, xend=max(SWBwide$abmis_p_1surv, na.rm=T), 
                              y=coef(lm5)[1], yend=coef(lm5)[1]+coef(lm5)[2]*max(SWBwide$abmis_p_1surv, na.rm=T)), 
                          size=1)
}   else{
  p1 <- p1 + geom_segment(aes(x=0, xend=max(SWBwide$abmis_p_1surv, na.rm=T), 
                              y=coef(lm5)[1], yend=coef(lm5)[1]+coef(lm5)[2]*max(SWBwide$abmis_p_1surv, na.rm=T)), 
                          size=1, linetype ="dashed", color="gray")
  
}
p1

p2 <- ggplot(SWBwide, aes(mismatch_p_2growth, cost_p_2growth)) +
  geom_point(size=2.5, alpha=0.5) +
  # geom_smooth(data = SWBwide[which(SWBwide$mismatch_p_2growth<0),], method = "lm", se=F) +
  # geom_smooth(data = SWBwide[which(SWBwide$mismatch_p_2growth>0),], method = "lm", se=F) +
  stat_poly_eq(data = SWBwide[which(SWBwide$mismatch_p_2growth<0),],
               aes(label = ..rr.label..), label.y.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(data = SWBwide[which(SWBwide$mismatch_p_2growth<0),],
                  method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, size = 4) +
  stat_poly_eq(data = SWBwide[which(SWBwide$mismatch_p_2growth>0),],
               aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.9, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(data = SWBwide[which(SWBwide$mismatch_p_2growth>0),],
                  method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.9, size = 4) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  scale_x_continuous(name ="Mismatch", 
                     breaks = seq(-40,40, by=20), 
                     limits=c(-max(SWBwide$abmis_p_2growth, na.rm=T), max(SWBwide$abmis_p_2growth, na.rm=T))) +
  scale_y_continuous(name = expression(paste("Growth cost, cm" ^ "2", sep="")),
                     breaks = seq(0, 50, by=10), 
                     limits=c(-4, max(SWBwide$cost_p_2growth, na.rm=T))) +
  ggtitle("c)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p2
if(summary(lm6)$coefficients[2,4]<=0.05){
  p2 <- p2 + geom_segment(aes(x=0, xend=(-max(SWBwide$abmis_p_2growth, na.rm=T)), 
                              y=coef(lm6)[1], yend=coef(lm6)[1]+coef(lm6)[2]*(-max(SWBwide$abmis_p_2growth, na.rm=T))), 
                          size=1) 
}   else{
  p2 <- p2 + geom_segment(aes(x=0, xend=(-max(SWBwide$abmis_p_2growth, na.rm=T)), 
                              y=coef(lm6)[1], yend=coef(lm6)[1]+coef(lm6)[2]*(-max(SWBwide$abmis_p_2growth, na.rm=T))), 
                          size=1, linetype ="dashed", color="gray")
  
}
p2
if(summary(lm7)$coefficients[2,4]<=0.05){
  p2 <- p2 + geom_segment(aes(x=0, xend=max(SWBwide$abmis_p_2growth, na.rm=T), 
                              y=coef(lm7)[1], yend=coef(lm7)[1]+coef(lm7)[2]*max(SWBwide$abmis_p_2growth, na.rm=T)), 
                          size=1)
}   else{
  p2 <- p2 + geom_segment(aes(x=0, xend=max(SWBwide$abmis_p_2growth, na.rm=T), 
                              y=coef(lm7)[1], yend=coef(lm7)[1]+coef(lm7)[2]*max(SWBwide$abmis_p_2growth, na.rm=T)), 
                          size=1, linetype ="dashed", color="gray")
  
}
p2

p3 <- ggplot(SWBwide, aes(mismatch_p_3rec, cost_p_3rec)) +
  geom_point(size=2.5, alpha=0.5) +
  # geom_smooth(data = SWBwide[which(SWBwide$mismatch_p_3rec<0),], method = "lm", se=F) +
  # geom_smooth(data = SWBwide[which(SWBwide$mismatch_p_3rec>0),], method = "lm", se=F) +
  stat_poly_eq(data = SWBwide[which(SWBwide$mismatch_p_3rec<0),],
               aes(label = ..rr.label..), label.y.npc = 0.85, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(data = SWBwide[which(SWBwide$mismatch_p_3rec<0),],
                  method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, size = 4) +
  stat_poly_eq(data = SWBwide[which(SWBwide$mismatch_p_3rec>0),],
               aes(label = ..rr.label..), label.y.npc = 0.85, label.x.npc = 0.9, size = 4, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(data = SWBwide[which(SWBwide$mismatch_p_3rec>0),],
                  method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1, label.x.npc = 0.9, size = 4) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  scale_x_continuous(name ="Mismatch", 
                     breaks = seq(-40,40, by=20), 
                     limits=c(-max(SWBwide$abmis_p_3rec, na.rm=T), max(SWBwide$abmis_p_3rec, na.rm=T))) +
  scale_y_continuous(name ="Recruitment cost, seedlings per plot", 
                     breaks = seq(0, 10, by=2), 
                     limits=c(-0.75, max(SWBwide$cost_p_3rec, na.rm=T))) +
  ggtitle("d)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p3  
if(summary(lm8)$coefficients[2,4]<=0.05){
  p3 <- p3 + geom_segment(aes(x=0, xend=(-max(SWBwide$abmis_p_3rec, na.rm=T)), 
                              y=coef(lm8)[1], yend=coef(lm8)[1]+coef(lm8)[2]*(-max(SWBwide$abmis_p_3rec, na.rm=T))), 
                          size=1) 
}   else{
  p3 <- p3 + geom_segment(aes(x=0, xend=(-max(SWBwide$abmis_p_3rec, na.rm=T)), 
                              y=coef(lm8)[1], yend=coef(lm8)[1]+coef(lm8)[2]*(-max(SWBwide$abmis_p_3rec, na.rm=T))), 
                          size=1, linetype ="dashed", color="gray")
}
p3

if(summary(lm9)$coefficients[2,4]<=0.05){
  p3 <- p3 + geom_segment(aes(x=0, xend=max(SWBwide$abmis_p_3rec, na.rm=T), 
                     y=coef(lm9)[1], yend=coef(lm9)[1]+coef(lm9)[2]*max(SWBwide$abmis_p_3rec, na.rm=T)), 
                 size=1) 
}   else{
  p3 <- p3 + geom_segment(aes(x=0, xend=max(SWBwide$abmis_p_3rec, na.rm=T), 
                     y=coef(lm9)[1], yend=coef(lm9)[1]+coef(lm9)[2]*max(SWBwide$abmis_p_3rec, na.rm=T)), 
                 size=1, linetype ="dashed", color="gray") 
}
p3

pdf("Output pdfs/jan19 abundance figures/Fig 3_abundance_simple.pdf", width = 6, height = 6)
grid.arrange(p0,p1, p2, p3, nrow=2)
dev.off()

# Percent reduction----
lm4 <- lm(scaled_cost_p_1surv ~ mismatch_p_1surv, data = SWBwide[which(SWBwide$mismatch_p_1surv<0),])
lm5 <- lm(scaled_cost_p_1surv ~ mismatch_p_1surv, data = SWBwide[which(SWBwide$mismatch_p_1surv>0),])
lm6 <- lm(scaled_cost_p_2growth ~ mismatch_p_2growth, data = SWBwide[which(SWBwide$mismatch_p_2growth<0),])
lm7 <- lm(scaled_cost_p_2growth ~ mismatch_p_2growth, data = SWBwide[which(SWBwide$mismatch_p_2growth>0),])
lm8 <- lm(scaled_cost_p_3rec ~ mismatch_p_3rec, data = SWBwide[which(SWBwide$mismatch_p_3rec<0),])
lm9 <- lm(scaled_cost_p_3rec ~ mismatch_p_3rec, data = SWBwide[which(SWBwide$mismatch_p_3rec>0),])
lm10 <- lm(scaled_cost_p_0abun ~ mismatch_p_0abun, data = SWBwide[which(SWBwide$mismatch_p_0abun<0),])
lm11 <- lm(scaled_cost_p_0abun ~ mismatch_p_0abun, data = SWBwide[which(SWBwide$mismatch_p_0abun>0),])

corr4 <- cor.test(SWBwide$mismatch_p_1surv[which(SWBwide$mismatch_p_1surv<0)], SWBwide$scaled_cost_p_1surv[which(SWBwide$mismatch_p_1surv<0)], method="kendall")
corr5 <- cor.test(SWBwide$mismatch_p_1surv[which(SWBwide$mismatch_p_1surv>0)], SWBwide$scaled_cost_p_1surv[which(SWBwide$mismatch_p_1surv>0)], method="kendall")
corr6 <- cor.test(SWBwide$mismatch_p_2growth[which(SWBwide$mismatch_p_2growth<0)], SWBwide$scaled_cost_p_2growth[which(SWBwide$mismatch_p_2growth<0)], method="kendall")
corr7 <- cor.test(SWBwide$mismatch_p_2growth[which(SWBwide$mismatch_p_2growth>0)], SWBwide$scaled_cost_p_2growth[which(SWBwide$mismatch_p_2growth>0)], method="kendall")
corr8 <- cor.test(SWBwide$mismatch_p_3rec[which(SWBwide$mismatch_p_3rec<0)], SWBwide$scaled_cost_p_3rec[which(SWBwide$mismatch_p_3rec<0)], method="kendall")
corr9 <- cor.test(SWBwide$mismatch_p_3rec[which(SWBwide$mismatch_p_3rec>0)], SWBwide$scaled_cost_p_3rec[which(SWBwide$mismatch_p_3rec>0)], method="kendall")
corr10 <- cor.test(SWBwide$mismatch_p_0abun[which(SWBwide$mismatch_p_0abun<0)], SWBwide$scaled_cost_p_0abun[which(SWBwide$mismatch_p_0abun<0)], method="kendall")
corr11 <- cor.test(SWBwide$mismatch_p_0abun[which(SWBwide$mismatch_p_0abun>0)], SWBwide$scaled_cost_p_0abun[which(SWBwide$mismatch_p_0abun>0)], method="kendall")

neg <- broom::glance(corr10)
pos <- broom::glance(corr11)
neglab <- substitute(paste(italic("r"[s]), " = ", est), list(est=signif(neg$estimate, digits=2)))
poslab <- substitute(paste(italic("r"[s]), " = ", est), list(est=signif(pos$estimate, digits=2)))

p0 <- ggplot(SWBwide, aes(mismatch_p_0abun, scaled_cost_p_0abun)) +
  geom_point(size=2.5, alpha=0.5) +
  annotate("text", x=-45, y=15, label = deparse(neglab), parse=TRUE) +
  annotate("text", x = -45, y = 5, label = paste("p = ", signif(neg$p.value, digits=2), sep="")) +
  annotate("text", x=45, y=15, label = deparse(poslab), parse=TRUE) +
  annotate("text", x = 45, y = 5, label = paste("p = ", signif(pos$p.value, digits=2), sep="")) +
  # geom_smooth(data = SWBwide[which(SWBwide$mismatch_p_0abun<0),], method = "lm", se=F, color="black") +
  # geom_smooth(data = SWBwide[which(SWBwide$mismatch_p_0abun>0),], method = "lm", se=F, color="black") +
  # stat_poly_eq(data = SWBwide[which(SWBwide$mismatch_p_0abun<0),],
  #              aes(label = ..rr.label..), label.y = 0.05, label.x = 0.1, size = 3.5,
  #              formula = y~x, parse = TRUE) +
  # stat_fit_glance(data = SWBwide[which(SWBwide$mismatch_p_0abun<0),],
  #                 method = 'spearman', method.args = list(formula = y~x),
  #                 geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")),
  #                 label.y = .13, label.x = 0.1, size = 3.5) +
  # stat_poly_eq(data = SWBwide[which(SWBwide$mismatch_p_0abun>0),],
  #              aes(label = ..rr.label..), label.y.npc = 0.05, label.x.npc = 0.9, size = 3.5, 
  #              formula = y~x, parse = TRUE) +
  # stat_fit_glance(data = SWBwide[which(SWBwide$mismatch_p_0abun>0),],
  #                 method = 'lm', method.args = list(formula = y~x),
  #                 geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
  #                 label.y.npc = 0.13, label.x.npc = 0.9, size = 3.5) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  scale_x_continuous(name ="Mismatch", 
                     breaks = seq(-80,80, by=20), 
                     limits=c(-max(SWBwide$abmis_p_0abun, na.rm=T), max(SWBwide$abmis_p_0abun, na.rm=T))
  ) +
  scale_y_continuous(name ="Density cost, % reduction", 
                     breaks = seq(0, 100, by=10),
                     limits=c(0,100)) +
  ggtitle("a)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p0

neg <- broom::glance(corr4)
pos <- broom::glance(corr5)
neglab <- substitute(paste(italic("r"[s]), " = ", est), list(est=signif(neg$estimate, digits=2)))
poslab <- substitute(paste(italic("r"[s]), " = ", est), list(est=signif(pos$estimate, digits=2)))

p1 <- ggplot(SWBwide, aes(mismatch_p_1surv, scaled_cost_p_1surv)) +
  geom_point(size=2.5, alpha=0.5) +
  annotate("text", x=-50, y=15, label = deparse(neglab), parse=TRUE) +
  annotate("text", x = -50, y = 5, label = paste("p = ", signif(neg$p.value, digits=2), sep="")) +
  annotate("text", x=50, y=15, label = deparse(poslab), parse=TRUE) +
  annotate("text", x = 50, y = 5, label = paste("p = ", signif(pos$p.value, digits=2), sep="")) +
  # geom_smooth(data = SWBwide[which(SWBwide$mismatch_p_1surv<0),], method = "lm", se=F, color="black") +
  # geom_smooth(data = SWBwide[which(SWBwide$mismatch_p_1surv>0),], method = "lm", se=F, color="black") +
  # stat_poly_eq(data = SWBwide[which(SWBwide$mismatch_p_1surv<0),],
  #              aes(label = ..rr.label..), label.y.npc = 0.05, label.x.npc = 0.1, size = 3.5, 
  #              formula = y~x, parse = TRUE) +
  # stat_fit_glance(data = SWBwide[which(SWBwide$mismatch_p_1surv<0),],
  #                 method = 'lm', method.args = list(formula = y~x),
  #                 geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
  #                 label.y.npc = .13, label.x.npc = 0.1, size = 3.5) +
  # stat_poly_eq(data = SWBwide[which(SWBwide$mismatch_p_1surv>0),],
  #              aes(label = ..rr.label..), label.y.npc = 0.05, label.x.npc = 0.9, size = 3.5, 
  #              formula = y~x, parse = TRUE) +
  # stat_fit_glance(data = SWBwide[which(SWBwide$mismatch_p_1surv>0),],
  #                 method = 'lm', method.args = list(formula = y~x),
  #                 geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
  #                 label.y.npc = .13, label.x.npc = 0.9, size = 3.5) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  scale_x_continuous(name ="Mismatch", 
                     breaks = seq(-80,80, by=20), 
                     limits=c(-max(SWBwide$abmis_p_1surv, na.rm=T), max(SWBwide$abmis_p_1surv, na.rm=T))
  ) +
  scale_y_continuous(name ="Survival cost, % reduction", 
                     breaks = seq(0, 100, by=10),
                     limits=c(0,100)) +
  ggtitle("b)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p1

neg <- broom::glance(corr6)
pos <- broom::glance(corr7)
neglab <- substitute(paste(italic("r"[s]), " = ", est), list(est=signif(neg$estimate, digits=2)))
poslab <- substitute(paste(italic("r"[s]), " = ", est), list(est=signif(pos$estimate, digits=2)))

p2 <- ggplot(SWBwide, aes(mismatch_p_2growth, scaled_cost_p_2growth)) +
  geom_point(size=2.5, alpha=0.5) +
  annotate("text", x=-40, y=15, label = deparse(neglab), parse=TRUE) +
  annotate("text", x = -40, y = 5, label = paste("p = ", signif(neg$p.value, digits=2), sep="")) +
  annotate("text", x=40, y=15, label = deparse(poslab), parse=TRUE) +
  annotate("text", x = 40, y = 5, label = paste("p = ", signif(pos$p.value, digits=2), sep="")) +
  # geom_smooth(data = SWBwide[which(SWBwide$mismatch_p_2growth<0),], method = "lm", se=F, color="black") +
  # geom_smooth(data = SWBwide[which(SWBwide$mismatch_p_2growth>0),], method = "lm", se=F, color="black") +
  # stat_poly_eq(data = SWBwide[which(SWBwide$mismatch_p_2growth<0),],
  #              aes(label = ..rr.label..), label.y.npc = 0.05, label.x.npc = 0.1, size = 3.5, 
  #              formula = y~x, parse = TRUE) +
  # stat_fit_glance(data = SWBwide[which(SWBwide$mismatch_p_2growth<0),],
  #                 method = 'lm', method.args = list(formula = y~x),
  #                 geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
  #                 label.y.npc = .13, label.x.npc = 0.1, size = 3.5) +
  # stat_poly_eq(data = SWBwide[which(SWBwide$mismatch_p_2growth>0),],
  #              aes(label = ..rr.label..), label.y.npc = 0.05, label.x.npc = 0.9, size = 3.5, 
  #              formula = y~x, parse = TRUE) +
  # stat_fit_glance(data = SWBwide[which(SWBwide$mismatch_p_2growth>0),],
  #                 method = 'lm', method.args = list(formula = y~x),
  #                 geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
  #                 label.y.npc = .13, label.x.npc = 0.9, size = 3.5) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  scale_x_continuous(name ="Mismatch", 
                     breaks = seq(-40,40, by=20), 
                     limits=c(-max(SWBwide$abmis_p_2growth, na.rm=T), max(SWBwide$abmis_p_2growth, na.rm=T))) +
  scale_y_continuous(name = "Growth cost, % reduction", 
                     breaks = seq(0, 100, by=10),
                     limits=c(0,100)) +
  ggtitle("c)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p2

neg <- broom::glance(corr8)
pos <- broom::glance(corr9)
neglab <- substitute(paste(italic("r"[s]), " = ", est), list(est=signif(neg$estimate, digits=2)))
poslab <- substitute(paste(italic("r"[s]), " = ", est), list(est=signif(pos$estimate, digits=2)))

p3 <- ggplot(SWBwide, aes(mismatch_p_3rec, scaled_cost_p_3rec)) +
  geom_point(size=2.5, alpha=0.5) +
  annotate("text", x=-40, y=15, label = deparse(neglab), parse=TRUE) +
  annotate("text", x = -40, y = 5, label = paste("p = ", signif(neg$p.value, digits=2), sep="")) +
  annotate("text", x=40, y=15, label = deparse(poslab), parse=TRUE) +
  annotate("text", x = 40, y = 5, label = paste("p = ", signif(pos$p.value, digits=2), sep="")) +
  # geom_smooth(data = SWBwide[which(SWBwide$mismatch_p_3rec<0),], method = "lm", se=F, color="black") +
  # geom_smooth(data = SWBwide[which(SWBwide$mismatch_p_3rec>0),], method = "lm", se=F, color="black") +
  # stat_poly_eq(data = SWBwide[which(SWBwide$mismatch_p_3rec<0),],
  #              aes(label = ..rr.label..), label.y.npc = 0.05, label.x.npc = 0.1, size = 3.5, 
  #              formula = y~x, parse = TRUE) +
  # stat_fit_glance(data = SWBwide[which(SWBwide$mismatch_p_3rec<0),],
  #                 method = 'lm', method.args = list(formula = y~x),
  #                 geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
  #                 label.y.npc = .13, label.x.npc = 0.1, size = 3.5) +
  # stat_poly_eq(data = SWBwide[which(SWBwide$mismatch_p_3rec>0),],
  #              aes(label = ..rr.label..), label.y.npc = 0.05, label.x.npc = 0.9, size = 3.5, 
  #              formula = y~x, parse = TRUE) +
  # stat_fit_glance(data = SWBwide[which(SWBwide$mismatch_p_3rec>0),],
  #                 method = 'lm', method.args = list(formula = y~x),
  #                 geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
  #                 label.y.npc = .13, label.x.npc = 0.9, size = 3.5) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  scale_x_continuous(name ="Mismatch", 
                     breaks = seq(-40,40, by=20), 
                     limits=c(-max(SWBwide$abmis_p_3rec, na.rm=T), max(SWBwide$abmis_p_3rec, na.rm=T))) +
  scale_y_continuous(name ="Recruitment cost, % reduction", 
                     breaks = seq(0, 100, by=10),
                     limits=c(0, 100)) +
  ggtitle("d)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p3  



pdf("Output pdfs/jan19 abundance figures/Fig 3_abundance_percent.pdf", width = 9, height = 7.5)
grid.arrange(p0,p1, p2, p3, nrow=2)
dev.off()

## SWB relationship----
p0 <- ggplot(filter(SWBs, vr=='0abun'), aes(med_plot_swb, mismatch_p)) +
  geom_point(size=2.5, alpha=0.5) +
  geom_smooth(method = "lm", se=F, color="black") +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.9, label.x.npc = 0.9, size = 3.5, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = .98, label.x.npc = 0.9, size = 3.5) +
  theme_classic() +
  scale_x_continuous(name ="Median SWB of sampled plots", 
                     breaks = seq(40, 100, by=20)) +
  scale_y_continuous(name ="Density mismatch", 
                     breaks = seq(-60,60, by=20)) +
  geom_hline(yintercept=0, linetype='dashed') +
  ggtitle("a)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p0

p1 <- ggplot(filter(SWBs, vr=='1surv'), aes(med_plot_swb, mismatch_p)) +
  geom_point(size=2.5, alpha=0.5) +
  geom_smooth(method = "lm", se=F, color="black") +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.9, label.x.npc = 0.9, size = 3.5, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = .98, label.x.npc = 0.9, size = 3.5) +
  theme_classic() +
  scale_x_continuous(name ="Median SWB of sampled plots", 
                     breaks = seq(40, 100, by=20)) +
  scale_y_continuous(name ="Survival mismatch", 
                     breaks = seq(-60,60, by=20)) +
  geom_hline(yintercept=0, linetype='dashed') +
  ggtitle("b)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p1

p2 <- ggplot(filter(SWBs, vr=='2growth'), aes(med_plot_swb, mismatch_p)) +
  geom_point(size=2.5, alpha=0.5) +
  geom_smooth(method = "lm", se=F, color="black") +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.9, label.x.npc = 0.9, size = 3.5, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = .98, label.x.npc = 0.9, size = 3.5) +
  theme_classic() +
  scale_x_continuous(name ="Median SWB of sampled plots", 
                     breaks = seq(40, 100, by=20)) +
  scale_y_continuous(name ="Growth mismatch", 
                     breaks = seq(-60,60, by=20)) +
  geom_hline(yintercept=0, linetype='dashed') +
  ggtitle("c)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p2  

p3 <- ggplot(filter(SWBs, vr=='3rec'), aes(med_plot_swb, mismatch_p)) +
  geom_point(size=2.5, alpha=0.5) +
  geom_smooth(method = "lm", se=F, color="black") +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.9, label.x.npc = 0.9, size = 3.5, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = .98, label.x.npc = 0.9, size = 3.5) +
  theme_classic() +
  scale_x_continuous(name ="Median SWB of sampled plots", 
                     breaks = seq(40, 100, by=20)) +
  scale_y_continuous(name ="Recruitment mismatch", 
                     breaks = seq(-60,60, by=20)) +
  geom_hline(yintercept=0, linetype='dashed') +
  ggtitle("d)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p3  

pdf("Output pdfs/jan19 abundance figures/SWB Mismatch.pdf", width = 8, height = 8)
grid.arrange(p0,p1, p2, p3, nrow=2)
dev.off()

## SWB relationship fullrange----

fullrange <- subset(SWBs, overlap_max>90)
p0 <- ggplot(filter(fullrange, vr=='0abun'), aes(med_plot_swb, mismatch_p)) +
  geom_point(size=2.5, alpha=0.5) +
  geom_smooth(method = "lm", se=F, color="gray", linetype='dashed') +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.9, label.x.npc = 0.9, size = 3.5, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = .98, label.x.npc = 0.9, size = 3.5) +
  theme_classic() +
  scale_x_continuous(name ="Median SWB of sampled plots", 
                     breaks = seq(40, 100, by=20)) +
  scale_y_continuous(name ="Density mismatch", 
                     breaks = seq(-60,60, by=20)) +
  geom_hline(yintercept=0, linetype='dashed') +
  ggtitle("a)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p0

p1 <- ggplot(filter(fullrange, vr=='1surv'), aes(med_plot_swb, mismatch_p)) +
  geom_point(size=2.5, alpha=0.5) +
  geom_smooth(method = "lm", se=F, color="gray", linetype='dashed') +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.9, label.x.npc = 0.9, size = 3.5, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = .98, label.x.npc = 0.9, size = 3.5) +
  theme_classic() +
  scale_x_continuous(name ="Median SWB of sampled plots", 
                     breaks = seq(40, 100, by=20)) +
  scale_y_continuous(name ="Survival mismatch", 
                     breaks = seq(-60,60, by=20)) +
  geom_hline(yintercept=0, linetype='dashed') +
  ggtitle("b)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p1

p2 <- ggplot(filter(fullrange, vr=='2growth'), aes(med_plot_swb, mismatch_p)) +
  geom_point(size=2.5, alpha=0.5) +
  geom_smooth(method = "lm", se=F, color="black") +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.9, label.x.npc = 0.9, size = 3.5, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = .98, label.x.npc = 0.9, size = 3.5) +
  theme_classic() +
  scale_x_continuous(name ="Median SWB of sampled plots", 
                     breaks = seq(40, 100, by=20)) +
  scale_y_continuous(name ="Growth mismatch", 
                     breaks = seq(-60,60, by=20)) +
  geom_hline(yintercept=0, linetype='dashed') +
  ggtitle("c)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p2  

p3 <- ggplot(filter(fullrange, vr=='3rec'), aes(med_plot_swb, mismatch_p)) +
  geom_point(size=2.5, alpha=0.5) +
  geom_smooth(method = "lm", se=F, color="black") +
  stat_poly_eq(aes(label = ..rr.label..), label.y.npc = 0.9, label.x.npc = 0.9, size = 3.5, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = .98, label.x.npc = 0.9, size = 3.5) +
  theme_classic() +
  scale_x_continuous(name ="Median SWB of sampled plots", 
                     breaks = seq(40, 100, by=20)) +
  scale_y_continuous(name ="Recruitment mismatch", 
                     breaks = seq(-60,60, by=20)) +
  geom_hline(yintercept=0, linetype='dashed') +
  ggtitle("d)") +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=13), plot.title = element_text(size=16, face="bold"))
p3  

pdf("Output pdfs/jan19 abundance figures/SWB Mismatch fullrange.pdf", width = 8, height = 8)
grid.arrange(p0,p1, p2, p3, nrow=2)
dev.off()



