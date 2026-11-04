## 5_Mismatch figures
## Figures of distribution-performance mismatch and demographic cost.
## Author: Teresa Bohner
## Project: where western trees grow vs. where they thrive (see README)

library(bbmle)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(ggalt)
library(tidyverse)
library(data.table)
library(doBy)
library(lme4)

source("R scripts/5_Mismatch trait data compilation.R")


## wide format----
SWBwide <- data.table::dcast(data.table(SWBs), 
                             gensp 
                             + hydrophyte
                             + totscore
                             + wood_density
                             + value_RootDepth
                             + sla_value
                             + DispMode
                             + logSeedMass
                             + value_Lifespan
                             + lat_extent.y~vr, 
                             value.var=c("mismatch_p",
                                         "abmis_p",
                                         "cost_p",
                                         "scaled_cost_p","sens"))
## ggplot labelling miscellainy----
vr_names <- list(
  "0abun"="Abundance",
  '1surv'="Survival",
  '2growth'="Growth",
  '3rec'="Recruitment"
)
vr_labeller <- function(variable,value){
  return(vr_names[value])
  
}

SWBs$star <- ""
SWBs$star[SWBs$novitals == 1] <- "A"
SWBs$star[SWBs$novitals == 2] <- "B"
SWBs$star[SWBs$novitals == 3] <- "C"
SWBs$star[SWBs$novitals == 4] <- ""


SWBs$star2 <- as.factor(ifelse(is.na(SWBs$abmis_p),0,1))
SWBs$cost2 <- ifelse(is.na(SWBs$cost_p), 0, SWBs$cost_p)

surv <- subset(SWBs, vr=="1surv")
grow <- subset(SWBs, vr=="2growth")
rec <- subset(SWBs, vr=="3rec")

ggplot(SWBs, aes(mismatch_p)) +
  geom_histogram(aes(fill=vr))+
  facet_wrap(~vr, ncol=1, scales="free_x", labeller=vr_labeller) +  
  geom_vline(aes(xintercept=0), linetype="dotted", size=1) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(strip.background = element_rect(fill="light gray", colour="black", size=1)) +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  labs(x="Mismatch", y="Number of Species")+
  theme(text = element_text(size=18), strip.text.x = element_text(size = 18, face = "bold"))


## Sensitivity plots----
pdf("Output pdfs/species sensitivity.pdf")
ggplot(SWBs, aes(x=reorder(gensp, -spid), y=scaled_sens)) +  coord_flip() +
  geom_bar(stat="identity")+ #, position=position_dodge()
  facet_grid(. ~vr, scales="free", labeller=vr_labeller) +
  xlab("Species\n") + ylab("Sensitivity (unity normalized)") +
  theme_classic() +
  theme(strip.background = element_rect(fill="light gray", colour="black", size=1)) +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(panel.grid.major.y = element_line(colour = "dark gray", size = 0.3))
ggplot(data=SWBs, aes(y=scaled_sens, x=reorder(gensp, -spid))) + 
  geom_bar(aes(fill=vr), stat="identity", color="white", size=0.1) +
  geom_hline(aes(yintercept=0), linetype="dotted") + 
  xlab("Species\n") + ylab("\n sum of sensitivity (unity normalized)") + coord_flip() + 
  labs(color="Demographic Rate") +
  scale_color_discrete(labels=c("Survival", "Growth", "Recruitment"), l=50) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(panel.grid.major.y = element_line(colour = "dark gray", size = 0.3)) +
  geom_text(aes(y=-0.1, label=star), size=2)
dev.off()

## Bring in abundance data to weight plots
abundance <- read.csv("Ch 1/data for stan/abundance.csv")
avgsens <- merge(abundance, SWBwide, by="gensp")
avgsens$sweighted <- avgsens$sens_1surv*avgsens$totba
avgsens$gweighted <- avgsens$sens_2growth*avgsens$totba
avgsens$rweighted <- avgsens$sens_3rec*avgsens$totba

# Aggregate species level weighted sensitivity to plot level
plotsens <- group_by(avgsens, plot.uniqueID) %>%
  summarize(lat=first(LAT),
            lon=first(LON),
            swb=first(ptcoeff),
            s_wavg = sum(sweighted, na.rm=T)/sum(totba, na.rm=T),
            g_wavg = sum(gweighted, na.rm=T)/sum(totba, na.rm=T),
            r_wavg = sum(rweighted, na.rm=T)/sum(totba, na.rm=T),
            dom = first(gensp[which(totba==max(totba))]))

# Bin data for plotting
plotsens <- transform(plotsens, bins=cut(swb, 95))
avg.plot <- summaryBy(swb + s_wavg + g_wavg + r_wavg ~ bins, data=plotsens)


pdf("Output pdfs/weighted sensitivity.pdf")
plot(plotsens$s_wavg~plotsens$swb, col=plotsens$dom
     , ylab="BA weighted plot sensitivity", xlab="SWB", main="Survival")
abline(lm(plotsens$s_wavg~plotsens$swb))
plot(s_wavg.mean ~ swb.mean, data=avg.plot, pch=19
     , ylab="BA weighted plot sensitivity", xlab="SWB", main="Survival")
abline(lm(plotsens$s_wavg~plotsens$swb), col="red")

plot(plotsens$g_wavg~plotsens$swb, col=plotsens$dom
     , ylab="BA weighted plot sensitivity", xlab="SWB", main="Growth")
abline(lm(plotsens$g_wavg~plotsens$swb))
plot(g_wavg.mean ~ swb.mean, data=avg.plot, pch=19
     , ylab="BA weighted plot sensitivity", xlab="SWB", main="Growth")
abline(lm(plotsens$g_wavg~plotsens$swb), col="red")

plot(plotsens$r_wavg~plotsens$swb, col=plotsens$dom
     , ylab="BA weighted plot sensitivity", xlab="SWB", main="Recruitment")
abline(lm(plotsens$r_wavg~plotsens$swb))
plot(r_wavg.mean ~ swb.mean, data=avg.plot, pch=19
     , ylab="BA weighted plot sensitivity", xlab="SWB", main="Recruitment")
abline(lm(plotsens$r_wavg~plotsens$swb), col="red")
dev.off()

SWBsummary <- group_by(SWBs, gensp) %>%
  mutate(sign = ifelse(mismatch_p>0, "pos",
                       ifelse(mismatch_p==0, "zer", "neg"))) %>%
  summarize(ntotal = length(which(!is.na(sign))),
            npos = length(which(sign=="pos")),
            nneg = length(which(sign=="neg")),
            nzer = length(which(sign=="zer")),
            p_pos = npos/ntotal,
            p_neg = nneg/ntotal)


## Occurrence figures ----
pdf("Output pdfs/all VR summary panels occurrence.pdf")

ggplot(data=SWBs, aes(x=cost_p, y=reorder(gensp, -spid))) +
  geom_segment(aes(yend=reorder(gensp, -spid)), xend=0, colour="black") +
  geom_point(aes(color=family)) +
  facet_grid(. ~vr, scales="free", labeller=vr_labeller) +
  # scale_color_gradient(low = "blue", high="red") +
  geom_vline(aes(xintercept=0), linetype="dotted") +
  ylab("Species\n") + xlab("\n demographic cost") + 
  labs(color="SWB Mismatch") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(strip.background = element_rect(fill="light gray", colour="black", size=1)) +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(panel.grid.major.y = element_line(colour = "dark gray", size = 0.3)) +
  labs(caption="(based on occurrence model)")

ggplot(data=SWBs, aes(y=mismatch_p, x=reorder(gensp, -spid), color=mismatch_p)) +
  facet_grid(. ~vr, scales="free", labeller=vr_labeller) +
  geom_point() +
  scale_colour_gradient(low = "blue", high="red") +
  geom_hline(aes(yintercept=0), linetype="dotted") +
  xlab("Species\n") + ylab("\n SWB mismatch") + coord_flip() +
  labs(color="SWB Mismatch") +
  theme_classic() +
  theme(strip.background = element_rect(fill="light gray", colour="black", size=1)) +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(panel.grid.major.y = element_line(colour = "dark gray", size = 0.3))

ggplot(data=SWBs, aes(y=mismatch_p, x=reorder(gensp, -spid), fill=mismatch_p)) +
  facet_grid(. ~vr, scales="free", labeller=vr_labeller) +
  geom_bar(stat="identity") +
  scale_fill_gradient(low = "blue", high="red") +
  geom_hline(aes(yintercept=0), linetype="dotted") +
  xlab("Species\n") + ylab("\n SWB mismatch") + coord_flip() +
  labs(color="SWB Mismatch") +
  theme_classic() +
  theme(strip.background = element_rect(fill="light gray", colour="black", size=1)) +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(panel.grid.major.y = element_line(colour = "dark gray", size = 0.3)) 

ggplot(data=SWBs, aes(y=cost_p, x=reorder(gensp, -spid), fill=mismatch_p)) +
  facet_grid(. ~vr, scales="free", labeller=vr_labeller) +
  geom_bar(stat="identity") +
  scale_fill_gradient(low = "blue", high="red") +
  geom_hline(aes(yintercept=0), linetype="dotted") +
  xlab("Species\n") + ylab("\n demographic cost") + coord_flip() +
  labs(color="SWB Mismatch") +
  theme_classic() +
  theme(strip.background = element_rect(fill="light gray", colour="black", size=1)) +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(panel.grid.major.y = element_line(colour = "dark gray", size = 0.3)) 

ggplot(data=SWBs, aes(y=cost_p, x=reorder(gensp, -spid), fill=scaled_cost_p)) +
  facet_grid(. ~vr, scales="free", labeller=vr_labeller) +
  geom_bar(stat="identity") +
  scale_fill_gradient(low = "blue", high="red") +
  geom_hline(aes(yintercept=0), linetype="dotted") +
  xlab("Species\n") + ylab("\n demographic cost (unity normalized") + coord_flip() +
  labs(color="SWB Mismatch") +
  theme_classic() +
  theme(strip.background = element_rect(fill="light gray", colour="black", size=1)) +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(panel.grid.major.y = element_line(colour = "dark gray", size = 0.3)) 

ggplot(data=SWBs, aes(y=mismatch_p, x=reorder(gensp, -spid), fill=family)) +
  facet_grid(. ~vr, scales="free", labeller=vr_labeller) +
  geom_bar(stat="identity") +
  # scale_fill_gradient(low = "blue", high="red") +
  geom_hline(aes(yintercept=0), linetype="dotted") +
  xlab("Species\n") + ylab("\n SWB mismatch") + coord_flip() +
  labs(color="SWB Mismatch") +
  theme_classic() +
  theme(strip.background = element_rect(fill="light gray", colour="black", size=1)) +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(panel.grid.major.y = element_line(colour = "dark gray", size = 0.3)) 

ggplot(data=SWBs, aes(y=mismatch_p, x=family, fill=family))+
  geom_boxplot() +
  facet_grid(. ~vr, scales="free", labeller=vr_labeller) +
  theme_classic() +
  theme(legend.position = "false") +
  geom_hline(aes(yintercept=0), linetype="dotted") +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))

ggplot(data=SWBs, aes(y=mismatch_p, x=family, fill=vr))+
  geom_boxplot(position="dodge") +
  theme_classic() +
  geom_hline(aes(yintercept=0), linetype="dotted") +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))

ggplot(data=SWBs, aes(y=mismatch_p, x=vr, fill=family))+
  geom_boxplot() +
  facet_wrap(~family)+
  theme_classic() +
  geom_hline(aes(yintercept=0), linetype="dotted") +
  theme(legend.position = "false")

ggplot(data=SWBs, aes(y=mismatch_p, x=reorder(gensp, -meanmis_p))) +
  geom_linerange(aes(ymin=minmis_p, ymax=maxmis_p, group=gensp)) +
  geom_point(aes(color=vr)) +
  geom_hline(aes(yintercept=0), linetype="dotted") +
  xlab("Species\n") + ylab("\n SWB mismatch") + coord_flip() +
  labs(color="Performance Rate") +
  scale_color_discrete(labels=c("Density","Survival", "Growth", "Recruitment"), l=50) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(panel.grid.major.y = element_line(colour = "dark gray", size = 0.3)) +
  geom_text(aes(y=-80, label=star), size=2)

ggplot(data=SWBs, aes(y=mismatch_p, x=reorder(gensp, -meanmis_p))) +
  geom_linerange(aes(ymin=minmis_p, ymax=maxmis_p, group=gensp)) +
  geom_point(aes(color=med_plot_swb)) +
  geom_hline(aes(yintercept=0), linetype="dotted") +
  xlab("Species\n") + ylab("\n SWB mismatch") + coord_flip() +
  labs(color="Median plot SWB") +
  scale_color_viridis(direction=-1) +
  # scale_color_discrete(labels=c("Density","Survival", "Growth", "Recruitment"), l=50) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  theme(panel.grid.major.y = element_line(colour = "dark gray", size = 0.3)) +
  geom_text(aes(y=-80, label=star), size=2)


ggplot(data=SWBs, aes(x=abmis_p, y=cost_p)) +
  geom_point()+ #aes(color=cost_p)
  # scale_colour_gradient(low = "blue", high="red") +
  facet_wrap(~vr, scales="free", labeller=vr_labeller)+
  # labs(color="Demographic \n Cost") +
  xlab("Absolute value of SWB mismatch") +
  ylab("Demographic cost")+
  theme_classic() +
  theme(strip.background = element_rect(fill="light gray", colour="black", size=1)) +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  labs(caption="(based on occurrence model)")

#bivariate mismatch_p
ggplot(data=SWBwide) +
  geom_point(aes(x=abmis_p_1surv, y=abmis_p_2growth)) +
  xlab("Absolute value of survival SWB mismatch") +
  ylab("Absolute value of growth SWB mismatch")+
  theme_classic() +
  theme(axis.line= element_line(size=0.5)) +
  labs(caption="(based on occurrence model)")

ggplot(data=SWBwide) +
  geom_point(aes(x=abmis_p_1surv, y=abmis_p_3rec)) +
  xlab("Absolute value of survival SWB mismatch") +
  ylab("Absolute value of recruitment SWB mismatch")+
  theme_classic() +
  theme(axis.line= element_line(size=0.5)) +
  labs(caption="(based on occurrence model)")

ggplot(data=SWBwide) +
  geom_point(aes(x=abmis_p_2growth, y=abmis_p_3rec)) +
  xlab("Absolute value of growth SWB mismatch") +
  ylab("Absolute value of recruitment SWB mismatch")+
  theme_classic() +
  theme(axis.line= element_line(size=0.5)) +
  labs(caption="(based on occurrence model)")

#bivariate cost_p
ggplot(data=SWBwide) +
  geom_point(aes(x=cost_p_1surv, y=cost_p_2growth), color="red") +
  xlab("Survival cost of  SWB mismatch") +
  ylab("Growth cost of SWB mismatch")+
  theme_classic() +
  theme(axis.line= element_line(size=0.5)) +
  labs(caption="(based on occurrence model)")

ggplot(data=SWBwide) +
  geom_point(aes(x=cost_p_1surv, y=cost_p_3rec), color="red") +
  xlab("Survival cost of  SWB mismatch") +
  ylab("Recruitment cost of SWB mismatch")+
  theme_classic() +
  theme(axis.line= element_line(size=0.5)) +
  labs(caption="(based on occurrence model)")

ggplot(data=SWBwide) +
  geom_point(aes(x=cost_p_2growth, y=cost_p_3rec), color="red") +
  xlab("Growth cost of  SWB mismatch") +
  ylab("Recruitment cost of SWB mismatch")+
  theme_classic() +
  theme(axis.line= element_line(size=0.5))+
  labs(caption="(based on occurrence model)")

dev.off()

ggplot(data=SWBs, aes(x=abmis_p, y=cost_p)) +
  geom_point()+ #aes(color=cost_n)
  # geom_smooth(method=lm) +
  # scale_colour_gradient(low = "blue", high="red") +
  facet_wrap(~vr, scales="free", labeller=vr_labeller)+
  # labs(color="Demographic \n Cost") +
  xlab("Absolute value of SWB mismatch") +
  ylab("Demographic cost")+
  theme_classic() +
  theme(strip.background = element_rect(fill="light gray", colour="black", size=1)) +
  theme(axis.line= element_line(size=0.5)) +
  labs(caption="(based on occurrence model)")

## Cost mismatch relationship----
## Spearman rank-correlation
corr_s <- cor.test(SWBwide$abmis_p_1surv, SWBwide$cost_p_1surv, method = 'spearman', exact=F)
corr_g <- cor.test(SWBwide$abmis_p_2growth, SWBwide$cost_p_2growth, method = 'spearman', exact=F)
corr_r <- cor.test(SWBwide$abmis_p_3rec, SWBwide$cost_p_3rec, method = 'spearman', exact=F)


pdf("Output pdfs/mismatch-cost scatter.pdf")

ggplot(data=SWBs, aes(x=abmis_p, y=cost_p, color=hydrophyte)) +
  geom_point()+ #aes(color=cost_n)
  # geom_smooth(method=lm) +
  # scale_colour_gradient(low = "blue", high="red") +
  facet_wrap(~vr, scales="free", labeller=vr_labeller)+
  # labs(color="Demographic \n Cost") +
  xlab("Absolute value of SWB mismatch") +
  ylab("Demographic cost")+
  theme_classic() +
  theme(strip.background = element_rect(fill="light gray", colour="black", size=1)) +
  theme(axis.line= element_line(size=0.5)) +
  labs(caption="(based on occurrence model)")

ggplot(data=SWBs, aes(x=abmis_p, y=scaled_cost_p)) +
  geom_point()+ #aes(color=cost_n)
  # geom_smooth(method=lm) +
  # scale_colour_gradient(low = "blue", high="red") +
  facet_wrap(~vr, scales="free", labeller=vr_labeller)+
  # labs(color="Demographic \n Cost") +
  xlab("Absolute value of SWB mismatch") +
  ylab("Demographic cost")+
  theme_classic() +
  theme(strip.background = element_rect(fill="light gray", colour="black", size=1)) +
  theme(axis.line= element_line(size=0.5)) +
  labs(caption="(based on occurrence model)")

ggplot(data=SWBwide, aes(x=abmis_p_1surv, y=cost_p_1surv, color=hydrophyte)) +
  geom_point(size=3)+ #aes(color=cost_p)
  # scale_colour_gradient(low = "blue", high="red") +
  # labs(color="Demographic \n Cost") +
  # geom_text(aes(y=cost_p+.01, label=gensp), size=4)+
  # geom_smooth(method = lm, se = FALSE, size=1, color="black")+
  geom_text_repel(aes(label=gensp)) +
  xlab("Absolute value of SWB mismatch") +
  ylab("Demographic cost")+
  ggtitle("Survival") +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))+
  theme(strip.background = element_rect(fill="light gray", colour="black", size=1)) +
  theme(text = element_text(size=18), strip.text.x = element_text(size = 18, face = "bold"))

ggplot(data=SWBwide, aes(x=abmis_p_2growth, y=cost_p_2growth)) +
  geom_point(size=3)+ #aes(color=cost_p)
  # scale_colour_gradient(low = "blue", high="red") +
  # labs(color="Demographic \n Cost") +
  # geom_text(aes(y=cost_p+.01, label=gensp), size=4)+
  # geom_smooth(method = lm, se = FALSE, size=1, color="black")+
  xlab("Absolute value of SWB mismatch") +
  ylab("Demographic cost")+
  ggtitle("Growth") +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))+
  theme(strip.background = element_rect(fill="light gray", colour="black", size=1)) +
  theme(text = element_text(size=18), strip.text.x = element_text(size = 18, face = "bold"))

ggplot(data=SWBwide, aes(x=abmis_p_3rec, y=cost_p_3rec)) +
  geom_point(size=3)+ #aes(color=cost_p)
  # scale_colour_gradient(low = "blue", high="red") +
  # labs(color="Demographic \n Cost") +
  # geom_text(aes(y=cost_p+.01, label=gensp), size=4)+
  # geom_smooth(method = lm, se = FALSE, size=1, color="black")+
  xlab("Absolute value of SWB mismatch") +
  ylab("Demographic cost")+
  ggtitle("Recruitment") +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA))+
  theme(strip.background = element_rect(fill="light gray", colour="black", size=1)) +
  theme(text = element_text(size=18), strip.text.x = element_text(size = 18, face = "bold"))
dev.off()

traits <- dplyr::select(SWBwide, c(abmis_p_2growth, cost_p_2growth)) %>% data.frame()
rownames(traits) <- SWBwide$gensp
droplist <- SWBwide$gensp[which(is.na(SWBwide$abmis_p_2growth))] 
tree1 <- drop.tip(tr$phylo, droplist)
traits1 <- traits[which(!rownames(traits)%in%droplist),] 
hist(traits$cost_p_2growth)
hist(log(traits$cost_p_2growth))


## Traits----
## Exploratory figures
fullrange <- subset(SWBs, overlap_max>90)
ggplot(fullrange, aes(x=med_plot_swb, y=mismatch_p)) +
  geom_text_repel(aes(label=gensp, color=family),
                  size=2,
                  segment.size = 0.2,
                  box.padding = 0.2,
                  segment.color = 'grey50') +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3,
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1.02, size = 3) +
  facet_wrap(~vr)

ggplot(filter(SWBs, vr=="0abun"), aes(x=med_plot_swb, y=meanmis_p)) +
  geom_text_repel(aes(label=gensp, color=family),
                  size=2,
                  segment.size = 0.2,
                  box.padding = 0.2,
                  segment.color = 'grey50') +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3,
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1.02, size = 3) 

ggplot(filter(SWBs, vr=="0abun"), aes(x=med_plot_temp, y=meanmis_p)) +
  geom_text_repel(aes(label=gensp, color=gymno2),
                  size=2,
                  segment.size = 0.2,
                  box.padding = 0.2,
                  segment.color = 'grey50') +
  geom_point(aes(color=gymno2)) +
  geom_smooth(method=lm, aes(color=gymno2)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3,
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1.02, size = 3) 

ggplot(filter(SWBs, vr=="0abun"), aes(x=med_plot_swb, y=med_plot_temp, color=meanmis_p)) +
  geom_text_repel(aes(label=gensp),
                  size=2,
                  segment.size = 0.2,
                  box.padding = 0.2,
                  segment.color = 'grey50') +
  geom_point() +
  scale_colour_gradient2() +
  # geom_smooth(method=lm, aes(color=gymno2)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted')

ggplot(filter(fullrange, !is.na(mismatch_p)), aes(x=med_plot_swb, y=med_plot_temp, color=mismatch_p)) +
  geom_text_repel(aes(label=gensp),
                  size=2,
                  segment.size = 0.2,
                  box.padding = 0.2,
                  segment.color = 'grey50') +
  geom_point() +
  scale_colour_gradient2(low = "dark blue", mid = "white",
                         high = "dark red") +
  facet_wrap(~vr) +
  # geom_smooth(method=lm, aes(color=gymno2)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted')

ggplot(filter(SWBs, vr=="0abun"), aes(x=Seed.spread.rate, y=meanmis_p)) +
  geom_boxplot()
# theme_classic() +
# theme(panel.border = element_rect(size=1, fill = NA)) +
geom_hline(aes(yintercept=0), linetype='dotted')


ggplot(fullrange, aes(x=med_plot_temp, y=mismatch_p)) +
  geom_text_repel(aes(label=gensp, color=family),
                  size=2,
                  segment.size = 0.2,
                  box.padding = 0.2,
                  segment.color = 'grey50') +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1.02, size = 3)

ggplot(fullrange, aes(mismatch_p)) +
  geom_histogram(aes(y = ..count..)) +
  geom_density(aes(y = ..density..*(60*6))) +
  facet_wrap(~vr) +
  geom_vline(xintercept=0, linetype="dashed")


pdf("Output pdfs/mismatch predictor scatterplots.pdf")
ggplot(SWBs, aes(x=family, y=mismatch_p, fill=vr)) +
  geom_boxplot() +
  theme_classic()+
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))

ggplot(SWBs, aes(x=med_plot_swb, y=mismatch_p)) +
  # geom_text_repel(aes(label=gensp, color=family),
  # size=2,
  # segment.size = 0.2,
  # box.padding = 0.2,
  # segment.color = 'grey50') +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3,
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1.05, size = 3) +
  facet_wrap(~vr)


ggplot(SWBs, aes(x=med_plot_temp, y=mismatch_p)) +
  # geom_text_repel(aes(label=gensp, color=family),
  #                 size=2,
  #                 segment.size = 0.2,
  #                 box.padding = 0.2,
  #                 segment.color = 'grey50') +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE) +
  stat_fit_glance(method = 'lm', method.args = list(formula = y~x),
                  geom = 'text', aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")), 
                  label.y.npc = 1.05, size = 3)

ggplot(filter(SWBs, !is.na(mismatch_p)), aes(x=med_plot_swb, y=med_plot_temp, color=mismatch_p)) +
  geom_text_repel(aes(label=gensp),
                  size=2,
                  segment.size = 0.2,
                  box.padding = 0.2,
                  segment.color = 'grey50') +
  geom_point() +
  scale_colour_gradient2(low = "dark blue", mid = "white",
                         high = "dark red") +
  facet_wrap(~vr) +
  # geom_smooth(method=lm, aes(color=gymno2)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted')

ggplot(filter(SWBs, !is.na(scaled_cost_p)), aes(x=med_plot_swb, y=mismatch_p, color=scaled_cost)) +
  # geom_text_repel(aes(label=gensp),
  #                 size=2,
  #                 segment.size = 0.2,
  #                 box.padding = 0.2,
  #                 segment.color = 'grey50') +
  geom_point() +
  scale_colour_gradient() +
  facet_wrap(~vr) +
  # geom_smooth(method=lm, aes(color=gymno2)) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted')

ggplot(SWBs, aes(x=med_plot_elev, y=mismatch_p)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE)

ggplot(SWBs, aes(x=med_range_lat, y=mismatch_p)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE)

ggplot(SWBs, aes(x=wood_density, y=mismatch_p)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE)

ggplot(SWBs, aes(x=Plant.tolerance.to.shade, y=mismatch_p)) +
  geom_boxplot(aes(fill=Plant.tolerance.to.shade)) +
  # geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') 

ggplot(SWBs, aes(x=Plant.tolerance.to.drought, y=mismatch_p)) +
  geom_boxplot(aes(fill=Plant.tolerance.to.drought)) +
  # geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted')

ggplot(SWBs, aes(x=Anaerobic.Tolerance, y=mismatch_p)) +
  geom_boxplot(aes(fill=Anaerobic.Tolerance)) +
  # geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') 

ggplot(SWBs, aes(x=sla_value, y=mismatch_p)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE)
ggplot(SWBs, aes(x=bien_value_leaf_nitrogen_content_per_leaf_dry_mass, y=mismatch_p)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE)
dev.off()

pdf("Output pdfs/mismatch magnitude predictor scatterplots.pdf")
ggplot(SWBs, aes(x=family, y=abmis_p, fill=vr)) +
  geom_boxplot() +
  theme_classic()+
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))

ggplot(SWBs, aes(x=med_plot_swb, y=abmis_p)) +
  geom_text_repel(aes(label=gensp, color=family),
                  size=2,
                  segment.size = 0.2,
                  box.padding = 0.2,
                  segment.color = 'grey50') +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE) +
  facet_wrap(~vr)


ggplot(SWBs, aes(x=med_plot_temp, y=abmis_p)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE)

ggplot(SWBs, aes(x=med_plot_elev, y=abmis_p)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE)

ggplot(SWBs, aes(x=med_range_lat, y=abmis_p)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE)

ggplot(SWBs, aes(x=wood_density, y=abmis_p)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE)

ggplot(SWBs, aes(x=Plant.tolerance.to.shade, y=abmis_p)) +
  geom_boxplot(aes(fill=Plant.tolerance.to.shade)) +
  # geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted')

ggplot(SWBs, aes(x=Plant.tolerance.to.drought, y=abmis_p)) +
  geom_boxplot(aes(fill=Plant.tolerance.to.drought)) +
  # geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted')

ggplot(SWBs, aes(x=Anaerobic.Tolerance, y=abmis_p)) +
  geom_boxplot(aes(fill=Anaerobic.Tolerance)) +
  # geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted')

ggplot(SWBs, aes(x=sla_value, y=abmis_p)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE)
ggplot(SWBs, aes(x=bien_value_leaf_nitrogen_content_per_leaf_dry_mass, y=abmis_p)) +
  geom_point(aes(color=family)) +
  geom_smooth(method=lm) +
  facet_wrap(~vr) +
  theme_classic() +
  theme(panel.border = element_rect(size=1, fill = NA)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  # theme(legend.position = "none") +
  stat_poly_eq(aes(label = ..rr.label..), size = 3, 
               formula = y~x, parse = TRUE)
dev.off()


