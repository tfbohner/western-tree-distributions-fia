## FIA cost and mismatch summaries 

library(tidyverse)
library(data.table)
library(reshape2)

data.path <- "Processed Data/Analysis Ready Data/"
table.path <- "Output tables/"

## Mismatch and cost data----
gSWBs <- read.csv(paste0(data.path, 'growswbs_temp.csv'))
sSWBs <- read.csv(paste0(data.path, 'survswbs_temp.csv'))
rSWBs <- read.csv(paste0(data.path, 'recswbs_temp.csv'))

nSWBs <- read.csv(paste0(data.path, 'abunswbs_temp.csv'))[, c(1:22, 31)]

## Remove species with unreliable model predictions here not later
# reclist <- c("Fraxinus latifolia","Larix lyallii","Pinus attenuata")
# rSWBs <- subset(rSWBs, !(gensp%in%reclist))


## Get columns ready for long-format join
sSWBs <- rename(sSWBs, "vr_swb"="opt_swb", "vr_hat_surv"="s_hat_opt", "vr_hat_pres"="s_hat_pres", 
                "vr_hat_wet"="s_hat_wet", "vr_hat_dry"="s_hat_dry", "vr_hat_opt"="s_hat_opt")
gSWBs <- rename(gSWBs, "vr_swb"="opt_swb", "vr_hat_grow"="g_hat_opt", "vr_hat_pres"="g_hat_pres",
                "vr_hat_wet"="g_hat_wet", "vr_hat_dry"="g_hat_dry", "vr_hat_opt"="g_hat_opt")
rSWBs <- rename(rSWBs, "vr_swb"="opt_swb", "vr_hat_rec"="r_hat_opt", "vr_hat_pres"="r_hat_pres", 
                "vr_hat_wet"="r_hat_wet", "vr_hat_dry"="r_hat_dry", "vr_hat_opt"="r_hat_opt")

nSWBs <- rename(nSWBs, "vr_swb"="opt_swb", "vr_hat_abun"="n_hat_opt", "vr_hat_pres"="n_hat_pres",
                "vr_hat_wet"="n_hat_wet", "vr_hat_dry"="n_hat_dry", "vr_hat_opt"="n_hat_opt")

## Adjust debated species names
# sSWBs$gensp <- sub("Lithocarpus","Notholithocarpus",sSWBs$gensp)
# sSWBs$gensp <- sub("Chamaecyparis nootkatensis","Callitropsis nootkatensis",sSWBs$gensp)
# sSWBs$gensp <- sub("Abies shastensis","Abies magnifica var. shastensis",sSWBs$gensp)
# gSWBs$gensp <- sub("Lithocarpus","Notholithocarpus",gSWBs$gensp)
# gSWBs$gensp <- sub("Chamaecyparis nootkatensis","Callitropsis nootkatensis",gSWBs$gensp)
# gSWBs$gensp <- sub("Abies shastensis","Abies magnifica var. shastensis",gSWBs$gensp)
# rSWBs$gensp <- sub("Lithocarpus","Notholithocarpus",rSWBs$gensp)
# rSWBs$gensp <- sub("Chamaecyparis nootkatensis","Callitropsis nootkatensis",rSWBs$gensp)
# rSWBs$gensp <- sub("Abies shastensis","Abies magnifica var. shastensis",rSWBs$gensp)
## Get individual vr datasets to have same spid as longest (surv) for ordering in plots
gSWBs <- dplyr::select(gSWBs,-spid) %>%
  merge(dplyr::select(nSWBs, spid, gensp), by="gensp", all=T)
rSWBs<- dplyr::select(rSWBs,-spid) %>%
  merge(dplyr::select(nSWBs, spid, gensp), by="gensp", all=T)
sSWBs<- dplyr::select(sSWBs,-spid) %>%
  merge(dplyr::select(nSWBs, spid, gensp), by="gensp", all=T)

## Labels...number makes order in ggplot easier
gSWBs$vr <- "2growth"
sSWBs$vr <- "1surv"
rSWBs$vr <- "3rec"
nSWBs$vr <- "0abun"

## Rank order for increasing mismatch (for plotting later)
spplist_p <- sSWBs[,c("gensp","mismatch_p")]
spplist_p <- spplist_p[order(spplist_p$mismatch_p),]
spplist_p$rankp <- seq(1:nrow(spplist_p))


## Bind all vital rates together
SWBs <- bind_rows(sSWBs, gSWBs, rSWBs, nSWBs)
SWBs$vr <- as.factor(SWBs$vr)

## Bring in mismatch rank
SWBs <- merge(SWBs, spplist_p[c(1,3)], by="gensp")

## Calculate absolute value of mismatch
SWBs$abmis_p <- abs(SWBs$mismatch_p)

SWBs <- group_by(SWBs, vr)%>%
  mutate(scaled_cost_p = ifelse(cost_p>range, 100, cost_p*100/range),
         scaled_sens = (sens - min(sens, na.rm=T))/(max(sens, na.rm=T)-min(sens, na.rm=T)))

## Calculate minimum, maximum, and range of mismatch and cost across vr
SWBs <- group_by(SWBs, gensp) %>%
  mutate(minmis_p = min(mismatch_p, na.rm=T),
         maxmis_p = max(mismatch_p, na.rm=T),
         meanmis_p = mean(mismatch_p, na.rm=T),
         diffmis_p = maxmis_p-minmis_p,
         mincost_p = min(scaled_cost_p, na.rm=T),
         maxcost_p = max(scaled_cost_p, na.rm=T),
         diffcost_p = maxcost_p-mincost_p,
         novitals = length(which(!is.na(cost_p))))

## get trait values----
species_traits <- read.csv(paste0(data.path, "all_traits.csv"))
SWBs <- merge(SWBs, species_traits, by="gensp", all.x=T)

## range overlap *** not sure where the code to make this csv is
ranges <- read.csv(paste0(table.path, "range comparison.csv"))
ranges <- group_by(ranges, FIA.sp) %>%
  summarize(overlap=max(overlap, NA.rm=T),
            overlap_bien=max(overlap_bien, NA.rm=T))%>%
  rename("gensp"="FIA.sp")

ranges$overlap_max <- ifelse(ranges$overlap_bien>ranges$overlap, ranges$overlap_bien, ranges$overlap)

SWBs <- merge(SWBs, ranges, by='gensp')
## save data
write.csv(SWBs, paste0(data.path, "all_mismatch_and_traits.csv"), row.names = F)


## Summary statistics----
meanmismatch_spp <- group_by(SWBs, gensp)%>%
  summarize(mean_mis_p = mean(mismatch_p, na.rm = T),
            sd_mis_p = sd(mismatch_p, na.rm=T),
            min_mis_p = min(mismatch_p, na.rm = T),
            max_mis_p = max(mismatch_p, na.rm = T),
            sign = ifelse(min_mis_p>0&max_mis_p>0, "allpos",
                          ifelse(min_mis_p<0&max_mis_p<0, "allneg", "mix")),
            no_vitals = length(which(!is.na(mismatch_p))))

mismatch_summary <- group_by(meanmismatch_spp, no_vitals) %>%
  summarize(pos = length(which(sign=="allpos")),
            neg = length(which(sign=="allneg")),
            mix = length(which(sign=="mix")))



meanmismatch_vr <- group_by(SWBs, vr)%>%
  summarize(mean_mis_p = mean(mismatch_p, na.rm = T),
            sd_mis_p = sd(mismatch_p, na.rm=T), 
            min_mis_p = min(mismatch_p, na.rm = T),
            max_mis_p = max(mismatch_p, na.rm = T),
            # pvalue = t.test(mismatch_p)$p.value
            med_mis_p = median(mismatch_p, na.rm=T),
            n_pos_mis = length(which(mismatch_p>0)),
            n_neg_mis = length(which(mismatch_p<0)),
            n_zero_mis = length(which(mismatch_p==0))
  )

abs_meanmismatch_vr <- group_by(SWBs, vr)%>%
  summarize(mean_mis_p = mean(abmis_p, na.rm = T),
            sd_mis_p = sd(abmis_p, na.rm=T), 
            min_mis_p = min(abmis_p, na.rm = T),
            max_mis_p = max(abmis_p, na.rm = T),
            # pvalue = t.test(mismatch_p)$p.value
            med_mis_p = median(abmis_p, na.rm=T)
  )


meancost_vr <- group_by(SWBs, vr)%>%
  summarize(mean_cost_p = mean(cost_p, na.rm = T),
            sd_cost_p = sd(cost_p, na.rm=T), 
            min_cost_p = min(cost_p, na.rm = T),
            max_cost_p = max(cost_p, na.rm = T)
  )

meancost_percent_vr <- group_by(SWBs, vr)%>%
  summarize(mean_cost_p = mean(scaled_cost_p, na.rm = T),
            med_cost_p = median(scaled_cost_p, na.rm = T),
            # sd_cost_p = sd(scaled_cost_p, na.rm=T), 
            min_cost_p = min(scaled_cost_p, na.rm = T),
            max_cost_p = max(scaled_cost_p, na.rm = T)
  )

## Export summary tables----
write.csv(meanmismatch_vr, paste0(table.path, "mismatch_vr.csv"), row.names=F)
write.csv(meanmismatch_spp, paste0(table.path, "mismatch_spp.csv"), row.names=F)
write.csv(abs_meanmismatch_vr, paste0(table.path, "abs_mismatch_vr.csv"), row.names=F)
write.csv(meancost_vr, paste0(table.path, "cost_vr.csv"), row.names=F)
write.csv(meancost_percent_vr, paste0(table.path, "percent_cost_vr.csv"), row.names=F)


