## 5_Beta summaries
## Summaries of regression coefficients joined to species traits.
## Author: Teresa Bohner
## Project: where western trees grow vs. where they thrive (see README)

library(tidyverse)
library(data.table)
library(reshape2)

data.path <- "Processed Data/Analysis Ready Data/"
table.path <- "Output tables/"

## Regression coefficients----
nbetas <- read.csv(paste0(data.path, "abun_betas_swb_temp.csv"))
nbetas <- filter(nbetas, gensp!= "Pseudotsuga macrocarpa")
nbetas$signif95 <- ifelse(nbetas$X2.5.>0&nbetas$X97.5.>0,1,
                          ifelse(nbetas$X2.5.<0&nbetas$X97.5.<0,1,0))
nbetas$signif90 <- ifelse(nbetas$X5.>0&nbetas$X95.>0,1,
                          ifelse(nbetas$X5.<0&nbetas$X95.<0,1,0))
nbetas$signif85 <- ifelse(nbetas$X7.5.>0&nbetas$X92.5.>0,1,
                          ifelse(nbetas$X5.<0&nbetas$X95.<0,1,0))
nbetas$signif <- as.factor(ifelse(nbetas$signif95==1, 95, 
                                  ifelse(nbetas$signif90==1, 90, 0)))
nbetas <- group_by(nbetas, param) %>%
  mutate(medparam=median(X50.[which(signif90==1)]))
nbetas$vr <- "0abun"

nbetas <- nbetas %>%
  mutate(
    G = str_sub(string = gensp, start = 0, end = 1),
    sp = str_split_fixed(string = gensp, pattern= " ", 2)[,2],
    G.sp = paste(G, sp, sep = ". "))

paramnames <- data.frame(c("b_swb", "b2_swb"),
                         c("b1","b2"),
                         c("1° SWB", "2° SWB"))
names(paramnames) <- c("param", "betas", "predictor")
paramnames <- paramnames[order(paramnames$betas),] 
paramnames2 <- paramnames[1:6,]

nbetas <- merge(nbetas, paramnames, by="param")
nbetas$sign <- ifelse(nbetas$X50.<0, "neg", 
                      ifelse(nbetas$X50.>0, "pos", "zero"))

nbetas2 <- data.table::dcast(data.table(nbetas), spid+gensp+vr~param, value.var=c("X50.", "signif90"))

nsignif <- group_by(nbetas, vr, param) %>%
  summarize(medparam = median(X50.[which(signif90==1)]),
            maxparam = max(X50.[which(signif90==1)]),
            minparam = min(X50.[which(signif90==1)]),
            nsignif = length(which(signif90==1)),
            npos = length(which(signif90==1&X50.>0)),
            nneg = length(which(signif90==1&X50.<0)))


gbetas <- read.csv(paste0(data.path, "grow_betas_swb_temp.csv"))
gbetas$signif95 <- ifelse(gbetas$X2.5.>0&gbetas$X97.5.>0,1,
                          ifelse(gbetas$X2.5.<0&gbetas$X97.5.<0,1,0))
gbetas$signif90 <- ifelse(gbetas$X5.>0&gbetas$X95.>0,1,
                          ifelse(gbetas$X5.<0&gbetas$X95.<0,1,0))
gbetas$signif85 <- ifelse(gbetas$X7.5.>0&gbetas$X92.5.>0,1,
                          ifelse(gbetas$X5.<0&gbetas$X95.<0,1,0))
gbetas$signif <- as.factor(ifelse(gbetas$signif95==1, 95, 
                                  ifelse(gbetas$signif90==1, 90, 0)))
gbetas <- group_by(gbetas, param) %>%
  mutate(medparam=median(X50.[which(signif90==1)]))
gbetas$vr <- "2growth"

gbetas <- gbetas %>%
  mutate(
    G = str_sub(string = gensp, start = 0, end = 1),
    sp = str_split_fixed(string = gensp, pattern= " ", 2)[,2],
    G.sp = paste(G, sp, sep = ". "))

paramnames <- data.frame(c("b_adcount", "b_ba", "b_dia", "b_swb", "b_temp", "b2_dia", "b2_swb", "b2_temp"),
                         c("b1","b2","b7","b3","b5","b8","b4","b6"),
                         c("Stem Count", "Mean (BA)", "1° Diameter", "1° SWB", "1° Temperature",
                           "2° Diameter", "2° SWB", "2° Temperature"))
names(paramnames) <- c("param", "betas", "predictor")
paramnames <- paramnames[order(paramnames$betas),] 
paramnames2 <- paramnames[1:6,]

gbetas <- merge(gbetas, paramnames, by="param")
gbetas$sign <- ifelse(gbetas$X50.<0, "neg", 
                      ifelse(gbetas$X50.>0, "pos", "zero"))

gbetas2 <- data.table::dcast(data.table(gbetas), spid+gensp+vr~param, value.var=c("X50.", "signif90"))

gsignif <- group_by(gbetas, vr, param) %>%
  summarize(medparam = median(X50.[which(signif90==1)]),
            maxparam = max(X50.[which(signif90==1)]),
            minparam = min(X50.[which(signif90==1)]),
            nsignif = length(which(signif90==1)),
            npos = length(which(signif90==1&X50.>0)),
            nneg = length(which(signif90==1&X50.<0)))


sbetas <- read.csv(paste0(data.path, "surv_betas_swb_temp.csv"))
sbetas$signif95 <- ifelse(sbetas$X2.5.>0&sbetas$X97.5.>0,1,
                          ifelse(sbetas$X2.5.<0&sbetas$X97.5.<0,1,0))
sbetas$signif90 <- ifelse(sbetas$X5.>0&sbetas$X95.>0,1,
                          ifelse(sbetas$X5.<0&sbetas$X95.<0,1,0))
sbetas$signif <- as.factor(ifelse(sbetas$signif95==1, 95, 
                                  ifelse(sbetas$signif90==1, 90, 0)))

sbetas <- group_by(sbetas, param) %>%
  mutate(medparam=median(X50.[which(signif90==1)]))
sbetas$vr <- "1surv"
sbetas <- sbetas %>%
  mutate(
    G = str_sub(string = gensp, start = 0, end = 1),
    sp = str_split_fixed(string = gensp, pattern= " ", 2)[,2],
    G.sp = paste(G, sp, sep = ". "))
sbetas <- merge(sbetas, paramnames, by="param") 
sbetas$sign <- ifelse(sbetas$X50.<0, "neg", 
                      ifelse(sbetas$X50.>0, "pos", "zero"))

sbetas2 <- data.table::dcast(data.table(sbetas), spid+gensp+vr~param, value.var=c("X50.", "signif90"))

ssignif <- group_by(sbetas, vr, param) %>%
  summarize(medparam = median(X50.[which(signif90==1)]),
            maxparam = max(X50.[which(signif90==1)]),
            minparam = min(X50.[which(signif90==1)]),
            nsignif = length(which(signif90==1)),
            npos = length(which(signif90==1&X50.>0)),
            nneg = length(which(signif90==1&X50.<0)))


rbetas <- read.csv(paste0(data.path, "rec_betas_swb_temp.csv"))
rbetas$signif95 <- ifelse(rbetas$X2.5.>0&rbetas$X97.5.>0,1,
                          ifelse(rbetas$X2.5.<0&rbetas$X97.5.<0,1,0))
rbetas$signif90 <- ifelse(rbetas$X5.>0&rbetas$X95.>0,1,
                          ifelse(rbetas$X5.<0&rbetas$X95.<0,1,0))
rbetas$signif <- as.factor(ifelse(rbetas$signif95==1, 95, 
                                  ifelse(rbetas$signif90==1, 90, 0)))
rbetas <- group_by(rbetas, param) %>%
  mutate(medparam=median(X50.[which(signif90==1)]))
rbetas$vr <- "3rec"
rbetas <- rbetas %>%
  mutate(
    G = str_sub(string = gensp, start = 0, end = 1),
    sp = str_split_fixed(string = gensp, pattern= " ", 2)[,2],
    G.sp = paste(G, sp, sep = ". "))
rbetas <- merge(rbetas, paramnames2, by="param")
rbetas$sign <- ifelse(rbetas$X50.<0, "neg", 
                      ifelse(rbetas$X50.>0, "pos", "zero"))

rbetas2 <- data.table::dcast(data.table(rbetas), spid+gensp+vr~param, value.var=c("X50.", "signif90"))

rsignif <- group_by(rbetas, vr, param) %>%
  summarize(medparam = median(X50.[which(signif90==1)]),
            maxparam = max(X50.[which(signif90==1)]),
            minparam = min(X50.[which(signif90==1)]),
            nsignif = length(which(signif90==1)),
            npos = length(which(signif90==1&X50.>0)),
            nneg = length(which(signif90==1&X50.<0)))

swbbetas <- filter(ssignif, param%in%c("b_swb", "b2_swb")) %>%
  rbind(filter(nsignif, param%in%c("b_swb", "b2_swb"))) %>%
  rbind(filter(gsignif, param%in%c("b_swb", "b2_swb"))) %>%
  rbind(filter(rsignif, param%in%c("b_swb", "b2_swb")))


allbetas <- bind_rows(list(gbetas, sbetas, rbetas, nbetas))

allbetas2 <- rbind(dplyr::select(sbetas2, -"spid"),
                   dplyr::select(gbetas2, -"spid")) %>%
  rbind(dplyr::select(rbetas2, -"spid"), fill=T) 

allbetas2$signif_90_swb <- ifelse(allbetas2$signif90_b_swb==1|allbetas2$signif90_b2_swb==1,1,0)
allbetas2$signif_90_temp <- ifelse(allbetas2$signif90_b_temp==1|allbetas2$signif90_b2_temp==1,1,0)
allbetas2$signif_90_dia <- ifelse(allbetas2$signif90_b_dia==1|allbetas2$signif90_b2_dia==1,1,0)
allbetas2$signif_90_adcount <- allbetas2$signif90_b_adcount
allbetas2$signif_90_ba <- allbetas2$signif90_b_ba

spp_sig <- group_by(allbetas2, gensp) %>%
  summarize(totvr = length(vr),
            totvr_dia = length(vr)-1,
            n_signif_swb = sum(signif_90_swb, na.rm=T),
            prop_signif_swb = n_signif_swb/totvr,

            n_signif_temp = sum(signif_90_temp, na.rm=T),
            prop_signif_temp = n_signif_temp/totvr,

            n_signif_dia = sum(signif_90_dia, na.rm=T),
            prop_signif_dia = n_signif_dia/totvr,

            n_signif_adcount = sum(signif_90_adcount, na.rm=T),
            prop_signif_adcount = n_signif_adcount/totvr,

            n_signif_ba = sum(signif_90_ba, na.rm=T),
            prop_signif_ba = n_signif_ba/totvr)

nbetas2$signif_90_swb <- ifelse(nbetas2$signif90_b_swb==1|nbetas2$signif90_b2_swb==1,1,0)
sum(nbetas2$signif_90_swb)

spp_sig$is_signif_swb = ifelse(spp_sig$n_signif_swb>0, 1, 0)
spp_sig$is_signif_temp = ifelse(spp_sig$n_signif_temp>0, 1, 0)
spp_sig$is_signif_dia = ifelse(spp_sig$n_signif_dia>0, 1, 0)
spp_sig$is_signif_adcount = ifelse(spp_sig$n_signif_adcount>0, 1, 0)
spp_sig$is_signif_ba = ifelse(spp_sig$n_signif_ba>0, 1, 0)

sum(spp_sig$is_signif_swb)
sum(spp_sig$prop_signif_swb[which(spp_sig$prop_signif_swb==1)])
sum(spp_sig$is_signif_temp)
sum(spp_sig$is_signif_dia)
sum(spp_sig$is_signif_adcount)
sum(spp_sig$is_signif_ba)

vr_sig <- group_by(allbetas2, vr) %>%
  summarize(totspp = length(signif90_b_swb),
            n_signif_swb = sum(signif_90_swb),
            prop_signif_swb = n_signif_swb/totspp,
            n_signif_temp = sum(signif_90_temp),
            prop_signif_temp = n_signif_temp/totspp,
            n_signif_dia = sum(signif_90_dia),
            prop_signif_dia = n_signif_dia/totspp,
            n_signif_adcount = sum(signif_90_adcount),
            prop_signif_adcount = n_signif_adcount/totspp,
            n_signif_ba = sum(signif_90_ba),
            prop_signif_ba = n_signif_ba/totspp)

stand_sig <- group_by(allbetas2, vr) %>%
  summarize(n_sig_count = length(which(signif_90_adcount==1)),
            n_sig_ba = length(which(signif_90_ba==1)),
            npos_sig_count = length(which(signif90_b_adcount==1&X50._b_adcount>0)),
            pos_per_count = npos_sig_count/n_sig_count,
            nneg_sig_count = length(which(signif90_b_adcount==1&X50._b_adcount<0)),
            neg_per_count = nneg_sig_count/n_sig_count,
            npos_sig_ba = length(which(signif90_b_ba==1&X50._b_ba>0)),
            pos_per_ba = npos_sig_ba/n_sig_ba,
            nneg_sig_ba = length(which(signif90_b_ba==1&X50._b_ba<0)),
            neg_per_ba = nneg_sig_ba/n_sig_ba)

allrates <- group_by(stand_sig) %>%
  summarize(pos_per_count = mean(pos_per_count),
            neg_per_count = mean(neg_per_count),
            pos_per_ba = mean(pos_per_ba),
            neg_per_ba = mean(neg_per_ba))

## Export regression coefficients----
write.csv(swbbetas, paste0(table.path, "swb_betas.csv"), row.names=F)


