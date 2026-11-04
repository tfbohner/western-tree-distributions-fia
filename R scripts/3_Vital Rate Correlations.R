## 3_Vital Rate Correlations
## Plot- and species-level correlations among average vital rates.
## Author: Teresa Bohner
## Project: where western trees grow vs. where they thrive (see README)

data.path <- "Processed Data/Analysis Ready Data/"

## Save data for analyses ----
occurrence <- read.csv(paste0(data.path, "occurrence.csv"))
abundance <- read.csv(paste0(data.path, "abundance.csv"))
trees <- read.csv(paste0(data.path, "trees.csv"))
surv <- read.csv(paste0(data.path, "surv.csv"))
recruits <- read.csv(paste0(data.path, "recruits.csv"))
spptable <- read.csv(paste0(data.path, "spptable.csv"))


avgdens <- abundance %>%
  group_by(plot.uniqueID, gensp) %>%
  summarize(avg_dens = mean(totabun, na.rm=T),
            avg_ba = mean(totba, na.rm=T)) %>%
  filter(gensp!="Pinus coulteri")
N_spp <- length(unique(avgdens$gensp))

trees$rel_bai <- trees$bai/trees$ba
avggrowth <- trees %>%
  group_by(plot.uniqueID, gensp) %>%
  summarize(avg_bai = mean(bai, na.rm=T),
            avg_rel_bai = mean(rel_bai, na.rm=T))

avgsurv <- surv %>%
  group_by(plot.uniqueID, gensp, spid) %>%
  summarize(avg_surv = mean(surv, na.rm=T))
N_spp <- length(unique(avgsurv$gensp))

recruits$pc_rec <- recruits$N.sdlgs/recruits$adcount.s  
avgrec <- recruits %>%
  group_by(plot.uniqueID, gensp) %>%
  summarize(avg_rec = mean(N.sdlgs, na.rm=T),
            avg_pc_rec = mean(pc_rec, na.rm=T))

all_rates <- avgdens %>%
  left_join(avggrowth, by=c("plot.uniqueID", "gensp")) %>%
  left_join(avgsurv, by=c("plot.uniqueID", "gensp")) %>%
  left_join(avgrec, by=c("plot.uniqueID", "gensp")) 

library(Hmisc)
cor_array <- array(data=NA, dim = c(4, 4, length(unique(all_rates$gensp))))
p_array <- array(data=NA, dim = c(4, 4, length(unique(all_rates$gensp))))
for(i in 1:N_spp) {
  data <- filter(all_rates, spid==i) %>% 
    dplyr::select(c(spid, gensp, plot.uniqueID, avg_dens, avg_bai, avg_surv, avg_rec))
  cor_array[,,i] <- rcorr(as.matrix(data[, 4:7]))$r
  p_array[,,i] <- rcorr(as.matrix(data[, 4:7]))$P
}

dimnames(cor_array) <- list(c("avg_dens", "avg_bai", "avg_surv", "avg_rec"), 
                            c("avg_dens", "avg_bai", "avg_surv", "avg_rec"),
                            unique(all_rates$gensp))
dimnames(p_array) <- list(c("avg_dens", "avg_bai", "avg_surv", "avg_rec"), 
                          c("avg_dens", "avg_bai", "avg_surv", "avg_rec"),
                          unique(all_rates$gensp))
for(i in unique(all_rates$gensp)) {
  data <- filter(all_rates, gensp==i) %>% 
    dplyr::select(c(spid, gensp, plot.uniqueID, avg_dens, avg_bai, avg_surv, avg_rec))
  cor_array[,,i] <- rcorr(as.matrix(data[, 4:7]))$r
  p_array[,,i] <- rcorr(as.matrix(data[, 4:7]))$P
}

detach("package:Hmisc", unload = T)


## some problems after here...addreess
library(corrplot)
pdf("Output pdfs/corrplots.pdf")
for(i in unique(all_rates$gensp)){
  data <- cor_array[,,i]
  rownames(data) <- c("Density", "Growth", "Survival", "Recruitment")
  colnames(data) <- rownames(data)
  corrplot.mixed(corr=data, upper="ellipse", title = spptable.n$gensp[which(spptable.n$gensp==i)], 
                 mar=c(0,0,2,0), na.label = "NA", tl.col = "black", lower.col="black")
}
dev.off()

spp_cor <- as.data.frame(cbind(NA, cor_array[1,2,], cor_array[1,3,], cor_array[1,4,],
                               NA, cor_array[2,3,], cor_array[2,4,],
                               NA, cor_array[3,4,],
                               NA))
names(spp_cor) <- c("dens-dens", "dens-growth", "dens-surv", "dens-rec", 
                    "growth-growth","growth-surv", "growth-rec", 
                    "surv-surv", "surv-rec",
                    "rec-rec")
spp_cor$spid <- c(1:59)
spp_cor <- left_join(spp_cor, spptable.n, by="spid")

spp_cor2 <- spp_cor %>%
  gather(key="pair", value="cor", -c(spid, gensp)) 
spp_cor2$row <- ifelse(str_detect(spp_cor2$pair, "dens-"), "0density",
                       ifelse(str_detect(spp_cor2$pair, "growth-"), "1growth",
                              ifelse(str_detect(spp_cor2$pair, "surv-"), "2survival", 
                                     ifelse(str_detect(spp_cor2$pair, "rec-"), "3recruitment", NA))))
spp_cor2$col <- ifelse(str_detect(spp_cor2$pair, "-dens"), "0density",
                       ifelse(str_detect(spp_cor2$pair, "-growth"), "1growth",
                              ifelse(str_detect(spp_cor2$pair, "-surv"), "2survival",
                                     ifelse(str_detect(spp_cor2$pair, "-rec"), "3recruitment", NA))))
spp_cor2$dum <- 1
labels <- c("0density" = "Density", "1growth" = "Growth", "2survival" = "Survival", "3recruitment" = "Recruitment")
pdf("Output pdfs/corrplot_summary.pdf")
ggplot(spp_cor2, aes(x=dum, y=cor)) +
  geom_violin(draw_quantiles=0.5) +
  xlim(c(0,2)) +
  theme(axis.title.y = element_text(size = 14),
        axis.text.y = element_text(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1, fill=NA),
        strip.background = element_rect(fill=NA),
        strip.text.x = element_text(angle=90, hjust = 0, size=12),
        strip.text.y = element_text(angle=0, hjust = 0, size=12)) +
  ylab(label = "correlation coefficient") +
  geom_hline(yintercept = 0, linetype="dashed") +
  facet_grid(row~col, labeller=labeller(row = labels, col = labels)) 
dev.off()


