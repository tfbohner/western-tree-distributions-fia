## Author: Teresa Bohner
## Purpose: Chapter 1 analyses: Data preparation for analyses
## Date Modified: 2-24-2019

## change tree.table 25 to be based off of abundance

library(tidyverse)
rm(list=ls(all=TRUE))


## read in processed FIA data

occurrence <- read.csv("Processed Data/Allplot.Occurdata.csv")
occurrence <- rename(occurrence, plot.uniqueID = fia.ID)
occurrence <- filter(occurrence, STDORGCD!=1)
abundance <- read.csv("Processed Data/Allspecies.abundance.csv")
abundance <- rename(abundance, plot.uniqueID = fia.ID)
abundance <- filter(abundance, STDORGCD!=1)
trees <- read.csv("Processed Data/All.Treedata.csv")
trees <- rename(trees, plot.uniqueID = fia.ID)
trees <- filter(trees, STDORGCD!=1)
surv <- trees
sdlg.data <- read.table("Processed Data/Seedling.data.csv",header=TRUE, sep=',') 
sdlg.data <- rename(sdlg.data, plot.uniqueID = fia.ID)
sdlg.data <- filter(sdlg.data, STDORGCD!=1)

## Get new plot.table for plot characteristics ----
plot.table <- group_by(occurrence, plot.uniqueID) %>%
  summarise(Lat = mean(LAT)
            ,Long = mean(LON)
            ,elev = mean(ELEV)
            ,ptcoeff=mean(ptcoeff)
            ,temp=mean(temp)
            ,origin=first(STDORGCD)
            ,physio_cd=first(PHYSCLCD)
            ,physio_cat = first(site_water))
length(which(plot.table$origin==1))
## 2667/21277 plots have artificial regeneration. Remove?

## Jeff does some plotting here from another file
plot(plot.table$ptcoeff~plot.table$Lat)
plot(plot.table$ptcoeff~plot.table$elev)
plot(plot.table$temp~plot.table$elev)
boxplot(plot.table$ptcoeff~plot.table$physio_cat, ylab='SWB', xlab='Physiographic code from FIA')

# trees <- merge(trees, wood[,c(3,4)], by.x = "gensp", by.y = "Binomial")
# names(trees)[43]<-"wood_dens"
## ba to basal area increment (still in m^2)
trees <- rename(trees, bai=BA.grwth)
sdlg.data$recPC <- sdlg.data$N.sdlgs/(sdlg.data$adcount.s + 1)

##incorporate single visit plots into survival
surv$surv <- ifelse(is.na(surv$surv)&surv$repmeas==1&surv$status==1, 1,
                    ifelse(is.na(surv$surv)&surv$repmeas==1&surv$status==0, 0, surv$surv))


## trees gets reduced from 1.1 million records to 240k in this subset below.  Why?
dim(trees)
dim(trees[!is.na(trees$bai),])  # bai is the big reduction, 
head(trees[is.na(trees$bai),])  # because All.Treedata.csv includes species without a baNext (single measurements) ??

## remove NAs
occurrence <- subset(occurrence, !is.na(pres)&!is.na(ptcoeff))
abundance <- subset(abundance, !is.na(liveabun)&!is.na(ptcoeff))
trees <- subset(trees, !is.na(bai)&!is.na(ptcoeff)&!is.na(adcount)&!is.na(avgba))
surv <- filter(surv, is.na(ptcoeff)==F, is.na(dia)==F, is.na(surv)==F, is.na(adcount)==F, is.na(avgba)==F)
R0 <- droplevels(filter(sdlg.data, is.na(adcount)==F, is.na(avgba)==F, is.na(ptcoeff)==F))

occurrence <- subset(occurrence, !is.na(pres)&!is.na(ptcoeff)&!is.na(temp))
abundance <- subset(abundance, !is.na(liveabun)&!is.na(ptcoeff)&!is.na(temp))
trees <- subset(trees, !is.na(bai)&!is.na(ptcoeff)&!is.na(adcount)&!is.na(avgba)&!is.na(temp))
surv <- filter(surv, is.na(ptcoeff)==F, is.na(dia)==F, is.na(surv)==F, is.na(adcount)==F, is.na(avgba)==F, is.na(temp)==F)
R0 <- droplevels(filter(sdlg.data, is.na(adcount)==F, is.na(avgba)==F, is.na(ptcoeff)==F,is.na(temp)==F))

trees <- subset(trees, alive==1&surv==1)

## Get species summaries -----
temp <- group_by(abundance, gensp)
tree.table1 <- summarise(temp,
                         Num.plots = length(which(totabun>0)))
sum(tree.table1$Num.plots)
temp <- group_by(occurrence, gensp)
tree.table2 <- summarise(temp,
                         Num.plots = sum(pres, na.rm=T))
sum(tree.table2$Num.plots)

temp <- group_by(trees, gensp)
tree.table3 <- summarise(temp,
                         Num.trees = length(unique(tree.uniqueID))
                         ,Num.plots = length(unique(plot.uniqueID)))
temp <- group_by(surv, gensp)
tree.table4 <- summarise(temp,
                         Num.trees = length(unique(tree.uniqueID))
                         ,Num.plots = length(unique(plot.uniqueID)))
tree.table5 <- group_by(R0, gensp) %>%
  summarise(Num.plots = sum(pres, na.rm=T))

# tree.table <- as.data.frame(tree.table3)
# tree.table$mean.density <- round(tree.table$Num.trees / tree.table$Num.plots, 1)
# tree.table <- tree.table[order(tree.table$Num.trees, decreasing=TRUE),]
# head(tree.table); dim(tree.table)
## do we want to take out any spp.??
## omit spp without distribution models
species.table25.n <- droplevels(subset(tree.table1, tree.table1$Num.plots>=25))
# species.table25.n1 <- species.table25.n[-grep("spp.",species.table25.n$gensp),]
species.table25.n <- species.table25.n[-grep("Tree",species.table25.n$gensp),]

dropped <- droplevels(subset(tree.table1, tree.table1$Num.plots<25))
dropped <- droplevels(tree.table1$gensp[which(tree.table1$Num.plots<25)])

gentree <- droplevels(tree.table1$gensp[grep("Tree",tree.table1$gensp)])
spp. <- droplevels(tree.table1$gensp[grep("spp.",tree.table1$gensp)])

spplist <- unique(species.table25.n$gensp)

species.table25.p <- droplevels(subset(tree.table2, tree.table2$Num.plots>=25)) %>%
  filter(gensp%in%spplist)
species.table25.g <- droplevels(subset(tree.table3, tree.table3$Num.plots>=25)) %>%
  filter(gensp%in%spplist)
species.table25.s <- droplevels(subset(tree.table4, tree.table4$Num.plots>=25)) %>%
  filter(gensp%in%spplist)
species.table25.r <- droplevels(subset(tree.table5, tree.table5$Num.plots>=25)) %>%
  filter(gensp%in%spplist)


## omit spp without enough data
spplist <- unique(species.table25.n$gensp)
abundance <- droplevels(filter(abundance, gensp%in%spplist))
abundance$spid <- as.numeric(abundance$gensp)
spptable.n <- summarize(group_by(abundance, spid), gensp=as.character(first(gensp)))

spplist <- unique(species.table25.p$gensp)
occurrence <- droplevels(filter(occurrence, gensp%in%spplist))
occurrence$spid <- as.numeric(occurrence$gensp)
spptable.p <- summarize(group_by(occurrence, spid), gensp=as.character(first(gensp)))

spplist <- unique(species.table25.g$gensp)
trees <- droplevels(filter(trees, gensp%in%spplist))
trees$spid <- as.numeric(trees$gensp)
spptable.g <- summarize(group_by(trees, spid), gensp=as.character(first(gensp)))

spplist <- unique(species.table25.s$gensp)
surv <- droplevels(filter(surv, gensp%in%spplist))
surv$spid <- as.numeric(surv$gensp)
spptable.s <- summarize(group_by(surv, spid), gensp=as.character(first(gensp)))

spplist <- unique(species.table25.r$gensp)
recruits <- droplevels(filter(R0, gensp%in%spplist))
recruits$spid <- as.numeric(recruits$gensp)
spptable.r <- summarize(group_by(recruits, spid), gensp=as.character(first(gensp)))

spptable <- merge(spptable.n, spptable.p, by='spid', all=T) %>%
  merge(spptable.g, by='spid', all=T) %>%
  merge(spptable.s, by='spid', all=T) %>%
  merge(spptable.r, by='spid', all=T)

names(spptable) <- c("spid","abundance","occurrence","growth","survival","recruitment")

data.path <- "Processed Data/Analysis Ready Data/"

## Save data for analyses ----
write.csv(occurrence, paste0(data.path, "occurrence.csv"), row.names = F)
write.csv(abundance, paste0(data.path, "abundance.csv"), row.names = F)
write.csv(trees, paste0(data.path, "trees.csv"), row.names = F)
write.csv(surv, paste0(data.path, "surv.csv"), row.names = F)
write.csv(recruits, paste0(data.path, "recruits.csv"), row.names = F)
write.csv(spptable, paste0(data.path, "spptable.csv"), row.names = F)




