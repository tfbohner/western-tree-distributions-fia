## Purpose: Assemble FIA data into workable format, assemble CSVs for analyses

# Major steps & filtering
#  FILTER: Drop records with contradictory data for trees (e.g. multiple azimuth, distances, spcd)
#  get current and next year data for growth and status for each repeat measurement
#  PAIRWISE data
    # lump in 'gensp25' all species in fewer than 25 plots as 'Other'
    # get numbers of adults of each species for each record, to enable pairwise analyses
    # records are 'target' species; in columns are 'neighbor' species
# FILTER: Subset tree dataset to those target species in >= 25 plots ----  
## Assemble seedling data

## Products:

# write.csv(CA.data, file='All.Treedata.csv')   FOR GROWTH AND SURVIVAL ANALYSES
# write.csv(sdlgs.together, file='Seedling.data.csv', row.names=FALSE)   FOR SEEDLING RECRUITMENT ANALYSIS

## FIA subplots = 24ft radius = 7.315 meter radius.   So, subplot = 168 m2 area.   This is approx area of a 13 x 13 meter square.   Small for trees.  
## Plot = 4 x 168 m2 = 672 m2.   =  approx 26 x 26 meter square area


rm(list=ls(all=TRUE))
library(reshape2)
library(dplyr)
library(data.table)
library(gdata)
library(raster)
library(rgdal)

data.path <- "Raw Data/CA 2/"

sp.refs <- read.table(paste(data.path,"REF_SPECIES.CSV",sep=''),header=T, sep=',') # read data file
sp.ref <- subset(sp.refs, select=c(SPCD, GENUS, SPECIES, COMMON_NAME))

catrees.raw <- read.table(paste(data.path,"CA_TREE.csv",sep=''),header=TRUE, sep=',') # read data file
casdlgs.raw <- read.table(paste(data.path,"CA_SEEDLING.csv",sep=''),header=TRUE, sep=',') # read data file
caplots.raw <- read.table(paste(data.path,"CA_PLOT.csv",sep=''),header=TRUE, sep=',') # read data file
catree.grm.raw <- read.table(paste(data.path,"CA_TREE_GRM_ESTN.csv",sep=''),header=TRUE, sep=',') # read data file
cacond.raw <- read.table(paste(data.path,"CA_COND.csv",sep=''),header=TRUE, sep=',') # read data file

data.path <- "Raw Data/OR/"
ortrees.raw <- read.table(paste(data.path,"OR_TREE.csv",sep=''),header=TRUE, sep=',') # read data file
orsdlgs.raw <- read.table(paste(data.path,"OR_SEEDLING.csv",sep=''),header=TRUE, sep=',') # read data file
orplots.raw <- read.table(paste(data.path,"OR_PLOT.csv",sep=''),header=TRUE, sep=',') # read data file
ortree.grm.raw <- read.table(paste(data.path,"OR_TREE_GRM_ESTN.csv",sep=''),header=TRUE, sep=',') # read data file
orcond.raw <- read.table(paste(data.path,"OR_COND.csv",sep=''),header=TRUE, sep=',') # read data file

data.path <- "Raw Data/WA/"

watrees.raw <- read.table(paste(data.path,"WA_TREE.csv",sep=''),header=TRUE, sep=',') # read data file
wasdlgs.raw <- read.table(paste(data.path,"WA_SEEDLING.csv",sep=''),header=TRUE, sep=',') # read data file
waplots.raw <- read.table(paste(data.path,"WA_PLOT.csv",sep=''),header=TRUE, sep=',') # read data file
watree.grm.raw <- read.table(paste(data.path,"WA_TREE_GRM_ESTN.csv",sep=''),header=TRUE, sep=',') # read data file
wacond.raw <- read.table(paste(data.path,"WA_COND.csv",sep=''),header=TRUE, sep=',') # read data file

trees.raw <- rbind(catrees.raw, ortrees.raw, watrees.raw)
sdlgs.raw <- rbind(casdlgs.raw, orsdlgs.raw, wasdlgs.raw)
plots.raw <- rbind(caplots.raw, orplots.raw, waplots.raw)
tree.grm.raw <- rbind(catree.grm.raw, ortree.grm.raw, watree.grm.raw)
cond.raw <- rbind(cacond.raw, orcond.raw, wacond.raw)

ptcoeff.fia <- read.xls("Raw Data/From FIA/alpha_extractvalu.xlsx")
temp.fia <- read.xls("Raw Data/From FIA/temp_vals_realloc.xlsx")
temp.fia$temp_rnd <- ifelse(temp.fia$temp_rnd==-9999, NA, temp.fia$temp_rnd)

# ptcoeff <- raster( "Raw Data/Global Soil-Water Balance/ALPHA/alpha/hdr.adf")
ptcoeff.fia$fia.ID <- with(ptcoeff.fia, paste(STATECD, COUNTYCD, PLOT_FIADB, sep='.'))
temp.fia$fia.ID <- with(temp.fia, paste(STATECD, COUNTYCD, PLOT_FIADB, sep='.'))
climate.fia <- merge(ptcoeff.fia, dplyr::select(temp.fia, c(fia.ID, temp_rnd)), by="fia.ID", all=T)
climate.fia <- rename(climate.fia, ptcoeff=RASTERVALU, temp=temp_rnd)


# t1 <- raster("Raw Data/wc2.0_30s_tavg/wc2.0_30s_tavg_01.tif")
# t2 <- raster("Raw Data/wc2.0_30s_tavg/wc2.0_30s_tavg_02.tif")
# t3 <- raster("Raw Data/wc2.0_30s_tavg/wc2.0_30s_tavg_03.tif")
# t4 <- raster("Raw Data/wc2.0_30s_tavg/wc2.0_30s_tavg_04.tif")
# t5 <- raster("Raw Data/wc2.0_30s_tavg/wc2.0_30s_tavg_05.tif")
# t6 <- raster("Raw Data/wc2.0_30s_tavg/wc2.0_30s_tavg_06.tif")
# t7 <- raster("Raw Data/wc2.0_30s_tavg/wc2.0_30s_tavg_07.tif")
# t8 <- raster("Raw Data/wc2.0_30s_tavg/wc2.0_30s_tavg_08.tif")
# t9 <- raster("Raw Data/wc2.0_30s_tavg/wc2.0_30s_tavg_09.tif")
# t10 <- raster("Raw Data/wc2.0_30s_tavg/wc2.0_30s_tavg_10.tif")
# t11 <- raster("Raw Data/wc2.0_30s_tavg/wc2.0_30s_tavg_11.tif")
# t12 <- raster("Raw Data/wc2.0_30s_tavg/wc2.0_30s_tavg_12.tif")
# 
# temp <- stack(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12)
# 
# states <- readOGR("/Users/diezlab/Google Drive File Stream/My Drive/GIS/shapefiles/cb_2016_us_state_500k/cb_2016_us_state_500k.shp")
# states <- readOGR("/Users/teresabohner/Google Drive File Stream/My Drive/GIS/shapefiles/cb_2016_us_state_500k/cb_2016_us_state_500k.shp")
# states1 <-  states[states$NAME=="California"|states$NAME=="Oregon"|states$NAME=="Washington",]
# states.f <- fortify(states1, regions="STUSPS")
# 
# temp1 <- crop(temp, states1)
# avgtemp <- mean(temp1)

# writeRaster(avgtemp, "Raw Data/wc2.0_30s_tavg/wc2.0_30s_tavg_annual_PNW.tif")

# plots.raw$LAT <- signif(plots.raw$LAT,3)
# plots.raw$LON <- signif(plots.raw$LON,4)

trees.raw$plot.uniqueID <- with(trees.raw, paste(STATECD, UNITCD, COUNTYCD, PLOT, sep='.'))
trees.raw$fia.ID <- with(trees.raw, paste(STATECD, COUNTYCD, PLOT, sep='.'))
sdlgs.raw$plot.uniqueID <- with(sdlgs.raw, paste(STATECD, UNITCD, COUNTYCD, PLOT, sep='.'))
sdlgs.raw$fia.ID <- with(sdlgs.raw, paste(STATECD, COUNTYCD, PLOT, sep='.'))
plots.raw$plot.uniqueID <- with(plots.raw, paste(STATECD, UNITCD, COUNTYCD, PLOT, sep='.'))
plots.raw$fia.ID <- with(plots.raw, paste(STATECD, COUNTYCD, PLOT, sep='.'))
cond.raw$plot.uniqueID <- with(cond.raw, paste(STATECD, UNITCD, COUNTYCD, PLOT, sep='.'))
cond.raw$fia.ID <- with(cond.raw, paste(STATECD, COUNTYCD, PLOT, sep='.'))

trees.sub <- subset(trees.raw, select=c(plot.uniqueID,SPCD,PLOT,SUBP,PLT_CN,INVYR,CONDID,TREE,STATUSCD,DIA,ACTUALHT,HT,AZIMUTH,DIST))
trees.sub2 <- subset(trees.raw, select=c(fia.ID,SPCD,PLOT,SUBP,PLT_CN,INVYR,CONDID,TREE,STATUSCD,DIA,ACTUALHT,HT,AZIMUTH,DIST))
trees.sub3 <- subset(trees.raw, select=c(fia.ID,STATECD,COUNTYCD,PLOT))
sdlgs.sub <- subset(sdlgs.raw, select=c(plot.uniqueID,SPCD,PLOT,SUBP,PLT_CN,CN,INVYR,TREECOUNT,TREECOUNT_CALC))
sdlgs.sub2 <- subset(sdlgs.raw, select=c(fia.ID,SPCD,PLOT,SUBP,PLT_CN,CN,INVYR,TREECOUNT,TREECOUNT_CALC))
plots.sub <- subset(plots.raw, select=c(plot.uniqueID,LAT,CN,LON,ELEV,MEASYEAR,MANUAL,ECOSUBCD))
plots.sub2 <- subset(plots.raw, select=c(fia.ID,LAT,CN,LON,ELEV,MEASYEAR,MANUAL,ECOSUBCD))
plots.sub3 <- subset(plots.raw, select=c(fia.ID,STATECD,COUNTYCD,PLOT,CN,MEASYEAR,MANUAL))
cond.sub <- subset(cond.raw, select=c(plot.uniqueID, CN, PLT_CN, INVYR, STDAGE, CONDPROP_UNADJ, TRTCD1, TRTCD2, TRTCD3, STDORGCD, PHYSCLCD, DSTRBCD1, AFFORESTATION_CD, PREV_AFFORESTATION_CD))
cond.sub2 <- subset(cond.raw, select=c(fia.ID, CN, PLT_CN, INVYR, STDAGE, CONDPROP_UNADJ, TRTCD1, TRTCD2, TRTCD3, STDORGCD, PHYSCLCD, DSTRBCD1, AFFORESTATION_CD, PREV_AFFORESTATION_CD))
# plots.for.fia.people <- subset(plots.raw, select=c(plot.uniqueID, STATECD, INVYR, UNITCD, COUNTYCD, PLOT, LAT, LON))
# plots.ptcoeff <- merge(plots.sub2, ptcoeff.fia, by="fia.ID", all.x=T)
plots.climate <- merge(plots.sub2, climate.fia, by="fia.ID", all.x=T)

# plots.for.john <- subset(plots.sub3, !fia.ID%in%ptcoeff.fia$fia.ID)
# plots.for.john2 <- subset(trees.sub3, !fia.ID%in%ptcoeff.fia$fia.ID)
# plots.for.john2 <- distinct(plots.for.john2, fia.ID, .keep_all=TRUE)
# 
# write.csv(plots.for.john[2:7], "plots with no match.csv", row.names = F)
# write.csv(plots.for.john[2:4], "plots with tree records no match.csv", row.names=F)
# plots.no.ptcoeff.fia <- subset(plots.ptcoeff, is.na(RASTERVALU))


# plots.ptcoeff <- merge(plots.sub2, ptcoeff.fia, by="fia.ID")
# plots.ptcoeff <- rename(plots.ptcoeff, ptcoeff=RASTERVALU)
# xy_co <- plots.ptcoeff[,c(4,2)]
# plots.ptcoeff$temp <- extract(avgtemp, xy_co) 

trees.sub <- subset(trees.sub, STATUSCD>0)   # 1=Live tree;  2=Dead tree;  3=Removed tree    0=not in current survey
trees.sub2 <- subset(trees.sub2, STATUSCD>0)

# xy_co <- plots.sub[,c(4,2)]
# plots.sub$ptcoeff <- extract(ptcoeff, xy_co)
# plots.sub$temp <- extract(avgtemp, xy_co)


# xy_co <- plots.for.fia.people[,c(8,7)]
# plots.for.fia.people$ptcoeff <- extract(ptcoeff, xy_co)

plot.cond.yrs <- merge(plots.climate, cond.sub2, by='fia.ID', all.x=T)
cond.short <- dplyr::select(cond.sub2, c('fia.ID', 'CONDPROP_UNADJ','STDORGCD', 'PHYSCLCD', 'TRTCD1', 'TRTCD2')) %>%
  arrange(desc(CONDPROP_UNADJ, fia.ID)) %>%
  group_by(fia.ID) %>% 
  summarize(STDORGCD = first(STDORGCD),
            PHYSCLCD = first(PHYSCLCD),
            TRTCD1 = first(TRTCD1),
            CONDPROP_UNADJ = first(CONDPROP_UNADJ))
plot.cond <- merge(plots.climate, cond.short, all.x=T)

plot.cond$site_water <- ifelse(is.na(plot.cond$PHYSCLCD), 'NA',
                               ifelse(plot.cond$PHYSCLCD<20, "xeric", 
                                      ifelse(plot.cond$PHYSCLCD <30, 'mesic',
                                             ifelse(plot.cond$PHYSCLCD<40, 'hydric'))))

plots.climate <- plot.cond

trees.sub <- merge(trees.sub2, plots.climate, by.x='PLT_CN', by.y='CN')
trees.sub <- filter(trees.sub, fia.ID.x==fia.ID.y)
trees.sub <- rename(trees.sub, fia.ID=fia.ID.x)
trees.sub$fia.ID.y=NULL

temp <- group_by(trees.sub, fia.ID)
temp1 <- summarise(temp,
                   N.NAs = length(which(is.na(MEASYEAR)==TRUE))
)
table(temp1$N.NAs)

sdlgs.sub <- merge(sdlgs.sub2, plots.climate, by.x='PLT_CN', by.y='CN', all.x=TRUE)
sdlgs.sub <- filter(sdlgs.sub, fia.ID.x==fia.ID.y)
sdlgs.sub <- rename(sdlgs.sub, fia.ID=fia.ID.x)
sdlgs.sub$fia.ID.y=NULL

plots.sub <- subset(plots.climate, duplicated(plots.climate$fia.ID)==FALSE, drop=TRUE)

write.table(plots.sub, "Processed Data/plots.sub.csv", sep=",", eol="\n", quote=FALSE, row.names=F)



######################## Adults----
trees <- merge(trees.sub, sp.ref, by="SPCD")
trees$gensp <- paste(trees$GENUS,trees$SPECIES,sep=' ')
tree.table <- as.data.frame(sort(table(trees$gensp), decreasing=TRUE))
colnames(tree.table) <- c("gensp","N")
tree.table

#  Make unique tree id
trees$tree.uniqueID <- NA
trees$tree.uniqueID <- as.factor(with(trees, paste(SPCD, fia.ID, SUBP, TREE, sep='.')))

#  Order by year
trees <- trees[order(trees$tree.uniqueID, trees$MEASYEAR),]

#  Total visits count column (for subsetting later)
##  dplyr code (faster:)
temp <- group_by(trees, tree.uniqueID)
trees <- mutate(temp, repmeas=length(unique(MEASYEAR)))

##  There is some weirdness in tree data that should, but doesn't stay consistent over the years
##  I found this out after melt/cast would return multiple rows for individuals while we only want 1
##  Problems were with PLT_CN (which I don't think we really need), SPCD, AZIMUTH, and DIST
inds <- group_by(trees, tree.uniqueID)
check <- summarise(inds,
                   spcd = length(unique(SPCD)),
                   plot = length(unique(PLOT)),
                   subp = length(unique(SUBP)),
                   tree = length(unique(TREE)),
                   azimuth = length(unique(AZIMUTH)),
                   d.az = max(AZIMUTH)-min(AZIMUTH),
                   dist = length(unique(DIST)),
                   d.dist = max(DIST)-min(DIST),
                   rm = length(unique(repmeas)))
length(which(check$plot>1))
length(which(check$spcd>1))
length(which(check$subp>1))
length(which(check$tree>1))
length(which(check$azimuth>1))
summary(check$d.az)
length(which(check$dist>1))
summary(check$d.dist)
length(which(check$rm>1))

#  Create list to drop id's from full dataset
#  Dropping:
#  When azimuth difference btw years is greater than 5 and less than 355 (4843)
#  When distance difference btw years is greater than 10 cm
# droptable <- subset(check, spcd>1 | (d.az>5 & d.az<355) | d.dist>10)
# droplist <- as.list(droptable$tree.uniqueID)

#  FILTER 1: Drop records with contradictory data for trees (e.g. multiple azimuth, distances) ----
# trees.tb <- subset(trees, !(tree.uniqueID %in% droplist), drop=TRUE )
# trees.tb <- trees.tb[,c(1:3, 5:22)]

#  We want a single value for azimuth and distance.  For now I will keep min values for ease
#  We might want to change this to an average (weird for azimuth becasue 0 and 360 are the same)
#  I thought it might be best to stick with the first measurement, but I'm not sure how to do this yet
# trees.tb <- ddply(trees.tb, .(tree.uniqueID), mutate,
#                   AZIMUTH = min(AZIMUTH, na.rm=T),
#                   DIST = min(DIST, na.rm=T))
# temp <- group_by(trees, tree.uniqueID)
# trees <- mutate(temp, 
#                    AZIMUTH = min(AZIMUTH, na.rm=TRUE),
#                    DIST = min(DIST, na.rm=TRUE))

#  Add the visit counter variable
#  Turn data.frame into a data.table
trees.dt <- data.table(trees)
#  Add index column (visit)
trees.dt[ , visit := 1:.N, by = "tree.uniqueID"]   # .N gives total number within each group
#  Turn data.table back into a data.frame for melt
trees <- data.frame(trees.dt)

# MOVED FILTER: subset to include only trees with multiple measurement years
#trees.tb <- subset(trees.tb, repmeas>1)  # do this later for the growth data only

trees$STATUSCD[trees$STATUSCD==2] <- 0   # dead trees

#basal area to m^2
trees$ba <- .00007854*((trees$DIA*2.54)^2)

## Presence/absence data ----
temp <- group_by(trees, fia.ID, visit, gensp)
temp <- summarise(temp,
                  liveabun=length(which(STATUSCD==1)),
                  totba=sum(ba, na.rm=T),
                  totabun = length(STATUSCD))


## need the absences
## get all unique plots-visits
temp1 <- group_by(trees, fia.ID, visit)
temp1 <- summarise(temp1, year=first(MEASYEAR))
temp2 <- rbind(temp1, temp1)
spplist <- unique(trees$gensp)

for(i in 1:(length(spplist)-2)){
  temp2<- rbind(temp1, temp2)
}

spp <- matrix(NA, nrow=nrow(temp1), ncol=length(spplist))
for(i in 1:length(spplist)){
  spp[,i] <-rep(spplist[i], nrow(temp1))
}

temp2$gensp <- unname(unmatrix(spp))

## get real data in allplots
allplots <-merge(temp2, temp, by=c("fia.ID","visit","gensp"), all=T)

## abundance to presence absence
allplots$liveabun[is.na(allplots$liveabun)] <- 0
allplots$totabun[is.na(allplots$totabun)] <- 0
allplots$pres <- ifelse(allplots$totabun==0,0,1)
allplots$fia.ID <- as.factor(allplots$fia.ID)
allplots$gensp <- as.factor(allplots$gensp)

length(allplots$pres==1)
length(allplots$totabun>0)

## Data from most recent survey
occurrence <- group_by(allplots, fia.ID, gensp) %>%
  summarize(pres=ceiling(mean(pres)))
occurrence <- merge(occurrence, plots.sub, by="fia.ID")
sum(occurrence$pres)
write.csv(occurrence, "Processed Data/Allplot.Occurdata.csv", row.names = F)

abundance <- subset(allplots, totabun>0) %>%
  group_by(fia.ID, gensp) %>%
  summarize(year = last(year),
            liveabun = last(liveabun),
            totabun = last(totabun),
            totba = last(totba))

abundance <- merge(abundance, plots.sub, by="fia.ID")

write.csv(abundance, "Processed Data/Allspecies.abundance.csv", row.names=F)

# demography data----
# get everything into year yearnext format

trees.1 <- subset(trees, repmeas==1)
trees.1$yearNext <- NA
trees.1$statusNext <- NA
trees.1$diaNext <- NA
trees.1$baNext<- NA
trees.1$htNext <- NA

trees.2 <- subset(trees, repmeas==2)
trees.2$yearNext <- NA
trees.2$statusNext <- NA
trees.2$diaNext <- NA
trees.2$baNext<- NA
trees.2$htNext <- NA

# for(i in 1:length(n.visits)){   
# Changed this - didn't need to CYCLE THROUGH THE 2ND VISIT HERE?   Already ALL NA'S RIGHT?
  trees.2$yearNext[trees.2$visit== 1] <- trees.2$MEASYEAR[trees.2$visit== 2]
  trees.2$statusNext[trees.2$visit== 1] <- trees.2$STATUSCD[trees.2$visit== 2]
  trees.2$diaNext[trees.2$visit== 1] <- trees.2$DIA[trees.2$visit== 2]
  trees.2$baNext[trees.2$visit== 1] <- trees.2$ba[trees.2$visit== 2]
  trees.2$htNext[trees.2$visit== 1] <- trees.2$HT[trees.2$visit== 2]
#}  

  ## This only works when the trees.2 are sorted by tree.uniqueID
  ## Check
  tt <- as.numeric(as.factor(trees.2$tree.uniqueID))
  table(tt - rep(seq(1:max(tt)),each=2))

  
## most recent data max(repmeas)=2
# trees.3 <- subset(trees, repmeas==3)
# trees.3$yearNext <- rep(NA, nrow(trees.3))
# trees.3$statusNext <- rep(NA, nrow(trees.3))
# trees.3$diaNext <- rep(NA, nrow(trees.3))
# trees.3$htNext <- rep(NA, nrow(trees.3))
# 
# n.visits <- c(1,2)
# for(i in 1:length(n.visits)){
#   trees.3$yearNext[trees.3$visit==n.visits[i]] = trees.3$MEASYEAR[trees.3$visit==n.visits[i]+1]
#   trees.3$statusNext[trees.3$visit==n.visits[i]] = trees.3$STATUSCD[trees.3$visit==n.visits[i]+1]
#   trees.3$diaNext[trees.3$visit==n.visits[i]] = trees.3$DIA[trees.3$visit==n.visits[i]+1]
#   trees.3$htNext[trees.3$visit==n.visits[i]] = trees.3$HT[trees.3$visit==n.visits[i]+1]
# }

## Check
# tt <- as.numeric(as.factor(trees.3$tree.uniqueID))
# table(tt - rep(seq(1:max(tt)),each=3))

trees.all <- rbind(trees.1, trees.2) #, trees.3

trees.all$interval <- trees.all$yearNext-trees.all$MEASYEAR


## Annual Basal Area growth rates ----  m^2 / year
trees.all$BA.grwth <- (trees.all$baNext - trees.all$ba) / trees.all$interval

trees.all <- rename(trees.all, year=MEASYEAR, plot=PLOT, subp=SUBP, tree=TREE, status=STATUSCD, dia=DIA, actualht=ACTUALHT, ht=HT, azimuth=AZIMUTH, dist=DIST, lat=LAT, lon=LON, elev=ELEV, genus=GENUS, species=SPECIES, commonName=COMMON_NAME)
trees.all$unique.subp <- paste(trees.all$fia.ID, trees.all$subp, sep='.')
trees.all$plotyear <- paste(trees.all$fia.ID, trees.all$year, sep='.')
trees.all$subplotyear <- paste(trees.all$unique.subp, trees.all$year, sep='.')


## for forest density meas (#ind's/area)
trees.all <- group_by(trees.all, fia.ID, year)%>%
  mutate(totba = sum(ba[status==1&dia>=5], na.rm=T),
         adcount = length(tree.uniqueID[status==1&dia>=5]),
         avgba = mean(ba[status==1&dia>=5], na.rm=T))

# densfactor <- group_by(trees.all,fia.ID, gensp, year)
# plot.densities <- summarize(densfactor,
#                             adba.s = sum(ba[status==1&dia>=5], na.rm=TRUE),
#                             smba.s = sum(ba[status==1&dia<5], na.rm=TRUE),
#                             adcount.s = length(tree.uniqueID[status==1&dia>=5]),
#                             smcount.s = length(tree.uniqueID[status==1 & dia<5]))
# ## Scale up small ind's and combine for density in factors (Accounting for different plot sizes for seedlings)
# # scaleupf <- (24^2*pi)/(6.8^2*pi)
# # plot.densities$totba.consp <- plot.densities$adba+plot.densities$smba*scaleupf
# # plot.densities$totcount.consp <- plot.densities$adcount+round(plot.densities$smcount*scaleupf)
# 
# plot.densities2 <- group_by(plot.densities, fia.ID, year) %>%
#   mutate(totba = sum(adba.s, na.rm=T),
#          adcount = sum(adcount.s, na.rm=T))
# 
# 
# trees.all <- merge(trees.all, plot.densities2, by=c('fia.ID', 'gensp','year'))

trees.all$plotyear <- paste(trees.all$fia.ID, trees.all$year, sep='.')


## Deleted species table and pairwise here


trees.all$alive <- NA; trees.all$surv <- NA
trees.all$alive <- ifelse(trees.all$status==1 , 1, 0)
trees.all$surv <- ifelse(trees.all$status==1 & trees.all$statusNext==1, 1, 
                         ifelse(trees.all$status==1 & trees.all$statusNext==0, 0, NA))  # ! Is this the right thing to do?



## Survival data ----
# ISSUE ----
#  Can we use non-repeat data for survival, based on FIA status code: 
#  
## For repmeas = 2, we want to remove those trees that go status 0 -> 1
## For repmeas = 2, we want to remove second visit for trees that go status 0 -> 0
## For repmeas = 3, we want to remove those trees that go status 1 -> 0 -> 1   OR 0 -> 0 -> 1 
## HOW TO DO THIS? For repmeas = 3, we want to remove the third year for trees that go status  1 -> 0 -> 0   

table(trees.all$status)

# See whether any trees have repeated 'dead' measures
table(trees.all$status, trees.all$statusNext,trees.all$repmeas)
# yes there are some trees with repeat measures recored at dead each time.  Take out the second and third ones:

temp <- group_by(trees.all, tree.uniqueID)
temp1 <- summarize(temp,
                   repmeas = mean(repmeas)
                   ,status.avg = signif(mean(status, na.rm=TRUE), 3)
                   ,status.len = length(status)
                   ,status.v1 = status[visit==1]
                   ,status.v2 = ifelse(min(repmeas)==2, status[visit==2], 999)
                   ,status.v3 = ifelse(min(repmeas)==3, status[visit==3], 999)
)
head(temp1)
summary(temp1)
remove.trees <- c(temp1$tree.uniqueID[temp1$repmeas==2 & temp1$status.v1==0 & temp1$status.v2==1]
                  ,temp1$tree.uniqueID[temp1$repmeas==3 & temp1$status.v1==0 & temp1$status.v2==1]
                  ,temp1$tree.uniqueID[temp1$repmeas==3 & temp1$status.v2==0 & temp1$status.v3==1]
                  ,temp1$tree.uniqueID[temp1$repmeas==3 & temp1$status.v1==1 & temp1$status.v2==0 & temp1$status.v3==1]) 

remove.v2 <- temp1$tree.uniqueID[temp1$repmeas==2 & temp1$status.v1==0 & temp1$status.v2==0 ] 
remove.v3 <- temp1$tree.uniqueID[temp1$repmeas==3 & temp1$status.v1==1 & temp1$status.v2==0 & temp1$status.v3==0] 

trees.all <- trees.all %>%
  filter(!tree.uniqueID %in% remove.trees) %>% 
  filter(!(tree.uniqueID %in% remove.v2 & visit==2)) %>%
  filter(!(tree.uniqueID %in% remove.v3 & visit==3))

table(trees.all$surv)
table(trees.all$repmeas)
table(trees.all$alive)

test <-trees.all[is.na(trees.all$surv)==TRUE & trees.all$repmeas==2,]



## Write All tree data - used for survival and growth ----
write.csv(trees.all, file='Processed Data/All.Treedata.csv', row.names = F)



################# Seedlings----
sdlgs.sub2 <- merge(sdlgs.sub, sp.ref, by='SPCD')
sdlgs.sub2$gensp <- paste(sdlgs.sub2$GENUS,sdlgs.sub2$SPECIES,sep=' ')
sdlgs.sub2 <- rename(sdlgs.sub2, year=MEASYEAR)
head(sdlgs.sub2)



# TREECOUNT_CALC
# Tree count used in calculations. This attribute is set either to COUNTCD, which was dropped in FIADB version 2.1, or TREECOUNT. When PLOT.MANUAL <2.0, the national core procedure was to record the actual seedling count up to six seedlings and then record 6+ if at least six seedlings were present. However, the following regions collected the actual seedling count when PLOT.MANUAL <2.0: Rocky Mountain Research Station (RMRS) and North Central Research Station (NCRS). If PLOT.MANUAL <2.0 and TREECOUNT is blank (null), then a value of 6 in TREECOUNT_CALC represents 6 or more seedlings.

## For manual < 2
# COUNTCD Seedling count code. Indicates the number of seedlings (DIA < 1.0 inches)
# present on the microplot. Conifer seedlings are at least 6 inches tall and
# hardwood seedlings are at least 12 inches tall. New in 1999.
# Code Description
# 1 1 seedling
# 2 2 seedlings
# 3 3 seedlings
# 4 4 seedlings
# 5 5 seedlings
# 6+ 6 or more seedlings
# 6+ indicated when TREECOUNT==0&TREECOUNT_CALC==6

temp<- mutate(sdlgs.sub2, 
              seedlings = ifelse(MANUAL >=2, TREECOUNT_CALC, ifelse(MANUAL<2&is.na(TREECOUNT)&TREECOUNT_CALC==6, 6, TREECOUNT_CALC))

)
head(temp)
sdlgs.sub1 <- temp
head(sdlgs.sub1)

## Merge species-specific data with seedling database ----     < LOSE LAT here >
adultdensity <- group_by(trees.all,fia.ID, gensp, year)%>%
  summarize(adcount.s = length(tree.uniqueID[status==1&dia>=5]))
adultdensity <- subset(adultdensity, adcount.s>0)

temp <- group_by(sdlgs.sub1, gensp, fia.ID, year) ## I think this is aggregating subplots into plots
temp1 <- summarise(temp,
                   N.sdlgs = sum(seedlings, na.rm=TRUE))
sdlgs.plot0 <- as.data.frame(temp1)
dim(sdlgs.plot0)
head(sdlgs.plot0)
sdlg.plots <- unique(sdlgs.plot0$fia.ID)

sdlgs.plot1 <- merge(adultdensity, sdlgs.plot0, by=c("fia.ID", "gensp", "year"), all=T )
sdlgs.plot1$plotyear <- paste(sdlgs.plot1$fia.ID, sdlgs.plot1$year, sep='.')


# should we only include adults for plots where there are seedling records?
trees.plotNEW <- group_by(trees.all, plotyear) %>%
  summarise(fia.ID = unique(fia.ID),
            ptcoeff = mean(ptcoeff, na.rm=T),
            temp = mean(temp, na.rm=T),
            elev = mean(elev, na.rm=T),
            STDORGCD = first(STDORGCD),
            PHYSCLCD = first(PHYSCLCD),
            site_water = first(site_water),
            totba = mean(totba, na.rm=T),
            adcount = mean(adcount,na.rm=T),
            avgba = mean(avgba, na.rm=T)) %>%
  as.data.frame()

trees.plotNEW$avgba[is.nan(trees.plotNEW$avgba)] <-0

head(trees.plotNEW); dim(trees.plotNEW)

sdlg.full<- merge(sdlgs.plot1, trees.plotNEW, by=c("plotyear",'fia.ID'), all.x=TRUE)
head(sdlg.full); dim(sdlg.full)

sdlg.full$N.sdlgs <- ifelse(is.na(sdlg.full$N.sdlgs), 0, sdlg.full$N.sdlgs)
length(which(sdlg.full$N.sdlgs==0))

r1 <- sdlg.full

r1$pres <- ifelse(r1$N.sdlgs==0, 0, 1)

length(which(r1$pres==1&r1$adcount==0))
# test <- subset(r1, adcount==0&pres==1)
# test2 <- merge(test, abundance, by=c("fia.ID", "gensp"))
# r1 <- subset(r1, adcount>0|pres==1)

presence <- subset(occurrence, pres==1)
presence <- rename(presence, ad.pres=pres)
r1 <- merge(r1, presence[1:3], by=c("fia.ID","gensp"), all.x=T)
length(which(r1$pres==1&r1$adcount==0))

# r2$adba.s <- ifelse(is.na(r2$adba.s),0,r2$adba.s)
# r2$smba.s <- ifelse(is.na(r2$smba.s),0,r2$smba.s)
# r2$adcount.s <- ifelse(is.na(r2$adcount.s),0,r2$adcount.s)
# r2$smcount.s <- ifelse(is.na(r2$smcount.s),0,r2$smcount.s)

length(which(is.na(r1$ptcoeff)&is.na(r1$temp)&is.na(r1$elev)&is.na(r1$totba)&r1$pres==1))

#there are 421 records without 
# get plot data for seedlings that don't have adults in the plot (had nothing to merge)
r2 <- merge(r1, plots.sub[c(1,2,4,5,12,13)], by="fia.ID")  # JD: i added LAT and LON to the plot variables to retain here
r2$ptcoeff <- ifelse(!is.na(r2$ptcoeff.x)&is.na(r2$ptcoeff.y), r2$ptcoeff.x,
                 ifelse(is.na(r2$ptcoeff.x)&!is.na(r2$ptcoeff.y), r2$ptcoeff.y, r2$ptcoeff.x))
r2$temp <- ifelse(!is.na(r2$temp.x)&is.na(r2$temp.y), r2$temp.x,
                     ifelse(is.na(r2$temp.x)&!is.na(r2$temp.y), r2$temp.y, r2$temp.x))
r2$elev <- ifelse(!is.na(r2$elev)&is.na(r2$ELEV), r2$elev,
                  ifelse(is.na(r2$elev)&!is.na(r2$ELEV), r2$ELEV, r2$elev))
## r2 <- r2[c(1:79,81:84,87)]  # should be able to delete this, since removing pairwise analysis 

length(which(is.na(r2$ptcoeff)&is.na(r2$temp)&is.na(r2$elev)&is.na(r2$totba)&r2$pres==1))

write.csv(r2, file='Processed Data/Seedling.data.csv', row.names=FALSE)




