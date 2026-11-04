## 3_Range Maps
## Check how much of each species' range (Little's range maps) the inventory plots cover.
## Author: Teresa Bohner
## Project: where western trees grow vs. where they thrive (see README)

library(data.table)
library(reshape2)
library(BIEN)
library(rgeos)
library(ggplot2)
library(tidyverse)

data.path <- "Processed Data/Analysis Ready Data/"
abundance <- read.csv(paste0(data.path, "abundance.csv"))

species_summary <- group_by(abundance, gensp)%>%
  summarize(med_plot_elev = median(ELEV, na.rm=T),
            med_plot_lat = median(LAT, na.rm=T),
            med_plot_swb = median(ptcoeff, na.rm=T),
            med_plot_temp = median(temp, na.rm=T),
            max_plot_lat = max(LAT, na.rm=T),
            min_plot_lat = min(LAT, na.rm=T),
            lat_extent = max_plot_lat - min_plot_lat)
species_summary$fia_gensp <- species_summary$gensp

alt_spp <- data.frame(c("Abies magnifica var. shastensis",
                        "Cupressus nootkatensis",
                        "Callitropsis nootkatensis",
                        "Notholithocarpus densiflorus"))
names(alt_spp) <- "gensp"
alt_spp$fia_gensp <- c("Abies shastensis",
                       "Chamaecyparis nootkatensis",
                       "Chamaecyparis nootkatensis",
                       "Lithocarpus densiflorus")

species_summary <- merge(species_summary, alt_spp, by=c("gensp", "fia_gensp"), all=T)
spplist <- as.character(unique(species_summary$gensp))

states <- readOGR("Raw Data/GIS/cb_2016_us_state_500k/cb_2016_us_state_500k.shp")
states <- readOGR("Raw Data/GIS/cb_2016_us_state_500k/cb_2016_us_state_500k.shp")
states1 <-  states[states$NAME=="California"|states$NAME=="Oregon"|states$NAME=="Washington",]
states.t <- spTransform(states1, CRS("+init=epsg:4326"))

states2 <- spTransform(states.t, 
                       CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96"))

spptable <- as.data.frame(cbind(spplist, 1:65)) 
names(spptable) <- c("gensp", "spid")
spptable$gensp <- as.character(spptable$gensp)
spptable$overlap <- NA
spptable$overlap_bien <-NA

## Little maps----
data.path2 <- "Raw Data/GIS/E. Little Maps/USTreeAtlas-master/"
mapspp <- read.csv(paste0(data.path2, "Little_datatable.csv"))
names(mapspp) <- c("gensp", "common", 'ref', 'map', 'filesp')
spptable <- merge(spptable, mapspp, by="gensp", all.x=T)
spptable$spid <- 1:65
spptable$filesp <- as.character(spptable$filesp)

for(i in 1:nrow(spptable)) {
  if(!is.na(spptable$filesp[i])) {
    dsn = paste("Raw Data/GIS/E. Little Maps/USTreeAtlas-master/SHP/", spptable$filesp[which(spptable$spid==i)], "/", sep='')
    shape <- readOGR(dsn = dsn , layer = spptable$filesp[which(spptable$spid==i)])
    proj4string(shape) <- CRS("+init=epsg:4267")
    shape2 <- spTransform(shape, 
                          CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96"))
    
    intersection<- gIntersection(shape2, states2)
    spptable$overlap[i] <- round(sum(gArea(intersection))/gArea(shape2)*100, 2)
  } else {
    spptable$overlap[i] <- NA
  }
  
}

plot(shape2, col="green")

for(i in 22:nrow(spptable)) {
  spp <- BIEN_ranges_load_species(species=spptable$gensp[which(spptable$spid==i)])
  spp2 <- spTransform(spp,
                      CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96"))
  # this is a well known R / GEOS hack (usually combined with the above) to
  # deal with "bad" polygons
  spp2 <- gBuffer(spp2, byid=TRUE, width=0)
  
  intersection<- gIntersection(spp2, states2)
  spptable$overlap_bien[i] <- round(sum(gArea(intersection))/gArea(spp2)*100, 2)
  
}

