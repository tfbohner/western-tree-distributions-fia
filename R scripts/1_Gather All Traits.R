## 1_Gather All Traits
## Compile functional traits (BIEN, TRY, USDA PLANTS, wood density, shade tolerance) and build the species phylogeny.
## Author: Teresa Bohner
## Project: where western trees grow vs. where they thrive (see README)

library(data.table)
library(reshape2)
library(BIEN)
library(rgeos)
library(taxize)
library(ape)
library(ggplot2)
library(tidyverse)

## Data from FIA plots----
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

## Latitudinal extent from BIEN----
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


occurrence <- BIEN_occurrence_species(species = spplist, native.status = T)
occurrence_sub <- filter(occurrence, is_introduced!=1)

lat <- group_by(occurrence_sub, scrubbed_species_binomial)%>%
  summarize(max_range_lat = max(latitude, na.rm=T),
            min_range_lat = min(latitude, na.rm=T),
            lat_extent = max_range_lat-min_range_lat,
            med_range_lat = median(latitude, na.rm=T))

lat <- subset(lat, scrubbed_species_binomial!="Sequoia sempervirens"|scrubbed_species_binomial!="Pinus albicaulis")

## Physical traits from BIEN----
traits <-BIEN_trait_species(species=spplist)
str(traits)
summary <- group_by(traits, scrubbed_species_binomial, trait_name, unit, project_pi)%>%
  summarize(mean_trait_value = ifelse(is.numeric(trait_value),mean(trait_value, na.rm=T),
                                      first(trait_value)))
tr1 <- group_by(summary, trait_name) %>%
  summarise(num.sp = length(unique(scrubbed_species_binomial))) %>%
  as.data.frame()
trait.table <- tr1[order(tr1$num.sp,decreasing=TRUE),]#[1:50,]

tr2 <- group_by(summary, scrubbed_species_binomial) %>%
  summarise(num.traits = length(unique(trait_name))) %>%
  as.data.frame()
trait.table2 <- tr2[order(tr2$num.traits,decreasing=TRUE),]#[1:50,]

interestingtraits <- c("seed mass", 
                       "leaf area per leaf dry mass",
                       "leaf nitrogen content per leaf dry mass",
                       "stem wood density")

bientraits <- subset(summary, trait_name%in%interestingtraits)
bientraits <- arrange(bientraits, scrubbed_species_binomial, trait_name, project_pi)
bientraits <- group_by(bientraits, scrubbed_species_binomial, trait_name, unit)%>%
  summarize(bien_value = first(mean_trait_value),
            bien_project_pi = first(project_pi))%>%
  data.table()

bientraits <- rename(bientraits, bien_unit=unit)
bientraits$trait_name <- str_replace_all(bientraits$trait_name, " ", "_")
bientraits_wide <- dcast.data.table(bientraits, scrubbed_species_binomial~trait_name, 
                                    value.var=c("bien_value", "bien_unit", "bien_project_pi"))
bientraits_all <- merge(lat, bientraits_wide, by="scrubbed_species_binomial", all=T)

## write bien traits data----
write.csv(bientraits_all, "Processed Data/Processed Traits/bien traits.csv", row.names=F)
bientraits_all <- read.csv("Processed Data/Processed Traits/bien traits.csv")


all_traits <- merge(species_summary, bientraits_all, by.x="gensp", by.y="scrubbed_species_binomial", all.x=T)

citingtraits <- BIEN_trait_traitbyspecies(species=spplist, trait = c("leaf area per leaf dry mass", "stem wood density"))

BIEN_metadata_citation(trait.dataframe=citingtraits)
temp_dir <- file.path(tempdir(), "Processed Data/Processed Traits/")
citation_info<-BIEN_metadata_citation(trait.dataframe = citingtraits,
                                      bibtex_file = file.path("Processed Data/Processed Traits/fia_bien_traits.bib"),
                                      acknowledgement_file=file.path("Processed Data/Processed Traits/fia_bie_traits_ack.txt"))

### USDA plants----

### USDA plants wetland indicator species
usda <- read.table("Raw Data/Traits/USDAwetlandReportDownload.txt", sep=",", header=T)

scoretable <- as.data.frame(cbind(c("UPL","FACU","FAC","FACW"), c(1:4)))
names(scoretable) <- c("Indicator", "Score")
usda2 <- usda %>% subset(Scientific.Name %in% spplist) %>%
  subset(Subregion=="Arid West"|Subregion=="Western Mountains, Valleys, and Coast")
usda2$Indicator <- usda2$Indicator %>% as.character()%>%
  str_trim( side = c("right")) %>%
  as.factor()
usda2 <- merge(usda2, scoretable, by="Indicator")


usda3 <- usda2 %>% 
  dplyr::select(c("Scientific.Name", "Indicator", "Subregion")) %>%
  spread(Subregion, Indicator) %>%
  merge(dplyr::select(usda2, c("Scientific.Name", "Score", "Subregion")) %>%
          spread(Subregion, Score), by="Scientific.Name") 
usda3$totscore <- as.numeric(usda3$`Western Mountains, Valleys, and Coast.y`) + as.numeric(usda3$`Arid West.y`)

all_traits <- merge(all_traits, 
                    dplyr::select(usda3, c("Scientific.Name", "totscore")), 
                    by.x="gensp", by.y="Scientific.Name", all=T)

usda4 <- read.csv("Raw Data/Traits/USDAplantsdownloadData.txt", sep=",", header=T)
usda4 <- usda4 %>% subset(Scientific.Name %in% spplist)
usda4$Height.at.Base.Age..Maximum..feet. <- NULL
usda4$Vegetative.Spread.Rate <- NULL

all_traits <- merge(all_traits, 
                    dplyr::select(usda4, -c("Accepted.Symbol", "Synonym.Symbol")),
                    by.x="gensp", by.y="Scientific.Name", all=T)

usda5 <- read.csv("Raw Data/Traits/USDAclasses.txt", sep=",", header=T)
usda5 <- usda5 %>% subset(Scientific.Name %in% spplist)
all_traits <- merge(all_traits, 
                    dplyr::select(usda5, -c("Accepted.Symbol", "Synonym.Symbol")), 
                    by.x="gensp", by.y="Scientific.Name", all=T)


## Wood density from worldagroforestry----
wood_dens2 <- read.delim("Raw Data/Traits/worldagro_wood density.txt")
wood_dens2$X <- NULL
wood_dens2 <- rename(wood_dens2, gensp=species)

all_traits <- merge(all_traits, wood_dens2[1:2], by="gensp", all.x=T)

## combine wood density into one column

all_traits$wood_density <- ifelse(is.na(all_traits$bien_value_stem_wood_density), all_traits$wd, all_traits$bien_value_stem_wood_density)
all_traits$wood_density <- as.numeric(all_traits$wood_density)
all_traits$wood_density_source <- ifelse(!is.na(all_traits$wood_density),
                                         ifelse(is.na(all_traits$bien_value_stem_wood_density), 
                                                "world agroforestry", as.character(all_traits$bien_project_pi_stem_wood_density)), NA)

## Traits from TRY----
# Chamaecyparis nootkatensis -> Callitropsis nootkatensis
# Lithocarpus densiflorus -> Notholithocarpus densiflorus 
# Abies shastensis is a variety of Abies magnifica

## This section for getting the format of species and traits to paste into TRY ----
sp.try <- read.table("Raw Data/Traits/TryAccSpecies.txt",sep='\t',header=TRUE)
head(sp.try)

temp <- merge(data.frame(spplist), sp.try, by.x='spplist', by.y='AccSpeciesName', all.x=TRUE)
head(temp); dim(temp)

filter(temp, is.na(AccSpeciesID)==TRUE)


## Now have TRY data - This section for reading it in etc ----

try.data.all <- read.delim("Raw Data/Traits/TRY/TRY data/4612.txt", header = T)

head(try.data.all); dim(try.data.all)

try.data <- dplyr::select(try.data.all, -c(Reference,Comment))
head(try.data)

length(unique(try.data$SpeciesName)) # 281 varieties of the  61 species
sort(unique(try.data$SpeciesName)) # varieties
length(unique(try.data$AccSpeciesName)) # Use AccSpeciesName !

## Which traits do we have for all species ----

tr1 <- group_by(try.data, TraitName) %>%
  summarise(num.sp = length(unique(AccSpeciesName))) %>%
  as.data.frame()
trait.table <- tr1[order(tr1$num.sp,decreasing=TRUE),]#[1:50,]
trait.table <- droplevels(filter(trait.table, num.sp>25))
dim(trait.table)
trait.table # Look through trait table for interesting ones

# Pick out potentially interesting ones
vars<-list(
  p1 <- "Seed dry mass"
  ,p2 <- "Plant tolerance to drought"
  ,p3 <- "Plant tolerance to shade"
  ,p4 <-"Dispersal unit type"
  ,p5 <-"Plant tolerance to fire" # 58 species
  ,p6 <-"Plant growth rate" # 57 species
  ,p7 <-"Plant tolerance to planting density" # 
  ,p8 <-"Root rooting depth" # 
  ,p9 <-"Seed spread rate" # 
  ,p10 <-"Shoot carbon/nitrogen (C/N) ratio" # 
  ,p11 <-"Dispersal syndrome" # 
  ,p12 <-"Plant mycorrhizal type" 
  ,p13 <- "Plant lifespan (longevity)"

)

try.sub <- subset(try.data, TraitName%in%vars)
try.sub.distinct <- distinct(try.sub, SpeciesName, AccSpeciesName, TraitName, DataName, OriglName, OrigValueStr, OrigUnitStr)
different.meas <- group_by(try.sub.distinct, AccSpeciesName,TraitName, DataName, OriglName, OrigValueStr, OrigUnitStr) %>%
  summarize() %>%
  group_by(TraitName, DataName, OriglName)%>%
  mutate(num.species.per.OriglName = length(unique(AccSpeciesName))) %>%
  group_by(TraitName, DataName, OriglName, OrigValueStr, OrigUnitStr)%>%
  summarize(num.species.per.OriglName=mean(num.species.per.OriglName),
            num.species.per.OrigValueStr = length(unique(AccSpeciesName)))

write.csv(different.meas,"Processed Data/Processed Traits/TraitSummaries.csv",row.names=FALSE)

vars2 <- list(p1 <- "Leaf area per leaf dry mass (SLA or 1/LMA): undefined if petiole and rachis are in- or excluded"
  ,p2 <- "Leaf area per leaf dry mass (SLA or 1/LMA): petiole and rachis excluded"
  ,p3 <- "Leaf area per leaf dry mass (specific leaf area, SLA)"
)
sla <- subset(try.data, TraitName%in%vars2)%>%
  distinct(SpeciesName, AccSpeciesName, TraitName, DataName, OriglName, OrigValueStr, OrigUnitStr) #%>%
sla$value <- as.numeric(as.character(sla$OrigValueStr))

try_sla_spp <- unique(sla$AccSpeciesName)%>%
  as.character()

all_sla_spp <- unique(all_traits$gensp[which(is.finite(all_traits$bien_value_leaf_area_per_leaf_dry_mass))])%>%
  as.character()

try_only <- setdiff(try_sla_spp, all_sla_spp)

sla2 <- subset(sla, AccSpeciesName%in%try_only) 

sla2$OrigUnitStr <- as.character(sla2$OrigUnitStr)
sla2$OrigUnitStr <- ifelse(sla2$OrigUnitStr=="cm2/g", "cm2 g-1", sla2$OrigUnitStr)
sla2$newvalue_sla <- ifelse(sla2$OrigUnitStr=="cm2 g-1", sla2$value*0.1, sla2$value)
sla2$newunit_sla <- "m2.kg-1"

sla3 <- group_by(sla2, AccSpeciesName, TraitName, DataName) %>%
  summarize(try_value_sla = mean(newvalue_sla, na.rm=T),
            try_unit_sla = first(newunit_sla))

all_traits <- merge(all_traits, sla3[c(1,4,5)], by.x="gensp", by.y="AccSpeciesName", all=T)
all_traits$sla_value <- ifelse(!is.na(all_traits$bien_value_leaf_area_per_leaf_dry_mass), 
                               all_traits$bien_value_leaf_area_per_leaf_dry_mass, all_traits$try_value_sla)
all_traits$sla_unit <- ifelse(!is.na(all_traits$bien_unit_leaf_area_per_leaf_dry_mass), 
                               all_traits$bien_unit_leaf_area_per_leaf_dry_mass, all_traits$try_unit_sla)


try.height <- subset(try.data, TraitName=="Plant height vegetative")%>%
  distinct(SpeciesName, AccSpeciesName, TraitName, DataName, OriglName, OrigValueStr, OrigUnitStr)
try.height$value <- as.numeric(as.character(try.height$OrigValueStr))

try.height.spp <- unique(try.height$AccSpeciesName)%>%
  as.character()

all.height.spp <- unique(all_traits$gensp[which(is.finite(all_traits$Height..Mature..feet.))])%>%
  as.character()

# More height measurements from TRY alone

try.height2 <- subset(try.height, OriglName=="MaximumHeightM")
height.2.spp <- unique(try.height2$AccSpeciesName)
try.height3 <- subset(try.height, OriglName=="MaximumHeightMaxM")
height.3.spp <- unique(try.height3$AccSpeciesName)

all.spp <- unique(species_summary$gensp)
setdiff(all.spp, height.2.spp)
setdiff(all.spp, height.3.spp)

intersect(setdiff(all.spp, height.2.spp), setdiff(all.spp, height.3.spp))
# only missing Abies shastensis

try.height4 <- rbind(try.height2, try.height3) %>%
  group_by(AccSpeciesName) %>%
  summarize(MaxHeightM = mean(value, na.rm=T))

all_traits <- merge(all_traits, try.height4, by.x="gensp", by.y="AccSpeciesName", all=T)
# check try and USDA value comparison
plot(all_traits$MaxHeightM, all_traits$Height..Mature..feet.)
abline(0,1)

try.wood3 <- group_by(try.wood2, AccSpeciesName, TraitName) %>%
  summarize(try_value_wd = mean(newvalue_wd, na.rm=T),
            try_unit_wd = first(newunit_wd))

try.wood <- subset(try.data, TraitName=="Stem dry mass per stem fresh volume (stem specific density, SSD, wood density)")%>%
  distinct(SpeciesName, AccSpeciesName, TraitName, DataName, OriglName, OrigValueStr, OrigUnitStr)
try.wood$value <- as.numeric(as.character(try.wood$OrigValueStr))

try_wood_spp <- unique(try.wood$AccSpeciesName)%>%
  as.character()

all_wood_spp <- unique(all_traits$gensp[which(is.finite(all_traits$wood_density))])%>%
  as.character()

try_only2 <- setdiff(try_wood_spp, all_wood_spp)

try.wood2 <- subset(try.wood, AccSpeciesName%in%try_only2)

try.wood2$OrigUnitStr <- as.character(try.wood2$OrigUnitStr)
try.wood2$newunit_wd <- "g.cm-3"
try.wood2$newvalue_wd <- try.wood2$value

try.wood3 <- group_by(try.wood2, AccSpeciesName, TraitName) %>%
  summarize(try_value_wd = mean(newvalue_wd, na.rm=T),
            try_unit_wd = first(newunit_wd))

all_traits <- merge(all_traits, try.wood3[c(1,3,4)], by.x="gensp", by.y="AccSpeciesName", all=T)
all_traits$wood_density <- ifelse(!is.na(all_traits$wood_density), 
                               all_traits$wood_density, all_traits$try_value_wd)
all_traits$wd_unit <- ifelse(!is.na(all_traits$wood_density), "g.cm-3", NA)
all_traits$wood_density_source <- ifelse(!is.na(all_traits$try_value_wd), "TRY", 
                                         all_traits$wood_density_source)


try.root <- subset(try.data, TraitName=="Root rooting depth")%>%
  distinct(SpeciesName, AccSpeciesName, TraitName, DataName, OriglName, OrigValueStr, OrigUnitStr)
try.root$value <- as.numeric(as.character(try.root$OrigValueStr))
try.root <- try.root[which(!is.na(try.root$value)),]
try.root.spp <- unique(try.root$AccSpeciesName)
try.root.sub <- subset(try.root, OriglName=="Root Depth, Minimum")
try.root.sub.spp <- unique(try.root.sub$AccSpeciesName)
usda.root.spp <- unique(all_traits$gensp[which(!is.na(all_traits$Root.Depth..Minimum..inches.))])

setdiff(usda.root.spp, try.root.sub.spp)

quant.names<-list(
  p1 <- "PlantLifeSpan"
  ,p2 <- "1000SeedWeight"
  ,p3 <- "Root Depth, Minimum")
  
clean.quant <- subset(try.sub.distinct, OriglName%in%quant.names) 
clean.quant$value <- as.numeric(as.character(clean.quant$OrigValueStr))
clean.quant$trait <- ifelse(clean.quant$TraitName=="Seed dry mass", "SeedMass",
                            ifelse(clean.quant$TraitName=="Root rooting depth", "RootDepth", "Lifespan"))
clean.quant <- clean.quant[which(!is.na(clean.quant$value)),] %>%
  group_by(AccSpeciesName, trait) %>%
  summarize(value=mean(value),
            unit=first(OrigUnitStr)) %>%
  data.table() %>%
  data.table::dcast(AccSpeciesName~trait, value.var = c("value","unit"))
  
write.csv(clean.quant,"Processed Data/Processed Traits/TRY_CleanQuantTraits.csv",row.names=FALSE)
all_traits <- merge(all_traits, clean.quant, by.x="gensp", by.y = "AccSpeciesName", all=T)

trait.names<-list(
  p1 <- "Drought Tolerance"
  ,p2 <- "TraitDroughtTolerance"
  ,p3 <- "Fire Resistant"
  ,p4 <-"Shade Tolerance"
  ,p5 <-"TraitShadeTolerance" 
  ,p6 <-"Growth Rate" 
  ,p7 <-"Seed Spread Rate" # 
  ,p8 <-"C:N Ratio" # 
)

clean.qual <- subset(try.sub.distinct, OriglName%in%trait.names)
clean.qual <- clean.qual[which(!is.na(clean.qual$OrigValueStr)),]
clean.qual2 <- rename(clean.qual, gensp=AccSpeciesName, trait=TraitName, value=OrigValueStr) %>%
  group_by(gensp, trait, value) %>%
  summarize()

unique(clean.qual2$value[which(clean.qual2$trait=="Plant tolerance to drought")])
drought.cat <- c("Low", "Medium", "High")
unique(clean.qual2$value[which(clean.qual2$trait=="Plant tolerance to shade")])
shade.cat <- c("Tolerant", "Intermediate", "Intolerant")

clean.qual3 <- clean.qual2 %>%
  filter(!(trait=="Plant tolerance to drought"&!value%in%drought.cat),
           !(trait=="Plant tolerance to shade"&!value%in%shade.cat)) %>%
  group_by(gensp, trait) %>%
  mutate(newvalue = paste(unique(value), collapse="/")) %>%
  group_by(gensp, trait, newvalue) %>%
  summarize() %>%
  data.table %>%
  data.table::dcast(gensp~trait, value.var="newvalue")

try.sub.distinct$TraitName <- as.character(try.sub.distinct$TraitName)
dispersal <- subset(try.sub.distinct, TraitName=="Dispersal syndrome")
dispersal$OrigValueStr <- tolower(dispersal$OrigValueStr)

dispersal2 <- filter(dispersal, !OriglName%in%c("dispersal type", "diaspore type", "Seed Dispersal", "Dispersal syndrome")) %>%
  as.tibble() 
dispersaltypes <- unique(dispersal2$OrigValueStr)

dispersal2$DispMode <- NA
dispersal2$DispMode[which(str_detect(dispersal2$OrigValueStr, "(animal|intentionally|bird)")==TRUE)] <- "animal"
dispersal2$DispMode[which(str_detect(dispersal2$OrigValueStr, "(wind)")==TRUE)] <- "wind"
dispersal2$DispMode[which(str_detect(dispersal2$OrigValueStr, "(unassisted|drop)")==TRUE)] <- "unassisted"

dispersal3 <- rename(dispersal2, gensp=AccSpeciesName, trait=TraitName) %>%
  group_by(gensp, DispMode) %>%
  summarize()
clean.traits.all <- merge(clean.qual3, dispersal3, by="gensp", all=T)

write.csv(clean.traits.all,"Processed Data/Processed Traits/TRY_CleanQualTraits.csv",row.names=FALSE)
all_traits <- merge(all_traits, clean.traits.all, by="gensp", all=T)

# ## Ordinal tolerance data from Niinemets and Valladares 2006
shade <- read.csv("Raw Data/Traits/Niinemets and Valladares 2006 shade tolerance.csv")
shade <- subset(shade, Data.set.1=="North America ")
shade <- rename(shade, gensp=Species)
shade$gensp <- str_sub(shade$gensp, 0, -2)
shade$gensp <- sub("Chrysoplepis","Chrysolepis",shade$gensp)
shade$gensp <- sub("Pinus contorta ssp. contorta","Pinus contorta",shade$gensp)
shade <- subset(shade, gensp%in%spplist)

setdiff(spplist, shade$gensp)

all_traits <- merge(all_traits, 
                    dplyr::select(shade, c("gensp", "Evergreen", "Gymnosperm")), by="gensp", all.x=T)

all_traits$Category <- as.character(all_traits$Category)
all_traits$Leaf.Retention <- as.character(all_traits$Leaf.Retention)
all_traits$Gymnosperm <- as.character(all_traits$Gymnosperm)
all_traits$Evergreen <- as.character(all_traits$Evergreen)

all_traits$Gymnosperm[which(all_traits$Gymnosperm=="Y ")] <- "Gymnosperm"
all_traits$Gymnosperm[which(all_traits$Gymnosperm=="N ")] <- "Dicot"
all_traits$Evergreen[which(all_traits$Evergreen=="Y ")] <- "Evergreen"
all_traits$Evergreen[which(all_traits$Evergreen=="N ")] <- "Deciduous"
all_traits$Leaf.Retention[which(all_traits$Leaf.Retention=="Yes")] <- "Evergreen"
all_traits$Leaf.Retention[which(all_traits$Leaf.Retention=="No")] <- "Deciduous"
all_traits$Leaf.Retention[which(all_traits$Leaf.Retention=="")] <- NA

all_traits$gymno2 <- ifelse(is.na(all_traits$Category), all_traits$Gymnosperm, all_traits$Category)
all_traits$evergreen2 <- ifelse(is.na(all_traits$Leaf.Retention), all_traits$Evergreen, all_traits$Leaf.Retention)


funfirst <- function(x) {
  if(all(is.na(x))){NA} else{
    first(na.omit(x))} 
  
}

weirdos <- subset(all_traits, fia_gensp=="Abies shastensis"|fia_gensp=="Lithocarpus densiflorus"|fia_gensp=="Chamaecyparis nootkatensis")
weirdos2 <- dplyr::select(weirdos, -gensp) %>%
  group_by(fia_gensp) %>%
  summarise_all(first)
weirdos2$gensp <- weirdos2$fia_gensp

all_traits <- subset(all_traits, fia_gensp!="Abies shastensis"&fia_gensp!="Lithocarpus densiflorus"&fia_gensp!="Chamaecyparis nootkatensis")%>%
  rbind(weirdos2)

##phylogenetic tree----
out <- classification(as.character(unique(all_traits$gensp)), db="ncbi")
tr <- class2tree(out)
plot(tr)
test <- tr$classification

families <- dplyr::select(tr$classification, c(species, family))
rownames(families) <- NULL
families <- families %>% distinct()
all_traits <- merge(all_traits, families, by.x="gensp", by.y="species", all.x=T)
names(all_traits)


## Cleaning any weirdness----
str(all_traits)
all_traits$bien_value_leaf_area_per_leaf_dry_mass <- as.numeric(all_traits$bien_value_leaf_area_per_leaf_dry_mass)
all_traits$logSeedMass <- log(all_traits$value_SeedMass)
all_traits$hydrophyte <- ifelse(!is.na(all_traits$totscore)&all_traits$totscore>=6, "yes", "no")

write.csv(all_traits, "Processed Data/Analysis Ready Data/all_traits.csv", row.names=F)


