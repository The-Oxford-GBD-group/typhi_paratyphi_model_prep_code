#Roughly prep the data for the ST-GPR
#use only the first location identified

library(raster)
library(maptools)
library(rgeos)
library(googlesheets)
rm(list = ls())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Merge on GBD location IDs  #
#    for GBD 2019               # 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

setwd("K:/AMR/Pathogens/typhi_paratyphi/model_prep")

new.data <- readRDS('clean_data/clean_data_crosswalked.rds')

#merge on the hospital location information
locations <- gs_read(ss = gs_title("Prevalence data - Typhi & Paratyphi"), ws = "facility lookup")

#split into those with coords, admin1's and nationals, and get coords
hosp_coords <- locations[!is.na(locations$lat),]
hosp_coords <- hosp_coords[c('facility_name',
                             'lat',
                             'long')]

#get the coords of hospitals with only admin information
hosp_adm <- locations[!is.na(locations$GAUL_ADM1) & is.na(locations$lat),]
hosp_adm <- hosp_adm[c("facility_name", 
                       "GAUL_ADM1")]

#read in admin1 shapefile
admin1 <- shapefile("K:/AMR/Shapefiles/admin2013_1.shp")
#get the centroids
admin1 <-  admin1[admin1$GAUL_CODE%in%hosp_adm$GAUL_ADM1,]
centroids <- gCentroid(admin1, byid = T)
centroids <- data.frame(centroids)
centroids$GAUL_ADM1 <- admin1$GAUL_CODE
hosp_adm <- merge(hosp_adm, centroids, by = 'GAUL_ADM1', all.x = T)
hosp_adm$GAUL_ADM1 <- NULL
colnames(hosp_adm) <- c('facility_name',
                        'long',
                        'lat')

hosp_coords <- rbind(hosp_coords, hosp_adm)

rm(centroids, locations, hosp_adm, admin1)

#merge on hospital lats and longs
mydata <- merge(new.data, hosp_coords, by.x = 'hospital_name1' , by.y = 'facility_name', all.x = T)

#combined coords of community and hosptials to one set of vars
mydata$lat[is.na(mydata$lat) & !is.na(mydata$latitude)] <- mydata$latitude[is.na(mydata$lat) & !is.na(mydata$latitude)]
mydata$long[is.na(mydata$long) & !is.na(mydata$longitude)] <- mydata$longitude[is.na(mydata$long) & !is.na(mydata$longitude)]

#get locations with admin units but no coords
ad1 <- mydata[is.na(mydata$lat) & mydata$admin_level == 1,]
ad2 <- mydata[is.na(mydata$lat) & mydata$admin_level == 2,]

#get the centroids for these locations
admin1 <- shapefile("K:/AMR/Shapefiles/admin2013_1.shp")
#get the centroids
admin1 <-  admin1[admin1$GAUL_CODE%in%ad1$gaul_code,]
ad1centroids <- gCentroid(admin1, byid = T)
ad1centroids <- data.frame(ad1centroids)
ad1centroids$gaul_code <- admin1$GAUL_CODE

#admin 2 shapefile
admin2 <- shapefile("K:/AMR/Shapefiles/admin2013_2.shp")
admin2 <-  admin2[admin2$GAUL_CODE%in%ad2$gaul_code,]
ad2centroids <- gCentroid(admin2, byid = T)
ad2centroids <- data.frame(ad2centroids)
ad2centroids$gaul_code <- admin2$GAUL_CODE

centroids <-  rbind(ad1centroids, ad2centroids)
centroids$gaul_code <- as.numeric(centroids$gaul_code)
rm(ad1, ad2, admin1, admin2, ad1centroids, ad2centroids)

#add the centroids onto the data
mydata <- merge(mydata, centroids, by = 'gaul_code', all.x = T)
mydata$long[is.na(mydata$long) & !is.na(mydata$x)] <- mydata$x[is.na(mydata$long) & !is.na(mydata$x)]
mydata$lat[is.na(mydata$lat) & !is.na(mydata$y)] <- mydata$y[is.na(mydata$lat) & !is.na(mydata$x)]
rm(centroids, new.data, hosp_coords)

#clean up columns
#get the IHME location IDs from the shapefile based on the coords
IHME_shp <- shapefile('K:/AMR/Shapefiles/GBD2019/GBD2019_analysis_final_loc_set_22.shp')

#get the IHME location IDs for each country without subnational info
IHME_shp_subnat <- IHME_shp[IHME_shp$level == 4,]
IHME_shp_nat <- IHME_shp[IHME_shp$level == 3,]

coords <- unique(mydata[c('long', 'lat')])
coords <- coords[!is.na(coords$lat),]
coords <- SpatialPoints(coords, CRS(projection(IHME_shp)))
coords$IHME_loc_4 <- over(coords, IHME_shp_subnat)$loc_id
coords$IHME_loc_3 <- over(coords, IHME_shp_nat)$loc_id
coords <- data.frame(coords)
coords$IHME_loc_4 <- as.numeric(as.character(coords$IHME_loc_4))
coords$IHME_loc_3 <- as.numeric(as.character(coords$IHME_loc_3))
coords$loc_id <- coords$IHME_loc_4
coords$loc_id[is.na(coords$loc_id)] <- coords$IHME_loc_3[is.na(coords$loc_id)]
coords$IHME_loc_3 <- NULL
coords$IHME_loc_4 <- NULL
coords$optional <- NULL

check <- coords[is.na(coords$loc_id),]
coords$loc_id[is.na(coords$loc_id)] <- 354
rm(check)

#get the country (iso code) and location id of the remaining national surveys
IHME_shp_nat <- data.frame(IHME_shp_nat)
countries <- unique(IHME_shp_nat[c('ihme_lc_id', 'loc_id')])
colnames(countries) <- c('iso3', 'loc_id')

#merge the location IDS onto the data
mydata <- merge(mydata, coords, by = c('lat', 'long'), all.x = T)

mydata <- merge(mydata, countries, by = 'iso3', all.x = T)
mydata$loc_id.x[is.na(mydata$loc_id.x)] <- mydata$loc_id.y[is.na(mydata$loc_id.x)]

check <- mydata[is.na(mydata$loc_id.x),]


rm(check, coords, countries, IHME_shp, IHME_shp_nat, IHME_shp_subnat)

#jitter 0's by 0.001 to allow variance calculations
mydata$prop_resistant <- mydata$number_resistant/mydata$no_examined
mydata$prop_resistant[mydata$prop_resistant==0] <- 0.001
mydata$prop_resistant[mydata$prop_resistant==1] <- 0.999

#calculate variance
mydata$variance <- (mydata$prop_resistant*(1-mydata$prop_resistant))/mydata$no_examined

#add other variables required
mydata$age_group_id <-  22
mydata$sex_id <- 3
mydata$measure <- 'continuous'
mydata$is_outlier <- 0

#rename some variables
colnames(mydata)[colnames(mydata)=='no_examined'] <- 'sample_size'
colnames(mydata)[colnames(mydata)=='source_id'] <- 'nid'
colnames(mydata)[colnames(mydata)=='prop_resistant'] <- 'data'
colnames(mydata)[colnames(mydata)=='mid_year'] <- 'year_id'
colnames(mydata)[colnames(mydata)=='loc_id.x'] <- 'location_id'

#subset to required variables
mydata <- mydata[c("nid",
                   "location_id",
                   "country",
                   "year_id",
                   "age_group_id",
                   "sex_id",
                   "measure",
                   "data",
                   "variance",
                   "sample_size",
                   "is_outlier",
                   'antimicrobial',
                   'species')]


#subset data for each drug-bug combo 
mdr_typhi <- mydata[mydata$antimicrobial == 'multi-drug' |
                     mydata$antimicrobial == 'multi-drug-adjusted',]
mdr_typhi <- mdr_typhi[mdr_typhi$species == 'Typhi',]


fqns_typhi <- mydata[grep('FQR', mydata$antimicrobial),]
fqns_typhi <- fqns_typhi[fqns_typhi$species == 'Typhi',]

mdr_paratyphi <- mydata[mydata$antimicrobial == 'multi-drug' |
                      mydata$antimicrobial == 'multi-drug-adjusted',]
mdr_paratyphi <- mdr_paratyphi[mdr_paratyphi$species == 'Paratyphi',]


fqns_paratyphi <- mydata[grep('FQR', mydata$antimicrobial),]
fqns_paratyphi <- fqns_paratyphi[fqns_paratyphi$species == 'Paratyphi',]

rm(mydata)

mdr_paratyphi$antimicrobial <- NULL
mdr_paratyphi$species <- NULL

mdr_typhi$antimicrobial <- NULL
mdr_typhi$species <- NULL

fqns_paratyphi$antimicrobial <- NULL
fqns_paratyphi$species <- NULL

fqns_typhi$antimicrobial <- NULL
fqns_typhi$species <- NULL

write.csv(mdr_typhi, 'K:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/mdr_typhi.csv', row.names = F)
write.csv(fqns_typhi, 'K:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/fqns_typhi.csv', row.names = F)
write.csv(mdr_paratyphi, 'K:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/mdr_paratyphi.csv', row.names = F)
write.csv(fqns_paratyphi, 'K:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/fqns_paratyphi.csv', row.names = F)
