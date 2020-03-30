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
admin1_master <- shapefile("K:/AMR/Shapefiles/admin2013_1.shp")

#get the centroids
admin1 <-  admin1_master[admin1_master$GAUL_CODE%in%hosp_adm$GAUL_ADM1,]
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

#get the centroids
admin1 <-  admin1_master[admin1_master$GAUL_CODE%in%ad1$gaul_code,]
ad1centroids <- gCentroid(admin1, byid = T)
ad1centroids <- data.frame(ad1centroids)
ad1centroids$gaul_code <- admin1$GAUL_CODE

#admin 2 shapefile
admin2_master <- shapefile("K:/AMR/Shapefiles/admin2013_2.shp")
admin2 <-  admin2_master[admin2_master$GAUL_CODE%in%ad2$gaul_code,]
ad2centroids <- gCentroid(admin2, byid = T)
ad2centroids <- data.frame(ad2centroids)
ad2centroids$gaul_code <- admin2$GAUL_CODE

centroids <-  rbind(ad1centroids, ad2centroids)
centroids$gaul_code <- as.numeric(centroids$gaul_code)
rm(ad1, ad2, ad1centroids, ad2centroids, admin1, admin2)

#add the centroids onto the data
mydata <- merge(mydata, centroids, by = 'gaul_code', all.x = T)
mydata$long[is.na(mydata$long) & !is.na(mydata$x)] <- mydata$x[is.na(mydata$long) & !is.na(mydata$x)]
mydata$lat[is.na(mydata$lat) & !is.na(mydata$y)] <- mydata$y[is.na(mydata$lat) & !is.na(mydata$x)]
rm(centroids, new.data, hosp_coords)
mydata$x <- NULL
mydata$y <- NULL

#get the admin1 and admin 2 polygons for each locations
coords <- unique(mydata[c('long', 'lat')])
coords <- coords[!is.na(coords$lat),]
coords <- SpatialPoints(coords, CRS(projection(admin1_master)))
coords$adm1 <- over(coords, admin1_master)$GAUL_CODE
coords$adm2 <- over(coords, admin2_master)$GAUL_CODE
coords <- data.frame(coords)
coords$adm1 <- as.numeric(as.character(coords$adm1))
coords$adm2 <- as.numeric(as.character(coords$adm2))
coords$optional <- NULL

check <- coords[is.na(coords$adm2),]
if(length(check$lat>0)){
  pts <- SpatialPoints(check[c('long', 'lat')], CRS(projection(admin1_master)))
  proj4string(admin1_master)
  # utmStr <- "+proj=utm +zone=%d +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
  # crs <- CRS(sprintf(utmStr, 32))
  # adm2UTM <- spTransform(admin2, crs)
  # ptsUTM <- spTransform(check, crs)
  
  
  ## Set up container for results
  n <- length(ptsUTM)
  nearest_adm2 <- character(n)
  nearest_adm1 <- character(n)
  
  ## For each point, find name of nearest polygon (in this case, Belgian cantons)
  for (i in seq_along(nearest_adm2)) {
    nearest_adm2[i] <- admin2_master$GAUL_CODE[which.min(gDistance(pts[i,], admin2_master, byid=TRUE))]
    nearest_adm1[i] <- admin2_master$PARENT_ID[which.min(gDistance(pts[i,], admin2_master, byid=TRUE))]
  }
  nearest_adm2 <- as.numeric(nearest_adm2)
  nearest_adm1 <- as.numeric(nearest_adm1)
  check$adm1 <- nearest_adm1
  check$adm2 <- nearest_adm2
  coords <- coords[!is.na(coords$adm2),]
  coords <- rbind(coords, check)
}

#merge the admin 1 and 2 codes onto the data and clean up

mydata <- merge(mydata, coords, by = c('lat', 'long'), all.x = T)
rm(check, pts, i, n, nearest_adm1, nearest_adm2, coords)

#calculate the rpoprotion of dr
mydata$rate <- mydata$number_resistant/mydata$no_examined

#calculate variance as the weight
mydata$weight <- (mydata$rate*(1-mydata$rate))/mydata$no_examined

mydata$shapefile <- NA
mydata$lat <- NA
mydata$long <- NA

mydata <- mydata[c("row_id",
                   "iso3",
                   "adm1",
                   "adm2",
                   'shapefile',
                   'lat',
                   'long',
                   "mid_year",
                   "species",
                   "antimicrobial",
                   "number_resistant",
                   "no_examined",
                   'rate',
                   'weight')]


#subset data for each drug-bug combo 
mdr_typhi <- mydata[mydata$antimicrobial == 'multi-drug' |
                      mydata$antimicrobial == 'multi-drug-adjusted',]
mdr_typhi <- mdr_typhi[mdr_typhi$species == 'Typhi',]

mdr_typhi$species <- NULL
mdr_typhi$antimicrobial <- NULL
colnames(mdr_typhi)[colnames(mdr_typhi)=='iso3'] <- 'country'
colnames(mdr_typhi)[colnames(mdr_typhi)=='mid_year'] <- 'year'
colnames(mdr_typhi)[colnames(mdr_typhi)=='number_resistant'] <- 'mdr_typhi'
colnames(mdr_typhi)[colnames(mdr_typhi)=='no_examined'] <- 'N'

adm1 <- mdr_typhi
adm1$shapefile <- 'admin2013_1'
colnames(adm1)[colnames(adm1)=='adm1'] <- 'location_code'
adm1$adm2 <-  NULL

adm2 <- mdr_typhi
adm2$shapefile <- 'admin2023_2'
colnames(adm2)[colnames(adm2)=='adm2'] <- 'location_code'
adm2$adm1 <-  NULL

write.csv(adm1, 'K:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/mdr_typhi_mbg_adm1.csv', row.names = F)
write.csv(adm2, 'K:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/mdr_typhi_mbg_adm2.csv', row.names = F)
