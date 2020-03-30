#Roughly prep the data for the ST-GPR
#use only the first location identified
rm(list = ls())
library(raster)
library(foreign)
library(maptools)
library(rgeos)
library(googlesheets)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Merge on GBD location IDs  #
#    for GBD 2019               # 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

setwd("K:/AMR/Pathogens/typhi_paratyphi/model_prep")

mydata <- readRDS('clean_data/clean_data_crosswalked.rds')

#merge on the hospital location information
locations <- gs_read(ss = gs_title("Prevalence data - Typhi & Paratyphi"), ws = "facility lookup")
locations <- data.frame(locations)

#get the hospital information
hosp_gauls <- locations[c('facility_name',
                             'GAUL_ADM1',
                             'GAUL_ADM2')]

hosp_gauls$GAUL_ADM1 <- as.numeric(hosp_gauls$GAUL_ADM1) 
hosp_gauls$GAUL_ADM2 <- as.numeric(hosp_gauls$GAUL_ADM2) 

hosp_gauls <- hosp_gauls[!is.na(hosp_gauls$GAUL_ADM1),]

#merge this onto the data (just look at hospital 1 for the moment)
mydata <- merge(mydata, hosp_gauls, by.x = 'hospital_name1', by.y = 'facility_name', all.x = T, all.y = F)

#replace the admin 1 and 2 columns for hospitals with those from community if avaialble
#long version as somethings not happy here
admin1s <- mydata[is.na(mydata$GAUL_ADM1),]
admin1s <- mydata[mydata$admin_level == 1,]
admin1s$GAUL_ADM1 <- admin1s$gaul_code

admin2s <- mydata[is.na(mydata$GAUL_ADM2),]
admin2s <- mydata[mydata$admin_level == 2,]
admin2s$GAUL_ADM2 <- admin2s$gaul_code
admin2s$GAUL_ADM1 <- NULL

#lookup admin1s 
lookup <- read.dbf('K:/AMR/Shapefiles/admin2013_2.dbf')
lookup <- lookup[c('GAUL_CODE', 'PARENT_ID')]
names(lookup) <- c('GAUL_ADM2', 'GAUL_ADM1')
admin2s <- merge(admin2s, lookup, by = 'GAUL_ADM2')

mydata <- mydata[!(mydata$row_id%in%admin1s$row_id),]
mydata <- mydata[!(mydata$row_id%in%admin2s$row_id),]

mydata <- rbind(mydata, admin1s, admin2s)
rm(lookup, admin1s, admin2s)

#get admin levels based on centroids of community coords
comm_coords <- mydata[!is.na(mydata$latitude),]
comm_coords <- comm_coords[is.na(comm_coords$GAUL_ADM1),]
pts <- comm_coords[, .(longitude, latitude)] 
admin2 <- shapefile('K:/AMR/Shapefiles/admin2013_2.shp')

pts <- SpatialPoints(pts, CRS(projection(admin2)))
comm_coords$GAUL_ADM2 <- over(pts, admin2)$GAUL_CODE
comm_coords$GAUL_ADM1 <- over(pts, admin2)$PARENT_ID
comm_coords$GAUL_ADM2 <- as.numeric(comm_coords$GAUL_ADM2)
comm_coords$GAUL_ADM1 <- as.numeric(comm_coords$GAUL_ADM1)

#merge with other data
mydata <- mydata[!(mydata$row_id%in%comm_coords$row_id),]
mydata <- rbind(mydata, comm_coords)

#look at data withough gaul codes
check <- mydata[is.na(mydata$GAUL_ADM1),]

rm(check, pts, admin2, comm_coords, locations, hosp_gauls)

#calculate the rpoprotion of dr
mydata$rate <- mydata$number_resistant/mydata$no_examined

#calculate variance as the weight
mydata$weight <- (mydata$rate*(1-mydata$rate))/mydata$no_examined

mydata <- mydata[,.(source_id,
                   iso3,
                   GAUL_ADM1,
                   GAUL_ADM2,
                   mid_year,
                   species,
                   antimicrobial,
                   number_resistant,
                   no_examined,
                   rate,
                   weight)]

colnames(mydata)[1] <- 'nid'

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
colnames(adm1)[colnames(adm1)=='GAUL_ADM1'] <- 'location_code'
adm1$GAUL_ADM2 <-  NULL
adm1 <- adm1[!is.na(adm1$location_code),]

adm2 <- mdr_typhi
adm2$shapefile <- 'admin2013_2'
colnames(adm2)[colnames(adm2)=='GAUL_ADM2'] <- 'location_code'
adm2$GAUL_ADM1 <-  NULL
adm2 <- adm2[!is.na(adm2$location_code),]

write.csv(adm1, 'K:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/mdr_typhi_mbg_adm1.csv', row.names = F)
write.csv(adm2, 'K:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/mdr_typhi_mbg_adm2.csv', row.names = F)
