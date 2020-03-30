#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Create data availability maps #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rm(list = ls())
library(raster)
library(googlesheets)
library(RColorBrewer)
library(sp)
library(sf)
library(ggplot2)

rm(list = ls())

#1. For MDR Typhi first
# create one map for all data
drugbug <- c('MDR_typhi', 'MDR_paratyphi', 'FQNS_typhi', 'FQNS_paratyphi')
adm0_master <- st_read('Z:/AMR/Shapefiles/admin2013_0.shp')
adm1_master <- st_read('Z:/AMR/Shapefiles/admin2013_1.shp')
adm2_master <- st_read('C:/Users/annieb/Desktop/admin2013_2.shp')

adm0_master <- st_simplify(adm0_master, dTolerance = 0.1, preserveTopology = T)
adm1_master <- st_simplify(adm1_master, dTolerance = 0.1, preserveTopology = T)
adm2_master <- st_simplify(adm2_master, dTolerance = 0.1, preserveTopology = T)

background <- adm0_master

for(db in drugbug){ 
  mydata <- readRDS(paste0('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/', db, '.rds')) 
  
  adm0_locs <- mydata[mydata$best_geo == 'national',]
  adm1_locs <- mydata[mydata$best_geo == 'admin1',]
  adm2_locs <- mydata[mydata$best_geo == 'admin2',]
  coords <- data.frame(mydata[mydata$best_geo == 'coords',])
  
  adm0 <- adm0_master[adm0_master$GAUL_CODE %in% adm0_locs$adm0,]
  adm1 <- adm1_master[adm1_master$GAUL_CODE %in% adm1_locs$adm1,]
  adm2 <- adm2_master[adm2_master$GAUL_CODE %in% adm2_locs$adm2,]
  coords <- unique(coords[c('long', 'lat')])
  coords <-  coords[!is.na(coords$lat),]
  
  # get the endmic area
  regions.lookup <- read.csv("Z:/AMR/Misc/GBD_regions//GBD regions lookup.csv")
  
  if(db == 'MDR_paratyphi' | db == 'FQNS_paratyphi'){
    regions.lookup <-  regions.lookup[regions.lookup$region == 'East Asia'|
                                        regions.lookup$region == 'North Africa & Middle East'|
                                        regions.lookup$region == 'South Asia'|
                                        regions.lookup$region == 'Southeast Asia'|
                                        regions.lookup$region == 'Oceania',]
  } else{
    regions.lookup <-  regions.lookup[regions.lookup$region == 'East Asia'|
                                        regions.lookup$region == 'North Africa & Middle East'|
                                        regions.lookup$region == 'South Asia'|
                                        regions.lookup$region == 'Southeast Asia'|
                                        regions.lookup$region =='Central Sub-Saharan Africa'|
                                        regions.lookup$region =='Western Sub-Saharan Africa'|
                                        regions.lookup$region =='Eastern Sub-Saharan Africa'|
                                        regions.lookup$region =='Southern Sub-Saharan Africa'|
                                        regions.lookup$region == 'Oceania',]
  } 
  
  non_endemic <- background[!(background$GAUL_CODE%in%regions.lookup$GAUL_CODE),]
  
  png(paste0('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/data_availability/maps/', db, '.png'),
      height = 15, width = 20, unit = 'cm', res = 300 )
  ggplot()+
    geom_sf(data = background, fill = '#bdbdbd',colour = 'black', size = 0.25)+
    geom_sf(data = adm0, fill = '#7fcdbb',colour = 'black', size = 0.25)+
    geom_sf(data = adm1, fill = '#31a354',colour = 'black', size = 0.25)+
    geom_sf(data = adm2, fill = '#2c7fb8',colour = 'black', size = 0.25)+
    geom_point(data = coords, aes(x = long, y = lat),colour = 'black', fill = '#c51b8a', pch = 21)+
    geom_sf(data = non_endemic, fill = '#636363',colour = 'black', size = 0.25)+
    theme_bw()+
    theme(line = element_blank(),
          axis.text = element_blank())+
    xlab("") + 
    ylab("")  
    
  dev.off()
} 
