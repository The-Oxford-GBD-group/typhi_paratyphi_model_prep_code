#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Create data availability map #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rm(list = ls())
library(raster)
library(googlesheets)
library(RColorBrewer)
library(sp)

setwd("K:/AMR/Pathogens/typhi_paratyphi/model_prep//")

# Load in the clean data
master.data <- readRDS("clean_data/clean_data_crosswalked.rds")

#restrict to data in the highly endemic regions
master.data <- master.data[master.data$super_region!="Latin America & Caribbean",]
master.data <- master.data[!(master.data$super_region=="Sub-Saharan Africa" & master.data$species == 'Paratyphi'),]

#clean up antibiotic group
master.data$antimicrobial[master.data$antimicrobial == 'multi-drug' | master.data$antimicrobial == 'multi-drug-adjusted'] <- 'MDR'
master.data$antimicrobial[grep('FQR', master.data$antimicrobial)] <- 'FQR'

#read in required files
#facility list
master.facilities <- gs_read(ss = gs_title("Prevalence data - Typhi & Paratyphi"), ws = "data")
names(master.facilities) <- tolower(names(master.facilities))
#location details
locations <- gs_read(ss = gs_title("Prevalence data - Typhi & Paratyphi"), ws = "facility lookup")
#admin 2 shapefile
admin2 <- shapefile("C:/Users/abrowne/Desktop/Shapefile directory/admin2013_2.shp")
#admin1 shapefile
admin1 <- shapefile("C:/Users/abrowne/Desktop/Shapefile directory/admin2013_1.shp")
#load in admin 0 as a base layer
admin0 <- shapefile("C:/Users/abrowne/desktop/Shapefile directory/admin2013_0.shp")
#load in the endemic area and plot over the base layer
regions.lookup <- read.csv("K:/AMR/Misc/GBD regions lookup.csv")

#Loop through & produce maps for All studies and MDR and FQR typhi and paratyphi
species <- c('All', 'Typhi', 'Paratyphi')
antimicrobial <- c('All', 'MDR', 'FQR')

for(sp in species){
  if(sp == 'All') {
    subset1 <- master.data
  } else {
    subset1 <- master.data[master.data$species == sp]
  }
  for(a in antimicrobial){
  if(a == 'All') {
    mydata <- subset1
  } else {
    mydata <- subset1[subset1$antimicrobial == a]
  }
  
  if(sp=='All' & a !="All"){
  }else if(sp!='All' & a =="All"){
  } else {   

# limit to facilities in the dataset
facilities <- master.facilities[master.facilities$source_id%in%mydata$source_id,]

#Get a list of the hospitals in the dataset
hospitals <- c(facilities$hospital_name1, facilities$hospital_name2, facilities$hospital_name3, facilities$hospital_name4, facilities$hospital_name5, facilities$hospital_name6, facilities$hospital_name7, facilities$hospital_name8, facilities$hospital_name9)
hospitals <- unique(hospitals)

#limit information with the hospitals to those in the dataset
hospitals <- locations[locations$facility_name%in%hospitals,]

#get health facility points
pnt_hospitals <- hospitals[!(is.na(hospitals$lat)),] 
pnt_hospitals <- pnt_hospitals[c("long", "lat")]

#get health facility admin units
admin1.hospitals <- unique(hospitals$GAUL_ADM1[is.na(hospitals$lat) & is.na(hospitals$GAUL_ADM2)])
admin2.hospitals <- unique(hospitals$GAUL_ADM2[is.na(hospitals$lat)])
rm(hospitals)

#get locations with polygons in dataset
comm_polys1 <-  unique(facilities$gaul_code[facilities$geometry_type == 'polygon' &
                                          !is.na(facilities$admin_level) & 
                                          facilities$setting == 'Community' &
                                          facilities$admin_level=='1'])

comm_polys2 <-  unique(facilities$gaul_code[facilities$geometry_type == 'polygon' &
                                              !is.na(facilities$admin_level) & 
                                              facilities$setting == 'Community' &
                                              facilities$admin_level=='2'])


hosp_polys1 <-  unique(facilities$gaul_code[facilities$geometry_type == 'polygon' &
                                              !is.na(facilities$admin_level) & 
                                              facilities$setting != 'Community' &
                                              facilities$admin_level=='1'])

hosp_polys2 <-  unique(facilities$gaul_code[facilities$geometry_type == 'polygon' &
                                              !is.na(facilities$admin_level) & 
                                              facilities$setting != 'Community' &
                                              facilities$admin_level=='2'])


admin1.hospitals <- c(admin1.hospitals, hosp_polys1)
admin2.hospitals <- c(admin2.hospitals, hosp_polys2)
admin2.hospitals <- admin2.hospitals[!is.na(admin2.hospitals)]
rm(hosp_polys1, hosp_polys2)

#get locations with points from the dataset
pnts <- facilities[facilities$geometry_type == 'point',]
comm_pnts <- pnts[pnts$setting == 'Community',]
hosp_pnts <- pnts[pnts$setting != 'Community',]

if(length(mydata$latitude[mydata$setting == 'Community'])>0){
  comm_pnts_lat <- c(comm_pnts$latitude, comm_pnts$latitude2)
  comm_pnts_long <- c(comm_pnts$longitude, comm_pnts$longitude2)
  comm_pnts_lat <- unique(comm_pnts_lat)
  comm_pnts_long <- unique(comm_pnts_long)
  comm_pnts_lat <- as.numeric(comm_pnts_lat[!is.na(comm_pnts_lat)])
  comm_pnts_long <- as.numeric(comm_pnts_long[!is.na(comm_pnts_long)])
  comm_pnts <- SpatialPoints(cbind(comm_pnts_long, comm_pnts_lat))
}

hosp_pnts_lat <- c(hosp_pnts$latitude, hosp_pnts$latitude2)
hosp_pnts_long <- c(hosp_pnts$longitude, hosp_pnts$longitude2)
hosp_pnts_lat <- unique(hosp_pnts_lat)
hosp_pnts_long <- unique(hosp_pnts_long)
hosp_pnts_lat <- hosp_pnts_lat[!is.na(hosp_pnts_lat)]
hosp_pnts_long <- hosp_pnts_long[!is.na(hosp_pnts_long)]
hosp_pnts_lat <- as.numeric(hosp_pnts_lat)
hosp_pnts_long <- as.numeric(hosp_pnts_long)
hosp_pnts <- cbind(hosp_pnts_long, hosp_pnts_lat)
colnames(hosp_pnts) <- c('long', 'lat')

hosp_pnts <- rbind(pnt_hospitals, hosp_pnts)
hosp_pnts <- SpatialPoints(hosp_pnts)

rm(pnts, comm_pnts_lat, comm_pnts_long, pnt_hospitals, hosp_pnts_lat, hosp_pnts_long)

#Restrict admin2 shapefile to polygons with community and hospital data
comm_polys2 <- admin2[admin2$GAUL_CODE%in%comm_polys2,]
admin2.hospitals <- admin2[admin2$GAUL_CODE%in%admin2.hospitals,]

#Restrict admin1 shapefile to hospital and lab polygons
admin1.hospitals <- admin1[admin1$GAUL_CODE%in%admin1.hospitals,]

#get the endemic regions for typhi and paratyphi
if(sp == 'Typhi' | sp == 'All'){
regions.lookup <-  regions.lookup[regions.lookup$region == 'East Asia'|
                                    regions.lookup$region == 'North Africa & Middle East'|
                                    regions.lookup$region == 'South Asia'|
                                    regions.lookup$region == 'Southeast Asia'|
                                    regions.lookup$region =='Central Sub-Saharan Africa'|
                                    regions.lookup$region =='Western Sub-Saharan Africa'|
                                    regions.lookup$region =='Eastern Sub-Saharan Africa'|
                                    regions.lookup$region == 'Oceania',]
} else{
  regions.lookup <-  regions.lookup[regions.lookup$region == 'East Asia'|
                                      regions.lookup$region == 'North Africa & Middle East'|
                                      regions.lookup$region == 'South Asia'|
                                      regions.lookup$region == 'Southeast Asia'|
                                      regions.lookup$region == 'Oceania',]
} 

endemic <- admin0[(admin0$GAUL_CODE%in%regions.lookup$GAUL_CODE),]

#Add all of these layers to the plot
if(sp == 'All'){
  jpeg(paste0("data_availability/maps/data_availability_map_all_studies.jpeg"),
       width = 30,
       height = 20,
       units = 'cm',
       res = 300)
}else{
jpeg(paste0("data_availability/maps/data_availability_map_", a, "_", sp, ".jpeg"),
     width = 30,
     height = 20,
     units = 'cm',
     res = 300)
}

plot(admin0,
     col = '#969696') 
     #density = 10, 
     #angle = 45
     
plot(endemic, 
     col = "#d9d9d9",
     add = T)

plot(admin1.hospitals,
     add = T,
     col ="#d95f0e")

plot(admin2.hospitals,
     add = T,
     col ="#d95f0e")

plot(comm_polys2,
     add = T,
     col = "#3182bd")

plot(hosp_pnts, 
     add = T, 
     col = "black",
     bg = "#d95f0e",
     pch = 21)

if(length(mydata$latitude[mydata$setting == 'Community'])>0){
plot(comm_pnts, 
     add = T,
     col = "black",
     bg = "#3182bd",
     pch = 21)
}

legend(-170, -40, 
       legend = c("Health facility/laboratory",
                  "Community study",
                  "Endemic area",
                  "Non endemic area"),
       fill = c("#d95f0e",
                "#3182bd",
                "#d9d9d9",
                '#969696'),
       density = c(NA,NA,NA,NA),
       cex = 0.8)
if(a == 'All' & sp == 'All'){
  title('Data availability map - all data')
}else{
  title(paste(a, sp, 'data availability map', sep = " "))
}

dev.off()
   
    rm(hosp_pnts, 
       comm_pnts, comm_polys2, comm_polys1,
       admin1.hospitals, admin2.hospitals,
       endemic, facilities)
    } 
  }
}
