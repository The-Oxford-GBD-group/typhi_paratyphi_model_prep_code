#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Merge on facility location info    #
# Select on location per record      #
# (smallest to incaptulate all data) #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rm(list =ls())
library(foreign)
library(googlesheets)
library(raster)

setwd("Z:/AMR/Pathogens/typhi_paratyphi/model_prep")

master_data <- readRDS('clean_data/clean_data.rds')
facilities <-gs_read(ss = gs_title("Prevalence data - Typhi & Paratyphi"), ws = "facility lookup")
facilities <- data.frame(facilities)

#1. Seperate out community data ####
table(master_data$geometry_type)
master_data$geometry_type[master_data$source_id == 4370 & master_data$site_name == 'Kasai'] <- 'national'
matched.data <- master_data[!is.na(master_data$latitude)| !is.na(master_data$gaul_code),]
national_data <- master_data[master_data$geometry_type == 'national',]
mydata <- master_data[which(is.na(master_data$latitude) & is.na(master_data$gaul_code) & master_data$geometry_type!='national'),]
  
#get the location infomation for the hosptials
facilities <- facilities[c("facility_name",
                           "GAUL_ADM2",
                           "GAUL_ADM1",
                           "lat",
                           "long")]

adm1 <- read.dbf('Z:/AMR/Shapefiles/admin2013_1.dbf')
adm2 <- read.dbf('Z:/AMR/Shapefiles/admin2013_2.dbf')
adm2_shp <- shapefile('C:/Users/annieb/Desktop/admin2013_2.shp')

#1. For the hospitals
#a. merge the locations onto each of the hosptitals
temp_data <-  mydata
temp_data <- merge(temp_data, facilities, by.x = c('hospital_name1'), by.y = c('facility_name'), all.x = T, all.y = F)
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM2'] <- 'GAUL_ADM2_1'
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM1'] <- 'GAUL_ADM1_1'
colnames(temp_data)[colnames(temp_data) == 'lat'] <- 'lat_1'
colnames(temp_data)[colnames(temp_data) == 'long'] <- 'long_1'

temp_data <- merge(temp_data, facilities, by.x = c('hospital_name2'), by.y = c('facility_name'), all.x = T, all.y = F)
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM2'] <- 'GAUL_ADM2_2'
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM1'] <- 'GAUL_ADM1_2'
colnames(temp_data)[colnames(temp_data) == 'lat'] <- 'lat_2'
colnames(temp_data)[colnames(temp_data) == 'long'] <- 'long_2'

temp_data <- merge(temp_data, facilities, by.x = c('hospital_name3'), by.y = c('facility_name'), all.x = T, all.y = F)
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM2'] <- 'GAUL_ADM2_3'
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM1'] <- 'GAUL_ADM1_3'
colnames(temp_data)[colnames(temp_data) == 'lat'] <- 'lat_3'
colnames(temp_data)[colnames(temp_data) == 'long'] <- 'long_3'

temp_data <- merge(temp_data, facilities, by.x = c('hospital_name4'), by.y = c('facility_name'), all.x = T, all.y = F)
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM2'] <- 'GAUL_ADM2_4'
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM1'] <- 'GAUL_ADM1_4'
colnames(temp_data)[colnames(temp_data) == 'lat'] <- 'lat_4'
colnames(temp_data)[colnames(temp_data) == 'long'] <- 'long_4'

temp_data <- merge(temp_data, facilities, by.x = c('hospital_name5'), by.y = c('facility_name'), all.x = T, all.y = F)
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM2'] <- 'GAUL_ADM2_5'
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM1'] <- 'GAUL_ADM1_5'
colnames(temp_data)[colnames(temp_data) == 'lat'] <- 'lat_5'
colnames(temp_data)[colnames(temp_data) == 'long'] <- 'long_5'

temp_data <- merge(temp_data, facilities, by.x = c('hospital_name6'), by.y = c('facility_name'), all.x = T, all.y = F)
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM2'] <- 'GAUL_ADM2_6'
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM1'] <- 'GAUL_ADM1_6'
colnames(temp_data)[colnames(temp_data) == 'lat'] <- 'lat_6'
colnames(temp_data)[colnames(temp_data) == 'long'] <- 'long_6'

temp_data <- merge(temp_data, facilities, by.x = c('hospital_name7'), by.y = c('facility_name'), all.x = T, all.y = F)
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM2'] <- 'GAUL_ADM2_7'
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM1'] <- 'GAUL_ADM1_7'
colnames(temp_data)[colnames(temp_data) == 'lat'] <- 'lat_7'
colnames(temp_data)[colnames(temp_data) == 'long'] <- 'long_7'

temp_data <- merge(temp_data, facilities, by.x = c('hospital_name8'), by.y = c('facility_name'), all.x = T, all.y = F)
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM2'] <- 'GAUL_ADM2_8'
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM1'] <- 'GAUL_ADM1_8'
colnames(temp_data)[colnames(temp_data) == 'lat'] <- 'lat_8'
colnames(temp_data)[colnames(temp_data) == 'long'] <- 'long_8'

temp_data <- merge(temp_data, facilities, by.x = c('hospital_name9'), by.y = c('facility_name'), all.x = T, all.y = F)
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM2'] <- 'GAUL_ADM2_9'
colnames(temp_data)[colnames(temp_data) == 'GAUL_ADM1'] <- 'GAUL_ADM1_9'
colnames(temp_data)[colnames(temp_data) == 'lat'] <- 'lat_9'
colnames(temp_data)[colnames(temp_data) == 'long'] <- 'long_9'

#Match each record to the locations which encompase all of the study locations
#add values in if there are not 9 hospitals
temp_data$GAUL_ADM1_2[is.na(temp_data$GAUL_ADM1_2)] <- temp_data$GAUL_ADM1_1[is.na(temp_data$GAUL_ADM1_2)]
temp_data$GAUL_ADM1_3[is.na(temp_data$GAUL_ADM1_3)] <- temp_data$GAUL_ADM1_1[is.na(temp_data$GAUL_ADM1_3)]
temp_data$GAUL_ADM1_4[is.na(temp_data$GAUL_ADM1_4)] <- temp_data$GAUL_ADM1_1[is.na(temp_data$GAUL_ADM1_4)]
temp_data$GAUL_ADM1_5[is.na(temp_data$GAUL_ADM1_5)] <- temp_data$GAUL_ADM1_1[is.na(temp_data$GAUL_ADM1_5)]
temp_data$GAUL_ADM1_6[is.na(temp_data$GAUL_ADM1_6)] <- temp_data$GAUL_ADM1_1[is.na(temp_data$GAUL_ADM1_6)]
temp_data$GAUL_ADM1_7[is.na(temp_data$GAUL_ADM1_7)] <- temp_data$GAUL_ADM1_1[is.na(temp_data$GAUL_ADM1_7)]
temp_data$GAUL_ADM1_8[is.na(temp_data$GAUL_ADM1_8)] <- temp_data$GAUL_ADM1_1[is.na(temp_data$GAUL_ADM1_8)]
temp_data$GAUL_ADM1_9[is.na(temp_data$GAUL_ADM1_9)] <- temp_data$GAUL_ADM1_1[is.na(temp_data$GAUL_ADM1_9)]

temp_data$GAUL_ADM2_2[is.na(temp_data$GAUL_ADM2_2)] <- temp_data$GAUL_ADM2_1[is.na(temp_data$GAUL_ADM2_2)]
temp_data$GAUL_ADM2_3[is.na(temp_data$GAUL_ADM2_3)] <- temp_data$GAUL_ADM2_1[is.na(temp_data$GAUL_ADM2_3)]
temp_data$GAUL_ADM2_4[is.na(temp_data$GAUL_ADM2_4)] <- temp_data$GAUL_ADM2_1[is.na(temp_data$GAUL_ADM2_4)]
temp_data$GAUL_ADM2_5[is.na(temp_data$GAUL_ADM2_5)] <- temp_data$GAUL_ADM2_1[is.na(temp_data$GAUL_ADM2_5)]
temp_data$GAUL_ADM2_6[is.na(temp_data$GAUL_ADM2_6)] <- temp_data$GAUL_ADM2_1[is.na(temp_data$GAUL_ADM2_6)]
temp_data$GAUL_ADM2_7[is.na(temp_data$GAUL_ADM2_7)] <- temp_data$GAUL_ADM2_1[is.na(temp_data$GAUL_ADM2_7)]
temp_data$GAUL_ADM2_8[is.na(temp_data$GAUL_ADM2_8)] <- temp_data$GAUL_ADM2_1[is.na(temp_data$GAUL_ADM2_8)]
temp_data$GAUL_ADM2_9[is.na(temp_data$GAUL_ADM2_9)] <- temp_data$GAUL_ADM2_1[is.na(temp_data$GAUL_ADM2_9)]

#get the admin 1, 2 and country to match all locations within a record
temp_data$final_adm2 <-apply(temp_data, 1, function(x) ifelse(x["GAUL_ADM2_1"]==x["GAUL_ADM2_2"]&
                                                        x["GAUL_ADM2_1"]==x["GAUL_ADM2_3"]&
                                                        x["GAUL_ADM2_1"]==x["GAUL_ADM2_4"]&
                                                        x["GAUL_ADM2_1"]==x["GAUL_ADM2_5"]&
                                                        x["GAUL_ADM2_1"]==x["GAUL_ADM2_6"]&
                                                        x["GAUL_ADM2_1"]==x["GAUL_ADM2_7"]&
                                                        x["GAUL_ADM2_1"]==x["GAUL_ADM2_8"]&
                                                        x["GAUL_ADM2_1"]==x["GAUL_ADM2_9"],
                                                      x["GAUL_ADM2_1"], NA))

temp_data$final_adm1 <-apply(temp_data, 1, function(x) ifelse(x["GAUL_ADM1_1"]==x["GAUL_ADM1_2"]&
                                                        x["GAUL_ADM1_1"]==x["GAUL_ADM1_3"]&
                                                        x["GAUL_ADM1_1"]==x["GAUL_ADM1_4"]&
                                                        x["GAUL_ADM1_1"]==x["GAUL_ADM1_5"]&
                                                        x["GAUL_ADM1_1"]==x["GAUL_ADM1_6"]&
                                                        x["GAUL_ADM1_1"]==x["GAUL_ADM1_7"]&
                                                        x["GAUL_ADM1_1"]==x["GAUL_ADM1_8"]&
                                                        x["GAUL_ADM1_1"]==x["GAUL_ADM1_9"],
                                                      x["GAUL_ADM1_1"], NA))


temp_data$final_adm0 <- lapply(temp_data$GAUL_ADM1_1, function(id, lookup_dataset = adm1){
  ad <- lookup_dataset$PARENT_ID[lookup_dataset$GAUL_CODE == as.numeric(id)]
  return(ad)
})

#if there is only one hospital and they have coordinates then assign these
temp_data$final_lat[!is.na(temp_data$lat_1) & is.na(temp_data$lat_2)] <- temp_data$lat_1[!is.na(temp_data$lat_1) & is.na(temp_data$lat_2)]
temp_data$final_long[!is.na(temp_data$long_1) & is.na(temp_data$long_2)] <- temp_data$long_1[!is.na(temp_data$long_1) & is.na(temp_data$long_2)]

temp_data$final_adm0 <- unlist(temp_data$final_adm0)
temp_data$final_adm1 <- temp_data$final_adm1
temp_data$final_adm2 <- temp_data$final_adm2
temp_data$final_lat <- temp_data$final_lat
temp_data$final_long <- temp_data$final_long

#Create a varaible concatinating variable names
# temp_data$hospital_name <- paste(temp_data$hospital_name1,
#                               temp_data$hospital_name2,
#                               temp_data$hospital_name3,
#                               temp_data$hospital_name4,
#                               temp_data$hospital_name5,
#                               temp_data$hospital_name6,
#                               temp_data$hospital_name7,
#                               temp_data$hospital_name8,
#                               temp_data$hospital_name9, sep = ';')
# 
# temp_data$hospital_name <-  gsub(';NA', '', temp_data$hospital_name)
temp_data <- temp_data[c("row_id",
                         "final_adm0",
                         "final_adm1",
                         "final_adm2" ,                                     
                         "final_lat",
                         "final_long")]

mydata <- merge(mydata, temp_data, by = 'row_id')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#get the lowest level greographies for non hospital data
matched.data$hospital_name <-  NA
matched.data$final_lat <-  NA
matched.data$final_long <-  NA
matched.data$final_adm1 <-  NA
matched.data$final_adm2 <-  NA
matched.data$final_adm0 <-  NA

#get the admin 2 locations of each set of coordinates
coords <- data.frame('long' = matched.data$longitude[!is.na(matched.data$longitude)],
                     'lat' = matched.data$latitude[!is.na(matched.data$longitude)])

coords <- SpatialPoints(coords, proj4string = CRS(projection(adm2_shp)))
matched.data$gaul_code[!is.na(matched.data$longitude)] <- over(coords, adm2_shp)$GAUL_CODE
matched.data$admin_level[!is.na(matched.data$longitude)] <- 2

coords <- data.frame('long' = matched.data$longitude2[!is.na(matched.data$longitude2)],
                     'lat' = matched.data$latitude2[!is.na(matched.data$longitude2)])

coords <- SpatialPoints(coords, proj4string = CRS(projection(adm2_shp)))
matched.data$gaul_code2[!is.na(matched.data$longitude2)] <- over(coords, adm2_shp)$GAUL_CODE
matched.data$admin_level2[!is.na(matched.data$longitude2)] <- 2

#if just one admin location assign those
matched.data$final_adm2[!is.na(matched.data$gaul_code) & is.na(matched.data$gaul_code2) & matched.data$admin_level == 2] <- 
  matched.data$gaul_code[!is.na(matched.data$gaul_code) & is.na(matched.data$gaul_code2) & matched.data$admin_level == 2]

matched.data$final_adm1[!is.na(matched.data$gaul_code) & is.na(matched.data$gaul_code2) & matched.data$admin_level == 1] <- 
  matched.data$gaul_code[!is.na(matched.data$gaul_code) & is.na(matched.data$gaul_code2) & matched.data$admin_level == 1]

#if multiple locations, select the smalest area which satisfies all locs
# for those with 2 locations
#if they (two sets of coordinates) in the same admin 2 
matched.data$final_adm2[!is.na(matched.data$gaul_code2) & is.na(matched.data$gaul_code3)]<-
  apply(matched.data[!is.na(matched.data$gaul_code2) & is.na(matched.data$gaul_code3),],1,function(x) ifelse(x["admin_level"]==2 &
                                                                                                               x["admin_level2"]==2 &
                                                                                                               x["gaul_code"] == x["gaul_code2"],
                                                                                                               x["gaul_code"] , NA)) 
#if 2 admin 2s in the same admin 1
matched.data$final_adm1[!is.na(matched.data$gaul_code2) & is.na(matched.data$gaul_code3)]<-
  apply(matched.data[!is.na(matched.data$gaul_code2) & is.na(matched.data$gaul_code3),],1,function(x) ifelse(x["admin_level"]==2 &
                                            x["admin_level2"]==2 &
                                            adm2$PARENT_ID[adm2$GAUL_CODE == x["gaul_code"]] == adm2$PARENT_ID[adm2$GAUL_CODE == x["gaul_code2"]],
                                            adm2$PARENT_ID[adm2$GAUL_CODE == x["gaul_code"]], NA)) 

#if two admin 2s in the same country
matched.data$final_adm0[!is.na(matched.data$gaul_code2) & is.na(matched.data$gaul_code3)]<-
  apply(matched.data[!is.na(matched.data$gaul_code2) & is.na(matched.data$gaul_code3),],1,function(x) ifelse(x["admin_level"]==2 &
                                                                                                               x["admin_level2"]==2 &
                                                                                                               adm1$PARENT_ID[adm1$GAUL_CODE==(adm2$PARENT_ID[adm2$GAUL_CODE == x["gaul_code"]])] == adm1$PARENT_ID[adm1$GAUL_CODE==(adm2$PARENT_ID[adm2$GAUL_CODE == x["gaul_code2"]])],
                                                                                                             adm1$PARENT_ID[(adm2$PARENT_ID[adm2$GAUL_CODE == x["gaul_code"]])], NA))
# NO ADMIN 1s FIT THIS ATM SO DONT RUN
# #if two admin 1s in the same country
# matched.data$final_adm0[!is.na(matched.data$gaul_code2) & is.na(matched.data$gaul_code3)]<-
#   apply(matched.data[!is.na(matched.data$gaul_code2) & is.na(matched.data$gaul_code3),],1,function(x) ifelse(x["admin_level"]==1 &
#                                                                                                                x["admin_level2"]==1 &
#                                                                                                                adm1$PARENT_ID[adm1$GAUL_CODE == x["gaul_code"]] == adm1$PARENT_ID[adm1$GAUL_CODE == x["gaul_code2"]],
#                                                                                                              adm1$PARENT_ID[adm1$GAUL_CODE == x["gaul_code"]], NA)) 


#for those with 3 locations
#if admin2s are all the same (i.e. there were 3 lots of coordinates from the same admin 2), assign this
matched.data$final_adm2[!is.na(matched.data$gaul_code3)]<-
  apply(matched.data[!is.na(matched.data$gaul_code3),],1,function(x) ifelse(x["admin_level"]==2 &
                                                                             x["admin_level2"]==2 &
                                                                             x["admin_level3"]==2 &
                                                                             x["gaul_code"] == x["gaul_code2"] &
                                                                             x["gaul_code"] == x["gaul_code3"],
                                                                             x["gaul_code"] , NA)) 

# #if the three admin 2 locations are in the same admin 1 assign this (none actually meet this requirement so skip for now)
# matched.data$final_adm1[!is.na(matched.data$gaul_code3)]<-
#   apply(matched.data[!is.na(matched.data$gaul_code3),],1,function(x) ifelse(x["admin_level"]==2 &
#                                                                               x["admin_level2"]==2 &
#                                                                               x["admin_level3"]==2 &
#                                                                               adm2$PARENT_ID[adm2$GAUL_CODE == x["gaul_code"]] == adm2$PARENT_ID[adm2$GAUL_CODE == x["gaul_code2"]]&
#                                                                               adm2$PARENT_ID[adm2$GAUL_CODE == x["gaul_code"]] == adm2$PARENT_ID[adm2$GAUL_CODE == x["gaul_code3"]],
#                                                                             adm2$PARENT_ID[adm2$GAUL_CODE == x["gaul_code"]], NA))

# # if they are admin 2s in differnt admin 1's assign the admin 0 (none of these at the moment so skip for now)
# matched.data$final_adm0[!is.na(matched.data$gaul_code3)]<-
#   apply(matched.data[!is.na(matched.data$gaul_code3),],1,function(x) ifelse(x["admin_level"]==2 &
#                                                                             x["admin_level2"]==2 &
#                                                                             x["admin_level3"]==2 &
#                                                                              adm1$PARENT_ID[adm1$GAUL_CODE==(adm2$PARENT_ID[adm2$GAUL_CODE == x["gaul_code"]])] == adm1$PARENT_ID[adm1$GAUL_CODE==(adm2$PARENT_ID[adm2$GAUL_CODE == x["gaul_code2"]])]&
#                                                                              adm1$PARENT_ID[adm1$GAUL_CODE==(adm2$PARENT_ID[adm2$GAUL_CODE == x["gaul_code"]])] == adm1$PARENT_ID[adm1$GAUL_CODE==(adm2$PARENT_ID[adm2$GAUL_CODE == x["gaul_code3"]])],
#                                                                             adm1$PARENT_ID[(adm2$PARENT_ID[adm2$GAUL_CODE == x["gaul_code"]])], NA)) 

#if they are admin 1s in the same country assign them the admin 0
matched.data$final_adm0[!is.na(matched.data$gaul_code3)]<-
  apply(matched.data[!is.na(matched.data$gaul_code3),],1,function(x) ifelse(x["admin_level"]==1 &
                                                                              x["admin_level2"]==1 &
                                                                              x["admin_level3"]==1 &
                                                                              adm1$PARENT_ID[adm1$GAUL_CODE == x["gaul_code"]] == adm1$PARENT_ID[adm1$GAUL_CODE == x["gaul_code2"]]&
                                                                              adm1$PARENT_ID[adm1$GAUL_CODE == x["gaul_code"]] == adm1$PARENT_ID[adm1$GAUL_CODE == x["gaul_code3"]],
                                                                              adm1$PARENT_ID[adm1$GAUL_CODE == x["gaul_code"]], NA)) 


#fill in any gaps
matched.data$final_adm1[is.na(matched.data$final_adm1) & !is.na(matched.data$final_adm2)] <- 
  lapply(matched.data$final_adm2[is.na(matched.data$final_adm1) & !is.na(matched.data$final_adm2)], 
         function(id, lookup_dataset = adm2){
         ad <- lookup_dataset$PARENT_ID[lookup_dataset$GAUL_CODE == as.numeric(id)]
      return(ad)
    })

matched.data$final_adm0[is.na(matched.data$final_adm0) & !is.na(matched.data$final_adm1)] <- 
  lapply(matched.data$final_adm1[is.na(matched.data$final_adm0) & !is.na(matched.data$final_adm1)], 
         function(id, lookup_dataset = adm1){
           ad <- lookup_dataset$PARENT_ID[lookup_dataset$GAUL_CODE == as.numeric(id)]
           return(ad)
         })

matched.data$final_adm0 <-  unlist(matched.data$final_adm0)
matched.data$final_adm1 <-  unlist(matched.data$final_adm1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Sort out the national only data #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
national_data$hospital_name <-  NA
national_data$final_lat <-  NA
national_data$final_long <-  NA
national_data$final_adm1 <-  NA
national_data$final_adm2 <-  NA
national_data$final_adm0 <-  NA

national_data$final_adm0 <- 
  lapply(national_data$iso3, 
         function(id, lookup_dataset = adm1){
           ad <- unique(lookup_dataset$PARENT_ID[lookup_dataset$COUNTRY_ID == id])
           return(ad)
         })
national_data$final_adm0[national_data$iso3=='IND'] <- 115
national_data$final_adm0[national_data$iso3=='PAK'] <- 188
national_data$final_adm0[national_data$iso3=='KEN'] <- 133
national_data$final_adm0 <-  unlist(national_data$final_adm0)

#join the datasets together
mydata <- rbind(mydata, matched.data)
mydata <- rbind(mydata, national_data)

#define best locations
mydata$best_geo[!is.na(mydata$final_lat)] <- 'coords'
mydata$best_geo[is.na(mydata$final_lat) & !is.na(mydata$final_adm2)] <- 'admin2'
mydata$best_geo[is.na(mydata$final_lat) & is.na(mydata$final_adm2) & !is.na(mydata$final_adm1)] <- 'admin1'
mydata$best_geo[is.na(mydata$final_lat) & is.na(mydata$final_adm2) & is.na(mydata$final_adm1) & !is.na(mydata$final_adm0)] <- 'national'

#get a final admin names
mydata$final_adm2_name <-  NA
mydata$final_adm1_name <- NA
adm2$NAME <-  as.character(adm2$NAME)
adm1$NAME <-  as.character(adm1$NAME)

mydata$final_adm2_name[!is.na(mydata$final_adm2)] <- 
  lapply(mydata$final_adm2[!is.na(mydata$final_adm2)], 
         function(id, lookup_dataset = adm2){
           ad <- lookup_dataset$NAME[lookup_dataset$GAUL_CODE == as.numeric(id)]
           return(ad)
         })

mydata$final_adm1_name[!is.na(mydata$final_adm1)] <- 
  lapply(mydata$final_adm1[!is.na(mydata$final_adm1)], 
         function(id, lookup_dataset = adm1){
           ad <- lookup_dataset$NAME[lookup_dataset$GAUL_CODE == as.numeric(id)]
           return(ad)
         })

mydata$final_adm1_name <-  unlist(mydata$final_adm1_name)
mydata$final_adm2_name <-  unlist(mydata$final_adm2_name)

#limit to required variables
names(mydata)
saveRDS(mydata, 'Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/cleaned_locations.RDS')

locs_by_row <- mydata[c("row_id",
                       "source_id",
                       "start_year",             
                       "end_year",
                       "mid_year",
                       "best_geo",
                       "setting",
                       "super_region",
                       "region",
                       "country",
                       "iso3",
                       "location_name",
                       "geometry_type",
                       "final_adm1_name",
                       "final_adm2_name",
                       "final_adm0",
                       "final_adm1",
                       "final_adm2",
                       "final_lat",
                       "final_long")]

write.csv(locs_by_row, 'Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/locations/locations_by_row.csv', row.names = F)

locs_study_site <- unique(mydata[c("source_id",
                                  "best_geo",
                                  "setting",
                                  "super_region",
                                  "region",
                                  "country",
                                  "iso3",
                                  "final_adm1_name",
                                  "final_adm2_name",
                                  # "site_name",
                                  # 'hospital_name',
                                  "location_name",
                                  "geometry_type",
                                  "final_adm0",
                                  "final_adm1",
                                  "final_adm2",
                                  "final_lat",
                                  "final_long")])

# locs_study_site <-  data.table(locs_study_site)
# locs_study_site$dups <- duplicated(locs_study_site, by = c('source_id', 'site_name', 'geometry_type'))
# table(locs_study_site$dups)
write.csv(locs_study_site, 'Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/locations/locations_by_study_site.csv', row.names = F)
