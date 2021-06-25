#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Merge on facility location info    #
# Select on location per record      #
# (smallest to incaptulate all data) #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rm(list =ls())
library(foreign)
library(googlesheets4)
library(raster)
library(data.table)
library(readxl)

setwd("Z:/AMR/Pathogens/typhi_paratyphi/model_prep")

master_data <- readRDS('clean_data/step1_cleaning.rds')

# fix the lists in the master data (messed up from you googlesheets function)
master_data$lower_age <- gsub('NULL', '999', master_data$lower_age) 
master_data$lower_age <-  gsub('999', NA, master_data$lower_age)
master_data$upper_age <- gsub('NULL', '999', master_data$upper_age) 
master_data$upper_age <-  gsub('999', NA, master_data$upper_age)
master_data$mean_age <- gsub('NULL', '999', master_data$mean_age) 
master_data$mean_age <-  gsub('999', NA, master_data$mean_age)
master_data$percent_male <- gsub('NULL', '999', master_data$percent_male) 
master_data$percent_male <-  gsub('999', NA, master_data$percent_male)
master_data$cultures_performed <- gsub('NULL', '-999', master_data$cultures_performed) 
master_data$cultures_performed <-  gsub('-999', NA, master_data$cultures_performed)

facilities <- read_sheet("https://docs.google.com/spreadsheets/d/109tfx8GOnephsbQf8PSrAtpa7vJQTfXIjKZRsFxYwGo/edit#gid=0", sheet = "facility lookup")
facilities$GAUL_ADM1 <- unlist(facilities$GAUL_ADM1)
facilities$lat <-  gsub('NULL', '999', facilities$lat)
facilities$lat <-  gsub('999', NA, facilities$lat)
facilities$lat <- as.numeric(facilities$lat)
# facilities <-gs_read(ss = gs_title("Prevalence data - Typhi & Paratyphi"), ws = "facility lookup")
# facilities <- data.frame(facilities)

#1. Seperate out community data ####
table(master_data$geometry_type)
master_data$geometry_type[master_data$source_id == 4370 & master_data$site_name == 'Kasai'] <- 'national'
matched.data <- master_data[!is.na(master_data$latitude)| !is.na(master_data$gaul_code),]
matched.data <-  matched.data[matched.data$geometry_type != 'national',]
national_data <- master_data[master_data$geometry_type == 'national',]
mydata <- master_data[which(is.na(master_data$latitude) & is.na(master_data$gaul_code) & master_data$geometry_type!='national'),]

#get the location infomation for the hosptials
facilities <- facilities[c("facility_name",
                           "GAUL_ADM2",
                           "GAUL_ADM1",
                           "lat",
                           "long")]

adm1 <- read.dbf('D:/Z_drive//Shapefiles/admin2013_1.dbf')
adm2 <- read.dbf('D:/Z_drive/Shapefiles/admin2013_2.dbf')
adm2_shp <- shapefile('D:/Z_drive/Shapefiles/admin2013_2.shp')

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
saveRDS(mydata, 'clean_data/cleaned_locations.RDS')

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

write.csv(locs_by_row, 'clean_data/locations/locations_by_row.csv', row.names = F, na = "")

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
write.csv(locs_study_site, 'clean_data/locations/locations_by_study_site.csv', row.names = F, na = "")

rm(adm1, adm2, adm2_shp, locs_study_site, locs_by_row, matched.data, national_data, temp_data, facilities)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Clean up the demographic info ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Get rid of the asterisks in age variables
mydata$lower_age <- gsub("[^0-9\\.]", "", mydata$lower_age) 
mydata$lower_age <- round(as.numeric(mydata$lower_age),0)
mydata$upper_age <- gsub("[^0-9\\.]", "", mydata$upper_age) 
mydata$upper_age <- round(as.numeric(mydata$upper_age),0)
mydata$mean_age <- gsub("[^0-9\\.]", "", mydata$mean_age) 
mydata$mean_age <-  as.numeric(mydata$mean_age)
mydata$mean_age <- round(mydata$mean_age, 1) 

#define the age range
mydata$age_range <- paste(mydata$lower_age, mydata$upper_age, sep = "-")
mydata$age_range[mydata$age_range == 'NA-NA'] <- 'No specified age restrictions'
mydata$age_range[mydata$age_range == '0-99'] <- 'No specified age restrictions'
mydata$child_adult <- 'No specified age restrictions/Adults and children'
mydata$child_adult[mydata$upper_age<=18] <- 'Children only'
mydata$child_adult[mydata$lower_age>=15] <- 'Adults only'

#clean inpatients and outpatients
mydata$inpatients_outpatients[mydata$source_id == 4526] <- 'outpatients'
mydata$inpatients_outpatients[mydata$inpatients_outpatients == 'community'] <- 'Community'
mydata$inpatients_outpatients[mydata$inpatients_outpatients == 'inpatients'] <- 'Inpatients'
mydata$inpatients_outpatients[mydata$inpatients_outpatients == 'outpatients'] <- 'Outpatients'
mydata$inpatients_outpatients[mydata$inpatients_outpatients == 'not specified'] <- 'Not specified'
mydata$inpatients_outpatients[mydata$inpatients_outpatients == 'inpatients & outpatients'] <- 'Inpatients & Outpatients'
mydata$inpatients_outpatients[mydata$inpatients_outpatients == 'outpatients & emergency dept'] <- 'Outpatients & Emergency department'

table(mydata$inpatients_outpatients)
table(mydata$age_range)
table(mydata$child_adult)

#remove asterisks for sex ratio
mydata$percent_male <- gsub("[^0-9\\.]", "", mydata$percent_male) 

#attach the notes for demographics
notes <- read.csv('clean_data/demographics/studies_demographics_notes.csv', stringsAsFactors = F)
notes <- notes[c(1, 4)]
mydata <- merge(mydata, notes, by = 'source_id', all.x = T)

#get the demographics per row
demographics_per_row <- mydata[c("row_id",
                                    "age_range" ,            
                                    "mean_age",
                                    "child_adult",
                                    "percent_male",
                                    "demographics_notes")]


# aggregate the demographics to the study level
temp_data <- data.table(mydata)
temp_data$percent_male <- as.numeric(temp_data$percent_male)
temp_data$number_male <- (temp_data$percent_male/100)*temp_data$no_examined 

#some changes to enable correct aggregation
temp_data$child_adult[temp_data$source_id == 465] <- 'No specified age restrictions/Adults and children'
temp_data$mean_age[temp_data$source_id == 4370 & temp_data$site_name == 'Bas-Congo'] <-  16
temp_data$mean_age[temp_data$source_id == 4370 & temp_data$site_name == 'Equateur'] <-  25.15
temp_data$mean_age[temp_data$source_id == 4370 & temp_data$site_name == 'Orientale'] <-  25.15
temp_data$mean_age[temp_data$source_id == 5136] <-  15.38
temp_data$mean_age[temp_data$source_id == 5137] <-  21.79
temp_data$mean_age[temp_data$source_id == 5134] <-  NA
temp_data$inpatients_outpatients[temp_data$source_id == 4370] <- 'Inpatients'
temp_data$mean_age[temp_data$source_id == 5148] <-  NA #(Find this out from the raw dataset later)
temp_data$child_adult[temp_data$source_id == 5148] <- 'No specified age restrictions/Adults and children'
temp_data$inpatients_outpatients[temp_data$source_id == 5146] <- 'Inpatients & Outpatients'



aggregated_demographics <- unique(temp_data[,.(source_id,
                                              location_name,
                                              species,
                                              mid_year,
                                              lower_age,
                                              upper_age,
                                              number_male,
                                              no_examined,
                                              mean_age,
                                              child_adult,
                                              inpatients_outpatients)])

all_demographics <- aggregated_demographics[, .(age_range = paste(min(lower_age), max(upper_age), sep = "-"),
                                                percent_male = round(sum(number_male)/sum(no_examined)*100,0)),
                                            by = c("source_id",
                                                   "location_name",
                                                   "species",
                                                   "inpatients_outpatients",
                                                   "child_adult",
                                                   "mean_age")]

# all_demographics$dups <- duplicated(all_demographics, by = c('source_id', 'species', 'location_name'))
# table(all_demographics$dups)

all_demographics$age_range[all_demographics$age_range == 'NA-NA'] <- 'No specified age restrictions'
all_demographics$age_range[all_demographics$age_range == '0-99'] <- 'No specified age restrictions'

#get rid of the notes which say aggregated by year as done this to the data now
notes$demographics_notes[grep('Demographics for all serotypes and years combined', notes$demographics_notes)] <- '*Demographics for all serotypes combined'
notes$demographics_notes[grep('Age demorgaphics are for all years combined; sex demographics are for all serotypes and years combined', notes$demographics_notes)] <- 
  '*Sex demographics for all serotypes combined'
notes <- notes[!grepl('Demographics for all years combined', notes$demographics_notes),]

all_demographics <- merge(all_demographics, notes, by = 'source_id', all.x = T) 

#remove mean ages which cannot be aggregated
aggregated_demographics$mean_age[aggregated_demographics$source_id == 207] <-  NA
aggregated_demographics$mean_age[aggregated_demographics$source_id == 883] <-  NA
aggregated_demographics$mean_age[aggregated_demographics$source_id == 959] <-  NA
aggregated_demographics$child_adult[aggregated_demographics$source_id == 1095] <-  'No specified age restrictions/Adults and children'
aggregated_demographics$inpatients_outpatients[aggregated_demographics$source_id == 2103] <-  'Not specified'
aggregated_demographics$child_adult[aggregated_demographics$source_id == 2103] <-  'No specified age restrictions/Adults and children'

combined_t_pt <- aggregated_demographics[, .(age_range = paste(min(lower_age), max(upper_age), sep = "-"),
                                             percent_male = round(sum(number_male)/sum(no_examined)*100,0)),
                                         by = c("source_id",
                                                "location_name",
                                                "inpatients_outpatients",
                                                "child_adult",
                                                "mean_age")]

# combined_t_pt$dups <- duplicated(combined_t_pt, by = c('source_id', 'location_name'))
# table(combined_t_pt$dups)

combined_t_pt$age_range[combined_t_pt$age_range == 'NA-NA'] <- 'No specified age restrictions'
combined_t_pt$age_range[combined_t_pt$age_range == '0-99'] <- 'No specified age restrictions'
notes <- notes[!grepl('serotype', notes$demographics_notes),]
combined_t_pt <- merge(combined_t_pt, notes, by = 'source_id', all.x = T) 

#Aggregate so one per study
aggregated_demographics$child_adult[aggregated_demographics$source_id == 2946] <-  'No specified age restrictions/Adults and children'
aggregated_demographics$child_adult[aggregated_demographics$source_id == 5065] <-  'No specified age restrictions/Adults and children'
aggregated_demographics$child_adult[aggregated_demographics$source_id == 4370] <-  'No specified age restrictions/Adults and children'
aggregated_demographics$mean_age[aggregated_demographics$source_id == 4370] <-  NA
aggregated_demographics$child_adult[aggregated_demographics$source_id == 2180] <-  'No specified age restrictions/Adults and children'
aggregated_demographics$mean_age[aggregated_demographics$source_id == 2180] <-  NA
aggregated_demographics$mean_age[aggregated_demographics$source_id == 2946] <-  NA

all_demogs <- aggregated_demographics[, .(age_range = paste(min(lower_age), max(upper_age), sep = "-"),
                                          percent_male = round(sum(number_male)/sum(no_examined)*100,0)),
                                      by = c("source_id",
                                             "inpatients_outpatients",
                                             "child_adult",
                                             "mean_age")]

# all_demogs$dups <- duplicated(all_demogs, by = c('source_id'))
# table(all_demogs$dups)

all_demogs$age_range[all_demogs$age_range == 'NA-NA'] <- 'No specified age restrictions'
all_demogs$age_range[all_demogs$age_range == '0-99'] <- 'No specified age restrictions'
all_demogs <- merge(all_demogs, notes, by = 'source_id', all.x = T) 

demographics_per_row <- mydata[c("row_id",
                                 "age_range" ,            
                                 "mean_age",
                                 "child_adult",
                                 "percent_male",
                                 "demographics_notes")]

#save the demographics datasets
write.csv(demographics_per_row, 'clean_data/demographics/demographics_per_row.csv', row.names = F, na = "")
write.csv(all_demographics, 'clean_data/demographics/demographics_by_study_site_serotype.csv', row.names = F, na = "")
write.csv(combined_t_pt, 'clean_data/demographics/demographics_by_study_site.csv', row.names = F, na = "")
write.csv(all_demogs, 'clean_data/demographics/demographics_by_study.csv', row.names = F, na = "")

rm(demographics_per_row, all_demographics, combined_t_pt, notes, temp_data, all_demogs, aggregated_demographics)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Clean the details of the testing ####
# 1. Cultures performed
# 2. Isolate_source
# 3. Testing method
# 4. Guidelines
# 5. Control strain
# 6. ISO accreditation
# 7. QA score
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mydata <- data.table(mydata)

#1. Cultures performed
#remove asterisks from cultures performed
mydata$cultures_performed <- gsub("[^0-9\\.]", "", mydata$cultures_performed)
mydata$cultures_performed <- as.numeric(mydata$cultures_performed)
N_cultures <- unique(mydata[,.(source_id, cultures_performed, location_name)])

# this study has the aggregated and individual results, keep just the agreggated number
N_cultures$cultures_performed[N_cultures$Source_id==2038 & N_cultures$cultures_performed!=9634] <- 0
N_cultures$cultures_performed[N_cultures$Source_id==2093 & N_cultures$cultures_performed!=6956] <- 0

# all other duplicates are due to number of cultures spread over years or sites so sum these
N_cultures_source <- N_cultures[, .(N_cultures = sum(cultures_performed)),
                         by = c('source_id')]

N_cultures_source_site <- N_cultures[, .(N_cultures = sum(cultures_performed)),
                                by = c('source_id', 'location_name')]


rm(N_cultures)

#2. Source of isolate
unique(mydata$diagnostic_method)
mydata$isolate_source <- 'Mixed sites'
mydata$isolate_source[mydata$diagnostic_method == "Blood culture"] <- 'Normally sterile sites'
mydata$isolate_source[mydata$diagnostic_method == "Blood and bone marrow culture"] <- 'Normally sterile sites'
mydata$isolate_source[mydata$diagnostic_method == "Blood and CSF culture"] <- 'Normally sterile sites'
mydata$isolate_source[mydata$diagnostic_method == "Bile culture"] <- 'Normally sterile sites'
mydata$isolate_source[mydata$diagnostic_method == "Stool culture"] <- 'Non sterile sites'
mydata$isolate_source[mydata$diagnostic_method == "Culture (unknown)"] <- 'Unknown'
mydata$isolate_source[mydata$diagnostic_method == "Unknown"] <- 'Unknown'
table(mydata$isolate_source)

isolate_source_site <- unique(mydata, by = c('source_id', 'location_name', 'isolate_source'))

isolate_source_site$dup <-  duplicated(isolate_source_site, by = c('source_id', 'location_name'))

for (i in isolate_source_site$source_id[isolate_source_site$dup == TRUE]){
  isolate_source_site$isolate_source[isolate_source_site$source_id == i] <-  'Mixed sites'
}

isolate_source_site <- unique(isolate_source_site, by = c('source_id', 'location_name', 'isolate_source'))
isolate_source_site <-  isolate_source_site[,.(source_id, location_name, isolate_source)]

isolate_source_study <- unique(isolate_source_site, by = c('source_id', 'isolate_source'))
isolate_source_study <-  isolate_source_study[,.(source_id, isolate_source)]

#3. AMR testing method
unique(mydata$amr_test)

mydata$amr_test[mydata$amr_test == 'disk-diffusion'] <- 'Disk-diffusion'
mydata$amr_test[mydata$amr_test == "DIsk-diffusion; VITEK 2"] <- "Disk-diffusion; VITEK 2" 
mydata$amr_test[mydata$amr_test == 'Etest'] <- 'E-test'
mydata$amr_test[mydata$amr_test == 'Plasmid resistance'] <- 'Unknown'
mydata$amr_test[mydata$amr_test == 'Disk-diffusion; unspecified MIC testing'] <- 'Disk-diffusion'

ast <- unique(mydata, by = c('source_id', 'location_name', 'amr_test'))
ast <-  ast[order(ast$amr_test)]

ast <- dcast(ast, source_id+location_name~amr_test, value.var = 'amr_test')
names(ast)

ast$`Disk-diffusion`[ast$`Disk-diffusion; E-test`=='Disk-diffusion; E-test'] <- 'Disk-diffusion'
ast$`E-test`[ast$`Disk-diffusion; E-test`=='Disk-diffusion; E-test'] <- 'E-test'
ast$`Disk-diffusion; E-test` <-  NULL

ast$`Disk-diffusion`[ast$`Disk-diffusion; Microdilution`=='Disk-diffusion; Microdilution'] <-  'Disk-diffusion'
ast$`Microdilution`[ast$`Disk-diffusion; Microdilution`=='Disk-diffusion; Microdilution'] <-  'Microdilution'
ast$`Disk-diffusion; Microdilution` <-  NULL

ast$`Disk-diffusion`[ast$`Disk-diffusion; Phoenix 100 automated system`=='Disk-diffusion; Phoenix 100 automated system'] <- 'Disk-diffusion'
ast$`Phoenix 100 automated system`[ast$`Disk-diffusion/Phoenix 100 automated system`=='Disk-diffusion; Phoenix 100 automated system'] <- 'Phoenix 100 automated system'
ast$`Disk-diffusion; Phoenix 100 automated system` <-  NULL

ast$`Disk-diffusion`[ast$`Disk-diffusion; VITEK 2` == 'Disk-diffusion; VITEK 2'] <- 'Disk-diffusion'
ast$`Automated methods`[ast$`Disk-diffusion; VITEK 2` == 'Disk-diffusion; VITEK 2'] <- 'Automated methods'
ast$`Disk-diffusion; VITEK 2` <-  NULL

ast$`Agar dilution`[ast$`Agar dilution; VITEK 2` == 'Agar dilution; VITEK 2'] <- 'Agar dilution'
ast$`VITEK 2`[ast$`Agar dilution; VITEK 2` == 'Agar dilution; VITEK 2'] <- 'VITEK 2'
ast$`Agar dilution; VITEK 2` <-  NULL

ast$Unknown <-  NULL

ast$`Dilution Assay` <- NA
ast$`Dilution Assay`[ast$`Microdilution` == 'Microdilution'] <- 'Microdilution'
ast$`Dilution Assay`[ast$`Agar dilution` == 'Agar dilution'] <- 'Microdilution'
ast$`Agar dilution` <-  NULL
ast$`Microdilution` <- NULL

ast$`Automated methods` <- NA
ast$`Automated methods`[ast$`Phoenix 100 automated system` == 'Phoenix 100 automated system'] <- 'Automated methods'
ast$`Automated methods`[ast$`Rapid ATB test` == 'Rapid ATB test'] <- 'Automated methods'
ast$`Automated methods`[ast$`VITEK 2` == 'VITEK 2'] <- 'Automated methods'
ast$`Phoenix 100 automated system` <- NULL
ast$`Rapid ATB test` <- NULL
ast$`VITEK 2` <- NULL

ast$ast <- paste(ast$`Disk-diffusion`, ast$`E-test`, ast$`Dilution Assay`, ast$`Automated methods`, sep = ' & ')
ast$ast <- gsub(' & NA', '', ast$ast)
ast$ast <- gsub('NA & ', '', ast$ast)

ast <-  ast[,.(source_id, location_name, ast)]
ast$ast[ast$ast=='NA'] <- 'Not specified'
table(ast$ast)

# ast$ast[ast$ast == 'Disk-diffusion & E-test'] <- 'Disk-diffusion & MIC determination'
# ast$ast[ast$ast == 'Disk-diffusion & Automated methods'] <- 'Disk-diffusion & MIC determination'
# ast$ast[ast$ast == 'Disk-diffusion & Microdilution'] <- 'Disk-diffusion & MIC determination'
# ast$ast[ast$ast == 'Disk-diffusion & E-test & Microdilution'] <- 'Disk-diffusion & MIC determination'
# 
# ast$ast[ast$ast == 'Microdilution & Automated methods'] <- 'Multiple MIC determination methods'
# ast$ast[ast$ast == 'E-test & Automated methods'] <- 'Multiple MIC determination methods'
# ast$ast[ast$ast == 'E-test & Microdilution'] <- 'Multiple MIC determination methods'


ast_study <- unique(mydata, by = c('source_id', 'amr_test'))
ast_study <-  ast_study[order(ast_study$amr_test)]

ast_study <- dcast(ast_study, source_id~amr_test, value.var = 'amr_test')
names(ast_study)

ast_study$`Disk-diffusion`[ast_study$`Disk-diffusion; E-test`=='Disk-diffusion; E-test'] <- 'Disk-diffusion'
ast_study$`E-test`[ast_study$`Disk-diffusion; E-test`=='Disk-diffusion; E-test'] <- 'E-test'
ast_study$`Disk-diffusion; E-test` <-  NULL

ast_study$`Disk-diffusion`[ast_study$`Disk-diffusion; Microdilution`=='Disk-diffusion; Microdilution'] <-  'Disk-diffusion'
ast_study$`Microdilution`[ast_study$`Disk-diffusion; Microdilution`=='Disk-diffusion; Microdilution'] <-  'Microdilution'
ast_study$`Disk-diffusion; Microdilution` <-  NULL

ast_study$`Disk-diffusion`[ast_study$`Disk-diffusion; Phoenix 100 automated system`=='Disk-diffusion; Phoenix 100 automated system'] <- 'Disk-diffusion'
ast_study$`Phoenix 100 automated system`[ast_study$`Disk-diffusion/Phoenix 100 automated system`=='Disk-diffusion; Phoenix 100 automated system'] <- 'Phoenix 100 automated system'
ast_study$`Disk-diffusion; Phoenix 100 automated system` <-  NULL

ast_study$`Disk-diffusion`[ast_study$`Disk-diffusion; VITEK 2` == 'Disk-diffusion; VITEK 2'] <- 'Disk-diffusion'
ast_study$`Automated methods`[ast_study$`Disk-diffusion; VITEK 2` == 'Disk-diffusion; VITEK 2'] <- 'Automated methods'
ast_study$`Disk-diffusion; VITEK 2` <-  NULL

ast_study$`Agar dilution`[ast_study$`Agar dilution; VITEK 2` == 'Agar dilution; VITEK 2'] <- 'Agar dilution'
ast_study$`VITEK 2`[ast_study$`Agar dilution; VITEK 2` == 'Agar dilution; VITEK 2'] <- 'VITEK 2'
ast_study$`Agar dilution; VITEK 2` <-  NULL

ast_study$Unknown <-  NULL

ast_study$`Dilution Assay` <- NA
ast_study$`Dilution Assay`[ast_study$`Microdilution` == 'Microdilution'] <- 'Microdilution'
ast_study$`Dilution Assay`[ast_study$`Agar dilution` == 'Agar dilution'] <- 'Microdilution'
ast_study$`Agar dilution` <-  NULL
ast_study$`Microdilution` <- NULL

ast_study$`Automated methods` <- NA
ast_study$`Automated methods`[ast_study$`Phoenix 100 automated system` == 'Phoenix 100 automated system'] <- 'Automated methods'
ast_study$`Automated methods`[ast_study$`Rapid ATB test` == 'Rapid ATB test'] <- 'Automated methods'
ast_study$`Automated methods`[ast_study$`VITEK 2` == 'VITEK 2'] <- 'Automated methods'
ast_study$`Phoenix 100 automated system` <- NULL
ast_study$`Rapid ATB test` <- NULL
ast_study$`VITEK 2` <- NULL

ast_study$ast <- paste(ast_study$`Disk-diffusion`, ast_study$`E-test`, ast_study$`Dilution Assay`, ast_study$`Automated methods`, sep = ' & ')
ast_study$ast <- gsub(' & NA', '', ast_study$ast)
ast_study$ast <- gsub('NA & ', '', ast_study$ast)

ast_study <-  ast_study[,.(source_id, ast)]
ast_study$ast[ast_study$ast=='NA'] <- 'Not specified'
table(ast_study$ast)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Resistanbce breakpoints ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
guidelines <- unique(mydata, by = c('source_id', 'location_name', 'resistance_breakpoints'))
guidelines <- guidelines[guidelines$antimicrobial!= "FQR - cipro & nalidixic acid",]
guidelines <- guidelines[guidelines$antimicrobial!=  "FQR - ciprofloxacin",]

length(unique(guidelines$source_id))
#remove imputed lines
guidelines$resistance_breakpoints[guidelines$source_id == 117] <- 'CLSI 2012; for azithromycin R \u2265 16 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 135] <- 'NCCLS 1997; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 137] <- 'Working party, 1991 & Ling et al. 1990; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 150] <- 'CLSI 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 151] <- 'CLSI 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 152] <- 'For AMP, CHL & SXT NCCLS 2000; For CIP & NAL CLSI 2005'
guidelines$resistance_breakpoints[guidelines$source_id == 207] <- 'CLSI 2012; for azithromicin R \u2264 18mm'
guidelines$resistance_breakpoints[guidelines$source_id == 235] <- 'BSAC (no year stated); for azithromicin R > 16 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 244] <- 'CLSI 2013'
guidelines$resistance_breakpoints[guidelines$source_id == 248] <- 'CLSI 2013; for azithromicin R > 16 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 267] <- 'For ciprofloxacin R \u2265 0.125 \u03BCg/ml; not stated for other antimicrobials'
guidelines$resistance_breakpoints[guidelines$source_id == 316] <- 'CLSI 2012; for azithromycin R \u2265 16 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 387] <- 'CLSI 2009'
guidelines$resistance_breakpoints[guidelines$source_id == 405] <- 'CLSI 2006; for ciprofloxacin R \u2265 0.19 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 394] <- 'CLSI 2012'
# guidelines$resistance_breakpoints[guidelines$source_id == 548] <- 'For ciprofloxacin R \u2265 0.125 \u03BCg/ml; not stated for other antimicrobials'
# guidelines$resistance_breakpoints[guidelines$source_id == 558] <- 'For ciprofloxacin R \u2265 0.125 \u03BCg/ml; not stated for other antimicrobials'
# guidelines$resistance_breakpoints[guidelines$source_id == 623] <- 'NCCLS (no year stated); for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 639] <- 'For cefepime R \u2264 16mm; for cefixime R \u2264 14mm'
# guidelines$resistance_breakpoints[guidelines$source_id == 693] <- 'For ciprofloxacin R \u2265 0.125 \u03BCg/ml; not stated for other antimicrobials'
guidelines$resistance_breakpoints[guidelines$source_id == 696] <- 'CLSI 2008'
guidelines$resistance_breakpoints[guidelines$source_id == 712] <- 'CLSI (no year stated); ciprofloxacin interpreted as R \u2265 0.125 \u03BCg/ml as MICs provided'
guidelines$resistance_breakpoints[guidelines$source_id == 732] <- 'For aztreonam R \u2264 27mm; ceftazidime R \u2264 22mm; for ceftriaxone R \u2264 25mm'
guidelines$resistance_breakpoints[guidelines$source_id == 737] <- 'NCCLS 2004; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 781] <- 'For ciprofloxacin R \u2265 0.125 \u03BCg/ml; not stated for other antimicrobials'
guidelines$resistance_breakpoints[guidelines$source_id == 787] <- 'CLSI 2005; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 797] <- 'NCCLS (no year stated); ciprofloxacin interpreted as R \u2265 0.125 \u03BCg/ml as MICs provided'
# guidelines$resistance_breakpoints[guidelines$source_id == 833] <- 'NCCLS 2002; ciprofloxacin interpreted as R \u2265 0.125 \u03BCg/ml as MICs provided'
guidelines$resistance_breakpoints[guidelines$source_id == 849] <- 'NCCLS (no year stated); decreased ciprofloxacin sensitivity as \u2265 0.064 \u03BCg/ml as MICs provided'
guidelines$resistance_breakpoints[guidelines$source_id == 850] <- 'CLSI 2003; for azithromycin R \u2265 32 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 852] <- 'CLSI 2009; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 855] <- 'CLSI 2005; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 891] <- 'CLSI 2012; for azithromycin R > 16 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 922] <- 'For ciprofloxacin R \u2265 0.5 \u03BCg/ml; For nalidixic acid R \u2265 32 \u03BCg/ml; not stated for other antimicrobials'
guidelines$resistance_breakpoints[guidelines$source_id == 923] <- 'CLSI 2013; for azithromycin R >16 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 929] <- 'CLSI 2007; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 934] <- 'CLSI 2004'
guidelines$resistance_breakpoints[guidelines$source_id == 1000] <- 'For ampicillin & norfloxacin R \u2264 17mm; for chloramphenicol R \u2264 18mm; for tetracycline R \u2264 19mm; not stated fpr cefotaxime'
guidelines$resistance_breakpoints[guidelines$source_id == 1062] <- 'CLSI 2005; for ciprofloxacin R \u2265 1 \u03BCg/ml; for gatifloxacin R \u2265 8 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 1064] <- 'NCCLS 2002; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1082] <- 'CLSI 2013; for azithrmoycin R > 16 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1095] <- 'CLSI 2007; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1103] <- 'CLSI 2014; EUCAST 2014 for various antimicrobials'
guidelines$resistance_breakpoints[guidelines$source_id == 1126] <- 'CLSI (no year stated); for ciprofloxacin R \u2265 1 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1137] <- 'CLSI 2009; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1223] <- 'CLSI 2008; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 1241] <- 'CLSI 2006; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 1243] <- 'NCCLS 2000; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 1277] <- 'CLSI 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 1334] <- 'NCCLS (no year stated); for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1363] <- 'NCCLS (no year stated); for ciprofloxacin R \u2265 0.25 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1412] <- 'CLSI 2006; For ciprofloxacin R \u2265 0.125 \u03BCg/ml; For azithromycin R \u2265 32 \u03BCg/ml; not stated for other antimicrobials'
guidelines$resistance_breakpoints[guidelines$source_id == 1424] <- 'CLSI 2011; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 1444] <- 'CLSI 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 1465] <- 'CLSI 2007; for azithrmoycin R \u2265 16 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 1466] <- 'CLSI 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 1467] <- 'NCCLS 1999; for ciprofloxacin R > 0.03 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1493] <- 'CLSI (no year stated); for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1494] <- 'EUCAST 2012; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1531] <- 'NCCLS (no year stated); for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1532] <- 'NCCLS 1997; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1533] <- 'NCCLS (no year stated); for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1612] <- 'CLSI 2006; For ciprofloxacin R \u2265 0.125 \u03BCg/ml; For azithromycin R \u2265 32 \u03BCg/ml; not stated for other antimicrobials'
guidelines$resistance_breakpoints[guidelines$source_id == 1694] <- 'NCCLS 2000'
guidelines$resistance_breakpoints[guidelines$source_id == 1698] <- 'CLSI (no year stated); for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1730] <- 'CLSI 2012; for ceftriaxone R >2 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1780] <- 'EUCAST; for azithromycin R >16 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1843] <- 'CLSI (no year stated); for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1861] <- 'NCCLS 1998; for azithromycin R \u2265 8 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1875] <- 'CLSI 2009; for azithromycin R >16 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1887] <- 'NCCLS (no year stated); for ciprofloxacin R \u2265 4 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1917] <- 'For ciprofloxacin R\u2265 0.125 \u03BCg/ml; for gentamicin R \u2265 8 \u03BCg/ml; not stated for other antimicrobials'
guidelines$resistance_breakpoints[guidelines$source_id == 1995] <- 'NCCLS 2001; for ciprofloxacin R \u2265 4 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 2038] <- 'CLSI 2012; EUCAST 2014 for ciprofloxacin'
guidelines$resistance_breakpoints[guidelines$source_id == 2069] <- 'NCCLS 1999; for ciprofloxacin R \u2265 4 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 2165] <- 'CLSI 2008; for ciprofloxacin R \u2265 4 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 2166] <- 'CLSI 2013; for azithromycin R > 16 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 2334] <- 'CLSI 2009; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 2376] <- 'NCCLS (no year stated); for ciprofloxacin R \u2265 4 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 2392] <- 'CLSI 2010; for azithromycin R \u2264 15mm'
guidelines$resistance_breakpoints[guidelines$source_id == 2401] <- 'NCCLS 2002; for ciprofloxacin and ofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 2824] <- 'CLSI 2009; CLSI 2012 for fluoroquinolones'
guidelines$resistance_breakpoints[guidelines$source_id == 2967] <- 'BSAC/EUCAST 2007-2016; for ciprofloxacin R \u2265 0.125 \u03BCg/ml; for ceftriaxone R \u2265 2 \u03BCg/ml; for azithromycin R \u2265 32 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 3012] <- 'CLSI 2012; for azithromycin BSAC 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 3113] <- 'CLSI 2008 and BSAC 2014 for various observations'
guidelines$resistance_breakpoints[guidelines$source_id == 3182] <- 'CLSI 2012; for gatifloxacin R > 1\u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 4050] <- 'CLSI 2012; CLSI 2013 for moxifloxacin and sparfloxacin'
guidelines$resistance_breakpoints[guidelines$source_id == 4229] <- 'CLSI 2012; for streptomycin R \u2265 32\u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 4363] <- 'CLSI 2012; for azithromycin BSAC 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 4370] <- 'CLSI (no year stated); for azithromycin BSAC (no year stated)'
guidelines$resistance_breakpoints[guidelines$source_id == 4490] <- 'CLSI (no year stated)'
guidelines$resistance_breakpoints[guidelines$source_id == 4515] <- 'CLSI 2013; for azithromycin EUCAST 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 4597] <- 'CLSI 2015-2016'
guidelines$resistance_breakpoints[guidelines$source_id == 5054] <- 'CLSI 2013; ciprofloxacin data are inconsistent with these guidelines and appears to be R  \u2265 0.5 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 5064] <- 'Pre 2012 CLSI (no year stated); 2012 and later CLSI 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 5079] <- 'CLSI 2016; for ciprofloxacin R \u2265 0.5 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 5110] <- 'CLSI 2012-2014; Pre 2012 interepreded using CLSI 2011; for azithrmoycin R \u2265 16 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 5124] <- 'CLSI 2012; CLSI (no year stated) for data pre 2012; for azithrmoycin R > 16 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 5148] <- 'CLSI 2000-2018'

guidelines_study <- unique(guidelines, by = c('source_id', 'resistance_breakpoints'))
guidelines_study <- guidelines_study[order(guidelines_study$source_id, guidelines_study$resistance_breakpoints),]
guidelines_study$dups <-  duplicated(guidelines_study, by = c('source_id'))
guidelines_study <- guidelines_study[!(guidelines_study$dups == TRUE & guidelines_study$resistance_breakpoints == 'not specified'),]
guidelines_study$dups <- duplicated(guidelines_study$source_id)
guidelines_study <- guidelines_study[!(guidelines_study$dups == TRUE & guidelines_study$resistance_breakpoints == 'Other'),]
guidelines_study <- unique(guidelines_study, by = c('source_id', 'resistance_breakpoints'))

guidelines_study$dups <- duplicated(guidelines_study$source_id)
unique(guidelines_study$source_id[guidelines_study$dups == TRUE])
guidelines_study <- guidelines_study[,.(source_id, resistance_breakpoints)]

guidelines_study_site <- unique(guidelines, by = c('source_id','location_name', 'resistance_breakpoints'))
guidelines_study_site <- guidelines_study_site[order(guidelines_study_site$source_id, guidelines_study_site$location_name, guidelines_study_site$resistance_breakpoints),]
guidelines_study_site$dups <-  duplicated(guidelines_study_site, by = c('source_id', 'location_name'))
guidelines_study_site <- guidelines_study_site[!(guidelines_study_site$dups == TRUE & guidelines_study_site$resistance_breakpoints == 'not specified'),]
guidelines_study_site$dups <- duplicated(guidelines_study_site$source_id)
guidelines_study_site <- guidelines_study_site[!(guidelines_study_site$dups == TRUE & guidelines_study_site$resistance_breakpoints == 'Other'),]
guidelines_study_site <- unique(guidelines_study_site, by = c('source_id', 'location_name', 'resistance_breakpoints'))

guidelines_study_site$dups <- duplicated(guidelines_study_site, by = c('source_id', 'location_name'))
unique(guidelines_study_site$source_id[guidelines_study_site$dups == TRUE])
guidelines_study_site <- guidelines_study_site[,.(source_id, location_name, resistance_breakpoints)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Control strain used ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
control_strain <- unique(mydata[,.(source_id, location_name, control_strain)])

control_strain$control_strain[control_strain$source_id==4370] <- 'E. coli ATCC 25922; K pneumoniae 700603; S aureus 29213'
control_strain$control_strain[control_strain$source_id==997] <- 'E coli NCTC 10418; S aureaus NCTC 6571'
control_strain$control_strain[control_strain$source_id==548] <- 'E coli ATCC 25922'
control_strain$control_strain[control_strain$source_id==558] <- 'E coli ATCC 25922'
control_strain$control_strain[control_strain$source_id==2999] <- 'E coli ATCC 25922; E faecalis ATCC 29212'
control_strain$control_strain[control_strain$source_id==2093] <- 'E coli NCTC 10418'

control_strain_study_site <- unique(control_strain[,.(source_id, location_name, control_strain)])
control_strain_study <- unique(control_strain[,.(source_id, control_strain)])

#~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. ISO accreditation ####
#~~~~~~~~~~~~~~~~~~~~~~~~~#
facilities <- read_sheet("https://docs.google.com/spreadsheets/d/109tfx8GOnephsbQf8PSrAtpa7vJQTfXIjKZRsFxYwGo/edit#gid=0", sheet = "facility lookup")

accredited <- facilities$facility_name[grep('YES', facilities$`ISO accreditation`)]

mydata$ISO_accredited <-  'No'
mydata$ISO_accredited[mydata$hospital_name1 %in% accredited |
                        mydata$hospital_name2 %in% accredited |
                        mydata$hospital_name3 %in% accredited |
                        mydata$hospital_name4 %in% accredited |
                        mydata$hospital_name5 %in% accredited |
                        mydata$hospital_name6 %in% accredited |
                        mydata$hospital_name7 %in% accredited |
                        mydata$hospital_name8 %in% accredited] <- 'Yes'

ISO_study_site <- unique(mydata[,.(source_id, location_name, ISO_accredited)])
ISO_study_site <- ISO_study_site[order(-ISO_study_site$ISO_accredited),]
ISO_study_site$dups <- duplicated(ISO_study_site, by = c('source_id', 'location_name'))
ISO_study_site <-  ISO_study_site[ISO_study_site$dups == FALSE]
ISO_study_site$dups <-  NULL

ISO_study <- unique(mydata[,.(source_id, ISO_accredited)])
ISO_study <- ISO_study[order(-ISO_study$ISO_accredited),]
ISO_study$dups <- duplicated(ISO_study$source_id)
ISO_study <-  ISO_study[ISO_study$dups == FALSE]
ISO_study$dups <-  NULL

rm(facilities, accredited)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Caclulate quality assessment score         #
# 1 = ISO 15819 accredited lab               #
# 2 = Not accredited but states diagnostic   #
#     method, AST method and guidelines used #
# 3 = Neither 1 nor 2                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mydata$QA <-  NA
mydata$QA[mydata$ISO_accredited == 'Yes'] <- 1
mydata$DM <-  'YES'
mydata$AST <- 'YES'
mydata$guidelines <- 'YES'

mydata$DM[mydata$diagnostic_method == 'Unknown' | mydata$diagnostic_method == 'Culture (Unknown)'] <- 'NO'
mydata$AST[mydata$amr_test == 'Unknown'] <- 'No'
mydata$guidelines[mydata$resistance_breakpoints == 'Not specified'] <- 'NO'

mydata$QA[is.na(mydata$QA) &
            mydata$DM == 'YES' &
            mydata$AST == 'YES' &
            mydata$guidelines == 'YES'] <- 2

mydata$QA[is.na(mydata$QA)] <- 3

QA_study_site <- unique(mydata[,.(source_id, location_name, QA)])
QA_study_site <- QA_study_site[order(QA_study_site$QA)]
QA_study_site$dups <- duplicated(QA_study_site, by = c('source_id', 'location_name'))
QA_study_site <- QA_study_site[QA_study_site$dups==FALSE]
QA_study_site$dups <-  NULL

QA_study <- unique(mydata[,.(source_id, QA)])
QA_study <- QA_study[order(QA_study$QA)]
QA_study$dups <- duplicated(QA_study$source_id)
QA_study <- QA_study[QA_study$dups==FALSE]
QA_study$dups <-  NULL


#merge the testing info together and save the dataset
testing_info_study <- merge(isolate_source_study, ast_study, by = 'source_id')
testing_info_study <-  merge(testing_info_study, N_cultures_source, by = 'source_id')
testing_info_study <-  merge(testing_info_study, guidelines_study, by = 'source_id')
testing_info_study <-  merge(testing_info_study, control_strain_study, by = 'source_id')
testing_info_study <-  merge(testing_info_study, ISO_study, by = 'source_id')
testing_info_study <-  merge(testing_info_study, QA_study, by = 'source_id')
testing_info_study$resistance_breakpoints[testing_info_study$resistance_breakpoints == 'not specified'] <-  'Not specified'
testing_info_study$control_strain[testing_info_study$control_strain == 'not specified'] <-  'Not specified'
testing_info_study$isolate_source[testing_info_study$isolate_source == 'Unknown'] <-  'Not specified'
testing_info_study$N_cultures[is.na(testing_info_study$N_cultures)] <-  'Not specified'

testing_info_study_site <- merge(isolate_source_site, ast, by = c('source_id', 'location_name'))
testing_info_study_site <-  merge(testing_info_study_site, N_cultures_source_site, by = c('source_id', 'location_name'))
testing_info_study_site <-  merge(testing_info_study_site, guidelines_study_site, by = c('source_id', 'location_name'))
testing_info_study_site <-  merge(testing_info_study_site, control_strain_study_site, by = c('source_id', 'location_name'))
testing_info_study_site <-  merge(testing_info_study_site, ISO_study_site, by = c('source_id', 'location_name'))
testing_info_study_site <-  merge(testing_info_study_site, QA_study_site, by = c('source_id', 'location_name'))
testing_info_study_site$resistance_breakpoints[testing_info_study_site$resistance_breakpoints == 'not specified'] <-  'Not specified'
testing_info_study_site$control_strain[testing_info_study_site$control_strain == 'not specified'] <-  'Not specified'
testing_info_study_site$isolate_source[testing_info_study_site$isolate_source == 'Unknown'] <-  'Not specified'
testing_info_study_site$N_cultures[is.na(testing_info_study_site$N_cultures)] <-  'Not specified'

write.csv(testing_info_study, 'clean_data/testing_info/testing_info_by_study.csv', row.names = F, na = "")
write.csv(testing_info_study_site, 'clean_data/testing_info/testing_info_by_study_site.csv', row.names = F, na = "")

#get the testing info per row
testing_info_per_row <- mydata[,.(row_id, isolate_source, amr_test, resistance_breakpoints, control_strain, ISO_accredited, QA)]
write.csv(testing_info_per_row, 'clean_data/testing_info/testing_info_per_row.csv', row.names = F, na = "")


rm(testing_info_per_row, testing_info_study, testing_info_study_site, 
   QA_study, QA_study_site,
   N_cultures_source, N_cultures_source_site,
   isolate_source_site, isolate_source_study, 
   ISO_study, ISO_study_site, 
   guidelines, guidelines_study, guidelines_study_site,
   control_strain, control_strain_study, control_strain_study_site,
   ast, ast_study)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Select rows for clean data set ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mydata <- mydata[,.(row_id,
                   source_id,
                   study_type,
                   super_region,
                   region,
                   country,
                   iso3,   
                   final_adm0, 
                   final_adm1_name, 
                   final_adm1,
                   final_adm2_name,  
                   final_adm2,
                   final_lat,
                   final_long,
                   best_geo,
                   location_name,
                   geometry_type,
                   setting,
                   start_year,
                   end_year,
                   mid_year,
                   age_range,
                   child_adult,
                   mean_age,
                   percent_male,
                   inpatients_outpatients,
                   demographics_notes,
                   cultures_performed,
                   diagnostic_method,
                   isolate_source,
                   amr_test,
                   resistance_breakpoints,
                   control_strain,
                   ISO_accredited,
                   QA,
                   species,
                   antimicrobial,
                   no_examined,
                   number_resistant,
                   percentage_resistant)] 


colnames(mydata)[colnames(mydata) == 'isolate_source'] <-  'isolate_source_clean'
colnames(mydata)[colnames(mydata) == 'diagnostic_method'] <- 'isolate_source'
colnames(mydata)[colnames(mydata) == 'amr_test'] <- 'AST_method'
colnames(mydata)[colnames(mydata) == 'species'] <- 'serotype'
colnames(mydata)[colnames(mydata) == 'no_examined'] <- 'sample_size'
colnames(mydata) <- gsub('final_', '', colnames(mydata))

saveRDS(mydata, 'clean_data/step2_cleaning.rds')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Get the study information ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mydata$study_type[mydata$source_id == 848] <- 'Prospective & retrospective'
study_info <- unique(mydata, by = c('source_id', 'study_type'))
study_info <- study_info[,.(source_id, study_type)]
citations <- read_excel('clean_data/study_info/citations.xlsx')
study_info <-  merge(study_info, citations, all.x = T)

study_info$publication_year <- as.numeric(unlist(lapply(study_info$author_year, function(x){strsplit(x, '[()]')[[1]][2]})))
write.csv(study_info, 'clean_data/study_info/study_info.csv', row.names = F)

rm(citations)

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# # Aggregate resistance data to one row per study site for data vis ####
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# dt <- mydata[, list(sample_size = sum(sample_size), 
#                     number_resistant = sum(number_resistant),
#                     percentage_resistant = round(sum(number_resistant)/sum(sample_size)*100,2),
#                     start_year = min(start_year),
#                     end_year = max(end_year)),
#              by = c("source_id",
#                     "location_name",
#                     "serotype",
#                     "antimicrobial")]
