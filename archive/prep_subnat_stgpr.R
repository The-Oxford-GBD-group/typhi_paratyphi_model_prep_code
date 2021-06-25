rm(list = ls())
library(foreign)
library(data.table)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Prep for ST-GPR subnat ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

locs <- read.dbf('Z:/AMR/Shapefiles/GBD2019/GBD2019_analysis_final.dbf')
locs <- locs[c('loc_id', "ihme_lc_id")]
adm1_lookup <- fread("Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/adm1_GBD_lookup2.csv")
KEN_ad1 <- adm1_lookup[iso3=='KEN',.(adm2, location_id)]  
adm1_lookup <- adm1_lookup[iso3!='KEN',.(adm1, location_id)]
dir.create('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/stgpr_files/')

for(db in c('MDR_typhi', 'MDR_paratyphi', 'FQNS_typhi', 'FQNS_paratyphi')){
  mydata <- readRDS(paste0('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/', db, '.rds'))
  
  mydata$adm1 <-  as.numeric(mydata$adm1)
  mydata$adm2 <-  as.numeric(mydata$adm2)
  
  mydata <-  merge(mydata, locs, by.x = c('iso3'), by.y = c('ihme_lc_id'))
  mydata$national_id <- mydata$loc_id
  mydata <-  merge(mydata, adm1_lookup, by = c('adm1'),all.x = T)
  mydata <-  merge(mydata, KEN_ad1, by = c('adm2'),all.x = T)
  
  mydata$loc_id[!is.na(mydata$location_id.y)] <- mydata$location_id.y[!is.na(mydata$location_id.y)]
  mydata$loc_id[!is.na(mydata$location_id.x)] <- mydata$location_id.x[!is.na(mydata$location_id.x)]
  
  mydata$location_id.y <-  NULL
  mydata$location_id.x <-  NULL
  
  #shift the proportions resistance away from 0 and 1
  mydata$percentage_resistant[mydata$percentage_resistant == 0] <- 0.1 
  mydata$percentage_resistant[mydata$percentage_resistant == 100] <- 99.9 
  
  mydata <- mydata[,.(nid = source_id,
                      super_region,
                      region,
                      country = iso3, 
                      national_id,
                      location_id = loc_id,
                      year_id = mid_year, 
                      age_group_id = 22,
                      sex_id = 3,
                      measure_id = 18,
                      is_outlier = 0,
                      sample_size,
                      val = percentage_resistant/100,
                      variance = ((percentage_resistant/100)*(1-(percentage_resistant/100)))/sample_size,
                      QA)]
  
  
  write.csv(mydata, paste0('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/stgpr_files/subnat_', db, '.csv'), row.names = F)
}