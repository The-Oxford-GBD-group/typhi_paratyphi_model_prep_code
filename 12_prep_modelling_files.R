library(data.table)
library(foreign)
rm(list = ls())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Clean up data for use in national STGPR, subnational STGPR, subnational ST-CAR ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
setwd('C:/Users/Annie/Documents/GRAM/typhi_paratyphi')
#read in the location files
locs <- read.dbf('C:/Users/Annie/Documents/GRAM/shapefiles/GBD2019_analysis_final.dbf')
locs <- locs[c('loc_id', "ihme_lc_id")]

adm1_lookup <- fread("model_prep/clean_data/adm1_GBD_lookup2.csv")
KEN_ad1 <- adm1_lookup[iso3=='KEN',.(adm2, location_id)]  
adm1_lookup <- adm1_lookup[iso3!='KEN',.(adm1, location_id)]

adm1_ids <- fread("covariates/all_admin1_typhi_covs.csv")
adm1_ids <- adm1_ids[,.(admin_code, adj_id, adj_id_sSA, adj_id_Asia)]
adm1_ids <- unique(adm1_ids)

#Merge on the locations for each dataset and clean up required variables 
for(db in c('MDR_typhi', 'MDR_paratyphi', 'FQNS_typhi', 'FQNS_paratyphi', 'ceph_typhi', 'ceph_paratyphi')){
  
  mydata <- readRDS(paste0('model_prep/clean_data/', db, '.rds'))
  mydata <- data.table(mydata)
  mydata <- mydata[,.(source_id, super_region, region, country = iso3, location_name, adm1 = as.numeric(adm1), adm2 = as.numeric(adm2), year_id = mid_year, 
                      number_resistant, sample_size, val = round(percentage_resistant/100,3), QA)]
  
  #merge on national location_id
  mydata <-  merge(mydata, locs, by.x = c('country'), by.y = c('ihme_lc_id'))
  mydata$national_id <- mydata$loc_id
  
  #merge on subnational GBD location IDs and create a composit of lowest level GBD locations ID
  mydata <-  merge(mydata, adm1_lookup, by = c('adm1'),all.x = T)
  mydata <-  merge(mydata, KEN_ad1, by = c('adm2'),all.x = T)
  mydata$loc_id[!is.na(mydata$location_id.y)] <- mydata$location_id.y[!is.na(mydata$location_id.y)]
  mydata$loc_id[!is.na(mydata$location_id.x)] <- mydata$location_id.x[!is.na(mydata$location_id.x)]
  mydata$location_id.y <-  NULL
  mydata$location_id.x <-  NULL
  
  #merge on adjacent admin1 IDs for ST-CAR subnational analysis
  mydata <- merge(mydata, adm1_ids, by.x = 'adm1', by.y = 'admin_code', all.x = T)
  
  #jitter 0's and 1's
  mydata$val[mydata$val ==0] <- 0.0001
  mydata$val[mydata$val ==1] <- 0.9999
  
  mydata <- mydata[,.(nid = source_id,
                      super_region,
                      region,
                      country,
                      location_name,
                      location_id = national_id,
                      subnat_location_id = loc_id,
                      adj_id, adj_id_sSA, adj_id_Asia,
                      year_id,
                      age_group_id = 22,
                      sex_id = 3,
                      measure_id = 18,
                      number_resistant,
                      sample_size,
                      val,
                      variance = (val*(1-val))/sample_size,
                      QA)]
  
  #restrict to required (endemic) locations
  mydata <- mydata[mydata$super_region == 'North Africa & Middle East'|
                     mydata$super_region == 'South Asia'|
                     mydata$super_region == 'Southeast Asia, East Asia & Oceania'|
                     mydata$super_region == 'Sub-Saharan Africa',]
  
  if(grepl('paratyphi', db)){
    mydata <- mydata[mydata$super_region == 'South Asia'|
                     mydata$super_region == 'Southeast Asia, East Asia & Oceania',]
  }
  
  #restrict to studies with sample_size >=5
  mydata <- mydata[mydata$sample_size >=5,]
  
  write.csv(mydata, paste0('model_prep/clean_data/', db, '_clean.csv'), row.names = F)
}
