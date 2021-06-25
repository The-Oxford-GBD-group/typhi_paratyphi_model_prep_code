rm(list = ls())
library(foreign)
library(data.table)
#~~~~~~~~~~~~~~~~~~~~#
# Prep for ST-GPR ####
#~~~~~~~~~~~~~~~~~~~~#


locs <- read.dbf('Z:/AMR/Shapefiles/GBD2019/GBD2019_analysis_final.dbf')
locs <- locs[c('loc_id', "ihme_lc_id")]

dir.create('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/stgpr_files/')

for(db in c('MDR_typhi', 'MDR_paratyphi', 'FQNS_typhi', 'FQNS_paratyphi')){

  mydata <- readRDS(paste0('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/', db, '.rds'))

  mydata <-  merge(mydata, locs, by.x = c('iso3'), by.y = c('ihme_lc_id'))
  
  #shift the proportions resistance away from 0 and 1
  mydata$percentage_resistant[mydata$percentage_resistant == 0] <- 0.1 
  mydata$percentage_resistant[mydata$percentage_resistant == 100] <- 99.9 
  
  mydata <- mydata[,.(nid = source_id,
                            super_region,
                            region,
                            country = iso3, 
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
  

  write.csv(mydata, paste0('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/stgpr_files/', db, '.csv'), row.names = F)
}