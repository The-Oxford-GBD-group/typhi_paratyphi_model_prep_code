#Clean up SEAP data
rm(list = ls())

mydata <- fread("Z:/AMR/Data/Received Data/SEAP 18.May.2021/SEAPData_6May2021.csv")

#clean up the year
mydata$AdmDt <- as.Date(mydata$AdmDt, format = "%d/%m/%Y")
mydata$year <- lubridate::year(mydata$AdmDt)

#create MDR
mydata$MDR <- 0
mydata$MDR[(mydata$ChloramphenicolSen == 'Resistant' | mydata$ChloramphenicolSen == 'Intermediate') &
             (mydata$AmpicillinSen == 'Resistant' | mydata$AmpicillinSen == 'Intermediate')&
             (mydata$CotrimoxazoleSen == 'Resistant' | mydata$AmpicillinSen == 'Intermediate')] <- 1

mydata$MDR[mydata$ChloramphenicolSen == "." | mydata$ChloramphenicolSen == 'Missing' | mydata$ChloramphenicolSen == 'Not performed'|
             mydata$AmpicillinSen == "." | mydata$AmpicillinSen == 'Missing' | mydata$AmpicillinSen == 'Not performed'|
             mydata$CotrimoxazoleSen == "." | mydata$CotrimoxazoleSen == 'Missing' | mydata$CotrimoxazoleSen == 'Not performed'] <- NA

mydata$fq <- NA
mydata$fq[mydata$CiprofloxacinSen == 'Sensitive'] <- 0
mydata$fq[mydata$CiprofloxacinSen == 'Resistant'|mydata$CiprofloxacinSen == 'Intermediate'] <- 1

mydata$ceph <- 0
mydata$ceph[mydata$CeftriaxoneSen == 'Resistant' | mydata$CeftriaxoneSen =='Intermediate' |
              mydata$CefiximeSen== 'Resistant' | mydata$CefiximeSen =='Intermediate'] <- 1

mydata$ceph[(mydata$CeftriaxoneSen == "." | mydata$CeftriaxoneSen == 'Missing' | mydata$CeftriaxoneSen == 'Not performed')&
             (mydata$CefiximeSen == "." | mydata$CefiximeSen == 'Missing' | mydata$CefiximeSen == 'Not performed')] <- NA

#remove missings
mydata <- mydata[!is.na(mydata$year) & mydata$year != 3000,]

# Clean up locations (all from each country are in the same admin 1)
mydata$Country[mydata$Country == 'Bangladesh'] <- 'BGD'
mydata$Country[mydata$Country == 'Nepal'] <- 'NPL'
mydata$Country[mydata$Country == 'Pakistan'] <- 'PAK'

mydata$adj_id[mydata$Country == 'BGD'] <- 823
mydata$adj_id[mydata$Country == 'PAK'] <- 902
mydata$adj_id[mydata$Country == 'NPL'] <- 892

mydata$location_name[mydata$Country == 'BGD'] <- 'Dhaka'
mydata$location_name[mydata$Country == 'PAK'] <- 'Karachi'
mydata$location_name[mydata$Country == 'NPL'] <- 'Kathmandu'

mydata$location_id[mydata$Country == 'BGD'] <- 161
mydata$location_id[mydata$Country == 'PAK'] <- 165
mydata$location_id[mydata$Country == 'NPL'] <- 164

#summarise and make the last datasets
MDR_Typhi <- mydata[Species =='typhi',.(number_resistant = sum(MDR, na.rm = T),
                                 sample_size = length(MDR[!is.na(MDR)])),
                    by = c('Country', 'location_name', 'location_id', 'adj_id', 'year')]

MDR_Paratyphi <- mydata[Species =='paratyphi',.(number_resistant = sum(MDR, na.rm = T),
                                        sample_size = length(MDR[!is.na(MDR)])),
                    by = c('Country', 'location_name','location_id','adj_id','year')]

FQ_Typhi <- mydata[Species =='typhi',.(number_resistant = sum(fq, na.rm = T),
                                        sample_size = length(fq[!is.na(fq)])),
                    by = c('Country', 'location_name', 'location_id', 'adj_id','year')]

FQ_Paratyphi <- mydata[Species =='paratyphi',.(number_resistant = sum(fq, na.rm = T),
                                                sample_size = length(fq[!is.na(fq)])),
                        by = c('Country', 'location_name', 'location_id', 'adj_id','year')]


ceph_Typhi <- mydata[Species =='typhi',.(number_resistant = sum(ceph, na.rm = T),
                                       sample_size = length(ceph[!is.na(ceph)])),
                   by = c('Country', 'location_name', 'location_id', 'adj_id','year')]

ceph_Paratyphi <- mydata[Species =='paratyphi',.(number_resistant = sum(ceph, na.rm = T),
                                               sample_size = length(ceph[!is.na(ceph)])),
                       by = c('Country', 'location_name', 'location_id', 'adj_id','year')]



MDR_Typhi <- MDR_Typhi[sample_size>=5 ,.(nid = 5151,
                                         super_region	= 'South Asia',
                                         region = 'South Asia',
                                         country = Country, 
                                         location_name,
                                         location_id,
                                         subnat_location_id = location_id,
                                         adj_id, adj_id_sSA = NA, adj_id_Asia = adj_id,
                                         year_id = year, 
                                         age_group_id = 22,
                                         sex_id = 3,
                                         measure_id = 18,
                                         number_resistant, sample_size, 
                                         val = number_resistant/sample_size, 
                                         variance = ((number_resistant/sample_size)*(1-(number_resistant/sample_size)))/sample_size, 
                                         QA = 1)]

MDR_Paratyphi <- MDR_Paratyphi[sample_size>=5 ,.(nid = 5151,
                                         super_region	= 'South Asia',
                                         region = 'South Asia',
                                         country = Country, 
                                         location_name,
                                         location_id,
                                         subnat_location_id = location_id,
                                         adj_id, adj_id_sSA = NA, adj_id_Asia = adj_id,
                                         year_id = year, 
                                         age_group_id = 22,
                                         sex_id = 3,
                                         measure_id = 18,
                                         number_resistant, sample_size, 
                                         val = number_resistant/sample_size, 
                                         variance = ((number_resistant/sample_size)*(1-(number_resistant/sample_size)))/sample_size, 
                                         QA = 1)]

FQ_Typhi <- FQ_Typhi[sample_size>=5 ,.(nid = 5151,
                                         super_region	= 'South Asia',
                                         region = 'South Asia',
                                         country = Country, 
                                         location_name,
                                         location_id,
                                         subnat_location_id = location_id,
                                         adj_id, adj_id_sSA = NA, adj_id_Asia = adj_id,
                                         year_id = year, 
                                         age_group_id = 22,
                                         sex_id = 3,
                                         measure_id = 18,
                                         number_resistant, sample_size, 
                                         val = number_resistant/sample_size, 
                                         variance = ((number_resistant/sample_size)*(1-(number_resistant/sample_size)))/sample_size, 
                                         QA = 1)]

FQ_Paratyphi <- FQ_Paratyphi[sample_size>=5 ,.(nid = 5151,
                                                 super_region	= 'South Asia',
                                                 region = 'South Asia',
                                                 country = Country, 
                                                 location_name,
                                                 location_id,
                                                 subnat_location_id = location_id,
                                                 adj_id, adj_id_sSA = NA, adj_id_Asia = adj_id,
                                                 year_id = year, 
                                                 age_group_id = 22,
                                                 sex_id = 3,
                                                 measure_id = 18,
                                                 number_resistant, sample_size, 
                                                 val = number_resistant/sample_size, 
                                                 variance = ((number_resistant/sample_size)*(1-(number_resistant/sample_size)))/sample_size, 
                                                 QA = 1)]

ceph_Typhi <- ceph_Typhi[sample_size>=5 ,.(nid = 5151,
                                       super_region	= 'South Asia',
                                       region = 'South Asia',
                                       country = Country, 
                                       location_name,
                                       location_id,
                                       subnat_location_id = location_id,
                                       adj_id, adj_id_sSA = NA, adj_id_Asia = adj_id,
                                       year_id = year, 
                                       age_group_id = 22,
                                       sex_id = 3,
                                       measure_id = 18,
                                       number_resistant, sample_size, 
                                       val = number_resistant/sample_size, 
                                       variance = ((number_resistant/sample_size)*(1-(number_resistant/sample_size)))/sample_size, 
                                       QA = 1)]

ceph_Paratyphi <- ceph_Paratyphi[sample_size>=5 ,.(nid = 5151,
                                               super_region	= 'South Asia',
                                               region = 'South Asia',
                                               country = Country, 
                                               location_name,
                                               location_id,
                                               subnat_location_id = location_id,
                                               adj_id, adj_id_sSA = NA, adj_id_Asia = adj_id,
                                               year_id = year, 
                                               age_group_id = 22,
                                               sex_id = 3,
                                               measure_id = 18,
                                               number_resistant, sample_size, 
                                               val = number_resistant/sample_size, 
                                               variance = ((number_resistant/sample_size)*(1-(number_resistant/sample_size)))/sample_size, 
                                               QA = 1)]

write.csv(MDR_Typhi, 'Z:/AMR/Pathogens/typhi_paratyphi/Data extraction/datasets/SEAP/MDR_Typhi.csv', row.names = F)
write.csv(MDR_Paratyphi, 'Z:/AMR/Pathogens/typhi_paratyphi/Data extraction/datasets/SEAP/MDR_Paratyphi.csv', row.names = F)
write.csv(FQ_Typhi, 'Z:/AMR/Pathogens/typhi_paratyphi/Data extraction/datasets/SEAP/FQ_Typhi.csv', row.names = F)
write.csv(FQ_Paratyphi, 'Z:/AMR/Pathogens/typhi_paratyphi/Data extraction/datasets/SEAP/FQ_Paratyphi.csv', row.names = F)
write.csv(ceph_Typhi, 'Z:/AMR/Pathogens/typhi_paratyphi/Data extraction/datasets/SEAP/ceph_Typhi.csv', row.names = F)
write.csv(ceph_Paratyphi, 'Z:/AMR/Pathogens/typhi_paratyphi/Data extraction/datasets/SEAP/ceph_Paratyphi.csv', row.names = F)
