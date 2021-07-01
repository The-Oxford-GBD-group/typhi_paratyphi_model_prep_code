#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Seperate out MDR FQNS Typhi and Paratyphi ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rm(list = ls())
setwd('C:/Users/Annie/Documents/GRAM/typhi_paratyphi/model_prep/')
mydata <- readRDS('clean_data/full_cleaned_dataset.rds')
mydata <- mydata[mydata$sample_size>=5]
table(mydata$antimicrobial)
MDR_typhi <- mydata[(mydata$antimicrobial == 'multi-drug' | mydata$antimicrobial == 'MDR - calculated') &
                      mydata$serotype == 'Typhi',]


MDR_paratyphi <- mydata[(mydata$antimicrobial == 'multi-drug' | mydata$antimicrobial == 'MDR - calculated') &
                      mydata$serotype == 'Paratyphi',]

FQNS_typhi <- mydata[(mydata$antimicrobial == 'FQR - cipro & nalidixic acid' | mydata$antimicrobial == 'FQR - ciprofloxacin' | mydata$antimicrobial == 'FQR - nalidixic acid') &
                      mydata$serotype == 'Typhi',]


FQNS_paratyphi <- mydata[(mydata$antimicrobial == 'FQR - cipro & nalidixic acid' | mydata$antimicrobial == 'FQR - ciprofloxacin'|mydata$antimicrobial == 'FQR - nalidixic acid') &
                          mydata$serotype == 'Paratyphi',]

#For 3GC take ceftriaxone if available, if not take cefixime
ceft_typhi <- mydata[mydata$antimicrobial == 'ceftriaxone' & mydata$serotype == 'Typhi',]
cefi_typhi <- mydata[mydata$antimicrobial == 'cefixime' & mydata$serotype == 'Typhi',]
cefi_typhi <- cefi_typhi[!(cefi_typhi$source_id %in% ceft_typhi$source_id),]
ceph_typhi <- rbind(ceft_typhi, cefi_typhi)
rm(ceft_typhi, cefi_typhi)

ceft_paratyphi <- mydata[mydata$antimicrobial == 'ceftriaxone' & mydata$serotype == 'Paratyphi',]
cefi_paratyphi <- mydata[mydata$antimicrobial == 'cefixime' & mydata$serotype == 'Paratyphi',]
cefi_paratyphi <- cefi_paratyphi[!(cefi_paratyphi$source_id %in% ceft_paratyphi$source_id),]
ceph_paratyphi <- rbind(ceft_paratyphi, cefi_paratyphi)
rm(ceft_paratyphi, cefi_paratyphi)

#check what locations are in each
table(MDR_typhi$best_geo)
table(FQNS_typhi$best_geo)
table(MDR_paratyphi$best_geo)
table(FQNS_paratyphi$best_geo)

#save files
saveRDS(MDR_typhi, 'clean_data/MDR_typhi.rds')
saveRDS(MDR_paratyphi, 'clean_data/MDR_paratyphi.rds')
saveRDS(FQNS_typhi, 'clean_data/FQNS_typhi.rds')
saveRDS(FQNS_paratyphi, 'clean_data/FQNS_paratyphi.rds')
saveRDS(ceph_typhi, 'clean_data/ceph_typhi.rds')
saveRDS(ceph_paratyphi, 'clean_data/ceph_paratyphi.rds')

#total number of studies
length(unique(c(MDR_typhi$source_id, MDR_paratyphi$source_id, FQNS_paratyphi$source_id, FQNS_typhi$source_id)))
