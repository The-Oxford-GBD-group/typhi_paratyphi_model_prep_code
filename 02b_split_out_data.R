#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Seperate out MDR FQNS Typhi and Paratyphi ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rm(list = ls())

mydata <- readRDS('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/full_cleaned_dataset.rds')

table(mydata$antimicrobial)
MDR_typhi <- mydata[(mydata$antimicrobial == 'multi-drug' | mydata$antimicrobial == 'MDR - calculated') &
                      mydata$serotype == 'Typhi',]


MDR_paratyphi <- mydata[(mydata$antimicrobial == 'multi-drug' | mydata$antimicrobial == 'MDR - calculated') &
                      mydata$serotype == 'Paratyphi',]

FQNS_typhi <- mydata[(mydata$antimicrobial == 'FQR - cipro & nalidixic acid' | mydata$antimicrobial == 'FQR - ciprofloxacin' | mydata$antimicrobial == 'FQR - nalidixic acid') &
                      mydata$serotype == 'Typhi',]


FQNS_paratyphi <- mydata[(mydata$antimicrobial == 'FQR - cipro & nalidixic acid' | mydata$antimicrobial == 'FQR - ciprofloxacin'|mydata$antimicrobial == 'FQR - nalidixic acid') &
                          mydata$serotype == 'Paratyphi',]

381+420+3
#check what locations are in each
table(MDR_typhi$best_geo)
table(FQNS_typhi$best_geo)
table(MDR_paratyphi$best_geo)
table(FQNS_paratyphi$best_geo)

#save files
saveRDS(MDR_typhi, 'Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/MDR_typhi.rds')
saveRDS(MDR_paratyphi, 'Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/MDR_paratyphi.rds')
saveRDS(FQNS_typhi, 'Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/FQNS_typhi.rds')
saveRDS(FQNS_paratyphi, 'Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/FQNS_paratyphi.rds')
