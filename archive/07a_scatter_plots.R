library(ggplot2)
library(ggforce)
rm(list = ls())

setwd('K:/AMR/Pathogens/typhi_paratyphi/model_prep/')
mydata <- readRDS('clean_data/clean_data_crosswalked.rds')

typhi <- mydata[mydata$species == 'Typhi',]
MDR_typhi <- typhi[typhi$antimicrobial=='multi-drug' | typhi$antimicrobial=='multi-drug-adjusted',]

pdf('K:/AMR/Pathogens/typhi_paratyphi/model_prep/input_data_plots/MDR_typhi.pdf',
    height = 17,
    width = 12)

for(i in 1:(ceiling((length(unique(MDR_typhi$iso3)))/2))){
  print(ggplot(MDR_typhi, aes(x = mid_year, y = percentage_resistant, label = source_id))+
    geom_point()+
    geom_text(aes(label = source_id), hjust = 0, vjust = 0)+
    facet_wrap_paginate(~iso3, nrow = 2, ncol = 1, page = i)
)}

dev.off()


FQR_typhi <- typhi[grep('FQR', typhi$antimicrobial),]
pdf('K:/AMR/Pathogens/typhi_paratyphi/model_prep/input_data_plots/FQR_typhi.pdf',
    height = 17,
    width = 12)

for(i in 1:(ceiling((length(unique(MDR_typhi$iso3)))/2))){
  print(ggplot(FQR_typhi, aes(x = mid_year, y = percentage_resistant, label = source_id))+
          geom_point()+
          geom_text(aes(label = source_id), hjust = 0, vjust = 0)+
          facet_wrap_paginate(~iso3, nrow = 2, ncol = 1, page = i)
  )}

dev.off()

paratyphi <- mydata[mydata$species == 'Paratyphi',]
MDR_paratyphi <- paratyphi[paratyphi$antimicrobial=='multi-drug' | paratyphi$antimicrobial=='multi-drug-adjusted',]

pdf('K:/AMR/Pathogens/typhi_paratyphi/model_prep/input_data_plots/MDR_paratyphi.pdf',
    height = 17,
    width = 12)

for(i in 1:(ceiling((length(unique(MDR_paratyphi$iso3)))/2))){
  print(ggplot(MDR_paratyphi, aes(x = mid_year, y = percentage_resistant, label = source_id))+
          geom_point()+
          geom_text(aes(label = source_id), hjust = 0, vjust = 0)+
          facet_wrap_paginate(~iso3, nrow = 2, ncol = 1, page = i)
  )}

dev.off()


FQR_paratyphi <- paratyphi[grep('FQR', paratyphi$antimicrobial),]
pdf('K:/AMR/Pathogens/typhi_paratyphi/model_prep/input_data_plots/FQR_paratyphi.pdf',
    height = 17,
    width = 12)

for(i in 1:(ceiling((length(unique(MDR_paratyphi$iso3)))/2))){
  print(ggplot(FQR_paratyphi, aes(x = mid_year, y = percentage_resistant, label = source_id))+
          geom_point()+
          geom_text(aes(label = source_id), hjust = 0, vjust = 0)+
          facet_wrap_paginate(~iso3, nrow = 2, ncol = 1, page = i)
  )}

dev.off()
