library(ggforce)
library(ggplot2)

rm(list = ls())

clean_data <- readRDS('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/full_cleaned_dataset.rds') 

# testing_info <- read.csv('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/testing_info/testing_info_by_study_site.csv', stringsAsFactors = F)
# demographics <- read.csv('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/demographics/demographics_by_study_site_serotype.csv', stringsAsFactors = F)
# names(demographics)[names(demographics)=='species'] <- 'serotype'
# study_info <- read.csv('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/study_info/study_info.csv', stringsAsFactors = F)


#restrict this to the resistance and location info then add on the aggregated testing and demographics
clean_data <- clean_data[,.(source_id,
                            location_name,
                            super_region,
                            region,
                            country,
                            iso3 ,
                            adm0,
                            adm1,
                            adm2,                  
                            lat,
                            long,
                            best_geo,
                            start_year,
                            end_year,
                            mid_year,
                            serotype,
                            antimicrobial,
                            sample_size,
                            number_resistant,
                            percentage_resistant,
                            QA,
                            isolate_source_clean,
                            resistance_breakpoints,
                            control_strain,
                            AST_method)]

colnames(clean_data)[colnames(clean_data)=='isolate_source_clean'] <- 'isolate_source'

# Set a sample size bins
clean_data$sample_size_bins <- NA
clean_data$sample_size_bins[clean_data$sample_size<50] <- '10-49'
clean_data$sample_size_bins[clean_data$sample_size>=50 &clean_data$sample_size<100 ] <- '50-99'
clean_data$sample_size_bins[clean_data$sample_size>=100 &clean_data$sample_size<500 ] <- '100-499'
clean_data$sample_size_bins[clean_data$sample_size>=500 ] <- '500+'
clean_data$sample_size_bins <-  as.factor(clean_data$sample_size_bins)
clean_data$sample_size_bins <- factor(clean_data$sample_size_bins, levels = c("10-49", "50-99", "100-499", "500+"))
clean_data <- clean_data[!is.na(clean_data$sample_size_bins)]

table(clean_data$AST_method)
clean_data$ast <- NA
clean_data$ast[clean_data$AST_method == 'Agar dilution'] <- 'MIC determination'
clean_data$ast[clean_data$AST_method == 'Agar dilution; VITEK 2'] <- 'MIC determination'
clean_data$ast[clean_data$AST_method == 'E-test'] <- 'MIC determination'
clean_data$ast[clean_data$AST_method == 'E-test; Agar dilution'] <- 'MIC determination'
clean_data$ast[clean_data$AST_method == 'Microdilution'] <- 'MIC determination'
clean_data$ast[clean_data$AST_method == 'Microscan system'] <- 'MIC determination'
clean_data$ast[clean_data$AST_method == 'Phoenix 100 automated system'] <- 'MIC determination'
clean_data$ast[clean_data$AST_method == 'Rapid ATB test'] <- 'MIC determination'
clean_data$ast[clean_data$AST_method == 'VITEK 2'] <- 'MIC determination'

clean_data$ast[clean_data$AST_method == 'Disk-diffusion'] <- 'Disk-diffusion'

clean_data$ast[clean_data$AST_method == 'Disk-diffusion; E-test; VITEK 2'] <- 'Disk-diffusion & MIC determination'
clean_data$ast[clean_data$AST_method == 'Disk-diffusion; E-test'] <- 'Disk-diffusion & MIC determination'
clean_data$ast[clean_data$AST_method == 'Disk-diffusion; Microdilution'] <- 'Disk-diffusion & MIC determination'

clean_data$ast[clean_data$AST_method == 'Unknown'] <- 'Not specified'
table(clean_data$ast)


#restrinct to anitbiotics required
MDR_typhi <- clean_data[(clean_data$antimicrobial == 'multi-drug' | clean_data$antimicrobial == 'MDR - calculated') &
                       clean_data$serotype == 'Typhi',]

MDR_paratyphi <- clean_data[(clean_data$antimicrobial == 'multi-drug' | clean_data$antimicrobial == 'MDR - calculated') &
                           clean_data$serotype == 'Paratyphi',]

FQNS_typhi <- clean_data[(clean_data$antimicrobial == 'FQR - cipro & nalidixic acid' | clean_data$antimicrobial == ' FQR - ciprofloxacin' | clean_data$antimicrobial == 'FQR - nalidixic acid') &
                        clean_data$serotype == 'Typhi',]

FQNS_paratyphi <- clean_data[(clean_data$antimicrobial == 'FQR - cipro & nalidixic acid' | clean_data$antimicrobial == ' FQR - ciprofloxacin'|clean_data$antimicrobial == 'FQR - nalidixic acid') &
                            clean_data$serotype == 'Paratyphi',]

for(st in c('MDR_typhi', 'MDR_paratyphi', 'FQNS_typhi', 'FQNS_paratyphi')){

  mydata <- get(st)
  pdf(paste0('C:/Users/annieb/Desktop/input_data_plots/', st,'/', st, '_by_isolate_source.pdf'),
    height = 17,
    width = 12)

  for(i in 1:(ceiling(length(unique(mydata$super_region))))){
    print(ggplot(mydata[mydata$super_region==unique(mydata$super_region)[i]], aes(x = mid_year, y = percentage_resistant, label = source_id))+
            geom_point(aes(colour = isolate_source, size = sample_size_bins))+
            # geom_text(aes(label = source_id), hjust = 0, vjust = 0)+
            facet_wrap_paginate(~iso3, page = i)
    )}
  
  dev.off()
  
  
  pdf(paste0('C:/Users/annieb/Desktop/input_data_plots/', st,'/', st, '_by_AST.pdf'),
      height = 17,
      width = 12)
  
  for(i in 1:(ceiling(length(unique(mydata$super_region))))){
    print(ggplot(mydata[mydata$super_region==unique(mydata$super_region)[i]], aes(x = mid_year, y = percentage_resistant, label = source_id))+
            geom_point(aes(colour = ast, size = sample_size_bins))+
            # geom_text(aes(label = source_id), hjust = 0, vjust = 0)+
            facet_wrap_paginate(~iso3, page = i)
    )}
  
  dev.off()
  
  mydata$QA <- as.character(mydata$QA)
  pdf(paste0('C:/Users/annieb/Desktop/input_data_plots/', st,'/', st, '_by_QA.pdf'),
      height = 17,
      width = 12)
  
  for(i in 1:(ceiling(length(unique(mydata$super_region))))){
    print(ggplot(mydata[mydata$super_region==unique(mydata$super_region)[i]], aes(x = mid_year, y = percentage_resistant, label = source_id))+
            geom_point(aes(colour = QA, size = sample_size_bins))+
            # geom_text(aes(label = source_id), hjust = 0, vjust = 0)+
            facet_wrap_paginate(~iso3, page = i)
    )}
  
  dev.off()
  
  mydata$control_strain[mydata$control_strain!= 'Not specified'] <- 'Stated'
  pdf(paste0('C:/Users/annieb/Desktop/input_data_plots/', st,'/', st, '_by_control_strain.pdf'),
      height = 17,
      width = 12)
  for(i in 1:(ceiling(length(unique(mydata$super_region))))){
    print(ggplot(mydata[mydata$super_region==unique(mydata$super_region)[i]], aes(x = mid_year, y = percentage_resistant, label = source_id))+
            geom_point(aes(colour = control_strain, size = sample_size_bins))+
            # geom_text(aes(label = source_id), hjust = 0, vjust = 0)+
            facet_wrap_paginate(~iso3, page = i)
    )}
  
  dev.off()
  
  mydata$resistance_breakpoints[mydata$resistance_breakpoints == 'Not specified'] <- 'No'
  mydata$resistance_breakpoints[mydata$resistance_breakpoints != 'No'] <- 'Yes'
  pdf(paste0('C:/Users/annieb/Desktop/input_data_plots/', st,'/', st, '_by_breakpoints.pdf'),
      height = 17,
      width = 12)
  for(i in 1:(ceiling(length(unique(mydata$super_region))))){
    print(ggplot(mydata[mydata$super_region==unique(mydata$super_region)[i]], aes(x = mid_year, y = percentage_resistant, label = source_id))+
            geom_point(aes(colour = resistance_breakpoints, size = sample_size_bins))+
            # geom_text(aes(label = source_id), hjust = 0, vjust = 0)+
            facet_wrap_paginate(~iso3, page = i)
    )}
  
  dev.off()
}
