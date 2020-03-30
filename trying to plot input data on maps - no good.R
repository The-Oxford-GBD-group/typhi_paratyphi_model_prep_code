library(data.table)
library(ggplot2)
library(viridis)
library(sf)
source("H:/Functions/round2.R")

adm0_master <- st_read('Z:/AMR/Shapefiles/admin2013_0.shp')
adm1_master <- st_read('Z:/AMR/Shapefiles/admin2013_1.shp')
adm2_master <- st_read('C:/Users/annieb/Desktop/admin2013_2.shp')

adm0_master <- st_simplify(adm0_master, dTolerance = 0.1, preserveTopology = T)
adm1_master <- st_simplify(adm1_master, dTolerance = 0.1, preserveTopology = T)
adm2_master <- st_simplify(adm2_master, dTolerance = 0.1, preserveTopology = T)

clean_data <- readRDS('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/full_cleaned_dataset.rds') 

testing_info <- read.csv('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/testing_info/testing_info_by_study_site.csv', stringsAsFactors = F)
demographics <- read.csv('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/demographics/demographics_by_study_site_serotype.csv', stringsAsFactors = F)
names(demographics)[names(demographics)=='species'] <- 'serotype'
study_info <- read.csv('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/study_info/study_info.csv', stringsAsFactors = F)

clean_data <- data.table(spare)

#restrict this to the resistance and location info then add on the aggregated testing and demographics
clean_data <- clean_data[,.(source_id,
                            location_name,
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
                            serotype,
                            antimicrobial,
                            sample_size,
                            number_resistant)]

clean_data <- merge(clean_data, demographics, by = c('source_id', 'location_name', 'serotype'), allow.cartesian = TRUE)
clean_data <- merge(clean_data, testing_info, by = c('source_id', 'location_name'), allow.cartesian = TRUE)
clean_data <- merge(clean_data, study_info, by = c('source_id'), allow.cartesian = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Aggregate resistance data to one row per study site for data vis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ag_data <- clean_data[, list(sample_size = sum(sample_size),
                              number_resistant = sum(number_resistant),
                              percentage_resistant = round(sum(number_resistant)/sum(sample_size)*100,2),
                              start_year = min(start_year),
                              end_year = max(end_year)),
             by = c("source_id",
                    "location_name",
                    "country",
                    "iso3" ,
                    "adm0",
                    "adm1",
                    "adm2",                  
                    "lat",
                    "long",
                    "best_geo",
                    "child_adult",
                    "isolate_source",
                    "ast",
                    "N_cultures",
                    "resistance_breakpoints",
                    "control_strain",        
                    "ISO_accredited",
                    "QA",
                    'serotype',
                    'antimicrobial',
                    'publication_year')]


#### calcaulte the mid year (need the publication year) and the time period
ag_data$start_year_temp <- ag_data$start_year
ag_data$end_year_temp <- ag_data$end_year

ag_data$start_year_temp[is.na(ag_data$start_year_temp)] <- ag_data$publication_year[is.na(ag_data$start_year_temp)]-median(ag_data$publication_year-ag_data$start_year, na.rm = T) 
ag_data$end_year_temp[is.na(ag_data$end_year_temp)] <- ag_data$publication_year[is.na(ag_data$end_year_temp)]-median(ag_data$publication_year-ag_data$end_year, na.rm = T) 

#mid year (round up if 0.5)
ag_data$mid_year <- round2(((ag_data$start_year_temp+ag_data$end_year_temp)/2),0)

ag_data$start_year_temp <- NULL
ag_data$end_year_temp <- NULL

#time periods
ag_data$time_period <- NA
ag_data$time_period[ag_data$mid_year>=1990 & ag_data$mid_year<1995] <-  "1990-1994"
ag_data$time_period[ag_data$mid_year>=1995 & ag_data$mid_year<2000] <-  "1995-1999"
ag_data$time_period[ag_data$mid_year>=2000 & ag_data$mid_year<2005] <-  "2000-2004"
ag_data$time_period[ag_data$mid_year>=2005 & ag_data$mid_year<2010] <-  "2005-2009"
ag_data$time_period[ag_data$mid_year>=2010 & ag_data$mid_year<2015] <-  "2010-2014"
ag_data$time_period[ag_data$mid_year>=2015] <-  "2015-2018"

# Set a sample size bins
ag_data$sample_size_bins <- NA
ag_data$sample_size_bins[ag_data$sample_size<50] <- '10-49'
ag_data$sample_size_bins[ag_data$sample_size>=50 &ag_data$sample_size<100 ] <- '50-99'
ag_data$sample_size_bins[ag_data$sample_size>=100 &ag_data$sample_size<500 ] <- '100-499'
ag_data$sample_size_bins[ag_data$sample_size>=500 ] <- '500+'
ag_data$sample_size_bins <-  as.factor(ag_data$sample_size_bins)
ag_data$sample_size_bins <- factor(ag_data$sample_size_bins, levels = c("10-49", "50-99", "100-499", "500+"))
ag_data <- ag_data[!is.na(ag_data$sample_size_bins)]

#restrinct to anitbiotics required
MDR_typhi <- ag_data[(ag_data$antimicrobial == 'multi-drug' | ag_data$antimicrobial == 'MDR - calculated') &
                          ag_data$serotype == 'Typhi',]


MDR_paratyphi <- ag_data[(ag_data$antimicrobial == 'multi-drug' | ag_data$antimicrobial == 'MDR - calculated') &
                              ag_data$serotype == 'Paratyphi',]

FQNS_typhi <- ag_data[(ag_data$antimicrobial == 'FQR - cipro & nalidixic acid' | ag_data$antimicrobial == ' FQR - ciprofloxacin' | ag_data$antimicrobial == 'FQR - nalidixic acid') &
                           ag_data$serotype == 'Typhi',]


FQNS_paratyphi <- ag_data[(ag_data$antimicrobial == 'FQR - cipro & nalidixic acid' | ag_data$antimicrobial == ' FQR - ciprofloxacin'|ag_data$antimicrobial == 'FQR - nalidixic acid') &
                               ag_data$serotype == 'Paratyphi',]


# Get the background
regions.lookup <- read.csv("Z:/AMR/Misc/GBD_regions//GBD regions lookup.csv")

paratyphi_endemic <-  regions.lookup[regions.lookup$region == 'East Asia'|
                                      regions.lookup$region == 'North Africa & Middle East'|
                                      regions.lookup$region == 'South Asia'|
                                      regions.lookup$region == 'Southeast Asia'|
                                      regions.lookup$region == 'Oceania',]

typhi_endemic <-  regions.lookup[regions.lookup$region == 'East Asia'|
                                      regions.lookup$region == 'North Africa & Middle East'|
                                      regions.lookup$region == 'South Asia'|
                                      regions.lookup$region == 'Southeast Asia'|
                                      regions.lookup$region =='Central Sub-Saharan Africa'|
                                      regions.lookup$region =='Western Sub-Saharan Africa'|
                                      regions.lookup$region =='Eastern Sub-Saharan Africa'|
                                      regions.lookup$region =='Southern Sub-Saharan Africa'|
                                      regions.lookup$region == 'Oceania',]


non_endemic_typhi <- adm0_master[!(adm0_master$GAUL_CODE%in%typhi_endemic$GAUL_CODE),]
non_endemic_paratyphi <- adm0_master[!(adm0_master$GAUL_CODE%in%paratyphi_endemic$GAUL_CODE),]


#get the admin files to plot
mdr_typhi_adm0 <- merge(adm0_master, MDR_typhi, by.x = c('GAUL_CODE'), by.y = c('adm0'), all.y = T) 
mdr_typhi_adm0 <- mdr_typhi_adm0[mdr_typhi_adm0$best_geo == 'national',] 

mdr_typhi_adm1 <- merge(adm1_master, MDR_typhi, by.x = c('GAUL_CODE'), by.y = c('adm1'), all.y = T) 
mdr_typhi_adm1 <- mdr_typhi_adm1[mdr_typhi_adm1$best_geo == 'admin1',] 

mdr_typhi_adm2 <- merge(adm2_master, MDR_typhi, by.x = c('GAUL_CODE'), by.y = c('adm2'), all.y = T) 
mdr_typhi_adm2 <- mdr_typhi_adm2[mdr_typhi_adm2$best_geo == 'admin2',] 

mdr_paratyphi_adm0 <- merge(adm0_master, MDR_paratyphi, by.x = c('GAUL_CODE'), by.y = c('adm0'), all.y = T) 
mdr_paratyphi_adm0 <- mdr_paratyphi_adm0[mdr_paratyphi_adm0$best_geo == 'national',] 

mdr_paratyphi_adm1 <- merge(adm1_master, MDR_paratyphi, by.x = c('GAUL_CODE'), by.y = c('adm1'), all.y = T) 
mdr_paratyphi_adm1 <- mdr_paratyphi_adm1[mdr_paratyphi_adm1$best_geo == 'admin1',] 

mdr_paratyphi_adm2 <- merge(adm2_master, MDR_paratyphi, by.x = c('GAUL_CODE'), by.y = c('adm2'), all.y = T) 
mdr_paratyphi_adm2 <- mdr_paratyphi_adm2[mdr_paratyphi_adm2$best_geo == 'admin2',] 

fqns_typhi_adm0 <- merge(adm0_master, FQNS_typhi, by.x = c('GAUL_CODE'), by.y = c('adm0'), all.y = T) 
fqns_typhi_adm0 <- fqns_typhi_adm0[fqns_typhi_adm0$best_geo == 'national',] 

fqns_typhi_adm1 <- merge(adm1_master, FQNS_typhi, by.x = c('GAUL_CODE'), by.y = c('adm1'), all.y = T) 
fqns_typhi_adm1 <- fqns_typhi_adm1[fqns_typhi_adm1$best_geo == 'admin1',] 

fqns_typhi_adm2 <- merge(adm2_master, FQNS_typhi, by.x = c('GAUL_CODE'), by.y = c('adm2'), all.y = T) 
fqns_typhi_adm2 <- fqns_typhi_adm2[fqns_typhi_adm2$best_geo == 'admin2',] 

fqns_paratyphi_adm0 <- merge(adm0_master, MDR_paratyphi, by.x = c('GAUL_CODE'), by.y = c('adm0'), all.y = T) 
fqns_paratyphi_adm0 <- fqns_paratyphi_adm0[fqns_paratyphi_adm0$best_geo == 'national',] 

fqns_paratyphi_adm1 <- merge(adm1_master, FQNS_paratyphi, by.x = c('GAUL_CODE'), by.y = c('adm1'), all.y = T) 
fqns_paratyphi_adm1 <- fqns_paratyphi_adm1[fqns_paratyphi_adm1$best_geo == 'admin1',] 

fqns_paratyphi_adm2 <- merge(adm2_master, FQNS_paratyphi, by.x = c('GAUL_CODE'), by.y = c('adm2'), all.y = T) 
fqns_paratyphi_adm2 <- fqns_paratyphi_adm2[fqns_paratyphi_adm2$best_geo == 'admin2',] 

png('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/input_data_plots/MDR_typhi.png',
    height = 30, width = 30, unit = 'cm', res = 300)
ggplot()+
  geom_sf(data = adm0_master, fill = '#bdbdbd',colour = 'black', size = 0.25)+
  # geom_sf(data = mdr_typhi_adm0, aes(fill = percentage_resistant) ,colour = 'black', size = 0.25)+
  # geom_sf(data = mdr_typhi_adm1, aes(fill = percentage_resistant) ,colour = 'black', size = 0.25)+
  # geom_sf(data = mdr_typhi_adm2, aes(fill = percentage_resistant), colour = 'black', size = 0.25)+
  geom_point(data = MDR_typhi, aes(x = long, y = lat, fill = percentage_resistant, size = sample_size_bins), colour = 'black', pch = 21)+
  geom_sf(data = non_endemic_typhi, fill = '#636363',colour = 'black', size = 0.25)+
  scale_fill_viridis(option='inferno', discrete = F, direction = -1)+
  scale_size_discrete(range = seq(3, 6))+
  facet_wrap(~time_period, ncol = 2)+
  theme_bw()+
  theme(line = element_blank(),
        axis.text = element_blank())+
  xlab("") + 
  ylab("")+
  labs(size = 'Sample size', fill = 'Resistance (%)', colour = 'Resistance (%)')+
  xlim(0, 180)

dev.off()



png('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/input_data_plots/MDR_paratyphi.png',
    height = 20, width = 20, unit = 'cm', res = 300)
ggplot()+
  geom_sf(data = adm0_master, fill = '#bdbdbd',colour = 'black', size = 0.25)+
  # geom_sf(data = mdr_paratyphi_adm0, aes(fill = percentage_resistant) ,colour = 'black', size = 0.25)+
  # geom_sf(data = mdr_paratyphi_adm1, aes(fill = percentage_resistant) ,colour = 'black', size = 0.25)+
  # geom_sf(data = mdr_paratyphi_adm2, aes(fill = percentage_resistant), colour = 'black', size = 0.25)+
  geom_point(data = MDR_paratyphi, aes(x = long, y = lat, fill = percentage_resistant, size = sample_size_bins), colour = 'black', pch = 21)+
  geom_sf(data = non_endemic_paratyphi, fill = '#636363',colour = 'black', size = 0.25)+
  scale_fill_viridis(option='inferno', discrete = F, direction = -1)+
  scale_size_discrete(range = seq(1, 4))+
  facet_wrap(~time_period, ncol = 2)+
  theme_bw()+
  theme(line = element_blank(),
        axis.text = element_blank())+
  xlab("") + 
  ylab("")+
  labs(size = 'Sample size', fill = 'Resistance (%)', colour = 'Resistance (%)')

dev.off()


png('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/input_data_plots/FQNS_typhi.png',
    height = 30, width = 30, unit = 'cm', res = 300)
ggplot()+
  geom_sf(data = adm0_master, fill = '#bdbdbd',colour = 'black', size = 0.25)+
  geom_sf(data = fqns_typhi_adm0, aes(fill = percentage_resistant) ,colour = 'black', size = 0.25)+
  geom_sf(data = fqns_typhi_adm1, aes(fill = percentage_resistant) ,colour = 'black', size = 0.25)+
  geom_sf(data = fqns_typhi_adm2, aes(fill = percentage_resistant), colour = 'black', size = 0.25)+
  geom_point(data = FQNS_typhi, aes(x = long, y = lat, fill = percentage_resistant, size = sample_size_bins), colour = 'black', pch = 21)+
  geom_sf(data = non_endemic_typhi, fill = '#636363',colour = 'black', size = 0.25)+
  scale_fill_viridis(option='inferno', discrete = F, direction = -1)+
  scale_size_discrete(range = seq(1, 4))+
  facet_wrap(~time_period, ncol = 2)+
  theme_bw()+
  theme(line = element_blank(),
        axis.text = element_blank())+
  xlab("") + 
  ylab("")+
  labs(size = 'Sample size', fill = 'Resistance (%)', colour = 'Resistance (%)')

dev.off()



png('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/input_data_plots/FQNS_paratyphi.png',
    height = 20, width = 20, unit = 'cm', res = 300)
ggplot()+
  geom_sf(data = adm0_master, fill = '#bdbdbd',colour = 'black', size = 0.25)+
  geom_sf(data = fqns_paratyphi_adm0, aes(fill = percentage_resistant) ,colour = 'black', size = 0.25)+
  geom_sf(data = fqns_paratyphi_adm1, aes(fill = percentage_resistant) ,colour = 'black', size = 0.25)+
  geom_sf(data = fqns_paratyphi_adm2, aes(fill = percentage_resistant), colour = 'black', size = 0.25)+
  geom_point(data = FQNS_paratyphi, aes(x = long, y = lat, fill = percentage_resistant, size = sample_size_bins), colour = 'black', pch = 21)+
  geom_sf(data = non_endemic_paratyphi, fill = '#636363',colour = 'black', size = 0.25)+
  scale_fill_viridis(option='inferno', discrete = F, direction = -1)+
  scale_size_discrete(range = seq(1, 4))+
  facet_wrap(~time_period, ncol = 2)+
  theme_bw()+
  theme(line = element_blank(),
        axis.text = element_blank())+
  xlab("") + 
  ylab("")+
  labs(size = 'Sample size', fill = 'Resistance (%)')

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
model1 <- lm(percentage_resistant~start_year, data = MDR_typhi)
summary(model1)

model2 <- lm(percentage_resistant~start_year+isolate_source, data = MDR_typhi)
summary(model2)

model3 <- lm(percentage_resistant~start_year+QA, data = MDR_typhi)
summary(model3)

#build a lasso regression to see which have any significant effect
MDR_typhi$sample_size_bins <- NA
MDR_typhi$sample_size_bins[MDR_typhi$sample_size<50] <- '10-49'
MDR_typhi$sample_size_bins[MDR_typhi$sample_size>=50 &MDR_typhi$sample_size<100 ] <- '50-99'
MDR_typhi$sample_size_bins[MDR_typhi$sample_size>=100 &MDR_typhi$sample_size<500 ] <- '100-499'
MDR_typhi$sample_size_bins[MDR_typhi$sample_size>=500 ] <- '500+'
MDR_typhi$sample_size_bins <-  as.factor(MDR_typhi$sample_size_bins)
MDR_typhi$sample_size_bins <- factor(MDR_typhi$sample_size_bins, levels = c("10-49", "50-99", "100-499", "500+"))
MDR_typhi <- MDR_typhi[!is.na(MDR_typhi$sample_size_bins)]

MDR_typhi$time_period <- NA
MDR_typhi$time_period[MDR_typhi$mid_year>=1990 & MDR_typhi$mid_year<1995] <-  "1990-1994"
MDR_typhi$time_period[MDR_typhi$mid_year>=1995 & MDR_typhi$mid_year<2000] <-  "1995-1999"
MDR_typhi$time_period[MDR_typhi$mid_year>=2000 & MDR_typhi$mid_year<2005] <-  "2000-2004"
MDR_typhi$time_period[MDR_typhi$mid_year>=2005 & MDR_typhi$mid_year<2010] <-  "2005-2009"
MDR_typhi$time_period[MDR_typhi$mid_year>=2010 & MDR_typhi$mid_year<2015] <-  "2010-2014"
MDR_typhi$time_period[MDR_typhi$mid_year>=2015] <-  "2015-2018"


ggplot()+
  geom_sf(data = IND, fill = '#bdbdbd',colour = 'black', size = 0.25)+
  geom_point(data = MDR_typhi, aes(x = long, y = lat, fill = percentage_resistant, size = sample_size_bins), colour = 'black', pch = 21)+
  scale_fill_viridis(option='inferno', discrete = F, direction = -1)+
  scale_size_discrete(range = seq(1, 4))+
  facet_wrap(~time_period, ncol = 2)+
  theme_bw()+
  theme(line = element_blank(),
        axis.text = element_blank())+
  xlab("") + 
  ylab("")+
  labs(size = 'Sample size', fill = 'Resistance (%)', colour = 'Resistance (%)')

