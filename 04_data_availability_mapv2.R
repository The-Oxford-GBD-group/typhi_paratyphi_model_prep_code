rm(list =ls())
library(sf)
library(ggplot2)
library(data.table)
library(gridExtra)
library(plyr)
library(cowplot)

# Set he directory to where youre input data is and where you want to save the map
setwd('Z:/AMR/Pathogens/data_availability_maps')

#read in your dataset
MDR_Typhi <- fread("Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/stgpr_files/MDR_typhi.csv")
FQNS_Typhi <- fread("Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/stgpr_files/FQNS_typhi.csv")
MDR_Paratyphi <- fread("Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/stgpr_files/MDR_paratyphi.csv")
FQNS_Paratyphi <- fread("Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/stgpr_files/FQNS_paratyphi.csv")
MDR_Typhi$indicator <- 'MDR_Typhi'
FQNS_Typhi$indicator <- 'FQNS_Typhi'
MDR_Paratyphi$indicator <- 'MDR_Paratyphi'
FQNS_Paratyphi$indicator <- 'FQNS_Paratyphi'

master.data <- rbind(MDR_Typhi, FQNS_Typhi, MDR_Paratyphi, FQNS_Paratyphi)
rm(MDR_Typhi, FQNS_Typhi, MDR_Paratyphi, FQNS_Paratyphi)
master.data <-  master.data[master.data$sample_size>=5,]

#remove non endemic countries
master.data <-  master.data[master.data$super_region != 'Central Europe, Eastern Europe & Central Asia',]
master.data <-  master.data[master.data$super_region != 'Latin America & Caribbean',]
master.data <-  master.data[master.data$super_region != 'High Income',]
master.data <-  master.data[master.data$country != 'ZAF',]


master.data <- master.data[!((master.data$indicator == 'FQNS_Paratyphi' | 
                         master.data$indicator == 'MDR_Paratyphi') & 
                        (master.data$super_region =='Sub-Saharan Africa'|
                           master.data$super_region =='North Africa & Middle East')),]


#get a summary of the amount of data
table(master.data$indicator)
totals <- master.data[,.(studies = .N,
                         isolates = sum(sample_size)),
                      by = c('indicator')]

length(unique(master.data$nid))
length(unique(master.data$country))
length(unique(master.data$year_id))

master.data$serotype <- NA
master.data$serotype[grepl('Typhi', master.data$indicator)] <- 'Typhi'
master.data$serotype[grepl('Paratyphi', master.data$indicator)] <- 'Paratyphi'
table(master.data$serotype)

#get approx numbers of isolates
isolates <- master.data[,.(isolates = max(sample_size)),
                        by = c('nid', 'country', 'year_id', 'serotype')]

sum(isolates$isolates)
sum(isolates$isolates[isolates$serotype == 'Typhi'])
sum(isolates$isolates[isolates$serotype == 'Paratyphi'])

#read in the shapefile and simplify
my_shp <- st_read('Z:/AMR/Shapefiles/GBD2020/GBD2020_mapping_final.shp')
my_shp <- st_simplify(my_shp, dTolerance = 0.1, preserveTopology = T)
my_shp <- my_shp[my_shp$level == 3,]

# For all years first 
#get the number of datapoints per country and group into something logical
mydata <- master.data[,.(data_points = .N),
                 by = c('location_id', 'indicator')]
hist(mydata$data_points)

mydata$data_points_group <- NA
mydata$data_points_group[mydata$data_points <=5] <- "1-5"
mydata$data_points_group[mydata$data_points>=6 & mydata$data_points <=10] <- "6-10"
mydata$data_points_group[mydata$data_points>=11 & mydata$data_points <=50] <- "11-50"
mydata$data_points_group[mydata$data_points>=51 & mydata$data_points <=100] <- "51-100"
mydata$data_points_group[mydata$data_points>100] <- ">100"
mydata$data_points_group <-  factor(mydata$data_points_group, levels = c('Non-endemic', 'No data', "1-5", "6-10", "11-50", "51-100", ">100"))

#create a shapefile with just the input countries
input_shp <-  merge(my_shp, mydata, by.x = 'loc_id', by.y = 'location_id')

# get the endemic countries
endemic_locs <- read.csv('Z:/AMR/Pathogens/typhi_paratyphi/misc/typhi_endemic_locations2.csv', stringsAsFactors = F)
endemic_locs  <- endemic_locs[c("loc_id", "loc_name", "spr_reg_id", "region_id", 'level', 'ihme_lc_id')]
endemic_locs <- endemic_locs[endemic_locs$region_id != 65,]
endemic_locs <- my_shp[my_shp$loc_id %in% endemic_locs$loc_id,]

paratyphi_endemic_locs <- endemic_locs[endemic_locs$spr_reg_id != 166,]
paratyphi_endemic_locs <- endemic_locs[endemic_locs$spr_reg_id != 137,]

# Western Sahara and French Guiana arent in the shapefile so get them from another
background <- st_read('Z:/AMR/Shapefiles/admin2013_0.shp')
background <- background[background$name == 'Western Sahara' | background$name == 'French Guiana',]
background <- st_simplify(background, dTolerance = 0.1, preserveTopology = T)

#plot out each of the maps
MDR_Typhi <- 
  ggplot()+
  geom_sf(data = my_shp, fill = '#737373',colour = 'black', size = 0.25)+
  geom_sf(data = background, fill = '#737373',colour = 'black', size = 0.25)+
  geom_sf(data = endemic_locs, fill = '#bdbdbd',colour = 'black', size = 0.25)+
  geom_sf(data = input_shp[input_shp$indicator == 'MDR_Typhi',], aes(fill = data_points_group),colour = 'black', size = 0.25)+
  scale_fill_manual(values = c('#ffffcc', '#c2e699', '#78c679', '#31a354', '#006837'),
                    labels = c("1-5","6-10", "11-50", "51-100", ">100"))+
  theme_bw()+
  theme(line = element_blank(),
        axis.text = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank())+
  xlim(-20,150)+
  ylim(-35,55)+
  guides(fill=FALSE)


MDR_Paratyphi <- 
  ggplot()+
  geom_sf(data = my_shp, fill = '#737373',colour = 'black', size = 0.25)+
  geom_sf(data = background, fill = '#737373',colour = 'black', size = 0.25)+
  geom_sf(data = paratyphi_endemic_locs, fill = '#bdbdbd',colour = 'black', size = 0.25)+
  geom_sf(data = input_shp[input_shp$indicator == 'MDR_Paratyphi' & input_shp$ihme_lc_id!='IRN',], aes(fill = data_points_group),colour = 'black', size = 0.25)+
  scale_fill_manual(values = c('#ffffcc', '#c2e699', '#78c679', '#31a354', '#006837'),
                    labels = c("1-5","6-10", "11-50", "51-100", ">100"))+
  theme_bw()+
  theme(line = element_blank(),
        axis.text = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank())+
  xlim(60, 150)+
  ylim(-12,55)+
  guides(fill=FALSE)

FQNS_Paratyphi <- 
  ggplot()+
  geom_sf(data = my_shp, fill = '#737373',colour = 'black', size = 0.25)+
  geom_sf(data = background, fill = '#737373',colour = 'black', size = 0.25)+
  geom_sf(data = paratyphi_endemic_locs, fill = '#bdbdbd',colour = 'black', size = 0.25)+
  geom_sf(data = input_shp[input_shp$indicator == 'FQNS_Paratyphi'& input_shp$ihme_lc_id!='IRN',], aes(fill = data_points_group),colour = 'black', size = 0.25)+
  scale_fill_manual(values = c('#ffffcc', '#c2e699', '#78c679', '#31a354', '#006837'),
                    labels = c("1-5","6-10", "11-50", "51-100", ">100"))+
  theme_bw()+
  theme(line = element_blank(),
        axis.text = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank())+
  xlim(60, 150)+
  ylim(-12,55)+
  guides(fill=FALSE)

#merge shapefiles together to get legend
FQ <- input_shp[input_shp$indicator == 'FQNS_Typhi',]
FQ$data_points <-  NULL
FQ$indicator <-  NULL

FQ_non_endemic <- my_shp[!(my_shp$loc_id %in% FQ$loc_id),]
FQ_non_endemic$data_points_group <- 'Non-endemic'

FQ_background <- endemic_locs[!(endemic_locs$loc_id %in% FQ$loc_id),]
FQ_background$data_points_group <- 'No data'

FQNS_shp <- rbind(FQ, FQ_non_endemic, FQ_background)

FQNS_Typhi <-
  ggplot()+
  geom_sf(data = background, fill = '#737373', colour = 'black', size = 0.25)+
  geom_sf(data = FQNS_shp, aes(fill = data_points_group),colour = 'black', size = 0.25)+
  scale_fill_manual(values = c('#737373', '#bdbdbd', '#ffffcc', '#c2e699', '#78c679', '#31a354', '#006837'),
                    labels = c('Non-endemic', 'No data', "1-5","6-10", "11-50", "51-100", ">100"))+
  theme_bw()+
  theme(line = element_blank(),
        axis.text = element_blank(),
        legend.position=c(.1,.22),
        legend.key.size = unit(0.2, "cm"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank())+
  xlim(-20,150)+
  ylim(-35,55)+
  labs(fill = 'Data points')

pdf('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/data_availability/data_availabilty_maps.pdf',
    height = 7, width = 10.8)
     plot_grid(MDR_Typhi, MDR_Paratyphi, FQNS_Typhi, FQNS_Paratyphi, nrow = 2 , rel_widths = c(2, 1.36, 2, 1.36), labels = c('a', 'b', 'c', 'd'))
dev.off()

png('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/data_availability/data_availabilty_maps.png',
    height = 7, width = 10.8, units = 'in', res = 300)
plot_grid(MDR_Typhi, MDR_Paratyphi, FQNS_Typhi, FQNS_Paratyphi, nrow = 2 , rel_widths = c(2, 1.36, 2, 1.36),
          labels = c('a', 'b', 'c', 'd'))
dev.off()



#plot out the number of studies
studies_by_loc <- unique(master.data[,.(nid, location_id)])
studies_by_loc <- studies_by_loc[,.(data_points = .N),
                                by = c('location_id')]
hist(studies_by_loc$data_points)
studies_by_loc$data_points_group <- NA
studies_by_loc$data_points_group[studies_by_loc$data_points <=5] <- "1-5"
studies_by_loc$data_points_group[studies_by_loc$data_points>=6 & studies_by_loc$data_points <=10] <- "6-10"
studies_by_loc$data_points_group[studies_by_loc$data_points>=11 & studies_by_loc$data_points <=50] <- "11-50"
studies_by_loc$data_points_group[studies_by_loc$data_points>=51 & studies_by_loc$data_points <=100] <- "51-100"
studies_by_loc$data_points_group[studies_by_loc$data_points>100] <- ">100"
studies_by_loc$data_points_group <-  factor(studies_by_loc$data_points_group, levels = c('Non-endemic', 'No data', "1-5", "6-10", "11-50", "51-100", ">100"))

studies_by_year <- unique(master.data[,.(nid, year_id)])

studies_by_loc <-  merge(my_shp, studies_by_loc, by.x = 'loc_id', by.y = 'location_id')

all_studies <- 
  ggplot()+
  geom_sf(data = my_shp, fill = '#737373',colour = 'black', size = 0.25)+
  geom_sf(data = background, fill = '#737373',colour = 'black', size = 0.25)+
  geom_sf(data = endemic_locs, fill = '#bdbdbd',colour = 'black', size = 0.25)+
  geom_sf(data = studies_by_loc, aes(fill = data_points_group),colour = 'black', size = 0.25)+
  scale_fill_manual(values = c('#ffffcc', '#c2e699', '#78c679', '#31a354', '#006837'),
                    labels = c("1-5","6-10", "11-50", "51-100", ">100"))+
  theme_bw()+
  theme(line = element_blank(),
        axis.text = element_blank(),
        legend.position=c(.1,.22),
        legend.key.size = unit(0.2, "cm"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank())+ 
  xlim(-20,150)+
  ylim(-35,55)+
  labs(fill = 'Data sources')

all_studies_bar <-
  ggplot()+
  geom_bar(data = studies_by_year, aes(x = year_id), fill = '#43a2ca', stat = "count")+
    theme_bw()+
    scale_x_continuous("Year", 
                       breaks = seq(1990, 2018, 1),
                       labels = seq(1990, 2018, 1))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          line = element_blank(),
          axis.title=element_text(size=8))+
    labs(y = 'Number of data sources')

png('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/data_availability/total_data_sources.png',
    height = 10, width = 15, units = 'cm', res = 300)
  all_studies
dev.off()

png('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/data_availability/total_data_sources_bar.png',
    height = 5, width = 15, units = 'cm', res = 300)
all_studies_bar
dev.off()
