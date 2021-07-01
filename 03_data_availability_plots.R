#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot out number of studies per country year #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(ggplot2)
rm(list = ls())

setwd("C:/Users/Annie/Documents/GRAM/typhi_paratyphi/model_prep//")

drugbug <- c('MDR_typhi', 'MDR_paratyphi', 'FQNS_typhi', 'FQNS_paratyphi', 'ceph_typhi', 'ceph_paratyphi', 'all_studies')

for(db in drugbug){
#read in cleaned csv
if(db == 'all_studies'){
        mydata <- readRDS('clean_data/full_cleaned_dataset.rds')
        mydata <- unique(mydata, by = c('super_region', 'region', 'country', 'mid_year', 'source_id'))
} else {
        mydata <- readRDS(paste0("clean_data/", db, ".rds"))
}

mydata <-  mydata[mydata$super_region != 'High Income',]
mydata$country[mydata$country=="Democratic Republic of the Congo"] <- "DRC"
mydata$country[mydata$country=="Malawi Mozambique border"] <- "Malawi/Mozambique"

#exclude africa from paratyphi
if(db == 'MDR_paratyphi' | db == 'FQNS_paratyphi'| db == 'ceph_paratyphi'){
        mydata <- mydata[mydata$super_region!='Sub-Saharan Africa',]
}

#Number of studies per year
studies <- mydata[,.(studies = .N), 
                   by = c('super_region', 'region', 'country', 'mid_year')]

studies$samplesize <- "1"
studies$samplesize[studies$studies>=2 & studies$studies<=5] <- "2-5"
studies$samplesize[studies$studies>=6 & studies$studies<=10] <- "6-10"
studies$samplesize[studies$studies>10] <- ">10"
studies$samplesize <- as.factor(studies$samplesize)
studies$samplesize <- factor(studies$samplesize, levels = c("1", "2-5", "6-10", ">10"))

#order countries
studies$country <- as.factor(studies$country)
studies$country <- factor(studies$country, levels=rev(levels(studies$country)))

png(paste0("data_availability/plots/", db, ".png"), 
     width =30, height = 30,  units = "cm", res = 300)

print(ggplot(studies, aes(x=mid_year, y=country, colour = super_region, order = factor(country)))+
        geom_vline(xintercept = c(1990, 1995, 2000, 2005, 2010, 2015), color = "black", size = 0.5)+
        geom_point(aes(size = factor(samplesize)))+
        # scale_colour_manual(values = c("#ACD788",
        #                                '#C7E2AD',
        #                                '#743128',
        #                                '#AE8275',
        #                                "#F1592A",
        #                                "#FEAE19",
        #                                "#FFCC7B",
        #                                "#95D6F6",
        #                                "#FBDD00",
        #                                "#ACD7D0",
        #                                "#83C5C1",
        #                                '#66B3AB',
        #                                "#12A89D"))+
        scale_x_continuous("Year", 
                           breaks = seq(1990, 2017, 5),
                           labels = c("1990", "1995", "2000", "2005", "2010", "2015"))+
        theme_bw()+
        ylab('Country')+
        labs(size = "Number of studies", colour = "GBD Region")+
        theme(plot.title = element_text(hjust = 0.5))+
        theme(legend.text = element_text(size = 8), legend.title = element_text(size = 8))+
        facet_wrap(~super_region, ncol = 2, scales='free')      
)
dev.off()
}
