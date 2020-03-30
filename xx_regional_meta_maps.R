rm(list = ls())
library(rgeos)
library(raster)
library(ggplot2)
library(maptools)
library(rgdal)
library(latticeExtra)
library(gridExtra)
library(data.table)
######################################################
#read in GBD dissolaved shapefile
regions_shp <- readOGR('K:/AMR/Misc/GBD_regions_global.shp')
regions_shp@data$region <-  as.character(regions_shp@data$region)

#merge on the data
mydata <- readRDS('K:/AMR/Pathogens/typhi_paratyphi/model_prep/meta_analysis/meta_analysis.RDS')

mydata <- mydata[c(1:5)]
mydata <- mydata[!is.na(mydata$pooled.prevalence),]
mydata <- mydata[mydata$region!='Andean Latin America',]
mydata <- mydata[mydata$region!='Central Asia',]
mydata <- mydata[mydata$region!='Southern Sub-Saharan Africa',]
mydata <- mydata[mydata$region!='Southern Latin America',]
mydata <- mydata[mydata$region!='High-income Asia Pacific',]
mydata <- mydata[mydata$region!='Central Latin America',]


mydata <- data.table(mydata)

MDR_typhi <- mydata[mydata$antimicrobial == 'MDR' &
                      mydata$serotype == 'Typhi',]

MDR_typhi <- dcast(MDR_typhi, region ~ time.period, value.var = 'pooled.prevalence')

FQNS_typhi <- mydata[mydata$antimicrobial == 'FQR' &
                      mydata$serotype == 'Typhi',]

FQNS_typhi <- dcast(FQNS_typhi, region ~ time.period, value.var = 'pooled.prevalence')


# Create non-endemic shapefile
nonend <-  regions_shp[regions_shp$region=='Andean Latin America' |
                         regions_shp$region=="Australasia"  |
                         regions_shp$region=="Caribbean"  |
                         regions_shp$region=="Central Asia"|
                         regions_shp$region==  "Central Europe" |
                         regions_shp$region == "Central Latin America"|
                         regions_shp$region ==   "Eastern Europe"| 
                         regions_shp$region == "High-income Asia Pacific"|
                         regions_shp$region == "High-income North America"|
                         regions_shp$region == 'Southern Latin America'|
                         regions_shp$region == 'Southern Sub-Saharan Africa'| 
                         regions_shp$region == 'Tropical Latin America'|
                         regions_shp$region ==   'Western Europe',]



#merge data onto the endmeic shapefile
colnames(MDR_typhi) <- c('region', 'MDR_1990', 'MDR_1995', 'MDR_2000', 'MDR_2005', 'MDR_2010', 'MDR_2015')
colnames(FQNS_typhi) <- c('region', 'FQNS_1990', 'FQNS_1995', 'FQNS_2000', 'FQNS_2005', 'FQNS_2010', 'FQNS_2015')

merged_shp <- merge(regions_shp, MDR_typhi, by = "region", all.x = F) 
merged_shp <- merge(merged_shp, FQNS_typhi, by = "region", all.x = F) 

#set up pdf and plot all time periods per drug/bug on one sheet.
jpeg('K:/AMR/Pathogens/typhi_paratyphi/model_prep/meta_maps/regional/MDR typhi maps.jpg',
     width = 21,
     height = 30,
     units = 'cm',
     res = 150)
grid.arrange(
  spplot(merged_shp, zcol = 'MDR_1990', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '1990-1994')+
    layer_(sp.polygons(merged_shp, fill = '#636363'))+
    layer_(sp.polygons(nonend, fill = '#bdbdbd')),
  
  spplot(merged_shp, zcol = 'MDR_1995', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '1995-1999')+
    layer_(sp.polygons(merged_shp, fill = '#636363'))+
    layer_(sp.polygons(nonend, fill = '#bdbdbd')),
  
  spplot(merged_shp, zcol = 'MDR_2000', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '2000-2004')+
    layer_(sp.polygons(merged_shp, fill = '#636363'))+
    layer_(sp.polygons(nonend, fill = '#bdbdbd')),
  
  spplot(merged_shp, zcol = 'MDR_2005', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '2005-2009')+
    layer_(sp.polygons(merged_shp, fill = '#636363'))+
    layer_(sp.polygons(nonend, fill = '#bdbdbd')),
  
  spplot(merged_shp, zcol = 'MDR_2010', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '2010-2014')+
    layer_(sp.polygons(merged_shp, fill = '#636363'))+
    layer_(sp.polygons(nonend, fill = '#bdbdbd')),
  
  ncol=1,
  top = 'MDR Typhi'
)
dev.off()

jpeg('K:/AMR/Pathogens/typhi_paratyphi/model_prep/meta_maps/regional/FQNS typhi maps.jpg',
     width = 21,
     height = 30,
     units = 'cm',
     res = 150)
grid.arrange(
  spplot(merged_shp, zcol = 'FQNS_1990', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '1990-1994')+
    layer_(sp.polygons(merged_shp, fill = '#636363'))+
    layer_(sp.polygons(nonend, fill = '#bdbdbd')),
  
  spplot(merged_shp, zcol = 'FQNS_1995', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '1995-1999')+
    layer_(sp.polygons(merged_shp, fill = '#636363'))+
    layer_(sp.polygons(nonend, fill = '#bdbdbd')),
  
  spplot(merged_shp, zcol = 'FQNS_2000', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '2000-2004')+
    layer_(sp.polygons(merged_shp, fill = '#636363'))+
    layer_(sp.polygons(nonend, fill = '#bdbdbd')),
  
  spplot(merged_shp, zcol = 'FQNS_2005', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '2005-2009')+
    layer_(sp.polygons(merged_shp, fill = '#636363'))+
    layer_(sp.polygons(nonend, fill = '#bdbdbd')),
  
  spplot(merged_shp, zcol = 'FQNS_2010', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '2010-2014')+
    layer_(sp.polygons(merged_shp, fill = '#636363'))+
    layer_(sp.polygons(nonend, fill = '#bdbdbd')),
  
  ncol=1,
  top = 'FQNS Typhi'
)
dev.off()

#also plot as single figs for gif

jpeg('K:/AMR/Pathogens/typhi_paratyphi/model_prep/meta_maps/regional/MDR typhi/1990-1994.jpg',
     width = 21,
     height = 6,
     units = 'cm',
     res = 150)
  spplot(merged_shp, zcol = 'MDR_1990', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '1990-1994')+
    layer_(sp.polygons(merged_shp, fill = '#636363'))+
    layer_(sp.polygons(nonend, fill = '#bdbdbd'))
dev.off()

jpeg('K:/AMR/Pathogens/typhi_paratyphi/model_prep/meta_maps/regional/MDR typhi/1995-1999.jpg',
     width = 21,
     height = 6,
     units = 'cm',
     res = 150)
  spplot(merged_shp, zcol = 'MDR_1995', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '1995-1999')+
    layer_(sp.polygons(merged_shp, fill = '#636363'))+
    layer_(sp.polygons(nonend, fill = '#bdbdbd'))
dev.off()
  
jpeg('K:/AMR/Pathogens/typhi_paratyphi/model_prep/meta_maps/regional/MDR typhi/2000-2004.jpg',
       width = 21,
       height = 6,
       units = 'cm',
       res = 150)
  spplot(merged_shp, zcol = 'MDR_2000', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '2000-2004')+
    layer_(sp.polygons(merged_shp, fill = '#636363'))+
    layer_(sp.polygons(nonend, fill = '#bdbdbd'))
dev.off()

jpeg('K:/AMR/Pathogens/typhi_paratyphi/model_prep/meta_maps/regional/MDR typhi/2005-2009.jpg',
       width = 21,
       height = 6,
       units = 'cm',
       res = 150)
  spplot(merged_shp, zcol = 'MDR_2005', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '2005-2009')+
    layer_(sp.polygons(merged_shp, fill = '#636363'))+
    layer_(sp.polygons(nonend, fill = '#bdbdbd'))
dev.off()
  
jpeg('K:/AMR/Pathogens/typhi_paratyphi/model_prep/meta_maps/regional/MDR typhi/2010-2014.jpg',
       width = 21,
       height = 6,
       units = 'cm',
       res = 150)
  spplot(merged_shp, zcol = 'MDR_2010', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '2010-2014')+
    layer_(sp.polygons(merged_shp, fill = '#636363'))+
    layer_(sp.polygons(nonend, fill = '#bdbdbd'))
dev.off()


#FQNS
jpeg('K:/AMR/Pathogens/typhi_paratyphi/model_prep/meta_maps/regional/FQNS typhi/1990-1994.jpg',
     width = 21,
     height = 6,
     units = 'cm',
     res = 150)
spplot(merged_shp, zcol = 'FQNS_1990', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '1990-1994')+
  layer_(sp.polygons(merged_shp, fill = '#636363'))+
  layer_(sp.polygons(nonend, fill = '#bdbdbd'))
dev.off()

jpeg('K:/AMR/Pathogens/typhi_paratyphi/model_prep/meta_maps/regional/FQNS typhi/1995-1999.jpg',
     width = 21,
     height = 6,
     units = 'cm',
     res = 150)
spplot(merged_shp, zcol = 'FQNS_1995', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '1995-1999')+
  layer_(sp.polygons(merged_shp, fill = '#636363'))+
  layer_(sp.polygons(nonend, fill = '#bdbdbd'))
dev.off()

jpeg('K:/AMR/Pathogens/typhi_paratyphi/model_prep/meta_maps/regional/FQNS typhi/2000-2004.jpg',
     width = 21,
     height = 6,
     units = 'cm',
     res = 150)
spplot(merged_shp, zcol = 'FQNS_2000', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '2000-2004')+
  layer_(sp.polygons(merged_shp, fill = '#636363'))+
  layer_(sp.polygons(nonend, fill = '#bdbdbd'))
dev.off()

jpeg('K:/AMR/Pathogens/typhi_paratyphi/model_prep/meta_maps/regional/FQNS typhi/2005-2009.jpg',
     width = 21,
     height = 6,
     units = 'cm',
     res = 150)
spplot(merged_shp, zcol = 'FQNS_2005', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '2005-2009')+
  layer_(sp.polygons(merged_shp, fill = '#636363'))+
  layer_(sp.polygons(nonend, fill = '#bdbdbd'))
dev.off()

jpeg('K:/AMR/Pathogens/typhi_paratyphi/model_prep/meta_maps/regional/FQNS typhi/2010-2014.jpg',
     width = 21,
     height = 6,
     units = 'cm',
     res = 150)
spplot(merged_shp, zcol = 'FQNS_2010', col.regions = rev(heat.colors(101)), at=seq(0, 100, length=100), main = '2010-2014')+
  layer_(sp.polygons(merged_shp, fill = '#636363'))+
  layer_(sp.polygons(nonend, fill = '#bdbdbd'))
dev.off()