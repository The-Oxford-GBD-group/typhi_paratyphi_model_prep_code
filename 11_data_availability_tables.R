#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Summarise the data and study characteristics for systematic review #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(doBy)
library(data.table)
rm(list = ls())

setwd("K:/AMR/Pathogens/typhi_paratyphi/model_prep/")

#read in cleaned csv
master.data <- readRDS("clean_data/clean_data_crosswalked.rds")

master.data$antimicrobial[grep("FQR", master.data$antimicrobial)] <- "FQR"
master.data$antimicrobial[master.data$antimicrobial=="multi-drug" | master.data$antimicrobial=="multi-drug-adjusted"] <- "MDR"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Table for number of studies all studies, all serotypes #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#Limit to required variables
master.data <- master.data[,.(no_examined,
                             number_resistant,
                             percentage_resistant,
                             antimicrobial,
                             species,
                             region,
                             country,
                             site_name,
                             source_id,
                             mid_year,
                             time_period)]

#Order by gbd region so all tables are in the same order
master.data <- master.data[order(master.data$region),]

antimicrobial.list <- c("All", "MDR", "FQR")
species.list <- c('All', "Typhi", "Paratyphi")

regions <- c("Central Europe",
             "Central Asia",
             "Eastern Europe",
             "Australasia",
             "Western Europe",
             "High-income Asia Pacific",
             "East Asia",
             "High-income North America",
             "Andean Latin America",
             "Tropical Latin America",
             "Southern Latin America",
             "Central Latin America",
             "Caribbean",
             "North Africa & Middle East",
             "South Asia",
             "Southeast Asia",
             "Oceania",
             "Central Sub-Saharan Africa",
             "Western Sub-Saharan Africa",
             "Southern Sub-Saharan Africa",
             "Eastern Sub-Saharan Africa")

for(a in antimicrobial.list){
  if(a == 'All'){
    subset1 <- master.data
  } else {
  subset1 <- master.data[master.data$antimicrobial== a,]
  }
  
  for(sp in species.list){
    if(sp == 'All'){
      mydata <- subset1
    } else {
    mydata <- subset1[subset1$species == sp,]
    }
    
    mydata <- data.table(mydata)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Get the number of studies per region/time period #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    my.studies <- unique(mydata, by = c("source_id", "region", "time_period"))
    total.studies <- length(unique(my.studies$source_id))
    
    my.studies <- my.studies[, .(studies = length(source_id)),
                             by = c("region", 
                                    "time_period")]
    
    my.studies <- dcast(my.studies, region ~ time_period, value.var = 'studies')
    
    colnames(my.studies) <- gsub("_", "-", colnames(my.studies))
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Get the number of isolates per region/time period #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Get the max number of isolates tested by each study location then sum them by region
      my.samples <- mydata[,.(no_examined = max(no_examined)),
                           by = c("region",
                                  'country',
                                  'site_name',
                                  "time_period",
                                  'source_id',
                                  'species')]
      
      my.samples <- my.samples[,.(no_examined = sum(no_examined)),
                               by = c("region", 
                                  "time_period")]
      
    my.samples <- dcast(my.samples, region ~ time_period, value.var = 'no_examined')
    
    colnames(my.samples) <- gsub("_", "-", colnames(my.samples))
    
    # Combine tables
    mytable = data.frame(regions)
    mytable$regions <- as.character(mytable$regions)
    colnames(mytable) <- "region"
    mytable$`1990-1994`  <- NA
    mytable$`1995-1999`  <- NA
    mytable$`2000-2004`  <- NA
    mytable$`2005-2009`  <- NA
    mytable$`2010-2014`  <- NA
    mytable$`2015-2018`  <- NA
    mytable <- mytable[order(mytable$region),]
  
    if(!("1990-1994"%in%colnames(my.studies))){
      my.studies$`1990-1994` <- NA
      my.samples$`1990-1994` <- NA
    }
    if(!("1995-1999"%in%colnames(my.studies))){
      my.studies$`1995-1999` <- NA
      my.samples$`1995-1999` <- NA
    }
    if(!("2000-2004"%in%colnames(my.studies))){
      my.studies$`2000-2004` <- NA
      my.samples$`2000-2004` <- NA
    }
    if(!("2005-2009"%in%colnames(my.studies))){
      my.studies$`2005-2009` <- NA
      my.samples$`2005-2009` <- NA
    }
    if(!("2010-2014"%in%colnames(my.studies))){
      my.studies$`2010-2014` <- NA
      my.samples$`2010-2014` <- NA
    }
    if(!("2015-2018"%in%colnames(my.studies))){
      my.studies$`2015-2018` <- NA
      my.samples$`2015-2018` <- NA
    }
    
    my.studies[is.na(my.studies)] <- 0
    my.samples[is.na(my.samples)] <- 0
    my.studies$region <- as.character(my.studies$region)
    my.samples$region <- as.character(my.samples$region)
    
    my.studies <- data.frame(my.studies)
    my.samples <- data.frame(my.samples)
    
    colnames(my.studies) <- c('region', "1990-1994", "1995-1999", "2000-2004", "2005-2009", "2010-2014", "2015-2018")
    colnames(my.samples) <- c('region', "1990-1994", "1995-1999", "2000-2004", "2005-2009", "2010-2014", "2015-2018")
    
    #ensure all regions are in the table
    for(j in 1:length(regions)){
      if(regions[j]%in%my.studies$region){
      } else{
        my.studies <- rbind(my.studies, c(regions[j], 0, 0, 0, 0, 0, 0))
        my.samples <- rbind(my.samples, c(regions[j], 0, 0, 0, 0, 0, 0))
        my.studies <- my.studies[order(my.studies$region),]
        my.samples <- my.samples[order(my.samples$region),]
        }
    }       

    colsnum <- c("1990-1994", "1995-1999", "2000-2004", "2005-2009", "2010-2014", "2015-2018")
    my.studies[colsnum] <- sapply(my.studies[colsnum], as.numeric)
    my.samples[colsnum] <- sapply(my.samples[colsnum], as.numeric)
    
    #get totals
    row.total.studies <- rowSums(my.studies[2:7])
    row.total.samples <- rowSums(my.samples[2:7])
    col.total.studies <- colSums(my.studies[,-1])
    col.total.samples <- colSums(my.samples[,-1])
    
    total.samples <- sum(col.total.samples)
    
    row.total.samples <- prettyNum(row.total.samples, big.mark = ",")
    col.total.samples <- prettyNum(col.total.samples, big.mark = ",")
    total.studies <- prettyNum(total.studies, big.mark = ",")
    total.samples <- prettyNum(total.samples, big.mark = ",")
    
    row.total <- paste0(row.total.studies, "(", row.total.samples, ')')
    col.total <- paste0(col.total.studies, "(", col.total.samples, ')')
    total.total <- paste0(total.studies, "(", total.samples, ')')
    
    #format with , 
    my.studies[colsnum] <- sapply(my.studies[colsnum], function(x)prettyNum(x, big.mark = ','))
    my.samples[colsnum] <- sapply(my.samples[colsnum], function(x)prettyNum(x, big.mark = ','))
    
    ## populate table with studies (samples)
    mytable$`1990-1994` <- paste0(my.studies$`1990-1994`, "(", my.samples$`1990-1994`, ")")
    mytable$`1995-1999` <- paste0(my.studies$`1995-1999`, "(", my.samples$`1995-1999`, ")")
    mytable$`2000-2004` <- paste0(my.studies$`2000-2004`, "(", my.samples$`2000-2004`, ")")
    mytable$`2005-2009` <- paste0(my.studies$`2005-2009`, "(", my.samples$`2005-2009`, ")")
    mytable$`2010-2014` <- paste0(my.studies$`2010-2014`, "(", my.samples$`2010-2014`, ")")
    mytable$`2015-2018` <- paste0(my.studies$`2015-2018`, "(", my.samples$`2015-2018`, ")")
    
    ##add totals
    mytable$Total <- row.total
    mytable <- rbind(mytable, c('Total', col.total, total.total))

    fileout <- paste0("data_availability/tables/", sp, "_", a, "_number_of_studies.csv")
    write.csv(mytable, fileout, row.names = F)
  
  }
}
