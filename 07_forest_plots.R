rm(list = ls())

#library(forestplot)
library(ggplot2)
library(metafor)
library(DescTools)
library(data.table)
library(readxl)

setwd("K:/AMR/Pathogens/typhi_paratyphi/model_prep/")
dir.create('forest_plots/MDR_Typhi/', showWarnings = F)
dir.create('forest_plots/FQR_Typhi/', showWarnings = F)
dir.create('forest_plots/MDR_Paratyphi/', showWarnings = F)
dir.create('forest_plots/FQR_Paratyphi/', showWarnings = F)

master.data <- readRDS("clean_data/clean_data_crosswalked.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Limit data to required variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

master.data <- master.data[,c("region",
                             "iso3",
                             "site_name",
                             "mid_year",
                             "time_period",
                             "species",
                             "antimicrobial",
                             "no_examined",
                             "number_resistant",
                             "percentage_resistant",
                             "source_id")]


master.data$antimicrobial[grep("FQR", master.data$antimicrobial)] <- "FQR"
master.data$antimicrobial[master.data$antimicrobial=="multi-drug"|master.data$antimicrobial=="multi-drug-adjusted"] <- "MDR"

# Clean and shorten some site names
master.data$site_name[is.na(master.data$site_name)] <- "Unspecified"
master.data$site_name[master.data$site_name==""] <- "Unspecified"
master.data$site_name[master.data$site_name=="Wards 29 & 30, Kolkata"] <- "Kolkata"
master.data$site_name[master.data$site_name=="Wards 29 & 30, Kolkata"] <- "Kolkata"
master.data$site_name[master.data$site_name=="Hijrat Colony, Sultanabad and Bilal Colony, Karachi"] <- "Karachi"
master.data$site_name[master.data$site_name=="Hijrat Colony & Sultanabad"] <- "Karachi"
master.data$site_name[master.data$site_name=="Kolkata, Delhi, Vellore, Karnataka & Nagpur"] <- "Multiple sites*"
master.data$site_name[master.data$site_name=="Nairobi, Embu & Nakuru"] <- "Multiple sites*"
master.data$site_name[master.data$site_name=="Ahmedabad, Vadodara & Bhavangar"] <- "Multiple sites*"
master.data$site_name[master.data$site_name=="Sharquia, Port Said, Quena and Sohag governates"] <- "Multiple sites*"
master.data$site_name[master.data$site_name=="Mingalartaungnyunt, Thingangyun, Tharkayta, Dalla, Kyauktada, Kyaukmyaung, Dawbone, Tarmwe, Hlaing, Hlaingtharyar, N. Okkalapa, Mingaladon, Laputta,"] <- "Yangon"
master.data$site_name[master.data$site_name=="Kampong Thom, Kampong Cham,Kampong Speu, Phnom Penh, Prey Veng, Kandal, Svay Rieng, Takeo, Kampot & Sihanoukville"] <- "Phnom Penh"
master.data$site_name[master.data$site_name=="Patan, Lalipur"] <- "Kathmandu"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Calculate 95% CIs using agresti-coull method #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#(method chosen instead of arcsince as the has upper CI <100 if % resistant == 100)
master.data<- escalc(measure = "PFT", xi = master.data$number_resistant, ni = master.data$no_examined, data = master.data, add = 0)
dat.back <- summary(master.data, transf=transf.ipft, ni=master.data$no_examined)

master.data$ci.lb <- dat.back$ci.lb*100
master.data$ci.ub <- dat.back$ci.ub*100

master.data <- master.data[order(master.data$species, master.data$antimicrobial, master.data$region, master.data$mid_year, master.data$iso3),]
master.data$result <- paste0(round(master.data$percentage_resistant,0), "[", round(master.data$ci.lb,0), "-", round(master.data$ci.ub,0), "]")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Read in the metanalysis results #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
MA <- readRDS('meta_analysis/meta_analysis.RDS')
MA <- MA[MA$N_studies > 0,]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Loop through all regions and create forest plots with median prevalences for 5 years #
# Not including south asia yet as need to split up further
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

antimicrobial.list <- c("MDR", "FQR")
species.list <- c("Typhi", "Paratyphi")

for(a in 1:length(antimicrobial.list)){
  subset1 <- master.data[master.data$antimicrobial==antimicrobial.list[a],]
  if(antimicrobial.list[a] == "FQR") {
    res.title = "Non-susceptibility (%)"
  } else{
    res.title = "Resistance (%)"
  }
  
  for(sp in 1:length(species.list)){
    subset2 <- subset1[subset1$species == species.list[sp],] 
    
    region.list <- unique(subset2$region)
    region.list <- region.list[region.list!="South Asia"]
    for(reg in 1:length(region.list)){
      mydata <- subset2[subset2$region == region.list[reg],]
      
      fileout <- paste0('forest_plots/',
                        antimicrobial.list[a],"_",
                        species.list[sp], "/",
                        region.list[reg], ".jpeg")
      
      fileout <- gsub(" ", "_", fileout)
      fileout <- gsub("&_", "", fileout)
      
      meta <- MA[(MA$antimicrobial==antimicrobial.list[a] &
                    MA$serotype == species.list[sp] &
                    MA$region == region.list[reg]),]
      
      timeperiods <- data.frame(table(mydata$time_period))
      tp <- timeperiods$Freq
      timeperiods <- as.character(timeperiods$Var1)
      k <- length(mydata$source_id)
      j <- length(timeperiods)
      gap <- -(7.8/(k+j-1))
      rows <- seq(7.8, 0, gap)
      
      
      for(i in 1:j){
        t <- tp[i]
        
        if(i==1){
          rowpos <- rows[1:t]
        } else {
          rowpos <- c(rowpos, rows[(length(rowpos)+i):(length(rowpos)+i+t-1)])
        }
      }
      
      polyrows <- rows[!(rows%in%rowpos)]
      meta$time.period <- gsub("_", "-", meta$time.period)
      #timeperiods <- paste0("Pooled prevalence [95% CI] ", timeperiods) 
      summarypolys <- paste0("Pooled prevalence [95% CI] for ", meta$time.period, '; I^2=', meta$I.2, '%; ', meta$p)
      summarypolys[meta$N_studies<3] <- paste0('No summary measure available; ', meta$time.period[meta$N_studies<3])
      
      jpeg(fileout,
           width =30, height = 15, units = "cm", res = 300)
      
      forest(x = mydata$percentage_resistant,
             ci.lb = mydata$ci.lb,
             ci.ub = mydata$ci.ub,
             xlim = c(-80, 100),
             ylim = c(0, 10),
             annotate = F,
             alim = c(0,250),
             at = seq(0,100,10),
             digits = c(0),
             xlab = ("Percentage resistant (%)"),
             ilab = cbind(mydata$source_id, mydata$iso3, mydata$site_name, mydata$mid_year, mydata$no_examined, mydata$result),
             ilab.pos = c(4, 4, 4, 4, 4, 2),
             ilab.xpos = c(-80, -55, -45, -20, -10, 120),
             cex = .5,
             refline = NA,
             rows = rowpos,
             slab = NA)
      
      addpoly.default(x  = meta$pooled.prevalence,
                      ci.lb = meta$lower.CI,
                      ci.ub = meta$upper.CI,
                      rows = polyrows,
                      digits = c(0),
                      mlab = summarypolys,
                      annotate = F,
                      cex = 0.5,
                      font = 2,
                      col = "blue")
      
      text(c(-80, -55, -45, -20, -10), 8.2, pos = 4, cex = 0.5, font = 2, c("Study", "Country", "Site", "Year", "N"))
      text(120, 8.4, pos = 2, cex = 0.5, font = 2, res.title)
      text(120, 8.2, pos = 2, cex = 0.5, font = 2, "[95% CI]")
      text(120, polyrows, 
           pos = 2, cex = 0.5, font = 2, 
           meta$`Pooled prevalence [95% CI]`)
      
      dev.off()
      
    }
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# For South Asia create seperate plots for each decade #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

southasia <- master.data[master.data$region == "South Asia",]

for(a in antimicrobial.list){
  subset1 <- southasia[southasia$antimicrobial==a,]
  if(a == "FQR") {
    res.title = "Non-susceptibility (%)"
  } else{
    res.title = "Resistance (%)"
  }
  
  for(sp in species.list){
    subset2 <- subset1[subset1$species == sp,] 
    # get the pooled prevalence
    meta <- MA[(MA$antimicrobial==a &
                  MA$serotype == sp &
                  MA$region == 'South Asia'),]
    
    for(saplot in 1:3){
      if(saplot == 1){
        mydata <-subset2[subset2$mid_year<2000,]
        meta_sub <- meta[meta$time.period=='1990_1994'|
                           meta$time.period=='1995_1999',]        
        
        fileout <- paste0('forest_plots/',
                          a,"_",
                          sp, "/",
                          "South_Asia_1990-1999.jpeg")
        
        
      } else if(saplot == 2){
        mydata <-subset2[subset2$mid_year>=2000 &
                           subset2$mid_year<2010,]
        meta_sub <- meta[meta$time.period=='2000_2004'|
                     meta$time.period=='2005_2009',]        
        fileout <- paste0('forest_plots/',
                          a,"_",
                          sp, "/",
                          "South_Asia_2000-2009.jpeg")
      } else {
        mydata <-subset2[subset2$mid_year>=2010,]
        meta_sub <- meta[meta$time.period=='2010_2014'|
                           meta$time.period=='2015_2017',]        
        
        fileout <- paste0('forest_plots/',
                          a,"_",
                          sp, "/",
                          "South_Asia_2010-2017.jpeg")
      }
      

      timeperiods <- data.frame(table(mydata$time_period))
      tp <- timeperiods$Freq
      timeperiods <- timeperiods$Var1
      k <- length(mydata$source_id)
      j <- length(timeperiods)
      gap <- -(7.8/(k+j-1))
      rows <- seq(7.8, 0, gap)

      for(i in 1:j){
        t <- tp[i]
        
        if(i==1){
          rowpos <- rows[1:t]
        } else {
          rowpos <- c(rowpos, rows[(length(rowpos)+i):(length(rowpos)+i+t-1)])
        }
      }
      
      polyrows <- rows[!(rows%in%rowpos)]
      meta_sub$time.period <- gsub('_', '-', meta_sub$time.period) 
      summarypolys <- paste0("Pooled prevalence [95% CI] for ", meta_sub$time.period, '; I^2=', meta_sub$I.2, '%; ', meta_sub$p)
      summarypolys[meta_sub$N_studies<3] <- paste0('No summary measure available; ', meta_sub$time.period[meta_sub$N_studies<3])

      jpeg(fileout,
           width =30, height = 30, units = "cm", res = 300)
      
      forest(x = mydata$percentage_resistant,
             ci.lb = mydata$ci.lb,
             ci.ub = mydata$ci.ub,
             xlim = c(-80, 100),
             ylim = c(0, 10),
             annotate = F,
             alim = c(0,250),
             at = seq(0,100,10),
             digits = c(0),
             xlab = ("Percentage resistant (%)"),
             ilab = cbind(mydata$source_id, mydata$iso3, mydata$site_name, mydata$mid_year, mydata$no_examined, mydata$result),
             ilab.pos = c(4, 4, 4, 4, 4, 2),
             ilab.xpos = c(-80, -55, -45, -20, -10, 120),
             cex = .5,
             refline = NA,
             rows = rowpos,
             slab = NA)
      
      addpoly.default(x  = meta_sub$pooled.prevalence,
                      ci.lb = meta_sub$lower.CI,
                      ci.ub = meta_sub$upper.CI,
                      rows = polyrows,
                      digits = c(0),
                      mlab = summarypolys,
                      annotate = F,
                      cex = 0.5,
                      font = 2,
                      col = "blue")
      
      text(c(-80, -55, -45, -20, -10), 8.2, pos = 4, cex = 0.5, font = 2, c("Study", "Country", "Site", "Year", "N"))
      text(120, 8.4, pos = 2, cex = 0.5, font = 2, res.title)
      text(120, 8.2, pos = 2, cex = 0.5, font = 2, "[95% CI]")
      text(120, polyrows, 
           pos = 2, cex = 0.5, font = 2, 
           meta_sub$`Pooled prevalence [95% CI]`)
      
      dev.off()
    }
  }
}


# For paratyphyi FQR as only have 1 study <2000

mydata <-southasia[southasia$mid_year<2009 & southasia$antimicrobial == 'FQR' & southasia$species == 'Paratyphi',]
meta_sub <- MA[MA$time.period=='1990_1994'|
               MA$time.period=='1995_1999'|
               MA$time.period=='2000_2004'|
               MA$time.period=='2005_2009',]        

meta_sub <-  meta_sub[meta_sub$region == 'South Asia'&
                        meta_sub$antimicrobial == 'FQR'&
                        meta_sub$serotype == 'Paratyphi',]

fileout <- paste0('forest_plots/FQR_Paratyphi/South_Asia_1995-2009.jpeg')

res.title = "Non-susceptibility (%)"

timeperiods <- data.frame(table(mydata$time_period))
tp <- timeperiods$Freq
timeperiods <- timeperiods$Var1
k <- length(mydata$source_id)
j <- length(timeperiods)
gap <- -(7.8/(k+j-1))
rows <- seq(7.8, 0, gap)

for(i in 1:j){
  t <- tp[i]
  
  if(i==1){
    rowpos <- rows[1:t]
  } else {
    rowpos <- c(rowpos, rows[(length(rowpos)+i):(length(rowpos)+i+t-1)])
  }
}

polyrows <- rows[!(rows%in%rowpos)]
meta_sub$time.period <- gsub('_', '-', meta_sub$time.period) 
#timeperiods <- paste0("Pooled prevalence [95% CI] ", timeperiods) 
summarypolys <- paste0("Pooled prevalence [95% CI] for ", meta_sub$time.period, '; I^2=', meta_sub$I.2, '%; ', meta_sub$p)
summarypolys[meta_sub$N_studies<3] <- paste0('No summary measure available; ', meta_sub$time.period[meta_sub$N_studies<3])


jpeg(fileout,
     width =30, height = 30, units = "cm", res = 300)

forest(x = mydata$percentage_resistant,
       ci.lb = mydata$ci.lb,
       ci.ub = mydata$ci.ub,
       xlim = c(-80, 100),
       ylim = c(0, 10),
       annotate = F,
       alim = c(0,250),
       at = seq(0,100,10),
       digits = c(0),
       xlab = ("Percentage resistant (%)"),
       ilab = cbind(mydata$source_id, mydata$iso3, mydata$site_name, mydata$mid_year, mydata$no_examined, mydata$result),
       ilab.pos = c(4, 4, 4, 4, 4, 2),
       ilab.xpos = c(-80, -55, -45, -20, -10, 120),
       cex = .5,
       refline = NA,
       rows = rowpos,
       slab = NA)

addpoly.default(x  = meta_sub$pooled.prevalence,
                ci.lb = meta_sub$lower.CI,
                ci.ub = meta_sub$upper.CI,
                rows = polyrows,
                digits = c(0),
                mlab = summarypolys,
                annotate = F,
                cex = 0.5,
                font = 2,
                col = "blue")

text(c(-80, -55, -45, -20, -10), 8.2, pos = 4, cex = 0.5, font = 2, c("Study", "Country", "Site", "Year", "N"))
text(120, 8.4, pos = 2, cex = 0.5, font = 2, res.title)
text(120, 8.2, pos = 2, cex = 0.5, font = 2, "[95% CI]")
text(120, polyrows, 
     pos = 2, cex = 0.5, font = 2, 
     meta_sub$`Pooled prevalence [95% CI]`)

dev.off()

