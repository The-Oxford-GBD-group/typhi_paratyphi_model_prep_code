#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# perform meta-analysis on subsets #
# of the data and compare for      #
# sensitivity analysis             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rm(list = ls())
library(metafor)
library(data.table)
library(ggplot2)
setwd("Z:/AMR/Pathogens/typhi_paratyphi/model_prep/")

master.data <- readRDS("clean_data/full_cleaned_dataset.rds")
master.data$time_period <-  NA
master.data$time_period[master.data$mid_year < 1995] <- "1990-1994"
master.data$time_period[master.data$mid_year>=1995 & master.data$mid_year < 1999] <- "1995-1999"
master.data$time_period[master.data$mid_year>=2000 & master.data$mid_year < 2004] <- "2000-2004"
master.data$time_period[master.data$mid_year>=2005 & master.data$mid_year < 2009] <- "2005-2009"
master.data$time_period[master.data$mid_year>=2010 & master.data$mid_year < 2014] <- "2010-2014"
master.data$time_period[master.data$mid_year>=2015] <- "2015-2018"
master.data <-  master.data[!is.na(master.data$time_period),]
master.data <- master.data[master.data$sample_size > 10, ]

#get the resistance breakpoints
master.data$resistance_breakpoints[master.data$resistance_breakpoints == 'not specified' | master.data$resistance_breakpoints == 'Not specified'] <- 'No'
master.data$resistance_breakpoints[master.data$resistance_breakpoints != 'No'] <- 'Yes'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Limit data to required variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
master.data <- master.data[,c("region",
                              "iso3",
                              "location_name",
                              "mid_year",
                              "time_period",
                              "serotype",
                              "antimicrobial",
                              "sample_size",
                              "number_resistant",
                              "percentage_resistant",
                              "source_id",
                              'resistance_breakpoints')]

master.data$antimicrobial[grep("FQR", master.data$antimicrobial)] <- "FQR"
master.data$antimicrobial[master.data$antimicrobial == "multi-drug" | master.data$antimicrobial == "MDR - calculated"] <- "MDR"

#Order by time location and serotype
master.data <- master.data[order(master.data$time_period, master.data$region, master.data$serotype, master.data$iso3, master.data$mid_year),]

#Transform the data using thr Freeman-Turkey double arcsine transformation
master.data <- master.data[master.data$number_resistant<master.data$sample_size,]
master.data<- escalc(measure = "PFT", xi = master.data$number_resistant, ni = master.data$sample_size, data = master.data, add = 0)

master.data <- master.data[master.data$resistance_breakpoints == 'Yes',]
  
#~~~~~~~~~~~~~~~#
# META-ANALYSIS #
#~~~~~~~~~~~~~~~#
  #Create empty dataframe for results (and one to bind on for each loop)
  
  results = data.frame(matrix(vector(), 0, 10,
                              dimnames=list(c(), c("serotype", "antimicrobial", "region", "time.period","pooled.prevalence", "lower.CI", "upper.CI", "I^2", "N_studies", "total_sampled"))),
                       stringsAsFactors=F)
  
  resultsA = data.frame(matrix(vector(), 1, 10,
                               dimnames=list(c(), c("serotype", "antimicrobial", "region", "time.period","pooled.prevalence", "lower.CI", "upper.CI", "I^2", "N_studies", "total_sampled"))),
                        stringsAsFactors=F)
  
  # Loop through each subgroup and perform metaanalysis
  
  antimicrobial.list <- c("MDR", "FQR")
  serotype.list <- c("Typhi", "Paratyphi")
  region.list <- unique(master.data$region)
  time.period.list <- unique(master.data$time_period)
  
  for(a in 1:length(antimicrobial.list)){
    subset1 <- master.data[master.data$antimicrobial==antimicrobial.list[a],]
    
    for(tp in 1:length(time.period.list)){
      subset2 <- subset1[subset1$time_period == time.period.list[tp],]
      
      for(sp in 1:length(serotype.list)){
        subset3 <- subset2[subset2$serotype == serotype.list[sp],] 
        
        for(rg in 1:length(region.list)){
          mydata <- subset3[subset3$region == region.list[rg],]
          
          if(length(mydata$antimicrobial)==0){
            ##Add results to the empty dataframe for each subgroup
            resultsA$time.period <- time.period.list[tp]
            resultsA$region <- region.list[rg]
            resultsA$antimicrobial <- antimicrobial.list[a]
            resultsA$serotype <- serotype.list[sp]
            resultsA$pooled.prevalence <- NA
            resultsA$lower.CI <- NA
            resultsA$upper.CI <- NA
            resultsA$I.2 <- NA
            resultsA$p <- NA
            resultsA$N_studies <- length(mydata$sample_size)
            resultsA$total_sampled <- sum(mydata$sample_size)
            results <- rbind(results, resultsA)
            
          } else if(length(mydata$antimicrobial)==1){
            
            #Back transform the data to get the  95% CIs for each study
            dat.back <- summary(mydata, transf=transf.ipft, ni=mydata$sample_size)
            
            ##Add results to the empty dataframe for each subgroup
            resultsA$time.period <- time.period.list[tp]
            resultsA$region <- region.list[rg]
            resultsA$antimicrobial <- antimicrobial.list[a]
            resultsA$serotype <- serotype.list[sp]
            resultsA$pooled.prevalence <- dat.back$yi
            resultsA$lower.CI <- dat.back$ci.lb
            resultsA$upper.CI <- dat.back$ci.ub
            resultsA$I.2 <- NA
            resultsA$p <- NA
            resultsA$N_studies <- length(mydata$sample_size)
            resultsA$total_sampled <- sum(mydata$sample_size)
            results <- rbind(results, resultsA)
            
          } else{
            
            #Back transform the data to get the  95% CIs for each study
            dat.back <- summary(mydata, transf=transf.ipft, ni=mydata$sample_size)
            
            #Perform meta-analysis (random effectds, weighted)
            res <- try(rma(yi, vi, method = "REML", data = mydata, weighted = TRUE))
            if(inherits(res, "try-error")){
              res <- rma(yi, vi, method = "REML", data = mydata, weighted = TRUE, control=list(stepadj=0.8))
            }
            
            #Tranform the pooled estimate back to prevalence using the harmonic mean of the sample size
            pred <- predict(res, transf=transf.ipft.hm, targs=list(ni=mydata$sample_size))
            
            ##Add results to the empty dataframe for each subgroup
            resultsA$time.period <- time.period.list[tp]
            resultsA$region <- region.list[rg]
            resultsA$antimicrobial <- antimicrobial.list[a]
            resultsA$serotype <- serotype.list[sp]
            resultsA$pooled.prevalence <- pred$pred
            resultsA$lower.CI <- pred$ci.lb
            resultsA$upper.CI <- pred$ci.ub
            resultsA$I.2 <- res$I2
            resultsA$p <- res$pval
            resultsA$N_studies <- res$k
            resultsA$total_sampled <- sum(dat.back$sample_size)
            results <- rbind(results, resultsA)
            
          } 
        }
      }
    }
  }
  
  #convert results to percentages and round to 2 decimal places
  results$pooled.prevalence <- round(results$pooled.prevalence*100,0)
  results$upper.CI <- round(results$upper.CI*100,0)
  results$lower.CI <- round(results$lower.CI*100,0)
  results$I.2 <- round(results$I.2,0)
  results$p <- round(results$p, 3)
  results$p[results$p<0.001] <- 'p<0.001'
  results$p[results$p!= 'p<0.001'&
              !is.na(results$p)] <- paste0('p=',results$p[results$p!= 'p<0.001'&
                                                            !is.na(results$p)])
  
  
  #Paste together prevalence and CIs for tables
  results$`Pooled prevalence [95% CI]`[results$N_studies>1] <- paste0(results$pooled.prevalence[results$N_studies>1], "% [", results$lower.CI[results$N_studies>1], "-", results$upper.CI[results$N_studies>1], "]")
  results$`Pooled prevalence [95% CI]`[results$N_studies==0] <- NA
  results$`Pooled prevalence [95% CI]`[results$N_studies==1] <- paste0(results$pooled.prevalence[results$N_studies==1], "% [", results$lower.CI[results$N_studies==1], "-", results$upper.CI[results$N_studies==1], "]*")
  
  #save as an RDS for forest plot
  saveRDS(results, 'sensitivity_analysis/RB/meta_analysis_RB_yes.RDS')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Stats test for significan differences between results #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rm(list = ls())

full_meta <- readRDS('meta_analysis/meta_analysis.RDS')
RB <- readRDS('sensitivity_analysis/RB/meta_analysis_RB_yes.RDS')
colnames(full_meta)[5] <- 'full_meta' 
colnames(RB)[5] <- 'RB_yes' 

all.data <- merge(full_meta, RB, by = c('serotype', 'antimicrobial', 'region', 'time.period'), all.x = T, all.y =T)

#T-test for all regions
t.test(all.data$full_meta, all.data$RB_yes)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Create tables for all regions #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
full_meta <- full_meta[full_meta$region =="Central Sub-Saharan Africa" |
                         full_meta$region == "East Asia" |
                         full_meta$region == "Eastern Sub-Saharan Africa" |
                         full_meta$region =="North Africa & Middle East"|
                         full_meta$region =="South Asia" |
                         full_meta$region == "Southeast Asia" |
                         full_meta$region == "Western Sub-Saharan Africa",]

fm <-  full_meta[c("serotype",
                   "antimicrobial",
                   "region",
                   'time.period',
                   "Pooled prevalence [95% CI]")]


RB <- RB[RB$region =="Central Sub-Saharan Africa" |
             RB$region == "East Asia" |
             RB$region == "Eastern Sub-Saharan Africa" |
             RB$region =="North Africa & Middle East"|
             RB$region =="South Asia" |
             RB$region == "Southeast Asia" |
             RB$region == "Western Sub-Saharan Africa",]

RB_sub <-  RB[c("serotype",
                  "antimicrobial",
                  "region",
                  'time.period',
                  "Pooled prevalence [95% CI]")]

colnames(fm)[5] <- 'Full model'
colnames(RB_sub)[5] <- 'Resistance breakpoints stated'

all.table <- merge(fm, RB_sub)

#get the t-test p values for all subgroups
t.test.table <- unique(all.table[c("serotype", "antimicrobial", 'region')])
t.test.table$p_value1 <- NA

serotype <- unique(t.test.table$serotype)
antimicrobial <- unique(t.test.table$antimicrobial)
region <- unique(t.test.table$region)

for(sp in serotype){
  subset1 <- all.data[all.data$serotype== sp,]
  
  for(a in antimicrobial){
    subset2 <- subset1[subset1$antimicrobial == a,]
    
    for(reg in region){
      subset3 <- subset2[subset2$region == reg,]
      
      if(length(subset3$full_meta[!is.na(subset3$full_meta)])>1 & 
         length(subset3$RB_yes[!is.na(subset3$RB_yes)])>1){
        t.result1 <- t.test(subset3$full_meta, subset3$RB_yes)
        t.result1 <- round(t.result1$p.value, 3)  
        t.test.table$p_value1[t.test.table$serotype == sp &
                                t.test.table$antimicrobial == a &
                                t.test.table$region == reg] <- t.result1
      } 
      
    }  
  }
}


all.table <- merge(all.table, t.test.table, by = c('serotype', 'antimicrobial', 'region'))
all.table$p_value1[all.table$time.period!="2015-2018"] <- NA

write.csv(all.table, 'sensitivity_analysis/RB/full_SA_table.csv', row.names = F, na = "")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# make plots for all subgroups #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rm(list =ls())
full_meta <- readRDS('meta_analysis/meta_analysis.RDS')
RB <- readRDS('sensitivity_analysis/RB/meta_analysis_RB_yes.RDS')

RB$subgroup <- 'Resistance breakpoints reported'
full_meta$subgroup <- 'All studies'

all.data <- rbind(full_meta, RB)
serotype <- unique(all.data$serotype)
antimicrobial <- unique(all.data$antimicrobial)
region = c("East Asia", "Southeast Asia", "South Asia",  "North Africa & Middle East",
           "Central Sub-Saharan Africa", "Eastern Sub-Saharan Africa",  "Western Sub-Saharan Africa")

for(sp in serotype){
  subset1 <- all.data[all.data$serotype== sp,]
  
  for(a in antimicrobial){
    subset2 <- subset1[subset1$antimicrobial == a,]
    
    for(reg in region){
      subset3 <- subset2[subset2$region == reg,]
      
      jpeg(paste0('sensitivity_analysis/RB/', a, "_", sp, "_", reg, '.jpeg'),
           width =20, height = 15,  units = "cm", res = 300)
      print(ggplot(subset3, aes(x=subgroup, y=pooled.prevalence, colour = subgroup)) +
              geom_pointrange(aes(ymin = lower.CI, ymax = upper.CI))+
              facet_wrap(~time.period)+
              theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank())+
              theme(legend.title=element_blank())+
              labs(y = paste0(a, ' (%)'))+
              ggtitle(paste(a, sp, reg, sep = " "))
      )
      dev.off()
    }
  }
}

