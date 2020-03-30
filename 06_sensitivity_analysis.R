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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Limit data to required variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
master.data <- master.data[,c("region",
                              "iso3",
                              "mid_year",
                              "time_period",
                              "serotype",
                              "antimicrobial",
                              "sample_size",
                              "number_resistant",
                              "percentage_resistant",
                              "source_id",
                              'QA')]

master.data$antimicrobial[grep("FQR", master.data$antimicrobial)] <- "FQR"
master.data$antimicrobial[master.data$antimicrobial == "multi-drug" | master.data$antimicrobial == "MDR - calculated"] <- "MDR"

#Order by time location and serotype
master.data <- master.data[order(master.data$time_period, master.data$region, master.data$serotype, master.data$iso3, master.data$mid_year),]

#Transform the data using thr Freeman-Turkey double arcsine transformation
master.data <- master.data[master.data$number_resistant<master.data$sample_size,]
master.data<- escalc(measure = "PFT", xi = master.data$number_resistant, ni = master.data$sample_size, data = master.data, add = 0)

for(i in 1:2) { 
  if(i==1){
    master.data <- master.data[master.data$QA <3,]
  }else{
    master.data <- master.data[master.data$QA == 1,]
  }
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
if(i==1){
  saveRDS(results, 'sensitivity_analysis/QA/meta_analysis_QA1_2.RDS')
}else{
  saveRDS(results, 'sensitivity_analysis/QA/meta_analysis_QA1.RDS')
}
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Stats test for significan differences between results #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rm(list = ls())

full_meta <- readRDS('meta_analysis/meta_analysis.RDS')
QA1 <- readRDS('sensitivity_analysis/QA/meta_analysis_QA1.RDS')
QA1_2 <- readRDS('sensitivity_analysis/QA/meta_analysis_QA1_2.RDS')
colnames(full_meta)[5] <- 'full_meta' 
colnames(QA1)[5] <- 'QA1' 
colnames(QA1_2)[5] <- 'QA1_2' 

all.data <- merge(full_meta, QA1, by = c('serotype', 'antimicrobial', 'region', 'time.period'), all.x = T, all.y =T)
all.data <- merge(all.data, QA1_2, by = c('serotype', 'antimicrobial', 'region', 'time.period'), all.x = T, all.y =T)

#T-test for all regions
t.test(all.data$full_meta, all.data$QA1)
t.test(all.data$full_meta, all.data$QA1_2)

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


QA1 <- QA1[QA1$region =="Central Sub-Saharan Africa" |
                         QA1$region == "East Asia" |
                         QA1$region == "Eastern Sub-Saharan Africa" |
                         QA1$region =="North Africa & Middle East"|
                         QA1$region =="South Asia" |
                         QA1$region == "Southeast Asia" |
                         QA1$region == "Western Sub-Saharan Africa",]

QA1_sub <-  QA1[c("serotype",
                   "antimicrobial",
                   "region",
                   'time.period',
                   "Pooled prevalence [95% CI]")]

QA1_2 <- QA1_2[QA1_2$region =="Central Sub-Saharan Africa" |
             QA1_2$region == "East Asia" |
             QA1_2$region == "Eastern Sub-Saharan Africa" |
             QA1_2$region =="North Africa & Middle East"|
             QA1_2$region =="South Asia" |
             QA1_2$region == "Southeast Asia" |
             QA1_2$region == "Western Sub-Saharan Africa",]

QA1_2_sub <-  QA1_2[c("serotype",
                  "antimicrobial",
                  "region",
                  'time.period',
                  "Pooled prevalence [95% CI]")]


colnames(fm)[5] <- 'Full model'
colnames(QA1_sub)[5] <- 'Quality level 1'
colnames(QA1_2_sub)[5] <- 'Quality level 1 & 2'


all.table <- merge(fm, QA1_sub)
all.table <- merge(all.table, QA1_2_sub)

#get the t-test p values for all subgroups
t.test.table <- unique(all.table[c("serotype", "antimicrobial", 'region')])
t.test.table$p_value1 <- NA
t.test.table$p_value2 <- NA

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
         length(subset3$QA1[!is.na(subset3$QA1)])>1){
      t.result1 <- t.test(subset3$full_meta, subset3$QA1)
      t.result1 <- round(t.result1$p.value, 3)  
      t.test.table$p_value1[t.test.table$serotype == sp &
                              t.test.table$antimicrobial == a &
                              t.test.table$region == reg] <- t.result1
      } 
      
      if(length(subset3$full_meta[!is.na(subset3$full_meta)])>1 & 
         length(subset3$QA1_2[!is.na(subset3$QA1_2)])>1){
      t.result2 <- t.test(subset3$full_meta, subset3$QA1_2)
      t.result2 <- round(t.result2$p.value,3)  
      t.test.table$p_value2[t.test.table$serotype == sp &
                              t.test.table$antimicrobial == a &
                              t.test.table$region == reg] <- t.result2
      }
    }  
  }
}


all.table <- merge(all.table, t.test.table)
all.table$p_value1[all.table$time.period!="2015-2018"] <- NA
all.table$p_value2[all.table$time.period!="2015-2018"] <- NA

write.csv(all.table, 'sensitivity_analysis/QA/full_SA_table.csv', row.names = F, na = "")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# make plots for all subgroups #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rm(list =ls())
full_meta <- readRDS('meta_analysis/meta_analysis.RDS')
QA1 <- readRDS('sensitivity_analysis/QA/meta_analysis_QA1.RDS')
QA1_2 <- readRDS('sensitivity_analysis/QA/meta_analysis_QA1_2.RDS')

QA1$subgroup <- 'Level 1 quality data'
QA1_2$subgroup <- 'Level 1 and 2 quality data'
full_meta$subgroup <- 'All studies'

all.data <- rbind(full_meta, QA1, QA1_2)
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
      
      jpeg(paste0('sensitivity_analysis/QA/', a, "_", sp, "_", reg, '.jpeg'),
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

      