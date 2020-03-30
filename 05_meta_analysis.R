#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Run a random effects metanalysis by GBD region                      #
# output tables with I^2 and pooled prevalence and median prevalence  # 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rm(list = ls())
library(metafor)
library(data.table)
setwd("Z:/AMR/Pathogens/typhi_paratyphi/model_prep")

master.data <- readRDS("clean_data/full_cleaned_dataset.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Limit data to required variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
master.data$time_period <-  NA
master.data$time_period[master.data$mid_year < 1995] <- "1990-1994"
master.data$time_period[master.data$mid_year>=1995 & master.data$mid_year < 1999] <- "1995-1999"
master.data$time_period[master.data$mid_year>=2000 & master.data$mid_year < 2004] <- "2000-2004"
master.data$time_period[master.data$mid_year>=2005 & master.data$mid_year < 2009] <- "2005-2009"
master.data$time_period[master.data$mid_year>=2010 & master.data$mid_year < 2014] <- "2010-2014"
master.data$time_period[master.data$mid_year>=2015] <- "2015-2018"
master.data <-  master.data[!is.na(master.data$time_period),]
master.data <- master.data[master.data$sample_size > 10, ]
master.data <- master.data[,c("region",
                             "iso3",
                             "mid_year",
                             "time_period",
                             "serotype",
                             "antimicrobial",
                             "sample_size",
                             "number_resistant",
                             "percentage_resistant",
                             "source_id")]

master.data$antimicrobial[grep("FQR", master.data$antimicrobial)] <- "FQR"
master.data$antimicrobial[master.data$antimicrobial == "multi-drug" | master.data$antimicrobial == "MDR - calculated"] <- "MDR"

#Order by time location and serotype
master.data <- master.data[order(master.data$time_period, master.data$region, master.data$serotype, master.data$iso3, master.data$mid_year),]

#Transform the data using thr Freeman-Turkey double arcsine transformation
master.data <- master.data[master.data$number_resistant<master.data$sample_size,]
master.data<- escalc(measure = "PFT", xi = master.data$number_resistant, ni = master.data$sample_size, data = master.data, add = 0)

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
  subset1 <- master.data[which(master.data$antimicrobial==antimicrobial.list[a]),]
  
  for(tp in 1:length(time.period.list)){
    subset2 <- subset1[which(subset1$time_period == time.period.list[tp]),]
    
    for(sp in 1:length(serotype.list)){
      # subset3 <- subset2[subset2$serotype == serotype.list[sp],] 
      subset3 <- subset2[grep(serotype.list[sp], subset2$serotype),] 
      
      for(rg in 1:length(region.list)){
        mydata <- subset3[which(subset3$region == region.list[rg]),]
        
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
saveRDS(results, 'meta_analysis/meta_analysis.RDS')

#Reshape tables to be wide by time-period
results <- data.table(results)
# Drop regions not in systematic review:
results <- results[results$region == "Eastern Sub-Saharan Africa" | 
                     results$region == "North Africa & Middle East" | 
                     results$region == "South Asia"| 
                     results$region == "Southeast Asia" | 
                     results$region =="East Asia" | 
                     results$region == "Western Sub-Saharan Africa" | 
                     results$region =="Central Sub-Saharan Africa",]
#rename sSA regions
results$region[results$region=='Central Sub-Saharan Africa'] <- 'Sub-Saharan Africa, Central'
results$region[results$region=='Eastern Sub-Saharan Africa'] <- 'Sub-Saharan Africa, Eastern'
results$region[results$region=='Western Sub-Saharan Africa'] <- 'Sub-Saharan Africa, Western'

metaanalysis <- results[, c("serotype",
                          "antimicrobial",
                          "region",                
                          "time.period",
                          "Pooled prevalence [95% CI]")]

metaanalysis <- dcast(metaanalysis, antimicrobial+serotype+region~time.period, value.var = 'Pooled prevalence [95% CI]')
metaanalysis <-  metaanalysis[order(-metaanalysis$serotype, -metaanalysis$antimicrobial, metaanalysis$region)]

I_2 <- results[,c("serotype",
                  "antimicrobial",
                  "region",                
                  "time.period",
                  "I.2")]

I_2 <- dcast(I_2, antimicrobial+serotype+region~time.period, value.var = 'I.2')
I_2 <-  I_2[order(-I_2$serotype, -I_2$antimicrobial, I_2$region)]

#save tables
write.csv(metaanalysis, "meta_analysis/pooled_prevalence.csv", row.names = F, na = '-')
write.csv(I_2, "meta_analysis/I_squared.csv", row.names = F, na = '-')


