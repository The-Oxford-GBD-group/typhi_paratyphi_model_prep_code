rm(list = ls())
library(lme4)
library(foreign)
library(ggplot2)
library(DDoutlier)
library(scales)
library(ggrepel)
library(stats)
library(ggforce)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Identify outliers in datasets using MAD and DDoutlier methods ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#read in data
setwd('C:/Users/Annie/Documents/GRAM/typhi_paratyphi')
mydata <- fread('model_prep/clean_data/ceph_paratyphi_clean.csv')
mydata$number_resistant <- round(mydata$number_resistant,0)

#match to national level covariates
covs <- read.csv('covariates/cleaned_covs.csv')

covs <- covs[names(covs) %in% c('location_id', 'year_id', 
                                "cv_mean_temperature",
                                "cv_sanitation_prop",
                                'cv_maternal_educ_yrs_pc',
                                "J01D")]

#centre scale covs
covs$year <- covs$year_id
covs[,3:7] <- scale(covs[,3:7])

#merge data and covs
mydata <- merge(mydata, covs, by = c('location_id', 'year_id'), all.x = T, all.y = F)

response = cbind(successes = mydata$number_resistant,
                 failures = mydata$sample_size)

model1 <- glmer(response ~ 1 + year +
               cv_mean_temperature+
               cv_sanitation_prop+
               cv_maternal_educ_yrs_pc+
               (1 |super_region / region / country), 
                data = mydata, 
                family = 'binomial',
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

summary(model1)

#predict out to all endmic locations
locs <- read.dbf('C:/Users/Annie/Documents/GRAM/shapefiles/typhi_endemic.dbf')
locs <- locs[c('loc_id', 'spr_reg_id', 'region_id', 'ihme_lc_id')]
locs[1:4] <- sapply(locs[1:4], as.character)
locs[1:3] <- sapply(locs[1:3], as.numeric)


covs <- merge(covs, locs, by.x = 'location_id', by.y = 'loc_id')
colnames(covs)[colnames(covs) == 'ihme_lc_id'] <- 'country'
colnames(covs)[colnames(covs) == 'region_id'] <- 'region'
colnames(covs)[colnames(covs) == 'spr_reg_id'] <- 'super_region'
covs$region[covs$region == 5] <- 'East Asia '
covs$region[covs$region == 9] <- 'Southeast Asia'
covs$region[covs$region == 159] <- 'South Asia'

covs$super_region[covs$super_region == 158] <- 'South Asia'
covs$super_region[covs$super_region == 4] <- 'Southeast Asia, East Asia & Oceania'

covs <- covs[covs$super_region == 'South Asia'|
               covs$super_region == 'Southeast Asia, East Asia & Oceania',]

covs$QA <- 1
covs$pred <- predict(model1, newdata = covs, type = 'response', allow.new.levels = TRUE)
summary(covs$pred)

covs <- merge(covs, mydata[,.(location_id, year_id, val, sample_size)], by = c('location_id', 'year_id'), all.x = T, all.y = T)

#determine the MAD
covs <- data.table(covs)
covs <- covs[!is.na(covs$super_region),]

MADs <-  covs[,.(upper_bound = pred + 3*mad(pred[!is.na(val)], val[!is.na(val)]),
                 lower_bound = pred - 3*mad(pred[!is.na(val)],val[!is.na(val)])),
              by = c('country')]

covs <- cbind(covs, MADs[,2:3])
covs$upper_bound[covs$upper_bound>1] <- 1
covs$lower_bound[covs$lower_bound<0] <- 0

pdf(paste0('model_prep/outliering/ceph_paratyphi/MAD2.pdf'),
    height = 17,
    width = 12) 

for(i in 1:(ceiling(length(unique(covs$super_region))))){ 
  print(ggplot(covs[covs$super_region==unique(covs$super_region)[i],])+ 
          geom_line(aes(x = year, y = pred))+
          geom_ribbon(aes(x = year, ymin = lower_bound, ymax = upper_bound, colour = 'red', fill = 'red', alpha  =0.5))+
          geom_point(aes(x = year, y = val))+
          facet_wrap_paginate(~country, page = i)
  )} 

dev.off()

#define outliers in the dataset
MADs <- covs[,.(country, year, lower_bound, upper_bound)]
MADs <-  unique(MADs)

mydata <- merge(mydata, MADs, by = c('country', 'year'))
mydata$MAD_outlier[mydata$val<mydata$lower_bound |mydata$val>mydata$upper_bound] <- 1
mydata$MAD_outlier[mydata$val>mydata$lower_bound & mydata$val<mydata$upper_bound] <- 0

#DD outliers ####
#rescale values
mydata$year.rescaled <- rescale(mydata$year_id)
mydata$location_id.rescaled <- rescale(mydata$location_id)
df <- mydata[,.(val, year.rescaled, location_id.rescaled)]

##Natural Neighbor (NAN) algorithm to return the self-adaptive neighborhood
K <- NAN(df, NaN_Edges=FALSE)$r

#calculate outlier scores for 14 algorithms
outlier_score_COF <- COF(dataset=df, k=K)
outlier_score_INFLO <- INFLO(dataset=df, k=K)
outlier_score_KDEOS <- KDEOS(dataset=df, k_min=K, k_max=15)
outlier_score_KNN_AGG <- KNN_AGG(dataset=df, k_min=K, k_max=15)
outlier_score_KNN_IN <- KNN_IN(dataset=df, k=K)
outlier_score_KNN_SUM <- KNN_SUM(dataset=df, k=K)
outlier_score_LDF <- LDF(dataset=df, k=K, h=2, c=1)$LDF
outlier_score_LDOF <- LDOF(dataset=df, k=K)
outlier_score_LOCI <- LOCI(dataset=df, alpha=0.5, nn=20, k=K)$norm_MDEF
outlier_score_LOF <- LOF(dataset=df, k=K)
outlier_score_LOOP <- LOOP(dataset=df, k=K, lambda=3)
outlier_score_NOF <- NOF(dataset=df)$NOF
outlier_score_RDOS <- RDOS(dataset=df, k=K, h=2)
outlier_score_RKOF_r <- RKOF(dataset=df, k = K, C = 1, alpha = 1, sigma2 = 1)

#Defined which are outliers in each algorithm
COF_outlier <- ifelse(outlier_score_COF>=quantile(outlier_score_COF, 0.95), 1,0) 
INFLO_outlier <- ifelse(outlier_score_INFLO>=quantile(outlier_score_INFLO, 0.95), 1,0) 
KDEOS_outlier <- ifelse(outlier_score_KDEOS>=quantile(outlier_score_KDEOS, 0.95), 1,0) 
KNN_AGG_outlier <- ifelse(outlier_score_KNN_AGG>=quantile(outlier_score_KNN_AGG, 0.95), 1,0) 
KNN_IN_outlier <- ifelse(outlier_score_KNN_IN>=quantile(outlier_score_KNN_IN, 0.95), 1,0) 
KNN_SUM_outlier <- ifelse(outlier_score_KNN_SUM>=quantile(outlier_score_KNN_SUM, 0.95), 1,0) 
LDF_outlier <- ifelse(outlier_score_LDF>=quantile(outlier_score_LDF, 0.95), 1,0) 
LDOF_outlier <- ifelse(outlier_score_LDOF>=quantile(outlier_score_LDOF, 0.95), 1,0) 
LOCI_outlier <- ifelse(outlier_score_LOCI>=quantile(outlier_score_LOCI, 0.95), 1,0) 
LOF_outlier <- ifelse(outlier_score_LOF>=quantile(outlier_score_LOF, 0.95), 1,0) 
LOOP_outlier <- ifelse(outlier_score_LOOP>=quantile(outlier_score_LOOP, 0.95), 1,0) 
NOF_outlier <- ifelse(outlier_score_NOF>=quantile(outlier_score_NOF, 0.95), 1,0) 
RDOS_outlier <- ifelse(outlier_score_RDOS>=quantile(outlier_score_RDOS, 0.95), 1,0) 
RKOF_r_outlier <- ifelse(outlier_score_RKOF_r>=quantile(outlier_score_RKOF_r, 0.95), 1,0) 

#get a comined DD outlier score
dd_outlier <- COF_outlier+INFLO_outlier+KDEOS_outlier+KNN_AGG_outlier+KNN_IN_outlier+KNN_SUM_outlier+
  LDF_outlier+LDOF_outlier+LOCI_outlier+LOF_outlier+LOOP_outlier+NOF_outlier+RDOS_outlier+RKOF_r_outlier+LOCI_outlier

#identify those which are outliers in at least 2 algorithms
mydata$dd_outlier <- ifelse(dd_outlier >=2,1,0)
mydata$outliers <- 'None'
mydata$outliers[mydata$dd_outlier==1] <- 'DD_outlier'
mydata$outliers[mydata$MAD_outlier==1] <- 'MAD_outlier'
mydata$outliers[mydata$dd_outlier==1 & mydata$MAD_outlier==1] <- 'Both'

mydata$outliers <- factor(mydata$outliers, levels = c('None', 'DD_outlier', 'MAD_outlier', 'Both'))

pdf('model_prep/outliering/ceph_paratyphi/DD_over2.pdf',
    height = 8.3, width = 11.7)
#plot out a page for each region
for(i in 1:length(unique(mydata$super_region))){
  subset <- mydata[mydata$super_region == unique(mydata$super_region)[i],]
  print(
    ggplot(subset)+
      geom_point(aes(x = year_id, y = val, colour = as.factor(outliers)))+
      facet_wrap(~country)+
      ylim(0,1)+
      ylab('Proportion DR')+
      theme_bw()+
      scale_x_continuous("Year", 
                         breaks = seq(1990, 2018, 5),
                         labels = c("1990", "1995", "2000", "2005", "2010", "2015"))
  )
}
dev.off()

#inspect this plot and select outliers
mydata <- mydata[,.(nid,
                    super_region,
                    region,
                    country,
                    location_name,
                    location_id,
                    subnat_location_id,
                    adj_id, adj_id_sSA, adj_id_Asia,
                    year_id,
                    age_group_id,
                    sex_id,
                    measure_id,
                    number_resistant,
                    sample_size,
                    val,
                    variance,
                    QA,
                    outliers)]

write.csv(mydata , "model_prep/outliering/ceph_paratyphi/inspect_outliers.csv", row.names = F)

#Inspect the output file and determine which of the MAD and DD outliers should be outliered
