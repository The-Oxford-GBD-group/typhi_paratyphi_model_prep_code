library(ggforce)
library(ggplot2)
source("H:/Functions/round2.R")

rm(list = ls())


clean_data <- readRDS('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/full_cleaned_dataset.rds') 
clean_data <- clean_data[clean_data$sample_size>=10,]
clean_data <- clean_data[clean_data$mid_year>=1990,]


testing_info <- read.csv('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/testing_info/testing_info_by_study_site.csv', stringsAsFactors = F)
demographics <- read.csv('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/demographics/demographics_by_study_site_serotype.csv', stringsAsFactors = F)
names(demographics)[names(demographics)=='species'] <- 'serotype'
study_info <- read.csv('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/study_info/study_info.csv', stringsAsFactors = F)


#restrict this to the resistance and location info then add on the aggregated testing and demographics
clean_data <- clean_data[,.(source_id,
                            location_name,
                            super_region,
                            region,
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
                            mid_year,
                            serotype,
                            antimicrobial,
                            sample_size,
                            number_resistant,
                            percentage_resistant)]

clean_data <- merge(clean_data, demographics, by = c('source_id', 'location_name', 'serotype'), allow.cartesian = TRUE)
clean_data <- merge(clean_data, testing_info, by = c('source_id', 'location_name'), allow.cartesian = TRUE)
clean_data <- merge(clean_data, study_info, by = c('source_id'), allow.cartesian = TRUE)

# Set a sample size bins
clean_data$sample_size_bins <- NA
clean_data$sample_size_bins[clean_data$sample_size<50] <- '10-49'
clean_data$sample_size_bins[clean_data$sample_size>=50 &clean_data$sample_size<100 ] <- '50-99'
clean_data$sample_size_bins[clean_data$sample_size>=100 &clean_data$sample_size<500 ] <- '100-499'
clean_data$sample_size_bins[clean_data$sample_size>=500 ] <- '500+'
clean_data$sample_size_bins <-  as.factor(clean_data$sample_size_bins)
clean_data$sample_size_bins <- factor(clean_data$sample_size_bins, levels = c("10-49", "50-99", "100-499", "500+"))
clean_data <- clean_data[!is.na(clean_data$sample_size_bins)]


#restrinct to anitbiotics required
MDR_typhi <- clean_data[(clean_data$antimicrobial == 'multi-drug' | clean_data$antimicrobial == 'MDR - calculated') &
                       clean_data$serotype == 'Typhi',]


MDR_paratyphi <- clean_data[(clean_data$antimicrobial == 'multi-drug' | clean_data$antimicrobial == 'MDR - calculated') &
                           clean_data$serotype == 'Paratyphi',]

FQNS_typhi <- clean_data[(clean_data$antimicrobial == 'FQR - cipro & nalidixic acid' | clean_data$antimicrobial == ' FQR - ciprofloxacin' | clean_data$antimicrobial == 'FQR - nalidixic acid') &
                        clean_data$serotype == 'Typhi',]



FQNS_paratyphi <- clean_data[(clean_data$antimicrobial == 'FQR - cipro & nalidixic acid' | clean_data$antimicrobial == ' FQR - ciprofloxacin'|clean_data$antimicrobial == 'FQR - nalidixic acid') &
                            clean_data$serotype == 'Paratyphi',]


table(MDR_typhi$ast)
MDR_typhi$ast[MDR_typhi$ast=='Disk-diffusion & Automated methods'] <- 'Disk-diffusion & MIC determination'
MDR_typhi$ast[MDR_typhi$ast=='Disk-diffusion & E-test'] <- 'Disk-diffusion & MIC determination'
MDR_typhi$ast[MDR_typhi$ast=='Disk-diffusion & E-test & Microdilution'] <- 'Disk-diffusion & MIC determination'
MDR_typhi$ast[MDR_typhi$ast=='Disk-diffusion & Microdilution'] <- 'Disk-diffusion & MIC determination'
MDR_typhi$ast[MDR_typhi$ast=='Disk-diffusion & E-test Disk-diffusion & MIC determination'] <- 'Disk-diffusion & MIC determination'

MDR_typhi$ast[MDR_typhi$ast=='E-test'] <- 'MIC determination'
MDR_typhi$ast[MDR_typhi$ast=='E-test & Automated methods'] <- 'MIC determination'
MDR_typhi$ast[MDR_typhi$ast=='E-test & Microdilution'] <- 'MIC determination'
MDR_typhi$ast[MDR_typhi$ast=='Microdilution & Automated methods'] <- 'MIC determination'
MDR_typhi$ast[MDR_typhi$ast=='Microdilution'] <- 'MIC determination'

MDR_typhi$control_strain[MDR_typhi$control_strain == 'Not specified'] <- 'No'
MDR_typhi$control_strain[MDR_typhi$control_strain != 'No'] <- 'Yes'

MDR_typhi$resistance_breakpoints[MDR_typhi$resistance_breakpoints == 'Not specified'] <- 'No'
MDR_typhi$resistance_breakpoints[MDR_typhi$resistance_breakpoints != 'No'] <- 'Yes'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# FIt stepwise linear regression to test for covariates to include ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
MDR_typhi_s_asia <- MDR_typhi[MDR_typhi$super_region == 'South Asia',]
MDR_typhi_se_asia <- MDR_typhi[MDR_typhi$super_region == 'Southeast Asia, East Asia & Oceania',]
MDR_typhi_sSA <- MDR_typhi[MDR_typhi$super_region == 'Sub-Saharan Africa',]

library(MASS)
#1. without ISO or MID YEAR

#a. All data
fit <- lm(percentage_resistant~isolate_source+ast+control_strain+resistance_breakpoints+QA+child_adult, data =MDR_typhi)
step <- stepAIC(fit, direction="both")
step$anova # display results

#b. South Asia
fit <- lm(percentage_resistant~isolate_source+ast+control_strain+resistance_breakpoints+QA++child_adult, data =MDR_typhi_s_asia)
step <- stepAIC(fit, direction="both")
step$anova # display results

#c. Southeast Asia...
fit <- lm(percentage_resistant~isolate_source+ast+control_strain+resistance_breakpoints+QA+child_adult, data =MDR_typhi_se_asia)
step <- stepAIC(fit, direction="both")
step$anova # display results

#d. sSA
fit <- lm(percentage_resistant~isolate_source+ast+control_strain+resistance_breakpoints+QA+child_adult, data =MDR_typhi_sSA)
step <- stepAIC(fit, direction="both")
step$anova # display results


#2. with MID YEAR

#a. All data
fit <- lm(percentage_resistant~mid_year+isolate_source+ast+control_strain+resistance_breakpoints+QA+child_adult, data =MDR_typhi)
step <- stepAIC(fit, direction="both")
step$anova # display results

#b. South Asia
fit <- lm(percentage_resistant~mid_year+isolate_source+ast+control_strain+resistance_breakpoints+QA+child_adult, data =MDR_typhi_s_asia)
step <- stepAIC(fit, direction="both")
step$anova # display results

#c. Southeast Asia...
fit <- lm(percentage_resistant~mid_year+isolate_source+ast+control_strain+resistance_breakpoints+QA+child_adult, data =MDR_typhi_se_asia)
step <- stepAIC(fit, direction="both")
step$anova # display results

#d. sSA
fit <- lm(percentage_resistant~mid_year+isolate_source+ast+control_strain+resistance_breakpoints+QA+child_adult, data =MDR_typhi_sSA)
step <- stepAIC(fit, direction="both")
step$anova # display results

#3. With ISO and MID YEAR

#a. All data
fit <- lm(percentage_resistant~iso3+mid_year+isolate_source+ast+control_strain+resistance_breakpoints+QA+child_adult, data =MDR_typhi)
step <- stepAIC(fit, direction="both")
step$anova # display results

#b. South Asia
fit <- lm(percentage_resistant~iso3+mid_year+isolate_source+ast+control_strain+resistance_breakpoints+QA+child_adult, data =MDR_typhi_s_asia)
step <- stepAIC(fit, direction="both")
step$anova # display results

#c. Southeast Asia...
fit <- lm(percentage_resistant~iso3+mid_year+isolate_source+ast+control_strain+resistance_breakpoints+QA+child_adult, data =MDR_typhi_se_asia)
step <- stepAIC(fit, direction="both")
step$anova # display results

#d. sSA
fit <- lm(percentage_resistant~iso3+mid_year+isolate_source+ast+control_strain+resistance_breakpoints+QA+child_adult, data =MDR_typhi_sSA)
step <- stepAIC(fit, direction="both")
step$anova # display results

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Fit Lasso to test for variables to include ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(glmnet)
#define covariates as a matrix
MDR_typhi$QA <- as.numeric(MDR_typhi$QA)
MDR_typhi$isolate_source <-  as.factor(MDR_typhi$isolate_source)
levels(MDR_typhi$isolate_source) <-  as.factor(MDR_typhi$isolate_source)
MDR_typhi$isolate_source <- factor(MDR_typhi$isolate_source, levels = c("Not specified", "Non sterile sites", "Mixed sites", "Normally sterile sites"))
MDR_typhi$isolate_source <- as.numeric(MDR_typhi$isolate_source)

MDR_typhi$control_strain[MDR_typhi$control_strain == 'No'] <- 0
MDR_typhi$control_strain[MDR_typhi$control_strain == 'Yes'] <- 1
MDR_typhi$control_strain <- as.numeric(MDR_typhi$control_strain)

MDR_typhi$ast <- as.factor(MDR_typhi$ast)
levels(MDR_typhi$ast) <-  as.factor(MDR_typhi$ast)
MDR_typhi$ast <- factor(MDR_typhi$ast, levels = c("Not specified", "Disk-diffusion", "Disk-diffusion & MIC determination", "MIC determination", "Automated methods"))
MDR_typhi$ast <-  as.numeric(MDR_typhi$ast)

MDR_typhi$resistance_breakpoints[MDR_typhi$resistance_breakpoints == 'Yes'] <- 1
MDR_typhi$resistance_breakpoints[MDR_typhi$resistance_breakpoints == 'No'] <- 0
MDR_typhi$resistance_breakpoints <- as.numeric(MDR_typhi$resistance_breakpoints)

MDR_typhi$iso3 <- as.factor(MDR_typhi$iso3)
MDR_typhi$iso3 <- as.numeric(MDR_typhi$iso3)

MDR_typhi$child_adult <-  as.factor(MDR_typhi$child_adult)
levels(MDR_typhi$child_adult)
MDR_typhi$child_adult <-  as.numeric(MDR_typhi$child_adult)

## shuffle the data into five random folds
MDR_typhi <- MDR_typhi[sample(nrow(MDR_typhi)),]
MDR_typhi[,fold_id := cut(seq(1,nrow(MDR_typhi)),breaks=5,labels=FALSE)]
MDR_typhi[, a_rowid := seq(1:nrow(MDR_typhi))]

#1. with vars but no iso3 or mid year
vars <- as.matrix(MDR_typhi[, .(isolate_source, resistance_breakpoints, ast, control_strain, QA, child_adult)], with = F)
colnames(vars) <- c('isolate_source', 'resistance_breakpoints', 'ast', 'control_strain', 'QA', 'child_adult')

response <- cbind(failures   = MDR_typhi$sample_size - MDR_typhi$number_resistant, 
                  successes = MDR_typhi$number_resistant)

#fit cross validated lasso to select lambda
cv_lasso = cv.glmnet(x = vars , y= response, family = 'binomial', alpha = 1, weights = MDR_typhi$sample_size, nfolds = 5, foldid = MDR_typhi$fold_id)

cv_lasso$lambda.1se
coef(cv_lasso, s = 0.05)

#2 Add mid-year
vars <- as.matrix(MDR_typhi[, .(mid_year, isolate_source, resistance_breakpoints, ast, control_strain, QA, child_adult)], with = F)
colnames(vars) <- c('mid_year', 'isolate_source', 'resistance_breakpoints', 'ast', 'control_strain', 'QA', 'child_adult')


response <- cbind(failures   = MDR_typhi$sample_size - MDR_typhi$number_resistant, 
                  successes = MDR_typhi$number_resistant)

#fit cross validated lasso to select lambda
cv_lasso = cv.glmnet(x = vars , y= response, family = 'binomial', alpha = 1, weights = MDR_typhi$sample_size, nfolds = 5, foldid = MDR_typhi$fold_id)

cv_lasso$lambda.1se

coef(cv_lasso, s = 0.05)

#2. Add ISO and mid_year
vars <- as.matrix(MDR_typhi[, .(iso3, mid_year, isolate_source, resistance_breakpoints, ast, control_strain, QA, child_adult)], with = F)
colnames(vars) <- c('iso3', 'mid_year', 'isolate_source', 'resistance_breakpoints', 'ast', 'control_strain', 'QA', 'child_adult')


response <- cbind(failures   = MDR_typhi$sample_size - MDR_typhi$number_resistant, 
                  successes = MDR_typhi$number_resistant)

#fit cross validated lasso to select lambda
cv_lasso = cv.glmnet(x = vars , y= response, family = 'binomial', alpha = 1, weights = MDR_typhi$sample_size, nfolds = 5, foldid = MDR_typhi$fold_id)

cv_lasso$lambda.1se

coef(cv_lasso, s = 0.05)


#repeat for SOuth asia only
MDR_typhi_s_asia <- MDR_typhi[MDR_typhi$super_region == 'South Asia',]
  
## shuffle the data into five random folds
MDR_typhi_s_asia <- MDR_typhi_s_asia[sample(nrow(MDR_typhi_s_asia)),]
MDR_typhi_s_asia[,fold_id := cut(seq(1,nrow(MDR_typhi_s_asia)),breaks=5,labels=FALSE)]
MDR_typhi_s_asia[, a_rowid := seq(1:nrow(MDR_typhi_s_asia))]

#1. with vars but no iso3 or mid year
vars <- as.matrix(MDR_typhi_s_asia[, .(isolate_source, resistance_breakpoints, ast, control_strain, QA, child_adult)], with = F)
colnames(vars) <- c('isolate_source', 'resistance_breakpoints', 'ast', 'control_strain', 'QA', 'child_adult')

response <- cbind(failures   = MDR_typhi_s_asia$sample_size - MDR_typhi_s_asia$number_resistant, 
                  successes = MDR_typhi_s_asia$number_resistant)

#fit cross validated lasso to select lambda
cv_lasso = cv.glmnet(x = vars , y= response, family = 'binomial', alpha = 1, weights = MDR_typhi_s_asia$sample_size, nfolds = 5, foldid = MDR_typhi_s_asia$fold_id)

cv_lasso$lambda.1se
coef(cv_lasso, s = 0.05)

#2 Add mid-year
vars <- as.matrix(MDR_typhi_s_asia[, .(mid_year, isolate_source, resistance_breakpoints, ast, control_strain, QA, child_adult)], with = F)
colnames(vars) <- c('mid_year', 'isolate_source', 'resistance_breakpoints', 'ast', 'control_strain', 'QA', 'child_adult')


response <- cbind(failures   = MDR_typhi_s_asia$sample_size - MDR_typhi_s_asia$number_resistant, 
                  successes = MDR_typhi_s_asia$number_resistant)

#fit cross validated lasso to select lambda
cv_lasso = cv.glmnet(x = vars , y= response, family = 'binomial', alpha = 1, weights = MDR_typhi_s_asia$sample_size, nfolds = 5, foldid = MDR_typhi_s_asia$fold_id)

cv_lasso$lambda.1se

coef(cv_lasso, s = 0.05)

#2. Add ISO and mid_year
vars <- as.matrix(MDR_typhi_s_asia[, .(iso3, mid_year, isolate_source, resistance_breakpoints, ast, control_strain, QA, child_adult)], with = F)
colnames(vars) <- c('iso3', 'mid_year', 'isolate_source', 'resistance_breakpoints', 'ast', 'control_strain', 'QA', 'child_adult')


response <- cbind(failures   = MDR_typhi_s_asia$sample_size - MDR_typhi_s_asia$number_resistant, 
                  successes = MDR_typhi_s_asia$number_resistant)

#fit cross validated lasso to select lambda
cv_lasso = cv.glmnet(x = vars , y= response, family = 'binomial', alpha = 1, weights = MDR_typhi_s_asia$sample_size, nfolds = 5, foldid = MDR_typhi_s_asia$fold_id)

cv_lasso$lambda.1se

coef(cv_lasso, s =0.05)


#repeat for Southeast asia... only
MDR_typhi_se_asia <- MDR_typhi[MDR_typhi$super_region == "Southeast Asia, East Asia & Oceania",]

## shuffle the data into five random folds
MDR_typhi_se_asia <- MDR_typhi_se_asia[sample(nrow(MDR_typhi_se_asia)),]
MDR_typhi_se_asia[,fold_id := cut(seq(1,nrow(MDR_typhi_se_asia)),breaks=5,labels=FALSE)]
MDR_typhi_se_asia[, a_rowid := seq(1:nrow(MDR_typhi_se_asia))]

#1. with vars but no iso3 or mid year
vars <- as.matrix(MDR_typhi_se_asia[, .(isolate_source, resistance_breakpoints, ast, control_strain, QA, child_adult)], with = F)
colnames(vars) <- c('isolate_source', 'resistance_breakpoints', 'ast', 'control_strain', 'QA', 'child_adult')

response <- cbind(failures   = MDR_typhi_se_asia$sample_size - MDR_typhi_se_asia$number_resistant, 
                  successes = MDR_typhi_se_asia$number_resistant)

#fit cross validated lasso to select lambda
cv_lasso = cv.glmnet(x = vars , y= response, family = 'binomial', alpha = 1, weights = MDR_typhi_se_asia$sample_size, nfolds = 5, foldid = MDR_typhi_se_asia$fold_id)

cv_lasso$lambda.1se
coef(cv_lasso, s = 0.05)

#2 Add mid-year
vars <- as.matrix(MDR_typhi_se_asia[, .(mid_year, isolate_source, resistance_breakpoints, ast, control_strain, QA, child_adult)], with = F)
colnames(vars) <- c('mid_year', 'isolate_source', 'resistance_breakpoints', 'ast', 'control_strain', 'QA', 'child_adult')


response <- cbind(failures   = MDR_typhi_se_asia$sample_size - MDR_typhi_se_asia$number_resistant, 
                  successes = MDR_typhi_se_asia$number_resistant)

#fit cross validated lasso to select lambda
cv_lasso = cv.glmnet(x = vars , y= response, family = 'binomial', alpha = 1, weights = MDR_typhi_se_asia$sample_size, nfolds = 5, foldid = MDR_typhi_se_asia$fold_id)

cv_lasso$lambda.1se

coef(cv_lasso, s = 'lambda.1se')

#3. Add ISO and mid_year
vars <- as.matrix(MDR_typhi_se_asia[, .(iso3, mid_year, isolate_source, resistance_breakpoints, ast, control_strain, QA, child_adult)], with = F)
colnames(vars) <- c('iso3', 'mid_year', 'isolate_source', 'resistance_breakpoints', 'ast', 'control_strain', 'QA', 'child_adult')


response <- cbind(failures   = MDR_typhi_se_asia$sample_size - MDR_typhi_se_asia$number_resistant, 
                  successes = MDR_typhi_se_asia$number_resistant)

#fit cross validated lasso to select lambda
cv_lasso = cv.glmnet(x = vars , y= response, family = 'binomial', alpha = 1, weights = MDR_typhi_se_asia$sample_size, nfolds = 5, foldid = MDR_typhi_se_asia$fold_id)

cv_lasso$lambda.1se

coef(cv_lasso, s =0.05)


#repeat for sSA only
MDR_typhi_sSA <- MDR_typhi[MDR_typhi$super_region == "Sub-Saharan Africa",]

## shuffle the data into five random folds
MDR_typhi_sSA <- MDR_typhi_sSA[sample(nrow(MDR_typhi_sSA)),]
MDR_typhi_sSA[,fold_id := cut(seq(1,nrow(MDR_typhi_sSA)),breaks=5,labels=FALSE)]
MDR_typhi_sSA[, a_rowid := seq(1:nrow(MDR_typhi_sSA))]

#1. with vars but no iso3 or mid year
vars <- as.matrix(MDR_typhi_sSA[, .(isolate_source, resistance_breakpoints, ast, control_strain, QA, child_adult)], with = F)
colnames(vars) <- c('isolate_source', 'resistance_breakpoints', 'ast', 'control_strain', 'QA', 'child_adult')

response <- cbind(failures   = MDR_typhi_sSA$sample_size - MDR_typhi_sSA$number_resistant, 
                  successes = MDR_typhi_sSA$number_resistant)

#fit cross validated lasso to select lambda
cv_lasso = cv.glmnet(x = vars , y= response, family = 'binomial', alpha = 1, weights = MDR_typhi_sSA$sample_size, nfolds = 5, foldid = MDR_typhi_sSA$fold_id)

cv_lasso$lambda.1se
coef(cv_lasso, s = 0.1)

#2 Add mid-year
vars <- as.matrix(MDR_typhi_sSA[, .(mid_year, isolate_source, resistance_breakpoints, ast, control_strain, QA, child_adult)], with = F)
colnames(vars) <- c('mid_year', 'isolate_source', 'resistance_breakpoints', 'ast', 'control_strain', 'QA', 'child_adult')


response <- cbind(failures   = MDR_typhi_sSA$sample_size - MDR_typhi_sSA$number_resistant, 
                  successes = MDR_typhi_sSA$number_resistant)

#fit cross validated lasso to select lambda
cv_lasso = cv.glmnet(x = vars , y= response, family = 'binomial', alpha = 1, weights = MDR_typhi_sSA$sample_size, nfolds = 5, foldid = MDR_typhi_sSA$fold_id)

cv_lasso$lambda.1se

coef(cv_lasso, s = 0.05)

#3. Add ISO and mid_year
vars <- as.matrix(MDR_typhi_sSA[, .(iso3, mid_year, isolate_source, resistance_breakpoints, ast, control_strain, QA, child_adult)], with = F)
colnames(vars) <- c('iso3', 'mid_year', 'isolate_source', 'resistance_breakpoints', 'ast', 'control_strain', 'QA', 'child_adult')


response <- cbind(failures   = MDR_typhi_sSA$sample_size - MDR_typhi_sSA$number_resistant, 
                  successes = MDR_typhi_sSA$number_resistant)

#fit cross validated lasso to select lambda
cv_lasso = cv.glmnet(x = vars , y= response, family = 'binomial', alpha = 1, weights = MDR_typhi_sSA$sample_size, nfolds = 5, foldid = MDR_typhi_sSA$fold_id)

cv_lasso$lambda.1se

coef(cv_lasso, s =0.05)

