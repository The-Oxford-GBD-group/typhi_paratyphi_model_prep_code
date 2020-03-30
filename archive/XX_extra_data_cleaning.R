#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Clean up specified datasets ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rm(list = ls())

master_data <- readRDS('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/clean_data.rds')
mydata <-  master_data

#1 Clean demographics
temp_data <- mydata 
temp_data$lower_age <- gsub("[^0-9\\.]", "", temp_data$lower_age) 
temp_data$lower_age <- round(as.numeric(temp_data$lower_age),0)
temp_data$upper_age <- gsub("[^0-9\\.]", "", temp_data$upper_age) 
temp_data$upper_age <- round(as.numeric(temp_data$upper_age),0)
temp_data$age_range <- paste(temp_data$lower_age, temp_data$upper_age, sep = "-")
temp_data$age_range[temp_data$age_range == 'NA-NA'] <- 'No specified age restrictions'
temp_data$age_range[temp_data$age_range == '0-99'] <- 'No specified age restrictions'
temp_data$mean_age <- gsub("[^0-9\\.]", "", temp_data$mean_age) 
temp_data$mean_age <-  as.numeric(temp_data$mean_age)
temp_data$mean_age <- round(temp_data$mean_age, 1) 
temp_data$child_adult <- 'No specified age restrictions/Adults and children'
temp_data$child_adult[temp_data$upper_age<=18] <- 'Children only'
temp_data$child_adult[temp_data$lower_age>=15] <- 'Adults only'

temp_data$inpatients_outpatients[temp_data$source_id == 4526] <- 'outpatients'
temp_data$inpatients_outpatients[temp_data$inpatients_outpatients == 'community'] <- 'Community'
temp_data$inpatients_outpatients[temp_data$inpatients_outpatients == 'inpatients'] <- 'Inpatients'
temp_data$inpatients_outpatients[temp_data$inpatients_outpatients == 'outpatients'] <- 'Outpatients'
temp_data$inpatients_outpatients[temp_data$inpatients_outpatients == 'not specified'] <- 'Not specified'
temp_data$inpatients_outpatients[temp_data$inpatients_outpatients == 'inpatients & outpatients'] <- 'Inpatients & Outpatients'
temp_data$inpatients_outpatients[temp_data$inpatients_outpatients == 'outpatients & emergency dept'] <- 'Outpatients & Emergency department'

table(temp_data$inpatients_outpatients)
table(temp_data$age_range)
table(temp_data$child_adult)

#remove asterisks
temp_data$percent_male <- gsub("[^0-9\\.]", "", temp_data$percent_male) 

#attach the notes for demographics
notes <- read.csv('Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/demographics/studies_demographics_notes.csv', stringsAsFactors = F)
notes <- notes[c(1, 4)]
temp_data <- merge(temp_data, notes, by = 'source_id', all.x = T)

#get the demographics per row
demographics_per_row <- temp_data[c("row_id",
                                    "age_range" ,            
                                    "mean_age",
                                    "child_adult",
                                    "percent_male",
                                    "demographics_notes")]

# aggregate the demographics to the study level
temp_data$percent_male <- as.numeric(temp_data$percent_male)
temp_data$number_male <- (temp_data$percent_male/100)*temp_data$no_examined 

#some changes to enable correct aggregation
temp_data$child_adult[temp_data$source_id == 465] <- 'No specified age restrictions/Adults and children'
temp_data$mean_age[temp_data$source_id == 4370 & temp_data$site_name == 'Bas-Congo'] <-  16
temp_data$mean_age[temp_data$source_id == 4370 & temp_data$site_name == 'Equateur'] <-  25.15
temp_data$mean_age[temp_data$source_id == 4370 & temp_data$site_name == 'Orientale'] <-  25.15
temp_data$mean_age[temp_data$source_id == 5136] <-  15.38
temp_data$mean_age[temp_data$source_id == 5137] <-  21.79
temp_data$mean_age[temp_data$source_id == 5134] <-  NA
temp_data$inpatients_outpatients[temp_data$source_id == 4370] <- 'Inpatients'
	
aggregated_demographics <- unique(temp_data[c('source_id',
                                              "site_name",
                                              "species",
                                              "mid_year",
                                              'lower_age',
                                               'upper_age',
                                               'number_male',
                                               'no_examined',
                                               "mean_age",
                                               "child_adult",
                                               "inpatients_outpatients")])

aggregated_demographics <-  data.table(aggregated_demographics)

all_demographics <- aggregated_demographics[, .(age_range = paste(min(lower_age), max(upper_age), sep = "-"),
                                                percent_male = round(sum(number_male)/sum(no_examined)*100,0)),
                                            by = c("source_id",
                                                   "site_name",
                                                   "species",
                                                   "inpatients_outpatients",
                                                   "child_adult",
                                                   "mean_age")]

# all_demographics$dups <- duplicated(all_demographics, by = c('source_id', 'species', 'site_name'))
# table(all_demographics$dups)

all_demographics$age_range[all_demographics$age_range == 'NA-NA'] <- 'No specified age restrictions'
all_demographics$age_range[all_demographics$age_range == '0-99'] <- 'No specified age restrictions'

#get rid of the notes which say aggregated by year as done this to the data now
notes$demographics_notes[grep('Demographics for all serotypes and years combined', notes$demographics_notes)] <- '*Demographics for all serotypes combined'
notes$demographics_notes[grep('Age demorgaphics are for all years combined; sex demographics are for all serotypes and years combined', notes$demographics_notes)] <- 
  '*Sex demographics for all serotypes combined'
notes <- notes[!grepl('Demographics for all years combined', notes$demographics_notes),]

all_demographics <- merge(all_demographics, notes, by = 'source_id', all.x = T) 

#remove mean ages which cannot be aggregated
aggregated_demographics$mean_age[aggregated_demographics$source_id == 207] <-  NA
aggregated_demographics$mean_age[aggregated_demographics$source_id == 883] <-  NA
aggregated_demographics$mean_age[aggregated_demographics$source_id == 959] <-  NA
aggregated_demographics$child_adult[aggregated_demographics$source_id == 1095] <-  'No specified age restrictions/Adults and children'
aggregated_demographics$inpatients_outpatients[aggregated_demographics$source_id == 2103] <-  'Not specified'

combined_t_pt <- aggregated_demographics[, .(age_range = paste(min(lower_age), max(upper_age), sep = "-"),
                                             percent_male = round(sum(number_male)/sum(no_examined)*100,0)),
                                         by = c("source_id",
                                                "site_name",
                                                "inpatients_outpatients",
                                                "child_adult",
                                                "mean_age")]


# combined_t_pt$dups <- duplicated(combined_t_pt, by = c('source_id', 'site_name'))
# table(combined_t_pt$dups)

combined_t_pt$age_range[combined_t_pt$age_range == 'NA-NA'] <- 'No specified age restrictions'
combined_t_pt$age_range[combined_t_pt$age_range == '0-99'] <- 'No specified age restrictions'
notes <- notes[!grepl('serotype', notes$demographics_notes),]
combined_t_pt <- merge(combined_t_pt, notes, by = 'source_id', all.x = T) 

#Aggregate so one per study
aggregated_demographics$child_adult[aggregated_demographics$source_id == 2946] <-  'No specified age restrictions/Adults and children'
aggregated_demographics$child_adult[aggregated_demographics$source_id == 5065] <-  'No specified age restrictions/Adults and children'
aggregated_demographics$child_adult[aggregated_demographics$source_id == 4370] <-  'No specified age restrictions/Adults and children'
aggregated_demographics$mean_age[aggregated_demographics$source_id == 4370] <-  NA
aggregated_demographics$child_adult[aggregated_demographics$source_id == 2180] <-  'No specified age restrictions/Adults and children'
aggregated_demographics$mean_age[aggregated_demographics$source_id == 2180] <-  NA

all_demogs <- aggregated_demographics[, .(age_range = paste(min(lower_age), max(upper_age), sep = "-"),
                                             percent_male = round(sum(number_male)/sum(no_examined)*100,0)),
                                         by = c("source_id",
                                                "inpatients_outpatients",
                                                "child_adult",
                                                "mean_age")]

# all_demogs$dups <- duplicated(all_demogs, by = c('source_id'))
# table(all_demogs$dups)

all_demogs$age_range[all_demogs$age_range == 'NA-NA'] <- 'No specified age restrictions'
all_demogs$age_range[all_demogs$age_range == '0-99'] <- 'No specified age restrictions'
all_demogs <- merge(all_demogs, notes, by = 'source_id', all.x = T) 

#save the demographics datasets
write.csv(demographics_per_row, 'Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/demographics/demographics_per_row.csv', row.names = F, na = "")
write.csv(all_demographics, 'Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/demographics/demographics_by_study_site_serotype.csv', row.names = F, na = "")
write.csv(combined_t_pt, 'Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/demographics/demographics_by_study_site.csv', row.names = F, na = "")
write.csv(all_demogs, 'Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/demographics/demographics_by_study.csv', row.names = F, na = "")

rm(demographics_per_row, all_demographics, combined_t_pt, notes, temp_data, all_demogs, aggregated_demographics)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Clean the details of the testing
# 1. Cultures performed
# 2. Isolate_source
# 3. Testing method
# 4. Guidelines
# 5. Control strain
# 6. ISO accreditation
# 7. QA score
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mydata <- data.table(mydata)

mydata$cultures_performed <- gsub("[^0-9\\.]", "", mydata$cultures_performed)
mydata$cultures_performed <- as.numeric(mydata$cultures_performed)
N_cultures <- unique(mydata[,.(source_id, cultures_performed)])

# the study has the aggregated and individual results, keep just the agreggated number
N_cultures$cultures_performed[N_cultures$Source_id==2038 & N_cultures$cultures_performed!=9634] <- 0
N_cultures$cultures_performed[N_cultures$Source_id==2093 & N_cultures$cultures_performed!=6956] <- 0

# all other duplicates are due to number of cultures spread over years or sites so sum these
N_cultures <- N_cultures[, .(N_cultures = sum(cultures_performed)),
                         by = c('source_id')]

#1.
unique(mydata$diagnostic_method)
mydata$isolate_source <- 'Mixed sites'
mydata$isolate_source[mydata$diagnostic_method == "Blood culture"] <- 'Normally sterile sites'
mydata$isolate_source[mydata$diagnostic_method == "Blood and bone marrow culture"] <- 'Normally sterile sites'
mydata$isolate_source[mydata$diagnostic_method == "Blood and CSF culture"] <- 'Normally sterile sites'
mydata$isolate_source[mydata$diagnostic_method == "Bile culture"] <- 'Normally sterile sites'
mydata$isolate_source[mydata$diagnostic_method == "Stool culture"] <- 'Non sterile sites'
mydata$isolate_source[mydata$diagnostic_method == "Culture (unknown)"] <- 'Unknown'
mydata$isolate_source[mydata$diagnostic_method == "Unknown"] <- 'Unknown'
table(mydata$isolate_source)

isolate_source <- unique(mydata, by = c('source_id', 'isolate_source'))
isolate_source$dup <-  duplicated(isolate_source$source_id)

for (i in isolate_source$source_id[isolate_source$dup == TRUE]){
  isolate_source$isolate_source[isolate_source$source_id == i] <-  'Mixed sites'
}
isolate_source <- unique(isolate_source, by = c('source_id', 'isolate_source'))
isolate_source <-  isolate_source[,.(source_id, isolate_source)]
#2.
unique(mydata$amr_test)

mydata$amr_test[mydata$amr_test == 'disk-diffusion'] <- 'Disk-diffusion'
mydata$amr_test[mydata$amr_test == "DIsk-diffusion; VITEK 2"] <- "Disk-diffusion; VITEK 2" 
mydata$amr_test[mydata$amr_test == 'Etest'] <- 'E-test'
mydata$amr_test[mydata$amr_test == 'Plasmid resistance'] <- 'Unknown'
mydata$amr_test[mydata$amr_test == 'Disk-diffusion; unspecified MIC testing'] <- 'Disk-diffusion'

ast <- unique(mydata, by = c('source_id', 'amr_test'))
ast <-  ast[order(ast$amr_test)]

ast <- dcast(ast, source_id~amr_test, value.var = 'amr_test')
names(ast)

ast$`Disk-diffusion`[ast$`Disk-diffusion; E-test`=='Disk-diffusion; E-test'] <- 'Disk-diffusion'
ast$`E-test`[ast$`Disk-diffusion; E-test`=='Disk-diffusion; E-test'] <- 'E-test'
ast$`Disk-diffusion; E-test` <-  NULL

ast$`Disk-diffusion`[ast$`Disk-diffusion; Microdilution`=='Disk-diffusion; Microdilution'] <-  'Disk-diffusion'
ast$`Microdilution`[ast$`Disk-diffusion; Microdilution`=='Disk-diffusion; Microdilution'] <-  'Microdilution'
ast$`Disk-diffusion; Microdilution` <-  NULL

ast$`Disk-diffusion`[ast$`Disk-diffusion; Phoenix 100 automated system`=='Disk-diffusion; Phoenix 100 automated system'] <- 'Disk-diffusion'
ast$`Phoenix 100 automated system`[ast$`Disk-diffusion/Phoenix 100 automated system`=='Disk-diffusion; Phoenix 100 automated system'] <- 'Phoenix 100 automated system'
ast$`Disk-diffusion; Phoenix 100 automated system` <-  NULL

ast$`Disk-diffusion`[ast$`Disk-diffusion; VITEK 2` == 'Disk-diffusion; VITEK 2'] <- 'Disk-diffusion'
ast$`Automated methods`[ast$`Disk-diffusion; VITEK 2` == 'Disk-diffusion; VITEK 2'] <- 'Automated methods'
ast$`Disk-diffusion; VITEK 2` <-  NULL

ast$`Agar dilution`[ast$`Agar dilution; VITEK 2` == 'Agar dilution; VITEK 2'] <- 'Agar dilution'
ast$`VITEK 2`[ast$`Agar dilution; VITEK 2` == 'Agar dilution; VITEK 2'] <- 'VITEK 2'
ast$`Agar dilution; VITEK 2` <-  NULL

ast$Unknown <-  NULL

ast$`Dilution Assay` <- NA
ast$`Dilution Assay`[ast$`Microdilution` == 'Microdilution'] <- 'Microdilution'
ast$`Dilution Assay`[ast$`Agar dilution` == 'Agar dilution'] <- 'Microdilution'
ast$`Agar dilution` <-  NULL
ast$`Microdilution` <- NULL

ast$`Automated methods` <- NA
ast$`Automated methods`[ast$`Phoenix 100 automated system` == 'Phoenix 100 automated system'] <- 'Automated methods'
ast$`Automated methods`[ast$`Rapid ATB test` == 'Rapid ATB test'] <- 'Automated methods'
ast$`Automated methods`[ast$`VITEK 2` == 'VITEK 2'] <- 'Automated methods'
ast$`Phoenix 100 automated system` <- NULL
ast$`Rapid ATB test` <- NULL
ast$`VITEK 2` <- NULL

ast$ast <- paste(ast$`Disk-diffusion`, ast$`E-test`, ast$`Dilution Assay`, ast$`Automated methods`, sep = ' & ')
ast$ast <- gsub(' & NA', '', ast$ast)
ast$ast <- gsub('NA & ', '', ast$ast)

ast <-  ast[,.(source_id, ast)]
ast$ast[ast$ast=='NA'] <- 'Not specified'
# ast[is.na(ast$ast)] <- 'Not specified'
table(ast$ast)

# ast$ast[ast$ast == 'Disk-diffusion & E-test'] <- 'Disk-diffusion & MIC determination'
# ast$ast[ast$ast == 'Disk-diffusion & Automated methods'] <- 'Disk-diffusion & MIC determination'
# ast$ast[ast$ast == 'Disk-diffusion & Microdilution'] <- 'Disk-diffusion & MIC determination'
# ast$ast[ast$ast == 'Disk-diffusion & E-test & Microdilution'] <- 'Disk-diffusion & MIC determination'
# 
# ast$ast[ast$ast == 'Microdilution & Automated methods'] <- 'Multiple MIC determination methods'
# ast$ast[ast$ast == 'E-test & Automated methods'] <- 'Multiple MIC determination methods'
# ast$ast[ast$ast == 'E-test & Microdilution'] <- 'Multiple MIC determination methods'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Resistanbce breakpoints ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
guidelines <- unique(mydata, by = c('source_id', 'resistance_breakpoints'))

#remove imputed lines
guidelines$resistance_breakpoints[guidelines$source_id == 117] <- 'CLSI 2012; for azithromycin R \u2265 16 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 135] <- 'NCCLS 1997; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 137] <- 'Working party, 1991 & Ling et al. 1990; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 150] <- 'CLSI 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 151] <- 'CLSI 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 152] <- 'For AMP, CHL & SXT NCCLS 2000; For CIP & NAL CLSI 2005'
guidelines$resistance_breakpoints[guidelines$source_id == 207] <- 'CLSI 2012; for azithromicin R \u2264 18mm'
guidelines$resistance_breakpoints[guidelines$source_id == 235] <- 'BSAC (no year stated); for azithromicin R > 16 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 244] <- 'CLSI 2013'
guidelines$resistance_breakpoints[guidelines$source_id == 248] <- 'CLSI 2013; for azithromicin R > 16 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 267] <- 'For ciprofloxacin R \u2265 0.125 \u03BCg/ml; not stated for other antimicrobials'
guidelines$resistance_breakpoints[guidelines$source_id == 316] <- 'CLSI 2012; for azithromycin R \u2265 16 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 387] <- 'CLSI 2009'
guidelines$resistance_breakpoints[guidelines$source_id == 405] <- 'CLSI 2006; for ciprofloxacin R \u2265 0.19 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 394] <- 'CLSI 2012'
# guidelines$resistance_breakpoints[guidelines$source_id == 548] <- 'For ciprofloxacin R \u2265 0.125 \u03BCg/ml; not stated for other antimicrobials'
# guidelines$resistance_breakpoints[guidelines$source_id == 558] <- 'For ciprofloxacin R \u2265 0.125 \u03BCg/ml; not stated for other antimicrobials'
# guidelines$resistance_breakpoints[guidelines$source_id == 623] <- 'NCCLS (no year stated); for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 639] <- 'For cefepime R \u2264 16mm; for cefixime R \u2264 14mm'
# guidelines$resistance_breakpoints[guidelines$source_id == 693] <- 'For ciprofloxacin R \u2265 0.125 \u03BCg/ml; not stated for other antimicrobials'
guidelines$resistance_breakpoints[guidelines$source_id == 696] <- 'CLSI 2008'
guidelines$resistance_breakpoints[guidelines$source_id == 712] <- 'CLSI (no year stated); ciprofloxacin interpreted as R \u2265 0.125 \u03BCg/ml as MICs provided'
guidelines$resistance_breakpoints[guidelines$source_id == 732] <- 'For aztreonam R \u2264 27mm; ceftazidime R \u2264 22mm; for ceftriaxone R \u2264 25mm'
guidelines$resistance_breakpoints[guidelines$source_id == 737] <- 'NCCLS 2004; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 781] <- 'For ciprofloxacin R \u2265 0.125 \u03BCg/ml; not stated for other antimicrobials'
guidelines$resistance_breakpoints[guidelines$source_id == 787] <- 'CLSI 2005; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 797] <- 'NCCLS (no year stated); ciprofloxacin interpreted as R \u2265 0.125 \u03BCg/ml as MICs provided'
# guidelines$resistance_breakpoints[guidelines$source_id == 833] <- 'NCCLS 2002; ciprofloxacin interpreted as R \u2265 0.125 \u03BCg/ml as MICs provided'
guidelines$resistance_breakpoints[guidelines$source_id == 849] <- 'NCCLS (no year stated); decreased ciprofloxacin sensitivity as \u2265 0.064 \u03BCg/ml as MICs provided'
guidelines$resistance_breakpoints[guidelines$source_id == 850] <- 'CLSI 2003; for azithromycin R \u2265 32 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 852] <- 'CLSI 2009; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 855] <- 'CLSI 2005; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 891] <- 'CLSI 2012; for azithromycin R > 16 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 922] <- 'For ciprofloxacin R \u2265 0.5 \u03BCg/ml; For nalidixic acid R \u2265 32 \u03BCg/ml; not stated for other antimicrobials'
guidelines$resistance_breakpoints[guidelines$source_id == 923] <- 'CLSI 2013; for azithromycin R >16 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 929] <- 'CLSI 2007; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 934] <- 'CLSI 2004'
guidelines$resistance_breakpoints[guidelines$source_id == 1000] <- 'For ampicillin & norfloxacin R \u2264 17mm; for chloramphenicol R \u2264 18mm; for tetracycline R \u2264 19mm; not stated fpr cefotaxime'
guidelines$resistance_breakpoints[guidelines$source_id == 1062] <- 'CLSI 2005; for ciprofloxacin R \u2265 1 \u03BCg/ml; for gatifloxacin R \u2265 8 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 1064] <- 'NCCLS 2002; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1082] <- 'CLSI 2013; for azithrmoycin R > 16 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1095] <- 'CLSI 2007; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1103] <- 'CLSI 2014; EUCAST 2014 for various antimicrobials'
guidelines$resistance_breakpoints[guidelines$source_id == 1126] <- 'CLSI (no year stated); for ciprofloxacin R \u2265 1 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1137] <- 'CLSI 2009; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1223] <- 'CLSI 2008; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 1241] <- 'CLSI 2006; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 1243] <- 'NCCLS 2000; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 1277] <- 'CLSI 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 1334] <- 'NCCLS (no year stated); for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1363] <- 'NCCLS (no year stated); for ciprofloxacin R \u2265 0.25 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1412] <- 'CLSI 2006; For ciprofloxacin R \u2265 0.125 \u03BCg/ml; For azithromycin R \u2265 32 \u03BCg/ml; not stated for other antimicrobials'
guidelines$resistance_breakpoints[guidelines$source_id == 1424] <- 'CLSI 2011; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 1444] <- 'CLSI 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 1465] <- 'CLSI 2007; for azithrmoycin R \u2265 16 \u03BCg/ml'
# guidelines$resistance_breakpoints[guidelines$source_id == 1466] <- 'CLSI 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 1467] <- 'NCCLS 1999; for ciprofloxacin R > 0.03 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1493] <- 'CLSI (no year stated); for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1494] <- 'EUCAST 2012; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1531] <- 'NCCLS (no year stated); for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1532] <- 'NCCLS 1997; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1533] <- 'NCCLS (no year stated); for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1612] <- 'CLSI 2006; For ciprofloxacin R \u2265 0.125 \u03BCg/ml; For azithromycin R \u2265 32 \u03BCg/ml; not stated for other antimicrobials'
guidelines$resistance_breakpoints[guidelines$source_id == 1694] <- 'NCCLS 2000'
guidelines$resistance_breakpoints[guidelines$source_id == 1698] <- 'CLSI (no year stated); for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1730] <- 'CLSI 2012; for ceftriaxone R >2 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1780] <- 'EUCAST; for azithromycin R >16 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1843] <- 'CLSI (no year stated); for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1861] <- 'NCCLS 1998; for azithromycin R \u2265 8 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1875] <- 'CLSI 2009; for azithromycin R >16 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1887] <- 'NCCLS (no year stated); for ciprofloxacin R \u2265 4 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 1917] <- 'For ciprofloxacin R\u2265 0.125 \u03BCg/ml; for gentamicin R \u2265 8 \u03BCg/ml; not stated for other antimicrobials'
guidelines$resistance_breakpoints[guidelines$source_id == 1995] <- 'NCCLS 2001; for ciprofloxacin R \u2265 4 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 2038] <- 'CLSI 2012; EUCAST 2014 for ciprofloxacin'
guidelines$resistance_breakpoints[guidelines$source_id == 2069] <- 'NCCLS 1999; for ciprofloxacin R \u2265 4 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 2165] <- 'CLSI 2008; for ciprofloxacin R \u2265 4 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 2166] <- 'CLSI 2013; for azithromycin R > 16 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 2334] <- 'CLSI 2009; for ciprofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 2376] <- 'NCCLS (no year stated); for ciprofloxacin R \u2265 4 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 2392] <- 'CLSI 2010; for azithromycin R \u2264 15mm'
guidelines$resistance_breakpoints[guidelines$source_id == 2401] <- 'NCCLS 2002; for ciprofloxacin and ofloxacin R \u2265 0.125 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 2824] <- 'CLSI 2009; CLSI 2012 for fluoroquinolones'
guidelines$resistance_breakpoints[guidelines$source_id == 2967] <- 'BSAC/EUCAST 2007-2016; for ciprofloxacin R \u2265 0.125 \u03BCg/ml; for ceftriaxone R \u2265 2 \u03BCg/ml; for azithromycin R \u2265 32 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 3012] <- 'CLSI 2012; for azithromycin BSAC 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 3113] <- 'CLSI 2008 and BSAC 2014 for various observations'
guidelines$resistance_breakpoints[guidelines$source_id == 3182] <- 'CLSI 2012; for gatifloxacin R > 1\u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 4050] <- 'CLSI 2012; CLSI 2013 for moxifloxacin and sparfloxacin'
guidelines$resistance_breakpoints[guidelines$source_id == 4229] <- 'CLSI 2012; for streptomycin R \u2265 32\u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 4363] <- 'CLSI 2012; for azithromycin BSAC 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 4370] <- 'CLSI (no year stated); for azithromycin BSAC (no year stated)'
guidelines$resistance_breakpoints[guidelines$source_id == 4490] <- 'CLSI (no year stated)'
guidelines$resistance_breakpoints[guidelines$source_id == 4515] <- 'CLSI 2013; for azithromycin EUCAST 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 4597] <- 'CLSI 2015-2016'
guidelines$resistance_breakpoints[guidelines$source_id == 5054] <- 'CLSI 2013; ciprofloxacin data are inconsistent with these guidelines and appears to be R  \u2265 0.5 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 5064] <- 'Pre 2012 CLSI (no year stated); 2012 and later CLSI 2012'
guidelines$resistance_breakpoints[guidelines$source_id == 5079] <- 'CLSI 2016; for ciprofloxacin R \u2265 0.5 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 5110] <- 'CLSI 2012-2014; Pre 2012 interepreded using CLSI 2011; for azithrmoycin R \u2265 16 \u03BCg/ml'
guidelines$resistance_breakpoints[guidelines$source_id == 5124] <- 'CLSI 2012; CLSI (no year stated) for data pre 2012; for azithrmoycin R > 16 \u03BCg/ml'


guidelines <- unique(guidelines, by = c('source_id', 'resistance_breakpoints'))
guidelines <- guidelines[order(guidelines$source_id, guidelines$resistance_breakpoints),]
guidelines$dups <-  duplicated(guidelines, by = c('source_id'))
guidelines$dups <- duplicated(guidelines$source_id)
guidelines <- guidelines[!(guidelines$dups == TRUE & guidelines$resistance_breakpoints == 'not specified'),]
guidelines$dups <- duplicated(guidelines$source_id)
guidelines <- guidelines[!(guidelines$dups == TRUE & guidelines$resistance_breakpoints == 'Other'),]
guidelines <- unique(guidelines, by = c('source_id', 'resistance_breakpoints'))

guidelines$dups <- duplicated(guidelines$source_id)
unique(guidelines$source_id[guidelines$dups == TRUE])

guidelines <- guidelines[,.(source_id, resistance_breakpoints)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Control strain used ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
control_strain <- unique(mydata[,.(source_id, control_strain)])
control_strain$dups <-  duplicated(control_strain$source_id)

control_strain$control_strain[control_strain$source_id==4370] <- 'E. coli ATCC 25922; K pneumoniae 700603; S aureus 29213'
control_strain$control_strain[control_strain$source_id==997] <- 'E coli NCTC 10418; S aureaus NCTC 6571'
control_strain$control_strain[control_strain$source_id==548] <- 'E coli ATCC 25922'
control_strain$control_strain[control_strain$source_id==558] <- 'E coli ATCC 25922'
control_strain$control_strain[control_strain$source_id==2999] <- 'E coli ATCC 25922; E faecalis ATCC 29212'
control_strain$control_strain[control_strain$source_id==2093] <- 'E coli NCTC 10418'

control_strain <- unique(control_strain[,.(source_id, control_strain)])

#~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. ISO accreditation ####
#~~~~~~~~~~~~~~~~~~~~~~~~~#
facilities <-gs_read(ss = gs_title("Prevalence data - Typhi & Paratyphi"), ws = "facility lookup")
accredited <- facilities$facility_name[grep('YES', facilities$`ISO accreditation`)]

mydata$ISO_accredited <-  'No'
mydata$ISO_accredited[mydata$hospital_name1 %in% accredited |
            mydata$hospital_name2 %in% accredited |
            mydata$hospital_name3 %in% accredited |
            mydata$hospital_name4 %in% accredited |
            mydata$hospital_name5 %in% accredited |
            mydata$hospital_name6 %in% accredited |
            mydata$hospital_name7 %in% accredited |
            mydata$hospital_name8 %in% accredited] <- 'Yes'


ISO <- unique(mydata[,.(source_id, ISO_accredited)])
ISO <- ISO[order(-ISO$ISO_accredited),]
ISO$dups <- duplicated(ISO$source_id)
ISO <-  ISO[ISO$dups == FALSE]
ISO$dups <-  NULL
rm(facilities, accredited)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Caclulate quality assessment score         #
# 1 = ISO 15819 accredited lab               #
# 2 = Not accredited but states diagnostic   #
#     method, AST method and guidelines used #
# 3 = Neither 1 nor 2                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mydata$QA <-  NA
mydata$QA[mydata$ISO_accredited == 'Yes'] <- 1
mydata$DM <-  'YES'
mydata$AST <- 'YES'
mydata$guidelines <- 'YES'

mydata$DM[mydata$diagnostic_method == 'Unknown' | mydata$diagnostic_method == 'Culture (Unknown)'] <- 'NO'
mydata$AST[mydata$amr_test == 'Unknown'] <- 'No'
mydata$guidelines[mydata$resistance_breakpoints == 'Not specified'] <- 'NO'

mydata$QA[is.na(mydata$QA) &
            mydata$DM == 'YES' &
            mydata$AST == 'YES' &
            mydata$guidelines == 'YES'] <- 2

mydata$QA[is.na(mydata$QA)] <- 3

QA_info <- unique(mydata[,.(source_id, QA)])
QA_info <- QA_info[order(QA_info$QA)]
QA_info$dups <- duplicated(QA_info$source_id)
QA_info <- QA_info[QA_info$dups==FALSE]
QA_info$dups <-  NULL


#merge the testing info together and save the dataset
testing_info <- merge(isolate_source, ast, by = 'source_id')
testing_info <-  merge(testing_info, N_cultures, by = 'source_id')
testing_info <-  merge(testing_info, guidelines, by = 'source_id')
testing_info <-  merge(testing_info, control_strain, by = 'source_id')
testing_info <-  merge(testing_info, ISO, by = 'source_id')
testing_info <-  merge(testing_info, QA_info, by = 'source_id')
testing_info$resistance_breakpoints[testing_info$resistance_breakpoints == 'not specified'] <-  'Not specified'
testing_info$control_strain[testing_info$control_strain == 'not specified'] <-  'Not specified'
testing_info$isolate_source[testing_info$isolate_source == 'Unknown'] <-  'Not specified'
testing_info$N_cultures[is.na(testing_info$N_cultures)] <-  'Not specified'

write.csv(testing_info, 'Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/testing_info/testing_info_by_study.csv', row.names = F)

#get the testing info per row
testing_info_per_row <- mydata[,.(row_id, isolate_source, amr_test, resistance_breakpoints, control_strain, ISO_accredited, QA)]
write.csv(testing_info_per_row, 'Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/testing_info/testing_info_per_row.csv', row.names = F)


location_info <- mydata[c("row_id",
                          "source_id",
                          "best_geo",
                          "setting",
                          "super_region",
                          "region",
                          "country",
                          "iso3",
                          "final_adm1_name",
                          "final_adm2_name",
                          "site_name",
                          "hospital_name",
                          "final_adm0",
                          "final_adm1",
                          "final_adm2",
                          "final_lat",
                          "final_long")]


resistance_data <- mydata[c("row_id",
                            "source_id",
                            "species",
                            "antimicrobial", 
                            "no_examined",
                            "number_resistant",
                            "percentage_resistant")] 

