#############################################################
### CLEAN DATA AND PREP FOR SYSTEMATIC REVIEW AND MAPPING ###
#############################################################
rm(list = ls())
library(doBy)
library(plyr)
library(googlesheets4)
library(data.table)

setwd("Z:/AMR/Pathogens/typhi_paratyphi/model_prep")

#Load required functions
source("H:/Functions/cbindfill.R")
source("H:/Functions/round2.R")

#Load in data from googlesheets
# master.data <-gs_read(ss = gs_title("Prevalence data - Typhi & Paratyphi"), ws = "data")
master.data <- read_sheet("https://docs.google.com/spreadsheets/d/109tfx8GOnephsbQf8PSrAtpa7vJQTfXIjKZRsFxYwGo/edit#gid=0", sheet = 'data')
mydata <- master.data

#change variable names to lowercase
names(mydata) <- tolower(names(mydata))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Remove isolates and change the inlcusion so that years with <10 isolates part of studies with >10 isolates will be included #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#remove the duplicates isolates
# duplicate.isolates <- gs_read(ss = gs_title("Tracking sheet"), ws = "Duplicated review isolates")
duplicate.isolates <- read_sheet("https://docs.google.com/spreadsheets/d/1cLQ1sNf74_9uW9EIR5sJFD8ypbUDKCePZ3vpjTHJFEA/edit#gid=208085225", sheet = 'Duplicated review isolates')
#surveillance.reports <- gs_read(ss = gs_title("Tracking sheet"), ws = "Surveillance reports")

mydata <- mydata[!(mydata$row_id%in%duplicate.isolates$row_id),]
rm(duplicate.isolates)

#we have been given the datasets so add source numbers for these and remove equating studies
# mydata$source_id[mydata$row_id>=15861 & mydata$row_id <=15892] <- 9999
# mydata$source_id[mydata$row_id>=15893 & mydata$row_id <=15906] <- 9998
# mydata$source_id[mydata$row_id>=15907 & mydata$row_id <=16032] <- 9997
mydata <- mydata[mydata$source_id!=1870,]
mydata <- mydata[mydata$source_id!=433,]

mydata <- mydata[mydata$source_id!=2823,] # another duplicated study

#redefine row_id so is definitly unique
mydata$row_id <- seq(1,length(mydata$row_id),1)

#add on GBD region
# merge on gbd regions (lower level)
regions <- read.csv("Z:/AMR/Misc/GBD_regions/GBD regions lookup.csv", stringsAsFactors = F)
regions <- unique(regions[c('super_region', 'region', 'country', 'iso3')])
mydata <- merge(mydata, regions, by = "country")
rm(regions)

#Clean up the species & drop paratyphi B and C records
mydata$species[which(mydata$species=="Salmonella paratyphi A")] <- "Paratyphi"
mydata$species[which(mydata$species=="Salmonella paratyphi")] <- "Paratyphi"
mydata$species[which(mydata$species=="Salmonella typhi")] <- "Typhi"
mydata$species[which(mydata$species=="Salmonella typhi & paratyphi")] <- "Typhi & Paratyphi"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Impute numbers from percentages if they are missing #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mydata$number_resistant[is.na(mydata$number_resistant) & !(is.na(mydata$percentage_resistant))] <- round((mydata$percentage_resistant[is.na(mydata$number_resistant) & !(is.na(mydata$percentage_resistant))]/100)*mydata$no_examined[is.na(mydata$number_resistant) & !(is.na(mydata$percentage_resistant))],0)
mydata$number_intermediate[is.na(mydata$number_intermediate) & !(is.na(mydata$percentage_intermediate))] <- round((mydata$percentage_intermediate[is.na(mydata$number_intermediate) & !(is.na(mydata$percentage_intermediate))]/100)*mydata$no_examined[is.na(mydata$number_intermediate) & !(is.na(mydata$percentage_intermediate))],0)
mydata$number_decreased_susceptibility[is.na(mydata$number_decreased_susceptibility) & !(is.na(mydata$percent_decreased_susceptibility))] <- round((mydata$percent_decreased_susceptibility[is.na(mydata$number_decreased_susceptibility) & !(is.na(mydata$percent_decreased_susceptibility))]/100)*mydata$no_examined[is.na(mydata$number_decreased_susceptibility) & !(is.na(mydata$percent_decreased_susceptibility))],0)
mydata$number_susceptible[is.na(mydata$number_susceptible) & !(is.na(mydata$percentage_susceptible))] <- round((mydata$percentage_susceptible[is.na(mydata$number_susceptible) & !(is.na(mydata$percentage_susceptible))]/100)*mydata$no_examined[is.na(mydata$number_susceptible) & !(is.na(mydata$percentage_susceptible))],0)

#repeat for the adjusted ciprofloxacin data
mydata$number_resistant_adj[is.na(mydata$number_resistant_adj) & !(is.na(mydata$percentage_resistant_adj))] <- round((mydata$percentage_resistant_adj[is.na(mydata$number_resistant_adj) & !(is.na(mydata$percentage_resistant_adj))]/100)*mydata$no_examined[is.na(mydata$number_resistant_adj) & !(is.na(mydata$percentage_resistant_adj))],0)
mydata$number_intermediate_adj[is.na(mydata$number_intermediate_adj) & !(is.na(mydata$percentage_intermediate_adj))] <- round((mydata$percentage_intermediate_adj[is.na(mydata$number_intermediate_adj) & !(is.na(mydata$percentage_intermediate_adj))]/100)*mydata$no_examined[is.na(mydata$number_intermediate_adj) & !(is.na(mydata$percentage_intermediate_adj))],0)
mydata$number_susceptible_adj[is.na(mydata$number_susceptible_adj) & !(is.na(mydata$percentage_susceptible_adj))] <- round((mydata$percentage_susceptible_adj[is.na(mydata$number_susceptible_adj) & !(is.na(mydata$percentage_susceptible_adj))]/100)*mydata$no_examined[is.na(mydata$number_susceptible_adj) & !(is.na(mydata$percentage_susceptible_adj))],0)
mydata$number_ri[is.na(mydata$number_ri) & !(is.na(mydata$percentage_ri))] <- round((mydata$percentage_ri[is.na(mydata$number_ri) & !(is.na(mydata$percentage_ri))]/100)*mydata$no_examined[is.na(mydata$number_ri) & !(is.na(mydata$percentage_ri))],0)

#remove percentage variables (will recalculate later)
mydata$percentage_resistant <- NULL
mydata$percentage_intermediate <- NULL
mydata$percent_decreased_susceptibility <- NULL
mydata$percentage_susceptible <- NULL
mydata$percentage_resistant_adj <- NULL
mydata$percentage_intermediate_adj <- NULL
mydata$percentage_susceptible_adj <- NULL

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Combined intermediate and resistant data #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mydata$NR <- rowSums(mydata[,c("number_resistant", "number_intermediate", "number_decreased_susceptibility")], na.rm=TRUE)
mydata$NR[is.na(mydata$number_resistant) & is.na(mydata$number_intermediate) & is.na(mydata$number_decreased_susceptibility)] <- NA

#replace missing numbers of R/S based ont he corresponding R/S values
mydata$NR[is.na(mydata$NR)& !(is.na(mydata$number_susceptible))] <- mydata$no_examined[is.na(mydata$NR)& !(is.na(mydata$number_susceptible))]-mydata$number_susceptible[is.na(mydata$NR)& !(is.na(mydata$number_susceptible))]
mydata$number_susceptible[is.na(mydata$number_susceptible) & !(is.na(mydata$NR))] <- mydata$no_examined[is.na(mydata$number_susceptible)&!(is.na(mydata$NR))]-mydata$NR[is.na(mydata$number_susceptible)&!(is.na(mydata$NR))]
mydata$number_resistant <- mydata$NR

#For the adjusted data the RI has been completed in the spreasheet so no corrections are required
mydata$number_intermediate <- NULL
mydata$percentage_intermediate <- NULL
mydata$number_decreased_susceptibility <- NULL
mydata$percentage_decreased_susceptibility <- NULL
mydata$NR <- NULL
mydata$number_intermediate_adj <- NULL
mydata$percentage_intermediate_adj <- NULL
mydata$number_resistant_adj <- NULL
mydata$percentager_resistant_adj <- NULL

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Calculate the FQNS                                  ####
# as ciprofloxacin based on the >=0.125ug/ml MIC /       #
# post 2012 CLSI guidelines or nalidixic acid resistance #        
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
cipro <- mydata[mydata$antimicrobial=="ciprofloxacin",]
cipro <- cipro[!is.na(cipro$number_ri),]
cipro1 <- cipro[cipro$species!="Typhi & Paratyphi",]
cipro2 <- cipro[cipro$species=="Typhi & Paratyphi", ] 

nal <-  mydata[mydata$antimicrobial=="nalidixic acid",]
nalid <- unique(nal$source_id[nal$species!="Typhi & Paratyphi"])

#drop ciprofloxacin records for combined serotypes if nalidixic acid are available for combined
cipro2 <- cipro2[which(!(cipro2$source_id%in%nalid)),]
cipro <- rbind(cipro1, cipro2)
rm(cipro1, cipro2)

#replace the resistance number and percentage with the adjusted ones
cipro$number_resistant <- cipro$number_ri
cipro$number_susceptible <- cipro$number_susceptible_adj

#identify nalidixic acid records for studies without ciprofloxacin records
cipro.id <- unique(cipro$source_id)

#remove 2796 from cipro.id and assign to antimicrobial cipro+nal as had 2 time periods - 1 with cipro results and one with nal 
cipro.id <- cipro.id[cipro.id!=2796]

nal <- nal[which(!(nal$source_id%in%cipro.id)),]
nal <-nal[nal$row_id != 2228,]

cipro$antimicrobial <- "FQR - ciprofloxacin"
cipro$resistance_breakpoints <- "R \u2265 0.125\u03BCg/ml"
nal$antimicrobial <- "FQR - nalidixic acid"

cipro.nal <- rbind(cipro, nal)

#clean up options for 2796
cipro.nal$antimicrobial[cipro.nal$source_id == 2796] <- "FQR - cipro & nalidixic acid"

cipro.nal$resistance_breakpoints[cipro.nal$source_id == 2796] <- "Not specified"

cipro.nal$row_id <- max(mydata$row_id)+1:nrow(cipro.nal)+max(mydata$row_id)
  
#join back onto the main dataset
# new.rows <- unique(cipro.nal$row_id)
# new.data <- mydata[which(!(mydata$row_id%in%new.rows)),]   # keep all rows for the moment 
# new.data <- mydata
mydata<- rbind(mydata, cipro.nal)

rm(cipro.nal, nalid, cipro.id, nal, cipro)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Clean up some aspects messing with the aggregation #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#1. Study 2093 has the individual year and total cipro results using different guidelines, remove combined data
mydata <- mydata[mydata$row_id!=2768,] 

#Study 21967 has year specific BSAC/EUCAST guidelines for each year so but comnbined years in. (only this study has these guidelines) 
mydata$resistance_breakpoints[grep("BSAC/EUCAST", mydata$resistance_breakpoints)] <- "BSAC/EUCAST 2007-2016"

#Study 1241 has a prospective and retrospective section
mydata$study_type[mydata$source_id==1241] <- "Prospective & Retrospective" 

#~~~~~~~~~~~~~~~~~~~#
# Calculate mid year#
#~~~~~~~~~~~~~~~~~~~#
#impute start and end years from publication date where not available and calculate mid-year
mydata$start_year_temp <- mydata$start_year
mydata$end_year_temp <- mydata$end_year

mydata$start_year_temp[is.na(mydata$start_year_temp)] <- mydata$publication_year[is.na(mydata$start_year_temp)]-median(mydata$publication_year-mydata$start_year, na.rm = T) 
mydata$end_year_temp[is.na(mydata$end_year_temp)] <- mydata$publication_year[is.na(mydata$end_year_temp)]-median(mydata$publication_year-mydata$end_year, na.rm = T) 

#mid year (round up if 0.5)
mydata$mid_year <- round2(((mydata$start_year_temp+mydata$end_year_temp)/2),0)

mydata$start_year_temp <- NULL
mydata$end_year_temp <- NULL

#~~~~~~~~~~~~~~~~~~~~~~~~~#
# recalculate percentages #
#~~~~~~~~~~~~~~~~~~~~~~~~~#
mydata$percentage_resistant <- round(mydata$number_resistant/mydata$no_examined*100,2)
mydata$percentage_susceptible <- round(mydata$number_susceptible/mydata$no_examined*100,2)

#Include multi-drug & tetracycline/streptpomycin as MDR (as these are on the same plasmid and generally co-occur)
mydata$antimicrobial[mydata$antimicrobial=="multi-drug & tetracycline"] <- "multi-drug"
mydata$antimicrobial[mydata$antimicrobial=="multi-drug & tetracycline & streptomycin"] <- "multi-drug"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Restrict to inclusion criteria #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

mydata <- mydata[which(mydata$mid_year>=1990),]

mydata <- mydata[mydata$no_examined>=5,] #do this later

mydata <- mydata[mydata$species=="Typhi" |
                   mydata$species=="Paratyphi",]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Create a location name based on hospitals and non hospital locations #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mydata$hospital_name <- paste(mydata$hospital_name1,
                                 mydata$hospital_name2,
                                 mydata$hospital_name3,
                                 mydata$hospital_name4,
                                 mydata$hospital_name5,
                                 mydata$hospital_name6,
                                 mydata$hospital_name7,
                                 mydata$hospital_name8,
                                 mydata$hospital_name9, sep = ';')

mydata$hospital_name <-  gsub(';NA', '', mydata$hospital_name)
mydata$hospital_name[mydata$hospital_name=="NA"] <- NA

mydata$location_name <- mydata$hospital_name
mydata$location_name[is.na(mydata$location_name)] <- mydata$site_name[is.na(mydata$location_name)]

#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Restrict and order data #
#~~~~~~~~~~~~~~~~~~~~~~~~~#

mydata <- mydata[c('row_id',
                   "super_region",
                   "region",
                   "iso3",
                   "country",
                   "site_name",
                   "hospital_name",
                   'location_name',
                   "geometry_type",
                   "setting",
                   "start_year",
                   "end_year",
                   "mid_year",
                   "lower_age",
                   "upper_age",
                   "mean_age",
                   "sd_age",
                   "percent_male",
                   "cultures_performed",
                   "species",
                   "diagnostic_method",
                   "amr_test",
                   "antimicrobial",
                   "no_examined",
                   "number_resistant",
                   "percentage_resistant",
                   "number_susceptible",
                   "percentage_susceptible",
                   "resistance_breakpoints",
                   "control_strain",
                   "inpatients_outpatients",
                   "hospital_name1",
                   "hospital_name2",
                   "hospital_name3",
                   "hospital_name4",
                   "hospital_name5",
                   "hospital_name6",                  
                   "hospital_name7",
                   "hospital_name8",
                   "hospital_name9",
                   "latitude",
                   "longitude",
                   "latitude2",
                   "longitude2",                      
                   "admin_level",
                   "gaul_code",
                   "admin_level2",
                   "gaul_code2",
                   "admin_level3",                    
                   "gaul_code3",
                   "source_id",
                   "study_type")]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Data checks and output the files #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#Data Checks
percent.over.100 <- mydata$source_id[mydata$percentage_resistant+mydata$percentage_susceptible>102|mydata$percentage_resistant+mydata$percentage_susceptible<98]
totalRS.over.tested <- mydata$source_id[mydata$number_resistant+mydata$number_susceptible>mydata$no_examined+2] #have +2 to account for rounding due to non-accurate percentages
totalRS.over.tested <- totalRS.over.tested[!(is.na(totalRS.over.tested))]
percent.over.100 <- percent.over.100[!is.na(percent.over.100)]

data.checks <- data.frame(cbind.fill(percent.over.100, totalRS.over.tested))
colnames(data.checks) <- c("dodgey_percentages", "dodgey_numbers")

# remove studies with known and accepted issues
data.checks$dodgey_percentages[data.checks$dodgey_percentages==689] <- NA
data.checks$dodgey_percentages[data.checks$dodgey_percentages==2848] <- NA
data.checks$dodgey_percentages[data.checks$dodgey_percentages==4519] <- NA
data.checks <- data.checks[!(is.na(data.checks$dodgey_numbers)) & !(is.na(data.checks$dodgey_percentages)),]

if(length(data.checks$dodgey_percentages)>0 | length(data.checks$dodgey_numbers)>0){
  message("Inconsistance numbers or percentages, check output")
  write.csv(data.checks, "clean_data/data_errors.csv", row.names = F, na = "")
}


#~~~~~~~~~~~~~~~~~~~~~#
# Any data exclusions #
#~~~~~~~~~~~~~~~~~~~~~#
mydata <- mydata[which(!(mydata$antimicrobial == 'multi-drug' & mydata$source_id == 2093)),]  # The MDR values in this paper un unrealistically low compared to the AMP/CHL/CoT values
mydata <- mydata[which(!(mydata$antimicrobial == 'multi-drug' & mydata$source_id == 4650)),]  # The MDR values in this paper un unrealistically low compared to the AMP/CHL/CoT values and some are higher
mydata <- mydata[which(mydata$source_id != 4651),]  # outbreak, appears unrepresentative, remove

#change the mid year for this study as has records from 04-04 for 1996-1999 and 04-10 for 2000 meaning there are 2 records with miod year == 2000
# as these studies start in the first half of the mid-year the mid year should = start year not end year.
mydata$mid_year[mydata$source_id == 455 & mydata$antimicrobial == 'multi-drug'] <- mydata$start_year[mydata$source_id == 455 & mydata$antimicrobial == 'multi-drug']

#~~~~~~~~~~~~#
# Save files #
#~~~~~~~~~~~~#

#save as an RDS 
saveRDS(mydata, "clean_data/step1_cleaning.rds")
included_studies <- unique(mydata$source_id)
#add 1870 back into this list as it is included but we know have an extended dataset for it (essentially replacing the details)
included_studies <- c(included_studies , 1870, 433)
write.csv(included_studies, 'clean_data/included_studies.csv', row.names = F)

