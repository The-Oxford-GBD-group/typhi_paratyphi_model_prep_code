#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Compare MDR and chloramphenicol data      #
# and look at the potential for a crosswalk #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rm(list = ls())

setwd("K:/AMR/Pathogens/typhi_paratyphi/model_prep")

clean_data <- readRDS('clean_data/clean_data.rds')

mydata <- clean_data[clean_data$antimicrobial == 'multi-drug' | clean_data$antimicrobial == 'chloramphenicol' ,]
mydata$antimicrobial[mydata$antimicrobial == 'multi-drug'] <- 'MDR'

mydata <- mydata[which(mydata$source_id != 455),]
mydata <- mydata[which(mydata$source_id != 3113),]

mydata <- mydata[c("source_id",
                   "iso3",
                   "site_name",
                   'species',
                   "mid_year",
                   'hospital_name1',
                   "antimicrobial",
                   "number_resistant",
                   "percentage_resistant",
                   "row_id")]



mydata <-  data.table(mydata)

mydata <- dcast(mydata, source_id+iso3+site_name+mid_year+species+hospital_name1 ~ antimicrobial, value.var = c("number_resistant", "percentage_resistant", "row_id"))

CHL <- mydata[!is.na(mydata$percentage_resistant_chloramphenicol) & is.na(mydata$percentage_resistant_MDR),]
model.data <- mydata[!is.na(mydata$percentage_resistant_chloramphenicol) & !is.na(mydata$percentage_resistant_MDR),]

scatter.smooth(model.data$percentage_resistant_chloramphenicol, model.data$percentage_resistant_MDR)
cor(model.data$percentage_resistant_chloramphenicol, model.data$percentage_resistant_MDR)
lm(formula= percentage_resistant_MDR ~ percentage_resistant_chloramphenicol, data = model.data)

#r = 0.91
#MDR = -1.537+(0.8475*chloramphenicol)

#Adjust the CHL to MDR and add onto cleaned dataset
CHL$mdr <- -1.537+(0.8479*CHL$percentage_resistant_chloramphenicol)
CHL$mdr[CHL$mdr<0] <- 0
CHL <- CHL[, .(row_id_chloramphenicol, mdr)]
colnames(CHL)[colnames(CHL)=='row_id_chloramphenicol'] <- 'row_id'

CHL <- merge(CHL, clean_data, by = 'row_id', all.x = T, all.y = F)
CHL$antimicrobial <- 'multi-drug-adjusted'
CHL$percentage_resistant <- CHL$mdr
CHL$number_resistant <- round((CHL$mdr/100)*CHL$no_examined,0)
CHL$percentage_susceptible <- 100-CHL$percentage_resistant
CHL$number_susceptible <- CHL$no_examined-CHL$number_resistant
CHL$mdr <- NULL
CHL$row_id <- seq(max(clean_data$row_id+1),max(clean_data$row_id)+length(CHL$row_id),1)

new.data <- rbind(clean_data, CHL)

saveRDS(new.data, 'clean_data/clean_data_crosswalked.rds')

#Repeat for ampicillin to check difference
mydata <- readRDS('clean_data/clean_data.rds')

mydata <- mydata[mydata$antimicrobial == 'multi-drug' | mydata$antimicrobial == 'ampicillin' |mydata$antimicrobial == 'ampicillin',]
mydata$antimicrobial[mydata$antimicrobial == 'multi-drug'] <- 'MDR'

mydata <- mydata[which(mydata$source_id != 455),]
mydata <- mydata[which(mydata$source_id != 3113),]

mydata <- mydata[c("source_id",
                   "iso3",
                   "site_name",
                   'species',
                   "mid_year",
                   'hospital_name1',
                   "antimicrobial",
                   "number_resistant",
                   "percentage_resistant")]

mydata <-  data.table(mydata)

mydata <- dcast(mydata, source_id+iso3+site_name+mid_year+species+hospital_name1 ~ antimicrobial, value.var = c("number_resistant", "percentage_resistant"))

mydata <- mydata[!is.na(mydata$percentage_resistant_ampicillin),]
mydata <- mydata[!is.na(mydata$percentage_resistant_MDR),]

scatter.smooth(mydata$percentage_resistant_ampicillin, mydata$percentage_resistant_MDR)
cor(mydata$percentage_resistant_ampicillin, mydata$percentage_resistant_MDR)
lm(formula= percentage_resistant_MDR ~ percentage_resistant_ampicillin, data = mydata)

#r = 0.808
#MDR = -2.6807+(0.7488*ampicillin)

