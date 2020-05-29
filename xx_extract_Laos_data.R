library(readxl)
library(data.table)
library(plyr)

master <- read_excel("Z:/AMR/Data/Laos/data/LOMWRU_BC_UIMahosot_2000to2019.xlsx")
master$Hospital <- 1

master2 <- read_excel("Z:/AMR/Data/Laos/data/LOMWRU_BC_UIOtherhospital_2009to2016.xlsx")

mydata <- rbind.fill(master, master2)

mydata <- mydata[c('UI_No',
                   'SpecimenDate',
                   "HospAdmissDate",
                   "OrganismGrown",
                   "Age_Year",
                   'Gender',
                   "Hospital",
                   'Ampicillin',
                   "Augmentin",
                   "Cefalothin", 
                   "Ceftazidime",
                   "Ceftriaxone",
                   "Ciprofloxacin",
                   "Chloramphenicol",
                   "Colistin",                  
                   "Doxycycline",          
                   "Erythromycin",           
                   "Gentamicin",           
                   "Cefpodoxime",
                   "Imipenem",                  
                   "Methicillin",
                   "NalidixicAcid",
                   "Ofloxacin",                 
                   "Oxacillin",
                   "Pencillin",
                   "TrimSulpha",                
                   "Vancomycin",          
                   "Tetracycline",
                   "Cefoxitin" )]


mydata <- mydata[mydata$OrganismGrown == 'Salmonella Typhi' |mydata$OrganismGrown == 'Salmonella Paratyphi A',]
mydata$SpecimenDate <- substr(mydata$SpecimenDate, 1, 4)
mydata$HospAdmissDate <- substr(mydata$HospAdmissDate, 1, 4)
mydata$year <- mydata$SpecimenDate
mydata$year[is.na(mydata$year)] <- mydata$HospAdmissDate[is.na(mydata$year)]

mydata <- mydata[!is.na(mydata$year),]
mydata <- data.table(mydata)
mydata$MDR <- NA
mydata$MDR[mydata$Ampicillin>=1 & mydata$Chloramphenicol >=1 & mydata$TrimSulpha >=1] <- 1

mydata$MDR[mydata$Ampicillin>=1 & mydata$Chloramphenicol >=1 & mydata$TrimSulpha ==0] <- 0
mydata$MDR[mydata$Ampicillin>=1 & mydata$Chloramphenicol ==0 & mydata$TrimSulpha >=1] <- 0
mydata$MDR[mydata$Ampicillin==0 & mydata$Chloramphenicol >=1 & mydata$TrimSulpha >=1] <- 0

mydata$MDR[mydata$Ampicillin>=1 & mydata$Chloramphenicol ==0 & mydata$TrimSulpha ==0] <- 0
mydata$MDR[mydata$Ampicillin==0 & mydata$Chloramphenicol >=1 & mydata$TrimSulpha ==0] <- 0
mydata$MDR[mydata$Ampicillin==0 & mydata$Chloramphenicol ==0 & mydata$TrimSulpha >=1] <- 0

mydata$MDR[mydata$Ampicillin==0 & mydata$Chloramphenicol ==0 & mydata$TrimSulpha ==0] <- 0

#reshape
mydata <- melt(mydata, id.vars = c('UI_No','Hospital', 'year', 'OrganismGrown', 'Age_Year', 'Gender'),
               measure.vars = names(mydata)[c(8:29, 31)],
               variable.name = 'antimicrobial',
               value.name = 'resistance')

mydata$resistance[mydata$resistance == 2] <- 1
mydata <- mydata[!is.na(mydata$resistance),]
mydata$Gender[mydata$Gender==2] <- 0 

#aggregate
by_hospital <- mydata[,.(min_age = min(Age_Year, na.rm = T),
                    max_age = max(Age_Year, na.rm = T),
                    mean_age = mean(Age_Year, na.rm = T),
                    sd_age = sd(Age_Year, na.rm = T),
                    median_age = median(Age_Year, na.rm = T),
                    lower_iq = quantile(Age_Year, .25, na.rm = T),
                    upper_iq = quantile(Age_Year, .75, na.rm = T),
                    percent_male = sum(Gender, na.rm = T)/length(Gender[!is.na(Gender)])*100,
                    no_examined = length(resistance),
                    number_resistant = sum(resistance),
                    percentage_resistant= sum(resistance)/length(resistance)*100),
                 by = c('OrganismGrown', 'Hospital', 'antimicrobial', 'year')]

total_laos <- mydata[,.(min_age = min(Age_Year, na.rm = T),
                    max_age = max(Age_Year, na.rm = T),
                    mean_age = mean(Age_Year, na.rm = T),
                    sd_age = sd(Age_Year, na.rm = T),
                    median_age = median(Age_Year, na.rm = T),
                    lower_iq = quantile(Age_Year, .25, na.rm = T),
                    upper_iq = quantile(Age_Year, .75, na.rm = T),
                    percent_male = sum(Gender, na.rm = T)/length(Gender[!is.na(Gender)])*100,
                    no_examined = length(resistance),
                    number_resistant = sum(resistance),
                    percentage_resistant= sum(resistance)/length(resistance)*100),
                 by = c('OrganismGrown', 'antimicrobial', 'year')]

by_hospital$percentage_resistant <- round(by_hospital$percentage_resistant,2)
total_laos$percentage_resistant <- round(total_laos$percentage_resistant,2)


write.csv(by_hospital, 'Z:/AMR/Pathogens/typhi_paratyphi/incoming_data/laos/LOMRU_clean_by_hospital.csv', row.names = F)
write.csv(total_laos, 'Z:/AMR/Pathogens/typhi_paratyphi/incoming_data/laos/LOMRU_clean_all_laos.csv', row.names = F)
