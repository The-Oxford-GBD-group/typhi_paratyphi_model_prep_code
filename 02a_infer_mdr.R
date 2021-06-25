#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Calculate MDR for those without AMR            ####
# If  2/3 of AMP, CHL SXT use the minimum resistant #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(data.table)
library(ggplot2)
rm(list = ls())

setwd("Z:/AMR/Pathogens/typhi_paratyphi/model_prep")

clean_data <- readRDS('clean_data/step2_cleaning.rds')

mydata <- clean_data[clean_data$antimicrobial == 'multi-drug' | clean_data$antimicrobial == 'chloramphenicol' | clean_data$antimicrobial == 'ampicillin' | clean_data$antimicrobial == 'co-trimoxazole',]
mydata$antimicrobial[mydata$antimicrobial == 'multi-drug'] <- 'MDR'

#this study has data based on both CLSI and BSAC guidelines seperately - just use the CLSI ones for this as that subset has more data
mydata <- mydata[!(mydata$source_id == 3113 & mydata$resistance_breakpoints == 'BSAC 2012'),]

mydata <- mydata[,.(row_id,
                    source_id,
                    location_name,
                    geometry_type,
                    serotype,
                    mid_year,
                    antimicrobial,
                    sample_size,
                    percentage_resistant)]
 
mydata <- dcast(mydata, source_id+location_name+geometry_type+mid_year+serotype ~ antimicrobial, value.var = c("sample_size", "percentage_resistant", "row_id"))


# Keep the data which has two of the three resistances for MDR
mydata$keep <- 0
mydata$keep[!is.na(mydata$percentage_resistant_chloramphenicol) &
             !is.na(mydata$percentage_resistant_ampicillin) &
              is.na(mydata$`percentage_resistant_co-trimoxazole`)] <- 1

mydata$keep[is.na(mydata$percentage_resistant_chloramphenicol) &
              !is.na(mydata$percentage_resistant_chloramphenicol) &
              !is.na(mydata$`percentage_resistant_co-trimoxazole`)] <- 1

mydata$keep[!is.na(mydata$percentage_resistant_chloramphenicol) &
              is.na(mydata$percentage_resistant_ampicillin) &
              !is.na(mydata$`percentage_resistant_co-trimoxazole`)] <- 1

mydata$keep[!is.na(mydata$percentage_resistant_chloramphenicol) &
              !is.na(mydata$percentage_resistant_ampicillin) &
              !is.na(mydata$`percentage_resistant_co-trimoxazole`)] <- 1

mydata <- mydata[mydata$keep==1,]

res <- mydata[,.(percentage_resistant_chloramphenicol, percentage_resistant_ampicillin, `percentage_resistant_co-trimoxazole`)]

mydata$inferred_mdr <- apply(res, 1, FUN=min, na.rm = T)

png('infered_MDR_check.png', height = 10, width = 10, units = 'cm', res = 300)
ggplot(mydata)+
  geom_point(aes(x = percentage_resistant_MDR, y = inferred_mdr))+
  xlim(0,100)+
  ylim(0,100)+
  xlab('True MDR')+
  ylab('Infered MDR')+
  theme_bw()
dev.off()

cor(mydata$percentage_resistant_MDR[!is.na(mydata$percentage_resistant_MDR)], mydata$inferred_mdr[!is.na(mydata$percentage_resistant_MDR)])^2

#check data points where the inferred value is higher than the actual MDR (allow 5% error)
check <- mydata[mydata$inferred_mdr+5<mydata$percentage_resistant_MDR]
check2 <- mydata[mydata$inferred_mdr>mydata$percentage_resistant_MDR+20]

# In theory this works - for how many rows can this help --> 445 records!!! :D
adjust_data <- mydata[is.na(mydata$percentage_resistant_MDR)]

#get the minimum number examined in AMP/CHL/SXT
res <- adjust_data[,.(sample_size_chloramphenicol, sample_size_ampicillin, `sample_size_co-trimoxazole`)]
adjust_data$sample_size <- apply(res, 1, FUN = min, na.rm = T)


#get the study information for the adjusted data
extra_data <- clean_data[clean_data$row_id %in% adjust_data$row_id_ampicillin |
                           clean_data$row_id %in% adjust_data$row_id_chloramphenicol|
                           clean_data$row_id %in% adjust_data$`row_id_co-trimoxazole`,]


adjust_data <- adjust_data[,.(source_id, 
                              location_name,
                              geometry_type,
                              serotype,
                              mid_year,
                              sample_size,
                              number_resistant = round((inferred_mdr/100)*sample_size,2),
                              percentage_resistant = inferred_mdr)]




extra_data$sample_size <-  NULL
extra_data$number_resistant <-  NULL
extra_data$percentage_resistant <-  NULL
extra_data$row_id <- NULL
extra_data$antimicrobial <- 'MDR - calculated' 
extra_data <-  unique(extra_data)

extra_data <-  merge(extra_data, adjust_data, by = c('source_id', 'location_name', 'geometry_type', 'serotype', 'mid_year'), all.x = F, all.y = T)

extra_data$row_id <- max(clean_data$row_id)+1:nrow(extra_data)+max(clean_data$row_id)

clean_data <- rbind(clean_data, extra_data)

saveRDS(clean_data, 'clean_data/full_cleaned_dataset.rds')

