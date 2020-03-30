#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Create rows for the fluoroquinolone insensitivity #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#create a rows for ciprofloxacin by current breakpoints or nalidixic acid if cipro is unavailable
#if cipro is only available for combined serotypes then included seperate NAR ones

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
cipro$resistance_breakpoints <- "R>=0.125ug/ml"
nal$antimicrobial <- "FQR - nalidixic acid"

cipro.nal <- rbind(cipro, nal)

#clean up options for 2796
cipro.nal$antimicrobial[cipro.nal$source_id == 2796] <- "FQR - cipro & nalidixic acid"

cipro.nal$resistance_breakpoints[cipro.nal$source_id == 2796] <- "Not specified"

#join back onto the main dataset
new.rows <- unique(cipro.nal$row_id)
# new.data <- mydata[which(!(mydata$row_id%in%new.rows)),]   # keep all rows for the moment 
new.data <- mydata
mydata<- rbind(new.data, cipro.nal)

rm(cipro.nal, new.data, new.rows, nalid, cipro.id, nal, cipro)
