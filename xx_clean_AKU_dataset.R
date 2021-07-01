library(data.table)
setwd('Z:/AMR/Pathogens/typhi_paratyphi/Data extraction/datasets/Aga Khan')
mydata_2007 <- fread("AGU_2007.csv")
mydata_2007 <- mydata_2007[,.(CITY, ORG, year = 2007, AMP, C, SXT, NAL = V39, CIP, CRO, CFM)]
mydata_2008 <- fread("AGU_2008.csv")
mydata_2008 <- mydata_2008[,.(CITY, ORG, year = 2008, AMP, C, SXT, NAL = V42, CIP, CRO, CFM)]
mydata_2009 <- fread("AGU_2009.csv")
mydata_2009 <- mydata_2009[,.(CITY, ORG, year = 2009, AMP, C, SXT, NAL = V38, CIP, CRO, CFM)]
mydata_2010 <- fread("AGU_2010.csv")
mydata_2010 <- mydata_2010[,.(CITY, ORG, year = 2010, AMP, C, SXT, NAL = V39, CIP, CRO, CFM)]
mydata_2011 <- fread("AGU_2011.csv")
mydata_2011 <- mydata_2011[,.(CITY, ORG, year = 2011, AMP, C, SXT, NAL = V51, CIP, CRO, CFM)]
mydata_2012 <- fread("AGU_2012.csv")
mydata_2012 <- mydata_2012[,.(CITY, ORG, year = 2012, AMP, C, SXT, NAL = V50, CIP, CRO, CFM)]
mydata_2013 <- fread("AGU_2013.csv")
mydata_2013 <- mydata_2013[,.(CITY, ORG, year = 2013, AMP, C, SXT, NAL = V51, CIP, CRO, CFM)]
mydata_2014 <- fread("AGU_2014.csv")
mydata_2014 <- mydata_2014[,.(CITY, ORG, year = 2014, AMP, C, SXT, NAL = V57, CIP, CRO, CFM)]
mydata_2015 <- fread("AGU_2015.csv")
mydata_2015 <- mydata_2015[,.(CITY = City, ORG, year = 2015, AMP, C, SXT, NAL = V51, CIP, CRO, CFM)]
mydata_2019 <- fread("AGU_2019.csv")
mydata_2019 <- mydata_2019[,.(CITY, ORG, year = 2019, AMP=AMPA, C = CA, SXT = SXTA, NAL = NAU, CIP=CIPA, CRO=CROA, CFM = CFMA)]

mydata <- rbindlist(list(mydata_2007, mydata_2008, mydata_2009, mydata_2010, mydata_2011, mydata_2012, 
                    mydata_2013, mydata_2014, mydata_2015, mydata_2019))

rm(mydata_2007, mydata_2008, mydata_2009, mydata_2010, mydata_2011, mydata_2012, 
   mydata_2013, mydata_2014, mydata_2015, mydata_2019)

mydata$MDR <- 0
mydata$MDR[mydata$AMP != 'S' & mydata$C != 'S' & mydata$SXT != 'S'] <- 1
mydata$MDR[is.na(mydata$AMP) |is.na(mydata$C) |is.na(mydata$SXT) ] <- NA
mydata$MDR[mydata$AMP=="" |mydata$C=="" |mydata$SXT=="" ] <- NA

mydata$NAL[mydata$NAL == 'S'] <- '0'
mydata$NAL[mydata$NAL == 'I' | mydata$NAL == 'R'] <- '1'
mydata$NAL <-  as.numeric(mydata$NAL)

mydata$CIP[mydata$CIP == 'S'] <- '0'
mydata$CIP[mydata$CIP == 'I' | mydata$CIP == 'R'] <- '1'
mydata$CIP <-  as.numeric(mydata$CIP)

mydata$CRO[mydata$CRO == 'S'] <- '0'
mydata$CRO[mydata$CRO == 'I' | mydata$CRO == 'R'] <- '1'
mydata$CRO <-  as.numeric(mydata$CRO)

mydata$CFM[mydata$CFM == 'S'] <- '0'
mydata$CFM[mydata$CFM == 'I' | mydata$CFM == 'R'] <- '1'
mydata$CFM <-  as.numeric(mydata$CFM)

mydata$ceph <- 0
mydata$ceph[mydata$CRO == 1 | mydata$CFM == 1] <- 1
mydata$ceph[is.na(mydata$CRO) & is.na(mydata$CFM)] <- NA

mydata <- mydata[mydata$ORG == 'ST' | mydata$ORG == 'SPA',]
mydata$CITY <- tolower(mydata$CITY)

locs <- read.csv('location_lookup.csv')
mydata <- merge(mydata, locs)

MDR_Typhi <- mydata[ORG =='ST',.(number_resistant = sum(MDR, na.rm = T),
                                 sample_size = length(MDR[!is.na(MDR)])),
                    by = c('location_name', 'adj_id', 'year')]

MDR_Typhi <- MDR_Typhi[sample_size>=5 & !is.na(adj_id),.(nid = 5150,
                                                         super_region	= 'South Asia',
                                                         region = 'South Asia',
                                                         country = 'PAK', 
                                                         location_name,
                                                         location_id = 165,
                                                         subnat_location_id = 165,
                                                         adj_id, adj_id_sSA = NA, adj_id_Asia = adj_id,
                                                         year_id = year, 
                                                         age_group_id = 22,
                                                         sex_id = 3,
                                                         measure_id = 18,
                                                         number_resistant, sample_size, 
                                                         val = number_resistant/sample_size, 
                                                         variance = ((number_resistant/sample_size)*(1-(number_resistant/sample_size)))/sample_size, 
                                                         QA = 1)]

write.csv(MDR_Typhi, 'MDR_Typhi.csv', row.names = F)

MDR_Paratyphi <- mydata[ORG =='SPA',.(number_resistant = sum(MDR, na.rm = T),
                                 sample_size = length(MDR[!is.na(MDR)])),
                    by = c('location_name', 'adj_id', 'year')]

MDR_Paratyphi <- MDR_Paratyphi[sample_size>=5 & !is.na(adj_id),.(nid = 5150,
                                                         super_region	= 'South Asia',
                                                         region = 'South Asia',
                                                         country = 'PAK', 
                                                         location_name,
                                                         location_id = 165,
                                                         subnat_location_id = 165,
                                                         adj_id, adj_id_sSA = NA, adj_id_Asia = adj_id,
                                                         year_id = year, 
                                                         age_group_id = 22,
                                                         sex_id = 3,
                                                         measure_id = 18,
                                                         number_resistant, sample_size, 
                                                         val = number_resistant/sample_size, 
                                                         variance = ((number_resistant/sample_size)*(1-(number_resistant/sample_size)))/sample_size, 
                                                         QA = 1)]
write.csv(MDR_Paratyphi, 'MDR_Paratyphi.csv', row.names = F)


FQ_Typhi <- mydata[ORG =='ST',.(number_resistant = sum(CIP, na.rm = T),
                                 sample_size = length(CIP[!is.na(CIP)])),
                    by = c('location_name', 'adj_id', 'year')]

FQ_Typhi <- FQ_Typhi[sample_size>=5 & !is.na(adj_id),.(nid = 5150,
                                                         super_region	= 'South Asia',
                                                         region = 'South Asia',
                                                         country = 'PAK', 
                                                         location_name,
                                                         location_id = 165,
                                                         subnat_location_id = 165,
                                                         adj_id, adj_id_sSA = NA, adj_id_Asia = adj_id,
                                                         year_id = year, 
                                                         age_group_id = 22,
                                                         sex_id = 3,
                                                         measure_id = 18,
                                                         number_resistant, sample_size, 
                                                         val = number_resistant/sample_size, 
                                                         variance = ((number_resistant/sample_size)*(1-(number_resistant/sample_size)))/sample_size, 
                                                         QA = 1)]

write.csv(FQ_Typhi, 'FQ_Typhi.csv', row.names = F)

FQ_Paratyphi <- mydata[ORG =='SPA',.(number_resistant = sum(CIP, na.rm = T),
                                      sample_size = length(CIP[!is.na(CIP)])),
                        by = c('location_name', 'adj_id', 'year')]

FQ_Paratyphi <- FQ_Paratyphi[sample_size>=5 & !is.na(adj_id),.(nid = 5150,
                                                                 super_region	= 'South Asia',
                                                                 region = 'South Asia',
                                                                 country = 'PAK', 
                                                                 location_name,
                                                                 location_id = 165,
                                                                 subnat_location_id = 165,
                                                                 adj_id, adj_id_sSA = NA, adj_id_Asia = adj_id,
                                                                 year_id = year, 
                                                                 age_group_id = 22,
                                                                 sex_id = 3,
                                                                 measure_id = 18,
                                                                 number_resistant, sample_size, 
                                                                 val = number_resistant/sample_size, 
                                                                 variance = ((number_resistant/sample_size)*(1-(number_resistant/sample_size)))/sample_size, 
                                                                 QA = 1)]
write.csv(FQ_Paratyphi, 'FQ_Paratyphi.csv', row.names = F)

ceph_Typhi <- mydata[ORG =='ST',.(number_resistant = sum(ceph, na.rm = T),
                                 sample_size = length(ceph[!is.na(ceph)])),
                    by = c('location_name', 'adj_id', 'year')]

ceph_Typhi <- ceph_Typhi[sample_size>=5 & !is.na(adj_id),.(nid = 5150,
                                                         super_region	= 'South Asia',
                                                         region = 'South Asia',
                                                         country = 'PAK', 
                                                         location_name,
                                                         location_id = 165,
                                                         subnat_location_id = 165,
                                                         adj_id, adj_id_sSA = NA, adj_id_Asia = adj_id,
                                                         year_id = year, 
                                                         age_group_id = 22,
                                                         sex_id = 3,
                                                         measure_id = 18,
                                                         number_resistant, sample_size, 
                                                         val = number_resistant/sample_size, 
                                                         variance = ((number_resistant/sample_size)*(1-(number_resistant/sample_size)))/sample_size, 
                                                         QA = 1)]

write.csv(ceph_Typhi, 'ceph_Typhi.csv', row.names = F)

ceph_Paratyphi <- mydata[ORG =='SPA',.(number_resistant = sum(ceph, na.rm = T),
                                      sample_size = length(ceph[!is.na(ceph)])),
                        by = c('location_name', 'adj_id', 'year')]

ceph_Paratyphi <- ceph_Paratyphi[sample_size>=5 & !is.na(adj_id),.(nid = 5150,
                                                                 super_region	= 'South Asia',
                                                                 region = 'South Asia',
                                                                 country = 'PAK', 
                                                                 location_name,
                                                                 location_id = 165,
                                                                 subnat_location_id = 165,
                                                                 adj_id, adj_id_sSA = NA, adj_id_Asia = adj_id,
                                                                 year_id = year, 
                                                                 age_group_id = 22,
                                                                 sex_id = 3,
                                                                 measure_id = 18,
                                                                 number_resistant, sample_size, 
                                                                 val = number_resistant/sample_size, 
                                                                 variance = ((number_resistant/sample_size)*(1-(number_resistant/sample_size)))/sample_size, 
                                                                 QA = 1)]
write.csv(ceph_Paratyphi, 'ceph_Paratyphi.csv', row.names = F)
