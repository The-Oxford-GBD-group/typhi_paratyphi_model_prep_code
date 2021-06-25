rm(list = ls())
library(data.table)
library(Hmisc)

mdr <- readRDS("Z:/AMR/Pathogens/typhi_paratyphi/model_prep/clean_data/MDR_typhi.rds")
mdr <- data.table(mdr)
# mdr <- mdr[mdr$antimicrobial == 'multi-drug',]
mdr <- mdr[mdr$country != 'Malawi Mozambique border',]
mdr <- mdr[mdr$sample_size>10,]

mdr$time_period <- NA
mdr$time_period[mdr$mid_year>=1990 & mdr$mid_year<=1994] <- "1990-1994"
mdr$time_period[mdr$mid_year>=1995 & mdr$mid_year<=1999] <- "1995-1999"
mdr$time_period[mdr$mid_year>=2000 & mdr$mid_year<=2004] <- "2000-2004"
mdr$time_period[mdr$mid_year>=2005 & mdr$mid_year<=2009] <- "2005-2009"
mdr$time_period[mdr$mid_year>=2010 & mdr$mid_year<=2014] <- "2010-2014"
mdr$time_period[mdr$mid_year>=2015 ] <- "2015-2018"

active <- mdr[mdr$setting == 'Community' | mdr$setting == 'Both',]
#aggregate up to study level as multiple records from the same study will bias the data
active <-  active[,.(source_id, country, time_period, sample_size, percentage_resistant)]


#add on one more community study which was just chloramphenicol
extra <- data.frame(4671, 'India', '1995-1999', 63, 29)
names(extra) <- names(active)
active <- rbind(active, extra)
rm(extra)

#calculate uncertainty intervals
active$variance <- (active$percentage_resistant*(100-active$percentage_resistant))/active$sample_size
active$lower_ci <- active$percentage_resistant-(1.96*sqrt(active$variance))
active$lower_ci[active$lower_ci<0] <- 0
active$upper_ci <- active$percentage_resistant+(1.96*sqrt(active$variance))
active$upper_ci[active$upper_ci>100] <- 100
active$variance <- NULL


medians <- mdr[,.(median = median(percentage_resistant),
              IQR_L = quantile(percentage_resistant, .25),
              IQR_U = quantile(percentage_resistant, .75),
              studies = length(row_id)),
           by = c('country', 'time_period')]


#get weighted mean and variance of values (weighted by sample size)
means <- mdr[,.(mean = weighted.mean(percentage_resistant, sample_size),
                mean_lower = weighted.mean(percentage_resistant, sample_size)-(1.96*sqrt(wtd.var(percentage_resistant, sample_size))),
                mean_upper = weighted.mean(percentage_resistant, sample_size)+(1.96*sqrt(wtd.var(percentage_resistant, sample_size)))),
               by = c('country', 'time_period')]

means$mean_lower[means$mean_lower<0] <- 0
means$mean_upper[means$mean_upper>0] <- 100

agg <- merge(means, medians, by = c('country', 'time_period'))

compare_data <- merge(active, agg, by = c('country', 'time_period'))
compare_data$difference <- compare_data$median-compare_data$percentage_resistant

#plot out data
png('Z:/AMR/Pathogens/typhi_paratyphi/bias_adjustment/MDR_Typhi_comparison.png',
    height = 20, width = 20, units = 'cm', res = 300)
ggplot(compare_data)+
  # geom_point(aes(x = median, y = percentage_resistance))
  geom_pointrange(aes(x=median, y = percentage_resistant, ymin = lower_ci, ymax = upper_ci))+
  geom_errorbarh(aes(y = percentage_resistant, xmin = IQR_L, xmax = IQR_U))+
  geom_abline(slope = 1, intercept = 0, colour = 'red')+
  theme_bw()+
  labs(x = 'Passive surveillance, median MDR (%)',
       y = 'Active surveillance, MDR (%)')+
  xlim(0,100)+
  ylim(0,100)
dev.off()

ggplot(compare_data)+
  geom_pointrange(aes(x=mean, y = percentage_resistant, ymin = lower_ci, ymax = upper_ci))+
  # geom_errorbarh(aes(y = percentage_resistant, xmin = mean_lower, xmax = mean_upper))+
  geom_abline(slope = 1, intercept = 0, colour = 'red')+
  theme_bw()+
  labs(x = 'Passive surveillance, mean resistance (%)',
       y = 'Active surveillance, resistance (%)')

#stats tests bewteen the data
t.test(compare_data$percentage_resistant, compare_data$median)

