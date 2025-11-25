# fusing output csvs of SDMs
library(tidyverse)
ver <- "VE8"
# list files
files <- list.files(paste0("/gpfs/gibbs/pi/jetz/from_loomis/data/results_30x30/birds/",
                           ver,"_30x30"), full.names=T,  recursive=F)
files <- files[grep(".csv", files)]
# loop files, combine
for (i in 1:length(files)){
  file <- read_csv(files[i])
  if(i>1){names(table) = names(file)}
  if(i==1){table <- file}else{table <- rbind(table, file)}
}
names(table)[71] = "evi_importance" # fix one name to remove seasonality for consistency

# check number of species
nrow(table[table$season=="summer",])/3
nrow(table[table$season=="winter",])/3
# save
write_csv(table, paste0("/gpfs/gibbs/pi/jetz/from_loomis/data/results_30x30/birds/sdm_birds_",ver,"_0.csv"))
