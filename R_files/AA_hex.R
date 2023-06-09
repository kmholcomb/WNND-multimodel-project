###
# Author: Karen Holcomb
# Last Updated: June 5, 2023
##
# Always Absent model as naive comparison model - always predict 0 cases
###

## Read in WNND case data
df_all <- readRDS("../data/WNND_hex_2005-2021.rds") #WNV cases per year per hex

## Always Absent model - always predict no cases
aa <- do.call(rbind.data.frame, lapply(2015:2021, function(y) { #for each prediction year
  
  probs <- sapply(1:length(unique(df_all$fips)), function(z) { #prob assigned to observed case count (0 if obs > 0 cases, 1 if obs = 0 cases)
    ifelse(df_all$tot_count[df_all$year == y & df_all$fips == unique(df_all$fips)[z]] == 0, 1, 0) 
  })
  data.frame(fips = unique(df_all$fips), year = y, score_aa = log(probs))
}))

saveRDS(aa %>% mutate(fips = as.numeric(as.character(fips))), "../data/AA_ls.rds") #save as R object so do not need to re-run
