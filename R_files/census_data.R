###
# Author: Karen Holcomb
# Last Updated: June 2, 2023
##
# Process census csv to get columns for modeling
# Census data from 2010 from: 
# https://www.census.gov/data/tables/time-series/demo/popest/2010s-counties-detail.html#par_textimage_1383669527
###

census <-read.csv('~/cc-est2018-alldata.csv') #census data
pop <- census[,c("STATE","COUNTY","STNAME","CTYNAME","YEAR","AGEGRP","TOT_POP")]
names(pop) <- tolower(names(pop))
pop$geoid <- paste0(sprintf("%02d", as.numeric(pop$state)), sprintf("%03d", as.numeric(pop$county))) #5 digit FIPS code, needed for merging with WNV case data
pop <- pop[pop$year == 1,] #just 2010 census data and no other year's estimates

oldpop <- aggregate(tot_pop~geoid, pop[pop$agegrp > 13,], "sum") #sum 65+ age groups
names(oldpop)[2] <- "pop65"
totpop <- pop[pop$agegrp == 0,-c(5:6)] #just total pop line (drop 'year','agegrp' col)

totpop <- merge(totpop, oldpop, by="geoid")
write.csv(totpop, "../data/census_data.csv")
