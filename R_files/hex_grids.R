###
# Author: Karen Holcomb
# Last Updated: June 2, 2023
##
# Make hexagonal grid across the contiguous US and assign to NOAA climate region
# NOAA climate regions: https://www.ncei.noaa.gov/access/monitoring/reference-maps/us-climate-regions
###

library(tidyr)
library(dplyr)
library(sf)
library(tigris)
options(tigris_use_cache = TRUE)

#counties shapefile
cts <- counties(state = unique(fips_codes$state[!fips_codes$state %in% c("AK","HI","PR","AS","GU","MP","UM","VI")]),
                        resolution = "20m", class="sf") %>% st_transform(., "EPSG:5070")

#states shapefile
sts <- states() %>% filter(., !(STUSPS %in% c("AK","HI","PR","AS","GU","MP","UM","VI"))) %>% st_transform(., "EPSG:5070")

#label states as climate region
sts$clim.regn <- ifelse(sts$STUSPS %in% c("IA","MI","MN","WI"), "UpperMidwest", 
                            ifelse(sts$STUSPS %in% c("IL","IN","KY","MO","OH","TN","WV"), "OhioValley",
                                   ifelse(sts$STUSPS %in% c("AL","FL","GA","NC","SC","VA"), "Southeast",
                                          ifelse(sts$STUSPS %in% c("MT","NE","ND","SD","WY"), "NRockiesandPlains",
                                                 ifelse(sts$STUSPS %in% c("AR","KS","LA","MS","OK","TX"),"South",
                                                        ifelse(sts$STUSPS %in% c("AZ","CO","NM","UT"), "Southwest",
                                                               ifelse(sts$STUSPS %in% c("ID","OR","WA"), "Northwest",
                                                                      ifelse(sts$STUSPS %in% c("CA","NV"), "West", 
                                                                             "Northeast"))))))))


#make hex grid across US
h <- st_make_grid(cts, cellsize = 2e5, square=FALSE) #cell size is diameter of hex (m)
h_cts <- h[cts] #crop to US outline

#which climate region each hex associated with
library(parallel) #if not using parallel, can use lapply instead of mclapply and remove the 'mc.cores' argument
h2 <- mclapply(1:length(h_cts), function(x) { #for each hex, for computation
  print(x) #monitor progress
  w_sts <- st_intersection(sts, h_cts[x]) %>%
    mutate(
      a_ovr = st_area(.), #area of overlap with state
      hnum = x) #which hex this intersection for

  w_sts %>% group_by(clim.regn) %>%
    summarise(tot = sum(a_ovr)) %>% slice_max(tot) %>%
    as.data.frame() %>% select(clim.regn) %>% #which clim region with largest area overlap (across all overlapping states in that region)
    st_set_geometry(., h_cts[x]) #replace geometry with hex geometry

}, mc.cores = 30)

h_clim <- do.call(rbind, h2)
h_clim$fips = 1:nrow(h_clim) #unique ID per hexagon

saveRDS(h_clim, "../data/h_clim_proj.rds") #R object of hexagonal grid
