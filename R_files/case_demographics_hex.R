###
# Author: Karen Holcomb
# Last Updated: June 2, 2023
##
# Process West Nile virus neuroinvasive disease (WNND) case data and census data from county to hexagon
# County-level data assigned to hexagons by overlapping county centroids
###

library(tidyr)
library(dplyr)
library(sf)
library(tigris)
options(tigris_use_cache = TRUE)

#2005-2021 WNV counts - researchers can request it from ArboNET by emailing dvbid2@cdc.gov
df_all <- read.csv("../data/NeuroWNV_by_county_2005-2021FULL.csv", header=T) %>%
  mutate(fips = sprintf("%05d", fips))

#counties shapefile
cts <- counties(state = unique(fips_codes$state[!fips_codes$state %in% c("AK","HI","PR","AS","GU","MP","UM","VI")]),
                resolution = "20m", class="sf") %>% st_transform(., "EPSG:5070")
cts_cents <- st_centroid(cts) #centroids for counties 

totpop <- read.csv("../data/census_data.csv", header=T) %>% #census data for 2010
  select(geoid, tot_pop, pop65) %>% mutate(geoid = sprintf("%05d", geoid)) %>% 
  left_join(., cts %>% as.data.frame() %>% select(GEOID, ALAND), by = c("geoid" = "GEOID"))

#make hex grid across US
h <- st_make_grid(cts, cellsize = 2e5, square=FALSE) #cell size is diameter of hex (m)
h_cts <- h[cts] #crop to US outline

#which county centroids in which hex
b <- st_intersects(h_cts, cts_cents) #list of which county centroid(s) in which hexagon

# calculate total counts, pop per hex
sums_df <- lapply(seq_along(b), function(x) {
  ovr_locs <- cts_cents$GEOID[b[[x]]] #GEOID/fips for any overlapping county centroids for that hexagon
  if(length(ovr_locs) > 0) {
    wnv = df_all[df_all$fips %in% ovr_locs,] #WNV data for overlap
    pop = totpop[totpop$geoid %in% ovr_locs,] #demographics data for overlap
    wnv %>% group_by(year) %>% 
      summarise(tot_count = sum(count)) %>% #sum of cases from overlapping counties per year
      mutate(totpop = sum(pop$tot_pop), #sum of total population from overlapping counties
             pop65 = sum(pop$pop65), #sum of pop > 65 from overlapping counties
             pdens = sum(pop$tot_pop)/(sum(pop$ALAND)/1e6), #total pop density (/km^2) (tot pop/area of overlapping counties)
             geometry = list(h_cts[[x]])) %>%
      st_as_sf(sf_column_name = "geometry", crs = st_crs(h_cts)) #add sfc polygon, sf object
  } else {
    data.frame(year = 2005:2021, tot_count = NA, totpop = NA, pop65 = NA, pdens = NA) %>% #NA counts if no centroid 
      mutate(geometry = list(h_cts[[x]])) %>%
      st_as_sf(sf_column_name = "geometry", crs = st_crs(h_cts)) #add sfc polygon, sf object
  }
})

sd_all <- do.call(rbind, sums_df) #all hexagons with data together
sd_all$fips <- rep(1:290, each=22) #add "fips" column for matching hex dat from climate and land use (fips == hex number)

sd_all <- sd_all %>% subset(!is.na(tot_count)) #remove hexagons w/o data (ones that no county centroid overlaps - internal or coastal)

saveRDS(sd_all, "../data/WNND_hex_2005-2021.rds") #R object with cases, pop per hex
