library(sf)
library(dplyr)

site_dat <- read.csv("FieldData/combined_lat_long.csv")
glimpse(site_dat)
date_dat <- read.csv("FieldData/all_site_data.csv") |> select(Site_ID, Date, Round_Number, CH4_avg, CO2_avg)
glimpse(date_dat)

join_dat <- left_join(site_dat, date_dat, by = join_by("Site_ID"))
glimpse(join_dat)

join_dat_sf <- st_as_sf(join_dat, coords = c("Longitude", "Latitude"))
plot(join_dat_sf["CH4_avg"])
st_write(join_dat_sf, "FieldData/NYS_GHG_Locs_Vect/nys_ghg_locs_vect.shp", append = FALSE)


clim_dat <- read.csv("FieldData/ERA5_Climate_Extraction.csv")
clim_dat |> select(where(is.numeric)) |> 
    cor() |> 
    corrplot::corrplot()
