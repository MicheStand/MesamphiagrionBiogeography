#01_polygons.R
#Create by Miguel Stand, 2026 - mstand20@gmail.com

#Load packages
library(sf)
library(dplyr)
library(ggplot2)

#Read occurence data
occ <- read.csv("../data/occurence.csv")

#Convert occ to spatial object
occ_sf <- st_as_sf(occ, coords=c("lon", "lat"), crs=4326)

#Create polygon per species
poligonos_sp <- occ_sf %>%  group_by(sp) %>%
  summarise(geometry=st_union(geometry), .groups="drop") %>%
  mutate(geometry=st_convex_hull(geometry))

#Check the structure
glimpse(poligonos_sp)

#Plot
ggplot() +
  geom_sf(data=poligonos_sp, aes(fill=sp), alpha=0.3) +
  geom_sf(data=occ_sf, aes(color=sp), size=1) +
  theme_classic()

#Export as shapefile
st_write(poligonos_sp, "../data/poligonos_spp.shp", delete_dsn = TRUE)

#Save RData
save.image(file = "01_polygons.RData")
#end :D