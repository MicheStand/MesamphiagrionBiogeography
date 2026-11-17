#02_phyloregion.R
#Create by Miguel Stand, 2026 - mstand20@gmail.com

#Load packages
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(phyloregion)
library(letsR)
library(trans4m8)
library(picante)
library(rnaturalearth)
library(elevatr)
library(ggnewscale)
sf_use_s2(FALSE)

#Load species polygons
maps <- read_sf("../data/")

#Time-calibrated phylogeny
phylo <- read.tree("../data/phylogeny.tre")

#Standardize taxonomy
maps <- maps %>% rename(BINOMIAL = sp) %>% mutate(BINOMIAL = str_replace_all(BINOMIAL, " ", "_"))

#Build a Presence–Absence Matrix (PAM)
#lets.presab() requires a SpatVector object
maps_vect <- terra::vect(maps)

#Overlay a 0.25 grid across the study area and record species occurrences within each grid cell
pam <- lets.presab(maps_vect, resol = 0.25, count = TRUE)

#Convert PAM into a community matrix
occurrences <- pam2points(pam, phylo_names = TRUE, gridid = FALSE)
names(occurrences) <- c("lon", "lat", "species")
community <- points2comm(dat = occurrences,
                         res = 0.5,
                         lon = "lon",
                         lat = "lat",
                         species = "species")
community_matrix <- community$comm_dat

#Match the phylogeny to the community matrix
#Remove species present in the phylogeny but absent from the community matrix
phylo_pruned <- keep.tip(phylo, unique(occurrences$species))

#Estimate Faith's Phylogenetic Diversity (PD)
#Faith's PD measures the total branch length connecting all species occurring within each community. Communities containing more distantly related species have larger PD values
#The first two columns of the PAM contain cell coordinates
#Save them before calculating PD
coordinates <- pam$Presence_and_Absence_Matrix[, 1:2]

#Community matrix containing only species presences
community_pa <- pam$Presence_and_Absence_Matrix[, -(1:2)]

#Ensure that the community matrix and phylogeny contain exactly the same species
matched <- match.phylo.comm(phy = phylo_pruned, comm = community_pa)

#Calculate Faith's Phylogenetic Diversity
pd_results <- pd(samp = matched$comm, tree = matched$phy, include.root = TRUE)

#Combine PD values with geographic coordinates
pd_df <- bind_cols(coordinates, pd_results$PD)
colnames(pd_df) <-c("x","y","pd")
head(pd_df)

#Convert PD values to a raster
pdraster <- rast(pd_df, type = "xyz", crs = crs(pam$Richness_Raster))
names(pdraster) <- "Faith_PD"
#Plot
ggplot() + 
  geom_spatraster(data = pdraster) +
  scale_fill_viridis_c(name = "Faith's PD") +
  coord_sf() +
  labs(title = "Phylogenetic Diversity",
       subtitle = "Faith's PD",
       x = "Longitude",
       y = "Latitude") +
  theme_bw(base_size = 14)

#Calculate phylogenetic beta diversity
#Pairwise phylogenetic turnover among communities
phylo_beta <- phylobeta(community_matrix, phylo_pruned)

#Simpson phylogenetic turnover. This metric represents pure turnover while minimizing richness effects
phylo_turnover <- phylo_beta$phylo.beta.sim

#Evaluate clustering algorithms
#Different hierarchical clustering algorithms preserve the dissimilarity matrix with different accuracy
best_linkage <- select_linkage(phylo_turnover)
print(best_linkage)

#Determine the optimal number of phylogenetic regions (k)
#Here the maximum number of regions is restricted to ten because BioGeoBEARS becomes computationally expensive with larger state spaces
optimal_regions <- optimal_phyloregion(phylo_turnover, method = "average", k = 10)

#Plot explained variance
optim_df <- as.data.frame(optimal_regions$df)
optimal_df <- tibble(k  = optimal_regions$optimal$k,
                     ev = optimal_regions$optimal$ev)

ggplot(optim_df, aes(k, ev)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_segment(data = optimal_df,
               aes(x = k, xend = k, y = 0, yend = ev),
               linewidth = 0.7, linetype = 2, colour = "firebrick") +
  geom_point(data = optimal_df,
             aes(x = k, y = ev), size = 5, shape = 21, fill = "firebrick") +
  labs(title = "Optimal Number of Phylogenetic Regions",
       subtitle = "Explained variance across candidate cluster numbers",
       x = "Number of regions (k)",
       y = "Explained variance") +
  theme_bw(base_size = 14)
#Delimit phylogenetic regions
optimal_regions$optimal$k
phylo_regions <- phyloregion(phylo_turnover, k = 4, method = "average", pol = community$map)

#NMDS visualization
#Visualize relationships among communities in reduced ordination space
plot_NMDS(phylo_regions, cex = 3)
text_NMDS(phylo_regions)

#Prepare regions for downstream analyses
phylo_regions$pol$cluster <- LETTERS[1:4] #the k selected
phylo_regions$pol$COLOURS <- c("#C44E52",
                               "#DDCC77",
                               "#7FA58A",
                               "#4C72B0")
plot(phylo_regions)

#Plot
#Countries
countries <- ne_countries(country = c("Colombia",
                                      "Ecuador",
                                      "Venezuela",
                                      "Peru",
                                      "Brazil",
                                      "Panama"),
                          scale = 10,
                          returnclass = "sf")
coastline <- ne_coastline(scale = 10, returnclass = "sf")
dem <- get_elev_raster(locations = countries, z = 8, clip = "locations")
dem <- rast(dem)
dem <- crop(dem, ext(-79, -69, -1, 11))

#Plot and export as .svg
svg(filename="../out/phyloregions.svg")
ggplot() +
  #DEM
  geom_spatraster(data = dem) +
  scale_fill_gradient(low = "grey95",
                      high = "grey30",
                      na.value = "white",
                      guide = "none") +
  new_scale_fill() +
  #Countries
  geom_sf(data = countries,
          fill = NA,
          colour = "black",
          linewidth = 0.3) +
  #Coastline
  geom_sf(data = coastline,
          colour = "black",
          linewidth = 0.35) +
  #Phyloregions
  geom_spatvector(data = phylo_regions$pol,
                  aes(fill = factor(cluster)),
                  colour = "black",
                  linewidth = 0.5,
                  alpha = 0.6) +
  scale_fill_manual(values = phylo_regions$pol$COLOURS) +
  coord_sf(xlim = c(-79, -69),
           ylim = c(-1, 11),
           expand = FALSE) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(fill = "Phyloregion")
dev.off()

#Calculate mean Faith's PD for each phylogenetic region
#Extract the average Faith's PD within each phyloregion polygon
pd_region <- terra::extract(x = pdraster,
                            y = phylo_regions$pol,
                            fun = mean,
                            na.rm = TRUE,
                            ID = FALSE)

#Attach the regional PD values to the polygon layer
phylo_regions$pol$Faith_PD <- pd_region$Faith_PD

ggplot() +
  geom_spatvector(data = phylo_regions$pol,
                  aes(fill = Faith_PD), colour = "black",
                  linewidth = 0.25) +
  scale_fill_viridis_c(option = "C", name = "Mean Faith's PD") +
  coord_sf() +
  labs(title = "Average Phylogenetic Diversity",
       subtitle = "Mean Faith's PD within each phylogenetic region") +
  theme_bw(base_size = 14)

#Assign each species to a phylogenetic region
region_df <- phylo_regions$region.df
species_regions <- sp2cluster(community_matrix, region_df, letters = TRUE)
head(species_regions)

#Calculate pairwise geographic distances among phylogenetic regions
region_distances <- cluster2dist(phylo_regions, BioGeoBEARS = TRUE)
region_distances
write.table(region_distances, "distances_k4.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#Export BioGeoBEARS geographic ranges
phylip_ranges <- data2phylip(species_regions)
head(phylip_ranges)
write.table(phylip_ranges, "phylip_k4.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Save RData
save.image(file = "02_phyloregion.RData")
#end :D