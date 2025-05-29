############################# GDI SCRIPTS ##################################
# =================== Karunarathne et al. ==================================
library(terra)
library(dplyr)
library(sf)
library(ggplot2)
library(scales) 

############################# META ANALYSIS #############################


# 1. plots to show species and markers for TDWG LEVE3_COD -------
pow_new_eu<-read.csv("data/EU_all_native_species_list.csv") # <-- complete European vascular plant list

all_meta_sum<-read.csv("data/marker_assigned_all_meta_data.csv") # <- marker assigned all meta data for European spcies

all_sp_mark_final<-data.table::fread("data/EU_all_species_marker_coordinate_final_list.txt",h=T) # <-- marker assigned geo-coded species list

wm<-vect("maps/TDWG/level3/level3.shp")

eu_range<-c(-10,45,35,85) # Europe range

all_sp_mark_final<-all_sp_mark_final[!is.na(all_sp_mark_final$Longitude),]
points_sf <- st_as_sf(all_sp_mark_final,
                      coords = c("Longitude", "Latitude"),
                      crs = 4326, # WGS 84
                      remove = FALSE)

eu_map_sf <- st_as_sf(eu_map)
sf::sf_use_s2(F) # disable s2 geometry
points_with_region <- st_join(points_sf, eu_map_sf["LEVEL3_COD"], left = FALSE)
summary_df <- points_with_region %>%
  st_drop_geometry() %>%
  group_by(LEVEL3_COD) %>%
  summarise(
    num_species = n_distinct(sp_name),
    num_acc = n(),
    coord_pri = sum(coord == "primary", na.rm = TRUE),
    coord_sec = sum(coord == "secondary", na.rm = TRUE)
  )
eu_map_sf <- eu_map_sf %>%
  left_join(summary_df, by = "LEVEL3_COD")

summary_l3 <- pow_new_eu %>%
  filter(!is.na(area_code_l3)) %>%
  group_by(area_code_l3) %>%
  summarise(num_species_l3 = n_distinct(sp_name)) %>%
  rename(LEVEL3_COD = area_code_l3)

# Join with eu_map
map_l3 <- eu_map_sf %>%
  left_join(summary_l3, by = "LEVEL3_COD")

summary_l2 <- pow_new_eu %>%
  filter(!is.na(region_code_l2)) %>%
  group_by(region_code_l2) %>%
  summarise(num_species_l2 = n_distinct(sp_name)) %>%
  rename(LEVEL2_COD = region_code_l2)

# Join with eu_map
map_l2 <- eu_map_sf %>%
  left_join(summary_l2, by = "LEVEL2_COD")


## adding accession counts ----------
r_int_raster<-rast("maps/all_sp_mark_final_joined_interpolated.tif")
r_int_raster_crop<-crop(r_int_raster,eu_range)
record_df <- as.data.frame(r_int_raster_crop, xy = TRUE, na.rm = TRUE)
names(record_df)[3] <- "record_count"


pdf("plots/meta_analysis_spp_acc_tdwg.pdf",w=8,h=5)
ggplot(data = eu_map_sf) +
  geom_sf(aes(fill = num_species)) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
  theme_minimal() +
  labs(title = "Number of Unique Species by Region", fill = "Species")

ggplot(data = eu_map_sf) +
  geom_sf(aes(fill = num_acc)) +
  scale_fill_viridis_c(option = "viridis", na.value = "grey90") +
  theme_minimal() +
  labs(title = "Number of Accessions by Region", fill = "Accessions")

# Plot: Species richness by LEVEL2_COD
ggplot(map_l2) +
  geom_sf(aes(fill = num_species_l2)) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
  theme_minimal() +
  labs(title = "Species Richness by LEVEL2 (region_code_l2)", fill = "Species")

# Plot: Species richness by LEVEL3_COD
ggplot(map_l3) +
  geom_sf(aes(fill = num_species_l3)) +
  scale_fill_viridis_c(option = "magma", na.value = "grey90") +
  theme_minimal() +
  labs(title = "Species Richness by LEVEL3 (area_code_l3)", fill = "Species")

ggplot() +
  geom_raster(data = record_df, aes(x = x, y = y, fill = record_count)) +
  scale_fill_viridis_c(option = "magma", na.value = "transparent") +
  geom_sf(data = eu_map_sf, fill = NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(-10, 70), ylim = c(20, 85), expand = FALSE) +
  theme_minimal() +
  labs(title = "Record Count per Grid Cell Overlaid on Regions",
       fill = "Record Count",
       x = "Longitude", y = "Latitude")

ggplot() +
  geom_raster(data = record_df, aes(x = x, y = y, fill = record_count)) +
  scale_fill_viridis_c(
    option = "magma",
    trans = "log10",                 # log scale for color mapping
    breaks = c(1, 10, 100, 1000),    # manual breaks on original scale
    labels = comma_format(),         
    na.value = "transparent"
  ) +
  geom_sf(data = eu_map_sf, fill = NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(-10, 70), ylim = c(20, 85), expand = FALSE) +
  theme_minimal() +
  labs(
    title = "Log-Scaled Record Count per Grid Cell",
    fill = "Record Count",
    x = "Longitude", y = "Latitude"
  )
dev.off()





