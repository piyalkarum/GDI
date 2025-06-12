###################### GDI SCRIPTS ############################
#================== Karunarathne et al. =======================
library(terra)
library(gSoup)

################ GDI calculation and raster creation ###########

# 1. Summary maps of GDI -------
wm<-vect("/maps/TDWG/level3/level3.shp")
eu_range<-c(-10,70,20,85)
eu<-crop(wm,eu_range)

# the 3-D matrix of extrapolated nucleotide diversity per species across European range
all_mats<-readRDS("data/extrapolated_matrix_array_coord_reassigned.rds")
# Calculate the element-wise mean across all matrices
mean_matrix <- apply(all_mats, c(1, 2), mean,na.rm=T)
# Calculate the element-wise sum across all matrices
sum_matrix <- apply(all_mats, c(1, 2), sum,na.rm=T)


# mapping and extrapolating the mean and sum matrices of GDI
# i. mean 
base_raster<-rast(nrow=260,ncol=320,extent=ext(-10 , 70 , 20 , 85),crs="+proj=longlat +datum=WGS84")
mean_values<-(as.vector(t(mean_matrix)))
base_raster<-setValues(base_raster,mean_values)

base_raster<-mask(base_raster,eu)
mt<-matrix(values(base_raster),nrow = nrow(base_raster),ncol = ncol(base_raster),byrow = T)
out_mat<-gSoup:::wind_intpol(mt,wind=10,width = 3,show_progress=T)
base_ras_int<-base_raster
values(base_ras_int)<-NA
values(base_ras_int)<-out_mat
base_ras_int<-mask(base_ras_int,eu)
writeRaster(base_ras_int,"maps/GenDiv_mean_smooth_EU_map_coor_reassigned.tif")

# ii. sum
base_raster_sum<-rast(nrow=260,ncol=320,extent=ext(-10 , 70 , 20 , 85),crs="+proj=longlat +datum=WGS84")
sum_values<-(as.vector(t(sum_matrix)))
base_raster_sum<-setValues(base_raster_sum,sum_values)

base_raster_sum<-mask(base_raster_sum,eu)
mt<-matrix(values(base_raster_sum),nrow = nrow(base_raster_sum),ncol = ncol(base_raster_sum),byrow = T)
out_mat<-gSoup:::wind_intpol(mt,wind=10,width = 3,show_progress=T)
base_ras_sum_int<-base_raster_sum
values(base_ras_sum_int)<-NA
values(base_ras_sum_int)<-out_mat
base_ras_sum_int<-mask(base_ras_sum_int,eu)
writeRaster(base_ras_sum_int,"maps/GenDiv_sum_smooth_EU_map_coord_reassigned.tif")



# 2. Correction Methods ---------------------
# >>>>>> DATA <<<<<<<<<<<<
eu_sp_rich<-rast("maps/species_richness_EU_smooth.tif")
sp_richness<-matrix(values(eu_sp_rich),nrow = nrow(eu_sp_rich),ncol = ncol(eu_sp_rich),byrow = T)
sp_richness[is.na(sp_richness)]<-0
all_mats<-readRDS("data/rasters/extrapolated_matrix_array_coord_reassigned.rds")




# 2.2. wGDI (WEIGHTED GDI) ---------------
eu_sp_rich<-rast("maps/species_richness_EU_smooth.tif")
sp_richness<-matrix(values(eu_sp_rich),nrow = nrow(eu_sp_rich),ncol = ncol(eu_sp_rich),byrow = T)
# matrices to store results
n_sampled <- matrix(0, nrow = 260, ncol = 320)
WMGD <- matrix(NA, nrow = 260, ncol = 320)
SP <- matrix(NA, nrow = 260, ncol = 320)
GDIw <- matrix(NA, nrow = 260, ncol = 320)
alpha<-0.5 # moderate correction

# Loop through each cell (site)
pb<-txtProgressBar(min=1,max=260,style = 3, width=50)
for (i in 1:260) {
  setTxtProgressBar(pb,i)
  for (j in 1:320) {
    diversity_values <- all_mats[i, j, ]
    # Remove NA values to consider only sampled species
    sampled_values <- diversity_values[!is.na(diversity_values)]
    n_sampled[i, j] <- length(sampled_values)
    N <- sp_richness[i, j]
    # Calculate WMGD: mean of sampled genetic diversities
    if (n_sampled[i, j] > 0) {
      WMGD[i, j] <- mean(sampled_values)
    } else {
      WMGD[i, j] <- NA
    }
    
    # Calculate SP: proportion of species sampled
    if(!is.na(N)){
      if (N > 0 ) {
        SP[i, j] <- n_sampled[i, j] / N
      } else {
        SP[i, j] <- NA
      }
    }
    
    # Calculate weighted GDI: WMGD * (SP)^alpha
    if (!is.na(WMGD[i, j]) && !is.na(SP[i, j])) {
      GDIw[i, j] <- WMGD[i, j] * (SP[i, j])^alpha
    } else {
      GDIw[i, j] <- NA
    }
  }
}

# Convert GDIw matrix to a raster for spatial representation
base_raster <- rast(nrows = 260, ncols = 320, extent = ext(-10, 70, 20, 85), crs = "+proj=longlat +datum=WGS84")
GDIw_raster <- setValues(base_raster, as.vector(t(GDIw)))
GDIw_raster<-mask(GDIw_raster,eu)

### apply smoothing to GDIw ------------
GDIw_mat<-matrix(values(GDIw_raster),nrow=nrow(GDIw_raster),ncol=ncol(GDIw_raster),byrow = T)
GDIw_mat[GDIw_mat==0]<-NA
out_mat<-gSoup:::wind_intpol(GDIw_mat,wind=10,width = 3,show_progress=T)
GDIw_int_raster<-GDIw_raster
values(GDIw_int_raster)<-NA
values(GDIw_int_raster)<-(out_mat)
GDIw_int_raster<-mask(GDIw_int_raster,eu)
writeRaster(GDIw_int_raster,"maps/weighted_genDivIndex_Interpolated_raster_jan31.tif")


# 3. rGDI (sample-based RAREFACTION) ------------
library(vegan)
cal_rGDI <- function(nuc_div, species_richness, num_samples = 10, n_reps = 1000) {
  nuc_div <- na.omit(nuc_div)
  if (length(nuc_div) == 0) return(NA)
  rarefaction_limit <- min(num_samples, species_richness, length(nuc_div))
  gdi_values <- replicate(n_reps, {
    sampled_species <- sample(nuc_div, size = rarefaction_limit, replace = FALSE)
    mean(sampled_species, na.rm = TRUE)
  })
  return(mean(gdi_values))
}

n_sampled <- matrix(0, nrow = 260, ncol = 320)
rGDI <- matrix(NA, nrow = 260, ncol = 320)
# Loop through each cell (site)
pb<-txtProgressBar(min=1,max=260,style = 3, width=50)
for (i in 1:260) {
  setTxtProgressBar(pb,i)
  for (j in 1:320) {
    diversity_values <- all_mats[i, j, ]
    # Remove NA values to consider only sampled species
    sampled_values <- diversity_values[!is.na(diversity_values)]
    n_sampled[i, j] <- length(sampled_values)
    N <- sp_richness[i, j]
    # Calculate WMGD: mean of sampled genetic diversities
    if (n_sampled[i, j] > 0 && N > 0) {
      rGDI[i, j] <- cal_rGDI(diversity_values,N)
    } else {
      rGDI[i, j] <- NA
    }
  }
}

## add to the map ----------------
wm<-vect("maps/TDWG/level3/level3.shp")
eu_range<-c(-10,70,20,85)
eu<-crop(wm,eu_range)
# Convert rgdi matrix to a raster for spatial representation
base_raster <- rast(nrows = 260, ncols = 320, extent = ext(-10, 70, 20, 85), crs = "+proj=longlat +datum=WGS84")
rgdi_raster <- setValues(base_raster, as.vector(t(rGDI)))
rgdi_raster<-mask(rgdi_raster,eu)

## apply smoothing to rGDI ------------
rgdi_mat<-matrix(values(rgdi_raster),nrow=nrow(rgdi_raster),ncol=ncol(rgdi_raster),byrow = T)
rgdi_mat[rgdi_mat==0]<-NA
out_mat<-gSoup:::wind_intpol(rgdi_mat,wind=10,width = 3,show_progress=T)
rgdi_int_raster<-rgdi_raster
values(rgdi_int_raster)<-NA
values(rgdi_int_raster)<-(out_mat)
rgdi_int_raster<-mask(rgdi_int_raster,eu)
writeRaster(rgdi_int_raster,"maps/rarefaction_genDivIndex_Interpolated_raster.tif")



# 4. sGDI (REGRESSION based correction) ---------------------------
# Compute mean across layers
n_sampled<-rast("maps/total_sampled_species_per_cell_coord_reassigned.tif")
n_sampled<-matrix(values(n_sampled),nrow = nrow(n_sampled),ncol = ncol(n_sampled),byrow = T)
eu_sp_rich<-rast("maps/species_richness_EU_smooth.tif")
sp_richness<-matrix(values(eu_sp_rich),nrow = nrow(eu_sp_rich),ncol = ncol(eu_sp_rich),byrow = T)
mean_map<-rast("maps/GenDiv_mean_smooth_EU_map_coor_reassigned.tif")
mean_map<-mask(mean_map,eu)
mean_matrix <- matrix(values(mean_map),nrow = nrow(mean_map),ncol = ncol(mean_map),byrow = T)

# Convert matrices to vectors
mean_vec <- as.vector(mean_matrix)
S_sampled <- as.vector(n_sampled)
S_total <- as.vector(sp_richness)
valid_indices <- complete.cases(mean_vec, S_sampled, S_total)
lm_gdi <- lm(mean_vec[valid_indices] ~ I(S_sampled[valid_indices] / (S_total[valid_indices] + 1e-6)))
# Extract residuals
GDI_corrected <- residuals(lm_gdi)

gdi_corrected_mat <- rep(NA,length(mean_vec))
gdi_corrected_mat[valid_indices] <- GDI_corrected
gdi_corrected_mat <- matrix(gdi_corrected_mat, nrow = 260, ncol = 320, byrow = F)

base_raster <- rast(nrows = 260, ncols = 320, extent = ext(-10, 70, 20, 85), crs = "+proj=longlat +datum=WGS84")
gdi_lm_raster <- setValues(base_raster, as.vector(t(gdi_corrected_mat)))
gdi_lm_raster <- mask(gdi_lm_raster, eu)

writeRaster(gdi_lm_raster,"maps/glmCorrected_genDivIndex_Interpolated_raster.tif")

## sampling map [number of species per site]
n_sampled<-readRDS("maps/total_sampled_species_per_cell_coord_reassigned.rda") # from 
base_raster <- rast(nrows = 260, ncols = 320, extent = ext(-10, 70, 20, 85), crs = "+proj=longlat +datum=WGS84")
n_samp_rast <- setValues(base_raster, as.vector(t(n_sampled)))
n_samp_rast <- mask(n_samp_rast, eu)

out_mat<-gSoup:::wind_intpol(n_sampled,wind=10,width = 3,show_progress=T)
base_raster <- rast(nrows = 260, ncols = 320, extent = ext(-10, 70, 20, 85), crs = "+proj=longlat +datum=WGS84")
n_samp_rast <- setValues(base_raster, as.vector(t(out_mat)))
n_samp_rast <- mask(n_samp_rast, eu)
writeRaster(n_samp_rast,"maps/Number_of_sp_sampled_EU.tif")



# PLOT ALL TOGETHER ------------
#+++++++++++++++++++++++++++++++++++++++
library(terra)
library(sf)
library(ggplot2)
library(dplyr)
library(grid)

wm<-vect("maps/TDWG/level3/level3.shp")
eu_range<-c(-10,70,20,85)
eu<-crop(wm,eu_range)
mean_map<-rast("maps/GenDiv_mean_smooth_EU_map_coor_reassigned.tif")
sum_gdi<-rast("maps/GenDiv_sum_smooth_EU_map_coord_reassigned.tif")
lm_cord<-rast("maps/glmCorrected_genDivIndex_Interpolated_raster.tif")
rare_cord<-rast("maps/rarefaction_genDivIndex_Interpolated_raster.tif")
weighted<-rast("maps/weighted_genDivIndex_Interpolated_raster_jan31.tif")
PD_map<-rast("maps/PhyD_map_Interpolated_raster_jan25.tif")
rosids<-rast("maps/GenDiv_mean_Superrosidae_coord_reassigned_smooth_EU_map.tif")
asterids<-rast("maps/GenDiv_mean_Superasteridae_coord_reassigned_smooth_EU_map.tif")
monocots<-rast("maps/GenDiv_mean_Monocot_coord_reassigned_smooth_EU_map.tif")
eu_sp_rich<-rast("maps/species_richness_EU_smooth.tif")
n_sample<-rast("maps/Number_of_sp_sampled_EU.tif")



pdf("plots/diversity_plots3_w_ggplot.pdf",h=5,w=6)
eu_range <- c(-10, 45, 35, 75)
# Load world map
wm <- vect("maps/TDWG/level3/level3.shp")
wm0<-st_as_sf(wm)
wm<-crop(wm,eu_range)
wm_sf <- st_as_sf(wm)
eu_range <- c(xmin = -10, xmax = 45, ymin = 35, ymax = 75)
eu_box <- st_polygon(list(rbind(
  c(eu_range["xmin"], eu_range["ymin"]),
  c(eu_range["xmin"], eu_range["ymax"]),
  c(eu_range["xmax"], eu_range["ymax"]),
  c(eu_range["xmax"], eu_range["ymin"]),
  c(eu_range["xmin"], eu_range["ymin"])
))) |> 
  st_sfc(crs = 4326) |> 
  st_sf()
p0 <- ggplot() +
  geom_sf(data = wm0, fill = "grey", color = NA, linewidth = 0.2) +
  geom_sf(data = eu_box, fill = "blue", color = "blue", alpha = 0.3) +  # <--- The box
  theme_minimal()

print(p0)

## 1. Rarefaction corrected meanGDI --------------
r <- rast("maps/rarefaction_genDivIndex_Interpolated_raster.tif")
r <- crop(r, ext(eu_range))
r <- project(r, "EPSG:4326")
df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
names(df)[3] <- "value"

# Calculate latitudinal density
lat_density <- df %>%
  filter(!is.na(value)) %>%
  group_by(y) %>%
  summarise(total = sum(value, na.rm = TRUE)) %>%
  mutate(norm_total = total / max(total, na.rm = TRUE))
x_base <- max(df$x, na.rm = TRUE)  # base longitude on right side of map
poly_right <- lat_density %>%
  mutate(x = x_base + norm_total * 10)  # control ribbon width
poly_left <- poly_right[rev(seq_len(nrow(poly_right))), ] %>%
  mutate(x = x_base)
lat_poly <- bind_rows(poly_right, poly_left)

cl1<-colorRamps::matlab.like2(10000)
# cl1<-rev(hcl.colors(10000,"Emrld"))
# Save the plot
p <- ggplot() +
  geom_sf(data = wm_sf, fill = "grey", color = NA, linewidth = 0.2) +
  geom_raster(data = df, aes(x = x, y = y, fill = value)) +
  geom_polygon(data = lat_poly, aes(x = x, y = y), fill = "grey20", alpha = 0.5,col=2) +
  scale_fill_gradientn(colors = cl1,na.value = NA,name = "GDI") +
  coord_sf(xlim = c(-15, 45), ylim = c(35, 75), expand = FALSE) +  
  labs(title = "Corrected Mean GDI", x = NULL, y = "Latitude") +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95),legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )

# Print it without clipping
gt <- ggplotGrob(p)
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.newpage()
grid.draw(gt)



## 2. Uncorrected meanGDI -------------------
r <- rast("maps/GenDiv_mean_smooth_EU_map_coor_reassigned.tif")
r <- crop(r, ext(eu_range))
r <- project(r, "EPSG:4326")
df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
names(df)[3] <- "value"

# Calculate latitudinal density
lat_density <- df %>%
  filter(!is.na(value)) %>%
  group_by(y) %>%
  summarise(total = sum(value, na.rm = TRUE)) %>%
  mutate(norm_total = total / max(total, na.rm = TRUE))
x_base <- max(df$x, na.rm = TRUE)  # base longitude on right side of map
poly_right <- lat_density %>%
  mutate(x = x_base + norm_total * 10)  # control ribbon width
poly_left <- poly_right[rev(seq_len(nrow(poly_right))), ] %>%
  mutate(x = x_base)
lat_poly <- bind_rows(poly_right, poly_left)

cl1<-hcl.colors(10000,"Inferno")
# Save the plot
p <- ggplot() +
  geom_sf(data = wm_sf, fill = "grey", color = NA, linewidth = 0.2) +
  geom_raster(data = df, aes(x = x, y = y, fill = value)) +
  geom_polygon(data = lat_poly, aes(x = x, y = y), fill = "grey20", alpha = 0.5,col=2) +
  scale_fill_gradientn(colors = cl1,na.value = NA,name = "GDI") +
  coord_sf(xlim = c(-15, 45), ylim = c(35, 75), expand = FALSE) +  
  labs(title = "Uncorrected Mean GDI", x = NULL, y = "Latitude") +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95),legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )

# Print it without clipping
gt <- ggplotGrob(p)
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.newpage()
grid.draw(gt)


## 3. Sum GDI ----------------------------------
r <- rast("maps/GenDiv_sum_smooth_EU_map_coord_reassigned.tif")
r <- crop(r, ext(eu_range))
r <- project(r, "EPSG:4326")
df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
names(df)[3] <- "value"

# Calculate latitudinal density
lat_density <- df %>%
  filter(!is.na(value)) %>%
  group_by(y) %>%
  summarise(total = sum(value, na.rm = TRUE)) %>%
  mutate(norm_total = total / max(total, na.rm = TRUE))
x_base <- max(df$x, na.rm = TRUE)  # base longitude on right side of map
poly_right <- lat_density %>%
  mutate(x = x_base + norm_total * 10)  # control ribbon width
poly_left <- poly_right[rev(seq_len(nrow(poly_right))), ] %>%
  mutate(x = x_base)
lat_poly <- bind_rows(poly_right, poly_left)

cl1<-(colorspace::terrain_hcl(10000))
# Save the plot
p <- ggplot() +
  geom_sf(data = wm_sf, fill = "grey", color = NA, linewidth = 0.2) +
  geom_raster(data = df, aes(x = x, y = y, fill = value)) +
  geom_polygon(data = lat_poly, aes(x = x, y = y), fill = "grey20", alpha = 0.5,col=2) +
  scale_fill_viridis_c(na.value = NA, name = "GDI") +
  coord_sf(xlim = c(-15, 45), ylim = c(35, 75), expand = FALSE) +  
  labs(title = "Sum GDI", x = NULL, y = "Latitude") +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95),legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )

# Print it without clipping
gt <- ggplotGrob(p)
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.newpage()
grid.draw(gt)

## 4. species richness ----------------------------------
r <- rast("maps/species_richness_EU_smooth.tif")
r <- crop(r, ext(eu_range))
r <- project(r, "EPSG:4326")
df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
names(df)[3] <- "value"

# Calculate latitudinal density
lat_density <- df %>%
  filter(!is.na(value)) %>%
  group_by(y) %>%
  summarise(total = sum(value, na.rm = TRUE)) %>%
  mutate(norm_total = total / max(total, na.rm = TRUE))
x_base <- max(df$x, na.rm = TRUE)  # base longitude on right side of map
poly_right <- lat_density %>%
  mutate(x = x_base + norm_total * 10)  # control ribbon width
poly_left <- poly_right[rev(seq_len(nrow(poly_right))), ] %>%
  mutate(x = x_base)
lat_poly <- bind_rows(poly_right, poly_left)

cl1<-(colorspace::terrain_hcl(10000))
# Save the plot
p <- ggplot() +
  geom_sf(data = wm_sf, fill = "grey", color = NA, linewidth = 0.2) +
  geom_raster(data = df, aes(x = x, y = y, fill = value)) +
  geom_polygon(data = lat_poly, aes(x = x, y = y), fill = "grey20", alpha = 0.5,col=2) +
  scale_fill_viridis_c(na.value = NA, name = "SR") +
  coord_sf(xlim = c(-15, 45), ylim = c(35, 75), expand = FALSE) +  
  labs(title = "Species Richness (no. species)", x = NULL, y = "Latitude") +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95),legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )

# Print it without clipping
gt <- ggplotGrob(p)
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.newpage()
grid.draw(gt)



## 5. Phylogenetic diversity ----------------------------------
r <- rast("maps/PhyD_map_Interpolated_raster_jan25.tif")
r <- crop(r, ext(eu_range))
r <- project(r, "EPSG:4326")
df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
names(df)[3] <- "value"

# Calculate latitudinal density
lat_density <- df %>%
  filter(!is.na(value)) %>%
  group_by(y) %>%
  summarise(total = sum(value, na.rm = TRUE)) %>%
  mutate(norm_total = total / max(total, na.rm = TRUE))
x_base <- max(df$x, na.rm = TRUE)  # base longitude on right side of map
poly_right <- lat_density %>%
  mutate(x = x_base + norm_total * 10)  # control ribbon width
poly_left <- poly_right[rev(seq_len(nrow(poly_right))), ] %>%
  mutate(x = x_base)
lat_poly <- bind_rows(poly_right, poly_left)

cl1<-rev(colorspace::terrain_hcl(10000))
# Save the plot
p <- ggplot() +
  geom_sf(data = wm_sf, fill = "grey", color = NA, linewidth = 0.2) +
  geom_raster(data = df, aes(x = x, y = y, fill = value)) +
  geom_polygon(data = lat_poly, aes(x = x, y = y), fill = "grey20", alpha = 0.5,col=2) +
  scale_fill_gradientn(colors = cl1,na.value = NA,name = "PD") +
  coord_sf(xlim = c(-15, 45), ylim = c(35, 75), expand = FALSE) +  
  labs(title = "Weighted Phylogenetic Diversity", x = NULL, y = "Latitude") +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95),legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )

# Print it without clipping
gt <- ggplotGrob(p)
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.newpage()
grid.draw(gt)

## 6. Sampling (no. species) ----------------------------------
r <- rast("maps/Number_of_sp_sampled_EU.tif")
r <- crop(r, ext(eu_range))
r <- project(r, "EPSG:4326")
df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
names(df)[3] <- "value"

# Calculate latitudinal density
lat_density <- df %>%
  filter(!is.na(value)) %>%
  group_by(y) %>%
  summarise(total = sum(value, na.rm = TRUE)) %>%
  mutate(norm_total = total / max(total, na.rm = TRUE))
x_base <- max(df$x, na.rm = TRUE)  # base longitude on right side of map
poly_right <- lat_density %>%
  mutate(x = x_base + norm_total * 10)  # control ribbon width
poly_left <- poly_right[rev(seq_len(nrow(poly_right))), ] %>%
  mutate(x = x_base)
lat_poly <- bind_rows(poly_right, poly_left)

cl1<-(hcl.colors(10000,"Inferno"))
# Save the plot
p <- ggplot() +
  geom_sf(data = wm_sf, fill = "grey", color = NA, linewidth = 0.2) +
  geom_raster(data = df, aes(x = x, y = y, fill = value)) +
  geom_polygon(data = lat_poly, aes(x = x, y = y), fill = "grey20", alpha = 0.5,col=2) +
  scale_fill_gradientn(colors = cl1,na.value = NA,name = "No. Species") +
  coord_sf(xlim = c(-15, 45), ylim = c(35, 75), expand = FALSE) +  
  labs(title = "Sampled Species", x = NULL, y = "Latitude") +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95),legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )

# Print it without clipping
gt <- ggplotGrob(p)
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.newpage()
grid.draw(gt)


## 7. Sampling (% samp/SR) ----------------------------------
eu_range <- c(-10, 45, 35, 72)
r0<-rast("maps/species_richness_EU_smooth.tif")
r0<-crop(r0,ext(eu_range))
r <- rast("maps/Number_of_sp_sampled_EU.tif")
r <- crop(r, ext(eu_range))
r<-r/r0
r <- project(r, "EPSG:4326")
df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
names(df)[3] <- "value"

# Calculate latitudinal density
lat_density <- df %>%
  filter(!is.na(value)) %>%
  group_by(y) %>%
  summarise(total = sum(value, na.rm = TRUE)) %>%
  mutate(norm_total = total / max(total, na.rm = TRUE))
x_base <- max(df$x, na.rm = TRUE)  # base longitude on right side of map
poly_right <- lat_density %>%
  mutate(x = x_base + norm_total * 10)  # control ribbon width
poly_left <- poly_right[rev(seq_len(nrow(poly_right))), ] %>%
  mutate(x = x_base)
lat_poly <- bind_rows(poly_right, poly_left)

cl1<-(hcl.colors(10000,"Inferno"))
# Save the plot
p <- ggplot() +
  geom_sf(data = wm_sf, fill = "grey", color = NA, linewidth = 0.2) +
  geom_raster(data = df, aes(x = x, y = y, fill = value)) +
  geom_polygon(data = lat_poly, aes(x = x, y = y), fill = "grey20", alpha = 0.5,col=2) +
  scale_fill_gradientn(colors = cl1,na.value = NA,name = "Sampling %") +
  coord_sf(xlim = c(-15, 45), ylim = c(35, 75), expand = FALSE) +  
  labs(title = "Sampled Species % (no. sp/SR)", x = NULL, y = "Latitude") +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95),legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )

# Print it without clipping
gt <- ggplotGrob(p)
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.newpage()
grid.draw(gt)

## 8. Weighted GDI ----------------------------------
eu_range <- c(-10, 45, 35, 75)
r <- rast("maps/weighted_genDivIndex_Interpolated_raster_jan31.tif")
r <- crop(r, ext(eu_range))
r <- project(r, "EPSG:4326")
df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
names(df)[3] <- "value"

# Calculate latitudinal density
lat_density <- df %>%
  filter(!is.na(value)) %>%
  group_by(y) %>%
  summarise(total = sum(value, na.rm = TRUE)) %>%
  mutate(norm_total = total / max(total, na.rm = TRUE))
x_base <- max(df$x, na.rm = TRUE)  # base longitude on right side of map
poly_right <- lat_density %>%
  mutate(x = x_base + norm_total * 10)  # control ribbon width
poly_left <- poly_right[rev(seq_len(nrow(poly_right))), ] %>%
  mutate(x = x_base)
lat_poly <- bind_rows(poly_right, poly_left)

cl1<-rev(hcl.colors(10000,"Emrld"))
# Save the plot
p <- ggplot() +
  geom_sf(data = wm_sf, fill = "grey", color = NA, linewidth = 0.2) +
  geom_raster(data = df, aes(x = x, y = y, fill = value)) +
  geom_polygon(data = lat_poly, aes(x = x, y = y), fill = "grey20", alpha = 0.5,col=2) +
  scale_fill_gradientn(colors = cl1,na.value = NA,name = "GDI") +
  coord_sf(xlim = c(-15, 45), ylim = c(35, 75), expand = FALSE) +  
  labs(title = "Weighted GDI (GDI/SR)", x = NULL, y = "Latitude") +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95),legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )

# Print it without clipping
gt <- ggplotGrob(p)
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.newpage()
grid.draw(gt)


# 9. GDI confidence ------------
eu_range <- c(-10, 45, 35, 75)
r <- rast("maps/EU_confidence_GDI_samp_vs_gdi.tif")
r <- rast("maps/EU_confidence_GDI_samp_vs_gdi_smoothed.tif")
r <- crop(r, ext(eu_range))
r <- project(r, "EPSG:4326")
df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
names(df)[3] <- "value"

# Calculate latitudinal density
lat_density <- df %>%
  filter(!is.na(value)) %>%
  group_by(y) %>%
  summarise(total = sum(value, na.rm = TRUE)) %>%
  mutate(norm_total = total / max(total, na.rm = TRUE))
x_base <- max(df$x, na.rm = TRUE)  # base longitude on right side of map
poly_right <- lat_density %>%
  mutate(x = x_base + norm_total * 10)  # control ribbon width
poly_left <- poly_right[rev(seq_len(nrow(poly_right))), ] %>%
  mutate(x = x_base)
lat_poly <- bind_rows(poly_right, poly_left)

cl1<-rev(hcl.colors(10000,"Emrld"))
# Save the plot
p <- ggplot() +
  geom_sf(data = wm_sf, fill = "grey", color = NA, linewidth = 0.2) +
  geom_raster(data = df, aes(x = x, y = y, fill = value)) +
  geom_polygon(data = lat_poly, aes(x = x, y = y), fill = "grey20", alpha = 0.5,col=2) +
  scale_fill_gradientn(colors = cl1,na.value = NA,name = "Confidence") +
  coord_sf(xlim = c(-15, 45), ylim = c(35, 75), expand = FALSE) +  
  labs(title = "GDI Reliabiligy Confidence", x = NULL, y = "Latitude") +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95),legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )

# Print it without clipping
gt <- ggplotGrob(p)
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.newpage()
grid.draw(gt)


# 10. Rosids -------------
eu_range <- c(-10, 45, 35, 75)
r <- rosids
r <- crop(r, ext(eu_range))
r <- project(r, "EPSG:4326")
df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
names(df)[3] <- "value"

# Calculate latitudinal density
lat_density <- df %>%
  filter(!is.na(value)) %>%
  group_by(y) %>%
  summarise(total = sum(value, na.rm = TRUE)) %>%
  mutate(norm_total = total / max(total, na.rm = TRUE))
x_base <- max(df$x, na.rm = TRUE)  # base longitude on right side of map
poly_right <- lat_density %>%
  mutate(x = x_base + norm_total * 10)  # control ribbon width
poly_left <- poly_right[rev(seq_len(nrow(poly_right))), ] %>%
  mutate(x = x_base)
lat_poly <- bind_rows(poly_right, poly_left)

library(scales)
colors <- c("#FFE4E1", "#F4A6C8", "#D33682")
cl1 <- colorRampPalette(colors)(10000)
# Save the plot
p <- ggplot() +
  geom_sf(data = wm_sf, fill = "grey", color = NA, linewidth = 0.2) +
  geom_raster(data = df, aes(x = x, y = y, fill = value)) +
  geom_polygon(data = lat_poly, aes(x = x, y = y), fill = "grey20", alpha = 0.5,col=2) +
  scale_fill_gradientn(colors = cl1,na.value = NA,name = "GDI") +
  coord_sf(xlim = c(-15, 45), ylim = c(35, 75), expand = FALSE) +  
  labs(title = "Mean GDI Superrosidae", x = NULL, y = "Latitude") +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95),legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )

# Print it without clipping
gt <- ggplotGrob(p)
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.newpage()
grid.draw(gt)


# 11. Asterids -------------
eu_range <- c(-10, 45, 35, 75)
r <- asterids
r <- crop(r, ext(eu_range))
r <- project(r, "EPSG:4326")
df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
names(df)[3] <- "value"

# Calculate latitudinal density
lat_density <- df %>%
  filter(!is.na(value)) %>%
  group_by(y) %>%
  summarise(total = sum(value, na.rm = TRUE)) %>%
  mutate(norm_total = total / max(total, na.rm = TRUE))
x_base <- max(df$x, na.rm = TRUE)  # base longitude on right side of map
poly_right <- lat_density %>%
  mutate(x = x_base + norm_total * 10)  # control ribbon width
poly_left <- poly_right[rev(seq_len(nrow(poly_right))), ] %>%
  mutate(x = x_base)
lat_poly <- bind_rows(poly_right, poly_left)

cl1<-rev(hcl.colors(10000,"Peach"))
# Save the plot
p <- ggplot() +
  geom_sf(data = wm_sf, fill = "grey", color = NA, linewidth = 0.2) +
  geom_raster(data = df, aes(x = x, y = y, fill = value)) +
  geom_polygon(data = lat_poly, aes(x = x, y = y), fill = "grey20", alpha = 0.5,col=2) +
  scale_fill_gradientn(colors = cl1,na.value = NA,name = "GDI") +
  coord_sf(xlim = c(-15, 45), ylim = c(35, 75), expand = FALSE) +  
  labs(title = "Mean GDI Superasteridae", x = NULL, y = "Latitude") +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95),legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )

# Print it without clipping
gt <- ggplotGrob(p)
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.newpage()
grid.draw(gt)



# 12. Monocots -------------
eu_range <- c(-10, 45, 35, 75)
r <- monocots
r <- crop(r, ext(eu_range))
r <- project(r, "EPSG:4326")
df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
names(df)[3] <- "value"

# Calculate latitudinal density
lat_density <- df %>%
  filter(!is.na(value)) %>%
  group_by(y) %>%
  summarise(total = sum(value, na.rm = TRUE)) %>%
  mutate(norm_total = total / max(total, na.rm = TRUE))
x_base <- max(df$x, na.rm = TRUE)  # base longitude on right side of map
poly_right <- lat_density %>%
  mutate(x = x_base + norm_total * 10)  # control ribbon width
poly_left <- poly_right[rev(seq_len(nrow(poly_right))), ] %>%
  mutate(x = x_base)
lat_poly <- bind_rows(poly_right, poly_left)

cl1<-rev(hcl.colors(10000,"Green-Yellow"))
# Save the plot
p <- ggplot() +
  geom_sf(data = wm_sf, fill = "grey", color = NA, linewidth = 0.2) +
  geom_raster(data = df, aes(x = x, y = y, fill = value)) +
  geom_polygon(data = lat_poly, aes(x = x, y = y), fill = "grey20", alpha = 0.5,col=2) +
  scale_fill_gradientn(colors = cl1,na.value = NA,name = "GDI") +
  coord_sf(xlim = c(-15, 45), ylim = c(35, 75), expand = FALSE) +  
  labs(title = "Mean GDI Monocots", x = NULL, y = "Latitude") +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95),legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )

# Print it without clipping
gt <- ggplotGrob(p)
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.newpage()
grid.draw(gt)

dev.off()








