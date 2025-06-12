######## SIMULATION OF SAMPLING VS SPECIES RICHNESS EFFECT ON GDI ##############
#===============================================================================
# This script runs a simulation for varying levels of sampling intensity vs species
# richness to calculate the expected GDI from true GDI. The accuracy is measured by 
# relative squared distance of observed values from the true values. The complete 
# methodology is explained in the Apendices of the manuscript.

library(pbmcapply)
set.seed(123)
n_sites <- 100000
species_richness <- round(rexp(n_sites, rate = 1/250))  
species_richness[species_richness < 1] <- 1
species_richness[species_richness > 3000] <- 3000

# Function to simulate realistic nucleotide diversity values (π)
simulate_pi <- function(N) {
  x <- rnorm(N, mean = 0.1, sd = 0.05)
  x[x < 0] <- 0
  x[x > 0.3] <- 0.3
  return(x)
}
# Create full dataset
sim_data <- do.call(rbind, lapply(1:n_sites, function(site_id) {
  N <- species_richness[site_id]
  data.frame(
    site = site_id,
    species = seq_len(N),
    diversity = simulate_pi(N)
  )
}))
# subsampling
subMeanDGI <- function(data, n_sampled, total_species) {
  if (n_sampled >= total_species) {
    return(mean(data$diversity, na.rm = TRUE))  # Perfect mean
  } else {
    sampled <- sample(data$diversity, n_sampled, replace = FALSE)
    return(mean(sampled, na.rm = TRUE))
  }
}
calc_site_stats <- function(site_id, sim_data, sampling_proportions) {
  site_data <- subset(sim_data, site == site_id)
  N <- nrow(site_data)
  true_mean <- mean(site_data$diversity)
  
  mGDI_values <- sapply(sampling_proportions, function(prop) {
    n_sampled <- floor(prop * N)
    n_sampled <- max(1, min(n_sampled, N))
    subMeanDGI(site_data, n_sampled, N)
  })
  
  data.frame(
    site = site_id,
    sp_richness = N,
    sampling_prop = sampling_proportions,
    mGDI = mGDI_values,
    true_mean = true_mean
  )
}
sampling_proportions <- seq(0.01, 1, by = 0.01)
results2 <- pbmclapply(seq_len(n_sites),
                       FUN = calc_site_stats,
                       sim_data = sim_data,
                       sampling_proportions = sampling_proportions,
                       mc.cores = parallel::detectCores() - 1)
# results_df <- do.call(rbind, results2)
saveRDS(results2, "simulations/GDI_simulation_results_v4.rds", compress = "xz")


########## smooth curve or reliability #############
library(dplyr)
library(ggplot2)
library(scales)

pdf("plots/GDI_sampling_vs_richness_simulation.pdf",h=4,w=5)
results2<-readRDS("simulations/GDI_simulation_results_v4.rds")
results_df<-do.call(rbind,results2)
results_df <- results_df %>%
  mutate(relative_error = abs(mGDI - true_mean) / true_mean)
results_df <- results_df %>%
  mutate(
    is_low_richness = sp_richness <= 2,
    reliable = ifelse(sp_richness > 2 & relative_error < 0.1, TRUE, FALSE)
  )

# filter extreme richness values
results_df_filtered <- results_df %>% filter(sp_richness >= 3, sp_richness <= 3000)

# Group by species richness and sampling proportion
threshold_curve <- results_df_filtered %>%
  group_by(sp_richness, sampling_prop) %>%
  summarise(prop_below_thresh = mean(relative_error < 0.1), .groups = "drop") %>%
  group_by(sp_richness) %>%
  filter(prop_below_thresh >= 0.95) %>%
  slice_min(sampling_prop, with_ties = FALSE)

# Smooth and plot
ggplot(threshold_curve, aes(x = sp_richness, y = sampling_prop)) +
  geom_point(size = 0.6, alpha = 0.6, color = "#1C86EE") +
  geom_smooth(method = "loess", span = 0.3, color = "#E74C3C", linewidth = 1.2) +
  # scale_x_continuous(trans = "log10") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Minimum Sampling Proportion Required for Reliable GDI",
    subtitle = "Curve shows minimum sampling proportion to keep error <10% in 95% of simulations",
    x = "Species Richness",
    y = "Minimum Sampling Proportion"
  ) +
  theme_minimal()

## >>>>>>>> gradient plot <<<<<<<<<<<<<<<
# with ggplot
results_df_dev <- left_join(results_df, threshold_curve, by = "sp_richness", suffix = c("", "_threshold"))
results_df_dev <- results_df_dev %>%
  mutate(
    deviation = sampling_prop - sampling_prop_threshold  # positive = good, negative = bad
  )
col<-rev(c("#1C86EE", "#00BFFF", "#DAA520", "#FF0000"))
ggplot(results_df_dev, aes(x = sp_richness, y = sampling_prop, color = deviation)) +
  geom_point(size = 0.5, alpha = 0.4) +
  scale_color_gradientn(
    colours = col, 
    limits = c(-0.99, 0.99),
    name = "Deviation\nfrom threshold"
  ) +
  scale_x_continuous(limits=c(0,1000)) +
  labs(
    title = "Deviation from Minimum Sampling Threshold",
    subtitle = "Negative values (red) = insufficient sampling; Positive (blue) = reliable estimation",
    x = "Species Richness",
    y = "Sampling Proportion"
  ) +
  theme_minimal()
dev.off()


# B. Predictive model for GDI reliability -----------
library(dplyr)
results_df_filtered <- results_df %>% filter(sp_richness >= 3, sp_richness <= 3000)
# Group by species richness and sampling proportion
threshold_curve <- results_df_filtered %>%
  group_by(sp_richness, sampling_prop) %>%
  summarise(prop_below_thresh = mean(relative_error < 0.1), .groups = "drop") %>%
  group_by(sp_richness) %>%
  filter(prop_below_thresh >= 0.95) %>%
  slice_min(sampling_prop, with_ties = FALSE)
threshold_model <- approxfun(
  x = threshold_curve$sp_richness,
  y = threshold_curve$sampling_prop,
  method = "linear",
  rule = 2  
)
confidence_model <- approxfun(
  x = threshold_curve$sp_richness,
  y = threshold_curve$prop_below_thresh,
  method = "linear",
  rule = 2
)

# Prediction function
gdi_reliability <- function(richness, sampling_prop) {
  threshold <- threshold_model(richness)
  confidence <- confidence_model(richness)
  
  if (sampling_prop >= threshold) {
    status <- "RELIABLE"
    prob <- confidence
  } else {
    status <- "UNRELIABLE"
    drop <- (threshold - sampling_prop) / threshold
    prob <- max(0, confidence * (1 - drop))  # Penalize linearly
  }
  
  return(list(
    status = status,
    threshold = round(threshold, 3),
    confidence = round(prob, 3)
  ))
}

gdi_reliability(100,0.5)

## apply this to out data -------
eu_sp_rich<-rast("maps/species_richness_EU_smooth.tif")
n_sample<-rast("maps/Number_of_sp_sampled_EU.tif")
s_proportion<-n_sample/eu_sp_rich

sprich<-values(eu_sp_rich)
sprich[is.na(sprich)]<-0
prop<-values(s_proportion)
prop[is.na(prop)]<-0

conf_vals<-NULL
for(i in seq_along(sprich)){
  conf_vals[i]<-gdi_reliability(sprich[i],prop[i])$confidence
}
base_raster <- rast(nrows = nrow(eu_sp_rich), ncols = ncol(eu_sp_rich), extent = ext(eu_sp_rich), crs = "+proj=longlat +datum=WGS84")
conf_raster <- setValues(base_raster, as.vector(conf_vals))
conf_raster<-mask(conf_raster,n_sample)
conf_raster<-crop(conf_raster,ext(c(-10, 45, 35, 75)))
writeRaster(conf_raster,"maps/EU_confidence_GDI_samp_vs_gdi.tif", overwrite=T)

## apply a smoothing --------------
conf_mat<-matrix(values(conf_raster),nrow=nrow(conf_raster),ncol=ncol(conf_raster),byrow = T)
conf_mat[conf_mat==0]<-NA
out_mat<-gSoup:::wind_intpol(conf_mat,wind=10,width = 3,show_progress=T)
conf_int_raster<-conf_raster
values(conf_int_raster)<-NA
values(conf_int_raster)<-(out_mat)
conf_int_raster<-mask(conf_int_raster,conf_raster)
writeRaster(conf_int_raster,"maps/EU_confidence_GDI_samp_vs_gdi_smoothed.tif", overwrite=T)

dev.off()


