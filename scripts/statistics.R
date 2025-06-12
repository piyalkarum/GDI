################### STATISTICS #############################
# ===========================================================


### CORRELATION PLOTS  -----------
library(terra)
library(tidyverse)
library(ggplot2)
eu_range <- c(-10, 45, 35, 75)
mean_map   <- rast("maps/GenDiv_mean_smooth_EU_map_coor_reassigned.tif")
rare_cord  <- rast("maps/glmCorrected_genDivIndex_Interpolated_raster.tif")
weighted   <- rast("maps/weighted_genDivIndex_Interpolated_raster_jan31.tif")
PD_map     <- rast("maps/PhyD_map_Interpolated_raster_jan25.tif")
sum_gdi    <- rast("maps/GenDiv_sum_smooth_EU_map_coord_reassigned.tif")
eu_sp_rich <- rast("maps/species_richness_EU_smooth.tif")
n_sample   <- rast("maps/Number_of_sp_sampled_EU.tif")
mean_map <-crop(mean_map,ext(eu_range))
rare_cord <-crop(rare_cord,ext(eu_range))
weighted <-crop(weighted,ext(eu_range))
PD_map <-crop(PD_map,ext(eu_range))
sum_gdi <-crop(sum_gdi,ext(eu_range))
eu_sp_rich <-crop(eu_sp_rich,ext(eu_range))
n_sample <-crop(n_sample,ext(eu_range))
sprop<-n_sample/eu_sp_rich

# Combine into a single stack
r_stack <- c(mean_map, rare_cord, weighted, PD_map, sum_gdi, eu_sp_rich, sprop)
names(r_stack) <- c("uncorrected_mean_GDI", "corrected_mean_GDI", "weighted_GDI", "PD", "sum_GDI", "species_richness", "n_sampled")
df <- as.data.frame(r_stack, xy = TRUE, na.rm = TRUE)
df_long <- df %>%
  pivot_longer(cols = c(uncorrected_mean_GDI, corrected_mean_GDI, weighted_GDI, PD, sum_GDI),
               names_to = "diversity_metric", values_to = "diversity_value")
# === Adjusted R² for species richness model ===
r2_rich <- df_long %>%
  group_by(diversity_metric) %>%
  summarise(
    adj_r2 = summary(lm(diversity_value ~ species_richness))$adj.r.squared
  )
# === Adjusted R² for n_sampled model ===
r2_sampled <- df_long %>%
  group_by(diversity_metric) %>%
  summarise(
    adj_r2 = summary(lm(diversity_value ~ n_sampled))$adj.r.squared
  )
# Merge with coordinates for annotation placement
r2_rich <- r2_rich %>% mutate(x = Inf, y = -Inf, label = paste0("Adj. R² = ", round(adj_r2, 3)))
r2_sampled <- r2_sampled %>% mutate(x = Inf, y = -Inf, label = paste0("Adj. R² = ", round(adj_r2, 3)))
# === Plot 1: vs species richness ===
p1 <- ggplot(df_long, aes(x = species_richness, y = diversity_value)) +
  geom_point(alpha = 0.1, size = 0.5) +
  facet_wrap(~ diversity_metric, scales = "free_y") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +  # 1:1 line
  geom_smooth(method = "lm", color = "#E74C3C", se = FALSE, linewidth = 1) +
  geom_text(data = r2_rich, aes(x = x, y = y, label = label), hjust = 1.05, vjust = -0.5, inherit.aes = FALSE, size = 3) +
  labs(
    title = "Diversity Metrics vs. Species Richness",
    x = "Species Richness",
    y = "Diversity Metric"
  ) +
  theme_minimal()
# === Plot 2: vs number of sampled species ===
df_long2 <- df_long %>% filter(diversity_metric != "PD")
r2_sampled2 <- r2_sampled %>% filter(diversity_metric != "PD")
p2 <- ggplot(df_long2, aes(x = n_sampled, y = diversity_value)) +
  geom_point(alpha = 0.1, size = 0.5) +
  facet_wrap(~ diversity_metric, scales = "free_y") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") + 
  geom_smooth(method = "lm", color = "#2C3E50", se = FALSE, linewidth = 1) +
  geom_text(data = r2_sampled2, aes(x = x, y = y, label = label), 
            hjust = 1.05, vjust = -0.5, inherit.aes = FALSE, size = 3) +
  scale_x_continuous(limits = c(0,0.4)) +
  labs(
    title = "Diversity Metrics vs. Sampled Species Proportion",
    x = "Sampled Species %",
    y = "Diversity Metric"
  ) +
  theme_minimal()
# Display
print(p1)
print(p2)
# Optional: save to file
ggsave("plots/diversity_vs_richness.pdf", p1, width = 8, height = 6)
ggsave("plots/diversity_vs_sampled_species.pdf", p2, width = 6, height = 6)

## cleaner with 1:1 line ----
library(dplyr)

# For species richness
limits_rich <- df_long %>%
  group_by(diversity_metric) %>%
  summarise(
    min_val = min(c(species_richness, diversity_value), na.rm = TRUE),
    max_val = max(c(species_richness, diversity_value), na.rm = TRUE)
  ) %>%
  mutate(min_val = 0, max_val = max_val * 1.1)  # add small buffer

# For n_sampled
limits_sampled <- df_long2 %>%
  group_by(diversity_metric) %>%
  summarise(
    min_val = min(c(n_sampled, diversity_value), na.rm = TRUE),
    max_val = max(c(n_sampled, diversity_value), na.rm = TRUE)
  ) %>%
  mutate(min_val = 0, max_val = max_val * 1.1)

library(ggplot2)
library(ggforce)

p1 <- ggplot(df_long, aes(x = species_richness, y = diversity_value)) +
  geom_point(alpha = 0.1, size = 0.5) +
  facet_wrap(~ diversity_metric, scales = "free", ncol = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_smooth(method = "lm", color = "#D95F02", se = FALSE, linewidth = 0.8) +
  geom_text(data = r2_rich, aes(x = x, y = y, label = label), hjust = 1.05, vjust = -0.5, inherit.aes = FALSE, size = 3.5) +
  labs(
    title = "Diversity Metrics vs. Species Richness",
    x = "Species Richness",
    y = "Diversity Metric"
  ) +
  theme_bw(base_size = 12)

p2 <- ggplot(df_long2, aes(x = n_sampled, y = diversity_value)) +
  geom_point(alpha = 0.1, size = 0.5) +
  facet_wrap(~ diversity_metric, scales = "free", ncol = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_smooth(method = "lm", color = "#1B9E77", se = FALSE, linewidth = 0.8) +
  geom_text(data = r2_sampled2, aes(x = x, y = y, label = label), hjust = 1.05, vjust = -0.5, inherit.aes = FALSE, size = 3.5) +
  scale_x_continuous(limits = c(0,0.4)) +
  labs(
    title = "Diversity Metrics vs. Sampled Species Proportion",
    x = "Sampled Species %",
    y = "Diversity Metric"
  ) +
  theme_bw(base_size = 12)

ggsave("plots/diversity_vs_richness_w1to1.pdf", p1, width = 8, height = 6)
ggsave("plots/diversity_vs_sampled_species_w1to1.pdf", p2, width = 8, height = 6)


## SUMMARY STATS ON SPECIES/FAMILIES/REGIONS -----------------
gdiv_summary<-read.table("/home/pkaru/piyal/projects/GeoGenoDiv/marvinK_eu_GGD/processed_data/gdiv_summary_europe_coordinate_reassigned.tsv",h=T)
rags<-as.numeric(gdiv_summary$gdi_max)-as.numeric(gdiv_summary$gdi_min)
rags[is.infinite(rags)]<-NA
rags_no0<-rags
rags_no0[rags_no0==0]<-NA
gdiv_summary[which.max(rags),] # 1. max range
gdiv_summary[which.min(rags_no0),] # 3. min range

max(as.numeric(gdiv_summary$gdi_max),na.rm = T) # 2. max gd
gdiv_summary[as.numeric(gdiv_summary$gdi_max)>0.15,] #  max gd species

# range of the maximum and minimum pi values per species
min_range<-gdiv_summary$gdi_min
min_range[is.infinite(min_range)]<-NA
max_range<-gdiv_summary$gdi_max
max_range[is.infinite(max_range)]<-NA
max_range[max_range==0]<-NA
range(max_range,na.rm=T)

## check maps of min and max species
plot(rast(fls[grep("Sisymbrium_officinale",fls)][1])) # max range
plot(rast(fls[grep("Aconitum_napellus",fls)][2]))

## family level variation
hi_fam<-gdiv_summary[gdiv_summary$gdi_max>0.05,]
par(mar=c(10,2,2,1))
barplot(sort(table(hi_fam$family),decreasing = T),col = hcl.colors(28), border = F,las=2, main="Species with pi>0.05 per family")
# low fams
low_fam<-gdiv_summary[gdiv_summary$gdi_max<0.05 & gdiv_summary$gdi_max>0,]
par(mar=c(10,4,2,1))
barplot(sort(table(low_fam$family),decreasing = T),col = hcl.colors(123), border = F,las=2, main="Species with pi<0.05 per family")






## META DATA SUMMARIES ----------------
# usable sequence counts Results 2.1 section
fiv_hi<-data.table::fread("data/sp_marker_sequence_counts.tsv", h=T)
fiv_hi$Marker<-gsub(".txt","",fiv_hi$Marker)
fiv_hi$Species<-paste0(fiv_hi$Genus,"_",fiv_hi$Species)

les1000<-fiv_hi$Sequence_Count[fiv_hi$Sequence_Count<1001 & fiv_hi$Sequence_Count>1]
hist(les1000,breaks=500)
mean(les1000)
mean(fiv_hi$Sequence_Count)


### REGIONS WITH HIGH GDI -----------------
wm<-vect("maps/TDWG/level3/level3.shp")
eu_range<-c(-10,70,20,85)
eu<-crop(wm,eu_range)
mean_map<-rast("maps/GenDiv_mean_smooth_EU_map_coor_reassigned.tif")
mean_map<-mask(mean_map,eu)

cl1<-colorRamps::matlab.like(length(unique(values(mean_map),na.rm=T)))
rng<-ext(-10,45,35,75)
plot(wm,border="grey85",col="grey99",ext=rng,main="MeanGDI")#,mar=c(3.1, 3.1, 2.1, 8.1)
plot(mean_map,add=T,smooth=T,col=cl1)


### sampling in major regions ---------------
wm<-vect("maps/TDWG/level3/level3.shp")
eu_range<-c(-10,70,20,85)
rng<-ext(-10,45,35,75)
n_sample<-rast("maps/Number_of_sp_sampled_EU.tif")
map4=n_sample
map4<-project(map4,crs(wm))
cl1<-(hcl.colors(length(unique(values(map4),na.rm=T)),"Inferno"))
plot(wm,border="grey85",col="grey99",ext=rng,main="Number of Sampled Species")#,mar=c(3.1, 3.1, 2.1, 8.1)
plot(map4,add=T,smooth=T,col=cl1)#

mmax_box<-click()
mmax_box<-vect(mmax_box)
mmax<-crop(map4,mmax_box)
range(values(mmax),na.rm = T)
mean(values(mmax),na.rm = T)


## calculate unique species from tdwg level2 regions
wm_l2<-vect("maps/TDWG/level2/level2.shp")
l2_eu<-crop(wm_l2,ext(-10 , 70 , 20 , 85))
plot(l2_eu)
text(l2_eu,labels=l2_eu$LEVEL2_NAM)

all_mats<-readRDS("data/extrapolated_matrix_array_coord_reassigned.rds")
n_species <- dim(all_mats)[3]
species_rasters <- rast(all_mats, extent = ext(-10, 70, 20, 85), crs = "EPSG:4326")
names(species_rasters) <- paste0("sp_", 1:n_species)
presence_rasters <- classify(species_rasters, cbind(NA, NA, 0), include.lowest = TRUE)
presence_rasters[!is.na(species_rasters)] <- 1
l2_eu <- project(l2_eu, crs(species_rasters))
extracted_vals <- extract(presence_rasters, l2_eu)

library(dplyr)
species_per_region <- extracted_vals %>%
  group_by(ID) %>%
  summarise(unique_species = sum(colSums(across(starts_with("sp_")), na.rm = TRUE) > 0))

l2_eu$unique_species <- species_per_region$unique_species[match(1:nrow(l2_eu), species_per_region$ID)]
plot(l2_eu)
text(l2_eu,labels=l2_eu$LEVEL2_COD)
l2df<-data.frame(l2_eu)
write.table(l2df,"stats/species_sampled_in_TDWG_L2.txt",row.names=F,quote=F,sep="\t")


## for level 3 
## calculate unique species from tdwg level2 regions
wm_l3<-vect("maps/TDWG/level3/level3.shp")
l3_eu<-crop(wm_l3,ext(-10 , 70 , 20 , 85))
plot(l3_eu)
text(l3_eu,labels=l3_eu$LEVEL3_NAM)

all_mats<-readRDS("maps/extrapolated_matrix_array_coord_reassigned.rds")
n_species <- dim(all_mats)[3]
species_rasters <- rast(all_mats, extent = ext(-10, 70, 20, 85), crs = "EPSG:4326")
names(species_rasters) <- paste0("sp_", 1:n_species)
presence_rasters <- classify(species_rasters, cbind(NA, NA, 0), include.lowest = TRUE)
presence_rasters[!is.na(species_rasters)] <- 1
l3_eu <- project(l3_eu, crs(species_rasters))
extracted_vals <- extract(presence_rasters, l3_eu)

library(dplyr)
species_per_region <- extracted_vals %>%
  group_by(ID) %>%
  summarise(unique_species = sum(colSums(across(starts_with("sp_")), na.rm = TRUE) > 0))

l3_eu$unique_species <- species_per_region$unique_species[match(1:nrow(l3_eu), species_per_region$ID)]
plot(l3_eu)
text(l3_eu,labels=l3_eu$LEVEL3_COD)
l3df<-data.frame(l3_eu)
write.table(l3df,"stats/species_sampled_in_TDWG_L3.txt",row.names=F,quote=F,sep="\t")



## Nucleotide diversity ranges of top 10 families ----------
all_sp_taxonomy<-read.csv("data/eu_new_taxonomy.csv")
length(unique(all_sp_taxonomy$sp_name))

top10fams<-sort(table(all_sp_taxonomy$family),decreasing=T)[1:10]
top10fams_names<-names(top10fams)

index<-all_sp_taxonomy$family%in%top10fams_names
top10fam_species<-all_sp_taxonomy[index,c(1,2,5:7,20,23,24)]
top10fam_species_with_marker<-top10fam_species[top10fam_species$marker>0,]

# species list used in the GDI
sp_pi_mat<-read.table("stats/sp_mar_for_matrix.tsv",h=T)
top10fam_spp_used_inGDI<-sp_pi_mat[match(top10fam_species$sp_name,sp_pi_mat$Species),]
top10fam_spp_used_inGDI<-data.frame(top10fam_species,top10fam_spp_used_inGDI)
top10fam_spp_used_inGDI<-top10fam_spp_used_inGDI[!is.na(top10fam_spp_used_inGDI$Species),]


# extract pi values and summarize them
fls<-list.files("data/pi_per_sp",full.names = T)
fams<-top10fams_names
fam_pi<-list()
for(i in seq_along(fams)){
  tm<-top10fam_spp_used_inGDI[top10fam_spp_used_inGDI$family==fams[i],"match_name"]
  pi_vals<-NULL
  for(j in seq_along(tm)){
    tryCatch({pi_fl<-fls[grep(tm[j],fls)]
    pi_tab<-data.table::fread(pi_fl)
    pi_vals<-c(pi_vals,pi_tab$gen_dist)},error=function(e){print(tm[j])})
  }
  fam_pi[[i]]<-pi_vals
}
names(fam_pi)<-fams
saveRDS(fam_pi,"data/Pi_per_fam_top10.rds")

# plot them 
library(ggplot2)
library(ggdist)
library(dplyr)
# Convert list to data frame
fam_pi_df <- do.call(rbind, lapply(names(fam_pi), function(fam) {
  data.frame(Family = fam, Pi = fam_pi[[fam]])
}))
fam_pi_df$Family <- factor(fam_pi_df$Family, levels = rev(fams))
fam_pi_df_filtered <- fam_pi_df %>% filter(Pi <= 0.1)
# plot only to 0.025 for clarity
ggplot(fam_pi_df_filtered, aes(x = Pi, y = Family, fill = Family)) +
  stat_halfeye(.width = 0.5, slab_alpha = 0.7) +
  geom_jitter(height = 0.01, width = 0, alpha = 0.3) +
  coord_cartesian(xlim = c(0, 0.025)) +
  theme_minimal() +
  labs(title = "Nucleotide diversity (pi) per family",
       x = expression(pi),
       y = "Family") +
  guides(fill = FALSE)


