############################# GDI SCRIPTS ##################################
# =================== Karunarathne et al. ==================================
library(terra)
library(dplyr)
library(sf)
library(ggplot2)
library(scales) 

############################# META ANALYSIS #############################
# 1. Search for genetic markers in accession definitions --------------
# Corrected vectors for markers and genomes
# Corrected vectors with matching lengths for markers and genomes
markers1 <- c("matK", "rbcL", "psbA-trnH", "trnL-trnF", "atpB-rbcL", "ycf1", 
              "trnK", "ndhF", "rpl16", "accD", "ndhJ", "rpoB", "rpoC1", "psbB", "psbC", 
              "atpF-atpH", "petA-psbJ", "trnS-trnG", "petD", "clpP", "rps16", "rps4", 
              "rpl20-rps12", "rpl14", "rbcL-atpB", "psbM-trnD", "psbB-psbH", "trnS-rpS4", 
              "GBSSI", "waxy", "PHYC", "LEAFY", "ncpGS", "G3pdh", "Adh", "CHS", "F3H", 
              "cox1", "cox2", "COI","cob", "nad1", "nad2", "nad4", "nad5", "atp1", "atp6", 
              "18S rRNA", "26S rRNA", "28S rRNA", "5.8S rRNA", "ITS", "ITS1", "ITS2", 
              "ETS", "WGS")

# Adding the missing genome entries for two markers
genomes <- c("Chloroplast", "Chloroplast", "Chloroplast", "Chloroplast", "Chloroplast", 
             "Chloroplast", "Chloroplast", "Chloroplast", "Chloroplast", "Chloroplast", 
             "Chloroplast", "Chloroplast", "Chloroplast", "Chloroplast", "Chloroplast", 
             "Chloroplast", "Chloroplast", "Chloroplast", "Chloroplast", "Chloroplast", 
             "Chloroplast", "Chloroplast", "Chloroplast", "Chloroplast", "Chloroplast", 
             "Chloroplast", "Chloroplast", "Chloroplast", "Chloroplast", "Nuclear", 
             "Nuclear", "Nuclear", "Nuclear", "Nuclear", "Nuclear", "Nuclear", "Nuclear", 
             "Mitochondrial", "Mitochondrial", "Mitochondrial", "Mitochondrial", "Mitochondrial",
             "Mitochondrial", "Mitochondrial", "Mitochondrial", "Mitochondrial", 
             "Nuclear", "Nuclear", "Nuclear", "Nuclear", "Nuclear", "Nuclear", "Nuclear", 
             "Nuclear", "Nuclear", "Nuclear")

# Combine the markers and genomes into a data frame
marker_genomes <- data.frame(Marker = markers1, Genome = genomes, stringsAsFactors = FALSE)
marker_genomes<-marker_genomes[order(marker_genomes$Genome),]


markers <- c("matK", "maturase K", "rbcL", "psbA-trnH", "ITS", "ITS1", "ITS2", "trnL-trnF", "atpB-rbcL", "COI", "ycf1", 
             "trnK", "ndhF", "rpl16", "accD", "ndhJ", "rpoB", "rpoC1", "psbB", "psbC", "atpF-atpH", "petA-psbJ", 
             "trnS-trnG", "petD", "clpP", "rps16", "rps4", "rpl20-rps12", "rpl14", "rbcL-atpB", "psbM-trnD", 
             "psbB-psbH", "trnS-rpS4", "GBSSI", "waxy", "PHYC", "LEAFY", "ncpGS", "G3pdh", "Adh", "CHS", "F3H", 
             "cox1", "cox2", "cob", "nad1", "nad2", "nad4", "nad5", "atp1", "atp6", "18S rRNA", "26S rRNA", "28S rRNA", 
             "5.8S rRNA", "ETS", "WGS", "whole genome")

# Function to find markers and other sequence types in definitions
find_sequences_in_def <- function(def, markers) {
  found_markers <- markers[sapply(markers, function(marker) grepl(marker, def, ignore.case = TRUE))]
  return(found_markers)
}


# Function to process each row
process_row <- function(row_data) {
  matched_sequences <- find_sequences_in_def(row_data$def, markers)
  matched_sequences <- unique(matched_sequences)
  if (length(matched_sequences) > 0) {
    row_results <- data.frame()
    for (sequence in matched_sequences) {
      row_results <- rbind(row_results, data.frame(row_data[-1], marker = sequence, stringsAsFactors = FALSE))
    }
    return(row_results)
  }
  return(NULL)
}

# Function to process a single file (input one file at a time)
process_file <- function(f,output_file=NULL) {
  if(is.character(f)){meta_set <- read.table(f, h = TRUE)} else {meta_set<-data.frame(f)}
  df <- data.frame(meta_set[, c("def", "sp_name", "accession", "date", "length", "country", "coord")], stringsAsFactors = FALSE)
  cl <- makeCluster(num_cores)
  clusterExport(cl, list("process_row", "find_sequences_in_def", "markers"))
  #nrow(df)
  results_list <- pblapply(1:nrow(df), function(i) {
    process_row(df[i, , drop = FALSE])  # Send only the specific row
  }, cl = cl)
  
  stopCluster(cl)  # Corrected this line to stop the cluster correctly.
  
  results <- do.call(rbind, results_list)
  # Save the results for the current file
  if(is.null(output_file)){output_file<-gsub(".txt", "", f)}
  output_file <- paste0(output_file, "_marker_assigned.rda")
  saveRDS(results, output_file, compress = "xz")
  return(paste0(basename(output_file), " DONE"))
}




# 2. plots to show species and markers for TDWG LEVE3_COD -------
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


# 2. EU taxonomic levels
eu_new<-read.csv("data/eu_new_taxonomy.csv")

## 2.1. Plot phylogenetic tree with meta data ---------------
library(V.PhyloMaker2)
library(ape)
library(dplyr)
library(tibble)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(RColorBrewer)

# download the phylogenetic tree using phylomaker
df<-read.table("data/EU_sp_taxonomy_for_phylo_tree.txt",h=T)
tree_result <- phylo.maker(sp.list = df)

# Extract the species tree
species_tree <- tree_result$scenario.3
species_tree$tip.label <- gsub(" ", "_", species_tree$tip.label)
df$species <- gsub(" ", "_", df$species)
species_families <- df[match(species_tree$tip.label, df$species), ]
# Select one representative species per family
family_reps <- species_families %>%
  distinct(family, .keep_all = TRUE) %>%
  pull(species)
# Prune tree to only representative species
family_tree <- keep.tip(species_tree, family_reps)
# Change tip labels to family names
species_info <- df[match(family_tree$tip.label, df$species), ]
family_tree$tip.label <- species_info$family
# ladderize and resolve polytomies
family_tree <- ladderize(multi2di(family_tree))

## plot the sequencing on to the tree 
family_stats <- eu_new %>%
  group_by(family) %>%
  summarise(
    total_species = n_distinct(sp_name),
    sequenced = n_distinct(sp_name[no_accession > 0]),
    geo_seq = n_distinct(sp_name[no_accession > 0 & (primary > 0 | secondary > 0)]),
    .groups = "drop"
  )

plot_data <- family_stats[match(family_tree$tip.label, family_stats$family), ]
plot_data <- as.data.frame(plot_data)
row.names(plot_data) <- family_tree$tip.label


# Prepare long format for bars
plot_data_long <- plot_data %>%
  select(family, total_species, sequenced, geo_seq) %>%
  pivot_longer(cols = c(total_species, sequenced, geo_seq), names_to = "category", values_to = "count")

# # Reorder category so bars are stacked bottom-up logically
plot_data_long$category <- factor(plot_data_long$category, levels = c("geo_seq", "sequenced", "total_species"))

plot_data_long$category <- factor(
  plot_data_long$category,
  levels = c("total_species", "sequenced", "geo_seq")
)

## add geom strips to show orders
strip_defs<-data.table::fread("data/meta_analysis_seq_geo_availability_phylo_tree_ordertable2.csv")

p <- ggtree(family_tree, layout = "circular")
p2<-p + geom_tiplab(size = 2, align = TRUE, offset = 0.5)
p <- p + geom_fruit(
  data = plot_data_long,
  geom = geom_col, 
  mapping = aes(y = family, x = count, fill = category),
  orientation = "y",
  offset = 0.1,      
  width = 1,       
  pwidth = 1      
) +
  scale_fill_manual(
    values = c("geo_seq" = "#F4DE35", "sequenced" = 4, "total_species" = 2),
    labels = c("Total", "Sequenced", "Geo + Seq")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, -.2))) +
  theme_void() +
  theme(legend.title = element_blank())

# Add strips one by one
for (i in 1:nrow(strip_defs)) {
  p <- p + geom_strip(
    taxa1 = strip_defs$taxa1[i],
    taxa2 = strip_defs$taxa2[i],
    label = strip_defs$order[i],
    barsize = 1,
    color = "grey40",
    offset = 0.5,
    angle = -90
  )
}

pdf("plots/meta_analysis_seq_geo_availability_phylo_tree_orders2.pdf",w=8,h=5)
p
p2
dev.off()












