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
eu_map<-crop(wm,eu_range)

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

### Data for the main text table 1 
sp_tdwg3<-as.data.frame(eu_map_sf)



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
strip_defs<-data.table::fread("data/meta_analysis_seq_geo_availability_phylo_tree_ordertable2.txt")

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




## Phylogenetic tree with pi heatmaps -----------------------------
library(data.table)

all_sp_taxonomy<-read.csv("data/eu_new_taxonomy1.csv")
length(unique(all_sp_taxonomy$sp_name))
# species list used in the GDI
sp_pi_mat<-read.table("stats/sp_mar_for_matrix.tsv",h=T)
sps<-sp_pi_mat$match_name
# extract pi values and summarize them
fls<-list.files("data/pi_per_sp",full.names = T)
# join all pi tables to one dataframe
all_pi<-NULL
pb<-txtProgressBar(max=length(sps),style = 3,width = 50)
for(i in seq_along(sps)){
  tryCatch({
    setTxtProgressBar(pb,i)
    pi_file<-fls[grep(sps[i],fls)]
    pi_tab<-data.table::fread(pi_file)
    all_pi<-rbind(all_pi,pi_tab)
  },error=function(e){print(e)})
}

# join taxonomy and pi tables
pi_tab <- as.data.table(all_pi,)
all_sp_taxonomy <- as.data.table(all_sp_taxonomy)
pi_tax <- merge(pi_tab, all_sp_taxonomy[, .(sp_name, family, order, clade)], 
                by = "sp_name", all.x = TRUE)
# Calculate mean π per species (if not already done)
pi_species <- pi_tax[!is.na(gen_dist), .(pi_mean = mean(gen_dist, na.rm = TRUE)), by = .(sp_name, family, order, clade)]


## fix missing classifications
pi_tax<-data.frame(pi_tax)
pi_tax<-pi_tax[!is.na(pi_tax$family),]
fams <- unique(na.omit(pi_tax[, c("family","order","clade")]))

for (i in seq_len(nrow(fams))) {
  fam <- as.character(fams$family[i])
  idx <- which(!is.na(pi_tax$family) & pi_tax$family == fam)  # no NAs in index
  pi_tax$order[idx] <- as.character(fams$order[i])
  pi_tax$clade[idx] <- as.character(fams$clade[i])
}

famNA<-unique(pi_tax[, c("family","order","clade")])

famNA[famNA$family=="Asteraceae",2]<-"Asterales"
famNA[famNA$family=="Pteridaceae",2]<-"Polypodiales"
famNA[famNA$family=="Viburnaceae",2]<-"Dipsacales"
famNA[famNA$family=="Fabaceae",2]<-"Fabales"
famNA[famNA$family=="Aspleniaceae",2]<-"Polypodiales"
famNA[famNA$family=="Ophioglossaceae",2]<-"Ophioglossales"
famNA[famNA$family=="Mazaceae",2]<-"Lamiales"
famNA[famNA$family=="Polypodiaceae",2]<-"Polypodiales"
famNA[famNA$family=="Equisetaceae",2]<-"Equisetales"
famNA[famNA$family=="Asphodelaceae",2]<-"Asparagales"
famNA[famNA$family=="Lycopodiaceae",2]<-"Lycopodiales"
famNA[famNA$family=="Psilotaceae",2]<-"Psilotales"
famNA[famNA$family=="Selaginellaceae",2]<-"Selaginellales"

famNA$clade_fix<-famNA$clade
famNA$clade_fix[famNA$order %in% c("Polypodiales","Ophioglossales","Psilotales","Equisetales")]<- "Ferns"
famNA$clade_fix[famNA$order %in% c("Lycopodiales","Selaginellales")]<- "Lycophytes"

for (i in seq_len(nrow(famNA))) {
  fam <- as.character(famNA$family[i])
  idx <- which(!is.na(pi_tax$family) & pi_tax$family == fam)  # no NAs in index
  pi_tax$order[idx] <- as.character(famNA$order[i])
  pi_tax$clade[idx] <- as.character(famNA$clade_fix[i])
}
# Fill missing clade values by matching on order
setDT(pi_tax)
order_clade_lookup <- unique(na.omit(pi_tax[, .(order, clade)]))
pi_tax[is.na(clade), clade := order_clade_lookup[.SD, on = "order", x.clade]]

write.csv(pi_tax,"/Volumes/PiyalKaru/geoGeno/GDI/data/all_sp_pi_tax_joined.csv",row.names = F)

## plotting =======================

library(ggplot2)
library(dplyr)
library(ggh4x) # for facet_nested

pi_tax<-read.csv("/Volumes/PiyalKaru/geoGeno/GDI/data/all_sp_pi_tax_joined.csv")
# Calculate summary statistics at Family level within taxonomic hierarchy
family_summary <- pi_tax %>%
  group_by(clade, order, family) %>%
  summarise(
    mean_pi = mean(gen_dist, na.rm = TRUE),
    sd_pi = sd(gen_dist, na.rm = TRUE),
    median_pi = median(gen_dist, na.rm = TRUE),
    n_species = n_distinct(sp_name),
    n_accessions = n(),
    .groups = "drop"
  ) %>%
  filter(n_species >= 3) %>% # Optional: filter families with few species
  arrange(clade, order, mean_pi)

# Create the pyramid plot
ggplot(family_summary, aes(x = mean_pi, y = reorder(family, mean_pi))) +
  geom_point(aes(size = n_species, color = order), alpha = 0.8) +
  geom_errorbarh(aes(xmin = mean_pi - sd_pi, xmax = mean_pi + sd_pi), 
                 height = 0.2, alpha = 0.6, linewidth = 0.5) +
  facet_nested(clade + order ~ ., scales = "free_y", space = "free") +
  labs(x = "Mean Nucleotide Diversity (π)", y = "Family",
       title = "Nucleotide Diversity Across European Plant Families",
       subtitle = "Error bars show ±1 SD, point size indicates number of species",
       size = "Number of\nSpecies", color = "Order") +
  scale_size_continuous(range = c(1, 6)) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )


# Calculate at Order level
order_summary <- pi_tax %>%
  group_by(clade, order) %>%
  summarise(
    mean_pi = mean(gen_dist, na.rm = TRUE),
    se_pi = sd(gen_dist, na.rm = TRUE) / sqrt(n()),
    n_families = n_distinct(family),
    n_species = n_distinct(sp_name),
    n_accessions = n(),
    .groups = "drop"
  ) %>%
  arrange(clade, mean_pi)


ggplot(order_summary, aes(x = mean_pi, y = reorder(order, mean_pi))) +
  geom_point(aes(size = n_species, color = clade), alpha = 0.8) +
  geom_errorbarh(aes(xmin = mean_pi - se_pi, xmax = mean_pi + se_pi),height = 0.2, alpha = 0.6, linewidth = 0.7) +
  facet_nested(clade ~ ., scales = "free_y", space = "free") +
  labs(x = "Mean Nucleotide Diversity (π)", y = "Order",
       title = "Nucleotide Diversity Across European Plant Orders",
       size = "Number of\nSpecies", color = "Major Clade") +
  scale_size_continuous(range = c(2, 8)) +
  theme_minimal()


## Phylogenetic tree with heatmaps

library(ggtree)
library(ggplot2)
library(dplyr)
library(tidytree)
library(V.PhyloMaker2)
library(tidyr)
library(ggnewscale)

pi_tax<-read.csv("/Volumes/PiyalKaru/geoGeno/GDI/data/all_sp_pi_tax_joined.csv")
# download the phylogenetic tree using phylomaker
sp_unique<-data.frame(unique(pi_tax[,c(1,11)]))
df<-data.frame(species=sp_unique$sp_name,genus=stringr::str_split_fixed(sp_unique$sp_name,"_",n=2)[,1],family=sp_unique$family)
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

tree<-family_tree
write.tree(tree,"family_tree.tre")


# Calculate stats
family_stats <- pi_tax %>%
  group_by(family) %>%
  summarise(
    Min = min(gen_dist, na.rm = TRUE),
    Mean = mean(gen_dist, na.rm = TRUE), 
    Max = max(gen_dist, na.rm = TRUE),
    N_Species = n_distinct(sp_name),
    .groups = "drop"
  )

# Create tree with very thin branches
p <- ggtree(tree, layout = "circular", size = 0.5)
tree_data <- p$data %>% filter(isTip) %>% select(label, x, y, angle)

# Make heatmap very prominent
max_tree_x <- max(tree_data$x)
heatmap_offset <- max_tree_x * 0.125
column_width <- max_tree_x * 0.1

# Prepare data
min_data <- family_stats %>% select(family, value = Min) %>% left_join(tree_data, by = c("family" = "label"))
mean_data <- family_stats %>% select(family, value = Mean) %>% left_join(tree_data, by = c("family" = "label"))
max_data <- family_stats %>% select(family, value = Max) %>% left_join(tree_data, by = c("family" = "label"))
species_data <- family_stats %>% select(family, value = N_Species) %>% left_join(tree_data, by = c("family" = "label"))

# Build plot with family labels
p<-p +
  # Heatmap columns (same as before)
  geom_tile(data = species_data, aes(x = x + heatmap_offset * 1, y = y, fill = value), width = column_width, height = 0.8) +
  scale_fill_viridis_c("N Species", option = "mako", guide = guide_colorbar(order = 1)) +
  new_scale_fill() +
  
  geom_tile(data = min_data, aes(x = x + heatmap_offset*2, y = y, fill = value), width = column_width, height = 0.8) +
  scale_fill_viridis_c("Min π", option = "viridis", guide = guide_colorbar(order = 2)) +
  
  new_scale_fill() +
  geom_tile(data = mean_data, aes(x = x + heatmap_offset * 3, y = y, fill = value), width = column_width, height = 0.8) +
  scale_fill_viridis_c("Mean π", option = "plasma", guide = guide_colorbar(order = 3)) +
  
  new_scale_fill() +
  geom_tile(data = max_data, aes(x = x + heatmap_offset * 4, y = y, fill = value), width = column_width, height = 0.8) +
  scale_fill_viridis_c("Max π", option = "inferno", guide = guide_colorbar(order = 4)) +
  
  # Family labels with circular layout adjustment
  geom_text(
    data = tree_data,
    aes(x = x + heatmap_offset * 4.8, y = y, label = label, angle = angle),
    size = 2.2,
    hjust = 0,
    nudge_x = column_width * 0.3,
    check_overlap = FALSE
  ) +
  
  # Expand limits for labels
  xlim(NA, max_tree_x + heatmap_offset * 13) +
  theme(legend.position = "right",
        legend.box = "vertical",
        plot.margin = margin(20, 20, 20, 20))


# Add stripes for the orders with curved labels perpendicular to the strip
strip_defs <- data.table::fread("data/meta_analysis_seq_geo_availability_phylo_tree_ordertable3.txt")

# Calculate midpoint positions and angles for each strip
strip_data <- data.frame()
for (i in 1:nrow(strip_defs)) {
  # Get the tip data for the two taxa that define the strip
  taxa1_data <- tree_data %>% filter(label == strip_defs$taxa1[i])
  taxa2_data <- tree_data %>% filter(label == strip_defs$taxa2[i])
  
  if (nrow(taxa1_data) > 0 && nrow(taxa2_data) > 0) {
    # Calculate midpoint angle (average of the two angles)
    mid_angle <- mean(c(taxa1_data$angle, taxa2_data$angle))
    
    # Calculate midpoint y position (average of the two y coordinates)
    mid_y <- mean(c(taxa1_data$y, taxa2_data$y))
    
    # Calculate appropriate x position (same offset as your strips)
    strip_x <- max_tree_x + heatmap_offset * 12
    
    strip_data <- rbind(strip_data, data.frame(
      order = strip_defs[i, 1],  # Assuming column 1 contains order names
      x = strip_x,
      y = mid_y,
      angle = mid_angle
    ))
  }
}

# Add the strips
for (i in 1:nrow(strip_defs)) {
  p <- p + geom_strip(
    taxa1 = strip_defs$taxa1[i],
    taxa2 = strip_defs$taxa2[i],
    barsize = 1,
    color = "grey40",
    offset = heatmap_offset * 12,
    angle = -90
  )
}

# Add the curved order labels
p <- p + geom_text(
  data = strip_data,
  aes(x = x, y = y, label = order, angle = angle),
  size = 2.5,
  hjust = 0,
  nudge_x = column_width * 0.5,  # Adjust this value to position labels outside the strips
  color = "grey40",
  fontface = "bold"
)






# ===================== TESTS ======================

# Add stripes for the orders with curved labels
strip_defs <- data.table::fread("data/meta_analysis_seq_geo_availability_phylo_tree_ordertable3.txt")

# Calculate midpoint positions and angles for each strip
strip_data <- data.frame()
for (i in 1:nrow(strip_defs)) {
  # Get the tip data for the two taxa that define the strip
  taxa1_data <- tree_data %>% filter(label == strip_defs$taxa1[i])
  taxa2_data <- tree_data %>% filter(label == strip_defs$taxa2[i])
  
  if (nrow(taxa1_data) > 0 && nrow(taxa2_data) > 0) {
    # Calculate midpoint angle (average of the two angles)
    mid_angle <- mean(c(taxa1_data$angle, taxa2_data$angle))
    
    # Calculate midpoint y position (average of the two y coordinates)
    mid_y <- mean(c(taxa1_data$y, taxa2_data$y))
    
    # Calculate appropriate x position (same offset as your strips)
    strip_x <- max_tree_x + heatmap_offset * 12
    
    strip_data <- rbind(strip_data, data.frame(
      order = strip_defs[i, 1],  # Assuming column 1 contains order names
      x = strip_x,
      y = mid_y,
      angle = mid_angle
    ))
  }
}

# Add the strips
for (i in 1:nrow(strip_defs)) {
  p <- p + geom_strip(
    taxa1 = strip_defs$taxa1[i],
    taxa2 = strip_defs$taxa2[i],
    barsize = 1,
    color = "grey40",
    offset = heatmap_offset * 12,
    angle = -90
  )
}

# Add the curved order labels - CORRECTED ANGLE
p <- p + geom_text(
  data = strip_data,
  aes(x = x, y = y, label = order, angle = angle - 90),  # Subtract 90 degrees to make text parallel to strip
  size = 2.5,
  hjust = 0.5,  # Center the text on the strip
  vjust = 0.5,  # Center vertically
  nudge_x = .5,  # No nudge since we're placing directly on the strip
  color = "grey40",
  fontface = "bold"
)







# Create tree with very thin branches
p <- ggtree(tree, layout = "circular", size = 0.5)
tree_data <- p$data %>% filter(isTip) %>% select(label, x, y, angle)

# Make heatmap very prominent
max_tree_x <- max(tree_data$x)
heatmap_offset <- max_tree_x * 0.125
column_width <- max_tree_x * 0.1

# Prepare data
min_data <- family_stats %>% select(family, value = Min) %>% left_join(tree_data, by = c("family" = "label"))
mean_data <- family_stats %>% select(family, value = Mean) %>% left_join(tree_data, by = c("family" = "label"))
max_data <- family_stats %>% select(family, value = Max) %>% left_join(tree_data, by = c("family" = "label"))
species_data <- family_stats %>% select(family, value = N_Species) %>% left_join(tree_data, by = c("family" = "label"))

# Build plot with family labels
p<-p +
  # Heatmap columns (same as before)
  geom_tile(data = species_data, aes(x = x + heatmap_offset * 1, y = y, fill = value), width = column_width, height = 0.8) +
  scale_fill_viridis_c("N Species", option = "mako", guide = guide_colorbar(order = 1)) +
  new_scale_fill() +
  
  geom_tile(data = min_data, aes(x = x + heatmap_offset*2, y = y, fill = value), width = column_width, height = 0.8) +
  scale_fill_viridis_c("Min π", option = "viridis", guide = guide_colorbar(order = 2)) +
  
  new_scale_fill() +
  geom_tile(data = mean_data, aes(x = x + heatmap_offset * 3, y = y, fill = value), width = column_width, height = 0.8) +
  scale_fill_viridis_c("Mean π", option = "plasma", guide = guide_colorbar(order = 3)) +
  
  new_scale_fill() +
  geom_tile(data = max_data, aes(x = x + heatmap_offset * 4, y = y, fill = value), width = column_width, height = 0.8) +
  scale_fill_viridis_c("Max π", option = "inferno", guide = guide_colorbar(order = 4)) +
  
  
  # Family labels with circular layout adjustment
  geom_text(
    data = tree_data,
    aes(x = x + heatmap_offset * 4.8, y = y, label = label, angle = angle),
    size = 2.2,
    hjust = 0,
    nudge_x = column_width * 0.3,
    check_overlap = FALSE
  ) +
  
  # Expand limits for labels
  xlim(NA, max_tree_x + heatmap_offset * 13) +
  theme(legend.position = "right",
        legend.box = "vertical",
        plot.margin = margin(20, 20, 20, 20))

