

library(data.table)
## Nucleotide diversity ranges of top 10 families ----------
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
p +
  # Heatmap columns (same as before)
  geom_tile(data = min_data, aes(x = x + heatmap_offset, y = y, fill = value), width = column_width, height = 0.8) +
  scale_fill_viridis_c("Min π", option = "viridis", guide = guide_colorbar(order = 1)) +
  
  new_scale_fill() +
  geom_tile(data = mean_data, aes(x = x + heatmap_offset * 2, y = y, fill = value), width = column_width, height = 0.8) +
  scale_fill_viridis_c("Mean π", option = "plasma", guide = guide_colorbar(order = 2)) +
  
  new_scale_fill() +
  geom_tile(data = max_data, aes(x = x + heatmap_offset * 3, y = y, fill = value), width = column_width, height = 0.8) +
  scale_fill_viridis_c("Max π", option = "inferno", guide = guide_colorbar(order = 3)) +
  
  new_scale_fill() +
  geom_tile(data = species_data, aes(x = x + heatmap_offset * 4, y = y, fill = value), width = column_width, height = 0.8) +
  scale_fill_viridis_c("N Species", option = "mako", guide = guide_colorbar(order = 4)) +
  
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
  xlim(NA, max_tree_x + heatmap_offset * 6) +
  theme(legend.position = "right",
        legend.box = "vertical",
        plot.margin = margin(20, 20, 20, 20))



# ===================== TESTS ======================