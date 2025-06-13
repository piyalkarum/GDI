
# Genetic Diversity Index (GDI) Project

This repository contains the code, data, and methodology used in the study assessing species richness, phylogenetic diversity, and spatial patterns of genetic diversity across European vascular plants.

## Overview of the Workflow

The workflow is structured into two main components:
1. **Data Acquisition and Meta-data Analysis**: Downloading, filtering, and processing occurrence records and genetic metadata.
2. **Diversity Analysis and Mapping**: Quantifying and visualizing species richness, phylogenetic diversity, and genetic diversity patterns.

An abstract overview of the key analytical steps is provided in the manuscript (main text Fig. 1A and 1B).

---

## 1. Data Acquisition and Meta-Data Analysis

### Occurrence Data
Species names were curated using:
- **Plants of the World Online (POWO)**  
- **World Checklist of Vascular Plants (WCVP)**
Files: `EU_all_native_species_list.csv`, `eu_new_taxonomy.csv`

A list of 26,242 accepted native European species was curated based on:
- Native and extant status
- Accepted taxonomy
- Geographical restriction to Europe

### Species Distributions
Species ranges were obtained from:
- **Daru (2024) dataset** (11,276 species)
- **Species-specific SDMs** for an additional 1,298 species
Files: `species_richness_EU_smooth.tif`

Environmental predictors were sourced from **WorldClim 2.0** at 1 arc-minute resolution.

### Genetic Metadata
Metadata was retrieved from **NCBI GenBank** using the `rentrez` R package for all 26,242 species.

### Filtering and Georeferencing
Accessions were categorized into:
- **Primary georeferenced** (coordinates available)
- **Secondary georeferenced** (descriptive location)

Secondary locations were geocoded using:
- `tmaptools`
- `tidygeocoder`

Files: `marker_assigned_all_meta_data.csv`

### Marker Assignment
Accessions were scanned for 50+ commonly used markers. One unique accession per species-marker combination was retained for those with at least 5 georeferenced records.
Files: `marker_assigned_all_meta_data.csv`

### Sequence Download
FASTA files were downloaded using **NCBI EDirect**, organized by species-marker.

---

## 2. Diversity Analyses

### Taxonomic Diversity
- Calculated from species distribution models
- Mapped at 0.25° resolution (~400 km² grid cells)
- 12,571 species used
Files: `species_richness_EU_smooth.tif`


### Phylogenetic Diversity
- Faith's PD using the `phyloregion` R package
- Phylogeny: Smith & Brown (2018)
- Weighted and standardized PD calculated
Files: `PhyD_map_Interpolated_raster_jan25.tif`

### Genetic Diversity and Genetic Diversity Index (GDI)

#### Nucleotide Diversity (π)
- Calculated per species-marker
- Sequences aligned with **MAFFT**
- Spatially linked to sampling locations
Files: `sequence_alignments/*.fasta`


#### Spatial Extrapolation
- Applied only within native SDM range
- Implemented using `wind_interpolate()` from the `gSoup` package

#### Multi-Species GDI Indices
1. **sGDI**: Sum of π across species
2. **mGDI**: Mean π per site
3. **cGDI**: Corrected indices to reduce bias
   - **wGDI**: Weighted by richness
   - **lGDI**: Regression-adjusted
   - **rGDI**: Sample-based rarefaction

All GDI values normalized between 0 (no diversity) and 1 (maximum observed).
Files: `pi_per_sp/*.txt`, `extrapolated_matrix_array_coord_reassigned.h5`

---

## 3. Validation of Correction Methods

### Bias Testing
- Correlation and permutation tests to detect sampling bias
- Simulation using synthetic data to test index behavior under varying richness and sampling
Files: `maps/*.tif`

---

## 4. Mapping and Visualization

- Raster processing with `terra` and `sf`
- Smoothing using a 50-cell moving window in `gSoup`
- Maps plotted with `ggplot2`
Files: `plots/*.pdf`

---

## Citation

If using this method or dataset, please cite:

- Karunarathne et al. (2025). Unlocking the Forgotten Dimension of Biodiversity: A Scalable Genetic Diversity Index for Multi-Species Analysis. *bioRxiv*.<https://doi.org/10.1101/2025.06.03.657643>
