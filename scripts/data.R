####################### DATA ###########################
# ======================================================



# A 3D matrix array of all nucleotide diversity values across species and 
# space laid on a 260 x 320 matrix over Europe within the bounding
# box of -10 , 70 , 20 , 85

all_mats<-readRDS("data/extrapolated_matrix_array_coord_reassigned.rds")

# save as an HDF5 file
library(rhdf5)
h5createFile("extrapolated_matrix_array_coord_reassigned.h5")
h5write(all_mats, "data/extrapolated_matrix_array_coord_reassigned.h5", "array_data")


# European species list with individual nucletide diversity and range
sp_pi_mat<-read.table("stats/sp_mar_for_matrix.tsv",h=T)


# All European vascular plant species screened for genetic data
# their distribution, and marker availability used in the study
all_sp_taxonomy<-read.csv("data/eu_new_taxonomy.csv")

