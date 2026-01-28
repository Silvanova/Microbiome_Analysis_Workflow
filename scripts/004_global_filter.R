#==============================================================================#
#                                                                              #
#  Global filtering: removing singletons, rare and/or scarce OTUs              #
#                                                                              #
#  Written by Karelle Rheault: karh@ign.ku.dk                                  #
#                                                                              #
#==============================================================================#

# Set working directory (modify the path below to match your own)
setwd("some_folder_path/my_working_directory")

# Load packages
source("scripts/001_required_packages.R")

# Load the cleaned OTU table from previous scripts
# "003_negative_control.R" or "003_positive_control.R"
load("RData/16S/otu_ctrl_neg_clean.RData")
load("RData/16S/otu_ctrl_pos_clean.RData")

# Singletons -------------------------------------------------------------------



# Rare OTUs --------------------------------------------------------------------



# Scarce OTUs ------------------------------------------------------------------



# Save cleaned OTU table -------------------------------------------------------

# OTU table with singletons filtered
save(otu_table_single, file = "RData/16S/otu_single_clean.RData")

# OTU table with rare taxa filtered
save(otu_table_rare, file = "RData/16S/otu_rare_clean.RData")

# OTU table with scarce taxa filtered
save(otu_table_scarce, file = "RData/16S/otu_scarce_clean.RData")
