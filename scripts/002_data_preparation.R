#==============================================================================#
#                                                                              #
#  Prepare OTU tables, manage auxiliary data, and convert raw files to .RData  #
#                                                                              #
#  Written by Karelle Rheault: karh@ign.ku.dk                                  #
#                                                                              #
#==============================================================================#

# Set working directory (modify the path below to match your own)
setwd("some_folder_path/my_working_directory")

# Load packages
source("001_required_packages.R")

# Load auxiliary data ----------------------------------------------------------

# Load an Excel sheet (.xlsx)
aux <- read.xlsx("../raw data/aux_data.xlsx")
