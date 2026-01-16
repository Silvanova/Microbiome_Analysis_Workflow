#==============================================================================#
#                                                                              #
#  Install and load packages                                                   #
#                                                                              #
#  Written by Karelle Rheault: karh@ign.ku.dk                                  #
#                                                                              #
#==============================================================================#

# installed from BiocManager ---------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# decontam
if (!require("decontam", quietly = TRUE))
  BiocManager::install("decontam")
library(decontam)

# phyloseq
if (!require("phyloseq", quietly = TRUE))
  BiocManager::install("phyloseq")
library(phyloseq)

# CRAN packages ----------------------------------------------------------------

if (!require("openxlsx", quietly = TRUE))
  install.packages("openxlsx")
library(openxlsx)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr) # always load dplyr last to avoid conflicts between packages

# Homemade functions -----------------------------------------------------------

