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

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
library(ggplot2)

if (!require("ggpubr", quietly = TRUE))
  install.packages("ggpubr")
library(ggpubr)

if (!require("vegan", quietly = TRUE))
  install.packages("vegan")
library(vegan)

if (!require("openxlsx", quietly = TRUE))
  install.packages("openxlsx")
library(openxlsx)

if (!require("plyr", quietly = TRUE))
  install.packages("plyr")
library(plyr)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr) # always load dplyr last to avoid conflicts between packages

# homemade functions -----------------------------------------------------------

remove_taxa <- function(badTaxa, physeq){
  allTaxa <- taxa_names(physeq)
  cleanTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(cleanTaxa, physeq))
}

