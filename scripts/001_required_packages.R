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

# installed from R Universe ----------------------------------------------------

if (!require("microViz", quietly = TRUE))
  install.packages(
    "microViz",
    repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
  )
library(microViz)

# CRAN packages ----------------------------------------------------------------

if (!require("Biostrings", quietly = TRUE))
  install.packages("Biostrings")
library(Biostrings)

if (!require("openxlsx", quietly = TRUE))
  install.packages("openxlsx")
library(openxlsx)

if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")
library(tidyr)

if (!require("plyr", quietly = TRUE))
  install.packages("plyr")
library(plyr)

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
library(ggplot2)

if (!require("ggpubr", quietly = TRUE))
  install.packages("ggpubr")
library(ggpubr)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr) # always load dplyr last to avoid conflicts between packages

# Homemade functions -----------------------------------------------------------

# Remove taxa from phyloseq objects
remove_taxa <- function(badTaxa, physeq){
  allTaxa <- taxa_names(physeq)
  cleanTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(cleanTaxa, physeq))
}
