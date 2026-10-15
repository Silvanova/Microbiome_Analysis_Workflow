#==============================================================================#
#                                                                              #
#  Standardization of count data: accounting for differing sequencing depths   #
#  using rarefaction or transformations                                        #
#                                                                              #
#  Written by Karelle Rheault: karh@ign.ku.dk                                  #
#                                                                              #
#==============================================================================#

# Set working directory (modify the path below to match your own)
setwd("some_folder_path/my_working_directory")

# Load packages and homemade functions
source("scripts/001_required_packages.R")

# Load the cleaned OTU table from previous scripts
load("RData/16S/otu_BLAST_clean.RData")

# Set seeds: for reproducibility, set seeds (random number generation) at the
# start of this script to always obtain the same results with rarefaction
set.seed(42)


# Rarefaction ------------------------------------------------------------------

# The most commonly used method for standardization of the sequencing depth
# across samples used to be rarefaction. This method resamples the OTU table
# such that all samples have the same library size.

## Rarefaction curves ####

# To assess whether the data are suitable for rarefaction,
# we can draw rarefaction curves.

# vegan::rarecurve() expects a samples x taxa matrix
# Transpose the OTU table if taxa_are_rows is TRUE
otu_mat <- as(otu_table(otu_BLAST_clean), "matrix")
if (taxa_are_rows(otu_BLAST_clean)) otu_mat <- t(otu_mat)

rarecurve(otu_mat, step = 100, label = FALSE)

# What is the smallest library size across samples? Is it sufficient?
min(sample_sums(otu_BLAST_clean))

# In this example, the smallest library size is 24663, which is sufficient to
# rarefy all samples at this depth. If some samples had a small library size,
# you could choose to rarefy to the size of the second, third, etc. sample
# until you find a sufficient depth.

# Rarefy to the library size of the sample with the lowest library size (default)
otu_table_rarefy <- rarefy_even_depth(otu_BLAST_clean, replace = F)

# Rarefy to a chosen library size
chosen_depth <- 24663
otu_table_rarefy <- rarefy_even_depth(otu_BLAST_clean, replace = F,
                                      sample.size = chosen_depth)

## Balanced library sizes ####

# To check if rarefaction is needed, we can compare library sizes between groups.

# Extract sample metadata
meta <- as(sample_data(otu_BLAST_clean), "data.frame")

# Add library size
meta$library_size <- sample_sums(otu_BLAST_clean)

# Plot results
ggplot(meta, aes(x = SiteID, y = library_size)) +
  geom_boxplot() +
  theme_bw()

# Test for homogeneity of variances
bartlett.test(library_size ~ SiteID, data = meta)

# Test for analysis of variance
kruskal.test(library_size ~ SiteID, data = meta)

# These results show that each group contain equal variance and their mean are
# not different. In addition, we previously found that sequencing depth is
# sufficient (>1000 reads per sample), which means that other compositional
# methods for normalization would be better suited than rarefaction.


# Other standardization methods ------------------------------------------------

# The advantage of other standardization methods over rarefaction is to not
# discard any reads. These other methods are more and more commonly used
# over rarefaction.


## Relative abundance transformations ####


### Simple relative abundance ####

# A simple transformation is the relative abundance. This method is presented
# here for demonstration purposes, but is not recommended in practice.

otu_table_rel <- transform_sample_counts(
  otu_BLAST_clean, function(x) x / sum(x))


### Hellinger transformation (sqrt of relative abundance) ####

# Hellinger transformation converts species abundances from absolute to relative
# values (i.e. standardizes the abundances to sample totals) and then square
# roots them. It is recommended to use this transformation with count data to
# account for the arch effect regularly observed in ecological count data.
# This transformation can be applied on both raw or rarefied data.

otu_table_hel <- transform_sample_counts(
  otu_BLAST_clean, function(x) sqrt(x / sum(x)))


## Ratio transformations ####

# Ratio transformations have the advantage of solving compositionality,
# an important issue when doing ordinations.

# Rarefaction and the Helligner transformation preserve abundance data, which
# is constrained by the sample total (the compositional "closure" problem).
# On the other hand, ratios are invariant to total counts and therefore capture
# the scale-independent relationships among taxa.

# Ratio transformations in compositional data analysis are log-ratio based.
# Because these transformations assume strictly positive data - because the
# log(0) is undefined - zero imputation is required with these transformations.
# Caveat: zero replacement changes the data structure and can influence
# downstream inference, especially when zeros are structural rather than
# below-detection limits.


### CLR (centered log-ratio) ####

# CLR transformation reduces data skewness and compositionality bias.
# The transformation requires zero imputation: (x + 1)

otu_table_clr <- transform_sample_counts(otu_BLAST_clean, function(x) {
  lx <- log(x + 1)
  lx - mean(lx)
})


### RCLR (robust CLR) ####

# RCLR is similar to regular CLR but allows data that contains zeros.
# However, this transformation is not a per-sample operation and cannot
# be performed using the phyloseq::transform_sample_counts function.
# For more information, look up the vegan::decostand function.

# Save OTU table using a new name
otu_table_rclr <- otu_BLAST_clean

# Extract OTU table from phyloseq object and convert to a matrix
mat <- as(otu_table(otu_table_rclr), "matrix")

# If taxa are rows in OTU table, the matrix must be transposed to be used
# with vegan::decostand (remember to transpose the matrix again before
# putting it back in the phyloseq object)
taxa_rows <- taxa_are_rows(otu_table_rclr)
if (taxa_rows) mat <- t(mat) # samples = rows for vegan

# rclr transformation
mat_rclr <- decostand(mat, "rclr", impute = TRUE)

# Put back row and column names (decostand removed them)
row.names(mat_rclr) <- row.names(mat)
colnames(mat_rclr) <- colnames(mat)

# transpose again and put it back into phyloseq object
if (taxa_rows) mat_rclr <- t(mat_rclr) # taxa = rows for phyloseq
otu_table(otu_table_rclr) <- otu_table(mat_rclr, taxa_are_rows = taxa_rows)


# Compare the methods using ordinations ----------------------------------------





