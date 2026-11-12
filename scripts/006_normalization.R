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

# Rarefy to a chosen library size:
# chosen_depth <- 24663
#  otu_table_rarefy <- rarefy_even_depth(otu_BLAST_clean, replace = F,
#                                      sample.size = chosen_depth)

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
# account for the arch effect regularly observed in ecological count data, and
# to make Euclidean-based methods (PCA, RDA) suitable for data with many zeros.
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
# Note: the choice of pseudocount (here 1) affects the result, especially for
# rare taxa.

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

# Each standardization method changes the data in a different way. To see what
# this means in practice, we can compare ordinations of the same samples
# calculated with each method. If the main patterns (e.g., the separation of
# the sites) are similar across methods, our conclusions do not depend on the
# choice of standardization. If they differ, we need to understand why before
# moving on to statistical analyses.

# A transformation and a distance measure always go together. Some
# combinations are standard because the distance is suited to the values
# produced by the transformation:
#   1. Rarefied counts         + Bray-Curtis distance
#   2. Relative abundance      + Bray-Curtis distance
#   3. Hellinger-transformed   + Euclidean distance (= Hellinger distance)
#   4. RCLR-transformed        + Euclidean distance (= robust Aitchison distance)

# Note on the Hellinger transformation: the Euclidean distance calculated on
# Hellinger-transformed data is called the Hellinger distance. This is the
# combination proposed by Legendre & Gallagher (2001), and the one implied when
# running a PCA or RDA on Hellinger-transformed data. The Bray-Curtis distance
# can also be calculated on Hellinger-transformed data, and you will find this
# combination in the literature, but it is less straightforward to interpret.

# Note on the Aitchison distance: the Aitchison distance is the Euclidean
# distance calculated on CLR-transformed data. When the data are RCLR-transformed
# instead, the Euclidean distance is called the robust Aitchison distance.
# The Bray-Curtis distance cannot be used on CLR or RCLR data, because these
# contain negative values.


## Calculate the distance matrices ####

# phyloseq::distance() calculates a distance matrix between samples from a
# phyloseq object. It uses vegan::vegdist() for most distance measures.
# The distance matrices are the basis of the ordinations, and we will also use
# them to compare the methods directly (see "Compare the distance matrices").

# 1. Rarefied counts + Bray-Curtis
dist_rare <- phyloseq::distance(otu_table_rarefy, method = "bray")

# 2. Relative abundance + Bray-Curtis
dist_rel <- phyloseq::distance(otu_table_rel, method = "bray")

# 3. Hellinger + Euclidean (= Hellinger distance)
dist_hel <- phyloseq::distance(otu_table_hel, method = "euclidean")

# 4. RCLR + Euclidean (= robust Aitchison distance)
dist_rclr <- phyloseq::distance(otu_table_rclr, method = "euclidean")

# Note: vegan::vegdist() can also calculate the Aitchison and robust Aitchison
# distances directly from the raw counts, without transforming the data first:
#   vegdist(otu_mat, method = "robust.aitchison")
# This gives a different result from the code above, because the zeros are
# handled differently: we imputed them with decostand() when creating
# otu_table_rclr (impute = TRUE is the default), whereas vegdist() does not
# impute them. Be consistent and report which approach you used.


## Calculate the NMDS ordinations ####

# NMDS (non-metric multidimensional scaling) represents the distances between
# samples in a small number of dimensions (here, k = 2). NMDS only uses the
# rank order of the distances, which makes it suitable for any distance measure.

# phyloseq::ordinate() calculates the distance matrix and runs the NMDS
# (vegan::metaMDS) in one step. The distance measure is given by its name,
# e.g., distance = "bray". Additional arguments are passed on to metaMDS():
#   - k: number of dimensions of the ordination
#   - trymax: maximum number of random starts used to find the best solution.
#     NMDS is an iterative method that starts from random configurations. The
#     seed set at the start of this script makes the results reproducible, as
#     long as the script is run from the start.
#   - autotransform = FALSE: by default, metaMDS() applies a square root and a
#     Wisconsin double standardization to count data. We turn this off, so that
#     each ordination is calculated on the data as we transformed them above.
#   - trace = FALSE: hides the progress messages

# Caveat: ordinate() also accepts a distance matrix that was already
# calculated (e.g., distance = dist_rare), but in that case it ignores all
# additional arguments and silently uses the default settings of metaMDS().

nmds_rare <- ordinate(otu_table_rarefy, method = "NMDS", distance = "bray",
                      k = 2, trymax = 100, autotransform = FALSE, trace = FALSE)
nmds_rel  <- ordinate(otu_table_rel, method = "NMDS", distance = "bray",
                      k = 2, trymax = 100, autotransform = FALSE, trace = FALSE)
nmds_hel  <- ordinate(otu_table_hel, method = "NMDS", distance = "euclidean",
                      k = 2, trymax = 100, autotransform = FALSE, trace = FALSE)
nmds_rclr <- ordinate(otu_table_rclr, method = "NMDS", distance = "euclidean",
                      k = 2, trymax = 100, autotransform = FALSE, trace = FALSE)

# Printing an NMDS object shows the stress and whether a convergent solution
# was found. If it says "Best solution was not repeated", increase trymax
# before interpreting the ordination.
nmds_rare

# The stress value indicates how well the 2D ordination represents the
# original distances. As a rule of thumb, stress < 0.1 is good, 0.1-0.2 is
# acceptable, and > 0.2 should be interpreted with caution.
nmds_stress <- c(Rarefied  = nmds_rare$stress,
                 RelAbund  = nmds_rel$stress,
                 Hellinger = nmds_hel$stress,
                 RCLR      = nmds_rclr$stress)
round(nmds_stress, 3)


## Plot the NMDS ordinations ####

# To compare the ordinations, we plot them side by side, with the samples
# coloured by site. The stress value is added to each plot title.

plot_nmds <- function(physeq, ord, title) {
  plot_ordination(physeq, ord, color = "SiteID") +
    geom_point(size = 2) +
    labs(title = title,
         subtitle = paste("Stress =", round(ord$stress, 3))) +
    theme_bw()
}

p_rare <- plot_nmds(otu_table_rarefy, nmds_rare, "Rarefied + Bray-Curtis")
p_rel  <- plot_nmds(otu_table_rel,    nmds_rel,  "Relative abundance + Bray-Curtis")
p_hel  <- plot_nmds(otu_table_hel,    nmds_hel,  "Hellinger + Euclidean")
p_rclr <- plot_nmds(otu_table_rclr,   nmds_rclr, "RCLR + Euclidean (robust Aitchison)")

ggarrange(p_rare, p_rel, p_hel, p_rclr,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")

# When comparing the plots, keep in mind that the orientation of NMDS axes is
# arbitrary: an ordination can be rotated or mirrored without changing its
# meaning. Focus on which samples group together and which groups are
# separated, not on the position of the groups in the plot.


## Compare the distance matrices ####

# The plots give a visual impression, but we can also measure how similar the
# methods are. The Mantel test calculates the correlation between two distance
# matrices. A correlation close to 1 means that both methods rank the
# distances between samples in the same way.

dist_list <- list(Rarefied  = dist_rare,
                  RelAbund  = dist_rel,
                  Hellinger = dist_hel,
                  RCLR      = dist_rclr)

# Create an empty matrix to store the Mantel correlations
mantel_r <- matrix(NA, nrow = length(dist_list), ncol = length(dist_list),
                   dimnames = list(names(dist_list), names(dist_list)))

# Calculate the Mantel correlation (Spearman) for each pair of methods
for (i in names(dist_list)) {
  for (j in names(dist_list)) {
    mantel_r[i, j] <- mantel(dist_list[[i]], dist_list[[j]],
                             method = "spearman", permutations = 0)$statistic
  }
}
round(mantel_r, 2)

# Note: permutations = 0 skips the permutation test, because we are only
# interested in the correlation here. Distances between the same samples
# are not independent, so we would not interpret a p-value in this context.

# Rarefaction and relative abundance usually give very similar results, as
# both are based on proportions. The RCLR transformation is usually the most
# different, because it gives more weight to rare taxa, whereas Bray-Curtis is
# dominated by the most abundant taxa.

# In this example, the rarefied, relative abundance and Hellinger distances are
# almost identical (r >= 0.99), and the four NMDS plots show a similar overall
# pattern: the agricultural sites and the grassland separate from the oak
# forests. The RCLR distances are less correlated with the others (r ~ 0.8),
# and the separation between the oak forest sites looks different.


## Check for library size effects ####

# If the data are not rarefied, differences in library size between samples
# can still affect the results. A simple way to check for this is to test
# whether library size explains a part of the variation in the ordinations.

# vegan::envfit() fits a variable onto an ordination and tests how much of the
# variation in sample positions it explains (r2). We use the original library
# sizes (before rarefaction) for all methods, so that we can check whether
# each method has removed the effect of sequencing depth. Library sizes are
# log-transformed because they vary over an order of magnitude.

# Library sizes are in the same sample order as the ordinations
log_lib <- data.frame(log_library_size = log10(sample_sums(otu_BLAST_clean)))

ord_list <- list(Rarefied  = nmds_rare,
                 RelAbund  = nmds_rel,
                 Hellinger = nmds_hel,
                 RCLR      = nmds_rclr)

# Fit library size onto each ordination and store the r2 and p-values
lib_fit <- t(sapply(ord_list, function(ord) {
  fit <- envfit(ord, log_lib, permutations = 999)
  c(r2 = unname(fit$vectors$r), p = unname(fit$vectors$pvals))
}))
round(lib_fit, 3)

# A low r2 and a non-significant p-value mean that library size does not
# explain the position of the samples in the ordination. If library size has a
# strong effect with the transformations, but not with rarefaction, consider
# using rarefied data. Note: if library size differs between groups (which we
# tested above), its effect cannot be separated from the effect of the groups.

# In this example, library size has no significant effect on the rarefied,
# relative abundance and Hellinger ordinations (r2 < 0.03), but it explains
# about 14% of the variation in the RCLR ordination (p < 0.01). A likely reason
# is that deeper sequenced samples detect more rare taxa (fewer zeros), and the
# RCLR transformation gives more weight to rare taxa. This shows why it is
# important to check for library size effects when the data are not rarefied.

# For a visual check, the ordinations can also be coloured by library size:
sample_data(otu_table_rclr)$log_library_size <- log_lib$log_library_size
plot_ordination(otu_table_rclr, nmds_rclr, color = "log_library_size") +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  labs(title = "RCLR + Euclidean, coloured by library size") +
  theme_bw()





