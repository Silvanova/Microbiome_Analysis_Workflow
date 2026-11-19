#==============================================================================#
#                                                                              #
#  Alpha diversity: rarefied and extrapolated richness, Shannon and Simpson    #
#  diversity, and comparisons between sites                                    #
#                                                                              #
#  Written by Davide Barsotti: db@ign.ku.dk                                    #
#  Modified from Sofia Gomes: s.i.gomes@biology.leidenuniv.nl                  #
#                                                                              #
#==============================================================================#

# Set working directory (modify the path below to match your own)
setwd("some_folder_path/my_working_directory")

# Load packages and homemade functions
source("scripts/001_required_packages.R")

# The figures of the ITS scripts are saved in a separate folder
dir.create("figures/ITS", showWarnings = FALSE)

# Load the cleaned OTU table from previous scripts
load("RData/16S/otu_BLAST_clean.RData")

# Set seeds: for reproducibility, set seeds (random number generation) at the
# start of this script to always obtain the same results with rarefaction
set.seed(100)


# Rarefaction ------------------------------------------------------------------

# Alpha diversity, and richness in particular, depends strongly on sequencing
# depth: samples with more reads inevitably contain more OTUs (see
# 006_normalization.R). Before comparing alpha diversity between samples, the
# data must therefore be standardized. Here we rarefy the data: all samples are
# randomly subsampled to the same number of reads.

# How to choose the rarefaction depth, and how many reads and OTUs are lost by
# rarefying, is explained in 006_normalization.R. In this example, the lowest
# sequencing depth (24663 reads) is high enough to keep all samples.

# Save the chosen rarefaction depth (you can change it, see Exercise 1)
rare_depth <- 24663

# phyloseq::rarefy_even_depth() subsamples all samples to the same depth.
#   - sample.size: the rarefaction depth. Samples with fewer reads are removed.
#   - replace = FALSE: subsampling without replacement (like drawing reads
#     out of the sample one at a time)
#   - rngseed: the seed for the random subsampling, for reproducibility
# OTUs that are no longer present in any sample after subsampling are removed
# automatically.
otu_table_rarefy <- rarefy_even_depth(otu_BLAST_clean, sample.size = rare_depth,
                                      replace = FALSE, rngseed = 100)

# Check the number of samples, reads per sample, and OTUs after rarefaction
nsamples(otu_table_rarefy) # 90 samples
sample_sums(otu_table_rarefy)
ntaxa(otu_table_rarefy)    # 7681 OTUs (8056 before rarefaction)

# Rarefying discards more than 80% of the reads, from the samples with a high
# sequencing depth. An alternative to rarefying is to extrapolate alpha
# diversity based on the accumulation of OTUs with increasing numbers of reads
# (see below).


# Interpolation and extrapolation of diversity with iNEXT ----------------------

# The R package iNEXT draws rarefaction (interpolation) and extrapolation
# curves, and estimates the diversity each sample would have if it was
# sequenced completely (asymptotic diversity). phyloseq has no equivalent
# function for this, so we extract the OTU table from the phyloseq object and
# use iNEXT directly.

# iNEXT uses Hill numbers of order q:
#   - q = 0: species richness (number of OTUs)
#   - q = 1: Shannon diversity, the exponential of the Shannon entropy
#   - q = 2: Simpson diversity, the inverse of the Simpson concentration
# Larger numbers always mean higher diversity. Hill numbers are expressed in
# "effective number of OTUs": the number of equally abundant OTUs that would
# give the same diversity. As q increases, rare OTUs get less weight.

## Run iNEXT for each sample ####

# iNEXT expects a data frame with OTUs as rows and samples as columns.
# Use the non-rarefied data: iNEXT does its own rarefaction and extrapolation.
otu_df <- as.data.frame(as(otu_table(otu_BLAST_clean), "matrix"))
if (!taxa_are_rows(otu_BLAST_clean)) otu_df <- as.data.frame(t(otu_df))

# Run iNEXT for each sample individually
#   - q: the orders of the Hill numbers to calculate
#   - datatype = "abundance": the data are read counts per OTU
#   - endpoint: the number of reads up to which the curves are drawn. Beyond
#     the library size of a sample, the curve is extrapolated. Extrapolation is
#     only reliable up to about twice the library size of the sample.
#   - knots: the number of points calculated along each curve. Fewer knots make
#     the calculation faster.
# NOTE: iNEXT takes a long time to run on a full dataset (possibly hours for
# the 90 samples of this example), because it calculates the curves and their
# confidence intervals (by bootstrapping) for each sample and each order q. Run this line before
# a long break or overnight, and save the result (see below) so that you do
# not have to run it again.
inext_res <- iNEXT(otu_df, q = c(0, 1, 2), datatype = "abundance",
                   endpoint = 2 * rare_depth, knots = 20)

# Save the iNEXT results, and load them next time instead of running iNEXT
save(inext_res, file = "RData/16S/inext_res.RData")
# load("RData/16S/inext_res.RData")

## Plot the curves ####

# Plot the rarefaction (solid line) and extrapolation (dashed line) curves of
# each sample, with one panel per order q. The vertical line shows the
# rarefaction depth.
p_inext <- ggiNEXT(inext_res, type = 1, facet.var = "Order.q", grey = TRUE) +
  geom_vline(xintercept = rare_depth) +
  theme(legend.position = "none")
p_inext
ggsave("figures/ITS/010_inext_curves.png", p_inext, width = 12, height = 5)

# For Shannon and Simpson diversity (q = 1 and q = 2), the curves reach a
# plateau before the rarefaction depth: the rarefied data capture the diversity
# of the samples. For richness (q = 0), the curves are still increasing, which
# means that more OTUs would be found with more sequencing. Richness is very
# sensitive to rare OTUs and to sequencing depth, and therefore requires the
# most caution.

## Estimated (asymptotic) diversity ####

# iNEXT returns a table with the observed and estimated (asymptotic) diversity
# of each sample, for each order q. For richness (q = 0), the estimator is the
# Chao1 estimator.
inext_div <- inext_res$AsyEst
head(inext_div)

# Note: Chao1 estimates the number of undetected OTUs from the number of OTUs
# that occur only once (singletons) or twice (doubletons) in a sample. In
# 004_global_filter.R, we removed singletons and rare OTUs (< 7 reads in the
# whole dataset). This removes many of the singletons and doubletons the
# estimator relies on, so the estimated richness is likely underestimated.
# Keep this in mind when interpreting estimated richness on filtered data.


# Rarefied vs estimated richness -----------------------------------------------

# What is the impact of using the estimated species richness instead of the
# rarefied species richness?

## Calculate richness ####

# Rarefied richness: phyloseq::estimate_richness() calculates several alpha
# diversity measures at once. "Observed" is the number of OTUs in each sample.
SR <- estimate_richness(otu_table_rarefy, measures = "Observed") %>%
  dplyr::rename(SRrare = Observed)

# estimate_richness() returns the samples in the same order as the phyloseq
# object, with sample names as row names
SR$Assemblage <- rownames(SR)

# Estimated richness from iNEXT
SR_est <- inext_div %>%
  filter(Diversity == "Species richness") %>%
  select(Assemblage, Estimator)

# Combine the two tables, keeping only the samples present in both
SR <- inner_join(SR, SR_est, by = "Assemblage")

# Add the sample data (site and land use)
meta <- as(sample_data(otu_table_rarefy), "data.frame")
SR$SiteID <- meta[SR$Assemblage, "SiteID"]
SR$Land_use <- meta[SR$Assemblage, "Land_use"]

## Correlation between the two measures ####

# R2 of the linear regression between the two measures
r2 <- summary(lm(Estimator ~ SRrare, data = SR))$r.squared

p_rare_vs_est <- ggplot(SR, aes(x = SRrare, y = Estimator)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  labs(x = "Rarefied richness (number of OTUs)",
       y = "Estimated richness (Chao1)",
       subtitle = paste("R2 =", round(r2, 3))) +
  theme_bw()
p_rare_vs_est
ggsave("figures/ITS/010_richness_rarefied_vs_estimated.png", p_rare_vs_est,
       width = 6, height = 5)

# Both measures are correlated. This is good, because we can use one or the
# other and obtain a similar general pattern. But does the choice matter for
# the comparison of richness between sites?

# phyloseq alternative: estimate_richness() can also calculate Chao1 directly
# from the non-rarefied data. It uses vegan::estimateR(), which is based on
# the same principle as iNEXT. Compare it with the iNEXT estimate:
#   chao1 <- estimate_richness(otu_BLAST_clean, measures = "Chao1")
#   plot(chao1[SR$Assemblage, "Chao1"], SR$Estimator)


# Compare richness between sites -----------------------------------------------

# We test if richness differs between sites using an analysis of variance
# (ANOVA), followed by a Tukey test to find which sites differ from each other.

# NOTE: Each site has 9 samples. The land use types are not replicated in the
# same way (6 oak forests, but only 1 site for each of the other land uses).
# Therefore, the samples of one site cannot be used as replicates of a land
# use, as this would be pseudoreplication. This is why we compare sites here,
# and colour the plots by land use.

## ANOVA and Tukey test ####

# Rarefied richness
M_rare <- aov(SRrare ~ SiteID, data = SR)
anova(M_rare)

# Estimated richness
M_est <- aov(Estimator ~ SiteID, data = SR)
anova(M_est)

# Check the assumptions of the ANOVA: the residuals should be normally
# distributed (Shapiro-Wilk test) and the variances equal between groups
# (Bartlett test). If not, consider a non-parametric test (Kruskal-Wallis).
shapiro.test(residuals(M_rare))
bartlett.test(SRrare ~ SiteID, data = SR)
shapiro.test(residuals(M_est))
bartlett.test(Estimator ~ SiteID, data = SR)

# Tukey test: compares all pairs of sites
T_rare <- TukeyHSD(M_rare)
T_est <- TukeyHSD(M_est)

# Adjusted p-values of all pairwise comparisons, for both measures
tuk <- data.frame(Tuk_SRrare = T_rare$SiteID[, "p adj"],
                  Tuk_Estimator = T_est$SiteID[, "p adj"])
tuk

# Compact letter display: sites that share a letter are not significantly
# different from each other
L_rare <- multcompLetters4(M_rare, T_rare)
L_est <- multcompLetters4(M_est, T_est)

letters_df <- data.frame(SiteID = names(L_rare$SiteID$Letters),
                         lettersSRrare = L_rare$SiteID$Letters)
letters_df$lettersEstimator <- L_est$SiteID$Letters[letters_df$SiteID]

# Find the position to write the letters above the boxplots
max_values <- SR %>%
  group_by(SiteID) %>%
  summarise(y_SRrare = max(SRrare) * 1.05, y_Estimator = max(Estimator) * 1.05)

label_data <- inner_join(letters_df, max_values, by = "SiteID")

## Plot ####

# Colours for each land use
land_use_col <- c("Conventional agriculture" = "bisque",
                  "Organic agriculture" = "bisque3",
                  "Grassland" = "orange",
                  "Oak" = "darkolivegreen2",
                  "Old Ref" = "darkolivegreen4")

p1 <- ggplot(SR, aes(x = SiteID, y = SRrare)) +
  geom_boxplot(aes(fill = Land_use), outlier.shape = NA) +
  geom_jitter(size = 1, width = 0.1, height = 0) +
  geom_text(data = label_data, aes(x = SiteID, y = y_SRrare, label = lettersSRrare)) +
  scale_fill_manual(values = land_use_col) +
  labs(y = "Rarefied richness (number of OTUs)", fill = "Land use") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p2 <- ggplot(SR, aes(x = SiteID, y = Estimator)) +
  geom_boxplot(aes(fill = Land_use), outlier.shape = NA) +
  geom_jitter(size = 1, width = 0.1, height = 0) +
  geom_text(data = label_data, aes(x = SiteID, y = y_Estimator, label = lettersEstimator)) +
  scale_fill_manual(values = land_use_col) +
  labs(y = "Estimated richness (Chao1)", fill = "Land use") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p1
ggsave("figures/ITS/010_richness_rarefied_per_site.png", p1, width = 8, height = 5)
p2
ggsave("figures/ITS/010_richness_estimated_per_site.png", p2, width = 8, height = 5)

# Plot both side by side
p_richness_sites <- ggarrange(p1, p2, nrow = 1, common.legend = TRUE, legend = "bottom")
p_richness_sites
ggsave("figures/ITS/010_richness_rarefied_vs_estimated_per_site.png", p_richness_sites,
       width = 12, height = 6, bg = "white")

# Compare the letters of both plots. The results are slightly different: be
# aware that the method chosen will have an impact on the results!


# Exercises --------------------------------------------------------------------

## Exercise 1 ####

# Repeat the calculation of richness, but choose another rarefaction depth
# (change rare_depth at the start of the "Rarefaction" section). Check whether
# it changes the comparison of richness between sites.

## Exercise 2 ####

# Calculate Shannon diversity for each site after rarefying, and test its
# statistical significance. Discuss the differences between richness and
# Shannon diversity.
# Hints:
#   - estimate_richness(otu_table_rarefy, measures = "Shannon") returns the
#     Shannon entropy (H'). Use exp() to convert it to Shannon diversity
#     (Hill number q = 1), in the same unit as richness (number of OTUs).
#   - The "InvSimpson" measure gives the Simpson diversity (q = 2).
#   - Reuse the code of the section "Compare richness between sites".
