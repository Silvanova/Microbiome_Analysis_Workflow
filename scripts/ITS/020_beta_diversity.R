#==============================================================================#
#                                                                              #
#  Beta diversity: distance matrices, PCoA, PERMANOVA and homogeneity of       #
#  group dispersions                                                           #
#                                                                              #
#  Modified from Sofia Gomes: s.i.gomes@biology.leidenuniv.nl                  #
#  Modified by Davide Barsotti: db@ign.ku.dk                                   #
#                                                                              #
#==============================================================================#

# Set working directory (modify the path below to match your own)
setwd("some_folder_path/my_working_directory")

# Load packages and homemade functions
source("scripts/001_required_packages.R")

# The figures of the ITS scripts are saved in a separate folder
dir.create("figures/ITS", showWarnings = FALSE)

# Load the cleaned OTU table from previous scripts
load("RData/ITS/otu_BLAST_clean.RData")

# Set seeds: for reproducibility, set seeds (random number generation) at the
# start of this script to always obtain the same results with rarefaction and
# permutation tests
set.seed(100)

# Order the land uses from agriculture to the old forest (for the plots and
# tables)
sample_data(otu_BLAST_clean)$Land_use <- factor(
  sample_data(otu_BLAST_clean)$Land_use,
  levels = c("Conventional agriculture", "Organic agriculture", "Grassland",
             "Oak", "Old Ref"))


# Data transformation ----------------------------------------------------------

# Beta diversity quantifies how different the taxonomic composition of two
# samples is. Metabarcoding data are compositional: the number of reads does
# not reflect the total abundance of the fungi in the soil, only the relative
# abundances within a sample can be compared. The data also contain a long
# tail of rare OTUs: a few OTUs are very abundant and most are rare.

## Rarefied data and Hellinger transformation ####

# As in 010_alpha_diversity.R, we start from the data rarefied to 1011 reads
# (same seed, so the same rarefied table).
otu_table_rarefy <- rarefy_even_depth(otu_BLAST_clean, sample.size = 1011,
                                      replace = FALSE, rngseed = 100)
nsamples(otu_table_rarefy) # 87 samples
ntaxa(otu_table_rarefy)    # 5153 OTUs

# Check: all samples should have 1011 reads
summary(sample_sums(otu_table_rarefy))

# The Hellinger transformation is the square root of the relative abundances
# (Legendre & Gallagher, 2001). Dividing by the sample total converts the
# counts to relative abundances, and the square root gives less weight to the
# few very abundant OTUs (it stabilises the variance).
otu_table_hel <- transform_sample_counts(otu_table_rarefy,
                                         function(x) sqrt(x / sum(x)))

## Non-rarefied data for the Aitchison distance ####

# The Aitchison distance is calculated on the non-rarefied data (all 89
# samples). vegan::vegdist() needs a matrix with the samples as rows, so we
# extract the OTU table from the phyloseq object.
otu_mat <- as(otu_table(otu_BLAST_clean), "matrix")
if (taxa_are_rows(otu_BLAST_clean)) otu_mat <- t(otu_mat) # samples = rows
dim(otu_mat) # 89 samples, 9036 OTUs

# Metadata of the samples, for the plots and tests below
meta_rare <- data.frame(sample_data(otu_table_rarefy)) # 87 samples
meta_all  <- data.frame(sample_data(otu_BLAST_clean))  # 89 samples


# Distance matrices ------------------------------------------------------------

# Several distance (or dissimilarity) metrics can be used to compare samples,
# and each of them emphasises different aspects of the communities:
#   - Bray-Curtis: abundance-based, dominated by the abundant OTUs. Here it is
#     calculated on the Hellinger-transformed rarefied data.
#   - Aitchison: the Euclidean distance between centred log-ratio (CLR)
#     transformed abundances. It is designed for compositional data, and
#     does not need rarefied data. The log of 0 is not defined, so a small
#     constant (pseudocount) is added to all counts.
#   - Jaccard: presence-absence only (binary = TRUE), the abundances are
#     ignored. Calculated on the rarefied data.
# Other metrics exist, e.g. UniFrac, which uses the phylogenetic relationships
# between OTUs (not available for ITS, which cannot be aligned across fungi).

# Bray-Curtis dissimilarities on the Hellinger-transformed data
dist_bray <- phyloseq::distance(otu_table_hel, method = "bray")
head(dist_bray)

# Aitchison distances on the non-rarefied data, with a pseudocount of 1
dist_ait <- vegdist(otu_mat, method = "aitchison", pseudocount = 1)
head(dist_ait)

# NOTE: the tutorial uses the "robust.aitchison" distance, which only uses the
# non-zero counts instead of adding a pseudocount (see 006_normalization.R).
# Both are valid. With a pseudocount, absent OTUs and OTUs with very few reads
# get similar values, which is important to interpret the PCoA below.

# Jaccard distances on the rarefied data (presence-absence)
dist_jac <- phyloseq::distance(otu_table_rarefy, method = "jaccard",
                               binary = TRUE)
head(dist_jac)


# Principal Coordinates Analysis (PCoA) ----------------------------------------

# A distance matrix contains one value for each pair of samples (87 samples:
# 3741 pairs), too many to look at. An ordination such as the Principal
# Coordinates Analysis (PCoA) summarises the distances in a few axes, so that
# similar samples are close to each other in the plot. phyloseq::ordinate()
# runs the PCoA with ape::pcoa().
pcoa_bray <- ordinate(otu_table_hel, method = "PCoA", distance = dist_bray)
pcoa_ait  <- ordinate(otu_BLAST_clean, method = "PCoA", distance = dist_ait)
pcoa_jac  <- ordinate(otu_table_rarefy, method = "PCoA", distance = dist_jac)

# Proportion of the variation explained by the first axes
head(pcoa_bray$values$Relative_eig)
head(pcoa_ait$values$Relative_eig)
head(pcoa_jac$values$Relative_eig)

## Plot ####

# Colours for each site and each land use
site_col <- c("Conv_agri" = "bisque", "Org_agri" = "bisque3",
              "Grass" = "orange",
              "Oak_20" = "darkolivegreen1", "Oak_28" = "darkolivegreen2",
              "Oak_33" = "darkolivegreen3", "Oak_44" = "palegreen3",
              "Oak_51" = "darkolivegreen4", "Oak_54" = "darkolivegreen",
              "Oak_ref" = "turquoise4")
land_use_col <- c("Conventional agriculture" = "bisque",
                  "Organic agriculture" = "bisque3",
                  "Grassland" = "orange",
                  "Oak" = "darkolivegreen2",
                  "Old Ref" = "turquoise4")

# The three PCoA plots use the same code, so we write it in a function.
# plot_pcoa() draws the samples on the first two axes, coloured by site, with
# an ellipse for each land use:
#   - physeq: the phyloseq object used for the distances (for the metadata)
#   - pcoa: the result of ordinate()
#   - title: the title of the plot
plot_pcoa <- function(physeq, pcoa, title) {
  # Sample scores on the first two axes, with the metadata of the samples
  scores <- plot_ordination(physeq, pcoa, justDF = TRUE)
  # Variation explained by the first two axes, for the axis labels
  var_expl <- round(100 * pcoa$values$Relative_eig[1:2], 1)

  ggplot(scores, aes(x = Axis.1, y = Axis.2)) +
    stat_ellipse(aes(fill = Land_use), geom = "polygon", alpha = 0.4) +
    geom_point(aes(colour = SiteID), size = 2.5) +
    scale_colour_manual(values = site_col, breaks = names(site_col)) +
    scale_fill_manual(values = land_use_col) +
    labs(title = title,
         x = paste0("Axis 1 (", var_expl[1], "%)"),
         y = paste0("Axis 2 (", var_expl[2], "%)"),
         colour = "Site", fill = "Land use") +
    coord_fixed() +
    theme_bw() +
    theme(panel.grid = element_blank())
}

p_pcoa_bray <- plot_pcoa(otu_table_hel, pcoa_bray,
                         "Bray-Curtis (Hellinger, rarefied)")
p_pcoa_bray
ggsave("figures/ITS/020_pcoa_bray.png", p_pcoa_bray, width = 8, height = 6)

p_pcoa_ait <- plot_pcoa(otu_BLAST_clean, pcoa_ait,
                        "Aitchison (CLR, non-rarefied)")
p_pcoa_ait
ggsave("figures/ITS/020_pcoa_aitchison.png", p_pcoa_ait, width = 8, height = 6)

p_pcoa_jac <- plot_pcoa(otu_table_rarefy, pcoa_jac,
                        "Jaccard (presence-absence, rarefied)")
p_pcoa_jac
ggsave("figures/ITS/020_pcoa_jaccard.png", p_pcoa_jac, width = 8, height = 6)

## Compare the distance metrics ####

p_pcoa_all <- ggarrange(p_pcoa_bray, p_pcoa_ait, p_pcoa_jac,
                        ncol = 3, common.legend = TRUE, legend = "right")
p_pcoa_all
ggsave("figures/ITS/020_pcoa_comparison.png", p_pcoa_all,
       width = 16, height = 6)

# In the three PCoAs, the first axis separates the agricultural sites from the
# grassland and the forests: this is the largest difference in the fungal
# communities. Conventional and organic agriculture overlap. With Bray-Curtis
# and Jaccard, the old reference forest (Oak_ref) is clearly separated from
# the oak plantations (only one sample of Oak_51 lies close to it). With
# Aitchison, the samples of the reference forest are spread along the second
# axis, and their ellipse overlaps with the plantations.
#
# The metrics therefore lead to different conclusions about how similar the
# plantations are to the reference forest. Which one to use? The Jaccard
# distances only use presence-absence, and show the same separation as
# Bray-Curtis: the plantations and the reference forest differ in which OTUs
# are present, not only in their abundances. With the Aitchison distance, a
# pseudocount is added to all counts before the log-ratio transformation, so
# absent OTUs and OTUs with very few reads get similar values. This can blur
# the differences caused by OTUs that are present in one group and absent in
# the other.


# PERMANOVA --------------------------------------------------------------------

# The PCoA shows the patterns, but does not test them. The permutational
# multivariate analysis of variance (PERMANOVA) tests whether the centroids
# (multivariate means) of groups of samples differ. vegan::adonis2() calculates
# a pseudo-F ratio, similar to the F ratio of an ANOVA, and compares it with
# the values obtained after randomly permuting the samples between groups.
# H0: the centroids of all groups are equivalent.

## Site and land use ####

# Test site and land use in the same model, each one after the other
# (by = "margin")
set.seed(100)
adonis2(dist_bray ~ SiteID + Land_use, data = meta_rare,
        permutations = 999, by = "margin")

# Land_use has 0 degrees of freedom and no F value. Why? Look at the metadata:
table(meta_rare$SiteID, meta_rare$Land_use)

# Each site belongs to a single land use: SiteID is nested within Land_use.
# Once the sites are in the model, land use cannot explain anything more:
# the two variables are redundant. We therefore test each variable alone.
set.seed(100)
adonis2(dist_bray ~ SiteID, data = meta_rare, permutations = 999)

set.seed(100)
adonis2(dist_bray ~ Land_use, data = meta_rare, permutations = 999)

# Both site and land use are significantly associated with the composition of
# the fungal communities (p = 0.001, the smallest p-value possible with 999
# permutations). Sites explain 43% of the variation (R2 = 0.43), land use 30%
# (R2 = 0.30). Sites explain more because they include the land use and the
# differences between the sites of the same land use (e.g., plantations of
# different ages).

# With by = "terms", the variables are tested in sequence: each term only gets
# the variation not explained by the terms before it. Land use first, then
# the differences between the sites of the same land use:
set.seed(100)
adonis2(dist_bray ~ Land_use + SiteID, data = meta_rare,
        permutations = 999, by = "terms")

# Land use explains 30% of the variation, and the sites explain another 13%
# within land uses: 30% + 13% = the 43% explained by the sites alone.

## Compare the distance metrics ####

# Is the effect of land use the same with the three distance metrics?
# The Aitchison distances include the two samples removed by rarefaction, so
# they use the metadata of all samples.
set.seed(100)
adonis2(dist_ait ~ Land_use, data = meta_all, permutations = 999)

set.seed(100)
adonis2(dist_jac ~ Land_use, data = meta_rare, permutations = 999)

# Land use is significant with the three metrics, and explains a similar part
# of the variation with Bray-Curtis (R2 = 0.30) and Aitchison (R2 = 0.28), but
# less with Jaccard (R2 = 0.19). Presence-absence gives the same weight to all
# OTUs: the many rare OTUs, found by chance in only a few samples, add
# variation that is not related to land use.

## Soil properties ####

# We can also test the effect of soil properties. Samples of the same land use
# are not independent. With the strata argument, the samples are only
# permuted within each land use (similar to a random factor in mixed models),
# so the test asks whether pH and phosphate explain the composition within
# land uses.
set.seed(100)
adonis2(dist_bray ~ pH + PO4, data = meta_rare, permutations = 999,
        by = "margin", strata = meta_rare$Land_use)

# Both pH (R2 = 0.13) and phosphate (R2 = 0.03) are significant (p = 0.001).
# NOTE: strata only changes how the samples are permuted, not the R2, which
# is still calculated on the total variation. pH differs strongly between
# land uses, so part of its R2 may also reflect the effect of land use.


# Homogeneity of group dispersions ---------------------------------------------

# Homogeneity of the dispersions is not an assumption of PERMANOVA, but it
# affects its results: a significant PERMANOVA can be caused by differences in
# the spread of the samples within groups (dispersion), not only by different
# centroids. vegan::betadisper(), a multivariate analogue of Levene's test,
# calculates the distance of each sample to the centroid of its group. We then
# test whether these distances differ between groups.

## Bray-Curtis ####

disp_bray <- betadisper(dist_bray, meta_rare$SiteID)
anova(disp_bray)

## Aitchison ####

disp_ait <- betadisper(dist_ait, meta_all$SiteID)
anova(disp_ait)

## Plot ####

# Distance to the centroid of each sample, for both metrics
disp_df <- rbind(
  data.frame(meta_rare[, c("SiteID", "Land_use")],
             distance = disp_bray$distances, metric = "Bray-Curtis"),
  data.frame(meta_all[, c("SiteID", "Land_use")],
             distance = disp_ait$distances, metric = "Aitchison"))
disp_df$metric <- factor(disp_df$metric, levels = c("Bray-Curtis", "Aitchison"))

# Order the sites by land use for the plot
disp_df$SiteID <- factor(disp_df$SiteID,
                         levels = unique(disp_df$SiteID[order(disp_df$Land_use,
                                                              disp_df$SiteID)]))

p_disp <- ggplot(disp_df, aes(x = SiteID, y = distance)) +
  geom_boxplot(aes(fill = Land_use), outlier.shape = NA) +
  geom_jitter(size = 1, width = 0.1, height = 0) +
  scale_fill_manual(values = land_use_col) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = NULL, y = "Distance to the site centroid", fill = "Land use") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_disp
ggsave("figures/ITS/020_dispersion.png", p_disp, width = 10, height = 5)

# The dispersions differ significantly between sites with both metrics
# (Bray-Curtis: F = 19.9, p < 0.001; Aitchison: F = 12.0, p < 0.001). The
# boxplots show that the agricultural samples have a lower dispersion than the
# other sites: their communities are more homogeneous (with Aitchison, also
# the youngest plantation, Oak_20). The reference forest
# has one of the highest dispersions. Part of the significant PERMANOVA may
# therefore be due to differences in dispersion, although the PCoA also shows
# a clear separation of the groups.

## Without the agricultural sites ####

# Do the other sites have similar dispersions? We remove the agricultural
# samples from the distance matrix and test again.
keep <- !meta_rare$Land_use %in% c("Conventional agriculture",
                                   "Organic agriculture")
dist_bray_noagri <- as.dist(as.matrix(dist_bray)[keep, keep])
disp_bray_noagri <- betadisper(dist_bray_noagri, meta_rare$SiteID[keep])
anova(disp_bray_noagri)

# Without the agricultural sites, the dispersions still differ between sites
# (F = 3.6, p = 0.003), but much less: the differences in dispersion are
# mainly due to the agricultural sites. Among the other sites, Oak_51 and
# Oak_ref are the most dispersed.


# Exercise ---------------------------------------------------------------------

# Is there an effect of age on the composition of the fungal communities in
# the forests (oak plantations and old reference forest)?
#
# Hints:
#   - Keep only the samples of the land uses "Oak" and "Old Ref", e.g. with
#     subset_samples() on otu_table_hel, then calculate the Bray-Curtis
#     distances again.
#   - The column Age of the metadata gives the age of the forest (years).
#   - Is age confounded with site? How many sites are there for each age?
