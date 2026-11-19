#==============================================================================#
#                                                                              #
#  Differential abundance: indicator species analysis and regression models    #
#  for multivariate abundances                                                 #
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

# Colours for each land use (as in 010_alpha_diversity.R and
# 020_beta_diversity.R)
land_use_col <- c("Conventional agriculture" = "bisque",
                  "Organic agriculture" = "bisque3",
                  "Grassland" = "orange",
                  "Oak" = "darkolivegreen2",
                  "Old Ref" = "turquoise4")


# Data preparation -------------------------------------------------------------

# Beta diversity (020_beta_diversity.R) showed that the fungal communities
# differ between land uses. Which OTUs are responsible for these differences?
# Several methods can identify OTUs associated with a group of samples. Here
# we use two of them:
#   - indicator species analysis (indicspecies package), a permutation test
#     of the association of each OTU with one group or a combination of
#     groups;
#   - regression models for multivariate abundances (mvabund package), which
#     fit a generalised linear model (GLM) to each OTU.

# As in 010_alpha_diversity.R and 020_beta_diversity.R, we start from the
# data rarefied to 1011 reads (same seed, so the same rarefied table).
otu_table_rarefy <- rarefy_even_depth(otu_BLAST_clean, sample.size = 1011,
                                      replace = FALSE, rngseed = 100)
nsamples(otu_table_rarefy) # 87 samples
ntaxa(otu_table_rarefy)    # 5153 OTUs

# indicspecies and mvabund do not work with phyloseq objects: they need the
# OTU table as a matrix with the samples as rows and the OTUs as columns.
# In our phyloseq object the OTUs are the rows (taxa_are_rows = TRUE), so the
# table is transposed.
taxa_are_rows(otu_table_rarefy)
otu_mat <- t(as(otu_table(otu_table_rarefy), "matrix"))
dim(otu_mat) # 87 samples x 5153 OTUs

# Metadata and taxonomy of the rarefied data
meta_rare <- data.frame(sample_data(otu_table_rarefy))
tax_rare  <- data.frame(tax_table(otu_table_rarefy))[, c("phylum", "class",
                                                         "order", "family",
                                                         "genus", "species")]

# Check that the samples are in the same order in the matrix and the metadata
all(rownames(otu_mat) == rownames(meta_rare)) # TRUE

# Number of samples per land use. The groups are unbalanced: 51 samples in
# the oak plantations (6 sites), 9 in each of the other land uses (1 site).
table(meta_rare$Land_use)


# Indicator species analysis ---------------------------------------------------

# Indicator species analysis (De Caceres & Legendre, 2009) tests, for each
# OTU, how strongly it is associated with a group of samples (here, a land
# use). indicspecies::multipatt() also tests combinations of groups: an OTU
# can be associated with several land uses together (e.g. both agricultural
# land uses), not only with one (De Caceres et al., 2010).
#
# Two association indices are commonly used:
#   - IndVal ("IndVal.g"): combines the specificity of the OTU (is it found
#     mostly in this group?) and its fidelity (is it found in most samples of
#     this group?). It only looks for OTUs that are associated with a group
#     (preference).
#   - Point-biserial correlation ("r.g"): the correlation between the
#     abundance of the OTU and the membership to a group. It can also be
#     negative (avoidance), and it considers the abundances outside the group.
#     As in the tutorial, we use "r.g"; the ".g" versions correct for the
#     different number of samples in each group, important here because the
#     oak group has many more samples.
#
# The p-values are obtained by permutations: the samples are randomly
# assigned to the groups 999 times. This takes about 30 seconds.
set.seed(100)
indval <- multipatt(otu_mat, meta_rare$Land_use, func = "r.g",
                    control = how(nperm = 999))

# The results are in indval$sign: for each OTU, the columns "s.<group>" show
# to which group(s) it is associated (1 = yes), "index" the number of the
# group combination, "stat" the value of the correlation and "p.value" the
# p-value of the permutation test.
head(indval$sign)

# The full summary is long (one line per indicator OTU, grouped by land use
# combination); summary() only shows the OTUs with p-value <= alpha.
summary(indval, alpha = 0.001)

## Multiple testing ####

# One test is done for each of the 5153 OTUs. With a threshold of 0.05,
# about 5% of the OTUs (~250) would be significant by chance only. The
# tutorial therefore uses a stricter threshold of p <= 0.001.
sum(indval$sign$p.value <= 0.05, na.rm = TRUE)  # 895 OTUs
sum(indval$sign$p.value <= 0.001, na.rm = TRUE) # 280 OTUs

# The p-values can also be corrected for multiple testing, e.g. with the
# false discovery rate (FDR) of Benjamini & Hochberg. With 999 permutations,
# the smallest possible p-value is 0.001: to obtain smaller corrected
# p-values, more permutations would be needed.
p_fdr <- p.adjust(indval$sign$p.value, method = "fdr")
sum(p_fdr <= 0.05, na.rm = TRUE)

# We keep the OTUs with p <= 0.001, as in the tutorial
indval_sig <- subset(indval$sign, p.value <= 0.001)
nrow(indval_sig) # 280 OTUs (tutorial: 287)

## Table of the indicator OTUs ####

# Name of the land use(s) each OTU is associated with: the columns of
# indval$comb are all the combinations of land uses, in the order of "index"
indval_sig$Land_uses <- colnames(indval$comb)[indval_sig$index]

# Add the taxonomy, and sort by land use combination and by strength of the
# association
indval_tab <- cbind(indval_sig[, c("Land_uses", "stat", "p.value")],
                    tax_rare[rownames(indval_sig), ])
indval_tab <- indval_tab[order(indval_sig$index, -indval_sig$stat), ]
head(indval_tab, 20)

# Number of indicator OTUs per land use or combination of land uses
sort(table(indval_tab$Land_uses), decreasing = TRUE)

# Most indicator OTUs are associated with a single land use: grassland (74),
# organic agriculture (59), conventional agriculture (57) and the old forest
# (55). Only 2 OTUs are indicators of the oak plantations alone. The six
# plantations differ from each other (020_beta_diversity.R), so few OTUs are
# found consistently in all of them. 20 OTUs are associated with both
# agricultural land uses together.

## Taxonomy of the indicator OTUs ####

# Number of indicator OTUs per land use (or combination) and phylum. The
# OTUs without a phylum are "Unidentified", and the phyla with few indicator
# OTUs are grouped as "Other" (the colour palette has 8 colours).
indval_tab$phylum_plot <- ifelse(is.na(indval_tab$phylum) |
                                   indval_tab$phylum == "NA",
                                 "Unidentified", indval_tab$phylum)
top_phyla <- names(sort(table(indval_tab$phylum_plot), decreasing = TRUE))[1:7]
indval_tab$phylum_plot[!indval_tab$phylum_plot %in% top_phyla] <- "Other"

indval_phylum <- count(indval_tab, Land_uses, phylum = phylum_plot)

# Order the land use combinations by their index (single land uses at the
# top of the plot), and put "Other" and "Unidentified" at the end of the
# legend
indval_phylum$Land_uses <- factor(
  indval_phylum$Land_uses,
  levels = rev(unique(indval_tab$Land_uses)))
indval_phylum$phylum <- factor(
  indval_phylum$phylum,
  levels = c(sort(setdiff(top_phyla, c("Other", "Unidentified"))),
             "Other", "Unidentified"))

p_indval_phylum <- ggplot(indval_phylum,
                          aes(x = n, y = Land_uses, fill = phylum)) +
  geom_col() +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Number of indicator OTUs", y = NULL, fill = "Phylum") +
  theme_bw()
p_indval_phylum
ggsave("figures/ITS/030_indicator_phylum.png", p_indval_phylum,
       width = 9, height = 5)

# Ascomycota are the most common indicators of the agricultural land uses and
# the grassland. Many indicators of the old forest belong to Rozellomycota
# (20 of 55 OTUs), a group of parasites of other fungi and protists, mostly
# known from environmental sequences.

## Plot of selected indicator OTUs ####

# We plot the abundance of the strongest indicator OTU of each single land
# use. The OTU names are long sequence identifiers, so the plot labels use
# the genus (or the lowest identified rank) instead.
top_otus <- indval_tab %>%
  tibble::rownames_to_column("OTU") %>%
  filter(Land_uses %in% levels(meta_rare$Land_use)) %>%
  group_by(Land_uses) %>%
  slice_max(stat, n = 1) %>%
  ungroup() %>%
  arrange(factor(Land_uses, levels = levels(meta_rare$Land_use)))
top_otus[, c("OTU", "Land_uses", "stat", "phylum", "genus")]

# psmelt() converts the phyloseq object into a long data frame (one row per
# OTU and sample), convenient for ggplot
top_df <- psmelt(prune_taxa(top_otus$OTU, otu_table_rarefy))

# Label of each OTU: land use it indicates, and lowest identified rank
top_otus$label <- paste0(
  top_otus$Land_uses, "\n",
  apply(top_otus[, c("phylum", "class", "order", "family", "genus")], 1,
        function(x) tail(x[!is.na(x) & x != "NA"], 1)))
top_df$label <- factor(top_otus$label[match(top_df$OTU, top_otus$OTU)],
                       levels = top_otus$label)

p_indval_otus <- ggplot(top_df, aes(x = Land_use, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(colour = Land_use), size = 2, width = 0.2, height = 0) +
  scale_colour_manual(values = land_use_col) +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Number of reads (rarefied)", colour = "Land use") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p_indval_otus
ggsave("figures/ITS/030_indicator_otus.png", p_indval_otus,
       width = 10, height = 6)

# The strongest indicators of the agricultural land uses and of the grassland
# are abundant in all samples of their land use and rare or absent elsewhere,
# e.g. the Mortierella OTU of organic agriculture (r = 0.97) has 7 to 18 reads
# in every organic agriculture sample, and almost none in the other land uses.
# The indicator of the oak plantations, a Penicillium OTU, has a weaker
# association (r = 0.65): it is absent from the open land uses, but its
# abundance varies a lot between the plantation samples, and it is also found
# in the old forest. The indicator of the old forest, a Solicoccozyma OTU, is
# present in all land uses but much more abundant in the old forest: the
# point-biserial correlation also detects differences in abundance, not only
# in presence.
#
# Indicator OTUs are useful to describe the groups, but remember that the
# taxonomy of many of them is incomplete (e.g. the OTU of conventional
# agriculture is only identified to the family Plectosphaerellaceae).


# Regression models for multivariate abundances --------------------------------

# Indicator species analysis uses permutation tests, which do not account for
# the strong mean-variance relationship of count data: abundant OTUs also
# vary much more between samples (overdispersion). An alternative is to fit
# a generalised linear model (GLM) to each OTU, with a distribution suited to
# overdispersed counts, such as the negative binomial. The mvabund package
# (Wang et al., 2012) fits these models with manyglm() and tests them with
# resampling (bootstrap).
#
# Resampling is slow with many OTUs. To keep the computation time short, we
# only use the OTUs present in at least 10% of the samples (9 of 87): rare
# OTUs, present in a few samples, carry little information for a GLM anyway.
prev <- colSums(otu_mat > 0)
otu_mat_prev <- otu_mat[, prev >= 9]
ncol(otu_mat_prev) # 402 OTUs

# Convert the matrix into an mvabund object
otu_mvabund <- mvabund(otu_mat_prev)

# Fit a negative binomial GLM to each OTU, with the land use as explanatory
# variable
glm_lu <- manyglm(otu_mvabund ~ Land_use, data = meta_rare,
                  family = "negative.binomial")

# Check the model: the residuals should be spread randomly around 0, without
# a fan shape (which would mean that the negative binomial distribution does
# not describe the mean-variance relationship well)
png("figures/ITS/030_manyglm_residuals.png", width = 7, height = 6,
    units = "in", res = 150)
plot(glm_lu)
dev.off()
plot(glm_lu)

# The residuals are spread evenly around 0, without a fan shape: the negative
# binomial distribution is suitable. The vertical lines on the left are OTUs
# that are absent from all samples of a land use (very low predicted
# abundance); they are common with sparse data such as OTU tables.

# Test the effect of the land use with anova(). nBoot is the number of
# bootstrap resamplings: we use 100 for a fast (about 1 minute) but less
# precise result. For a publication, use e.g. 999. p.uni = "adjusted" also
# tests each OTU separately, with p-values adjusted for multiple testing.
set.seed(100)
anova_lu <- anova(glm_lu, nBoot = 100, p.uni = "adjusted")

# Multivariate test: does the land use affect the community as a whole?
anova_lu$table

# Univariate tests: which OTUs differ between land uses? mvabund changes the
# OTU names that start with a number (it adds an "X"), so we take the names
# from our matrix instead (same order).
p_uni <- anova_lu$uni.p["Land_use", ]
otus_glm_sig <- colnames(otu_mat_prev)[p_uni <= 0.05]
length(otus_glm_sig)

# The land use has a significant effect on the community as a whole
# (p = 0.01, the smallest p-value possible with 100 resamplings: 1/101).
# 180 of the 402 OTUs differ significantly between land uses (adjusted
# p <= 0.05).

## Coefficients ####

# The coefficients of the significant OTUs (one column per OTU). The
# reference (intercept) is the first level of Land_use, conventional
# agriculture. A negative coefficient means that the OTU is less abundant in
# that land use than in conventional agriculture, a positive one that it is
# more abundant.
coefs <- coef(glm_lu)[, p_uni <= 0.05]
colnames(coefs) <- otus_glm_sig

# We look at the 6 most abundant significant OTUs
top_glm <- names(sort(colSums(otu_mat_prev[, otus_glm_sig]),
                      decreasing = TRUE))[1:6]
tax_rare[top_glm, c("phylum", "order", "genus")]
round(coefs[, top_glm], 2)

# How to interpret the coefficients? manyglm uses a negative binomial
# distribution with a log link: the coefficients are on a natural log scale.
# To obtain the fold change, take the exponential of the coefficient. For
# example, a coefficient of -1.38 gives exp(-1.38) = 0.25: the OTU is 4
# times less abundant (75% reduction) than in the reference land use. Very
# large coefficients (e.g. -15 or less, +15 or more) mean that the OTU is
# absent from a land use: the log of 0 is minus infinity, and the model
# returns a very large number instead.
round(exp(coefs[, top_glm]), 2)

# Examples:
#   - The Plectosphaerellaceae OTU (the strongest indicator of conventional
#     agriculture) has on average exp(5.36) = 213 reads in conventional
#     agriculture. It is 0.39 times as abundant in organic agriculture (61%
#     less), and 0.01 to 0.03 times as abundant in the other land uses.
#   - The Solicoccozyma OTU (indicator of the old forest) is 6.45 times more
#     abundant in the old forest than in conventional agriculture.
#   - The three Cortinarius OTUs (ectomycorrhizal fungi) are absent from the
#     agricultural land uses and the grassland: the intercept is -14.82 and
#     the coefficients of the forests are very large (+15 to +19). They depend
#     on the oak trees. One of them is found only in the plantations.

## Comparison with the indicator species analysis ####

# How many of the OTUs selected by manyglm are also indicators (p <= 0.001)?
# Only the prevalent OTUs were tested with manyglm, so we compare within
# these 402 OTUs.
indval_prev <- intersect(rownames(indval_sig), colnames(otu_mat_prev))
length(indval_prev)                        # indicators among the 402 OTUs
length(intersect(indval_prev, otus_glm_sig)) # selected by both methods

# Among the 402 prevalent OTUs, 143 are indicators (p <= 0.001), and 118 of
# them are also selected by manyglm. manyglm selects 62 more OTUs: it tests
# any difference in abundance between land uses, while the indicator species
# analysis looks for OTUs associated with a specific group (or combination
# of groups). The two methods largely agree on the OTUs that differ most
# strongly between land uses.

# NOTE: another popular method is DESeq2 (Love et al., 2014), which also fits
# negative binomial GLMs, on the raw (non-rarefied) counts, with its own
# normalisation. phyloseq objects can be converted with
# phyloseq::phyloseq_to_deseq2(). DESeq2 (Bioconductor) is not installed by
# 001_required_packages.R.
