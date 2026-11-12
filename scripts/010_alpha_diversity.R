#==============================================================================#
#                                                                              #
#  Alpha diversity: richness, Shannon and Simpson diversity, rarefaction and   #
#  extrapolation curves, and comparisons between sites                         #
#                                                                              #
#  Written by Davide Barsotti: db@ign.ku.dk                                    #
#  Modified from Sofia Gomes: s.i.gomes@biology.leidenuniv.nl                  #
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

# Alpha diversity measures, and richness in particular, depend strongly on
# sequencing depth: the more reads a sample has, the more OTUs are detected.
# For this reason, rarefaction is recommended before comparing alpha diversity
# between samples, even when other standardization methods are used for beta
# diversity (see 006_normalization.R).

# What is the smallest library size across samples?
min(sample_sums(otu_BLAST_clean)) # 24663 reads

# Rarefy to the library size of the sample with the lowest library size (default)
otu_table_rarefy <- rarefy_even_depth(otu_BLAST_clean, replace = F)

# Check the number of reads per sample: all samples should have the same depth
sample_sums(otu_table_rarefy)

# Save the rarefaction depth to use it in the plots below
rare_depth <- min(sample_sums(otu_table_rarefy))

# How many OTUs remained after rarefaction?
ntaxa(otu_BLAST_clean)
ntaxa(otu_table_rarefy)

# OTUs that only had a few reads in the samples with large library sizes can be
# lost during rarefaction (rarefy_even_depth removes them automatically).


# Rarefaction and extrapolation ------------------------------------------------

# Rarefaction discards reads. To check whether the rarefaction depth is
# sufficient to capture the diversity of the samples, we can draw rarefaction
# and extrapolation curves with the package iNEXT. These curves show how
# diversity increases with the number of reads, and estimate how it would
# continue to increase with more reads (extrapolation).

# iNEXT uses Hill numbers of order q:
#   q = 0: species richness (number of OTUs)
#   q = 1: Shannon diversity (the exponential of the Shannon entropy)
#   q = 2: Simpson diversity (the inverse of the Simpson concentration)
# Larger numbers always mean higher diversity. Hill numbers are expressed in
# "effective number of OTUs": the number of equally abundant OTUs that would
# give the same diversity. As q increases, rare OTUs get less weight.

# iNEXT expects a data frame with OTUs as rows and samples as columns.
# Use the non-rarefied count data: iNEXT does its own rarefaction.
otu_df <- as.data.frame(as(otu_table(otu_BLAST_clean), "matrix"))
if (!taxa_are_rows(otu_BLAST_clean)) otu_df <- as.data.frame(t(otu_df))

# Run iNEXT for each sample individually
#   - endpoint: the largest number of reads to extrapolate to. Extrapolation is
#     only reliable up to about twice the library size of the sample.
#   - knots: the number of points calculated along each curve. Fewer knots make
#     the calculation faster.
inext_res <- iNEXT(otu_df, q = c(0, 1, 2), datatype = "abundance",
                   endpoint = 2 * rare_depth, knots = 20)

# It takes a while to run...

# Plot the rarefaction (solid line) and extrapolation (dashed line) curves.
# The vertical line shows the rarefaction depth.
ggiNEXT(inext_res, type = 1, facet.var = "Order.q", grey = TRUE) +
  geom_vline(xintercept = rare_depth) +
  theme(legend.position = "none")

# For Shannon and Simpson diversity (q = 1 and q = 2), the curves reach an
# asymptote (plateau) well before the rarefaction depth: the rarefied data
# capture the diversity of the samples. For richness (q = 0), the curves keep
# increasing, which means that more OTUs would be found with more sequencing.
# Richness is very sensitive to rare OTUs and to sequencing depth, and is
# therefore the measure that requires the most caution.

# The sample completeness (sample coverage) gives the proportion of the
# community represented by the OTUs detected in each sample. A coverage close
# to 1 means that very few new OTUs would be found with more sequencing.
ggiNEXT(inext_res, type = 2, grey = TRUE) +
  geom_vline(xintercept = rare_depth) +
  theme(legend.position = "none")

# iNEXT also returns the observed and estimated (asymptotic) diversity of each
# sample. The estimated richness is the Chao1 estimator.
inext_div <- inext_res$AsyEst
head(inext_div)


# Alpha diversity measures -----------------------------------------------------

# Alpha diversity can be represented by multiple measures, e.g., richness,
# Shannon, Simpson, Chao1, etc. Here, we calculate them on the rarefied data.

# phyloseq::estimate_richness() calculates several measures at once:
#   - Observed: the number of OTUs (richness)
#   - Shannon: the Shannon entropy (H')
#   - InvSimpson: the inverse Simpson index (= Hill number q = 2)
alpha_div <- estimate_richness(otu_table_rarefy,
                               measures = c("Observed", "Shannon", "InvSimpson"))

# Convert the Shannon entropy to Shannon diversity (= Hill number q = 1), so
# that the three measures are expressed in the same unit (number of OTUs)
alpha_div$Shannon_div <- exp(alpha_div$Shannon)

# Add the estimated richness (Chao1) from iNEXT, calculated on the
# non-rarefied data
rich_est <- inext_div %>%
  filter(Diversity == "Species richness") %>%
  select(Assemblage, Estimator) %>%
  dplyr::rename(Chao1 = Estimator)

alpha_div$Chao1 <- rich_est$Chao1[match(rownames(alpha_div), rich_est$Assemblage)]

# Combine with the sample data
# Note: estimate_richness() returns the samples in the same order as the
# phyloseq object.
alpha_div <- cbind.data.frame(sample_data(otu_table_rarefy)) %>%
  cbind(alpha_div)

# Look at the data
head(alpha_div[, c("SiteID", "Land_use", "Observed", "Chao1",
                   "Shannon_div", "InvSimpson")])

# The measures are often correlated with each other
pairs(alpha_div[, c("Observed", "Chao1", "Shannon_div", "InvSimpson")])


# Compare alpha diversity between sites ----------------------------------------

# We now test if alpha diversity differs between sites using an analysis of
# variance (ANOVA), followed by a Tukey test to find which sites differ from
# each other. The significance letters from the Tukey test (compact letter
# display) are added to the plots: sites that share a letter are not
# significantly different.

# NOTE: Each site has 9 samples. The land use types are not replicated in the
# same way (6 oak forests, but only 1 site for each of the other land uses).
# Therefore, the samples of one site cannot be used as replicates of a land
# use, as this would be pseudoreplication. This is why we compare sites here,
# and colour the plots by land use.

# Order the sites by land use for the plots
alpha_div$SiteID <- factor(alpha_div$SiteID,
                           levels = unique(alpha_div$SiteID[order(alpha_div$Land_use)]))

## Richness ####

# ANOVA
M_rich <- aov(Observed ~ SiteID, data = alpha_div)
anova(M_rich)

# Check the assumptions of the ANOVA: the residuals should be normally
# distributed (Shapiro-Wilk test) and the variances equal between groups
# (Bartlett test). If not, consider a non-parametric test (Kruskal-Wallis).
shapiro.test(residuals(M_rich))
bartlett.test(Observed ~ SiteID, data = alpha_div)

# Tukey test
M_rich_tukey <- TukeyHSD(M_rich)

# Compact letter display
cld_rich <- multcompLetters4(M_rich, M_rich_tukey)

# Table with the mean per site and the significance letters
dt_rich <- group_by(alpha_div, SiteID, Land_use) %>%
  summarise(mean = mean(Observed), sd = sd(Observed), .groups = "drop")
dt_rich$cld <- cld_rich$SiteID$Letters[as.character(dt_rich$SiteID)]

# Plot barplot with the significance letters on top of the bars
ggplot(alpha_div, aes(x = SiteID, y = Observed)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = Land_use)) +
  geom_errorbar(data = dt_rich, aes(y = mean, ymin = mean - sd, ymax = mean + sd),
                width = 0.2) +
  geom_jitter(size = 2, alpha = 0.5, width = 0.1) +
  geom_text(data = dt_rich, aes(label = cld, y = max(alpha_div$Observed) * 1.08)) +
  ylab("Bacterial richness (number of OTUs)") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## Estimated richness (Chao1) ####

# We saw that the richness curves did not reach an asymptote. Let's see if the
# estimated richness shows the same pattern as the observed richness.

# ANOVA
M_chao <- aov(Chao1 ~ SiteID, data = alpha_div)
anova(M_chao)

# Check the assumptions
shapiro.test(residuals(M_chao))
bartlett.test(Chao1 ~ SiteID, data = alpha_div)

# Tukey test and compact letter display
M_chao_tukey <- TukeyHSD(M_chao)
cld_chao <- multcompLetters4(M_chao, M_chao_tukey)

# Table with the mean per site and the significance letters
dt_chao <- group_by(alpha_div, SiteID, Land_use) %>%
  summarise(mean = mean(Chao1), sd = sd(Chao1), .groups = "drop")
dt_chao$cld <- cld_chao$SiteID$Letters[as.character(dt_chao$SiteID)]

# Plot barplot with the significance letters on top of the bars
ggplot(alpha_div, aes(x = SiteID, y = Chao1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = Land_use)) +
  geom_errorbar(data = dt_chao, aes(y = mean, ymin = mean - sd, ymax = mean + sd),
                width = 0.2) +
  geom_jitter(size = 2, alpha = 0.5, width = 0.1) +
  geom_text(data = dt_chao, aes(label = cld, y = max(alpha_div$Chao1) * 1.08)) +
  ylab("Estimated bacterial richness (Chao1)") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## Shannon diversity ####

# ANOVA
M_shan <- aov(Shannon_div ~ SiteID, data = alpha_div)
anova(M_shan)

# Check the assumptions
shapiro.test(residuals(M_shan))
bartlett.test(Shannon_div ~ SiteID, data = alpha_div)

# Tukey test and compact letter display
M_shan_tukey <- TukeyHSD(M_shan)
cld_shan <- multcompLetters4(M_shan, M_shan_tukey)

# Table with the mean per site and the significance letters
dt_shan <- group_by(alpha_div, SiteID, Land_use) %>%
  summarise(mean = mean(Shannon_div), sd = sd(Shannon_div), .groups = "drop")
dt_shan$cld <- cld_shan$SiteID$Letters[as.character(dt_shan$SiteID)]

# Plot barplot with the significance letters on top of the bars
ggplot(alpha_div, aes(x = SiteID, y = Shannon_div)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = Land_use)) +
  geom_errorbar(data = dt_shan, aes(y = mean, ymin = mean - sd, ymax = mean + sd),
                width = 0.2) +
  geom_jitter(size = 2, alpha = 0.5, width = 0.1) +
  geom_text(data = dt_shan, aes(label = cld, y = max(alpha_div$Shannon_div) * 1.08)) +
  ylab("Bacterial Shannon diversity (q = 1)") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## Simpson diversity ####

# ANOVA
M_simp <- aov(InvSimpson ~ SiteID, data = alpha_div)
anova(M_simp)

# Check the assumptions
shapiro.test(residuals(M_simp))
bartlett.test(InvSimpson ~ SiteID, data = alpha_div)

# Tukey test and compact letter display
M_simp_tukey <- TukeyHSD(M_simp)
cld_simp <- multcompLetters4(M_simp, M_simp_tukey)

# Table with the mean per site and the significance letters
dt_simp <- group_by(alpha_div, SiteID, Land_use) %>%
  summarise(mean = mean(InvSimpson), sd = sd(InvSimpson), .groups = "drop")
dt_simp$cld <- cld_simp$SiteID$Letters[as.character(dt_simp$SiteID)]

# Plot barplot with the significance letters on top of the bars
ggplot(alpha_div, aes(x = SiteID, y = InvSimpson)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = Land_use)) +
  geom_errorbar(data = dt_simp, aes(y = mean, ymin = mean - sd, ymax = mean + sd),
                width = 0.2) +
  geom_jitter(size = 2, alpha = 0.5, width = 0.1) +
  geom_text(data = dt_simp, aes(label = cld, y = max(alpha_div$InvSimpson) * 1.08)) +
  ylab("Bacterial Simpson diversity (q = 2)") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Compare the patterns of the four measures: richness gives the same weight to
# all OTUs, whereas Shannon and Simpson diversity give more weight to the
# abundant OTUs. If a site has a high richness but a low Simpson diversity, its
# community contains many rare OTUs but is dominated by a few abundant ones.


# Save alpha diversity table ---------------------------------------------------

save(alpha_div, file = "RData/16S/alpha_div.RData")
