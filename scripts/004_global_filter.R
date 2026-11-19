#==============================================================================#
#                                                                              #
#  Global filtering: removing singletons, rare and/or scarce OTUs              #
#                                                                              #
#  Written by Karelle Rheault: karh@ign.ku.dk                                  #
#                                                                              #
#==============================================================================#

# Set working directory (modify the path below to match your own)
setwd("some_folder_path/my_working_directory")

# Load packages and homemade functions
source("scripts/001_required_packages.R")

# Load the cleaned OTU table from previous scripts
load("RData/ITS/otu_ctrl_clean.RData")

# Taxonomy filter --------------------------------------------------------------

# This first simple step serves to preserve only the relevant OTUs that were
# targeted with each primer sets: i.e., Bacteria for 16S, Fungi for ITS and
# Metazoa for CO1.

# Select organisms that should be included in the analysis
otu_table_tax <- subset_taxa(otu_ctrl_clean, kingdom == "Fungi")

# Count OTUs before and after filtering
ntaxa(otu_ctrl_clean) # 16687
ntaxa(otu_table_tax) # 9043 (~46% of the OTUs were not Fungi)

# The ITS primers also amplify many other eukaryotes (e.g., Alveolata, Rhizaria,
# Amoebozoa, Metazoa and plants), as well as the synthetic spike-in of the
# positive controls (kingdom "Unispike1"). These are all removed here.

# How many reads were removed?
sum(sample_sums(otu_ctrl_clean)) - sum(sample_sums(otu_table_tax)) # 69358 reads
sum(sample_sums(otu_table_tax)) # 360336 reads belong to Fungi

# Note that each genetic marker needs its own taxonomy filter. For example, the
# 16S primers also amplify chloroplasts and mitochondria, which contain DNA
# similar to bacteria. In the databases used to BLAST the 16S marker,
# chloroplasts are assigned at the order level and mitochondria at the family
# level, both within the kingdom Bacteria, and must be removed separately:
#   subset_taxa(otu_table_tax, order != "Chloroplast" & family != "Mitochondria")

# Explore the abundance and prevalence of OTUs ---------------------------------

# Sequencing datasets are characterized by a long tail of rare OTUs, i.e., OTUs
# represented by a low number of reads. Some of them may be products of
# sequencing errors that were not removed by the bioinformatic pipeline. Are
# they real? Is their presence meaningful in our analysis? It is often hard to
# decide which cutoff to use for an OTU to be considered "real", and the choice
# is often a bit arbitrary. It is therefore important to explore the data
# before deciding.

## Number of reads per OTU ####

# Some OTUs have 0 reads at this stage: they were only present in the control
# samples, which were removed in 003_negative_control.R. Exclude them from the
# exploration.
otu_table_explore <- prune_taxa(taxa_sums(otu_table_tax) > 0, otu_table_tax)

# Total number of reads per OTU
reads_otu <- data.frame(reads = taxa_sums(otu_table_explore))

min(reads_otu$reads) # 1 read
max(reads_otu$reads) # 16793 reads

# Plot the distribution of the number of reads per OTU (log scale)
p_reads_otu <- ggplot(reads_otu, aes(x = reads)) +
  geom_histogram(bins = 20, fill = "steelblue") +
  scale_x_log10(labels = scales::comma) +
  labs(x = "Number of reads per OTU (log10 scale)", y = "Number of OTUs") +
  theme_minimal()
p_reads_otu
ggsave("figures/004_reads_per_otu_histogram.png", p_reads_otu, width = 7, height = 5)

# How many OTUs have more than 1000 reads?
sum(reads_otu$reads > 1000) # 58 OTUs

# Most OTUs have a low number of reads, and only a few OTUs are very abundant.

## Prevalence of OTUs ####

# The prevalence of an OTU is the number of samples in which it is present.

# Convert the OTU table to presence-absence (1 = present, 0 = absent)
otu_table_pa <- transform_sample_counts(otu_table_explore, function(x) 1 * (x > 0))

# Number of OTUs per sample
sample_sums(otu_table_pa)

# Number of samples in which each OTU is present (prevalence)
prev_otu <- data.frame(prevalence = taxa_sums(otu_table_pa))

# How many OTUs are present in only one or two samples?
sum(prev_otu$prevalence == 1) # 4843 OTUs
sum(prev_otu$prevalence == 2) # 1367 OTUs

# Plot the distribution of the prevalence of OTUs
p_prev_otu <- ggplot(prev_otu, aes(x = prevalence)) +
  geom_histogram(bins = 60, fill = "steelblue") +
  labs(x = "Prevalence of OTUs (number of samples)", y = "Number of OTUs") +
  theme_minimal()
p_prev_otu
ggsave("figures/004_otu_prevalence_histogram.png", p_prev_otu, width = 7, height = 5)

# Most OTUs occur in a very low number of samples. Filtering OTUs based on
# their prevalence would therefore lead to the loss of a large number of OTUs.

# What if we remove the OTUs present in less than 10 samples?
prev_otu10 <- filter(prev_otu, prevalence >= 10)

p_prev_otu10 <- ggplot(prev_otu10, aes(x = prevalence)) +
  geom_histogram(bins = 60, fill = "steelblue") +
  labs(x = "Prevalence of OTUs (number of samples)", y = "Number of OTUs") +
  theme_minimal()
p_prev_otu10
ggsave("figures/004_otu_prevalence_histogram_min10.png", p_prev_otu10,
       width = 7, height = 5)

# Even after removing the OTUs present in less than 10 samples, the same
# pattern is observed: OTUs with a low prevalence are the most frequent. It is
# hard to choose a cutoff based on prevalence only. Let's look at filters
# based on the abundance of OTUs instead (see below).


# Singletons -------------------------------------------------------------------

# Remove singletons
otu_table_single <- prune_taxa(taxa_sums(otu_table_tax) > 1, otu_table_tax)

# Count OTUs after filtering
ntaxa(otu_table_single) # 5825 (~36% of the fungal OTUs were singletons)

# In this dataset, a large proportion of the OTUs are singletons (OTUs with a
# single read in the whole dataset). Together, however, they represent less
# than 1% of the reads.

# Rare OTUs --------------------------------------------------------------------

# OTUs are considered rare if their abundance is small compared to other OTUs.
# They may be present in many samples, but if their total sum across samples is
# smaller than k% of the average sequencing depth (relative abundance), or
# smaller than k reads (absolute abundance) then they are considered rare.

# NOTE: When calculating the average sequencing depth, ignore controls and
# blanks, as they artificially lower the average.

## Relative abundance ####

# Calculate the average sequencing depth
seq_depth_avg <- mean(sample_sums(otu_table_tax)) # 4048.7 reads

# Calculate thresholds based on k% of the average sequencing depth
# For the sake of this exercise, multiple thresholds between 0.05% and 1% are
# tested to compare their effects on various community analyses.
thresh05 <- ceiling(seq_depth_avg*0.05/100) # 3 reads
thresh1 <- ceiling(seq_depth_avg*0.1/100) # 5 reads
thresh10 <- ceiling(seq_depth_avg*1/100) # 41 reads

# Note: because the sequencing depth of this ITS dataset is low (~4000 reads
# per sample), lower thresholds (e.g., 0.005% or 0.01%) correspond to 1 read,
# i.e., they are identical to removing singletons. With a deeper sequencing
# (e.g., ~140000 reads per sample, typical of 16S datasets), the same
# percentages correspond to more reads (7 and 14 reads, respectively).

# Remove all OTUs that have a total sum lower than the threshold
otu_table_rare05 <- otu_table_tax %>% prune_taxa(taxa_sums(.) > thresh05, .)
otu_table_rare1 <- otu_table_tax %>% prune_taxa(taxa_sums(.) > thresh1, .)
otu_table_rare10 <- otu_table_tax %>% prune_taxa(taxa_sums(.) > thresh10, .)

# Compare the number of OTU in each filtering step
ntaxa(otu_ctrl_clean) #16687 (before filtering)
ntaxa(otu_table_tax) #9043 (~54% of OTUs are Fungi)
ntaxa(otu_table_single) #5825 (~36% singletons)
ntaxa(otu_table_rare05) #3597 (~60% of OTUs have 3 reads or less; 0.05%)
ntaxa(otu_table_rare1) #2763 (~69% of OTUs have 5 reads or less; 0.1%)
ntaxa(otu_table_rare10) #787 (~91% of OTUs have 41 reads or less; 1%)

### Compare thresholds ####

# To visualize the effect of these different thresholds, we can create a data
# frame per filtered OTU table that summarizes the sequencing depth per sample
# and the Observed OTU (richness).

# Non-filtered OTU table
seq_depth_raw <- cbind.data.frame(sample_data(otu_table_tax)) %>%
  cbind(estimate_richness(otu_table_tax, measures = "Observed")) %>%
  mutate(LibrarySize = sample_sums(otu_table_tax)) %>%
  arrange(LibrarySize) %>% mutate(Index = row_number()) %>%
  mutate(seq_depth = "raw")

# Singletons
seq_depth_single <- cbind.data.frame(sample_data(otu_table_single)) %>%
  cbind(estimate_richness(otu_table_single, measures = "Observed")) %>%
  mutate(LibrarySize = sample_sums(otu_table_single)) %>%
  arrange(LibrarySize) %>% mutate(Index = row_number()) %>%
  mutate(seq_depth = "singletons")

# Rare 0.05% (<=3 reads)
seq_depth_rare05 <- cbind.data.frame(sample_data(otu_table_rare05)) %>%
  cbind(estimate_richness(otu_table_rare05, measures = "Observed")) %>%
  mutate(LibrarySize = sample_sums(otu_table_rare05)) %>%
  arrange(LibrarySize) %>% mutate(Index = row_number()) %>%
  mutate(seq_depth = "rare 0.05%")

# Rare 0.1% (<=5 reads)
seq_depth_rare1 <- cbind.data.frame(sample_data(otu_table_rare1)) %>%
  cbind(estimate_richness(otu_table_rare1, measures = "Observed")) %>%
  mutate(LibrarySize = sample_sums(otu_table_rare1)) %>%
  arrange(LibrarySize) %>% mutate(Index = row_number()) %>%
  mutate(seq_depth = "rare 0.1%")

# Rare 1% (<=41 reads)
seq_depth_rare10 <- cbind.data.frame(sample_data(otu_table_rare10)) %>%
  cbind(estimate_richness(otu_table_rare10, measures = "Observed")) %>%
  mutate(LibrarySize = sample_sums(otu_table_rare10)) %>%
  arrange(LibrarySize) %>% mutate(Index = row_number()) %>%
  mutate(seq_depth = "rare 1%")

# Combine data frames
seq_depth_all <- seq_depth_raw %>% rbind(seq_depth_single) %>%
  rbind(seq_depth_rare05) %>% rbind(seq_depth_rare1) %>%
  rbind(seq_depth_rare10)

#### Effect on LibrarySize ####

# Plot the library size of all samples
p_libsize <- ggplot(data = seq_depth_all, aes(x = Index, y = LibrarySize, color = seq_depth)) +
  geom_point() + geom_line()
p_libsize
ggsave("figures/004_filter_library_size.png", p_libsize, width = 7, height = 5)

# Zoom in on the large library size part
p_libsize_zoom <- seq_depth_all %>% filter(Index > 80) %>%
  ggplot(aes(x = Index, y = LibrarySize, color = seq_depth)) +
  geom_point() + geom_line()
p_libsize_zoom
ggsave("figures/004_filter_library_size_zoom.png", p_libsize_zoom, width = 7, height = 5)

# With this analysis, we can see that removing singletons barely affects
# LibrarySize (pink vs blue), even though we lost ~36% of the OTUs.

# In addition, the filtering of the rare OTUs at a threshold of 0.1% or less
# (5 reads or less) changes very little in the LibrarySize, even though we lost
# ~69% of OTUs. Only the threshold of 1% removes a substantial number of reads.

#### Effect on OTU richness ####

# To confirm our choice of the best filtering method, let's compare the effect
# of filtering on OTU richness.

# Plot the OTU richness of all samples
p_richness <- ggplot(data = seq_depth_all, aes(x = Index, y = Observed, color = seq_depth)) +
  geom_point() + geom_line()
p_richness
ggsave("figures/004_filter_richness.png", p_richness, width = 7, height = 5)

# This plot shows that every filtering step lowers the OTU richness of the
# samples, especially in the samples with a large library size. The thresholds
# up to 0.1% return the same pattern of OTU richness as a function of library
# size, whereas the threshold of 1% flattens the differences between samples.

#### Effect on beta diversity ####

# To confirm our choice of filtering further, we can compare how the different
# filtering methods affect beta diversity

# Non-filtered OTU table
beta_div_raw <- otu_table_tax %>%
  tax_transform("identity", rank = "taxa_OTU") %>%
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "Land_use", auto_caption = NA) +
  stat_ellipse(aes(colour = Land_use))

# Singletons
beta_div_single <- otu_table_single %>%
  tax_transform("identity", rank = "taxa_OTU") %>%
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "Land_use", auto_caption = NA) +
  stat_ellipse(aes(colour = Land_use))

# Rare 0.05%
beta_div_rare05 <- otu_table_rare05 %>%
  tax_transform("identity", rank = "taxa_OTU") %>%
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "Land_use", auto_caption = NA) +
  stat_ellipse(aes(colour = Land_use))

# Rare 0.1%
beta_div_rare1 <- otu_table_rare1 %>%
  tax_transform("identity", rank = "taxa_OTU") %>%
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "Land_use", auto_caption = NA) +
  stat_ellipse(aes(colour = Land_use))

# Rare 1%
beta_div_rare10 <- otu_table_rare10 %>%
  tax_transform("identity", rank = "taxa_OTU") %>%
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "Land_use", auto_caption = NA) +
  stat_ellipse(aes(colour = Land_use))

# Visualize all plots together
p_beta_div_all <- ggarrange(
  beta_div_raw, beta_div_single, beta_div_rare05, beta_div_rare1,
  beta_div_rare10,
  labels = c(
    "Non-filtered (9043 OTUs)", "Singletons (5825 OTUs)",
    "0.05% (<=3 reads; 3597 OTUs)", "0.1% (<=5 reads; 2763 OTUs)",
    "1% (<=41 reads; 787 OTUs)"),
  ncol = 3, nrow = 2, common.legend = TRUE,
  hjust = 0, label.x = 0.2, font.label = list(size = 10))
p_beta_div_all
ggsave("figures/004_filter_thresholds_pcoa.png", p_beta_div_all,
       width = 14, height = 9, bg = "white")

# With this analysis, we can see that the overall patterns are the same with
# all filtering thresholds: the agricultural sites and the grassland separate
# from the oak forests. Note that the orientation of the axes is arbitrary: in
# the plot without singletons, the second axis is simply flipped.

# Overall, the three plots showed that rare OTUs were mostly found in samples
# with a large library size, and that filtering out rare OTUs does not affect
# the main patterns. Only a threshold as high as 1% considerably affected the
# richness of the samples.

## Choose a filtering threshold ####

# Is it then necessary to filter rare OTUs? It is hard to choose a cutoff, and
# the choice is often arbitrary. Removing rare OTUs mainly lowers the richness
# of the samples, and has little effect on the community composition.

# In this workflow, we keep all fungal OTUs (no abundance filter), for two
# reasons:
#   1. Rare OTUs, and singletons in particular, are needed to estimate the
#      number of undetected OTUs in each sample (extrapolated richness, see
#      010_alpha_diversity.R). Removing them would bias these estimates.
#   2. The sequencing depth of this dataset is low, so that even the lowest
#      thresholds remove a large proportion of the OTUs.
# If you decide to filter rare OTUs in your own data, report the threshold you
# used, and check that it does not change your conclusions.

# Keep all fungal OTUs, but remove the OTUs with 0 reads (only present in the
# controls, see above)
otu_table_rare <- prune_taxa(taxa_sums(otu_table_tax) > 0, otu_table_tax)
ntaxa(otu_table_rare) # 9039 OTUs

# To use a filtering threshold instead, replace the line above with e.g.:
#   otu_table_rare <- otu_table_single   # remove singletons
#   otu_table_rare <- otu_table_rare05   # remove OTUs with 3 reads or less


# Save cleaned OTU table -------------------------------------------------------

# OTU table with taxonomy filtered
save(otu_table_tax, file = "RData/ITS/otu_taxa_clean.RData")

# OTU table with singletons filtered
save(otu_table_single, file = "RData/ITS/otu_single_clean.RData")

# OTU table with rare taxa filtered
save(otu_table_rare, file = "RData/ITS/otu_rare_clean.RData")

