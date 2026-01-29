#==============================================================================#
#                                                                              #
#  Global filtering: removing singletons, rare and/or scarce OTUs              #
#                                                                              #
#  Written by Karelle Rheault: karh@ign.ku.dk                                  #
#                                                                              #
#==============================================================================#

# Set working directory (modify the path below to match your own)
setwd("some_folder_path/my_working_directory")

# Load packages
source("scripts/001_required_packages.R")

# Load the cleaned OTU table from previous scripts
# "003_negative_control.R" or "003_positive_control.R"
load("RData/16S/otu_ctrl_neg_clean.RData")
load("RData/16S/otu_ctrl_pos_clean.RData")

# Rename your OTU table to otu_table_raw to match this workflow
otu_table_raw <- otu_ctrl_pos_clean

# Singletons -------------------------------------------------------------------

# Count OTUs before filtering
ntaxa(otu_table_raw) # 25390

# Remove singletons
otu_table_single <- prune_taxa(taxa_sums(otu_table_raw) > 1, otu_table_raw)

# Count OTUs after filtering
ntaxa(otu_table_single) # 24923 (~2% of the OTUs were singletons)

# Rare OTUs --------------------------------------------------------------------

# OTUs are considered rare if their abundance is small compared to other OTUs.
# They may be present in many samples, but if their total sum across samples is
# smaller than k% of the average sequencing depth (relative abundance), or
# smaller than k reads (absolute abundance) then they are considered rare.

# NOTE: When calculating the average sequencing depth, ignore controls and
# blanks, as they artificially lower the average.

## Relative abundance ####

# Calculate the average sequencing depth
seq_depth_avg <- mean(sample_sums(otu_table_raw)) # 150745.2 reads

# Calculate thresholds based on k% of the average sequencing depth
# For the sake of this exercise, multiple thresholds between 0.005% and 1% are
# tested to compare their effects on various community analyses.
thresh005 <- ceiling(seq_depth_avg*0.005/100) # 8 reads
thresh01 <- ceiling(seq_depth_avg*0.01/100) # 16 reads
thresh05 <- ceiling(seq_depth_avg*0.05/100) # 76 reads
thresh1 <- ceiling(seq_depth_avg*0.1/100) # 151 reads
thresh10 <- ceiling(seq_depth_avg*1/100) # 1508 reads

# Remove all OTUs that have a total sum lower than the threshold
otu_table_rare005 <- otu_table_raw %>% prune_taxa(taxa_sums(.) > thresh005, .)
otu_table_rare01 <- otu_table_raw %>% prune_taxa(taxa_sums(.) > thresh01, .)
otu_table_rare05 <- otu_table_raw %>% prune_taxa(taxa_sums(.) > thresh05, .)
otu_table_rare1 <- otu_table_raw %>% prune_taxa(taxa_sums(.) > thresh1, .)
otu_table_rare10 <- otu_table_raw %>% prune_taxa(taxa_sums(.) > thresh10, .)

# Compare the number of OTU in each filtering step
ntaxa(otu_table_raw) #25390 (before filtering)
ntaxa(otu_table_single) #24923 (~2% singletons)
ntaxa(otu_table_rare005) #8638 (~66% of OTUs have less than 8 reads; 0.005%)
ntaxa(otu_table_rare01) #6631 (~74% of OTUs have less than 16 reads; 0.01%)
ntaxa(otu_table_rare05) #3783 (~85% of OTUs have less than 76 reads; 0.05%)
ntaxa(otu_table_rare1) #2929 (~88% of OTUs have less than 151 reads; 0.1%)
ntaxa(otu_table_rare10) #1042 (~96% of OTUs have less than 1508 reads; 1%)

### Compare thresholds ####

# To visualize the effect of these different thresholds, we can create a data
# frame per filtered OTU table that summarizes the sequencing depth per sample
# and the Observed OTU (richness).

# Non-filtered OTU table
seq_depth_raw <- cbind.data.frame(sample_data(otu_table_raw)) %>%
  cbind(OTU_rich = estimate_richness(otu_table_raw, measures = "Observed")) %>%
  mutate(LibrarySize = sample_sums(otu_table_raw)) %>%
  arrange(LibrarySize) %>% mutate(Index = row_number()) %>%
  mutate(seq_depth = "raw")

# Singletons
seq_depth_single <- cbind.data.frame(sample_data(otu_table_single)) %>%
  cbind(OTU_rich = estimate_richness(otu_table_single, measures = "Observed")) %>%
  mutate(LibrarySize = sample_sums(otu_table_single)) %>%
  arrange(LibrarySize) %>% mutate(Index = row_number()) %>%
  mutate(seq_depth = "singletons")

# Rare 0.005% (<8 reads)
seq_depth_rare005 <- cbind.data.frame(sample_data(otu_table_rare005)) %>%
  cbind(OTU_rich = estimate_richness(otu_table_rare005, measures = "Observed")) %>%
  mutate(LibrarySize = sample_sums(otu_table_rare005)) %>%
  arrange(LibrarySize) %>% mutate(Index = row_number()) %>%
  mutate(seq_depth = "rare 0.005%")

# Rare 0.01% (<16 reads)
seq_depth_rare01 <- cbind.data.frame(sample_data(otu_table_rare01)) %>%
  cbind(OTU_rich = estimate_richness(otu_table_rare01, measures = "Observed")) %>%
  mutate(LibrarySize = sample_sums(otu_table_rare01)) %>%
  arrange(LibrarySize) %>% mutate(Index = row_number()) %>%
  mutate(seq_depth = "rare 0.01%")

# Rare 0.05% (<76 reads)
seq_depth_rare05 <- cbind.data.frame(sample_data(otu_table_rare05)) %>%
  cbind(OTU_rich = estimate_richness(otu_table_rare05, measures = "Observed")) %>%
  mutate(LibrarySize = sample_sums(otu_table_rare05)) %>%
  arrange(LibrarySize) %>% mutate(Index = row_number()) %>%
  mutate(seq_depth = "rare 0.05%")

# Rare 0.1% (<151 reads)
seq_depth_rare1 <- cbind.data.frame(sample_data(otu_table_rare1)) %>%
  cbind(OTU_rich = estimate_richness(otu_table_rare1, measures = "Observed")) %>%
  mutate(LibrarySize = sample_sums(otu_table_rare1)) %>%
  arrange(LibrarySize) %>% mutate(Index = row_number()) %>%
  mutate(seq_depth = "rare 0.1%")

# Rare 1% (<1508 reads)
seq_depth_rare10 <- cbind.data.frame(sample_data(otu_table_rare10)) %>%
  cbind(OTU_rich = estimate_richness(otu_table_rare10, measures = "Observed")) %>%
  mutate(LibrarySize = sample_sums(otu_table_rare10)) %>%
  arrange(LibrarySize) %>% mutate(Index = row_number()) %>%
  mutate(seq_depth = "rare 1%")

# Combine data frames
seq_depth_all <- seq_depth_raw %>% rbind(seq_depth_single) %>%
  rbind(seq_depth_rare005) %>% rbind(seq_depth_rare01) %>%
  rbind(seq_depth_rare05) %>% rbind(seq_depth_rare1) %>%
  rbind(seq_depth_rare10)

#### Effect on LibrarySize ####

# Plot the library size of all samples
ggplot(data = seq_depth_all, aes(x = Index, y = LibrarySize, color = seq_depth)) +
  geom_point() + geom_line()

# Zoom in on the large library size part
seq_depth_all %>% filter(Index > 85) %>%
  ggplot(aes(x = Index, y = LibrarySize, color = seq_depth)) +
  geom_point() + geom_line()

# With this analysis, we can see that the first filtering step (singletons)
# does not affect LibrarySize, as the points completely overlap (purple vs pink).

# In addition, the filtering of the rare OTUs at a threshold of 0.01% or less
# changes very little in the LibrarySize, even though we lost ~73% of OTUs

#### Effect on OTU richness ####

# To confirm our choice of the best filtering method, let's compare the effect
# of filtering on OTU richness.

# Plot the OTU richness of all samples
ggplot(data = seq_depth_all, aes(x = Index, y = Observed, color = seq_depth)) +
  geom_point() + geom_line()

# This plot shows that these different filtering thresholds return the same
# pattern of OTU richness as a function of library size.

#### Effect on beta diversity ####

# To confirm our choice of filtering further, we can compare how the different
# filtering methods affect beta diversity

# Non-filtered OTU table
beta_div_raw <- otu_table_raw %>%
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

# Rare 0.005%
beta_div_rare005 <- otu_table_rare005 %>%
  tax_transform("identity", rank = "taxa_OTU") %>%
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "Land_use", auto_caption = NA) +
  stat_ellipse(aes(colour = Land_use))

# Rare 0.01%
beta_div_rare01 <- otu_table_rare01 %>%
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
ggarrange(beta_div_raw, beta_div_single, beta_div_rare005, beta_div_rare01,
          beta_div_rare05, beta_div_rare1, beta_div_rare10,
          labels = c("Non-filtered (25390 OTUs)", "Singletons (24923 OTUs)",
                     "0.005% (<8 reads; 9147 OTUs)", "0.01% (<16 reads; 6796 OTUs)",
                     "0.05% (<76 reads; 3798 OTUs)", "0.1% (<151 reads; 2934 OTUs)",
                     "1% (<1508 reads; 1042 OTUs)"),
          ncol = 3, nrow = 3, common.legend = TRUE,
          hjust = 0, label.x = 0.2, font.label = list(size = 10))

# With this analysis, we can see that removing rare OTUs reduces the variance
# within treatment a little bit, but does not affect overall patterns.

# Overall, all three plots showed that rare OTUs were mostly found in samples
# with large library size, and filtering out rare OTUs does not affect patterns.
# Only a threshold as high as 1% considerably affected patterns.

# Select the best filtering threshold
otu_table_rare <- otu_table_rare005

# Scarce OTUs ------------------------------------------------------------------

# After filtering out rare OTUs, we can look at scarce OTUs (low-prevalence).
# Any OTU that occur in less than 3 samples (number of replicates) could be
# considered scarce, but it is important to very that they are not abundant
# before removing them.

# Convert OTU table to presence-absence data and relative abundance
otu_table_rare.pa <- transform_sample_counts(otu_table_rare, function(x) 1*(x>0))

# Calculate prevalence and total read counts of each OTUs
prevalence <- data.frame(prevalence = taxa_sums(otu_table_rare.pa),
                         abundance = taxa_sums(otu_table_rare)) %>%
  # Filter OTUs with low prevalence (occur in less than 3 samples)
  # 3 samples is the number of replicates used in this experiment
  filter(prevalence < 3)

# 2177 OTUs were found in less than 3 samples. However, important OTUs could be
# removed this way. Some OTUs may be found in few samples, but still be abundant.
# It is important to verify their abundance before removing them.

# We already removed anything with less than 8 reads total, but we could increase
# this threshold to remove scarce OTUs. To decide whether an OTU is important or
# not, one can use the same method as used for filtering out the rare OTUs.

# Identify scarce OTUs filtered out for each threshold
scarce01 <- prevalence %>% filter(abundance < thresh01) #1243 OTUs
scarce05 <- prevalence %>% filter(abundance < thresh05) #2117 OTUs
scarce1 <- prevalence %>% filter(abundance < thresh1) #2152 OTUs
scarce10 <- prevalence %>% filter(abundance < thresh10) #2176 OTUs

# Remove all scarce OTUs that have a total sum lower than the threshold
otu_table_scarce01 <- remove_taxa(rownames(scarce01), otu_table_rare)
otu_table_scarce05 <- remove_taxa(rownames(scarce05), otu_table_rare)
otu_table_scarce1 <- remove_taxa(rownames(scarce1), otu_table_rare)
otu_table_scarce10 <- remove_taxa(rownames(scarce10), otu_table_rare)

# Compare the number of OTU in each filtering step
ntaxa(otu_table_raw) #25390 (before filtering)
ntaxa(otu_table_single) #24923 (~2% singletons)
ntaxa(otu_table_rare) #8638 (~66% of OTUs have less than 8 reads; 0.005%)
ntaxa(otu_table_scarce01) #7395 (+ ~4.9% of scarce OTUs with <16 reads; 0.01%)
ntaxa(otu_table_scarce05) #6521 (+ ~8.3% of scarce OTUs with <76 reads; 0.05%)
ntaxa(otu_table_scarce1) #6486 (+ ~8.5% of scarce OTUs with <151 reads; 0.1%)
ntaxa(otu_table_scarce10) #6462 (+ ~8.6% of scarce OTUs with <1508 reads; 1%)

# Compared to previous methods that were filtering out rare OTUs based on read
# counts alone, the second highest threshold (0.01% abundance / <16 reads total)
# removed ~74% of OTUs and we could barely see any significant effect on library
# size, alpha diversity and beta diversity. Here, even at the highest threshold,
# only a total of ~75% of OTUs are removed. We can therefore assume that the
# patterns would remain intact (or comparable to the effects seen at a threshold
# of 0.01% / <16 reads).

# Furthermore, only one OTU had a higher read count than the highest threshold
# at 1% (1508 reads), which had 1894 reads (~1.3%) across two samples. Therefore,
# all identified scarce OTUs can be considered rare and removing them could
# reduce noise in the data. For the sake of staying conservative, we choose to
# remove only scarce OTUs with an abundance <0.01% (16 reads).

otu_table_scarce <- otu_table_scarce01

# Save cleaned OTU table -------------------------------------------------------

# OTU table with singletons filtered
save(otu_table_single, file = "RData/16S/otu_single_clean.RData")

# OTU table with rare taxa filtered
save(otu_table_rare, file = "RData/16S/otu_rare_clean.RData")

# OTU table with scarce taxa filtered
save(otu_table_scarce, file = "RData/16S/otu_scarce_clean.RData")
