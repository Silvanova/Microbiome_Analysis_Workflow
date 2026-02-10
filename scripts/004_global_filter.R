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
load("RData/16S/otu_ctrl_clean.RData")

# Taxonomy filter --------------------------------------------------------------

# This first simple step serves to preserve only the relevant OTUs that were
# targeted with each primer sets: i.e., Bacteria for 16S, Fungi for ITS and
# Metazoa for CO1.

# Select organisms that should be included in the analysis
otu_table_tax <- subset_taxa(otu_ctrl_clean, kingdom == "Bacteria")

# Remove organisms that should not be included in the analysis
otu_table_tax <- subset_taxa(otu_table_tax, order != "Chloroplast" &
                               family != "Mitochondria")

# Count OTUs before and after filtering
ntaxa(otu_ctrl_clean) # 25390
ntaxa(otu_table_tax) # 22046 (~13% of the OTUs were not Bacteria)

# Singletons -------------------------------------------------------------------

# Remove singletons
otu_table_single <- prune_taxa(taxa_sums(otu_table_tax) > 1, otu_table_tax)

# Count OTUs after filtering
ntaxa(otu_table_single) # 21632 (~2% of the remaining OTUs were singletons)

# Rare OTUs --------------------------------------------------------------------

# OTUs are considered rare if their abundance is small compared to other OTUs.
# They may be present in many samples, but if their total sum across samples is
# smaller than k% of the average sequencing depth (relative abundance), or
# smaller than k reads (absolute abundance) then they are considered rare.

# NOTE: When calculating the average sequencing depth, ignore controls and
# blanks, as they artificially lower the average.

## Relative abundance ####

# Calculate the average sequencing depth
seq_depth_avg <- mean(sample_sums(otu_table_tax)) # 139853.2 reads

# Calculate thresholds based on k% of the average sequencing depth
# For the sake of this exercise, multiple thresholds between 0.005% and 1% are
# tested to compare their effects on various community analyses.
thresh005 <- ceiling(seq_depth_avg*0.005/100) # 7 reads
thresh01 <- ceiling(seq_depth_avg*0.01/100) # 14 reads
thresh05 <- ceiling(seq_depth_avg*0.05/100) # 70 reads
thresh1 <- ceiling(seq_depth_avg*0.1/100) # 140 reads
thresh10 <- ceiling(seq_depth_avg*1/100) # 1399 reads

# Remove all OTUs that have a total sum lower than the threshold
otu_table_rare005 <- otu_table_tax %>% prune_taxa(taxa_sums(.) > thresh005, .)
otu_table_rare01 <- otu_table_tax %>% prune_taxa(taxa_sums(.) > thresh01, .)
otu_table_rare05 <- otu_table_tax %>% prune_taxa(taxa_sums(.) > thresh05, .)
otu_table_rare1 <- otu_table_tax %>% prune_taxa(taxa_sums(.) > thresh1, .)
otu_table_rare10 <- otu_table_tax %>% prune_taxa(taxa_sums(.) > thresh10, .)

# Compare the number of OTU in each filtering step
ntaxa(otu_ctrl_clean) #25390 (before filtering)
ntaxa(otu_table_tax) #22046 (~87% of OTUs are Bacteria)
ntaxa(otu_table_single) #21632 (~2% singletons)
ntaxa(otu_table_rare005) #8064 (~63% of OTUs have less than 7 reads; 0.005%)
ntaxa(otu_table_rare01) #6183 (~72% of OTUs have less than 14 reads; 0.01%)
ntaxa(otu_table_rare05) #3467 (~84% of OTUs have less than 70 reads; 0.05%)
ntaxa(otu_table_rare1) #2710 (~88% of OTUs have less than 140 reads; 0.1%)
ntaxa(otu_table_rare10) #979 (~96% of OTUs have less than 1399 reads; 1%)

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

# Rare 0.005% (<7 reads)
seq_depth_rare005 <- cbind.data.frame(sample_data(otu_table_rare005)) %>%
  cbind(estimate_richness(otu_table_rare005, measures = "Observed")) %>%
  mutate(LibrarySize = sample_sums(otu_table_rare005)) %>%
  arrange(LibrarySize) %>% mutate(Index = row_number()) %>%
  mutate(seq_depth = "rare 0.005%")

# Rare 0.01% (<14 reads)
seq_depth_rare01 <- cbind.data.frame(sample_data(otu_table_rare01)) %>%
  cbind(estimate_richness(otu_table_rare01, measures = "Observed")) %>%
  mutate(LibrarySize = sample_sums(otu_table_rare01)) %>%
  arrange(LibrarySize) %>% mutate(Index = row_number()) %>%
  mutate(seq_depth = "rare 0.01%")

# Rare 0.05% (<70 reads)
seq_depth_rare05 <- cbind.data.frame(sample_data(otu_table_rare05)) %>%
  cbind(estimate_richness(otu_table_rare05, measures = "Observed")) %>%
  mutate(LibrarySize = sample_sums(otu_table_rare05)) %>%
  arrange(LibrarySize) %>% mutate(Index = row_number()) %>%
  mutate(seq_depth = "rare 0.05%")

# Rare 0.1% (<140 reads)
seq_depth_rare1 <- cbind.data.frame(sample_data(otu_table_rare1)) %>%
  cbind(estimate_richness(otu_table_rare1, measures = "Observed")) %>%
  mutate(LibrarySize = sample_sums(otu_table_rare1)) %>%
  arrange(LibrarySize) %>% mutate(Index = row_number()) %>%
  mutate(seq_depth = "rare 0.1%")

# Rare 1% (<1399 reads)
seq_depth_rare10 <- cbind.data.frame(sample_data(otu_table_rare10)) %>%
  cbind(estimate_richness(otu_table_rare10, measures = "Observed")) %>%
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
# changes very little in the LibrarySize, even though we lost ~72% of OTUs

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
          labels = c(
            "Non-filtered (22046 OTUs)", "Singletons (21632 OTUs)",
            "0.005% (<7 reads; 8064 OTUs)", "0.01% (<14 reads; 6183 OTUs)",
            "0.05% (<70 reads; 3467 OTUs)", "0.1% (<140 reads; 2710 OTUs)",
            "1% (<1399 reads; 979 OTUs)"),
          ncol = 3, nrow = 3, common.legend = TRUE,
          hjust = 0, label.x = 0.2, font.label = list(size = 10))

# With this analysis, we can see that removing rare OTUs reduces the variance
# within treatment a little bit, but does not affect overall patterns.

# Overall, all three plots showed that rare OTUs were mostly found in samples
# with large library size, and filtering out rare OTUs does not affect patterns.
# Only a threshold as high as 1% considerably affected patterns.

# Select the best filtering threshold
otu_table_rare <- otu_table_rare005


# Save cleaned OTU table -------------------------------------------------------

# OTU table with taxonomy filtered
save(otu_table_tax, file = "RData/16S/otu_taxa_clean.RData")

# OTU table with singletons filtered
save(otu_table_single, file = "RData/16S/otu_single_clean.RData")

# OTU table with rare taxa filtered
save(otu_table_rare, file = "RData/16S/otu_rare_clean.RData")

