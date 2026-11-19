#==============================================================================#
#                                                                              #
#  Dealing with positive controls (synthetic community)                        #
#                                                                              #
#  Modified from Karelle Rheault: karh@ign.ku.dk                               #
#  Modified by Davide Barsotti: db@ign.ku.dk                                   #
#                                                                              #
#==============================================================================#

# Set working directory (modify the path below to match your own)
setwd("some_folder_path/my_working_directory")

# Load packages and homemade functions
source("scripts/001_required_packages.R")

# The figures of the ITS scripts are saved in a separate folder
dir.create("figures/ITS", showWarnings = FALSE)

# Load raw OTU table
load("RData/ITS/otu_table_raw.RData")

# Inspect samples --------------------------------------------------------------

# Convert OTU table from count data to presence-absence
otu_table.pa <- transform_sample_counts(otu_table_raw, function(x) 1*(x>0))

# Convert OTU table from count data to relative abundance
otu_table.rel <- transform_sample_counts(otu_table_raw, function(x) x/sum(x))

## Identify synthetic community ####

# Identify OTUs present in positive controls
otu_table.pa.pos <- prune_samples(
  sample_data(otu_table.pa)$sample_type == "positive control", otu_table.pa)

table.pa.pos <- data.frame(pa.pos = taxa_sums(otu_table.pa.pos))

otu_pos <- filter(table.pa.pos, pa.pos > 0) %>% rownames()

# Filter OTU table to keep only OTUs found in positive controls
otu_table_pos <- prune_taxa(otu_pos, otu_table.rel)

# Visualize the abundance of these OTUs in all samples
p_pos_otus <- plot_bar(otu_table_pos, fill = "kingdom") +
  facet_grid(~sample_type, scales = "free", space = "free")
p_pos_otus
ggsave("figures/ITS/003_pos_control_otus_relative_abundance.png", p_pos_otus,
       width = 12, height = 6)

## Compare community profiles ####

# If one suspects issues with sequencing, positive controls can be used for
# visual inspection, by making sure the taxonomic profiles match expectations.

# Keep only positive controls for visualization
otu_table.rel.pos <- prune_samples(
  sample_data(otu_table.rel)$sample_type == "positive control", otu_table.rel)

# Visualize the community profiles of these samples
p_pos_profiles <- plot_bar(otu_table.rel.pos, fill = "kingdom")
p_pos_profiles
ggsave("figures/ITS/003_pos_control_profiles.png", p_pos_profiles, width = 7, height = 5)

# Clean and save OTU table -----------------------------------------------------

# In this example, the positive controls contain a synthetic ITS sequence
# (a "spike-in", assigned to the kingdom "Unispike1"). Nearly all reads in the
# positive controls belong to the spike-in, as expected, plus a few reads (1-2)
# of soil fungi that most likely come from the true samples. The spike-in is
# absent from the true samples, so there was no contamination from the positive
# controls into the samples.

# No obvious issues were detected in the positive controls, and the spike-in
# will be automatically removed when performing taxonomic filtering (keeping
# only Fungi, see 004_global_filter.R). Therefore, positive controls can simply
# be removed.

# Remove positive controls
otu_ctrl_pos_clean <- prune_samples(
  sample_data(otu_table_raw)$sample_type != "positive control", otu_table_raw)

# Save cleaned OTU table
save(otu_ctrl_pos_clean, file = "RData/ITS/otu_ctrl_pos_clean.RData")


