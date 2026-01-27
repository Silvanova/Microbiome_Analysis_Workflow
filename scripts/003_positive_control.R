#==============================================================================#
#                                                                              #
#  Dealing with positive controls (synthetic community)                        #
#                                                                              #
#  Written by Karelle Rheault: karh@ign.ku.dk                                  #
#                                                                              #
#==============================================================================#

# Set working directory (modify the path below to match your own)
setwd("some_folder_path/my_working_directory")

# Load packages
source("scripts/001_required_packages.R")

# Either load the cleaned OTU table from previous script "003_negative_control.R"
load("RData/16S/otu_ctrl_neg_clean.RData")

# Or load raw OTU table if there were no negative controls to deal with
load("RData/16S/otu_table_raw.RData")

# Inspect samples --------------------------------------------------------------

# Convert OTU table from count data to presence-absence
otu_table.pa <- transform_sample_counts(otu_ctrl_neg_clean, function(x) 1*(x>0))

# Convert OTU table from count data to relative abundance
otu_table.rel <- transform_sample_counts(otu_ctrl_neg_clean, function(x) x/sum(x))

## Identify synthetic community ####

# Identify OTUs present in positive controls
otu_table.pa.pos <- prune_samples(
  sample_data(otu_table.pa)$sample_type == "positive control", otu_table.pa)

table.pa.pos <- data.frame(pa.pos = taxa_sums(otu_table.pa.pos))

otu_pos <- filter(table.pa.pos, pa.pos > 0) %>% rownames()

# Filter OTU table to keep only OTUs found in positive controls
otu_table_pos <- prune_taxa(otu_pos, otu_table.rel)

# Visualize the abundance of these OTUs in all samples
plot_bar(otu_table_pos, fill = "kingdom") +
  facet_grid(~sample_type, scales = "free", space = "free")

## Compare community profiles ####

# If one suspects issues with sequencing, positive controls can be used for
# visual inspection, by making sure the taxonomic profiles match expectations.

# Keep only positive controls for visualization
otu_table.rel.pos <- prune_samples(
  sample_data(otu_table.rel)$sample_type == "positive control", otu_table.rel)

# Remove taxa with less than 1% abundance (for clarity of the plot)
otu_table.rel.pos <- filter_taxa(otu_table.rel.pos, function(x) sum(x) > .01, TRUE)

# Visualize the community profiles of these samples
plot_bar(otu_table.rel.pos, fill = "phylum")

# Clean and save OTU table -----------------------------------------------------

# In this example, no obvious issues were detected in positive controls and
# the synthetic community will be automatically removed when performing
# taxonomic filtering. Therefore, positive controls can simply be removed.

# Remove positive controls
otu_ctrl_pos_clean <- prune_samples(
  sample_data(otu_ctrl_neg_clean)$sample_type != "positive control", otu_ctrl_neg_clean)

# Save cleaned OTU table
save(otu_ctrl_pos_clean, file = "RData/16S/otu_ctrl_pos_clean.RData")


