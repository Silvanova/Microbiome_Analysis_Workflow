#==============================================================================#
#                                                                              #
#  Dealing with negative controls (and blanks)                                 #
#                                                                              #
#  Written by Karelle Rheault: karh@ign.ku.dk                                  #
#                                                                              #
#==============================================================================#

# Set working directory (modify the path below to match your own)
setwd("some_folder_path/my_working_directory")

# Load packages
source("scripts/001_required_packages.R")

# Load files
load("RData/16S/otu_table_raw.RData")

# Before using the package decontam to remove contaminants, we must first
# inspect samples (and controls) to determine if there are any visible problems.

# Pre inspection of samples ----------------------------------------------------

## Library size ----------------------------------------------------------------

# Inspect library size (sequencing depth)
# Negative controls and blanks should have few sequences
seq_depth <- cbind.data.frame(sample_data(otu_table_raw)) %>%
  mutate(LibrarySize = sample_sums(otu_table_raw)) %>%
  # Sort the data set by LirarySize and create an Index to match the sorted data
  arrange(LibrarySize) %>% mutate(Index = row_number())

# Visualize LibrarySize as a function of Index
ggplot(data = seq_depth, aes(x = Index, y = LibrarySize, color = sample_type)) +
  geom_point()

# Zoom in on the lower part of the plot
ggplot(data = seq_depth, aes(x = Index, y = LibrarySize, color = sample_type)) +
  geom_point() + ylim(0,27000) + xlim(0,11)

# We can see that, compared to the rest of the samples, all controls have a
# small library size. Nevertheless, all negative controls have a library size
# larger than zero. Therefore, further analyses are required to know if these
# OTUs are contaminants.

# Furthermore, one sample has a very small library size compared to other
# samples (and compared to controls). This sample should be further inspected
# to determine if the microbial community in this sample resemble that of other
# samples (i.e., replicates). If not, it should be discarded.

### Inspect small library size sample ####

# Identify sample with a small library size (X29239)
subset(seq_depth, sample_type == "sample" & LibrarySize < 10000)

# Add a column to sample data in phyloseq object to identify this sample
sample_data(otu_table_raw)$low_depth <- sample_data(otu_table_raw)$Soil_ID=="X29239"

# Plot an ordination of all samples to compare the low seq sample to replicates
otu_table_raw %>%
  tax_transform("identity", rank = "taxa_OTU") %>%
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "SiteID", shape = "low_depth", auto_caption = NA) +
  stat_ellipse(aes(colour = SiteID))

# This sample (shown with a triangle) is located far away from samples of the
# same site (Oak_51, in light blue), and close to the controls (NA, in gray)
# demonstrating that it harbors a different community composition than samples
# from the same site. If this difference cannot be explained, this sample
# should be removed.

## Prevalence of OTUs in negative controls -------------------------------------

# For the following analyses, remove positive controls from the data set
otu_table_raw.noPos <- prune_samples(
  sample_data(otu_table_raw)$sample_type != "positive control", otu_table_raw)

# Visualize the prevalence of OTUs in samples and controls

# Make phyloseq objects of presence-absence in controls and true samples
otu_table.pa <- transform_sample_counts(otu_table_raw.noPos, function(x) 1*(x>0))
otu_table.pa.neg <- prune_samples(
  sample_data(otu_table.pa)$sample_type == "negative control", otu_table.pa)
otu_table.pa.sam <- prune_samples(
  sample_data(otu_table.pa)$sample_type == "sample", otu_table.pa)

# Make data.frame of prevalence in controls and true samples
otu_table.pa <- data.frame(pa.sam = taxa_sums(otu_table.pa.sam),
                           pa.neg = taxa_sums(otu_table.pa.neg))

# Plot the prevalence of OTU in negative controls vs true samples
ggplot(data=otu_table.pa, aes(x=pa.neg, y=pa.sam)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
# Many OTUs are present in at least one negative control
# ('Prevalence (Negative Controls)' >= 1).

# We must find out if there was a contamination in our samples from the method,
# or if the negative samples have been contaminated by the other samples.

# For further analyses, convert OTU table from count data to relative abundance
otu_table_rel <- transform_sample_counts(otu_table_raw.noPos, function(x) x/sum(x))

# Identify OTUs present in negative controls
otu_neg <- filter(otu_table.pa, pa.neg > 0) %>% rownames()

# Filter OTU table to keep only OTUs found in negative controls
otu_table_neg <- prune_taxa(otu_neg, otu_table_rel)

# Visualize the abundance of these OTUs in all samples
plot_bar(otu_table_neg, fill = "phylum") +
  facet_grid(~sample_type, scales = "free", space = "free")

# A commonly recommended practice with negative controls is to remove any
# OTU found in negative controls from the rest of the data set. However, as
# you can see in this demonstration, the OTUs found in the negative controls
# are found in all samples and represent ~25% of the community in all samples.
# Therefore, it is most likely that the OTUs found in the negative controls
# are present due to a contamination from the rest of the samples and not
# the other way around.

# Another way to confirm this, is to look at the count data instead of the
# relative abundance:

# Filter OTU table to keep only OTUs found in negative controls
otu_table_neg <- prune_taxa(otu_neg, otu_table_raw.noPos)

# Visualize the abundance of these OTUs in all samples
plot_bar(otu_table_neg, fill = "phylum") +
  facet_grid(~sample_type, scales = "free", space = "free")

# Here we can see that the abundance of these contaminants in the
# negative controls is very small compared to true samples.

# R package decontam -----------------------------------------------------------

# A way to confirm our observations is to use the package decontam, which
# uses simple statistical tools to identify contaminant automatically.

# Visit their website for more information.
# https://bioconductor.org/packages/devel/bioc/vignettes/decontam/inst/doc/decontam_intro.html

# Use the function isContaminant with the prevalence method
# to find which OTUs are contaminants

# First, add a TRUE/FALSE column "is.neg" to the sample_data
# to identify negative controls
sample_data(otu_table_raw.noPos)$is.neg <-
  sample_data(otu_table_raw.noPos)$sample_type == "negative control"

contam.prev <- isContaminant(otu_table_raw.noPos, method="prevalence",
                             neg="is.neg", threshold = 0.1)

table(contam.prev$contaminant)

# The method found 21 contaminant sequences at a threshold of 0.1 (default).

# Increase the threshold to try and find more contaminants. This more aggressive
# threshold will identify as contaminants all sequences that are are more
# prevalent in negative controls than in true samples
contam.prev <- isContaminant(otu_table_raw.noPos, method="prevalence",
                             neg="is.neg", threshold = 0.5)

table(contam.prev$contaminant)

# The method found 71 contaminant sequences at threshold = 0.5.

## Inspect contaminants ####

# Visualize the prevalence of these potential contaminants in samples and
# controls using the presence-absence table created earlier
otu_table.pa <- mutate(otu_table.pa, contaminant = contam.prev$contaminant)

# Plot the prevalence of OTU in negative controls vs true samples
ggplot(data=otu_table.pa, aes(x=pa.neg, y=pa.sam, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
# The method identified OTUs that were mostly present in negative controls
# and few true samples

# Filter OTU table to keep only OTUs suspected of being contaminants
contaminants <- filter(contam.prev, contaminant == TRUE) %>% rownames()
otu_table_contam <- prune_taxa(contaminants, otu_table_rel)

# Visualize the abundance of these OTUs in all samples
plot_bar(otu_table_contam, fill = "OTU") +
  facet_grid(~sample_type, scales = "free", space = "free") +
  theme(legend.position = "none")

# After inspection of these OTUs, it is safe to consider them as contaminants
# and remove them from the data set. The other OTUs present in the negative
# controls were most likely coming from true samples and are therefore not
# considered contaminants.

# Clean and save OTU table -----------------------------------------------------

# Remove samples with low sequencing depth (sampleID = X29239)
otu_table_raw <- prune_samples(
  sample_names(otu_table_raw) != "X29239", otu_table_raw)

# Remove contaminants (note that remove_taxa is a homemade function)
otu_ctrl_clean <- remove_taxa(contaminants, otu_table_raw)

# Remove all controls
otu_ctrl_clean <- prune_samples(
  sample_data(otu_ctrl_clean)$sample_type != "positive control", otu_ctrl_clean)

otu_ctrl_clean <- prune_samples(
  sample_data(otu_ctrl_clean)$sample_type != "negative control", otu_ctrl_clean)

# Save cleaned OTU table
save(otu_ctrl_clean, file = "RData/16S/otu_ctrl_clean.RData")
