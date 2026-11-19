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
load("RData/ITS/otu_table_raw.RData")

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
p_libsize <- ggplot(data = seq_depth, aes(x = Index, y = LibrarySize, color = sample_type)) +
  geom_point()
p_libsize

# Save the plot in the figures folder. ggsave() saves a ggplot object; the file
# type is given by the extension (.png, .pdf, ...) and the size is in inches.
ggsave("figures/003_neg_library_size.png", p_libsize, width = 7, height = 5)

# Zoom in on the lower part of the plot
p_libsize_zoom <- ggplot(data = seq_depth, aes(x = Index, y = LibrarySize, color = sample_type)) +
  geom_point() + ylim(0,2000) + xlim(0,12)
p_libsize_zoom
ggsave("figures/003_neg_library_size_zoom.png", p_libsize_zoom, width = 7, height = 5)

# We can see that, compared to the rest of the samples, all controls have a
# small library size. The negative controls contain very few reads (1 to 19
# reads), but they are not empty. Therefore, further analyses are required to
# know if these OTUs are contaminants.

# Furthermore, two samples have a very small library size compared to other
# samples (fewer reads than the positive controls). These samples should be
# further inspected to determine if their fungal community resembles that of
# other samples from the same site (i.e., replicates). If not, they should be
# discarded.

### Inspect small library size samples ####

# Identify samples with a small library size (SNF4__29363 and SNF4__29347)
subset(seq_depth, sample_type == "sample" & LibrarySize < 1000)

# Add a column to sample data in phyloseq object to identify these samples
sample_data(otu_table_raw)$low_depth <-
  sample_names(otu_table_raw) %in% c("SNF4__29363", "SNF4__29347")

# Plot an ordination of all samples to compare the low depth samples to replicates
p_low_depth_pcoa <- otu_table_raw %>%
  tax_transform("identity", rank = "taxa_OTU") %>%
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "SiteID", shape = "low_depth", auto_caption = NA) +
  stat_ellipse(aes(colour = SiteID))
p_low_depth_pcoa
ggsave("figures/003_neg_low_depth_sample_pcoa.png", p_low_depth_pcoa, width = 7, height = 5)

# These samples (shown with triangles) are located among the other samples of
# their site (Oak_20 and Oak_33) and the other oak forests. Their community
# composition is therefore similar to that of their replicates, and there is no
# reason to remove them at this stage. Note that the controls (NA, in gray) are
# located close to the oak forest samples: the few reads in the controls most
# likely come from the samples (see below).

# However, their library size is very small: they will most likely be removed
# later, when choosing a rarefaction depth (see 006_normalization.R).

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
p_prev <- ggplot(data=otu_table.pa, aes(x=pa.neg, y=pa.sam)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
p_prev
ggsave("figures/003_neg_prevalence_controls_vs_samples.png", p_prev, width = 7, height = 5)
# 37 OTUs are present in at least one negative control
# ('Prevalence (Negative Controls)' >= 1). Most of them are present in only
# one control, but in many true samples.

# We must find out if there was a contamination in our samples from the method,
# or if the negative samples have been contaminated by the other samples.

# For further analyses, convert OTU table from count data to relative abundance
otu_table_rel <- transform_sample_counts(otu_table_raw.noPos, function(x) x/sum(x))

# Identify OTUs present in negative controls
otu_neg <- filter(otu_table.pa, pa.neg > 0) %>% rownames()

# Filter OTU table to keep only OTUs found in negative controls
otu_table_neg <- prune_taxa(otu_neg, otu_table_rel)

# Visualize the abundance of these OTUs in all samples
p_neg_rel <- plot_bar(otu_table_neg, fill = "phylum") +
  facet_grid(~sample_type, scales = "free", space = "free")
p_neg_rel
ggsave("figures/003_neg_control_otus_relative_abundance.png", p_neg_rel,
       width = 12, height = 6)

# A commonly recommended practice with negative controls is to remove any
# OTU found in negative controls from the rest of the data set. However, as
# you can see in this demonstration, the OTUs found in the negative controls
# are found in all samples, and represent on average ~10% (and up to ~37%) of
# the community in the samples. Removing them would remove a large part of the
# fungal community. It is most likely that the OTUs found in the negative
# controls are present due to a contamination from the samples (e.g., during
# pipetting) and not the other way around.

# Another way to confirm this, is to look at the count data instead of the
# relative abundance:

# Filter OTU table to keep only OTUs found in negative controls
otu_table_neg <- prune_taxa(otu_neg, otu_table_raw.noPos)

# Visualize the abundance of these OTUs in all samples
p_neg_counts <- plot_bar(otu_table_neg, fill = "phylum") +
  facet_grid(~sample_type, scales = "free", space = "free")
p_neg_counts
ggsave("figures/003_neg_control_otus_counts.png", p_neg_counts,
       width = 12, height = 6)

# Here we can see that the abundance of these OTUs in the negative controls
# (1 to 19 reads) is very small compared to true samples (often hundreds of
# reads).

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

contam.prev01 <- isContaminant(otu_table_raw.noPos, method="prevalence",
                               neg="is.neg", threshold = 0.1)

table(contam.prev01$contaminant)

# The method found 4 contaminant sequences at a threshold of 0.1 (default).

# Increase the threshold to try and find more contaminants. This more aggressive
# threshold will identify as contaminants all sequences that are more
# prevalent in negative controls than in true samples
contam.prev05 <- isContaminant(otu_table_raw.noPos, method="prevalence",
                               neg="is.neg", threshold = 0.5)

table(contam.prev05$contaminant)

# The method found 17 contaminant sequences at threshold = 0.5.

## Inspect contaminants ####

# Visualize the prevalence of these potential contaminants in samples and
# controls using the presence-absence table created earlier
otu_table.pa <- mutate(otu_table.pa,
                       contaminant = case_when(
                         contam.prev01$contaminant ~ "threshold 0.1",
                         contam.prev05$contaminant ~ "threshold 0.5 only",
                         .default = "not a contaminant"))

# Plot the prevalence of OTU in negative controls vs true samples
p_prev_contam <- ggplot(data=otu_table.pa, aes(x=pa.neg, y=pa.sam, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
p_prev_contam
ggsave("figures/003_neg_prevalence_contaminants.png", p_prev_contam, width = 7, height = 5)

# The negative controls contain so few reads that every OTU is present in at
# most 2 controls. With the threshold of 0.5, decontam also flags OTUs that are
# present in many true samples. Let's look at the taxonomy of the flagged OTUs
# and at their number of reads in the true samples.

# Number of reads of each OTU in the true samples
otu_table_sam <- prune_samples(
  sample_data(otu_table_raw.noPos)$sample_type == "sample", otu_table_raw.noPos)

contam_check <- data.frame(
  tax_table(otu_table_raw.noPos)[, c("kingdom", "phylum", "genus")],
  reads_samples = taxa_sums(otu_table_sam),
  prevalence_samples = otu_table.pa$pa.sam,
  prevalence_controls = otu_table.pa$pa.neg,
  contaminant = otu_table.pa$contaminant) %>%
  filter(contaminant != "not a contaminant") %>%
  arrange(contaminant, desc(reads_samples))
contam_check

# The 4 OTUs flagged with the default threshold (0.1) have very few reads in
# the true samples (1 to 11 reads). One of them is Malassezia, a yeast living
# on human skin and a well-known laboratory contaminant. These are most likely
# true contaminants.

# On the other hand, many of the additional OTUs flagged with the threshold of
# 0.5 are abundant in the true samples (up to ~600 reads), and some belong to
# typical soil and forest fungi, such as the ectomycorrhizal genera Cortinarius,
# Tomentella and Sebacina. These OTUs were most likely carried over from the
# samples into the negative controls, and removing them would remove real
# members of the fungal community. We therefore use the default threshold.

contaminants <- filter(contam.prev01, contaminant == TRUE) %>% rownames()

# Filter OTU table to keep only OTUs suspected of being contaminants
otu_table_contam <- prune_taxa(contaminants, otu_table_rel)

# Visualize the abundance of these OTUs in all samples
p_contam <- plot_bar(otu_table_contam, fill = "OTU") +
  facet_grid(~sample_type, scales = "free", space = "free") +
  theme(legend.position = "none")
p_contam
ggsave("figures/003_neg_contaminants_relative_abundance.png", p_contam,
       width = 12, height = 6)

# These OTUs make up a large part of the reads in the negative controls, and a
# negligible part of the reads in the true samples. It is safe to consider them
# as contaminants and remove them from the data set. The other OTUs present in
# the negative controls were most likely coming from true samples and are
# therefore not considered contaminants.

# Note: with low-biomass controls such as these (few reads), decontam has
# little information to work with. Always inspect the flagged OTUs before
# removing them, as we did here.

# Clean and save OTU table -----------------------------------------------------

# The two samples with a low sequencing depth (SNF4__29363 and SNF4__29347)
# are kept for now, as their community resembles that of their replicates.

# Remove contaminants (note that remove_taxa is a homemade function)
otu_ctrl_clean <- remove_taxa(contaminants, otu_table_raw)

# Remove all controls
otu_ctrl_clean <- prune_samples(
  sample_data(otu_ctrl_clean)$sample_type != "positive control", otu_ctrl_clean)

otu_ctrl_clean <- prune_samples(
  sample_data(otu_ctrl_clean)$sample_type != "negative control", otu_ctrl_clean)

# Save cleaned OTU table
save(otu_ctrl_clean, file = "RData/ITS/otu_ctrl_clean.RData")
