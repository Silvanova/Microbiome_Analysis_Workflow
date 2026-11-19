#==============================================================================#
#                                                                              #
#  Explore the data: OTU table, taxonomy table and metadata                    #
#                                                                              #
#  Modified from Sofia Gomes: s.i.gomes@biology.leidenuniv.nl                  #
#  Modified by Davide Barsotti: db@ign.ku.dk                                   #
#                                                                              #
#==============================================================================#

# Set working directory (modify the path below to match your own)
setwd("some_folder_path/my_working_directory")

# Load packages and homemade functions
source("scripts/001_required_packages.R")
dir.create("figures/ITS", showWarnings = FALSE, recursive = TRUE)

# Load files
load("RData/ITS/otu_table_raw.RData")
load("RData/auxfile.RData")

# Before cleaning the data, we take a first look at what we have. A phyloseq
# object combines the three tables we need to start working:
#   - the OTU table (number of reads of each OTU in each sample)
#   - the taxonomy table (taxonomic assignation and BLAST results of each OTU)
#   - the sample data, or metadata (information about each sample)
otu_table_raw


# Inspect the three tables -----------------------------------------------------

## OTU table ####

# OTUs are in rows and samples in columns. Most values are 0: each OTU is only
# found in a few samples.
otu_table(otu_table_raw)[1:10, 1:3]

## Taxonomy table ####

# The first columns contain the taxonomic assignation, from kingdom to species.
tax_table(otu_table_raw)[1:3, 1:7]

# The other columns contain the BLAST results for each OTU
colnames(tax_table(otu_table_raw))

# Meaning of the BLAST columns:
#   qlen:     query sequence length (length of the OTU sequence)
#   slen:     subject sequence length (length of the matched reference sequence)
#   qstart:   start position of the alignment on the query sequence
#   qend:     end position of the alignment on the query sequence
#   sstart:   start position of the alignment on the subject (reference) sequence
#   send:     end position of the alignment on the subject (reference) sequence
#   evalue:   expect value. The number of hits expected by chance in a database
#             of this size (lower values mean more significant matches)
#   length:   alignment length
#   nident:   number of identical matches in the alignment
#   mismatch: number of mismatches in the alignment
#   gapopen:  number of gap openings
#   gaps:     total number of gaps (gap openings plus gap extensions)
#   sstrand:  strand of the subject sequence (plus or minus)
#   qcovs:    query coverage per subject (percentage of the query sequence
#             covered by the alignment)
#   pident:   percentage of identical matches (nident relative to length)
# These columns are used later to check the quality of the taxonomic
# assignations (script 005_BLAST_filter.R).

## Metadata ####

# The sample data contains all the information associated with the samples
# that we can use to answer ecological questions: site, land use, age of the
# plantation, and soil physicochemical properties.
head(sample_data(otu_table_raw))


# Match samples and OTUs across tables -----------------------------------------

# When working with separate tables, we must make sure that the OTUs are in the
# same order in the OTU table and in the taxonomy table, and that the samples
# are in the same order in the OTU table and in the metadata. phyloseq does
# this for us: when the object is created (script 002_data_preparation.R), the
# tables are matched by their names, and samples or OTUs missing from one of
# the tables are dropped.

# Check that the tables match
identical(taxa_names(otu_table(otu_table_raw)), taxa_names(tax_table(otu_table_raw)))
identical(sample_names(otu_table(otu_table_raw)), sample_names(sample_data(otu_table_raw)))

# Samples of the auxiliary data that are not in the phyloseq object (they were
# not sequenced, or were removed during bioinformatic processing)
setdiff(aux$Library_ID, sample_names(otu_table_raw))

# The bioinformatics pipeline gave a long string (a hash of the sequence) as
# name of each OTU. We can replace them with shorter names. phyloseq renames
# the OTUs in all tables at once. Here we do it on a copy of the object only,
# to keep the original names, which are used in the next scripts.
otu_renamed <- otu_table_raw
taxa_names(otu_renamed) <- paste0("OTU", seq_len(ntaxa(otu_renamed)))
tax_table(otu_renamed)[1:6, 1:7]


# Let's look at the data! ------------------------------------------------------

## Samples and controls ####

# How many samples and how many OTUs?
nsamples(otu_table_raw) # 97
ntaxa(otu_table_raw)    # 16691

# The data set contains both samples and controls (see "Dealing with controls")
table(sample_data(otu_table_raw)$sample_type)

# Keep only the samples to explore the data
otu_samples <- subset_samples(otu_table_raw, sample_type == "sample")

# Remove the OTUs only found in the controls
otu_samples <- prune_taxa(taxa_sums(otu_samples) > 0, otu_samples)

nsamples(otu_samples) # 89
ntaxa(otu_samples)    # 16681

## Reads ####

# How many reads in total in the data set?
sum(sample_sums(otu_samples)) # 429715

# How many reads per sample (library size, or sequencing depth)?
sample_sums(otu_samples)
summary(sample_sums(otu_samples))

# Visualize the distribution of library sizes
p_libsize <- data.frame(LibrarySize = sample_sums(otu_samples)) %>%
  ggplot(aes(x = LibrarySize)) +
  geom_histogram(bins = 30, fill = "grey60", color = "black") +
  labs(x = "Number of reads per sample", y = "Number of samples")
p_libsize
ggsave("figures/ITS/002b_library_size.png", p_libsize, width = 6, height = 4)

# Library sizes vary a lot between samples, from 135 to 11053 reads. This must
# be taken into account before comparing samples (see "Rarefaction and
# normalization").

## Metadata ####

# Extract the metadata as a data frame to work with dplyr
meta <- data.frame(sample_data(otu_samples))

# How many samples per land use?
table(meta$Land_use)

# What is the age distribution of the sites?
table(meta$Age)

# Summaries of multiple variables together
meta %>% dplyr::count(Land_use, SiteID, Age)

# The data set contains oak plantations of different ages (20 to 54 years),
# old oak forest (231 years, reference), and agricultural fields and a
# grassland (age 0 and 44).

## Taxonomy ####

# Extract the taxonomy table as a data frame
tax <- data.frame(tax_table(otu_samples))

# How many OTUs per kingdom?
tax %>% dplyr::count(kingdom, sort = TRUE)

# The ITS primers also amplify many organisms other than fungi (protists,
# animals, plants). Note also that some names are written in different ways
# (e.g., "Viridiplantae" and "viridiplantae"). Only the fungi are kept in the
# next steps (see "Global filtering").

# Visualize the number of OTUs per kingdom
p_kingdom <- tax %>%
  dplyr::count(kingdom) %>%
  ggplot(aes(x = n, y = reorder(kingdom, n))) +
  geom_col() +
  labs(x = "Number of OTUs", y = "Kingdom")
p_kingdom
ggsave("figures/ITS/002b_otus_per_kingdom.png", p_kingdom, width = 6, height = 5)

# Keep only fungi to look at the fungal orders
tax_fungi <- tax %>% filter(kingdom == "Fungi")
nrow(tax_fungi) # 9042

# How many fungal OTUs per order?
tax_fungi %>% dplyr::count(order, sort = TRUE)

# Many fungal OTUs are not assigned to an order: the reference databases are
# still incomplete for soil fungi.
