#==============================================================================#
#                                                                              #
#  BLAST quality check: cleaning OTU table based on BLAST accuracy             #
#                                                                              #
#  Written by Karelle Rheault: karh@ign.ku.dk                                  #
#                                                                              #
#==============================================================================#

# Set working directory (modify the path below to match your own)
setwd("some_folder_path/my_working_directory")

# Load packages and homemade functions
source("scripts/001_required_packages.R")

# Load the cleaned OTU table from previous scripts
load("RData/ITS/otu_rare_clean.RData")

# BLAST filtering --------------------------------------------------------------

# Filter out remaining artefacts: short sequences (for ITS, qlen < 100) or OTUs
# with a query cover lower than 15% (qcovs < 15).

# Extract BLAST results from OTU table
BLAST_rare <- as.data.frame(otu_table_rare@tax_table) %>%
  # Convert column class automatically
  type.convert(as.is = TRUE)

# Identify artefacts based on BLAST results
otu_qlen <- filter(BLAST_rare, qlen < 100)[["taxa_OTU"]] # no sequence identified
otu_qcovs <- filter(BLAST_rare, qcovs < 15)[["taxa_OTU"]] # 3 sequences identified

# Inspect the sequences with a low query cover
filter(BLAST_rare, taxa_OTU %in% otu_qcovs) %>%
  select(phylum, genus, qlen, qcovs, pident, evalue)
taxa_sums(otu_table_rare)[otu_qcovs]

# These 3 OTUs are much longer (1147 to 1638 bp) than the expected length of
# the amplicon (~280 bp; see the table of primer sets on the website), only a
# small part of them matches a reference sequence, and they are represented by
# only 1 or 2 reads. They are most likely artefacts (e.g., chimeras) and can be
# removed.

# Filter artefacts
BLAST_filter <- filter(BLAST_rare, !taxa_OTU %in% c(otu_qlen, otu_qcovs))
otu_BLAST_filter <- remove_taxa(c(otu_qlen, otu_qcovs), otu_table_rare)

# Count OTUs before and after filtering
ntaxa(otu_table_rare) # 9039
ntaxa(otu_BLAST_filter) # 9036

# In this demonstration, most artefacts were removed in prior steps
# (bioinformatic processing).

# BLAST cleaning ---------------------------------------------------------------

## Unknown organisms ####

# One last step of filtering can be done based on BLAST accuracy:
# 1. OTUs with a pident < 70% cannot be identified to the kingdom level
# 2. OTUs with a evalue > e-20 cannot be identified to the kingdom level
# 3. OTUs with a evalue between e-20 and e-50 are not certain
# 4. OTUs with a qcovs < 50% are not certain

# Identify unknown OTUs
BLAST_unknown <- filter(BLAST_filter, pident < 70 | evalue > 1e-50 | qcovs < 50)
nrow(BLAST_unknown) # 1078 OTUs

# How many reads do they represent?
sum(taxa_sums(otu_BLAST_filter)[BLAST_unknown$taxa_OTU]) # 5474 reads (~1.5%)

# To which phyla were they assigned?
sort(table(BLAST_unknown$phylum), decreasing = TRUE)

# In a smaller dataset, these OTUs could be investigated one by one: a manual
# check against the 10 best BLAST matches allows for an accurate assignment,
# using the same database as the one used to assign the taxonomy of the OTU
# table (for ITS: UNITE, https://unite.ut.ee/analysis.php). OTUs that cannot be
# confidently identified as Fungi could then be removed.

# Here, more than 1000 OTUs fall in this category, which makes a manual check
# impractical. Many of them belong to groups of fungi that are poorly
# represented in reference databases, such as Rozellomycota and Zoopagomycota,
# which explains why their BLAST matches are weak. Together, they represent
# only ~1.5% of the reads. We therefore keep these OTUs, but we do not trust
# their taxonomy below the kingdom level: their assignation at all ranks below
# kingdom will be removed (see below).

# Identify OTUs that should only be assigned to the kingdom level
otu_BLAST_assign <- BLAST_unknown[["taxa_OTU"]]

# Visualize the proportion of reads belonging to these OTUs in each sample
unknown_prop <- data.frame(
  sample_data(otu_BLAST_filter),
  prop_unknown = sample_sums(prune_taxa(otu_BLAST_assign, otu_BLAST_filter)) /
    sample_sums(otu_BLAST_filter) * 100)

p_unknown <- ggplot(unknown_prop, aes(x = SiteID, y = prop_unknown)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  labs(x = "Site", y = "Reads of OTUs with uncertain BLAST matches (%)") +
  theme_bw()
p_unknown
ggsave("figures/005_BLAST_unknown_otus.png", p_unknown, width = 7, height = 5)

# In all samples, these OTUs represent a small proportion of the reads (less
# than 5%).

## Rank-level confidence ####

# Now that we are confident about our OTU selection, it is time to clean the
# taxonomy assignation. Rank-level confidence is based on pident.
# See Table 3 online ()
# for a summary of thresholds for each target gene and taxonomic rank found in the literature.
# Any taxonomic assignation below these thresholds is changed to NA.

# For ITS, we use the rules of thumb of Tedersoo et al. (2014, 2019, 2021):
# species 98%, genus 90%, family 85%, order 80%, class 75% and phylum 75%.

# In the ITS reference database, ranks without assignation are written as ".".
# Replace them with "NA", so that all missing assignations are written the
# same way.

# Clean the OTU table based on taxonomy thresholds
# (remove assignation if confidence is too low).
BLAST_clean <- BLAST_filter %>%
  # Replace missing assignations (".") with "NA"
  mutate(across(phylum:species, ~ if_else(.x == ".", "NA", .x))) %>%
  mutate(species = if_else(pident < 98, "NA", species),
         genus = if_else(pident < 90, "NA", genus),
         family = if_else(pident < 85, "NA", family),
         order = if_else(pident < 80, "NA", order),
         class = if_else(pident < 75, "NA", class),
         phylum = if_else(pident < 75, "NA", phylum)) %>%
  # Remove uncultured assignation
  mutate(across(phylum:species, ~ gsub(".*uncultured.*", "NA", .x))) %>%
  # Remove metagenome assignation
  mutate(across(phylum:species, ~ gsub(".*metagenome.*", "NA", .x))) %>%
  # Remove assignation below kingdom for previously identified OTUs
  mutate(across(phylum:species, ~ if_else(taxa_OTU %in% otu_BLAST_assign, "NA", .x)))

# How many OTUs are assigned to each rank, before and after cleaning?
colSums(BLAST_filter[, c("phylum", "class", "order", "family", "genus", "species")] != ".")
colSums(BLAST_clean[, c("phylum", "class", "order", "family", "genus", "species")] != "NA")

# Replace taxonomy table in the phyloseq object
otu_BLAST_clean <- otu_BLAST_filter
otu_BLAST_clean@tax_table <- tax_table(as.matrix(BLAST_clean))

# Save cleaned OTU table -------------------------------------------------------

# OTU table with cleaned taxonomy, based on BLAST accuracy
save(otu_BLAST_clean, file = "RData/ITS/otu_BLAST_clean.RData")
