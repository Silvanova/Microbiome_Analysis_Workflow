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
load("RData/16S/otu_rare_clean.RData")

# BLAST filtering --------------------------------------------------------------

# Filter out remaining artefacts: short sequences or (e.g., qlen < 250)
# or OTUs with a query cover lower than 15% (qcovs < 15).

# Extract BLAST results from OTU table
BLAST_rare <- as.data.frame(otu_table_rare@tax_table) %>%
  # Convert column class automatically
  type.convert(as.is = TRUE)

# Identify synthetic sequences based on BLAST results
otu_qlen <- filter(BLAST_rare, qlen < 250)[["taxa_OTU"]] # no sequence identified
otu_qcovs <- filter(BLAST_rare, qcovs < 15)[["taxa_OTU"]] # no sequence identified

# Filter synthetic sequences
BLAST_filter <- filter(BLAST_rare, !taxa_OTU %in% c(otu_qlen, otu_qcovs))
otu_BLAST_filter <- remove_taxa(c(otu_qlen, otu_qcovs), otu_table_rare)

# Count OTUs before and after filtering
ntaxa(otu_table_rare) # 8064
ntaxa(otu_BLAST_filter) # 8064

# In this demonstration, all artefacts were removed in prior steps
# (bioinformatics processing or singletons/rare filtering)

# BLAST cleaning ---------------------------------------------------------------

## Unknown organisms ####

# One last step of filtering can be done based on BLAST accuracy:
# 1. OTUs with a pident < 70% cannot be identified to the kingdom level
# 2. OTUs with a evalue > e-20 cannot be identified to the kingdom level
# 3. OTUs with a evalue between e-20 and e-50 are not certain
# 4. OTUs with a qcovs < 50% are not certain

# Identify unknown OTUs
BLAST_unknown <- filter(BLAST_filter, pident < 70 | evalue > 1e-50 | qcovs < 50)
View(BLAST_unknown)

# These OTUs should be investigated. A manual check against the 10 best BLAST
# matches allow for accurate assignment https://blast.ncbi.nlm.nih.gov/Blast.cgi

# After manual BLAST of these 44 OTUs, 8 could not be confidently
# identified as Bacteria. Before removing them, let's inspect them.

# Identify OTUs to remove after manual BLASTn
otu_BLAST_remove <- c("1d27c83d9692aab8c80bf374e70e6bd36ab107db",
                      "97fd840ad5b43b3e40763cb0f286c7f12482b5aa",
                      "a247774bebc7f5d0d55919174fa28717fde309c0",
                      "b32bc6dc7e56daabd7cbb16d6e64d75f5bbf1be7",
                      "c2273056078eeba2d2a2efc1cb351e5290d7bb91",
                      "fc6386fb30284b857e60f04b2101e43d368f60a1",
                      "501c4d4b723af4468417d493ca4fadbb1d6bfe4c",
                      "8767164c7b5703151beba4829286c69a6124f4e4")

# Identify OTUs that should only be assigned to the kingdom level
otu_BLAST_assign <- setdiff(BLAST_unknown[["taxa_OTU"]], otu_BLAST_remove)

# Filter OTU table to keep manually BLASTed OTUs to remove
otu_table_BLAST_remove <- prune_taxa(otu_BLAST_remove, otu_BLAST_filter)

# Visualize the abundance and prevalence of these OTUs across samples
plot_bar(otu_table_BLAST_remove, fill = "OTU")

# All OTUs identified are both scarce (low prevalence) and rare (low abundance),
# thus they can be removed.

# Remove remaining non-bacterial OTUs
BLAST_filter2 <- filter(BLAST_filter, !taxa_OTU %in% otu_BLAST_remove)
otu_BLAST_filter2 <- remove_taxa(otu_BLAST_remove, otu_BLAST_filter)

## Rank-level confidence ####

# Now that we are confident about our OTU selection, it is time to clean the
# taxonomy assignation. Rank-level confidence is based on pident.
# See Table 3 online ()
# for a summary of thresholds for each target gene and taxonomic rank found in the literature.
# Any taxonomic assignation below these thresholds is changed to NA.
# Uncultured and metagenome assignations are also changed to NA.

# Clean the OTU table based on taxonomy thresholds
# (remove assignation if confidence is too low).
BLAST_clean <- BLAST_filter2 %>%
  mutate(species = if_else(pident < 98.7, "NA", species),
         genus = if_else(pident < 94.5, "NA", genus),
         family = if_else(pident < 86.5, "NA", family),
         order = if_else(pident < 82.0, "NA", order),
         class = if_else(pident < 78.5, "NA", class),
         phylum = if_else(pident < 75.0, "NA", phylum)) %>%
  # Remove uncultured assignation
  mutate(across(phylum:species, ~ gsub(".*uncultured.*", "NA", .x))) %>%
  # Remove metagenome assignation
  mutate(across(phylum:species, ~ gsub(".*metagenome.*", "NA", .x))) %>%
  # Remove assignation above kingdom for previously identified OTUs
  mutate(across(phylum:species, ~ if_else(taxa_OTU %in% otu_BLAST_assign, "NA", .x)))

# Replace taxonomy table in the phyloseq object
otu_BLAST_clean <- otu_BLAST_filter2
otu_BLAST_clean@tax_table <- tax_table(as.matrix(BLAST_clean))

# Save cleaned OTU table -------------------------------------------------------

# OTU table with cleaned taxonomy, based on BLAST accuracy
save(otu_BLAST_clean, file = "RData/16S/otu_BLAST_clean.RData")

